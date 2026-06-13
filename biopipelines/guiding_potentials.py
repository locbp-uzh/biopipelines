# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Typed builders for RFdiffusion guiding potentials.

A plain contig constrains only what is *fixed*; it does not tell the model that
a newly generated region should *contact* a ligand, so the new segment often
drifts away from the bound molecule. Guiding potentials bias the diffusion
trajectory toward contact.

Hand-writing the hydra string
(``'["type:ligand_ncontacts,weight:3,r_0:8,d_0:4"]'``) is error-prone: the exact
``type:`` names, parameter keys, comma/colon syntax and bracket/quoting all have
to be right. These builders produce that string from named, validated arguments
with sensible defaults, so a tool can accept either a builder (or a list of
them) or the raw string.

Usage (via the nested alias on the tool)::

    from biopipelines import RFdiffusionAllAtom as RFAA

    RFAA(
        pdb=poses, ligand=Ligand(code="LIG"), contigs="40-60,A84-182",
        guiding_potentials=RFAA.GuidingPotential.ligand_ncontacts(
            weight=3, r_0=8, d_0=4),
        guide_decay="quadratic",
    )

``render_guiding_potentials(...)`` turns a builder / list / raw string into the
single hydra value the wrapper emits.
"""

from typing import Any, Dict, List, Optional, Union


def _fmt(v: Any) -> str:
    """Format a parameter value for the hydra token (drop trailing .0 on ints)."""
    if isinstance(v, bool):
        return "True" if v else "False"
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


# Shell-dangerous chars rejected in a per-param value or a per-token string
# (those get wrapped in double quotes when emitted, so a stray " breaks it).
_DANGEROUS_CHARS = "\"'`$\\"

# A pre-formed raw string is assumed already bracketed and double-quoted
# (``["type:...","type:..."]``), so its own quotes are structural and allowed.
# Only command-substitution / escape chars can corrupt the emitted shell arg.
_RAW_DANGEROUS_CHARS = "`$\\"

# Non-numeric potential params: bool flags and the free-form contact string.
# Everything else in _ALLOWED is a numeric weight/distance (binderlen a
# positive int), checked as such in GuidingPotential.__init__.
_BOOL_PARAMS = {"olig_intra_all", "olig_inter_all"}
_STRING_PARAMS = {"olig_custom_contact"}


def _reject_dangerous(s: str, context: str, chars: str = _DANGEROUS_CHARS) -> None:
    bad = set(chars) & set(s)
    if bad:
        raise ValueError(
            f"{context} contains {sorted(bad)}, which would corrupt the "
            "guiding-potential token or the generated shell argument."
        )


class GuidingPotential:
    """One RFdiffusion guiding-potential node (``type:...,key:val,...``).

    Construct via the typed classmethods (``substrate_contacts``, ``monomer_ROG``,
    …) rather than directly, so parameter names/types are checked. Each builder
    keeps only the parameters the user set (plus required ones), in a stable
    order, and serializes to the comma-joined ``type:...`` token RFdiffusion
    expects.

    The implemented potentials differ per RFdiffusion variant, so each tool
    exposes its OWN nested ``Tool.GuidingPotential`` subclass that whitelists
    only the types that tool supports (``_SUPPORTED``). Calling a factory for an
    unsupported type raises rather than emitting a token the model rejects at
    init. ``GuidingPotential`` itself (``_SUPPORTED = None``) accepts every known
    type and is not bound to any tool.
    """

    # type -> allowed parameter keys (order preserved for stable token output).
    # This is the full registry; per-tool subclasses narrow it via _SUPPORTED.
    _ALLOWED: Dict[str, List[str]] = {
        "ligand_ncontacts": ["weight", "r_0", "d_0"],
        "substrate_contacts": ["weight", "s", "r_0", "d_0", "rep_r_0", "rep_s", "rep_r_min"],
        "monomer_ROG": ["weight", "min_dist"],
        "monomer_contacts": ["weight", "r_0", "d_0"],
        "olig_contacts": ["weight_intra", "weight_inter", "r_0", "d_0",
                          "olig_intra_all", "olig_inter_all", "olig_custom_contact"],
        "binder_ROG": ["binderlen", "weight", "min_dist"],
        "binder_ncontacts": ["binderlen", "weight", "r_0", "d_0"],
        "interface_ncontacts": ["binderlen", "weight", "r_0", "d_0"],
    }

    # Some keys are NOT token params: RFdiffusion's PotentialManager reads them
    # from the top-level potentials config (e.g. potentials.olig_intra_all), not
    # from the `type:...,key:val` token. Putting them in the token makes the
    # per-node float parser choke (float('True')). Map type -> {param: hydra key}
    # so they are emitted as standalone `potentials.<key>=` overrides instead.
    _STANDALONE: Dict[str, Dict[str, str]] = {
        "olig_contacts": {
            "olig_intra_all": "potentials.olig_intra_all",
            "olig_inter_all": "potentials.olig_inter_all",
            "olig_custom_contact": "potentials.olig_custom_contact",
        },
    }

    # Per-tool whitelist of supported potential types. None = accept all (the
    # bare base class). Subclasses set this to the types their model implements.
    _SUPPORTED: Optional[set] = None

    def __init__(self, ptype: str, params: Dict[str, Any]):
        if ptype not in self._ALLOWED:
            raise ValueError(
                f"Unknown guiding-potential type {ptype!r}. Known types: "
                f"{', '.join(sorted(self._ALLOWED))}."
            )
        if self._SUPPORTED is not None and ptype not in self._SUPPORTED:
            raise ValueError(
                f"{type(self).__name__} does not support the {ptype!r} potential. "
                f"This variant implements: {', '.join(sorted(self._SUPPORTED))}."
            )
        allowed = self._ALLOWED[ptype]
        standalone_keys = self._STANDALONE.get(ptype, {})
        for k, v in params.items():
            if k not in allowed:
                raise ValueError(
                    f"{ptype}: unknown parameter {k!r}. Allowed: {', '.join(allowed)}."
                )
            # Type/range checks by parameter kind. Catch wrong types here rather
            # than after a slow job hits a Hydra/model error.
            if k in _BOOL_PARAMS:
                if not isinstance(v, bool):
                    raise ValueError(f"{ptype}.{k} must be a bool, got {type(v).__name__}")
            elif k in _STRING_PARAMS:
                if not isinstance(v, str):
                    raise ValueError(f"{ptype}.{k} must be a string, got {type(v).__name__}")
            elif k == "binderlen":
                if not isinstance(v, int) or isinstance(v, bool) or v <= 0:
                    raise ValueError(f"{ptype}.{k} must be a positive int, got {v!r}")
            else:  # numeric weights / distances
                if not isinstance(v, (int, float)) or isinstance(v, bool):
                    raise ValueError(f"{ptype}.{k} must be a number, got {type(v).__name__}")
                if v < 0:
                    raise ValueError(f"{ptype}.{k} must be non-negative, got {v}")
            # Reject anything that would break the hydra token / shell arg. The
            # comma/colon/bracket check is skipped for olig_custom_contact, whose
            # syntax legitimately uses them (e.g. "A&B,A!C") — but shell-dangerous
            # chars are still forbidden everywhere.
            tok = _fmt(v)
            structural = ",[] :"
            bad = set(_DANGEROUS_CHARS) & set(tok)
            if k != "olig_custom_contact":
                bad |= set(structural) & set(tok)
            if bad:
                raise ValueError(
                    f"{ptype}.{k}={v!r} contains {sorted(bad)}, which would corrupt the "
                    "guiding-potential token or the generated shell argument."
                )
        # Token params (everything except the standalone keys), in canonical order.
        self.ptype = ptype
        self.params = {k: params[k] for k in allowed
                       if k in params and k not in standalone_keys}
        # Standalone overrides: {hydra key: value} emitted as separate potentials.*=
        self.standalone = {standalone_keys[k]: params[k]
                           for k in allowed if k in params and k in standalone_keys}

    def to_token(self) -> str:
        """Serialize to ``type:<ptype>,<key>:<val>,...`` (no surrounding quotes)."""
        parts = [f"type:{self.ptype}"]
        parts += [f"{k}:{_fmt(v)}" for k, v in self.params.items()]
        return ",".join(parts)

    def __repr__(self) -> str:
        return f"GuidingPotential({self.to_token()!r})"

    # --- Typed factory methods (one per RFdiffusion potential type) ---------

    @classmethod
    def ligand_ncontacts(cls, weight: float = 1, r_0: float = 8,
                         d_0: float = 4) -> "GuidingPotential":
        """Pull the design into contact with the ligand (maximize protein-ligand
        contacts). This is the ONLY potential RFdiffusion-AllAtom implements; use it
        to steer an all-atom design to wrap the bound ligand. weight 1-10; r_0 ~8 Å
        (contact switching distance); d_0 ~4 Å (always-in-contact distance)."""
        return cls("ligand_ncontacts", _drop_none(dict(weight=weight, r_0=r_0, d_0=d_0)))

    @classmethod
    def substrate_contacts(cls, weight: float = 1, s: float = None, r_0: float = None,
                           d_0: float = None, rep_r_0: float = None, rep_s: float = None,
                           rep_r_min: float = None) -> "GuidingPotential":
        """Pull the design around a ligand/substrate (attractive + repulsive terms).

        Requires ``potentials.substrate=<ligand code>`` (the tool fills this from
        the ligand automatically). weight 1-10; s 0.1-2; r_0 6-12 Å; d_0 2-6 Å;
        rep_r_0 2-6 Å; rep_s 1-10; rep_r_min 1-3 Å.
        """
        return cls("substrate_contacts", _drop_none(dict(
            weight=weight, s=s, r_0=r_0, d_0=d_0,
            rep_r_0=rep_r_0, rep_s=rep_s, rep_r_min=rep_r_min,
        )))

    @classmethod
    def monomer_ROG(cls, weight: float = 1, min_dist: float = None) -> "GuidingPotential":
        """Encourage a compact monomer (minimize radius of gyration). weight 1-10; min_dist 10-20 Å."""
        return cls("monomer_ROG", _drop_none(dict(weight=weight, min_dist=min_dist)))

    @classmethod
    def monomer_contacts(cls, weight: float = 1, r_0: float = None,
                         d_0: float = None) -> "GuidingPotential":
        """Encourage intra-monomer contacts. weight 1-10; r_0 6-12 Å; d_0 2-6 Å."""
        return cls("monomer_contacts", _drop_none(dict(weight=weight, r_0=r_0, d_0=d_0)))

    @classmethod
    def binder_ROG(cls, binderlen: int, weight: float = 1,
                   min_dist: float = None) -> "GuidingPotential":
        """Compact binder (radius of gyration over the first ``binderlen`` residues)."""
        return cls("binder_ROG", _drop_none(dict(binderlen=binderlen, weight=weight, min_dist=min_dist)))

    @classmethod
    def binder_ncontacts(cls, binderlen: int, weight: float = 1, r_0: float = None,
                         d_0: float = None) -> "GuidingPotential":
        """Maximize intra-binder contacts. weight 1-10; r_0 6-12 Å; d_0 2-6 Å."""
        return cls("binder_ncontacts", _drop_none(dict(binderlen=binderlen, weight=weight, r_0=r_0, d_0=d_0)))

    @classmethod
    def interface_ncontacts(cls, binderlen: int, weight: float = 1, r_0: float = None,
                            d_0: float = None) -> "GuidingPotential":
        """Maximize binder-target interface contacts. weight 1-10; r_0 6-12 Å; d_0 2-6 Å."""
        return cls("interface_ncontacts", _drop_none(dict(binderlen=binderlen, weight=weight, r_0=r_0, d_0=d_0)))

    @classmethod
    def olig_contacts(cls, weight_intra: float = 1, weight_inter: float = 1,
                      r_0: float = None, d_0: float = None,
                      olig_intra_all: bool = None, olig_inter_all: bool = None,
                      olig_custom_contact: str = None) -> "GuidingPotential":
        """Inter/intra-chain contacts for symmetric oligomers. weight_* 0.1-2."""
        return cls("olig_contacts", _drop_none(dict(
            weight_intra=weight_intra, weight_inter=weight_inter, r_0=r_0, d_0=d_0,
            olig_intra_all=olig_intra_all, olig_inter_all=olig_inter_all,
            olig_custom_contact=olig_custom_contact,
        )))


def _drop_none(d: Dict[str, Any]) -> Dict[str, Any]:
    return {k: v for k, v in d.items() if v is not None}


def render_guiding_potentials(
    potentials: Union[str, GuidingPotential, List[Union[str, GuidingPotential]], None]
) -> Optional[str]:
    """Render a builder / list of builders / raw string into the hydra value.

    - None -> None.
    - A raw string is returned unchanged (back-compat; assumed already bracketed).
    - A GuidingPotential or list thereof -> ``["type:...","type:..."]`` with each
      node double-quoted (RFdiffusion's expected list form). Strings mixed into a
      list are taken as already-formed tokens.
    """
    if potentials is None:
        return None
    if isinstance(potentials, str):
        _reject_dangerous(potentials, "guiding_potentials raw string", _RAW_DANGEROUS_CHARS)
        return potentials
    nodes = potentials if isinstance(potentials, (list, tuple)) else [potentials]
    tokens = []
    for n in nodes:
        if isinstance(n, GuidingPotential):
            tokens.append(n.to_token())
        elif isinstance(n, str):
            token = n.strip().strip('"')
            _reject_dangerous(token, "guiding_potentials raw token")
            tokens.append(token)
        else:
            raise ValueError(
                f"guiding_potentials entries must be GuidingPotential or str, got {type(n).__name__}"
            )
    inner = ",".join(f'"{t}"' for t in tokens)
    return f"[{inner}]"


def guiding_potential_overrides(
    potentials: Union[str, GuidingPotential, List[Union[str, GuidingPotential]], None]
) -> List[str]:
    """Standalone ``potentials.<key>=<val>`` overrides for the given potentials.

    Some builder params (e.g. olig_contacts' olig_intra_all/olig_inter_all/
    olig_custom_contact) are read by RFdiffusion from the top-level potentials
    config, not from the ``type:...`` token. They are collected here and emitted
    as separate hydra overrides alongside ``potentials.guiding_potentials``. A
    raw-string potentials value carries no structured info, so it yields none.
    """
    if potentials is None or isinstance(potentials, str):
        return []
    nodes = potentials if isinstance(potentials, (list, tuple)) else [potentials]
    overrides = []
    for n in nodes:
        if isinstance(n, GuidingPotential):
            overrides += [f"{k}={_fmt(v)}" for k, v in n.standalone.items()]
    return overrides


class RFdiffusionGuidingPotential(GuidingPotential):
    """Potentials implemented by base RFdiffusion (scripts/run_inference.py).

    Bound to ``RFdiffusion.GuidingPotential``. Covers monomer/oligomer/binder
    contact and radius-of-gyration potentials; the symmetric-oligomer workflow
    uses ``olig_contacts``.
    """

    _SUPPORTED = {
        "monomer_ROG", "monomer_contacts", "olig_contacts",
        "binder_ROG", "binder_ncontacts", "interface_ncontacts",
        "substrate_contacts",
    }


class RFdiffusionAllAtomGuidingPotential(GuidingPotential):
    """Potentials implemented by RFdiffusion-AllAtom (only ``ligand_ncontacts``).

    Bound to ``RFdiffusionAllAtom.GuidingPotential``. Its
    ``potentials/potentials.py`` whitelists a single type; anything else trips an
    AssertionError at model init, so the others are rejected up front here.
    """

    _SUPPORTED = {"ligand_ncontacts"}
