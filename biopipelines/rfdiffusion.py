# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
RFdiffusion configuration for protein backbone generation.

RFdiffusion is a diffusion-based generative model for designing protein backbones.
It supports unconditional generation, motif scaffolding, binder design, and partial
diffusion. See https://github.com/RosettaCommons/RFdiffusion for full documentation.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import Resolve, TableReference
    from .combinatorics import generate_multiplied_ids_pattern
    from .guiding_potentials import (RFdiffusionGuidingPotential,
                                      render_guiding_potentials, guiding_potential_overrides)
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import Resolve, TableReference
    from combinatorics import generate_multiplied_ids_pattern
    from guiding_potentials import (RFdiffusionGuidingPotential,
                                    render_guiding_potentials, guiding_potential_overrides)


# Point groups accepted by inference.symmetry. Cyclic/dihedral take a degree
# (c4, d2, …); the polyhedral groups are fixed names.
_SYMMETRY_NAMED = {"tetrahedral", "octahedral", "icosahedral"}

# Decay schedules accepted by potentials.guide_decay.
_GUIDE_DECAYS = {"constant", "linear", "quadratic", "cubic"}


def _normalize_selection_arg(name, value):
    """Normalize a per-PDB selection argument to a (kind, token) pair.

    RFdiffusion's ``contigs`` / ``inpaint`` / ``inpaint_str`` each accept either
    a plain string — broadcast to every input PDB — or a per-PDB column
    reference (a ``TableReference`` from ``tool.tables.X.col``, or a
    ``(TableInfo, "col")`` tuple). A reference is resolved per-id at runtime via
    ``Resolve.table_column``; a literal is emitted directly.

    Returns ``("literal", str)`` or ``("table", "TABLE_REFERENCE:path:col")``.
    A literal is the only kind validated by ``_validate_freeform_string`` (a
    table reference reaches bash only through a resolved subshell, not as raw
    interpolation).
    """
    if isinstance(value, TableReference):
        return ("table", str(value))
    # (TableInfo, "column") tuple — mirror the LigandMPNN/ThermoMPNN convention.
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], str):
        return ("table", str(TableReference(value[0].info.path, value[1])))
    if isinstance(value, str):
        return ("literal", value)
    raise ValueError(
        f"{name} must be a string or a (table, column) reference, got {type(value)}"
    )


def _selection_cli_arg(kind_token) -> str:
    """CLI arg for the contigs resolver: the token, or '-' when empty.

    ``kind_token`` is the ``(kind, token)`` pair from
    ``_normalize_selection_arg``. Shared by all three RFdiffusion tools.
    """
    kind, token = kind_token
    if kind == "literal" and not token:
        return "-"
    return token


def _normalize_symmetry(value):
    """Normalize a symmetry spec to the token inference.symmetry expects.

    Accepts a point-group name (``tetrahedral`` / ``octahedral`` / ``icosahedral``)
    or a cyclic/dihedral group as ``c<N>`` / ``d<N>`` (also spelled ``cyclic_N`` /
    ``C4`` etc.). Returns the lowercase token RFdiffusion wants (``c4``, ``d2``,
    ``tetrahedral``). Raises ValueError on anything else.
    """
    s = str(value).strip().lower()
    if s in _SYMMETRY_NAMED:
        return s
    s = s.replace("cyclic_", "c").replace("dihedral_", "d")
    if (s.startswith("c") or s.startswith("d")) and s[1:].isdigit() and int(s[1:]) >= 1:
        return s
    raise ValueError(
        f"symmetry must be one of {sorted(_SYMMETRY_NAMED)}, or a cyclic/dihedral "
        f"group like 'c4' / 'd2', got {value!r}"
    )


class RFdiffusion(BaseConfig):
    """
    Configuration for RFdiffusion protein backbone generation.

    RFdiffusion generates novel protein backbones using a denoising diffusion
    probabilistic model conditioned on structural motifs, partial structures,
    or protein-protein interaction targets.

    Main workflows:
        - Unconditional generation: produce backbones of a given length range
          without any input structure (contigs only, no pdb).
        - Motif scaffolding: hold a fixed structural motif in place and design
          the surrounding scaffold (pdb + contigs specifying the motif).
        - Binder design: design a new protein that binds to a target (pdb of
          the target + contigs describing the binder length), optionally steered
          to specific interface residues with ``hotspot_res``.
        - Partial diffusion: apply limited noise to an existing structure and
          re-diffuse, producing near-neighbour variants (partial_steps).
        - Symmetric oligomers: build a point-group-symmetric assembly from a
          single contigs length (``symmetry`` + optional ``olig_contacts``
          guiding potential).
        - Fold conditioning: condition diffusion on a target topology / scaffold
          set (``scaffold_dir`` or ``target_ss`` + ``target_adj``).
        - Cyclic peptides: close the backbone into a macrocycle (``cyclic``).

    Model checkpoints (see WEIGHTS / DEFAULT_WEIGHTS):
        By default only Base and Complex_base are downloaded (~2 GB). Add others
        via RFdiffusion.install(weights=[...]) when needed.

    Installation:
        with Pipeline(...):
            RFdiffusion.install()                         # default weights
            RFdiffusion.install(weights=["Base",
                                         "InpaintSeq"])   # custom selection
            rfd = RFdiffusion(contigs="50-100", ...)

    Reference:
        Watson et al. (2023) De novo design of protein structure and function
        with RFdiffusion. Nature 620, 1089-1100.
        https://github.com/RosettaCommons/RFdiffusion
    """

    TOOL_NAME = "RFdiffusion"
    TOOL_VERSION = "1.1"

    # Typed builder for guiding potentials, e.g.
    #   RFdiffusion.GuidingPotential.olig_contacts(weight_intra=1, weight_inter=0.1)
    # Restricted to the types base RFdiffusion implements (see the subclass).
    GuidingPotential = RFdiffusionGuidingPotential

    # Mapping of weight name -> (url_hash, filename)
    WEIGHTS = {
        "Base":              ("6f5902ac237024bdd0c176cb93063dc4", "Base_ckpt.pt"),
        "Complex_base":      ("e29311f6f1bf1af907f9ef9f44b8328b", "Complex_base_ckpt.pt"),
        "Complex_Fold_base": ("60f09a193fb5e5ccdc4980417708dbab", "Complex_Fold_base_ckpt.pt"),
        "InpaintSeq":        ("74f51cfb8b440f50d70878e05361d8f0", "InpaintSeq_ckpt.pt"),
        "InpaintSeq_Fold":   ("76d00716416567174cdb7ca96e208296", "InpaintSeq_Fold_ckpt.pt"),
        "ActiveSite":        ("5532d2e1f3a4738decd58b19d633b3c3", "ActiveSite_ckpt.pt"),
        "Base_epoch8":       ("12fc204edeae5b57713c5ad7dcb97d39", "Base_epoch8_ckpt.pt"),
        "RF_structure_prediction": ("1befcb9b28e2f778f53d47f18b7597fa", "RF_structure_prediction_weights.pt"),
    }
    # Weights downloaded by default — covers the two main workflows:
    #   Base          → unconditional generation, motif scaffolding
    #   Complex_base  → binder design / protein-protein interaction
    # Add more via RFdiffusion.install(weights=[...]):
    #   Complex_Fold_base  → binder design with fold/topology conditioning
    #   InpaintSeq         → auto-selected when contigmap.inpaint_seq is used
    #   InpaintSeq_Fold    → inpaint + fold conditioning together
    #   ActiveSite         → small active-site motifs (requires ckpt_override_path)
    #   Base_epoch8        → alternative base checkpoint
    #   RF_structure_prediction → legacy structure prediction (unrelated to diffusion)
    DEFAULT_WEIGHTS = ["Base", "Complex_base"]

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False,
                        weights=None, **kwargs):
        """
        Generate the bash installation script for RFdiffusion.

        Clones the RFdiffusion repository, downloads the requested model
        checkpoints, and creates (or reuses) the SE3nv conda/mamba environment.

        Args:
            folders: Resolved pipeline folder paths (must contain "RFdiffusion"
                     and "biopipelines" keys).
            env_manager: "mamba", "conda", or "micromamba".
            force_reinstall: If True, skip the already-installed early-exit check.
            weights: List of checkpoint names to download. Defaults to
                     DEFAULT_WEIGHTS (["Base", "Complex_base"]).
                     Valid names are the keys of WEIGHTS.
        """
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("RFdiffusion", "")
        parent_dir = os.path.dirname(repo_dir)

        if weights is None:
            weights = cls.DEFAULT_WEIGHTS
        invalid = [w for w in weights if w not in cls.WEIGHTS]
        if invalid:
            raise ValueError(
                f"Unknown weight(s): {invalid}. "
                f"Valid options are: {list(cls.WEIGHTS.keys())}"
            )
        wget_lines = "\n".join(
            f"wget -nc http://files.ipd.uw.edu/pub/RFdiffusion/{cls.WEIGHTS[w][0]}/{cls.WEIGHTS[w][1]}"
            f" || echo \"WARNING: failed to download {cls.WEIGHTS[w][1]}, skipping\""
            for w in weights
        )

        # Skip only when the env exists AND `import rfdiffusion` actually
        # works inside it. Just the env name being present is not enough:
        # ProteinMPNN.install() may have created a leaner SE3nv (PyTorch
        # only) before us, in which case we must layer DGL, SE3Transformer,
        # and the rfdiffusion package on top instead of skipping.
        # The import probe is the load-bearing check (env + package); add the
        # weights-file test on top so we don't skip when checkpoints are missing.
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}/models" ] && [ -f "{repo_dir}/models/Base_ckpt.pt" ] \\
   && {env_manager} run -n SE3nv python -c "import rfdiffusion" >/dev/null 2>&1; then
    echo "RFdiffusion already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("SE3nv", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("SE3nv", env_manager, biopipelines)
        return f"""echo "=== Installing RFdiffusion ==="
{skip}mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/RosettaCommons/RFdiffusion.git
fi
cd {repo_dir}

# Download model weights
mkdir -p models && cd models
{wget_lines}
cd ..

# Create SE3nv environment from BioPipelines specification
{remove_block}
{env_block}
if [ $? -ne 0 ]; then
    echo "WARNING: BioPipelines SE3nv env creation failed. Trying official RFdiffusion environment..."
    {env_manager} env create -f env/SE3nv.yml
    if [ $? -ne 0 ]; then
        echo "ERROR: SE3nv environment creation failed with both methods."
        echo "This is likely a CUDA version mismatch for your system."
        exit 1
    fi
fi

# Install DGL from pre-built wheel (special --find-links syntax)
{env_manager} run -n SE3nv pip install dgl -f https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html --no-deps

# Install SE3Transformer and RFdiffusion (editable)
cd env/SE3Transformer
{env_manager} run -n SE3nv pip install --no-deps .
cd ../..
{env_manager} run -n SE3nv pip install -e .

# Verify installation
if {env_manager} run -n SE3nv python -c "import rfdiffusion" >/dev/null 2>&1 || [ -f "{repo_dir}/models/Base_ckpt.pt" ]; then
    touch "$INSTALL_SUCCESS"
    echo "=== RFdiffusion installation complete ==="
else
    echo "ERROR: RFdiffusion verification failed"
    exit 1
fi
"""

    # Lazy path descriptors
    #   main_table        — TableInfo CSV describing each design's provenance
    #                       (fixed/designed regions, pLDDT). Lives under tables/.
    #   table_py_file     — pipe script that builds main_table from .trb files.
    #   inference_py_file — RFdiffusion's own CLI entry point.
    #   pdb_ds_json       — serialized input DataStream (config-time artifact).
    #   update_map_py     — pipe script that rewrites structures_map.csv to
    #                       match actual files on disk after the run.
    main_table = Path(lambda self: self.table_path("structures"))
    table_py_file = Path(lambda self: self.pipe_script_path("pipe_rfdiffusion_table.py"))
    inference_py_file = Path(lambda self: os.path.join(self.folders["RFdiffusion"], "scripts", "run_inference.py"))
    pdb_ds_json = Path(lambda self: self.configuration_path("input_structures.json"))
    update_map_py = Path(lambda self: self.pipe_script_path("pipe_update_structures_map.py"))
    #   contigs_json       — per-PDB resolved {contigs, inpaint, inpaint_str},
    #                        written at runtime by pipe_rfdiffusion_contigs.py.
    #   contigs_py         — that resolver (one pass over all ids → JSON).
    #   contigs_reader_py  — per-id reader the bash loop calls for each PDB.
    contigs_json = Path(lambda self: self.configuration_path("contig_options.json"))
    contig_args_json = Path(lambda self: self.configuration_path("contig_args.json"))
    contigs_py = Path(lambda self: self.pipe_script_path("pipe_rfdiffusion_contigs.py"))
    contigs_reader_py = Path(lambda self: self.pipe_script_path("resolve_rfdiffusion_contigs.py"))

    def __init__(self,
                 contigs: Union[str, TableReference, "tuple"] = "",
                 pdb: Optional[Union[DataStream, StandardizedOutput]] = None,
                 inpaint: Union[str, TableReference, "tuple"] = "",
                 inpaint_str: Union[str, TableReference, "tuple"] = "",
                 num_designs: int = 1,
                 active_site: bool = False,
                 steps: int = 50,
                 partial_steps: int = 0,
                 reproducible: bool = False,
                 design_startnum: int = 1,
                 hotspot_res: Optional[List[str]] = None,
                 symmetry: Optional[str] = None,
                 guiding_potentials: Union[str, "GuidingPotential", List, None] = None,
                 guide_scale: Optional[float] = None,
                 guide_decay: Optional[str] = None,
                 cyclic: bool = False,
                 cyc_chains: Optional[str] = None,
                 provide_seq: Optional[str] = None,
                 noise_scale_ca: Optional[float] = None,
                 noise_scale_frame: Optional[float] = None,
                 inpaint_str_helix: Optional[str] = None,
                 inpaint_str_strand: Optional[str] = None,
                 scaffold_dir: Optional[str] = None,
                 target_ss: Optional[str] = None,
                 target_adj: Optional[str] = None,
                 target_path: Optional[str] = None,
                 mask_loops: Optional[bool] = None,
                 sampled_insertion: Optional[int] = None,
                 sampled_N: Optional[int] = None,
                 sampled_C: Optional[int] = None,
                 **kwargs):
        """
        Initialize RFdiffusion configuration.

        Args:
            contigs: Contig map describing which residues are fixed and which
                     are to be generated. Residues from the input PDB are
                     specified as chain+range (e.g. "A1-100"); new backbone
                     segments are specified as length ranges (e.g. "50-100").
                     Multiple segments are separated by "/" (e.g. "A1-50/30-50/A60-100").
                     For unconditional generation (no pdb) use a length range alone
                     (e.g. "100-200"). Either a plain string (broadcast to every
                     input PDB) or a per-PDB column reference
                     (``tool.tables.X.col`` / ``(TableInfo, "col")``) resolved by
                     id match at runtime — use the latter to give each input PDB
                     its own contigs. A column reference requires ``pdb``.
            pdb: Optional input structure(s). Required for motif scaffolding,
                 binder design, and partial diffusion. Accepts a DataStream or a
                 StandardizedOutput. When it carries multiple structures,
                 RFdiffusion runs once per input PDB (``num_designs`` designs
                 each); output ids are ``<pdb_id>_<n>``.
            inpaint: Residues whose sequence should be masked during diffusion
                     (same chain+range format as contigs). When set, the
                     InpaintSeq checkpoint is used automatically — ensure it
                     was downloaded at install time. Accepts the same
                     string-or-column-reference forms as ``contigs``.
            inpaint_str: Residues whose secondary structure should be masked
                     during diffusion (same chain+range format as contigs).
                     Wired to upstream ``contigmap.inpaint_str``. Empty string
                     disables (default). The overall length stays governed by
                     the ``contigs`` argument (use a range like ``"50-100"``
                     to control it). Accepts the same string-or-column-reference
                     forms as ``contigs``.
            num_designs: Number of independent backbone designs to generate.
            active_site: If True, use the ActiveSite checkpoint instead of
                         Base. Intended for scaffolding very small functional
                         motifs (< ~10 residues). Requires the ActiveSite
                         weight to have been downloaded at install time.
            steps: Number of denoising diffusion steps (default 50). Fewer
                   steps are faster but may reduce quality.
            partial_steps: Number of partial diffusion steps. When > 0, the
                           input structure is noised for this many steps and
                           then re-diffused, producing near-neighbour variants.
                           Must be less than steps.
            reproducible: If True, use deterministic (fixed-seed) sampling so
                          the same contigs always produce the same outputs.
            design_startnum: Integer appended to the pipeline name to number
                             output files (e.g. design_startnum=1 → name_1.pdb).
                             Useful when continuing a previous run.
            hotspot_res: Interface residues on the target the binder should
                         engage, e.g. ``["A30", "A33", "A34"]``
                         (``ppi.hotspot_res``). Binder design only — requires
                         ``pdb``.
            symmetry: Point group for symmetric-oligomer generation. A named
                      group (``"tetrahedral"``, ``"octahedral"``,
                      ``"icosahedral"``) or a cyclic/dihedral group
                      (``"c4"``, ``"d2"``, …). Switches RFdiffusion to its
                      ``symmetry`` config (``--config-name symmetry
                      inference.symmetry=...``). Symmetric runs are
                      unconditional, so ``pdb`` must not be set; ``contigs`` is
                      the total length of the assembly (e.g. ``"360"`` for a
                      tetrahedral 360-mer). Pair with an ``olig_contacts``
                      guiding potential to enforce inter/intra-chain contacts.
            guiding_potentials: Guiding potential(s) biasing the diffusion
                      trajectory (``potentials.guiding_potentials``). Use the
                      typed builder, e.g.
                      ``RFdiffusion.GuidingPotential.olig_contacts(weight_intra=1,
                      weight_inter=0.1, olig_intra_all=True, olig_inter_all=True)``;
                      a raw hydra string or a list of builders is also accepted.
            guide_scale: Global multiplier on the guiding potential
                      (``potentials.guide_scale``).
            guide_decay: How the potential decays over the trajectory
                      (``potentials.guide_decay``): ``"constant"``, ``"linear"``,
                      ``"quadratic"``, or ``"cubic"``.
            cyclic: If True, close the designed backbone into a macrocycle
                    (``inference.cyclic``). Use with ``cyc_chains``.
            cyc_chains: Chain letter(s) to cyclize (``inference.cyc_chains``,
                        e.g. ``"a"``). Only meaningful when ``cyclic=True``.
            provide_seq: Residue range(s) whose sequence is kept fixed during
                         partial diffusion (``contigmap.provide_seq``, e.g.
                         ``"100-119"``). Requires ``partial_steps > 0``.
            noise_scale_ca: Translational noise scale (``denoiser.noise_scale_ca``).
                            Values < 1 reduce diversity but improve quality
                            (binder design often uses ~0).
            noise_scale_frame: Rotational noise scale
                               (``denoiser.noise_scale_frame``). Same trade-off
                               as ``noise_scale_ca``.
            inpaint_str_helix: Residues whose secondary structure is masked and
                               forced to helix (``contigmap.inpaint_str_helix``).
            inpaint_str_strand: Residues whose secondary structure is masked and
                                forced to strand (``contigmap.inpaint_str_strand``).
            scaffold_dir: Directory of scaffold ``_ss.pt`` / ``_adj.pt`` files for
                          fold conditioning (``scaffoldguided.scaffold_dir``).
                          Setting it turns on ``scaffoldguided.scaffoldguided``.
                          Fold conditioning runs from the scaffold set rather than
                          ``contigs``.
            target_ss: Precomputed target secondary-structure tensor file
                       (``scaffoldguided.target_ss``, a ``.pt`` produced by
                       upstream ``helper_scripts/make_secstruc_adj.py`` — a small
                       per-target input, not a model weight). Implies
                       ``scaffoldguided.target_pdb=True``.
            target_adj: Precomputed target block-adjacency tensor file
                        (``scaffoldguided.target_adj``, the companion ``.pt``).
            target_path: Target structure for fold-conditioned binder design
                         (``scaffoldguided.target_path``). Used with
                         ``target_ss`` / ``target_adj``.
            mask_loops: If False, keep input loops fixed; if True, allow loops to
                        be remodelled (``scaffoldguided.mask_loops``).
            sampled_insertion: Max residues to insert when sampling scaffold
                               length variants (``scaffoldguided.sampled_insertion``).
            sampled_N: Max residues to sample at the N terminus
                       (``scaffoldguided.sampled_N``).
            sampled_C: Max residues to sample at the C terminus
                       (``scaffoldguided.sampled_C``).

        Output:
            Streams: structures (.pdb)
            Tables:
                structures: id | pdb | fixed | designed | source_fixed | plddt_mean | status
        """
        # Resolve optional pdb input — store stream for runtime resolution.
        # All input PDBs are iterated at execution time (see generate_script);
        # never index ids[0] at config time, which breaks under lazy IDs.
        self.pdb_stream: Optional[DataStream] = None
        if pdb is not None:
            if isinstance(pdb, StandardizedOutput):
                self.pdb_stream = pdb.streams.structures
            elif isinstance(pdb, DataStream):
                self.pdb_stream = pdb
            else:
                raise ValueError(f"pdb must be DataStream or StandardizedOutput, got {type(pdb)}")

        # Normalize the per-PDB selection args to (kind, token) pairs. A literal
        # is broadcast to every PDB; a table reference is resolved per-id at
        # runtime. Store both the original (for display) and the normalized form.
        self.contigs = contigs
        self.inpaint = inpaint
        self.inpaint_str = inpaint_str
        self._contigs_arg = _normalize_selection_arg("contigs", contigs)
        self._inpaint_arg = _normalize_selection_arg("inpaint", inpaint)
        self._inpaint_str_arg = _normalize_selection_arg("inpaint_str", inpaint_str)
        self.num_designs = num_designs
        self.active_site = active_site
        self.steps = steps
        self.partial_steps = partial_steps
        self.reproducible = reproducible
        self.design_startnum = design_startnum

        self.hotspot_res = hotspot_res or []
        self.symmetry = _normalize_symmetry(symmetry) if symmetry is not None else None
        # Accept a GuidingPotential builder (or list), or a raw hydra string;
        # render to the single bracketed string the shell emits. Some builder
        # params (olig_*) are top-level potentials.* overrides, not token params.
        self.guiding_potentials = render_guiding_potentials(guiding_potentials)
        self._guiding_potential_overrides = guiding_potential_overrides(guiding_potentials)
        self.guide_scale = guide_scale
        self.guide_decay = guide_decay
        self.cyclic = cyclic
        self.cyc_chains = cyc_chains
        self.provide_seq = provide_seq
        self.noise_scale_ca = noise_scale_ca
        self.noise_scale_frame = noise_scale_frame
        self.inpaint_str_helix = inpaint_str_helix
        self.inpaint_str_strand = inpaint_str_strand
        # Fold/scaffold conditioning (scaffoldguided.*). scaffold_dir or a
        # target_ss/target_adj pair turns the mode on.
        self.scaffold_dir = scaffold_dir
        self.target_ss = target_ss
        self.target_adj = target_adj
        self.target_path = target_path
        self.mask_loops = mask_loops
        self.sampled_insertion = sampled_insertion
        self.sampled_N = sampled_N
        self.sampled_C = sampled_C
        self._scaffold_guided = bool(scaffold_dir or target_ss or target_adj)

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate RFdiffusion-specific parameters."""
        contigs_kind, contigs_token = self._contigs_arg
        # Fold conditioning runs from a scaffold set (scaffoldguided.*), not from
        # contigs; in that mode contigs is optional. Every other mode needs it.
        if not self._scaffold_guided and contigs_kind == "literal" and not contigs_token:
            raise ValueError("contigs parameter is required")

        if self.num_designs <= 0:
            raise ValueError("num_designs must be positive")

        if self.steps <= 0:
            raise ValueError("steps must be positive")

        if self.partial_steps < 0:
            raise ValueError("partial_steps cannot be negative")

        # Symmetric generation uses RFdiffusion's `symmetry` config and builds
        # the assembly from a contigs length — it is unconditional, so a pdb
        # would be ignored. Reject the combination rather than silently drop it.
        if self.symmetry is not None and self.pdb_stream is not None:
            raise ValueError(
                "symmetry generates an unconditional symmetric oligomer and "
                "cannot be combined with a pdb input"
            )

        # Hotspots only mean something against a target structure.
        if self.hotspot_res and self.pdb_stream is None:
            raise ValueError("hotspot_res requires a pdb input (the binder target)")

        # provide_seq fixes sequence during partial diffusion.
        if self.provide_seq and self.partial_steps <= 0:
            raise ValueError("provide_seq requires partial_steps > 0")

        if self.cyc_chains is not None and not self.cyclic:
            raise ValueError("cyc_chains requires cyclic=True")

        # Fold conditioning is incompatible with symmetric generation (distinct
        # hydra configs).
        if self._scaffold_guided and self.symmetry is not None:
            raise ValueError("scaffold conditioning cannot be combined with symmetry")

        # Boundary checks on the typed params — catch invalid values here rather
        # than after a slow job lands on a Hydra/model error.
        if self.guide_decay is not None and self.guide_decay not in _GUIDE_DECAYS:
            raise ValueError(
                f"guide_decay must be one of {sorted(_GUIDE_DECAYS)}, got {self.guide_decay!r}"
            )
        for name in ("guide_scale", "noise_scale_ca", "noise_scale_frame"):
            v = getattr(self, name)
            if v is not None:
                if not isinstance(v, (int, float)) or isinstance(v, bool):
                    raise ValueError(f"{name} must be a number, got {type(v).__name__}")
                if v < 0:
                    raise ValueError(f"{name} must be non-negative, got {v}")
        for name in ("sampled_insertion", "sampled_N", "sampled_C"):
            v = getattr(self, name)
            if v is not None:
                if not isinstance(v, int) or isinstance(v, bool):
                    raise ValueError(f"{name} must be an int, got {type(v).__name__}")
                if v < 0:
                    raise ValueError(f"{name} must be non-negative, got {v}")

        # A column reference is keyed by input-PDB id, so it needs PDBs to
        # match against.
        for name, (kind, _) in (("contigs", self._contigs_arg),
                                 ("inpaint", self._inpaint_arg),
                                 ("inpaint_str", self._inpaint_str_arg)):
            if kind == "table" and self.pdb_stream is None:
                raise ValueError(
                    f"{name} given as a table column reference requires a pdb input"
                )

        # Only literal selections are interpolated raw into bash; table
        # references reach the shell through a resolved subshell.
        for name, (kind, token) in (("contigs", self._contigs_arg),
                                     ("inpaint", self._inpaint_arg),
                                     ("inpaint_str", self._inpaint_str_arg)):
            if kind == "literal":
                _validate_freeform_string(name, token)

        # Free-form strings reaching bash as raw interpolation.
        for i, res in enumerate(self.hotspot_res):
            _validate_freeform_string(f"hotspot_res[{i}]", res)
        for name in ("cyc_chains", "provide_seq", "guide_decay",
                     "inpaint_str_helix", "inpaint_str_strand",
                     "scaffold_dir", "target_ss", "target_adj", "target_path"):
            _validate_freeform_string(name, getattr(self, name))

        # guiding_potentials is a structured hydra value (commas, brackets,
        # colons, double quotes) emitted single-quoted, so the generic freeform
        # check would reject it. Forbid only what breaks single quoting.
        if self.guiding_potentials is not None:
            if not isinstance(self.guiding_potentials, str):
                raise ValueError("guiding_potentials must be a string")
            bad = set("'`$\\") & set(self.guiding_potentials)
            if bad:
                raise ValueError(
                    f"guiding_potentials contains {sorted(bad)}, which would break the "
                    "single-quoted shell argument."
                )

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get RFdiffusion configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"CONTIGS: {self.contigs}",
            f"NUM DESIGNS: {self.num_designs}",
            f"ACTIVE SITE: {self.active_site}",
            f"STEPS: {self.steps}"
        ])

        if self.pdb_stream:
            config_lines.append(f"PDB: {', '.join(self.pdb_stream.ids)}")
        if self._inpaint_arg[1]:
            config_lines.append(f"INPAINT: {self.inpaint}")
        if self._inpaint_str_arg[1]:
            config_lines.append(f"INPAINT_STR: {self.inpaint_str}")
        if self.partial_steps > 0:
            config_lines.append(f"PARTIAL STEPS: {self.partial_steps}")
        if self.reproducible:
            config_lines.append(f"REPRODUCIBLE: {self.reproducible}")
        if self.symmetry is not None:
            config_lines.append(f"SYMMETRY: {self.symmetry}")
        if self.hotspot_res:
            config_lines.append(f"HOTSPOTS: {','.join(self.hotspot_res)}")
        if self.guiding_potentials is not None:
            config_lines.append(f"GUIDING POTENTIALS: {self.guiding_potentials}")
        if self.cyclic:
            config_lines.append(f"CYCLIC: {self.cyclic}")
        if self._scaffold_guided:
            config_lines.append("FOLD CONDITIONING: scaffoldguided")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate RFdiffusion execution script."""
        # Serialize input DataStream to JSON for runtime file resolution.
        # The configuration/ folder was already created by the pipeline
        # after get_output_files() returned.
        if self.pdb_stream:
            self.pdb_stream.save_json(self.pdb_ds_json)

        script_content = "#!/bin/bash\n"
        script_content += "# RFdiffusion execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        # e3nn 0.3.3 uses torch.load() without weights_only=False,
        # which fails on PyTorch 2.6+ where the default flipped to True
        script_content += "export TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1\n" 
        script_content += self._generate_script_run_rfdiffusion()
        script_content += self._generate_script_create_table()
        script_content += self._generate_script_update_structures_map()
        script_content += self.generate_completion_check_footer()
        return script_content

    def _common_rfd_options(self) -> List[str]:
        """Inference args shared across runs (no contigs/pdb/prefix).

        Returned as a list of hydra tokens; each is emitted as one double-quoted
        bash-array element (see _generate_script_run_rfdiffusion), so values that
        carry brackets, commas, or colons (hotspot_res, guiding_potentials) need
        no extra quoting and are never re-parsed by `eval`.
        """
        opts = [
            f"inference.num_designs={self.num_designs}",
            f"inference.deterministic={self.reproducible}",
            f"inference.design_startnum={self.design_startnum}",
        ]
        if self.steps != 50:
            opts.append(f"diffuser.T={self.steps}")
        if self.partial_steps > 0:
            opts.append(f"diffuser.partial_T={self.partial_steps}")
        if self.active_site:
            opts.append("inference.ckpt_override_path=models/ActiveSite_ckpt.pt")

        if self.symmetry is not None:
            opts.append(f"inference.symmetry={self.symmetry}")
        if self.hotspot_res:
            opts.append(f"ppi.hotspot_res=[{','.join(self.hotspot_res)}]")

        if self.guiding_potentials is not None:
            opts.append(f"potentials.guiding_potentials={self.guiding_potentials}")
        # olig_intra_all / olig_inter_all / olig_custom_contact are top-level
        # potentials.* keys, not token params (RFdiffusion floats the token).
        opts += self._guiding_potential_overrides
        if self.guide_scale is not None:
            opts.append(f"potentials.guide_scale={self.guide_scale}")
        if self.guide_decay is not None:
            opts.append(f"potentials.guide_decay={self.guide_decay}")

        if self.cyclic:
            opts.append("inference.cyclic=True")
        if self.cyc_chains is not None:
            opts.append(f"inference.cyc_chains={self.cyc_chains}")

        if self.provide_seq:
            opts.append(f"contigmap.provide_seq=[{self.provide_seq}]")
        if self.inpaint_str_helix:
            opts.append(f"contigmap.inpaint_str_helix=[{self.inpaint_str_helix}]")
        if self.inpaint_str_strand:
            opts.append(f"contigmap.inpaint_str_strand=[{self.inpaint_str_strand}]")

        if self.noise_scale_ca is not None:
            opts.append(f"denoiser.noise_scale_ca={self.noise_scale_ca}")
        if self.noise_scale_frame is not None:
            opts.append(f"denoiser.noise_scale_frame={self.noise_scale_frame}")

        opts += self._scaffold_options()
        return opts

    def _scaffold_options(self) -> List[str]:
        """Fold-conditioning args (scaffoldguided.*); empty when unused."""
        if not self._scaffold_guided:
            return []
        opts = ["scaffoldguided.scaffoldguided=True"]
        if self.scaffold_dir is not None:
            opts.append(f"scaffoldguided.scaffold_dir={self.scaffold_dir}")
        # target_ss/target_adj describe a target topology -> target_pdb=True;
        # scaffold_dir alone (no target) keeps target_pdb=False.
        if self.target_ss is not None or self.target_adj is not None:
            opts.append("scaffoldguided.target_pdb=True")
            if self.target_path is not None:
                opts.append(f"scaffoldguided.target_path={self.target_path}")
            if self.target_ss is not None:
                opts.append(f"scaffoldguided.target_ss={self.target_ss}")
            if self.target_adj is not None:
                opts.append(f"scaffoldguided.target_adj={self.target_adj}")
        else:
            opts.append("scaffoldguided.target_pdb=False")
        if self.mask_loops is not None:
            opts.append(f"scaffoldguided.mask_loops={self.mask_loops}")
        if self.sampled_insertion is not None:
            opts.append(f"scaffoldguided.sampled_insertion={self.sampled_insertion}")
        if self.sampled_N is not None:
            opts.append(f"scaffoldguided.sampled_N={self.sampled_N}")
        if self.sampled_C is not None:
            opts.append(f"scaffoldguided.sampled_C={self.sampled_C}")
        return opts

    def _config_name_args(self) -> List[str]:
        """Leading --config-name args for modes that need a non-default config."""
        if self.symmetry is not None:
            return ["--config-name", "symmetry"]
        return []

    def _generate_script_run_rfdiffusion(self) -> str:
        """Generate the RFdiffusion execution part of the script.

        With a PDB input, RFdiffusion runs once per input structure inside a
        bash loop over the input ids (resolved at runtime — never ids[0] at
        config time, which breaks under lazy IDs). The per-PDB contigs /
        inpaint / inpaint_str selections (literal broadcast, or a table-column
        reference) are pre-resolved in one pass into contig_options.json by
        pipe_rfdiffusion_contigs.py, and the loop reads each id's values from
        it. Without a PDB it is a single unconditional run named after the
        pipeline.
        """
        common_args = self._common_rfd_options()
        # Bash-array elements: double-quote each so $VAR expands while brackets,
        # commas, and literal " inside a value (guiding_potentials) survive
        # without `eval` re-parsing them.
        common_array = ' '.join('"' + a.replace('"', '\\"') + '"' for a in common_args)
        config_name_array = ' '.join('"' + a + '"' for a in self._config_name_args())
        structures_dir = self.stream_folder("structures")
        prefix = os.path.join(structures_dir, self.pipeline_name)

        if self._scaffold_guided:
            # Fold conditioning: no contigs, no input-pdb loop — the run is
            # driven entirely by the scaffoldguided.* options in `common`.
            return f"""echo "Starting RFdiffusion (fold conditioning)"
echo "Output folder: {self.output_folder}"

cd {self.folders["RFdiffusion"]}
RFD_OPTIONS=(
    {config_name_array}
    "inference.output_prefix={prefix}"
    {common_array}
)
{self.container_prefix()}python {self.inference_py_file} "${{RFD_OPTIONS[@]}}"

"""

        if not self.pdb_stream:
            # Unconditional / symmetric generation — single run, contigs is a
            # literal (validate_params already rejected a table ref without pdb).
            # Inner single quotes keep each bracket entry a Hydra STRING — without
            # them a numeric-only length (e.g. [90]) is parsed as an int and
            # RFdiffusion's .strip() on it raises.
            contigs = self._contigs_arg[1]
            head = [f"contigmap.contigs=['{contigs}']"]
            if self._inpaint_arg[1]:
                head.append(f"contigmap.inpaint_seq=['{self._inpaint_arg[1]}']")
            if self._inpaint_str_arg[1]:
                head.append(f"contigmap.inpaint_str=['{self._inpaint_str_arg[1]}']")
            head_array = ' '.join('"' + a + '"' for a in head)
            mode = "symmetric" if self.symmetry is not None else "unconditional"
            return f"""echo "Starting RFdiffusion ({mode})"
echo "Output folder: {self.output_folder}"

cd {self.folders["RFdiffusion"]}
RFD_OPTIONS=(
    {config_name_array}
    {head_array}
    "inference.output_prefix={prefix}"
    {common_array}
)
{self.container_prefix()}python {self.inference_py_file} "${{RFD_OPTIONS[@]}}"

"""

        # PDB-conditioned: pre-resolve per-PDB selections (one pass over all
        # ids) into a JSON, then loop over the input ids at runtime.
        contigs_arg = _selection_cli_arg(self._contigs_arg)
        inpaint_arg = _selection_cli_arg(self._inpaint_arg)
        inpaint_str_arg = _selection_cli_arg(self._inpaint_str_arg)

        with open(self.contig_args_json, "w") as f:
            json.dump({
                "structures_json": str(self.pdb_ds_json),
                "contigs": contigs_arg,
                "inpaint": inpaint_arg,
                "inpaint_str": inpaint_str_arg,
                "output_json": str(self.contigs_json),
            }, f, indent=2)

        return f"""echo "Resolving per-PDB contig options"
python {self.contigs_py} "{self.contig_args_json}"

echo "Starting RFdiffusion"
echo "Output folder: {self.output_folder}"
cd {self.folders["RFdiffusion"]}

for STRUCT_ID in {Resolve.stream_ids(self.pdb_ds_json)}; do
    INPUT_PDB={Resolve.stream_item(self.pdb_ds_json, '$STRUCT_ID')}
    OPTS=$(python {self.contigs_reader_py} "{self.contigs_json}" "$STRUCT_ID")
    CONTIGS=$(echo "$OPTS" | sed -n '1p')
    INPAINT_SEL=$(echo "$OPTS" | sed -n '2p')
    INPAINT_STR_SEL=$(echo "$OPTS" | sed -n '3p')
    RFD_OPTIONS=(
        {config_name_array}
        "contigmap.contigs=['$CONTIGS']"
        "inference.input_pdb=$INPUT_PDB"
        "inference.output_prefix={structures_dir}/$STRUCT_ID"
        {common_array}
    )
    if [ -n "$INPAINT_SEL" ]; then
        RFD_OPTIONS+=("contigmap.inpaint_seq=['$INPAINT_SEL']")
    fi
    if [ -n "$INPAINT_STR_SEL" ]; then
        RFD_OPTIONS+=("contigmap.inpaint_str=['$INPAINT_STR_SEL']")
    fi
    echo "Running RFdiffusion for $STRUCT_ID"
    {self.container_prefix()}python {self.inference_py_file} "${{RFD_OPTIONS[@]}}"
done

"""

    def _generate_script_create_table(self) -> str:
        """Generate the table creation part of the script."""
        # pipe_rfdiffusion_table.py scans the structures/ stream folder for every
        # produced .pdb/.trb pair — one prefix per input PDB in the multi-PDB
        # case, one for the unconditional case.
        structures_dir = self.stream_folder("structures")
        return f"""echo "Creating results table"
python {self.table_py_file} "{structures_dir}" "{self.main_table}"

"""

    def _generate_script_update_structures_map(self) -> str:
        """Generate script to write structures_map.csv from the actual runtime PDBs."""
        structures_map = self.stream_map_path("structures")
        structures_dir = self.stream_folder("structures")
        # Each design id is "<pdb_id>_<n>"; its parent PDB is the id minus the
        # trailing "_<n>" suffix. Derive the `structures.id` provenance column
        # from that suffix at runtime (works for one or many parent PDBs). With
        # no PDB input there is no parent, so no provenance column.
        prov_arg = ' --provenance-from-suffix "structures.id"' if self.pdb_stream else ""
        return f"""echo "Writing structures map from actual output files"
python {self.update_map_py} --structures-map "{structures_map}" --output-folder "{structures_dir}"{prov_arg}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after RFdiffusion execution."""
        # Pattern-based IDs — PDBs land under <output_folder>/structures/.
        # The per-design map_table is written at runtime by
        # _generate_script_update_structures_map(); here we only declare the
        # stream and its map_table path.
        start = self.design_startnum
        end = self.design_startnum + self.num_designs - 1
        suffix_pattern = f"<{start}..{end}>"
        if self.pdb_stream:
            # One <pdb_id>_<n> fan-out per input PDB. Keep parent ids compact /
            # lazy-safe via the shared multiplier helper.
            structure_ids = generate_multiplied_ids_pattern(
                self.pdb_stream.ids, suffix_pattern,
                input_stream_name="structures"
            )
        else:
            structure_ids = [f"{self.pipeline_name}_{suffix_pattern}"]
        structures_dir = self.stream_folder("structures")
        file_template = [os.path.join(structures_dir, "<id>.pdb")]
        structures_map = self.stream_map_path("structures")

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=file_template,
            map_table=structures_map,
            format="pdb"
        )

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.main_table,
                columns=["id", "pdb", "fixed", "designed", "source_fixed", "plddt_mean", "status"],
                description="RFdiffusion structure generation results"
            )
        }

        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "rfd_params": {
                "pdb_input_ids": self.pdb_stream.ids if self.pdb_stream else None,
                "contigs": str(self.contigs),
                "inpaint": str(self.inpaint),
                "inpaint_str": str(self.inpaint_str),
                "num_designs": self.num_designs,
                "active_site": self.active_site,
                "steps": self.steps,
                "partial_steps": self.partial_steps,
                "reproducible": self.reproducible,
                "design_startnum": self.design_startnum,
                "hotspot_res": self.hotspot_res,
                "symmetry": self.symmetry,
                "guiding_potentials": self.guiding_potentials,
                "guide_scale": self.guide_scale,
                "guide_decay": self.guide_decay,
                "cyclic": self.cyclic,
                "cyc_chains": self.cyc_chains,
                "provide_seq": self.provide_seq,
                "noise_scale_ca": self.noise_scale_ca,
                "noise_scale_frame": self.noise_scale_frame,
                "inpaint_str_helix": self.inpaint_str_helix,
                "inpaint_str_strand": self.inpaint_str_strand,
                "scaffold_dir": self.scaffold_dir,
                "target_ss": self.target_ss,
                "target_adj": self.target_adj,
                "target_path": self.target_path,
                "mask_loops": self.mask_loops,
                "sampled_insertion": self.sampled_insertion,
                "sampled_N": self.sampled_N,
                "sampled_C": self.sampled_C,
            }
        })
        return base_dict
