# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Shared runtime helper for resolving a ligand residue code from a compounds
stream.

Tools that consume a ligand `code` accept any compounds stream (Ligand, Boltz2,
Load, custom), so the residue-code contract (1-5 alphanumeric, extended CCD) is
re-checked here at the consumer boundary rather than trusting the producer. By
default a single distinct code is required — a tool that genuinely supports
several ligands passes ``allow_multiple=True`` and gets the full set back.
"""

import os
import re

# Extended CCD codes are 1-5 alphanumeric (matches Ligand's _CCD_CODE_RE).
_CODE_RE = re.compile(r'^[A-Za-z0-9]{1,5}$')


def _read_codes(ligand_json: str):
    # Imported lazily so the pure helpers (auth_ligand_field) stay importable in
    # the pipe-script standalone mode without dragging in the datastream chain.
    try:
        from .biopipelines_io import load_datastream
    except ImportError:
        from biopipelines_io import load_datastream  # type: ignore
    ds = load_datastream(ligand_json)
    # Read codes straight from the rows the map_table actually carries. A
    # compounds stream may declare ids that an upstream filter dropped (the
    # static declaration outlives the runtime survivors); walking ids_expanded
    # would KeyError on those. The distinct codes present are all this needs.
    map_data = ds._get_map_data()
    if map_data is None or "code" not in map_data.columns:
        raise ValueError("ligand compounds stream has no `code` column")
    codes = []
    seen = set()
    for raw in map_data["code"].tolist():
        c = "" if raw is None else str(raw).strip()
        if c and c.lower() != "nan" and c not in seen:
            seen.add(c)
            codes.append(c)
    if not codes:
        raise ValueError("ligand compounds stream `code` column is empty")
    for c in codes:
        if not _CODE_RE.match(c):
            raise ValueError(f"invalid ligand code {c!r}: a residue code must be 1-5 alphanumeric characters")
    return codes


def resolve_ligand_smiles(ligand_json: str):
    """Return the ligand SMILES from a compounds-stream JSON, or None if absent.

    Unlike the `code`, SMILES is optional (code-only ligands have none). Returns
    the first non-empty `smiles` value; used as a bond-order template by tools
    that need correct ligand chemistry from coordinates (PoseBusters, etc.).
    """
    try:
        from .biopipelines_io import load_datastream, iterate_values
    except ImportError:
        from biopipelines_io import load_datastream, iterate_values  # type: ignore
    try:
        ds = load_datastream(ligand_json)
        for _cid, values in iterate_values(ds, columns=["smiles"]):
            s = str(values.get("smiles", "") or "").strip()
            if s:
                return s
    except Exception:
        pass
    return None


def resolve_ligand_code(ligand_json: str) -> str:
    """Resolve the single ligand `code` from a compounds-stream JSON.

    Errors if the stream carries more than one distinct code — use
    ``resolve_ligand_codes`` for tools that support several ligands.
    """
    codes = _read_codes(ligand_json)
    if len(codes) > 1:
        raise ValueError(
            f"ligand compounds stream has {len(codes)} distinct codes ({codes}); expected one")
    return codes[0]


def resolve_ligand_codes(ligand_json: str):
    """Resolve all distinct ligand `code`s from a compounds-stream JSON.

    For tools that genuinely operate on several ligands at once (e.g. PLIP's
    multi-ligand filter). Each code is validated against the residue-code shape.
    """
    return _read_codes(ligand_json)


def _clean(row, key) -> str:
    """A compounds-row field, with serialized-NaN ('', 'nan', 'none') as empty."""
    v = str(row.get(key, "") or "").strip()
    return "" if v.lower() in ("", "nan", "none") else v


def ligand_chemistry(row) -> dict:
    """The chemistry representations a compounds row actually carries.

    Returns ``{"ccd": <code or "">, "smiles": <smiles or "">}`` with serialized
    NaN normalized to empty. The presence test for any chemistry-consuming tool:
    a row has usable chemistry iff at least one value is non-empty.
    """
    return {"ccd": _clean(row, "ccd"), "smiles": _clean(row, "smiles")}


def auth_ligand_field(row) -> tuple:
    """Pick whether a ligand should be expressed by CCD code or SMILES.

    `row` is a compounds-stream row (dict-like with `ccd`, `smiles`, `source`).
    A CCD code names a specific authoritative chemistry; SMILES is preferred only
    when the row's SMILES is known to differ from that CCD. RCSB-fetched ligands
    and streams marked ``format=ccd`` carry both from the same entry, so they
    correspond (use the CCD).

    Raises with an actionable message when the row carries no chemistry — a
    code-only ``Ligand(code=...)`` is a structural HETATM label, not chemistry;
    fetch the CCD with ``Ligand("<code>")`` or pass ``smiles=`` instead.
    """
    chem = ligand_chemistry(row)
    ccd, smiles = chem["ccd"], chem["smiles"]
    source = str(row.get("source", "") or "").strip().lower()
    fmt = str(row.get("format", "") or "").strip().lower()
    if ccd and smiles:
        return ("ccd", ccd) if source == "rcsb" or fmt == "ccd" else ("smiles", smiles)
    if ccd:
        return ("ccd", ccd)
    if smiles:
        return ("smiles", smiles)
    code = str(row.get("code", "") or "").strip()
    raise ValueError(
        f"ligand {code or '?'!r} has no chemistry (no smiles, no ccd): a code-only "
        f"Ligand(code=...) is a structural HETATM label, not a molecule. To give a "
        f"tool real chemistry, fetch the CCD with Ligand({code or '<code>'!r}) or pass "
        f"Ligand(smiles=...).")


def templated_ligand_mol(coord_path: str, smiles: str = None):
    """Read a coordinate ligand (PDB/mol2) and return a sanitized RDKit mol with
    correct bond orders, KEEPING its coordinates. THE single shared bond-order
    implementation — used by write_ligand_sdf and by tools that need the mol
    object directly (e.g. ProLIF) rather than an SDF on disk.

    A bare coordinate ligand lacks bond orders, which breaks RDKit sanitization
    and downstream graph builders. When a `smiles` template is available, assign
    bond orders from it (matching by heavy-atom skeleton, coords preserved),
    trying the bare template first then AddHs-both as a fallback. Without a
    template, fall back to RDKit's own perception.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    ext = os.path.splitext(coord_path)[1].lower()
    if ext == ".pdb":
        mol = Chem.MolFromPDBFile(coord_path, removeHs=False, sanitize=False)
    elif ext == ".mol2":
        mol = Chem.MolFromMol2File(coord_path, removeHs=False, sanitize=False)
    else:
        raise RuntimeError(f"unsupported ligand coordinate format: {ext}")
    if mol is None:
        raise RuntimeError(f"RDKit could not read ligand {coord_path}")

    if smiles:
        template = Chem.MolFromSmiles(smiles)
        if template is None:
            raise ValueError(f"RDKit failed to parse template SMILES: {smiles!r}")
        # Crystal/extracted PDBs usually strip hydrogens, so try the bare template
        # first; fall back to AddHs on both sides for ligands that DO carry Hs.
        assigned = None
        errs = []
        for tmpl, lig in ((template, mol),
                          (Chem.AddHs(template), Chem.AddHs(mol, addCoords=True))):
            try:
                assigned = AllChem.AssignBondOrdersFromTemplate(tmpl, lig)
                break
            except Exception as exc:  # ValueError("No matching found"), etc.
                errs.append(exc)
        if assigned is None:
            # Keep both attempts' errors — for large/dye-sized ligands the
            # bare-template failure is often the more informative one.
            detail = "; ".join(f"[{label}] {e}" for label, e in zip(("bare", "AddHs"), errs))
            raise RuntimeError(
                f"AssignBondOrdersFromTemplate failed for {coord_path}: {detail}. "
                f"Check that the SMILES matches the ligand's skeleton.")
        mol = assigned

    Chem.SanitizeMol(mol)
    return mol


def write_ligand_sdf(src_path: str, dst_sdf: str, smiles: str = None) -> None:
    """Write a chemically valid SDF from a coordinate ligand, KEEPING its
    coordinates. THE single shared file→SDF implementation — OpenBabel routes
    through it, and the scorer/staging tools that consume a ligand SDF rely on it
    instead of each reimplementing the conversion.

    `.sdf` / `.mol` inputs already carry bond orders, so they are copied as-is;
    `.pdb` / `.mol2` are bond-order-restored via templated_ligand_mol.
    """
    import shutil
    from rdkit import Chem

    ext = os.path.splitext(src_path)[1].lower()
    if ext in (".sdf", ".mol"):
        shutil.copyfile(src_path, dst_sdf)
        return
    mol = templated_ligand_mol(src_path, smiles)
    writer = Chem.SDWriter(dst_sdf)
    writer.write(mol)
    writer.close()


def posed_ligand_mol(coord_path: str, smiles: str = None):
    """Read a posed ligand from any supported coordinate format and return a
    sanitized, H-complete RDKit mol keeping its coordinates. Extends
    templated_ligand_mol to sdf/mol (which already carry bond orders) and adds
    the hydrogens a force field needs, placed on the existing heavy-atom frame.
    """
    from rdkit import Chem

    ext = os.path.splitext(coord_path)[1].lower()
    if ext in (".sdf", ".mol"):
        mol = Chem.MolFromMolFile(coord_path, removeHs=False)
        if mol is None:
            raise RuntimeError(f"RDKit could not read ligand {coord_path}")
        if smiles:
            from rdkit.Chem import AllChem
            template = Chem.MolFromSmiles(smiles)
            if template is None:
                raise ValueError(f"RDKit failed to parse template SMILES: {smiles!r}")
            mol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    else:
        mol = templated_ligand_mol(coord_path, smiles)

    return Chem.AddHs(mol, addCoords=True)


# The rotatable-bond definition used to pick which torsions to restrain.
_ROTATABLE_SMARTS = "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"


def _rotatable_bonds(mol):
    from rdkit import Chem
    patt = Chem.MolFromSmarts(_ROTATABLE_SMARTS)
    return [tuple(m) for m in mol.GetSubstructMatches(patt)]


def _named_bonds(mol, atom_name_pairs):
    """Resolve (atom_name, atom_name) pairs against the mol's PDB atom names."""
    by_name = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is not None:
            by_name.setdefault(info.GetName().strip(), atom.GetIdx())

    bonds = []
    for a_name, b_name in atom_name_pairs:
        if a_name not in by_name or b_name not in by_name:
            raise ValueError(
                f"restrain_bonds atom name not found in structure: {a_name!r}-{b_name!r} "
                f"(available: {sorted(by_name)})")
        i, j = by_name[a_name], by_name[b_name]
        if mol.GetBondBetweenAtoms(i, j) is None:
            raise ValueError(f"restrain_bonds pair {a_name!r}-{b_name!r} is not a bond")
        bonds.append((i, j))
    return bonds


def _torsion_for_bond(mol, i, j):
    """Pick a heavy-atom torsion (h, i, j, k) spanning the bond i-j."""
    def _neighbor(center, exclude):
        cands = [n.GetIdx() for n in mol.GetAtomWithIdx(center).GetNeighbors()
                 if n.GetIdx() != exclude and n.GetAtomicNum() > 1]
        if not cands:
            cands = [n.GetIdx() for n in mol.GetAtomWithIdx(center).GetNeighbors()
                     if n.GetIdx() != exclude]
        return cands[0] if cands else None

    h = _neighbor(i, j)
    k = _neighbor(j, i)
    if h is None or k is None:
        return None
    return (h, i, j, k)


def _force_field(mol, ff: str):
    """Build a force field for `mol`, returning (ff_object_factory, engine_name).

    The factory is called to (re)build a field on the same mol — needed because a
    restrained minimisation must be scored on a field WITHOUT the restraint term.
    """
    from rdkit.Chem import AllChem

    if ff in ("auto", "mmff"):
        props = AllChem.MMFFGetMoleculeProperties(mol)
        if props is not None:
            return (lambda m: AllChem.MMFFGetMoleculeForceField(m, AllChem.MMFFGetMoleculeProperties(m)),
                    "MMFF")
        if ff == "mmff":
            raise RuntimeError("MMFF cannot type this molecule (MMFFGetMoleculeProperties returned None)")

    field = AllChem.UFFGetMoleculeForceField(mol)
    if field is None:
        raise RuntimeError("UFF cannot type this molecule")
    return (lambda m: AllChem.UFFGetMoleculeForceField(m), "UFF")


def conformer_strain(mol, restrain_bonds=None, ff: str = "auto", max_iters: int = 2000):
    """Torsional strain of a posed conformer, in kcal/mol.

    The reference state is a TORSION-RESTRAINED minimum, not the pose itself.
    Predicted/docked coordinates carry bond lengths and angles that differ from
    force-field ideals; scoring the raw pose charges every ligand the same large
    constant for that mismatch (tens of kcal/mol) and buries the conformational
    signal. Restraining the rotatable torsions at their posed values lets bonds
    and angles relax while the conformation is held, so the difference against a
    free local minimisation from the same coordinates isolates torsional strain.
    The relaxed reference is deliberately the NEAREST local minimum — a global
    conformer search would add a constant bulk offset to every pose.

    Returns (e_pose, e_relaxed, strain, ff_engine).
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolTransforms

    make_ff, engine = _force_field(mol, ff)

    bonds = _named_bonds(mol, restrain_bonds) if restrain_bonds else _rotatable_bonds(mol)

    restrained = Chem.Mol(mol)
    field = make_ff(restrained)
    conf = restrained.GetConformer()
    for i, j in bonds:
        torsion = _torsion_for_bond(restrained, i, j)
        if torsion is None:
            continue
        angle = rdMolTransforms.GetDihedralDeg(conf, *torsion)
        if engine == "MMFF":
            field.MMFFAddTorsionConstraint(*torsion, False, angle, angle, 1000.0)
        else:
            field.UFFAddTorsionConstraint(*torsion, False, angle, angle, 1000.0)
    _minimize(field, max_iters, "restrained")
    # Score without the restraint term: rebuild a clean field on the minimised mol.
    e_pose = make_ff(restrained).CalcEnergy()

    relaxed = Chem.Mol(mol)
    relaxed_field = make_ff(relaxed)
    _minimize(relaxed_field, max_iters, "free")
    e_relaxed = relaxed_field.CalcEnergy()

    return float(e_pose), float(e_relaxed), float(e_pose - e_relaxed), engine


def _minimize(field, max_iters: int, which: str) -> None:
    """Minimize to convergence. RDKit returns non-zero when it hits max_iters, and the
    energy at that point is not a minimum — so the difference of the two would not be a
    strain. Raise rather than emit an unreliable score."""
    if field.Minimize(maxIts=max_iters) != 0:
        raise RuntimeError(
            f"{which} minimisation did not converge in {max_iters} iterations; "
            f"strain would be unreliable (raise max_iters if the molecule is large)")
