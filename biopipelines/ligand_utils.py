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


def auth_ligand_field(row) -> tuple:
    """Pick whether a ligand should be expressed by CCD code or SMILES.

    `row` is a compounds-stream row (dict-like with `ccd`, `smiles`, `source`).
    A CCD code names a specific authoritative chemistry; SMILES is preferred only
    when the row's SMILES is known to differ from that CCD. RCSB-fetched ligands
    and streams marked ``format=ccd`` carry both from the same entry, so they
    correspond (use the CCD).
    """
    ccd = str(row.get("ccd", "") or "").strip()
    smiles = str(row.get("smiles", "") or "").strip()
    source = str(row.get("source", "") or "").strip().lower()
    fmt = str(row.get("format", "") or "").strip().lower()
    if ccd and smiles:
        return ("ccd", ccd) if source == "rcsb" or fmt == "ccd" else ("smiles", smiles)
    if ccd:
        return ("ccd", ccd)
    if smiles:
        return ("smiles", smiles)
    raise ValueError(f"ligand row has neither ccd nor smiles: {dict(row)}")


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
