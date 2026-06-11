# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Boltz2 configuration for protein-ligand complex prediction.

Handles apo and holo structure prediction with MSA caching,
ligand binding affinity calculation, and comprehensive analysis.
"""

import os
import sys
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import generate_combinatorics_config, get_mode, predict_output_ids, predict_output_ids_with_provenance, Bundle, Each
    from .datastream_resolver import resolve_input_to_datastream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_combinatorics_config, get_mode, predict_output_ids, predict_output_ids_with_provenance, Bundle, Each
    from datastream_resolver import resolve_input_to_datastream


class Boltz2(BaseConfig):
    """
    Boltz2 configuration for protein-ligand complex prediction.

    Predicts both apo (protein-only) and holo (protein-ligand) structures
    with automatic MSA management and comprehensive scoring.
    """

    TOOL_NAME = "Boltz2"
    TOOL_VERSION = "1.0"

    @classmethod
    def predict_atom_names(cls, ligand, *, path=None, size=(500, 500), scale=3):
        """Predict and depict the atom names Boltz2 will assign to a ligand.

        Boltz names ligand atoms deterministically from chemistry alone, so the
        names used in ligand-atom constraint tokens (contacts / covalent /
        metal_coord, e.g. ``["B", "C8"]``) can be looked up before any prediction
        runs — no GPU, no Pipeline context.

        Two naming regimes, matching Boltz's parser:

        - SMILES ligands (raw SMILES, ``Ligand(smiles=...)``, a PubChem name/CID/
          CAS lookup, ``CompoundLibrary``): ``AddHs`` then
          ``element.upper() + str(CanonicalRankAtoms + 1)``. Hydrogens take part
          in the ranking (so heavy-atom indices are NOT 1..N) but are hidden from
          the picture.
        - CCD ligands (``Ligand("ATP")`` and any RCSB CCD code): atom names come
          from the CCD definition (``_chem_comp_atom.atom_id`` in the RCSB CIF).
          PDB output collapses the residue NAME to ``LIG`` but leaves atom names
          intact.

        Args:
            ligand: a SMILES string, a Ligand / CompoundLibrary instance, or the
                    StandardizedOutput one returns when constructed inside a
                    Pipeline (unwrapped to its tool). No Pipeline context required.
            path: PNG output path. Defaults to ``./<id>_atoms.png`` in the current
                  directory (``./ligand_atoms.png`` for a bare SMILES or a multi-
                  compound grid). The PNG is always written.
            size: per-compound depiction size in pixels.
            scale: supersampling factor; each depiction is rasterized at
                   ``size * scale`` for a crisper, higher-resolution image.

        Returns:
            ``Boltz2.AtomNameResult`` — renders inline in Jupyter, and exposes
            ``.names`` (``{id: {atom_name: serial}}``), ``.path``, ``.regime``.

        Network is used only to resolve a name/CID/CAS or CCD lookup to chemistry;
        a raw SMILES or ``Ligand(smiles=...)`` is fully offline.
        """
        compounds = cls._atom_name_compounds(ligand)

        names: Dict[str, Dict[str, int]] = {}
        regime: Dict[str, str] = {}
        images = []
        legends = []
        for cid, reg, payload in compounds:
            if reg == "smiles":
                img, name_map = cls._draw_smiles_atom_names(payload, size, scale)
            else:
                img, name_map = cls._draw_ccd_atom_names(payload, size, scale)
            names[cid] = name_map
            regime[cid] = reg
            images.append(img)
            legends.append(cid)

        if len(images) == 1:
            final_img = images[0]
            default_name = f"{legends[0]}_atoms.png"
        else:
            final_img = cls._atom_name_grid(images, legends)
            default_name = "ligand_atoms.png"

        out_path = path or os.path.join(os.getcwd(), default_name)
        final_img.save(out_path)
        return cls.AtomNameResult(names=names, image=final_img, path=out_path, regime=regime)

    class AtomNameResult:
        """Result of ``Boltz2.predict_atom_names``.

        Attributes:
            names: ``{compound_id: {atom_name: heavy-atom serial}}`` — the token to
                   put in a Boltz constraint.
            image: a PIL image (single compound) or a stitched grid (multiple).
                   Renders inline in Jupyter via ``_repr_png_``.
            path: the PNG written to disk.
            regime: ``{compound_id: "smiles" | "ccd"}`` — which scheme was used.
        """

        def __init__(self, names, image, path, regime):
            self.names = names
            self.image = image
            self.path = path
            self.regime = regime

        def _repr_png_(self):
            if self.image is not None and hasattr(self.image, "_repr_png_"):
                return self.image._repr_png_()
            return None

        def __repr__(self):
            per = ", ".join(f"{cid}: {len(n)} heavy atoms ({self.regime[cid]})"
                            for cid, n in self.names.items())
            return f"AtomNameResult({per}; png={self.path})"

    @staticmethod
    def _import_pipe_ligand():
        """Import the runtime SMILES/CIF resolvers from pipe_scripts/pipe_ligand.py."""
        pipe_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pipe_scripts")
        if pipe_dir not in sys.path:
            sys.path.insert(0, pipe_dir)
        import pipe_ligand
        return pipe_ligand

    @classmethod
    def _atom_name_compounds(cls, ligand):
        """Normalize ligand input to a list of (id, regime, payload) tuples.

        regime is "smiles" or "ccd"; payload is the SMILES or the CCD code.
        Accepts a SMILES string, a Ligand/CompoundLibrary instance, or the
        StandardizedOutput one returns inside a pipeline (incl. the never-closing
        Colab/Jupyter context). Network lookups resolve a name/CID/CAS to SMILES.
        """
        if isinstance(ligand, str):
            return [("ligand", "smiles", ligand)]

        # Unwrap a StandardizedOutput (.tool) or ToolOutput (.config) to the instance.
        tool_obj = ligand
        if hasattr(ligand, "tool") and hasattr(getattr(ligand, "tool"), "TOOL_NAME"):
            tool_obj = ligand.tool
        elif hasattr(ligand, "config") and hasattr(getattr(ligand, "config"), "TOOL_NAME"):
            tool_obj = ligand.config
        if not hasattr(tool_obj, "TOOL_NAME"):
            raise ValueError(
                "predict_atom_names accepts a SMILES string, a Ligand/"
                "CompoundLibrary instance, or the StandardizedOutput from one; "
                f"got {type(ligand).__name__}.")

        tool = tool_obj.TOOL_NAME
        if tool == "CompoundLibrary":
            comps = getattr(tool_obj, "expanded_compounds", None)
            ids = getattr(tool_obj, "compound_ids", None)
            if not comps or not ids:
                raise ValueError(
                    "CompoundLibrary has no expanded compounds yet. Dict and CDXML "
                    "libraries expand on construction; a CSV-file library only "
                    "expands inside a Pipeline. Build it from a dict/CDXML, or pass "
                    "the SMILES directly.")
            return [(cid, "smiles", c["smiles"]) for cid, c in zip(ids, comps)]

        if tool != "Ligand":
            raise ValueError(
                f"predict_atom_names expects a Ligand or CompoundLibrary, got {tool}.")

        if getattr(tool_obj, "code_only", False):
            raise ValueError(
                "code-only Ligand(code=...) carries no chemistry, so atom names "
                "cannot be predicted. Pass the SMILES or a CCD/name lookup instead.")
        if getattr(tool_obj, "_structures_only", False):
            raise ValueError(
                "Ligand(structures=...) takes atom names from the bound structure, "
                "not from chemistry. Inspect the extracted PDB after running instead.")

        pipe = cls._import_pipe_ligand()
        out = []
        smiles_vals = list(getattr(tool_obj, "smiles_values", []) or [])
        lookup_vals = list(getattr(tool_obj, "lookup_values", []) or [])
        ids = list(getattr(tool_obj, "custom_ids", []) or [])

        # custom_ids ordering is lookup_values first, then smiles_values (Ligand.__init__).
        id_iter = iter(ids)
        for lk in lookup_vals:
            cid = next(id_iter, lk)
            lk_type = tool_obj._detect_lookup_type(lk)
            if lk_type == "ccd":
                out.append((cid, "ccd", lk))
            else:
                smi = pipe.fetch_smiles_from_pubchem(lk, lk_type)
                if not smi:
                    raise ValueError(f"Could not resolve SMILES for lookup {lk!r} from PubChem.")
                out.append((cid, "smiles", smi))
        for smi in smiles_vals:
            cid = next(id_iter, smi)
            out.append((cid, "smiles", smi))

        if not out:
            raise ValueError("Ligand has no resolvable chemistry to draw.")
        return out

    @staticmethod
    def _smiles_atom_names(smiles):
        """Reproduce Boltz's SMILES atom-naming verbatim.

        Returns (rdkit_atom_idx, element, atom_name) for every atom incl. H — the
        ranking must include hydrogens to match Boltz even though they are hidden.
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = AllChem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
        mol = AllChem.AddHs(mol)

        canonical_order = AllChem.CanonicalRankAtoms(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

        out = []
        for atom, can_idx in zip(mol.GetAtoms(), canonical_order):
            atom_name = atom.GetSymbol().upper() + str(can_idx + 1)
            if len(atom_name) > 4:
                raise ValueError(
                    f"{smiles} has an atom with a name longer than 4 characters: "
                    f"{atom_name}. Boltz rejects this ligand.")
            out.append((atom.GetIdx(), atom.GetSymbol(), atom_name))
        return out

    @classmethod
    def _draw_smiles_atom_names(cls, smiles, size, scale=3):
        """2D depiction of a SMILES, heavy atoms annotated by Boltz atom names.

        Returns (PIL image, {atom_name: heavy_serial}).
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
        from PIL import Image
        import io

        named = cls._smiles_atom_names(smiles)
        name_by_idx = {idx: name for idx, _el, name in named}

        mol = AllChem.MolFromSmiles(smiles)
        mol = AllChem.AddHs(mol)
        names_map = {}
        serial = 0
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "H":
                continue
            name = name_by_idx[atom.GetIdx()]
            atom.SetProp("atomNote", name)
            names_map[name] = serial
            serial += 1
        mol = Chem.RemoveHs(mol)
        AllChem.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DCairo(size[0] * scale, size[1] * scale)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return Image.open(io.BytesIO(drawer.GetDrawingText())), names_map

    @staticmethod
    def _read_cif_loop(text, prefix):
        try:
            from .pdb_parser import read_cif_loop
        except ImportError:
            from pdb_parser import read_cif_loop
        return read_cif_loop(text, prefix)

    @classmethod
    def _fetch_ccd_cif(cls, ccd_code):
        """Download an RCSB chem-comp CIF; raise a clear error if absent."""
        import requests

        code = ccd_code.strip().upper()
        url = f"https://files.rcsb.org/ligands/download/{code}.cif"
        headers = {"User-Agent": "BioPipelines-Ligand/1.0 (https://github.com/locbp-uzh/biopipelines)"}
        resp = requests.get(url, headers=headers, timeout=15)
        if resp.status_code == 404:
            raise ValueError(f"No CCD entry {code!r} at RCSB (cannot read atom names).")
        resp.raise_for_status()
        return resp.text

    # mmCIF bond value_order -> RDKit bond type.
    _CIF_BOND_ORDER = {"SING": 1.0, "DOUB": 2.0, "TRIP": 3.0, "QUAD": 4.0, "AROM": 1.5}

    @classmethod
    def _build_ccd_mol(cls, ccd_code):
        """Build an RDKit mol straight from the CCD's own atoms+bonds+coords.

        Names attach to atoms by construction (no SMILES re-ordering), so the
        depiction is guaranteed to label the right atom. Returns (mol, names_map)
        where mol carries heavy atoms only, each with its CCD ``atom_id`` set as
        ``atomNote``, laid out from the ideal 2-D projection of the CCD coords.
        """
        from rdkit import Chem
        from rdkit.Geometry import Point3D

        text = cls._fetch_ccd_cif(ccd_code)
        atoms = cls._read_cif_loop(text, "_chem_comp_atom.")
        bonds = cls._read_cif_loop(text, "_chem_comp_bond.")
        if not atoms:
            raise ValueError(f"Could not parse atoms from CCD CIF for {ccd_code!r}.")

        # Heavy atoms only, in CIF order.
        heavy = [a for a in atoms if a.get("type_symbol", "").upper() != "H"]
        idx_by_id = {}
        rw = Chem.RWMol()
        conf_xyz = []
        for a in heavy:
            aid = a["atom_id"]
            el = a["type_symbol"]
            atom = Chem.Atom(el.capitalize())
            charge = a.get("charge", "0")
            try:
                atom.SetFormalCharge(int(float(charge)))
            except (TypeError, ValueError):
                pass
            atom.SetProp("atomNote", aid)
            new_idx = rw.AddAtom(atom)
            idx_by_id[aid] = new_idx
            x = float(a.get("pdbx_model_Cartn_x_ideal", a.get("model_Cartn_x", 0.0)) or 0.0)
            y = float(a.get("pdbx_model_Cartn_y_ideal", a.get("model_Cartn_y", 0.0)) or 0.0)
            z = float(a.get("pdbx_model_Cartn_z_ideal", a.get("model_Cartn_z", 0.0)) or 0.0)
            conf_xyz.append((x, y, z))

        for b in bonds:
            a1, a2 = b.get("atom_id_1"), b.get("atom_id_2")
            if a1 not in idx_by_id or a2 not in idx_by_id:
                continue  # bond to a hydrogen (dropped) — skip
            order = cls._CIF_BOND_ORDER.get(b.get("value_order", "SING"), 1.0)
            bt = Chem.BondType.AROMATIC if order == 1.5 else Chem.BondType.values[int(order)]
            rw.AddBond(idx_by_id[a1], idx_by_id[a2], bt)

        mol = rw.GetMol()
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y, z) in enumerate(conf_xyz):
            conf.SetAtomPosition(i, Point3D(x, y, z))
        mol.AddConformer(conf, assignId=True)
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            # Layout/labels don't need a fully sanitized mol; keep going.
            pass

        names_map = {aid: i for aid, i in idx_by_id.items()}
        return mol, names_map

    @classmethod
    def _draw_ccd_atom_names(cls, ccd_code, size, scale=3):
        """Depiction for a CCD ligand, annotated by the CCD's own atom names."""
        from rdkit.Chem.Draw import rdMolDraw2D
        from rdkit.Chem import rdDepictor
        from PIL import Image
        import io

        mol, names_map = cls._build_ccd_mol(ccd_code)
        # Flatten the CCD's 3-D ideal coords to a clean 2-D layout for display.
        rdDepictor.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DCairo(size[0] * scale, size[1] * scale)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return Image.open(io.BytesIO(drawer.GetDrawingText())), names_map

    @staticmethod
    def _atom_name_grid(images, legends):
        """Stitch per-compound images side by side with id legends underneath."""
        from PIL import Image, ImageDraw, ImageFont

        max_h = max(im.height for im in images)
        font_px = max(14, max_h // 22)
        try:
            font = ImageFont.truetype("DejaVuSans.ttf", font_px)
        except OSError:
            font = ImageFont.load_default()
        pad = max(8, max_h // 50)
        legend_h = font_px + pad
        w = sum(im.width for im in images) + pad * (len(images) + 1)
        h = max_h + legend_h + pad * 2
        canvas = Image.new("RGB", (w, h), "white")
        draw = ImageDraw.Draw(canvas)
        x = pad
        for im, legend in zip(images, legends):
            canvas.paste(im, (x, pad))
            draw.text((x, im.height + pad), str(legend), fill="black", font=font)
            x += im.width + pad
        return canvas

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("Boltz2Env", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "Boltz2 already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("Boltz2Env", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("Boltz2Env", env_manager, biopipelines)
        return f"""echo "=== Installing Boltz2 ==="
{skip}{remove_block}
{env_block}

# Verify installation
if {env_manager} run -n Boltz2Env python -c "import boltz" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== Boltz2 installation complete ==="
else
    echo "ERROR: Boltz2 verification failed (cannot import boltz)"
    exit 1
fi
"""

    # Path descriptors - lazy evaluation after output_folder is set.
    #
    # New layout:
    #   configuration/  — input YAMLs for boltz predict, combinatorics JSON,
    #                     sequence_ids.csv, the ligands-from-SMILES helper CSV.
    #   execution/      — raw boltz_results_* dumps go here.
    #   structures/     — predicted PDB/CIF + structures_map.csv.
    #   sequences/      — content-bearing stream CSV (sequences.csv acts as map).
    #   msas/           — MSA files (.csv or .a3m) + msas_map.csv.
    #   tables/         — standalone TableInfo CSVs: confidence, affinity, missing.
    #   _extras/        — scores_info.txt and other ancillary dumps.
    config_files_dir = Path(lambda self: self.configuration_path("config_files"))

    # Configuration-time artifacts
    queries_csv = Path(lambda self: self.configuration_path(f"{self.get_effective_job_name() or 'prediction'}_queries.csv"))
    queries_fasta = Path(lambda self: self.configuration_path(f"{self.get_effective_job_name() or 'prediction'}_queries.fasta"))
    expanded_library_csv = Path(lambda self: self.configuration_path("expanded_smiles_library.csv"))
    fasta_files_list_file = Path(lambda self: self.configuration_path(".input_fasta_files.txt"))
    sequence_ids_file = Path(lambda self: self.configuration_path("sequence_ids.csv"))
    combinatorics_config_file = Path(lambda self: self.configuration_path("combinatorics_config.json"))

    # Raw boltz dumps
    prediction_folder = Path(lambda self: self.execution_folder)
    msa_cache_folder = Path(lambda self: self.stream_folder("msas"))

    # Stream maps (lineage)
    structures_map_csv = Path(lambda self: self.stream_map_path("structures"))
    msas_csv = Path(lambda self: self.stream_map_path("msas"))
    # Emitted compounds map: same ids/SMILES as the input ligands, but with the
    # `code` column overwritten with the residue code Boltz assigns (CCD code
    # for CCD ligands, "LIG" for SMILES ligands). Written at runtime.
    compounds_map_csv = Path(lambda self: self.stream_map_path("compounds"))

    # Content-bearing stream (sequences.csv IS the content table)
    sequences_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))

    # Standalone TableInfo CSVs
    confidence_csv = Path(lambda self: self.table_path("confidence"))
    affinity_csv = Path(lambda self: self.table_path("affinity"))
    missing_csv = Path(lambda self: self.table_path("missing"))

    # Extras (informational)
    scores_info_file = Path(lambda self: os.path.join(self.extras_folder, "scores_info.txt"))

    # Helper script paths
    boltz_config_unified_py = Path(lambda self: self.pipe_script_path("pipe_boltz_config_unified.py"))
    boltz_postprocessing_py = Path(lambda self: self.pipe_script_path("pipe_boltz_postprocessing.py"))
    boltz_msa_copy_py = Path(lambda self: self.pipe_script_path("pipe_boltz_msa_copy.py"))
    boltz_compounds_py = Path(lambda self: self.pipe_script_path("pipe_boltz_compounds.py"))

    def __init__(self,
                 # Primary input parameters
                 config: Optional[str] = None,
                 proteins: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ssDNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 dsDNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ssRNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 dsRNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ligands: Optional[Union[DataStream, StandardizedOutput]] = None,
                 msas: Optional[StandardizedOutput] = None,
                 # Core prediction parameters
                 affinity: bool = True,
                 output_format: str = "pdb",
                 # Advanced prediction parameters
                 recycling_steps: Optional[int] = None,
                 diffusion_samples: Optional[int] = None,
                 top_only: bool = True,
                 use_potentials: bool = False,
                 # Template parameters
                 template: Optional[str] = None,
                 template_chain_ids: Optional[List[str]] = None,
                 template_force: bool = True,
                 template_threshold: float = 5.0,
                 # Pocket constraint parameters
                 pocket_residues: Optional[List[int]] = None,
                 pocket_max_distance: float = 7.0,
                 pocket_force: bool = True,
                 # Glycosylation parameters
                 glycosylation: Optional[Dict[str, List[int]]] = None,
                 # Covalent linkage parameters
                 covalent_linkage: Optional[Dict[str, Any]] = None,
                 # Contact constraint parameters
                 contacts: Optional[List[Dict[str, Any]]] = None,
                 # Disulfide bond constraints
                 disulfide_bonds: Optional[List[Dict[str, Any]]] = None,
                 # Metal coordination bond constraints
                 metal_coord: Optional[List[Dict[str, Any]]] = None,
                 **kwargs):
        """
        Initialize Boltz2 configuration.

        Args:
            config: Direct YAML configuration string
            proteins: Protein sequences as DataStream or StandardizedOutput
            ssDNA: Single-stranded DNA sequences (one chain per sequence)
            dsDNA: Double-stranded DNA sequences (two chains per sequence, reverse complement auto-generated)
            ssRNA: Single-stranded RNA sequences (one chain per sequence)
            dsRNA: Double-stranded RNA sequences (two chains per sequence, reverse complement auto-generated)
            ligands: DataStream or StandardizedOutput with a compounds stream (e.g. Ligand, CompoundLibrary)
            msas: Precomputed MSAs (StandardizedOutput with an msas table, e.g. from
                MMseqs2 or a previous Boltz2 run). When provided, the public MSA
                server is not queried. When omitted, MSAs are generated via the
                public MSA server.
            affinity: Whether to calculate binding affinity
            output_format: Output format ("pdb" or "mmcif")
            recycling_steps: Number of recycling steps
            diffusion_samples: Number of diffusion samples
            use_potentials: Enable potentials for improved structure prediction
            template: Path to PDB template file
            template_chain_ids: List of chain IDs to apply template to
            template_force: Whether to force template usage
            template_threshold: RMSD threshold for template matching
            pocket_residues: List of residue positions defining binding pocket
            pocket_max_distance: Maximum distance for pocket constraint
            pocket_force: Whether to force pocket constraint
            glycosylation: Dict mapping chain IDs to Asn positions for N-glycosylation
            covalent_linkage: Dict specifying covalent attachment (chain, position,
                      protein_atom, ligand_atom). `position` may be a residue index
                      or a selection string resolved per-structure against the chain's
                      sequence (e.g. "C in ILIPCH" finds the catalytic Cys regardless
                      of renumbering); it must resolve to exactly one residue.
            contacts: List of contact constraints, each a dict with token1, token2,
                      optional max_distance (4-20A, default 6.0), optional force (bool).
                      Tokens are [chain_id, residue_index] for proteins or
                      [chain_id, atom_name] for ligands (e.g. ["B", "C8"]). Use
                      Boltz2.predict_atom_names(<ligand>) to see the exact ligand
                      atom names before writing the constraint.
            disulfide_bonds: List of cysteine-cysteine bond constraints. Each entry is a
                      dict with token1=[chain, residue] and token2=[chain, residue]. Atom
                      names default to SG/SG and are filled in automatically.
            metal_coord: List of metal-coordination bond constraints. Each entry is a dict
                      with atom1=[chain, residue, atom] and atom2=[chain, residue, atom]
                      (the ligand_residue defaults to 1 for single-residue metal entities).
            **kwargs: Additional parameters

        Output:
            Streams: structures (.pdb/.mmcif), sequences (.csv), compounds (.csv), msas (.csv/.a3m)
            Tables:
                structures: id | file
                confidence: id | input_file | confidence_score | ptm | iptm | complex_plddt | complex_iplddt
                sequences: id | sequence
                msas: id | sequences.id | sequence | msa_file
                affinity: id | input_file | affinity_pred_value | affinity_probability_binary
                compounds: id | format | code | smiles | ccd  (code = residue code Boltz assigned)
                missing: id | removed_by | kind | cause
        """
        self.config = config
        self.msas = msas

        # Store raw inputs for combinatorics config generation
        self.proteins = proteins
        self.ssDNA = ssDNA
        self.dsDNA = dsDNA
        self.ssRNA = ssRNA
        self.dsRNA = dsRNA
        self.ligands = ligands

        # Combinatorics modes
        self._proteins_mode = get_mode(proteins)
        self._ssDNA_mode = get_mode(ssDNA)
        self._dsDNA_mode = get_mode(dsDNA)
        self._ssRNA_mode = get_mode(ssRNA)
        self._dsRNA_mode = get_mode(dsRNA)
        self._ligands_mode = get_mode(ligands)

        # Resolve sequence inputs to DataStreams using shared utility
        self.proteins_stream: Optional[DataStream] = resolve_input_to_datastream(proteins, fallback_stream="sequences")
        self.ssDNA_stream: Optional[DataStream] = resolve_input_to_datastream(ssDNA, fallback_stream="sequences")
        self.dsDNA_stream: Optional[DataStream] = resolve_input_to_datastream(dsDNA, fallback_stream="sequences")
        self.ssRNA_stream: Optional[DataStream] = resolve_input_to_datastream(ssRNA, fallback_stream="sequences")
        self.dsRNA_stream: Optional[DataStream] = resolve_input_to_datastream(dsRNA, fallback_stream="sequences")

        # Resolve ligand input to a compounds stream.
        self.ligands_stream: Optional[DataStream] = None
        if ligands is not None:
            self.ligands_stream = resolve_input_to_datastream(ligands, fallback_stream="compounds")

        # Override affinity to False if no ligands
        if self.ligands is None and affinity:
            print("Warning: No ligands detected, setting affinity=False")
            self.affinity = False
        else:
            self.affinity = affinity

        self.output_format = output_format
        self.recycling_steps = recycling_steps
        self.diffusion_samples = diffusion_samples
        # When False, every diffusion sample is surfaced as a separate structure
        # (<id>_1..N) instead of just the top model (<id>) — e.g. to feed a pose
        # ensemble into downstream design. Default True keeps the single-best output.
        self.top_only = top_only
        self.use_potentials = use_potentials

        # Template parameters
        self.template = template
        self.template_chain_ids = template_chain_ids
        self.template_force = template_force
        self.template_threshold = template_threshold

        # Pocket constraint parameters
        self.pocket_residues = pocket_residues
        self.pocket_max_distance = pocket_max_distance
        self.pocket_force = pocket_force

        # Glycosylation, covalent linkage, and contacts
        self.glycosylation = glycosylation
        self.covalent_linkage = covalent_linkage
        self.contacts = contacts
        self.disulfide_bonds = disulfide_bonds
        self.metal_coord = metal_coord

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate Boltz2-specific parameters."""
        # Must have some form of input
        has_sequence_input = any([
            self.proteins_stream is not None,
            self.ssDNA_stream is not None,
            self.dsDNA_stream is not None,
            self.ssRNA_stream is not None,
            self.dsRNA_stream is not None,
        ])
        has_input = self.config is not None or has_sequence_input
        if not has_input:
            raise ValueError("Either config or at least one sequence parameter (proteins/ssDNA/dsDNA/ssRNA/dsRNA) is required")

        # Cannot specify both config and sequence parameters
        if self.config is not None and has_sequence_input:
            raise ValueError("Cannot specify both config and sequence parameters (proteins/ssDNA/dsDNA/ssRNA/dsRNA)")

        # Validate enum values
        if self.output_format not in ["pdb", "mmcif"]:
            raise ValueError("output_format must be 'pdb' or 'mmcif'")

        _validate_freeform_string("template", self.template)
        _validate_freeform_string("config", self.config)

        if self.recycling_steps is not None and (not isinstance(self.recycling_steps, int) or self.recycling_steps < 1):
            raise ValueError("recycling_steps must be a positive integer")

        if self.diffusion_samples is not None and (not isinstance(self.diffusion_samples, int) or self.diffusion_samples < 1):
            raise ValueError("diffusion_samples must be a positive integer")
        if not isinstance(self.top_only, bool):
            raise ValueError("top_only must be a bool")

        if self.contacts is not None:
            if not isinstance(self.contacts, list):
                raise ValueError("contacts must be a list of dicts")
            for i, contact in enumerate(self.contacts):
                if not isinstance(contact, dict):
                    raise ValueError(f"contacts[{i}] must be a dict")
                for key in ('token1', 'token2'):
                    if key not in contact:
                        raise ValueError(f"contacts[{i}] missing required field '{key}'")
                    token = contact[key]
                    if not isinstance(token, list) or len(token) not in (2, 3):
                        raise ValueError(f"contacts[{i}]['{key}'] must be a list of 2 or 3 elements [chain, residue_or_atom_name] or [chain, residue, atom]")
                if 'max_distance' in contact:
                    md = contact['max_distance']
                    if not isinstance(md, (int, float)) or md < 4 or md > 20:
                        raise ValueError(f"contacts[{i}]['max_distance'] must be a number between 4 and 20")

        if self.disulfide_bonds is not None:
            if not isinstance(self.disulfide_bonds, list):
                raise ValueError("disulfide_bonds must be a list of dicts")
            for i, bond in enumerate(self.disulfide_bonds):
                if not isinstance(bond, dict):
                    raise ValueError(f"disulfide_bonds[{i}] must be a dict")
                for key in ('token1', 'token2'):
                    if key not in bond:
                        raise ValueError(f"disulfide_bonds[{i}] missing required field '{key}'")
                    token = bond[key]
                    if not isinstance(token, list) or len(token) != 2:
                        raise ValueError(f"disulfide_bonds[{i}]['{key}'] must be [chain, residue]")

        if self.metal_coord is not None:
            if not isinstance(self.metal_coord, list):
                raise ValueError("metal_coord must be a list of dicts")
            for i, entry in enumerate(self.metal_coord):
                if not isinstance(entry, dict):
                    raise ValueError(f"metal_coord[{i}] must be a dict")
                for key in ('atom1', 'atom2'):
                    if key not in entry:
                        raise ValueError(f"metal_coord[{i}] missing required field '{key}'")
                    atom = entry[key]
                    if not isinstance(atom, list) or len(atom) != 3:
                        raise ValueError(f"metal_coord[{i}]['{key}'] must be [chain, residue, atom]")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files."""
        self.folders = pipeline_folders

    def _build_combinatorics_kwargs(self) -> Dict:
        """Build combinatorics kwargs for all active axes."""
        kwargs = {}
        if self.proteins is not None:
            kwargs['proteins'] = (self.proteins, "sequences", "protein")
        if self.ssDNA is not None:
            kwargs['ssDNA'] = (self.ssDNA, "sequences", "ssdna")
        if self.dsDNA is not None:
            kwargs['dsDNA'] = (self.dsDNA, "sequences", "dsdna")
        if self.ssRNA is not None:
            kwargs['ssRNA'] = (self.ssRNA, "sequences", "ssrna")
        if self.dsRNA is not None:
            kwargs['dsRNA'] = (self.dsRNA, "sequences", "dsrna")
        if self.ligands is not None:
            kwargs['ligands'] = (self.ligands, "compounds", "ligand")
        return kwargs

    def _write_combinatorics_config(self) -> str:
        """Write combinatorics config file at configuration time."""
        config_path = self.combinatorics_config_file
        generate_combinatorics_config(config_path, **self._build_combinatorics_kwargs())
        return config_path

    def _generate_extra_config_params(self) -> str:
        """Generate extra CLI parameters for pipe_boltz_config_unified.py."""
        extra_params = []

        if self.template:
            extra_params.append(f'--template "{self.template}"')
            if self.template_chain_ids:
                extra_params.append(f'--template-chains "{",".join(self.template_chain_ids)}"')
            if self.template_force:
                extra_params.append('--template-force')
            extra_params.append(f'--template-threshold {self.template_threshold}')

        if self.pocket_residues:
            extra_params.append(f'--pocket-residues "{self.pocket_residues}"')
            extra_params.append(f'--pocket-max-distance {self.pocket_max_distance}')
            if self.pocket_force:
                extra_params.append('--pocket-force')

        if self.glycosylation:
            glyco_json = json.dumps(self.glycosylation)
            extra_params.append(f"--glycosylation '{glyco_json}'")

        if self.covalent_linkage:
            covalent_json = json.dumps(self.covalent_linkage)
            extra_params.append(f"--covalent-linkage '{covalent_json}'")

        if self.contacts:
            contacts_json = json.dumps(self.contacts)
            extra_params.append(f"--contacts '{contacts_json}'")

        if self.disulfide_bonds:
            disulfide_json = json.dumps(self.disulfide_bonds)
            extra_params.append(f"--disulfide-bonds '{disulfide_json}'")

        if self.metal_coord:
            metal_json = json.dumps(self.metal_coord)
            extra_params.append(f"--metal-coord '{metal_json}'")

        return " ".join(extra_params)

    def _supplied_msa_format(self) -> str:
        """Format of the upstream msas stream (csv/a3m); defaults to csv."""
        try:
            return self.msas.streams.msas.format or "csv"
        except Exception:
            return "csv"

    def _get_msa_table_flag(self) -> str:
        """Get the MSA table flag for the config generator."""
        if not self.msas:
            return ""

        if hasattr(self.msas, 'tables'):
            if hasattr(self.msas.tables, '_tables') and 'msas' in self.msas.tables._tables:
                msa_table_path = self.msas.tables._tables['msas'].info.path
                return f'--msa-table "{msa_table_path}"'
            elif hasattr(self.msas.tables, 'msas'):
                msa_table = self.msas.tables.msas
                if hasattr(msa_table, 'info'):
                    return f'--msa-table "{msa_table.info.path}"'

        return ""

    def generate_script(self, script_path: str) -> str:
        """Generate bash script for Boltz2 execution."""
        boltz_cache_folder = self.folders["BoltzCache"]
        # Query the public MSA server only when no MSAs were supplied.
        msa_option = "" if self.msas is not None else " --use_msa_server"

        script_content = "#!/bin/bash\n"
        script_content += "# Boltz2 execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        # Fix NVRTC builtins path for CUDA 13.0 (Colab compatibility)
        script_content += """# NVRTC builtins fix for CUDA 13.0
CU13_LIB=$(find "$CONDA_PREFIX" -path "*/nvidia/cu13/lib" -type d 2>/dev/null | head -1)
if [ -n "$CU13_LIB" ] && [ -d "$CU13_LIB" ]; then
    export LD_LIBRARY_PATH="$CU13_LIB:${LD_LIBRARY_PATH:-}"
    echo "Added $CU13_LIB to LD_LIBRARY_PATH for NVRTC compatibility"
fi

"""

        # execution/, configuration/, and stream folders are auto-created
        # by the pipeline. msas/ gets populated at runtime below.

        # Handle MSA recycling if provided
        if self.msas:
            script_content += self._generate_msa_recycling_section()

        # Write combinatorics config (into configuration/)
        combinatorics_config_path = self._write_combinatorics_config()

        effective_job_name = self.get_effective_job_name() or "prediction"
        config_file_path = self.configuration_path(f"{effective_job_name}.yaml")

        if self.config:
            # Direct YAML configuration — write at configuration time.
            with open(config_file_path, 'w') as f:
                f.write(self.config)
            script_content += f"""
echo "Using direct YAML configuration: {config_file_path}"

"""
            uses_unified_config = False
        else:
            # Use unified config generation
            script_content += self._generate_unified_config_section(combinatorics_config_path, effective_job_name)
            uses_unified_config = True

        # Build Boltz2 options — boltz writes boltz_results_* under
        # execution_folder (raw model dumps).
        boltz_options = f"--cache {boltz_cache_folder} --out_dir {self.execution_folder}{msa_option} --output_format {self.output_format}"

        if self.recycling_steps is not None:
            boltz_options += f" --recycling_steps {self.recycling_steps}"

        if self.diffusion_samples is not None:
            boltz_options += f" --diffusion_samples {self.diffusion_samples}"

        if self.use_potentials:
            boltz_options += " --use_potentials"

        # Run Boltz2 prediction (wrapped in container_prefix when configured).
        # boltz predict accepts a directory and iterates its .yaml/.fasta
        # files in a single process, avoiding the per-file startup cost of
        # a bash loop.
        cp = self.container_prefix()
        predict_input = self.config_files_dir if uses_unified_config else config_file_path
        script_content += f"""
echo "Running Boltz2 prediction"
{cp}boltz predict {predict_input} {boltz_options}

"""

        # Post-process results
        script_content += self._generate_postprocess_section()

        # Emit the compounds stream with the residue code Boltz assigned
        script_content += self._generate_compounds_section()

        # Propagate missing table
        script_content += self._generate_missing_table_propagation()

        # Completion check
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_unified_config_section(self, combinatorics_config_path: str, config_name: str) -> str:
        """Generate script section for unified config generation."""
        script = f"""
echo "Generating Boltz2 configurations using unified config generator"
mkdir -p {self.config_files_dir}

"""
        # Build command for unified config generator
        cmd_parts = [
            f'python {self.boltz_config_unified_py}',
            f'--combinatorics-config "{combinatorics_config_path}"',
            f'--output-dir "{self.config_files_dir}"',
            f'--sequence-ids "{self.sequence_ids_file}"',
        ]

        msa_table_flag = self._get_msa_table_flag()
        if msa_table_flag:
            cmd_parts.append(msa_table_flag)

        if self.affinity:
            cmd_parts.append('--affinity')

        cmd_parts.append(f'--sequences-csv "{self.sequences_csv}"')

        extra_params = self._generate_extra_config_params()
        if extra_params:
            cmd_parts.append(extra_params)

        script += '# Generate config files\n'
        script += ' \\\n    '.join(cmd_parts) + '\n\n'

        return script

    def _generate_msa_recycling_section(self) -> str:
        """Generate script section for MSA recycling."""
        if not self.msas:
            return ""

        # Get MSA table path from msas input
        msa_table_path = None
        if hasattr(self.msas, 'tables'):
            if hasattr(self.msas.tables, '_tables') and 'msas' in self.msas.tables._tables:
                msa_table_path = self.msas.tables._tables['msas'].info.path
            elif hasattr(self.msas.tables, 'msas'):
                msa_table = self.msas.tables.msas
                if hasattr(msa_table, 'info'):
                    msa_table_path = msa_table.info.path

        if msa_table_path:
            return f"""
echo "Recycling MSAs from previous prediction"
python {self.boltz_msa_copy_py} \\
    --msa-table "{msa_table_path}" \\
    --output-folder "{self.msa_cache_folder}"

"""
        return ""

    def _generate_postprocess_section(self) -> str:
        """Generate script section for post-processing.

        Boltz writes raw ``boltz_results_*`` dumps into execution/. The post-
        processing script copies/renames structures into structures/, MSAs
        into msas/, and writes the standalone confidence/affinity tables
        into tables/. Each destination is passed explicitly so the script
        doesn't have to know the layout convention.
        """
        structures_dir = self.stream_folder("structures")
        return f"""
echo "Post-processing Boltz2 results"
python {self.boltz_postprocessing_py} {self.execution_folder} {self.output_folder} {self.sequence_ids_file} \\
    --structures-map {self.structures_map_csv} \\
    --structures-folder "{structures_dir}" \\
    --msas-folder "{self.msa_cache_folder}" \\
    --confidence-csv "{self.confidence_csv}" \\
    --affinity-csv "{self.affinity_csv}" \\
    --sequences-csv "{self.sequences_csv}" \\
    --msas-csv "{self.msas_csv}" \\
    --scores-info "{self.scores_info_file}"{"" if self.top_only else " --all-samples"}

if [ $? -ne 0 ]; then
    echo "Error: Post-processing failed"
    exit 1
fi

echo "Post-processing completed"

"""

    def _generate_compounds_section(self) -> str:
        """Emit the compounds stream map with the residue code Boltz assigned.

        Boltz writes SMILES ligands into the complex with the residue name
        ``LIG`` and CCD ligands with their CCD code. The pipe script reads the
        input ligand metadata (id, format, smiles, ccd) and writes a fresh
        compounds map that preserves id/smiles but sets ``code`` accordingly,
        so downstream HETATM-selector tools read the right code automatically.
        """
        if not self.ligands_stream:
            return ""
        return f"""
echo "Emitting compounds stream with assigned residue codes"
python {self.boltz_compounds_py} \\
    --combinatorics-config "{self.combinatorics_config_file}" \\
    --output-format "{self.output_format}" \\
    --output-compounds "{self.compounds_map_csv}"

"""

    def _missing_input_sources(self):
        """All id-bearing input axes that may carry an upstream `missing` table."""
        return (
            self.proteins, self.proteins_stream,
            self.ligands, self.ligands_stream,
            self.ssDNA, self.ssDNA_stream,
            self.dsDNA, self.dsDNA_stream,
            self.ssRNA, self.ssRNA_stream,
            self.dsRNA, self.dsRNA_stream,
        )

    def _generate_missing_table_propagation(self) -> str:
        """Propagate the union of upstream `missing` manifests into tables/missing.csv."""
        return self.generate_missing_propagation(
            *self._missing_input_sources(), missing_csv=self.missing_csv
        )


    def _predict_sequence_ids(self) -> List[str]:
        """Predict sequence IDs from input sources using combinatorics module."""
        return predict_output_ids(
            **self._build_combinatorics_kwargs()
        )

    def _predict_sequence_ids_with_provenance(self):
        """Predict sequence IDs and provenance from input sources."""
        return predict_output_ids_with_provenance(
            **self._build_combinatorics_kwargs()
        )

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after Boltz2 execution."""
        predicted_ids, _provenance = self._predict_sequence_ids_with_provenance()

        # Structure files land inside the structures/ stream folder.
        # The per-design map_table (with {stream}.id provenance) is written at
        # runtime by pipe_boltz_postprocessing.py from the combinatorics config;
        # here we only declare the stream and its map_table path.
        structure_ext = ".pdb" if self.output_format == "pdb" else ".cif"
        structure_files = [self.stream_path("structures", f"<id>{structure_ext}")]

        # With top_only=False the postprocess surfaces EVERY diffusion sample as
        # <id>_1..K (for any K, including K=1), so the declared structure ids must
        # mirror that exactly or the completion check mismatches. Suffix whenever
        # top_only is off; K = diffusion_samples (Boltz's effective default is 1).
        # Use the compact lazy-id range pattern (<1..K>) rather than materializing.
        # Other streams (msas/sequences/compounds) stay per-complex on the base ids.
        structure_ids = predicted_ids
        if not self.top_only:
            k = self.diffusion_samples or 1
            structure_ids = [f"{pid}_<1..{k}>" for pid in predicted_ids]

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=structure_files,
            map_table=self.structures_map_csv,
            format=self.output_format
        )

        # Supplied MSAs keep their upstream format; public-server MSAs are csv.
        msa_fmt = self._supplied_msa_format() if self.msas is not None else "csv"
        msa_files = [self.stream_path("msas", f"<id>.{msa_fmt}")]

        msas = DataStream(
            name="msas",
            ids=predicted_ids,
            files=msa_files,
            map_table=self.msas_csv,
            format=msa_fmt
        )

        # Sequences DataStream (from input)
        sequences = DataStream(
            name="sequences",
            ids=predicted_ids,
            files=[],  # Sequences are value-based, stored in map_table
            map_table=self.sequences_csv,
            format="csv"
        )

        # Compounds DataStream — emit a FRESH value-based csv carrying the same
        # ids/SMILES as the input ligands but with `code` overwritten with the
        # residue code Boltz assigns (written at runtime by pipe_boltz_compounds).
        if self.ligands_stream:
            compounds = DataStream(
                name="compounds",
                ids=self.ligands_stream.ids,
                files=[],
                map_table=self.compounds_map_csv,
                format="csv",
            )
        else:
            compounds = None

        # Tables
        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_map_csv,
                columns=["id", "file"],
                description="Boltz2 predicted structures"
            ),
            "confidence": TableInfo(
                name="confidence",
                path=self.confidence_csv,
                columns=["id", "input_file", "confidence_score", "ptm", "iptm", "complex_plddt", "complex_iplddt"],
                description="Boltz2 confidence scores"
            ),
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence"],
                description="Input protein sequences"
            ),
            "msas": TableInfo(
                name="msas",
                path=self.msas_csv,
                columns=["id", "sequences.id", "sequence", "msa_file"],
                description="MSA files for recycling"
            )
        }

        if self.affinity:
            tables["affinity"] = TableInfo(
                name="affinity",
                path=self.affinity_csv,
                columns=["id", "input_file", "affinity_pred_value", "affinity_probability_binary"],
                description="Boltz2 affinity predictions"
            )

        if self.ligands_stream:
            tables["compounds"] = TableInfo(
                name="compounds",
                path=self.compounds_map_csv,
                columns=["id", "format", "code", "smiles", "ccd"],
                description="Ligand compounds with the residue code Boltz assigned"
            )

        # Declare a `missing` table whenever any input axis carries one, so the
        # completion check excuses upstream-filtered ids instead of failing.
        if self._collect_upstream_missing_paths(*self._missing_input_sources()):
            tables["missing"] = self.missing_table_info(self.missing_csv)

        return {
            "structures": structures,
            "sequences": sequences,
            "compounds": compounds,
            "msas": msas,
            "tables": tables,
            "output_folder": self.output_folder,
            "rendering_parameters": {
                "structures": {
                    "color_by": "plddt",
                    "plddt_upper": 100,
                }
            }
        }

    def get_config_display(self) -> List[str]:
        """Get configuration display lines for pipeline output."""
        config_lines = super().get_config_display()

        if self.ssDNA_stream:
            config_lines.append(f"ssDNA: {len(self.ssDNA_stream)} sequences")
        if self.dsDNA_stream:
            config_lines.append(f"dsDNA: {len(self.dsDNA_stream)} sequences")
        if self.ssRNA_stream:
            config_lines.append(f"ssRNA: {len(self.ssRNA_stream)} sequences")
        if self.dsRNA_stream:
            config_lines.append(f"dsRNA: {len(self.dsRNA_stream)} sequences")
        if self.ligands_stream:
            config_lines.append(f"Ligands: {len(self.ligands_stream)} compounds")

        config_lines.extend([
            f"Output format: {self.output_format}",
            f"MSA source: {'supplied' if self.msas is not None else 'public server'}",
            f"Affinity calculation: {self.affinity}"
        ])

        if self.recycling_steps is not None:
            config_lines.append(f"Recycling steps: {self.recycling_steps}")

        if self.diffusion_samples is not None:
            config_lines.append(f"Diffusion samples: {self.diffusion_samples}")

        if self.use_potentials:
            config_lines.append(f"Use potentials: {self.use_potentials}")

        if self.template:
            config_lines.append(f"Template: {os.path.basename(self.template)}")
            if self.template_chain_ids:
                config_lines.append(f"Template chains: {', '.join(self.template_chain_ids)}")

        if self.pocket_residues:
            config_lines.append(f"Pocket residues: {self.pocket_residues}")

        if self.glycosylation:
            config_lines.append(f"Glycosylation: {self.glycosylation}")

        if self.covalent_linkage:
            config_lines.append(f"Covalent linkage: {self.covalent_linkage}")

        if self.contacts:
            config_lines.append(f"Contact constraints: {len(self.contacts)}")

        if self.disulfide_bonds:
            config_lines.append(f"Disulfide bonds: {len(self.disulfide_bonds)}")

        if self.metal_coord:
            config_lines.append(f"Metal coordination bonds: {len(self.metal_coord)}")

        return config_lines

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "boltz2_params": {
                "config": self.config,
                "proteins": str(self.proteins) if self.proteins else None,
                "ssDNA": str(self.ssDNA) if self.ssDNA else None,
                "dsDNA": str(self.dsDNA) if self.dsDNA else None,
                "ssRNA": str(self.ssRNA) if self.ssRNA else None,
                "dsRNA": str(self.dsRNA) if self.dsRNA else None,
                "ligands": str(self.ligands) if self.ligands else None,
                "msas": str(self.msas) if self.msas else None,
                "affinity": self.affinity,
                "output_format": self.output_format,
                "recycling_steps": self.recycling_steps,
                "diffusion_samples": self.diffusion_samples,
                "use_potentials": self.use_potentials,
                "template": self.template,
                "template_chain_ids": self.template_chain_ids,
                "template_force": self.template_force,
                "template_threshold": self.template_threshold,
                "pocket_residues": self.pocket_residues,
                "pocket_max_distance": self.pocket_max_distance,
                "pocket_force": self.pocket_force,
                "glycosylation": self.glycosylation,
                "covalent_linkage": self.covalent_linkage,
                "contacts": self.contacts,
                "disulfide_bonds": self.disulfide_bonds,
                "metal_coord": self.metal_coord,
            }
        })
        return base_dict
