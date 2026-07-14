# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""RDKit tool: cheminformatics descriptors from SMILES and conformer strain from poses.

Reads SMILES from a compounds map_table and writes a wide table of per-compound RDKit descriptors (MW, logP, TPSA, HBA, HBD, rotatable bonds, QED, fraction sp3, ...). With a `structures` input it also scores each posed ligand's internal (torsional) strain against a torsion-restrained reference state. Runs in the `biopipelines` env — rdkit is already pinned there.
"""

import json
import os
from typing import Dict, List, Any, Union, Optional, Tuple

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import TableReference
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import TableReference


# Default descriptor set (RDKit Descriptors function names). Users can pass a
# subset or extend via the `descriptors` param. Names match `Descriptors.<name>`.
DEFAULT_DESCRIPTORS = [
    "MolWt", "MolLogP", "TPSA", "NumHAcceptors", "NumHDonors",
    "NumRotatableBonds", "NumAromaticRings", "FractionCSP3", "HeavyAtomCount",
    "RingCount", "qed",
]

FF_CHOICES = ("auto", "mmff", "uff")

STRAIN_COLUMNS = ["id", "smiles", "e_pose", "e_relaxed", "strain", "ff_engine"]


def _normalize_smiles_template(smiles):
    """Normalize the bond-order template to a literal SMILES or a TableReference."""
    if smiles is None or isinstance(smiles, (str, TableReference)):
        return smiles
    if isinstance(smiles, tuple) and len(smiles) == 2:
        table, column = smiles
        if not isinstance(column, str):
            raise ValueError(
                f"smiles tuple must be (TableInfo|path, column_name), got column {type(column).__name__}"
            )
        if isinstance(table, TableInfo):
            return TableReference(table.info.path, column)
        if isinstance(table, str):
            return TableReference(table, column)
        raise ValueError(
            f"smiles tuple's first element must be a TableInfo or path string, got {type(table).__name__}"
        )
    raise ValueError(
        f"smiles must be a SMILES string, a TableReference (e.g. tool.tables.compounds.smiles), "
        f"or a (TableInfo, column) tuple; got {type(smiles).__name__}"
    )


class RDKit(BaseConfig):
    """
    RDKit: cheminformatics descriptors from SMILES, conformer strain from poses.

    Inputs:
        compounds: SMILES read from the compounds map_table (`smiles` column).
            Produces the `descriptors` table.
        structures: posed ligand coordinate files, one per id (e.g. the structures
            stream of `Ligand(code=..., structures=poses)`). Produces the `strain`
            table.
        descriptors: list of RDKit descriptor function names to compute
            (default: MolWt, MolLogP, TPSA, NumHAcceptors, NumHDonors,
            NumRotatableBonds, NumAromaticRings, FractionCSP3, HeavyAtomCount,
            RingCount, qed).
        morgan_fp: if True, append a Morgan fingerprint (radius=2, n_bits=1024)
            as columns fp_0 .. fp_1023.
        smiles: bond-order template for the posed coordinates — a literal SMILES,
            or a (TableInfo, column) reference for a per-id template. Defaults to
            the `smiles` column of the compounds stream when one is supplied, whose             compound is matched to each pose through the structures map's `compounds.id` provenance (so a `prot+lig` complex resolves to `lig`).
            Bond-order templating is mandatory for the strain path.
        restrain_bonds: explicit torsions to restrain, as (atom_name, atom_name)
            pairs naming the two central atoms of each rotatable bond (PDB atom
            names). Default None restrains every rotatable bond.
        ff: force field — "auto" (MMFF, falling back to UFF where MMFF cannot
            type the molecule), "mmff", or "uff".

    Outputs:
        Streams: (none)
        Tables:
            descriptors: id | smiles | <descriptor columns> [| fp_0 .. fp_1023]
                (only when `compounds` is supplied)
            strain: id | smiles | e_pose | e_relaxed | strain | ff_engine
                (only when `structures` is supplied)
    """

    TOOL_NAME = "RDKit"
    TOOL_VERSION = "1.1"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== RDKit ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== RDKit ready ==="
"""

    descriptors_csv = Path(lambda self: self.table_path("descriptors"))
    strain_csv = Path(lambda self: self.table_path("strain"))
    config_json = Path(lambda self: self.configuration_path("rdkit_config.json"))
    compounds_json = Path(lambda self: self.configuration_path("compounds_ds.json"))
    structures_json = Path(lambda self: self.configuration_path("structures_ds.json"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_rdkit_descriptors.py"))

    def __init__(self,
                 compounds: Union[DataStream, StandardizedOutput] = None,
                 structures: Union[DataStream, StandardizedOutput] = None,
                 descriptors: List[str] = None,
                 morgan_fp: bool = False,
                 smiles: Union[str, Tuple[TableInfo, str], 'TableReference'] = None,
                 restrain_bonds: List[Tuple[str, str]] = None,
                 ff: str = "auto",
                 **kwargs):
        self.compounds_stream: Optional[DataStream] = None
        if compounds is not None:
            if isinstance(compounds, StandardizedOutput):
                self.compounds_stream = compounds.streams.compounds
            elif isinstance(compounds, DataStream):
                self.compounds_stream = compounds
            else:
                raise ValueError(f"compounds must be DataStream or StandardizedOutput, got {type(compounds)}")

        self.structures_stream: Optional[DataStream] = None
        if structures is not None:
            if isinstance(structures, StandardizedOutput):
                self.structures_stream = structures.streams.structures
            elif isinstance(structures, DataStream):
                self.structures_stream = structures
            else:
                raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        self.descriptors = list(descriptors) if descriptors else list(DEFAULT_DESCRIPTORS)
        self.morgan_fp = bool(morgan_fp)
        self.smiles = _normalize_smiles_template(smiles)
        self.restrain_bonds = [tuple(b) for b in restrain_bonds] if restrain_bonds else None
        self.ff = ff
        super().__init__(**kwargs)

    def validate_params(self):
        if self.compounds_stream is None and self.structures_stream is None:
            raise ValueError("RDKit requires compounds= (descriptors) and/or structures= (strain)")
        if self.compounds_stream is not None and len(self.compounds_stream) == 0:
            raise ValueError("compounds must not be empty")
        if self.structures_stream is not None and len(self.structures_stream) == 0:
            raise ValueError("structures must not be empty")

        if not self.descriptors:
            raise ValueError("descriptors must be a non-empty list")
        for i, d in enumerate(self.descriptors):
            if not isinstance(d, str):
                raise ValueError(f"descriptors[{i}] must be a string (RDKit function name)")
            _validate_freeform_string(f"descriptors[{i}]", d)

        if self.ff not in FF_CHOICES:
            raise ValueError(f"ff must be one of {FF_CHOICES}, got {self.ff!r}")

        if self.structures_stream is None:
            return

        # Bond-order templating is mandatory for strain: coordinate-only perception
        # mis-assigns conjugated and charged systems.
        if self.smiles is None and self.compounds_stream is None:
            raise ValueError(
                "structures= requires a bond-order template: pass smiles=<SMILES>, "
                "smiles=(table, 'column'), or a compounds= stream carrying a 'smiles' column"
            )
        if isinstance(self.smiles, str):
            _validate_freeform_string("smiles", self.smiles)

        if self.restrain_bonds is not None:
            if not self.restrain_bonds:
                raise ValueError("restrain_bonds must be a non-empty list, or None to restrain all rotatable bonds")
            for i, bond in enumerate(self.restrain_bonds):
                if len(bond) != 2 or not all(isinstance(a, str) for a in bond):
                    raise ValueError(f"restrain_bonds[{i}] must be a pair of atom-name strings, got {bond!r}")
                for j, atom in enumerate(bond):
                    _validate_freeform_string(f"restrain_bonds[{i}][{j}]", atom)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        if self.compounds_stream is not None:
            lines.append(f"COMPOUNDS: {len(self.compounds_stream)} molecules")
            lines.append(f"DESCRIPTORS: {len(self.descriptors)} ({', '.join(self.descriptors[:6])}{'...' if len(self.descriptors) > 6 else ''})")
            if self.morgan_fp:
                lines.append("MORGAN FP: 1024 bits, radius 2")
        if self.structures_stream is not None:
            lines.append(f"STRAIN: {len(self.structures_stream)} pose(s), ff={self.ff}")
            if self.restrain_bonds:
                lines.append(f"RESTRAINED BONDS: {', '.join('-'.join(b) for b in self.restrain_bonds)}")
        return lines

    def generate_script(self, script_path: str) -> str:
        cfg = {
            "descriptors": self.descriptors,
            "morgan_fp": self.morgan_fp,
        }
        if self.compounds_stream is not None:
            self.compounds_stream.save_json(self.compounds_json)
            cfg["compounds_json"] = self.compounds_json
            cfg["descriptors_csv"] = self.descriptors_csv
        if self.structures_stream is not None:
            self.structures_stream.save_json(self.structures_json)
            cfg["structures_json"] = self.structures_json
            cfg["strain_csv"] = self.strain_csv
            cfg["smiles"] = str(self.smiles) if self.smiles is not None else None
            cfg["restrain_bonds"] = [list(b) for b in self.restrain_bonds] if self.restrain_bonds else None
            cfg["ff"] = self.ff

        os.makedirs(os.path.dirname(self.config_json), exist_ok=True)
        with open(self.config_json, "w", encoding="utf-8") as f:
            json.dump(cfg, f, indent=2)

        what = []
        if self.compounds_stream is not None:
            what.append(f"descriptors for {len(self.compounds_stream)} compound(s)")
        if self.structures_stream is not None:
            what.append(f"strain for {len(self.structures_stream)} pose(s)")

        script = "#!/bin/bash\n"
        script += "# RDKit descriptors script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Computing RDKit {' and '.join(what)}"
python "{self.helper_py}" --config "{self.config_json}"
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        tables = {}
        if self.compounds_stream is not None:
            columns = ["id", "smiles"] + list(self.descriptors)
            if self.morgan_fp:
                columns += [f"fp_{i}" for i in range(1024)]
            tables["descriptors"] = TableInfo(
                name="descriptors",
                path=self.descriptors_csv,
                columns=columns,
                description="RDKit physicochemical descriptors per compound",
            )
        if self.structures_stream is not None:
            tables["strain"] = TableInfo(
                name="strain",
                path=self.strain_csv,
                columns=list(STRAIN_COLUMNS),
                description="Force-field conformer strain per posed ligand (kcal/mol)",
            )
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
