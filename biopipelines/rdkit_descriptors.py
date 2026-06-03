# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""RDKit tool: compute cheminformatics descriptors for a compounds stream.

Reads SMILES from the compounds map_table and writes a wide table of per-compound
RDKit descriptors (MW, logP, TPSA, HBA, HBD, rotatable bonds, QED, fraction
sp3, ...). Runs in the `biopipelines` env — rdkit is already pinned there.
"""

import json
import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream


# Default descriptor set (RDKit Descriptors function names). Users can pass a
# subset or extend via the `descriptors` param. Names match `Descriptors.<name>`.
DEFAULT_DESCRIPTORS = [
    "MolWt", "MolLogP", "TPSA", "NumHAcceptors", "NumHDonors",
    "NumRotatableBonds", "NumAromaticRings", "FractionCSP3", "HeavyAtomCount",
    "RingCount", "qed",
]


class RDKit(BaseConfig):
    """
    RDKit: compute cheminformatics descriptors for a compounds stream.

    Inputs:
        compounds: SMILES read from the compounds map_table (`smiles` column).
        descriptors: list of RDKit descriptor function names to compute
            (default: MolWt, MolLogP, TPSA, NumHAcceptors, NumHDonors,
            NumRotatableBonds, NumAromaticRings, FractionCSP3, HeavyAtomCount,
            RingCount, qed).
        morgan_fp: if True, append a Morgan fingerprint (radius=2, n_bits=1024)
            as columns fp_0 .. fp_1023.

    Outputs:
        Streams: (none)
        Tables:
            descriptors: id | smiles | <descriptor columns> [| fp_0 .. fp_1023]
    """

    TOOL_NAME = "RDKit"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== RDKit ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== RDKit ready ==="
"""

    descriptors_csv = Path(lambda self: self.table_path("descriptors"))
    config_json = Path(lambda self: self.configuration_path("rdkit_config.json"))
    compounds_json = Path(lambda self: self.configuration_path("compounds_ds.json"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_rdkit_descriptors.py"))

    def __init__(self,
                 compounds: Union[DataStream, StandardizedOutput],
                 descriptors: List[str] = None,
                 morgan_fp: bool = False,
                 **kwargs):
        if isinstance(compounds, StandardizedOutput):
            self.compounds_stream: DataStream = compounds.streams.compounds
        elif isinstance(compounds, DataStream):
            self.compounds_stream = compounds
        else:
            raise ValueError(f"compounds must be DataStream or StandardizedOutput, got {type(compounds)}")

        self.descriptors = list(descriptors) if descriptors else list(DEFAULT_DESCRIPTORS)
        self.morgan_fp = bool(morgan_fp)
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.compounds_stream or len(self.compounds_stream) == 0:
            raise ValueError("compounds parameter is required and must not be empty")
        if not self.descriptors:
            raise ValueError("descriptors must be a non-empty list")
        for i, d in enumerate(self.descriptors):
            if not isinstance(d, str):
                raise ValueError(f"descriptors[{i}] must be a string (RDKit function name)")
            _validate_freeform_string(f"descriptors[{i}]", d)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"COMPOUNDS: {len(self.compounds_stream)} molecules")
        lines.append(f"DESCRIPTORS: {len(self.descriptors)} ({', '.join(self.descriptors[:6])}{'...' if len(self.descriptors) > 6 else ''})")
        if self.morgan_fp:
            lines.append("MORGAN FP: 1024 bits, radius 2")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.compounds_stream.save_json(self.compounds_json)
        os.makedirs(os.path.dirname(self.config_json), exist_ok=True)
        with open(self.config_json, "w", encoding="utf-8") as f:
            json.dump({
                "compounds_json": self.compounds_json,
                "output_csv": self.descriptors_csv,
                "descriptors": self.descriptors,
                "morgan_fp": self.morgan_fp,
            }, f, indent=2)

        script = "#!/bin/bash\n"
        script += "# RDKit descriptors script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Computing RDKit descriptors for {len(self.compounds_stream)} compound(s)"
python "{self.helper_py}" --config "{self.config_json}"
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        columns = ["id", "smiles"] + list(self.descriptors)
        if self.morgan_fp:
            columns += [f"fp_{i}" for i in range(1024)]
        tables = {
            "descriptors": TableInfo(
                name="descriptors",
                path=self.descriptors_csv,
                columns=columns,
                description="RDKit physicochemical descriptors per compound",
            )
        }
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
