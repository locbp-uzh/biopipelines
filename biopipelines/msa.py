# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
MSA tool for converting Multiple Sequence Alignment files between formats.

Converts MSA files from CSV (Boltz2 public server) to A3M (AlphaFold/ColabFold)
format and vice versa, enabling MSA recycling between prediction tools.
"""

import os
import json
from typing import Dict, List, Any

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream, create_map_table
    from .datastream_resolver import resolve_input_to_datastream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream, create_map_table
    from datastream_resolver import resolve_input_to_datastream


class MSA(BaseConfig):
    """
    Pipeline tool for converting MSA files between formats.

    Converts Multiple Sequence Alignment files between CSV (Boltz2 public server)
    and A3M (AlphaFold/ColabFold) formats, enabling MSA recycling between tools.

    Note: a CSV holds only `key, sequence` — Boltz rejects any other column — so the
    `>` headers carrying the UniRef identifiers ColabFold pairs chains on are kept in
    an `.a3m_meta` sidecar written beside it. convert="a3m" restores them when that
    sidecar is present, and falls back to synthesised indices when it is not, which is
    the case for any CSV produced before the sidecar existed. Such an A3M folds as a
    monomer but cannot pair complex chains, so prefer MMseqs2(output_format="a3m")
    when the alignment is destined for AlphaFold.

    Example:
        af_result = AlphaFold(proteins=seq)
        csv_msas = MSA(af_result, convert="csv")
        boltz = Boltz2(proteins=seq, ligands=lig, msas=csv_msas)
    """

    TOOL_NAME = "MSA"
    TOOL_VERSION = "1.6"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== MSA ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== MSA ready ==="
"""

    # Lazy path descriptors — MSA files live inside the msas/ stream folder
    # alongside their co-located map_table.
    output_msas_folder = Path(lambda self: self.stream_folder("msas"))
    output_msas_csv = Path(lambda self: self.stream_map_path("msas"))
    config_json = Path(lambda self: self.configuration_path("msa_config.json"))
    msa_py = Path(lambda self: self.pipe_script_path("pipe_msa.py"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))

    def __init__(self,
                 msas: StandardizedOutput,
                 convert: str,
                 **kwargs):
        """
        Initialize MSA conversion tool.

        Args:
            msas: Tool output with an msas stream (e.g., from Boltz2, AlphaFold, MMseqs2)
            convert: Target format - "a3m" or "csv"

        Output:
            Streams: msas (.a3m or .csv)
            Tables:
                msas: id | sequence | file
                      One row per input row, under the input row's id; `sequences.id` is carried only when the input had one that its id does not already spell.
                missing: id | removed_by | kind | cause
                      A `failure` row for an id the input table lists whose alignment file is gone or whose `file` cell is empty, plus the upstream manifest forwarded.
        """
        self.msas_input = msas
        self.msas_stream = resolve_input_to_datastream(msas, fallback_stream="msas")
        self.convert = convert

        # Get the MSA table path from upstream
        self.input_msa_table = msas.tables.msas.info.path

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate MSA conversion parameters."""
        if self.convert not in ("a3m", "csv"):
            raise ValueError(f"convert must be 'a3m' or 'csv', got '{self.convert}'")

        if not self.msas_stream or len(self.msas_stream) == 0:
            raise ValueError("Input must have a non-empty msas stream")

        input_format = self.msas_stream.format
        if input_format == self.convert:
            print(f"  Warning: Input MSAs are already in {self.convert} format. "
                  f"Files will be copied without conversion.")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"INPUT FORMAT: {self.msas_stream.format}",
            f"TARGET FORMAT: {self.convert}",
            f"MSA COUNT: {len(self.msas_stream)}"
        ])
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate MSA conversion script."""
        config_data = {
            "input_msa_table": self.input_msa_table,
            "convert": self.convert,
            "output_folder": self.output_msas_folder,
            "output_msas_csv": self.output_msas_csv,
            "local_missing_csv": self.local_missing_csv,
        }

        with open(self.config_json, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# MSA format conversion script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""echo "Converting MSA files to {self.convert} format"
echo "Input table: {self.input_msa_table}"
echo "Output folder: {self.output_msas_folder}"

python "{self.msa_py}" --config "{self.config_json}"

"""
        script_content += self.generate_missing_propagation(
            self.msas_input, missing_csv=self.missing_csv, local_missing=self.local_missing_csv)
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after MSA conversion."""
        # One row per input row, under the input row's id.
        msa_ids = list(self.msas_stream.ids)

        ext = ".a3m" if self.convert == "a3m" else ".csv"
        msa_files = [os.path.join(self.output_msas_folder, f"<id>{ext}")]

        msas = DataStream(
            name="msas",
            ids=msa_ids,
            files=msa_files,
            map_table=self.output_msas_csv,
            format=self.convert
        )

        tables = {
            "msas": TableInfo(
                name="msas",
                path=self.output_msas_csv,
                columns=["id", "sequence", "file"],
                description="Converted MSA files"
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }

        return {
            "msas": msas,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "msa_params": {
                "convert": self.convert,
                "input_format": self.msas_stream.format,
                "input_msa_table": self.input_msa_table,
                "msa_count": len(self.msas_stream)
            }
        })
        return base_dict
