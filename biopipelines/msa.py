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

    Note: AlphaFold → CSV → Boltz2 recycling works. Boltz2 → A3M → AlphaFold
    does not — ColabFold ignores converted A3M files and re-queries the MSA server.

    Example:
        af_result = AlphaFold(proteins=seq)
        csv_msas = MSA(af_result, convert="csv")
        boltz = Boltz2(proteins=seq, ligands=lig, msas=csv_msas)
    """

    TOOL_NAME = "MSA"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== MSA ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== MSA ready ==="
"""

    # Lazy path descriptors
    output_msas_folder = Path(lambda self: os.path.join(self.output_folder, "MSAs"))
    output_msas_csv = Path(lambda self: os.path.join(self.output_folder, "msas.csv"))
    config_json = Path(lambda self: os.path.join(self.output_folder, "msa_config.json"))
    msa_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_msa.py"))

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
                msas: id | sequences.id | sequence | msa_file
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
        os.makedirs(self.output_folder, exist_ok=True)

        config_data = {
            "input_msa_table": self.input_msa_table,
            "convert": self.convert,
            "output_folder": self.output_msas_folder,
            "output_msas_csv": self.output_msas_csv,
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
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after MSA conversion."""
        sequence_ids = list(self.msas_stream.ids)

        ext = ".a3m" if self.convert == "a3m" else ".csv"
        msa_files = [os.path.join(self.output_msas_folder, f"<id>{ext}")]

        msas = DataStream(
            name="msas",
            ids=sequence_ids,
            files=msa_files,
            map_table=self.output_msas_csv,
            format=self.convert
        )

        tables = {
            "msas": TableInfo(
                name="msas",
                path=self.output_msas_csv,
                columns=["id", "sequences.id", "sequence", "msa_file"],
                description="Converted MSA files",
                count=len(sequence_ids)
            )
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
