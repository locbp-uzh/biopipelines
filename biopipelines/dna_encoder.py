# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
DNAEncoder tool for reverse translation of protein sequences to DNA with codon optimization.

Takes protein sequences as input and generates DNA sequences optimized for specific organisms
using organism-specific codon usage tables (CoCoPUTs). Outputs both CSV tables and Excel files
with color-coded codon frequencies.
"""

import os
from typing import Dict, List, Any, Union, Optional

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


class DNAEncoder(BaseConfig):
    """
    Pipeline tool for reverse translating protein sequences to optimized DNA sequences.

    Uses organism-specific codon usage tables from CoCoPUTs (HIVE, updated April 2024)
    to generate DNA sequences with optimal codon frequencies.

    Supports:
    - Single or multiple organisms (EC: E. coli, SC: S. cerevisiae, HS: H. sapiens)
    - Thresholded weighted sampling (recommended for synthesis): samples codons ≥10‰ frequency
    - CSV output with DNA sequences
    - Excel output with color-coded codon frequencies

    Citation: Please cite CoCoPUTs (HIVE) when using this tool.
    """

    # Tool identification
    TOOL_NAME = "DNAEncoder"
    TOOL_VERSION = "1.2"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== DNAEncoder ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== DNAEncoder ready ==="
"""

    # Content-bearing: the stream CSV doubles as the map_table.
    sequences_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))
    dna_excel = Path(lambda self: os.path.join(self.extras_folder, "encoded_sequences.xlsx"))
    info_txt = Path(lambda self: os.path.join(self.extras_folder, "dna_info.txt"))
    config_file = Path(lambda self: self.configuration_path("dna_encoder_config.json"))
    sequences_json = Path(lambda self: self.configuration_path(".input_sequences.json"))
    sequences_csv_path = Path(lambda self: self.configuration_path(".input_sequences.csv"))
    encoder_py = Path(lambda self: self.pipe_script_path("pipe_dna_encoder.py"))

    def __init__(self,
                 sequences: Union[DataStream, StandardizedOutput],
                 organism: str = "EC",
                 exclude_sites: Optional[List[str]] = None,
                 **kwargs):
        """
        Initialize DNA encoder tool.

        Args:
            sequences: Input protein sequences as DataStream or StandardizedOutput
            organism: Target organism(s) for codon optimization. Options:
                     - "EC" (Escherichia coli)
                     - "SC" (Saccharomyces cerevisiae)
                     - "HS" (Homo sapiens)
                     - "EC&HS" (optimized for both E. coli and human)
                     - "EC&SC" (optimized for both E. coli and yeast)
                     - "HS&SC" (optimized for both human and yeast)
                     - "EC&HS&SC" (optimized for all three organisms)
                     With more than one organism it is more likely if not inevitable to have rare codons.
            exclude_sites: Sequences to keep out of the coding DNA, e.g. restriction
                     sites ["GAATTC", "GGATCC"]. IUPAC codes are accepted. Each site is
                     excluded on both strands, since a site on the reverse strand cuts
                     too. Codons are re-drawn (with backtracking) until the site is gone,
                     so the result stays within the normal codon-usage thresholds.
            **kwargs: Additional parameters

        Output:
            Streams: sequences (content-bearing; map_table is the table below)
            Tables:
                sequences: id | sequence (DNA) | protein_sequence | organism | method
        """
        # Resolve input to DataStream
        if isinstance(sequences, StandardizedOutput):
            self.sequences_stream: DataStream = sequences.streams.sequences
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
        else:
            raise ValueError(f"sequences must be DataStream or StandardizedOutput, got {type(sequences)}")

        self.organism = organism
        self.exclude_sites = list(exclude_sites) if exclude_sites else []

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate DNAEncoder parameters."""
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")

        valid_iupac = set("ACGTRYSWKMBDHVN")
        for site in self.exclude_sites:
            clean = str(site).strip().upper()
            if not clean:
                raise ValueError("exclude_sites entries must not be empty")
            bad = set(clean) - valid_iupac
            if bad:
                raise ValueError(
                    f"exclude_sites entry '{site}' has non-IUPAC characters: {sorted(bad)}. "
                    f"Give the recognition sequence, e.g. 'GAATTC' for EcoRI.")
            if len(clean) < 4:
                raise ValueError(
                    f"exclude_sites entry '{site}' is only {len(clean)} bases. Sites shorter "
                    f"than 4 occur by chance every few hundred bp and would over-constrain "
                    f"the encoding.")

        # Validate organism parameter
        valid_organisms = ["EC", "SC", "HS"]
        organism_parts = self.organism.split("&")

        for org in organism_parts:
            if org.strip() not in valid_organisms:
                raise ValueError(
                    f"Invalid organism '{org}'. Must be one of: {valid_organisms} "
                    f"or combinations like 'EC&HS', 'EC&SC', 'HS&SC', 'EC&HS&SC'"
                )

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"SEQUENCES: {len(self.sequences_stream)} sequences",
            f"ORGANISM: {self.organism}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate script to perform DNA encoding."""
        self.sequences_stream.save_json(self.sequences_json)

        script_content = "#!/bin/bash\n"
        script_content += "# DNAEncoder execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.generate_filtered_map_table_block(
            self.sequences_json, self.sequences_csv_path, required_columns=["id", "sequence"]
        )
        script_content += self.activate_environment()
        script_content += self.generate_script_run_encoding()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_encoding(self) -> str:
        """Generate the DNA encoding part of the script."""
        import json

        config_data = {
            "sequences_csv": self.sequences_csv_path,
            "organism": self.organism,
            "exclude_sites": self.exclude_sites,
            "dna_output": self.sequences_csv,
            "excel_output": self.dna_excel,
            "info_output": self.info_txt
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Encoding protein sequences to DNA"
echo "Target organism(s): {self.organism}"
echo "Output folder: {self.output_folder}"

python "{self.encoder_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after DNA encoding."""
        # DNA sequences inherit IDs from input sequences
        sequence_ids = list(self.sequences_stream.ids)

        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[],
            map_table=self.sequences_csv,
            format="csv"
        )

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence", "protein_sequence", "organism", "method"],
                description="DNA sequences with thresholded weighted codon optimization"
            )
        }

        return {
            "sequences": sequences,
            "tables": tables,
            "output_folder": self.output_folder,
            "excel": self.dna_excel,
            "info": self.info_txt
        }
