"""
DNAEncoder tool for reverse translation of protein sequences to DNA with codon optimization.

Takes protein sequences as input and generates DNA sequences optimized for specific organisms
using organism-specific codon usage tables (CoCoPUTs). Outputs both CSV tables and Excel files
with color-coded codon frequencies.
"""

import os
from typing import Dict, List, Any, Union

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

    # Lazy path descriptors
    dna_csv = Path(lambda self: os.path.join(self.output_folder, "dna.csv"))
    dna_excel = Path(lambda self: os.path.join(self.output_folder, "dna_sequences.xlsx"))
    info_txt = Path(lambda self: os.path.join(self.output_folder, "dna_info.txt"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "dna_encoder_config.json"))
    encoder_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_dna_encoder.py"))

    def __init__(self,
                 sequences: Union[DataStream, StandardizedOutput],
                 organism: str = "EC",
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
            **kwargs: Additional parameters

        Examples:
            # Encode for E. coli
            dna = DNAEncoder(sequences=lmpnn, organism="EC")

            # Encode for both E. coli and human (conservative approach)
            dna = DNAEncoder(sequences=designed_sequences, organism="EC&HS")
        """
        # Resolve input to DataStream
        if isinstance(sequences, StandardizedOutput):
            self.sequences_stream: DataStream = sequences.streams.sequences
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
        else:
            raise ValueError(f"sequences must be DataStream or StandardizedOutput, got {type(sequences)}")

        self.organism = organism

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate DNAEncoder parameters."""
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")

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
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# DNAEncoder execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_encoding()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_encoding(self) -> str:
        """Generate the DNA encoding part of the script."""
        import json

        config_data = {
            "sequences_csv": self.sequences_stream.map_table,
            "organism": self.organism,
            "dna_output": self.dna_csv,
            "excel_output": self.dna_excel,
            "info_output": self.info_txt
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Encoding protein sequences to DNA"
echo "Input sequences: {self.sequences_stream.map_table}"
echo "Target organism(s): {self.organism}"
echo "Output folder: {self.output_folder}"

python "{self.encoder_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after DNA encoding."""
        # DNA sequences inherit IDs from input sequences
        sequence_ids = self.sequences_stream.ids.copy()

        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[],
            map_table=self.dna_csv,
            format="dna"
        )

        tables = {
            "dna": TableInfo(
                name="dna",
                path=self.dna_csv,
                columns=["id", "protein_sequence", "dna_sequence", "organism", "method"],
                description="DNA sequences with thresholded weighted codon optimization",
                count=len(sequence_ids)
            )
        }

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": sequences,
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder,
            "excel": self.dna_excel,
            "info": self.info_txt
        }
