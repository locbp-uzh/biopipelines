"""
DNAEncoder tool for reverse translation of protein sequences to DNA with codon optimization.

Takes protein sequences as input and generates DNA sequences optimized for specific organisms
using organism-specific codon usage tables (CoCoPUTs). Outputs both CSV tables and Excel files
with color-coded codon frequencies.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class DNAEncoder(BaseConfig):
    """
    Pipeline tool for reverse translating protein sequences to optimized DNA sequences.

    Uses organism-specific codon usage tables from CoCoPUTs (HIVE, updated April 2024)
    to generate DNA sequences with optimal codon frequencies.

    Supports:
    - Single or multiple organisms (EC: E. coli, SC: S. cerevisiae, HS: H. sapiens)
    - Thresholded weighted sampling (recommended): samples codons ≥10‰ frequency
    - CSV output with DNA sequences
    - Excel output with color-coded codon frequencies

    Citation: Please cite CoCoPUTs (HIVE) when using this tool.
    """

    # Tool identification
    TOOL_NAME = "DNAEncoder"
    

    def __init__(self,
                 sequences: Union[ToolOutput, StandardizedOutput],
                 organism: str = "EC",
                 **kwargs):
        """
        Initialize DNA encoder tool.

        Args:
            sequences: Input protein sequences (ToolOutput or StandardizedOutput with sequences table)
            organism: Target organism(s) for codon optimization. Options:
                     - "EC" (Escherichia coli)
                     - "SC" (Saccharomyces cerevisiae)
                     - "HS" (Homo sapiens)
                     - "EC&HS" (optimized for both E. coli and human)
                     - "EC&SC" (optimized for both E. coli and yeast)
                     - "HS&SC" (optimized for both human and yeast)
                     - "EC&HS&SC" (optimized for all three organisms)
            **kwargs: Additional parameters

        Examples:
            # Encode for E. coli
            dna = pipeline.add(DNAEncoder(
                sequences=lmpnn,
                organism="EC"
            ))

            # Encode for both E. coli and human (conservative approach)
            dna = pipeline.add(DNAEncoder(
                sequences=designed_sequences,
                organism="EC&HS"
            ))
        """
        self.sequences_input = sequences
        self.organism = organism

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        if hasattr(sequences, 'config'):
            self.dependencies.append(sequences.config)

    def validate_params(self):
        """Validate DNAEncoder parameters."""
        if not isinstance(self.sequences_input, (ToolOutput, StandardizedOutput)):
            raise ValueError("sequences must be a ToolOutput or StandardizedOutput object")

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
        """Configure input sequences from previous tools."""
        self.folders = pipeline_folders

        # Extract sequences path
        self.sequences_path = self._extract_sequences_path(self.sequences_input)

    def _extract_sequences_path(self, input_obj: Union[ToolOutput, StandardizedOutput]) -> str:
        """Extract sequences path from ToolOutput or StandardizedOutput."""
        # Try to get sequences table from input
        if hasattr(input_obj, 'tables'):
            tables = input_obj.tables

            # Check for sequences table
            if hasattr(tables, 'sequences'):
                if isinstance(tables.sequences, str):
                    return tables.sequences
                else:
                    if hasattr(tables.sequences, 'path'):
                        return tables.sequences.path
                    else:
                        raise ValueError("Cannot extract path of sequences table")
            elif hasattr(tables, '_tables'):
                # Standard BioPipelines format
                for name, info in tables._tables.items():
                    if 'sequence' in name.lower():
                        return info.path
            elif isinstance(tables, dict):
                # Dict format
                for name, info in tables.items():
                    if 'sequence' in name.lower():
                        if isinstance(info, str):
                            return info
                        elif isinstance(info, dict) and 'path' in info:
                            return info['path']
                        elif hasattr(info, 'path'):
                            return info.path

        # Fallback: try to predict sequences file in output folder
        if hasattr(input_obj, 'output_folder'):
            predicted_files = [
                os.path.join(input_obj.output_folder, 'sequences.csv'),
                os.path.join(input_obj.output_folder, 'queries.csv'),
                os.path.join(input_obj.output_folder, 'table.csv')
            ]
            return predicted_files[0]  # Return first prediction

        raise ValueError("Could not extract sequences path from input")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"SEQUENCES: {type(self.sequences_input).__name__}",
            f"ORGANISM: {self.organism}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate script to perform DNA encoding.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        runtime_folder = os.path.dirname(script_path)
        output_folder = self.output_folder
        os.makedirs(output_folder, exist_ok=True)

        # Output files
        dna_csv = os.path.join(output_folder, "dna.csv")
        dna_excel = os.path.join(output_folder, "dna_sequences.xlsx")
        info_txt = os.path.join(output_folder, "dna_info.txt")

        # Create config file for DNA encoder
        config_file = os.path.join(output_folder, "dna_encoder_config.json")
        config_data = {
            "sequences_csv": self.sequences_path,
            "organism": self.organism,
            "dna_output": dna_csv,
            "excel_output": dna_excel,
            "info_output": info_txt
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# DNAEncoder execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Encoding protein sequences to DNA"
echo "Input sequences: {self.sequences_path}"
echo "Target organism(s): {self.organism}"
echo "Output folder: {output_folder}"

# Run Python DNA encoder script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_dna_encoder.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully encoded sequences to DNA"
    echo "DNA sequences (CSV): {dna_csv}"
    echo "DNA sequences (Excel): {dna_excel}"
    echo "Encoding info: {info_txt}"
else
    echo "Error: Failed to encode sequences to DNA"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after DNA encoding.

        Returns:
            Dictionary with output file paths and table information
        """
        dna_csv = os.path.join(self.output_folder, "dna.csv")
        dna_excel = os.path.join(self.output_folder, "dna_sequences.xlsx")
        info_txt = os.path.join(self.output_folder, "dna_info.txt")

        # Create standardized output
        output = {
            "output_folder": self.output_folder,
            "tables": {
                "dna": TableInfo(
                    name="dna",
                    path=dna_csv,
                    columns=["id", "protein_sequence", "dna_sequence", "organism", "method"],
                    description="DNA sequences with thresholded weighted codon optimization",
                    count=0  # Will be determined at runtime
                )
            },
            "files": {
                "excel": dna_excel,
                "info": info_txt
            }
        }

        return output
