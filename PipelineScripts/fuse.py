# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Fuse configuration for protein fusion sequence generation.

Creates fusion sequences by linking multiple proteins with flexible linkers
of varying lengths, generating all possible combinations for downstream
folding and analysis.
"""

import os
import json
from typing import Dict, List, Any, Union
from itertools import product

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


class Fuse(BaseConfig):
    """
    Configuration for protein fusion sequence generation.

    Generates all combinations of protein fusions with variable linker lengths,
    producing meaningful sequence IDs and structured output for downstream
    folding tools like AlphaFold.
    """

    TOOL_NAME = "Fuse"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba"):
        return """echo "=== Fuse ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== Fuse ready ==="
"""

    # Lazy path descriptors
    queries_csv = Path(lambda self: os.path.join(self.output_folder, f"{self._get_job_base()}_queries.csv"))
    queries_fasta = Path(lambda self: os.path.join(self.output_folder, f"{self._get_job_base()}_queries.fasta"))
    fuse_queries_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_fuse_queries.py"))

    def __init__(self,
                 proteins: Union[List[Union[str, DataStream, StandardizedOutput]], DataStream, StandardizedOutput],
                 name: str = "",
                 linker: str = "GGGGSGGGGSGGGGSGGGGS",
                 linker_lengths: List[str] = None,
                 **kwargs):
        """
        Initialize Fuse configuration.

        Args:
            proteins: List of protein sequences/structures/DataStreams/StandardizedOutputs,
                     or a single DataStream/StandardizedOutput containing multiple sequences
            name: Job name for output files
            linker: Linker sequence to use (will be truncated to specified lengths)
            linker_lengths: List of length ranges for each junction (e.g., ["1-6", "1-6"])
            **kwargs: Additional parameters
        """
        # Resolve input proteins to a list of strings (sequences or file paths)
        self.input_proteins = self._resolve_proteins(proteins)

        # Validate we have enough proteins for fusion
        if len(self.input_proteins) < 2:
            raise ValueError("Fuse requires at least 2 proteins for fusion")

        # Store Fuse-specific parameters
        self.name = name
        self.linker = linker

        # Set default linker_lengths based on number of proteins if not provided
        if linker_lengths is None:
            expected_junctions = len(self.input_proteins) - 1
            self.linker_lengths = ["1-6"] * expected_junctions
        else:
            self.linker_lengths = linker_lengths

        # Validate linker_lengths matches protein count
        expected_junctions = len(self.input_proteins) - 1
        if len(self.linker_lengths) != expected_junctions:
            raise ValueError(f"linker_lengths must have {expected_junctions} entries for {len(self.input_proteins)} proteins")

        # Initialize base class
        super().__init__(**kwargs)

    def _resolve_proteins(self, proteins) -> List[str]:
        """
        Resolve protein input to a list of strings (sequences or file paths).

        Args:
            proteins: Various input formats

        Returns:
            List of protein strings (sequences or file paths)
        """
        if isinstance(proteins, StandardizedOutput):
            # Extract sequences from StandardizedOutput
            if proteins.streams.sequences and len(proteins.streams.sequences) > 0:
                return proteins.streams.sequences.files
            elif proteins.streams.structures and len(proteins.streams.structures) > 0:
                return proteins.streams.structures.files
            else:
                raise ValueError("StandardizedOutput has no sequences or structures")

        elif isinstance(proteins, DataStream):
            # Extract files from DataStream
            if len(proteins) == 0:
                raise ValueError("DataStream is empty")
            return proteins.files

        elif isinstance(proteins, list):
            # Process each item in the list
            resolved = []
            for item in proteins:
                if isinstance(item, StandardizedOutput):
                    # Extract first file from StandardizedOutput
                    if item.streams.structures and len(item.streams.structures) > 0:
                        resolved.append(item.streams.structures.files[0])
                    elif item.streams.sequences and len(item.streams.sequences) > 0:
                        resolved.append(item.streams.sequences.files[0])
                    else:
                        raise ValueError("StandardizedOutput has no structures or sequences")
                elif isinstance(item, DataStream):
                    # Extract first file from DataStream
                    if len(item) == 0:
                        raise ValueError("DataStream is empty")
                    resolved.append(item.files[0])
                elif isinstance(item, str):
                    resolved.append(item)
                else:
                    raise ValueError(f"Unsupported protein type in list: {type(item)}")
            return resolved
        else:
            raise ValueError(f"proteins must be List, DataStream, or StandardizedOutput, got {type(proteins)}")

    def _get_job_base(self) -> str:
        """Get job base name for file naming."""
        return self.name if self.name else "fuse"

    def validate_params(self):
        """Validate Fuse-specific parameters."""
        if not self.input_proteins:
            raise ValueError("input_proteins parameter is required")

        if len(self.input_proteins) < 2:
            raise ValueError("Fuse requires at least 2 proteins")

        if not self.linker:
            raise ValueError("linker sequence is required")

        if not self.linker_lengths:
            raise ValueError("linker_lengths is required")

        # Validate linker_lengths format (basic check)
        for length_spec in self.linker_lengths:
            if not isinstance(length_spec, str) or not any(c in length_spec for c in '0123456789'):
                raise ValueError(f"Invalid linker_lengths format: {length_spec}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input proteins."""
        self.folders = pipeline_folders

    def _predict_sequence_ids(self) -> List[str]:
        """
        Predict sequence IDs that will be generated by fusion combinations.

        IDs follow the format: {name}_{length1}_{length2}_{...}
        e.g., sensor_2_3, sensor_2_4, sensor_3_3, etc.

        This matches the output format of pipe_fuse_queries.py.
        """
        # Parse linker lengths to get all combinations
        def parse_length_spec(spec: str) -> List[int]:
            """Parse length specification like '1-6' or '3+5-7'."""
            lengths = []
            if '+' in spec:
                for part in spec.split('+'):
                    lengths.extend(parse_length_spec(part.strip()))
            elif '-' in spec and not spec.startswith('-'):
                if spec.count('-') == 1:
                    start, end = map(int, spec.split('-'))
                    lengths.extend(range(start, end + 1))
                else:
                    # Handle negative values or complex ranges
                    lengths.append(int(spec))
            else:
                lengths.append(int(spec))
            return lengths

        # Get all length combinations
        length_ranges = [parse_length_spec(spec) for spec in self.linker_lengths]

        # Determine base name
        name_base = self.name if self.name else "fused"

        # Generate all combinations using itertools.product
        sequence_ids = []
        for length_combo in product(*length_ranges):
            lengths_str = "_".join(str(l) for l in length_combo)
            seq_id = f"{name_base}_{lengths_str}"
            sequence_ids.append(seq_id)

        return sequence_ids

    def get_config_display(self) -> List[str]:
        """Get Fuse configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"PROTEINS: {len(self.input_proteins)} proteins",
            f"LINKER: {self.linker[:20]}{'...' if len(self.linker) > 20 else ''}",
            f"LINKER_LENGTHS: {', '.join(self.linker_lengths)}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate Fuse execution script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        # Build proteins string for fuse_queries.py
        proteins_str = ";".join(self.input_proteins)

        # Format linker lengths for fuse_queries.py (needs L...L wrapping)
        linker_lengths_str = f"L{';'.join(self.linker_lengths)}L"

        job_base = self._get_job_base()

        script_content = "#!/bin/bash\n"
        script_content += "# Fuse execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""echo "Using direct protein input"
PROTEINS_STR="{proteins_str}"
echo "Proteins: $PROTEINS_STR"

echo "Generating fusion sequence combinations"
echo "Linker: {self.linker}"
echo "Linker lengths: {linker_lengths_str}"

# Call fuse_queries.py to generate all fusion combinations
python {self.fuse_queries_py} "{job_base}" "$PROTEINS_STR" "{self.linker}" "{linker_lengths_str}" "{self.queries_csv}" "{self.queries_fasta}"

# Check if files were created successfully
if [ ! -f "{self.queries_csv}" ]; then
    echo "ERROR: Failed to generate queries CSV file"
    exit 1
fi

if [ ! -f "{self.queries_fasta}" ]; then
    echo "ERROR: Failed to generate queries FASTA file"
    exit 1
fi

echo "Successfully generated fusion sequences:"
echo "CSV file: {self.queries_csv}"
echo "FASTA file: {self.queries_fasta}"

# Show summary
NUM_SEQUENCES=$(tail -n +2 "{self.queries_csv}" | wc -l)
echo "Generated $NUM_SEQUENCES fusion sequence combinations"

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """
        Get expected output files after Fuse execution.

        Returns:
            Dictionary with DataStream objects and tables
        """
        # Predict sequence IDs for downstream tools
        sequence_ids = self._predict_sequence_ids()

        # Build position column names dynamically based on number of proteins
        # For n proteins: D1, L1, D2, L2, ..., D(n-1), L(n-1), Dn
        num_proteins = len(self.input_proteins)
        position_columns = []
        for i in range(1, num_proteins + 1):
            position_columns.append(f"D{i}")
            if i < num_proteins:
                position_columns.append(f"L{i}")

        # Organize tables by content type
        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.queries_csv,
                columns=["id", "sequence", "lengths"] + position_columns,
                description="Fusion sequences with domain/linker positions in PyMOL selection format",
                count=len(sequence_ids)
            )
        }

        # Create sequences DataStream
        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[self.queries_fasta],
            map_table=self.queries_csv,
            format="fasta"
        )

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": sequences,
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including Fuse-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "fuse_params": {
                "name": self.name,
                "linker": self.linker,
                "linker_lengths": self.linker_lengths,
                "num_proteins": len(self.input_proteins)
            }
        })
        return base_dict
