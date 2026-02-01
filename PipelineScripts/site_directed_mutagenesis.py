"""
Site-Directed Mutagenesis (SDM) tool for generating amino acid substitutions.

Creates systematic amino acid substitutions at specified positions using
various class-based strategies for comprehensive mutagenesis studies.
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


class SDM(BaseConfig):
    """
    Site-Directed Mutagenesis tool for systematic amino acid substitutions.

    Generates amino acid substitutions at specified positions using various
    class-based strategies including saturation mutagenesis and targeted
    amino acid class substitutions.
    """

    TOOL_NAME = "SDM"

    AMINO_ACID_CLASSES = {
        "saturation": "ACDEFGHIKLMNPQRSTVWY",  # All 20 amino acids
        "hydrophobic": "AFILMVWY",
        "hydrophilic": "DEHKNQRST",
        "charged": "DEHKR",
        "polar": "CDEHKNQRSTY",
        "nonpolar": "AFGILMPVW",
        "aromatic": "FHWY",
        "aliphatic": "AGILV",
        "positive": "HKR",
        "negative": "DE"
    }

    # Lazy path descriptors
    sequences_csv = Path(lambda self: os.path.join(self.output_folder, "sequences.csv"))
    missing_sequences_csv = Path(lambda self: os.path.join(self.output_folder, "missing_sequences.csv"))
    sdm_helper_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_site_directed_mutagenesis.py"))

    def __init__(self,
                 original: Union[str, DataStream, StandardizedOutput],
                 position: int,
                 mode: str = "saturation",
                 include_original: bool = False,
                 exclude: str = "",
                 prefix: str = "",
                 **kwargs):
        """
        Initialize Site-Directed Mutagenesis tool.

        Args:
            original: Input sequence - can be:
                - DataStream or StandardizedOutput containing sequences
                - Direct protein sequence string
            position: Target position for mutagenesis (1-indexed)
            mode: Mutagenesis strategy - one of:
                - "saturation": All amino acids except original
                - "hydrophobic": A, F, I, L, M, V, W, Y
                - "hydrophilic": D, E, H, K, N, Q, R, S, T
                - "charged": D, E, H, K, R
                - "polar": C, D, E, H, K, N, Q, R, S, T, Y
                - "nonpolar": A, F, G, I, L, M, P, V, W
                - "aromatic": F, H, W, Y
                - "aliphatic": A, G, I, L, V
                - "positive": H, K, R
                - "negative": D, E
            include_original: Whether to include original amino acid (default: False)
            exclude: Amino acids to exclude as single string (e.g., "CEFGL")
            prefix: Prefix for sequence IDs (used only for string sequence input)
            **kwargs: Additional parameters

        Examples:
            # Saturation mutagenesis from tool output
            sdm = pipeline.add(SDM(original=boltz_output, position=167, mode="saturation"))

            # Hydrophobic substitutions from sequence string
            sdm = pipeline.add(SDM(original="MKLLVV...", position=5,
                                  mode="hydrophobic", prefix="test"))

            # Charged amino acids excluding histidine
            sdm = pipeline.add(SDM(original=structure, position=175,
                                  mode="charged", exclude="H"))
        """
        # Store SDM-specific parameters
        self.position = position
        self.mode = mode
        self.include_original = include_original
        self.exclude = exclude.upper()
        self.prefix = prefix

        # Handle different input types
        self.input_sequence = None
        self.input_sequence_id = None
        self.sequences_stream = None

        if isinstance(original, StandardizedOutput):
            if original.sequences and len(original.sequences) > 0:
                self.sequences_stream = original.sequences
                self.input_sequence_id = original.sequences.ids[0] if original.sequences.ids else "sequence"
            else:
                raise ValueError("StandardizedOutput has no sequences")
        elif isinstance(original, DataStream):
            if len(original) == 0:
                raise ValueError("DataStream is empty")
            self.sequences_stream = original
            self.input_sequence_id = original.ids[0] if original.ids else "sequence"
        elif isinstance(original, str):
            # Direct sequence string
            self.input_sequence = original.upper()
            self.input_sequence_id = prefix if prefix else "sequence"
        else:
            raise ValueError(f"original must be DataStream, StandardizedOutput, or string, got {type(original)}")

        # Initialize base class
        super().__init__(**kwargs)

    def validate_params(self):
        """Validate SDM-specific parameters."""
        # Validate mode
        if self.mode not in self.AMINO_ACID_CLASSES:
            valid_modes = list(self.AMINO_ACID_CLASSES.keys())
            raise ValueError(f"Invalid mode '{self.mode}'. Valid modes: {valid_modes}")

        # Validate position
        if self.position <= 0:
            raise ValueError("Position must be positive (1-indexed)")

        # Validate exclude string
        if self.exclude:
            valid_aas = set("ACDEFGHIKLMNPQRSTVWY")
            invalid_aas = set(self.exclude) - valid_aas
            if invalid_aas:
                raise ValueError(f"Invalid amino acids in exclude: {invalid_aas}")

        # Validate sequence length if we have the sequence
        if self.input_sequence and self.position > len(self.input_sequence):
            raise ValueError(f"Position {self.position} exceeds sequence length {len(self.input_sequence)}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and dependencies."""
        self.folders = pipeline_folders

    def generate_script(self, script_path: str) -> str:
        """Generate SDM execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# Site-Directed Mutagenesis execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        # Determine sequence source
        if self.input_sequence:
            # Direct sequence provided
            sequence_param = f'"{self.input_sequence}"'
            sequence_id_param = f'"{self.input_sequence_id}"'
            sequence_source = "direct"
        else:
            # Sequence from DataStream file
            sequence_param = f'"{self.sequences_stream.files[0]}"'
            sequence_id_param = f'"{self.input_sequence_id}"'
            sequence_source = "file"

        script_content += self._generate_script_run_sdm(sequence_source, sequence_param, sequence_id_param)
        script_content += self._generate_script_missing_sequences()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_sdm(self, sequence_source: str, sequence_param: str, sequence_id_param: str) -> str:
        """Generate the SDM execution part of the script."""
        return f"""echo "Running Site-Directed Mutagenesis"
echo "Position: {self.position}"
echo "Mode: {self.mode}"
echo "Include original: {self.include_original}"
echo "Exclude: {self.exclude}"

# Run SDM generation
python {self.sdm_helper_py} \\
    --sequence-source {sequence_source} \\
    --sequence {sequence_param} \\
    --sequence-id {sequence_id_param} \\
    --position {self.position} \\
    --mode {self.mode} \\
    --include-original {str(self.include_original).lower()} \\
    --exclude "{self.exclude}" \\
    --output "{self.sequences_csv}"

echo "SDM completed successfully"
echo "Generated sequences saved to: {self.sequences_csv}"

"""

    def _generate_script_missing_sequences(self) -> str:
        """Generate the missing sequences creation part of the script."""
        return f"""# Generate missing sequences CSV if original is excluded
if [ "{str(self.include_original).lower()}" = "false" ]; then
    echo "Creating missing sequences table..."
    python -c "
import pandas as pd
import sys

# Get original sequence and ID from sequences CSV
try:
    sequences_df = pd.read_csv('{self.sequences_csv}')
    if len(sequences_df) > 0:
        first_row = sequences_df.iloc[0]
        original_sequence = first_row['sequence']
        base_id = first_row['id'].rsplit('_', 1)[0]  # Remove _position_aa suffix
        original_aa = first_row['original_aa']

        # Create missing sequence entry
        original_id = f'{{base_id}}_{self.position}{{original_aa}}'
        missing_data = [{{
            'id': original_id,
            'sequence': original_sequence[:({self.position}-1)] + original_aa + original_sequence[{self.position}:],
            'reason': 'Original amino acid excluded from mutagenesis'
        }}]

        missing_df = pd.DataFrame(missing_data)
        missing_df.to_csv('{self.missing_sequences_csv}', index=False)
        print(f'Created missing sequences CSV: {self.missing_sequences_csv}')
    else:
        print('Warning: No sequences found in main CSV')
except Exception as e:
    print(f'Error creating missing sequences CSV: {{e}}')
"
fi

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after SDM execution."""
        # Calculate expected number of mutants
        amino_acids = self.AMINO_ACID_CLASSES[self.mode]

        # Remove excluded amino acids
        if self.exclude:
            amino_acids = ''.join(aa for aa in amino_acids if aa not in self.exclude)

        # Determine original amino acid and remove it unless include_original=True
        if self.input_sequence:
            original_aa = self.input_sequence[self.position - 1]
            if not self.include_original and original_aa in amino_acids:
                amino_acids = amino_acids.replace(original_aa, '')

        num_mutants = len(amino_acids)

        # Generate sequence IDs
        base_id = self.input_sequence_id
        if self.input_sequence:
            original_aa = self.input_sequence[self.position - 1]
        else:
            original_aa = "X"  # Will be determined at runtime

        sequence_ids = []
        for aa in amino_acids:
            mutant_id = f"{base_id}_{self.position}{aa}"
            sequence_ids.append(mutant_id)

        # Prepare tables
        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence", "mutation", "position", "original_aa", "new_aa"],
                description=f"Site-directed mutants at position {self.position} using {self.mode} mode",
                count=num_mutants
            )
        }

        # Add missing sequences table if original is excluded
        if not self.include_original:
            tables["missing_sequences"] = TableInfo(
                name="missing_sequences",
                path=self.missing_sequences_csv,
                columns=["id", "sequence", "reason"],
                description="Sequences excluded from mutagenesis (original amino acid)",
                count=1
            )

        # Create sequences DataStream
        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[self.sequences_csv],
            map_table=self.sequences_csv,
            format="csv"
        )

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": sequences,
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

    def get_config_display(self) -> List[str]:
        """Get configuration display lines for pipeline output."""
        config_lines = super().get_config_display()

        if self.input_sequence:
            config_lines.append(f"INPUT: Direct sequence ({len(self.input_sequence)} aa)")
        else:
            config_lines.append(f"INPUT: {len(self.sequences_stream)} sequences from DataStream")

        config_lines.extend([
            f"POSITION: {self.position}",
            f"MODE: {self.mode}",
            f"INCLUDE ORIGINAL: {self.include_original}",
            f"EXCLUDE: {self.exclude if self.exclude else 'None'}"
        ])

        return config_lines

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including SDM-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "sdm_params": {
                "position": self.position,
                "mode": self.mode,
                "include_original": self.include_original,
                "exclude": self.exclude,
                "prefix": self.prefix,
                "input_sequence": self.input_sequence,
                "input_sequence_id": self.input_sequence_id
            }
        })
        return base_dict
