# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Mutagenesis (Mutagenesis) tool for generating amino acid substitutions.

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


class Mutagenesis(BaseConfig):
    """
    Mutagenesis tool for systematic amino acid substitutions.

    Generates amino acid substitutions at specified positions using various
    class-based strategies including saturation mutagenesis and targeted
    amino acid class substitutions.
    """

    TOOL_NAME = "Mutagenesis"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Mutagenesis (Mutagenesis) ==="
echo "Requires MutationEnv (installed with MutationProfiler.install())"
echo "No additional installation needed."
echo "=== Mutagenesis ready ==="
"""

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
    mutagenesis_helper_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_mutagenesis.py"))

    def __init__(self,
                 original: Union[str, DataStream, StandardizedOutput],
                 position: int,
                 mutate_to: str = "",
                 mode: str = "specific",
                 include_original: bool = False,
                 exclude: str = "",
                 prefix: str = "",
                 **kwargs):
        """
        Initialize Mutagenesis tool.

        Args:
            original: Input sequence - can be:
                - DataStream or StandardizedOutput containing sequences
                - Direct protein sequence string
            position: Target position for mutagenesis (1-indexed)
            mutate_to: Target amino acid(s) for "specific" mode, as single letter
                codes (e.g., "A" for alanine, "AV" for alanine and valine).
                Required when mode is "specific".
            mode: Mutagenesis strategy - one of:
                - "specific": Only the amino acid(s) given in `mutate_to` (default)
                - "saturation": All 20 amino acids
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
            # Convert position 50 to alanine
            sdm = Mutagenesis(original=structure, position=50, mutate_to="A")

            # Try alanine and valine at position 50
            sdm = Mutagenesis(original=structure, position=50, mutate_to="AV")

            # Saturation mutagenesis
            sdm = Mutagenesis(original=boltz_output, position=167, mode="saturation")

            # Charged amino acids excluding histidine
            sdm = Mutagenesis(original=structure, position=175,
                              mode="charged", exclude="H")
        """
        # Store Mutagenesis-specific parameters
        self.position = position
        self.mutate_to = mutate_to.upper()
        self.mode = mode
        self.include_original = include_original
        self.exclude = exclude.upper()
        self.prefix = prefix

        # Handle different input types
        self.input_sequence = None
        self.input_sequence_id = None
        self.sequences_stream = None

        if isinstance(original, StandardizedOutput):
            if original.streams.sequences and len(original.streams.sequences) > 0:
                self.sequences_stream = original.streams.sequences
                self.input_sequence_id = original.streams.sequences.ids[0] if original.streams.sequences.ids else "sequence"
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
        """Validate Mutagenesis-specific parameters."""
        valid_aas = set("ACDEFGHIKLMNPQRSTVWY")

        # Validate mode
        valid_modes = ["specific"] + list(self.AMINO_ACID_CLASSES.keys())
        if self.mode not in valid_modes:
            raise ValueError(f"Invalid mode '{self.mode}'. Valid modes: {valid_modes}")

        # Validate specific mode requires mutate_to
        if self.mode == "specific":
            if not self.mutate_to:
                raise ValueError("mutate_to is required when mode is 'specific'. "
                                 "Provide target amino acid(s) (e.g., mutate_to='A') "
                                 "or use a different mode (e.g., mode='saturation').")
            invalid_aas_target = set(self.mutate_to) - valid_aas
            if invalid_aas_target:
                raise ValueError(f"Invalid amino acids in mutate_to: {invalid_aas_target}")

        # Validate position
        if self.position <= 0:
            raise ValueError("Position must be positive (1-indexed)")

        # Validate exclude string
        if self.exclude:
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
        """Generate Mutagenesis execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# Mutagenesis execution script\n"
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
        """Generate the Mutagenesis execution part of the script."""
        mutate_to_line = f'    --mutate-to "{self.mutate_to}" \\\n' if self.mutate_to else ""
        return f"""echo "Running Mutagenesis"
echo "Position: {self.position}"
echo "Mode: {self.mode}"
{"echo " + '"' + "Mutate to: " + self.mutate_to + '"' if self.mutate_to else ""}
echo "Include original: {self.include_original}"
echo "Exclude: {self.exclude}"

# Run Mutagenesis generation
python {self.mutagenesis_helper_py} \\
    --sequence-source {sequence_source} \\
    --sequence {sequence_param} \\
    --sequence-id {sequence_id_param} \\
    --position {self.position} \\
    --mode {self.mode} \\
{mutate_to_line}    --include-original {str(self.include_original).lower()} \\
    --exclude "{self.exclude}" \\
    --output "{self.sequences_csv}"

echo "Mutagenesis completed successfully"
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
        """Get expected output files after Mutagenesis execution."""
        # Calculate expected number of mutants
        if self.mode == "specific":
            amino_acids = self.mutate_to
        else:
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
                description=f"mutants at position {self.position} using {self.mode} mode",
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

        config_lines.append(f"POSITION: {self.position}")
        config_lines.append(f"MODE: {self.mode}")
        if self.mutate_to:
            config_lines.append(f"MUTATE TO: {self.mutate_to}")
        config_lines.append(f"INCLUDE ORIGINAL: {self.include_original}")
        config_lines.append(f"EXCLUDE: {self.exclude if self.exclude else 'None'}")

        return config_lines

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including Mutagenesis-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "sdm_params": {
                "position": self.position,
                "mutate_to": self.mutate_to,
                "mode": self.mode,
                "include_original": self.include_original,
                "exclude": self.exclude,
                "prefix": self.prefix,
                "input_sequence": self.input_sequence,
                "input_sequence_id": self.input_sequence_id
            }
        })
        return base_dict
