"""
Site-Directed Mutagenesis (SDM) tool for generating amino acid substitutions.

Creates systematic amino acid substitutions at specified positions using
various class-based strategies for comprehensive mutagenesis studies.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


class SDM(BaseConfig):
    """
    Site-Directed Mutagenesis tool for systematic amino acid substitutions.

    Generates amino acid substitutions at specified positions using various
    class-based strategies including saturation mutagenesis and targeted
    amino acid class substitutions.
    """

    # Tool identification
    TOOL_NAME = "SDM"
    DEFAULT_ENV = "ProteinEnv"
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "4GB", "time": "1:00:00"}

    # Amino acid classifications
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

    def __init__(self,
                 original: Union[str, ToolOutput, StandardizedOutput],
                 position: int,
                 mode: str = "saturation",
                 include_original: bool = False,
                 exclude: str = "",
                 prefix: str = "",
                 **kwargs):
        """
        Initialize Site-Directed Mutagenesis tool.

        Args:
            original: Input structure/sequence - can be:
                - ToolOutput or StandardizedOutput (with sequence_id)
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
        self.original = original
        self.position = position
        self.mode = mode
        self.include_original = include_original
        self.exclude = exclude.upper()  # Ensure uppercase for consistency
        self.prefix = prefix

        # Handle different input types
        self.input_is_tool_output = False
        self.input_sequence = None
        self.input_sequence_id = None

        if isinstance(original, (ToolOutput, StandardizedOutput)):
            self.input_is_tool_output = True
            if isinstance(original, ToolOutput):
                # Add dependency for ToolOutput
                self.dependencies = [original.config]
            # Sequence and ID will be extracted in configure_inputs
        elif isinstance(original, str):
            # Direct sequence string
            self.input_sequence = original.upper()
            self.input_sequence_id = prefix if prefix else "sequence"
        else:
            raise ValueError(f"Invalid original input type: {type(original)}. "
                           "Must be ToolOutput, StandardizedOutput, or string sequence.")

        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths
        self.sequences_csv = None
        self.sdm_helper_py = None

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
        self._setup_file_paths()

        # Extract sequence information from tool outputs
        if self.input_is_tool_output:
            if isinstance(self.original, StandardizedOutput):
                # StandardizedOutput case
                if hasattr(self.original, 'sequences') and self.original.sequences:
                    # We'll read the sequence from the file at runtime
                    self.input_sources = {"sequences": self.original.sequences[0]}
                    if hasattr(self.original, 'sequence_ids') and self.original.sequence_ids:
                        self.input_sequence_id = self.original.sequence_ids[0]
                    else:
                        self.input_sequence_id = "sequence"
                else:
                    raise ValueError("No sequences found in StandardizedOutput")
            else:  # ToolOutput
                # ToolOutput case
                sequences = self.original.get_output_files("sequences")
                if sequences:
                    self.input_sources = {"sequences": sequences[0]}
                    # Try to get sequence IDs
                    sequence_ids = self.original.get_output_files("sequence_ids")
                    if sequence_ids:
                        self.input_sequence_id = sequence_ids[0] if isinstance(sequence_ids[0], str) else "sequence"
                    else:
                        self.input_sequence_id = "sequence"
                else:
                    raise ValueError("No sequences found in ToolOutput")
        else:
            # Direct sequence - no additional configuration needed
            self.input_sources = {}

    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        # Output files
        self.sequences_csv = os.path.join(self.output_folder, "sequences.csv")
        self.missing_sequences_csv = os.path.join(self.output_folder, "missing_sequences.csv")

        # Helper script paths
        if hasattr(self, 'folders') and self.folders:
            self.sdm_helper_py = os.path.join(self.folders["HelpScripts"], "pipe_site_directed_mutagenesis.py")
        else:
            self.sdm_helper_py = None

    def generate_script(self, script_path: str) -> str:
        """Generate SDM execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# Site-Directed Mutagenesis execution script\n"
        script_content += "# Generated by BioPipelines pipeline system\n\n"
        script_content += self.generate_completion_check_header()

        # Determine sequence source
        if self.input_sequence:
            # Direct sequence provided
            sequence_param = f'"{self.input_sequence}"'
            sequence_id_param = f'"{self.input_sequence_id}"'
            sequence_source = "direct"
        else:
            # Sequence from file
            sequence_param = f'"{list(self.input_sources.values())[0]}"'
            sequence_id_param = f'"{self.input_sequence_id}"'
            sequence_source = "file"

        script_content += f"""echo "Running Site-Directed Mutagenesis"
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

# Generate missing sequences CSV if original is excluded
if [ "{str(self.include_original).lower()}" = "false" ]; then
    echo "Creating missing sequences datasheet..."
    python -c "
import pandas as pd
import sys

# Get original sequence and ID from sequences CSV
try:
    sequences_df = pd.read_csv('{self.sequences_csv}')
    if len(sequences_df) > 0:
        first_row = sequences_df.iloc[0]
        original_sequence = first_row['sequence']
        base_id = first_row['id'].rsplit('_', 1)[0]  # Remove _{self.position}{aa} suffix
        original_aa = first_row['original_aa']

        # Create missing sequence entry
        original_id = f"{base_id}_{self.position}{original_aa}"
        missing_data = [{
            'id': original_id,
            'sequence': original_sequence[:(self.position-1)] + original_aa + original_sequence[self.position:],
            'reason': 'Original amino acid excluded from mutagenesis'
        }]

        missing_df = pd.DataFrame(missing_data)
        missing_df.to_csv(self.missing_sequences_csv, index=False)
        print(f'Created missing sequences CSV: {self.missing_sequences_csv}')
    else:
        print('Warning: No sequences found in main CSV')
except Exception as e:
    print(f'Error creating missing sequences CSV: {e}')
"
fi

"""
        script_content += self.generate_completion_check_footer()

        return script_content

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
        if self.input_sequence:
            base_id = self.input_sequence_id
            original_aa = self.input_sequence[self.position - 1]
        else:
            base_id = self.input_sequence_id
            original_aa = "X"  # Will be determined at runtime

        sequence_ids = []
        for aa in amino_acids:
            mutant_id = f"{base_id}_{self.position}{aa}"
            sequence_ids.append(mutant_id)

        # Prepare datasheets
        datasheets = {
            "sequences": DatasheetInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence", "mutation", "position", "original_aa", "new_aa"],
                description=f"Site-directed mutants at position {self.position} using {self.mode} mode",
                count=num_mutants
            )
        }

        # Add missing sequences datasheet if original is excluded
        if not self.include_original:
            datasheets["missing_sequences"] = DatasheetInfo(
                name="missing_sequences",
                path=self.missing_sequences_csv,
                columns=["id", "sequence", "reason"],
                description="Sequences excluded from mutagenesis (original amino acid)",
                count=1
            )

        return {
            "sequences": [self.sequences_csv],
            "sequence_ids": sequence_ids,
            "structures": [],  # SDM doesn't produce structures
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "datasheets": datasheets,
            "output_folder": self.output_folder
        }

    def get_config_display(self) -> List[str]:
        """Get configuration display lines for pipeline output."""
        config_lines = super().get_config_display()

        if self.input_sequence:
            config_lines.append(f"Input: Direct sequence ({len(self.input_sequence)} aa)")
        else:
            config_lines.append(f"Input: {type(self.original).__name__}")

        config_lines.extend([
            f"Position: {self.position}",
            f"Mode: {self.mode}",
            f"Include original: {self.include_original}",
            f"Exclude: {self.exclude if self.exclude else 'None'}"
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