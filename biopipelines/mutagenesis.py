# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
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
    from .combinatorics import generate_multiplied_ids
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_multiplied_ids


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
    missing_csv = Path(lambda self: os.path.join(self.output_folder, "missing.csv"))
    mutagenesis_helper_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_mutagenesis.py"))

    def __init__(self,
                 original: Union[DataStream, StandardizedOutput],
                 position: int,
                 mutate_to: str = "",
                 mode: str = "specific",
                 include_original: bool = False,
                 exclude: str = "",
                 **kwargs):
        """
        Initialize Mutagenesis tool.

        Args:
            original: Input sequence as DataStream or StandardizedOutput
                      (e.g., from Sequence() entity or upstream tool)
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
            **kwargs: Additional parameters

        Examples:
            # Convert position 50 to alanine
            seq = Sequence("MKTVRQ...", ids="my_protein")
            sdm = Mutagenesis(original=seq, position=50, mutate_to="A")

            # Saturation mutagenesis
            sdm = Mutagenesis(original=seq, position=167, mode="saturation")

            # Charged amino acids excluding histidine
            sdm = Mutagenesis(original=seq, position=175,
                              mode="charged", exclude="H")
        """
        # Store Mutagenesis-specific parameters
        self.position = position
        self.mutate_to = mutate_to.upper()
        self.mode = mode
        self.include_original = include_original
        self.exclude = exclude.upper()

        # Handle input types
        self.sequences_stream = None

        if isinstance(original, StandardizedOutput):
            if original.streams.sequences and len(original.streams.sequences) > 0:
                self.sequences_stream = original.streams.sequences
            else:
                raise ValueError("StandardizedOutput has no sequences")
        elif isinstance(original, DataStream):
            if len(original) == 0:
                raise ValueError("DataStream is empty")
            self.sequences_stream = original
        else:
            raise ValueError(f"original must be DataStream or StandardizedOutput, got {type(original)}")

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

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and dependencies."""
        self.folders = pipeline_folders

    def generate_script(self, script_path: str) -> str:
        """Generate Mutagenesis execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# Mutagenesis execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        script_content += self._generate_script_run_sdm()
        script_content += self._generate_script_missing_sequences()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_sdm(self) -> str:
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
    --sequences "{self.sequences_stream.map_table}" \\
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
        step_tool_name = os.path.basename(self.output_folder)
        return f"""# Generate missing.csv if original is excluded
if [ "{str(self.include_original).lower()}" = "false" ]; then
    echo "Creating missing table..."
    python -c "
import pandas as pd

try:
    sequences_df = pd.read_csv('{self.sequences_csv}')
    if len(sequences_df) > 0:
        first_row = sequences_df.iloc[0]
        base_id = first_row['id'].rsplit('_', 1)[0]
        original_aa = first_row['original_aa']

        original_id = f'{{base_id}}_{self.position}{{original_aa}}'
        missing_data = [{{
            'id': original_id,
            'removed_by': '{step_tool_name}',
            'cause': f'Original amino acid ({{original_aa}}) excluded from mutagenesis at position {self.position}'
        }}]

        missing_df = pd.DataFrame(missing_data)
        missing_df.to_csv('{self.missing_csv}', index=False)
        print(f'Created missing.csv: {self.missing_csv}')
    else:
        print('Warning: No sequences found in main CSV')
        pd.DataFrame(columns=['id', 'removed_by', 'cause']).to_csv('{self.missing_csv}', index=False)
except Exception as e:
    print(f'Error creating missing.csv: {{e}}')
"
else
    echo "id,removed_by,cause" > "{self.missing_csv}"
fi

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after Mutagenesis execution."""
        # Calculate expected number of mutants (approximate — original_aa unknown at configuration time)
        if self.mode == "specific":
            amino_acids = self.mutate_to
        else:
            amino_acids = self.AMINO_ACID_CLASSES[self.mode]

        # Remove excluded amino acids
        if self.exclude:
            amino_acids = ''.join(aa for aa in amino_acids if aa not in self.exclude)

        # Multiply by number of input sequences
        num_input = len(self.sequences_stream.ids)
        num_mutants_per_input = len(amino_acids)
        num_mutants = num_mutants_per_input * num_input

        # Generate sequence IDs (approximate — original_aa resolved at runtime)
        suffixes = [f"{self.position}{aa}" for aa in amino_acids]
        sequence_ids, provenance = generate_multiplied_ids(
            self.sequences_stream.ids, suffixes,
            input_stream_name="original"
        )

        # Prepare tables
        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "original.id", "sequence", "mutations", "mutation_positions", "original_aa", "new_aa"],
                description=f"mutants at position {self.position} using {self.mode} mode",
                count=num_mutants
            )
        }

        # Add missing table if original is excluded
        if not self.include_original:
            tables["missing"] = TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "cause"],
                description="IDs removed (original amino acid excluded from mutagenesis)",
                count=1
            )

        # Create sequences DataStream (id-value stream, files empty)
        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[],
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
                "exclude": self.exclude
            }
        })
        return base_dict
