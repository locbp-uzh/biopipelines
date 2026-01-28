"""
GhostFold configuration for database-free protein structure prediction.

GhostFold generates synthetic, structure-aware MSAs from single sequences using
ProstT5, then uses ColabFold for structure prediction. This eliminates the need
for large sequence databases while maintaining prediction accuracy.

Reference: https://github.com/brineylab/ghostfold
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class GhostFold(BaseConfig):
    """
    Configuration for GhostFold structure prediction.

    GhostFold performs database-free protein folding by generating synthetic
    structure-aware MSAs using ProstT5, then running ColabFold for structure
    prediction.

    Modes:
    - Full pipeline (default): Generate MSAs + predict structures
    - MSA-only: Only generate synthetic MSAs for use with other tools
    """

    TOOL_NAME = "GhostFold"
    

    def __init__(self, sequences: Union[str, List[str], ToolOutput, Dict[str, Any]] = None,
                 name: str = "",
                 msa_only: bool = False,
                 num_recycle: int = 10,
                 num_models: int = 5,
                 num_seeds: int = 5,
                 subsample: bool = False,
                 mask_msa: Optional[float] = None,
                 **kwargs):
        """
        Initialize GhostFold configuration.

        Args:
            sequences: Input sequences - can be FASTA file, CSV file, list, ToolOutput, or dict
            name: Job name for output files
            msa_only: If True, only generate synthetic MSAs without structure prediction
            num_recycle: Number of recycling iterations for ColabFold (default 10)
            num_models: Number of AlphaFold2 models to use (default 5)
            num_seeds: Number of random seeds for diversity (default 5)
            subsample: Enable multi-level MSA subsampling (tests max-seq: 16, 32, 64, 128)
            mask_msa: Optional MSA masking fraction (0.0-1.0, e.g., 0.15 for 15%)
            **kwargs: Additional parameters
        """
        # Handle different input formats
        if sequences is not None:
            if isinstance(sequences, StandardizedOutput):
                self.input_sequences = sequences.sequences
                self.input_tables = sequences.tables
                self.input_is_tool_output = False
                self.standardized_input = sequences
            elif isinstance(sequences, ToolOutput):
                self.input_sequences = sequences
                self.input_tables = sequences.get_output_files("tables")
                self.input_is_tool_output = True
                self.standardized_input = None
            elif isinstance(sequences, dict):
                self.input_sequences = sequences.get('sequences', [])
                self.input_tables = sequences.get('tables', {})
                self.input_is_tool_output = False
                self.standardized_input = None
            else:
                self.input_sequences = sequences
                self.input_tables = {}
                self.input_is_tool_output = isinstance(sequences, ToolOutput)
                self.standardized_input = None
        else:
            raise ValueError("sequences parameter is required")

        # Store GhostFold-specific parameters
        self.name = name or kwargs.get('job_name', '')
        self.msa_only = msa_only
        self.num_recycle = num_recycle
        self.num_models = num_models
        self.num_seeds = num_seeds
        self.subsample = subsample
        self.mask_msa = mask_msa

        # Initialize base class
        super().__init__(**kwargs)

        # Initialize file paths
        self._initialize_file_paths()

    def validate_params(self):
        """Validate GhostFold-specific parameters."""
        if not self.input_sequences:
            raise ValueError("input sequences parameter is required")

        if self.num_recycle < 1:
            raise ValueError("num_recycle must be at least 1")

        if self.num_models < 1 or self.num_models > 5:
            raise ValueError("num_models must be between 1 and 5")

        if self.num_seeds < 1:
            raise ValueError("num_seeds must be at least 1")

        if self.mask_msa is not None:
            if self.mask_msa < 0.0 or self.mask_msa > 1.0:
                raise ValueError("mask_msa must be between 0.0 and 1.0")

    def _initialize_file_paths(self):
        """Initialize common file paths used throughout the class."""
        self.queries_fasta = None
        self.queries_csv = None
        self.msa_folder = None
        self.folding_folder = None

        # Helper script paths
        self.ghostfold_sh = None
        self.csv_to_fasta_py = None
        self.ghostfold_confidence_py = None

    def _setup_file_paths(self):
        """Set up all file paths after output_folder is known."""
        job_base = self.name or self._extract_job_name()

        # Core input/output files
        self.queries_fasta = os.path.join(self.output_folder, f"{job_base}_queries.fasta")
        self.queries_csv = os.path.join(self.output_folder, f"{job_base}_queries.csv")
        self.msa_folder = os.path.join(self.output_folder, "MSAs")
        self.folding_folder = os.path.join(self.output_folder, "Folding")

        # Helper script paths
        if hasattr(self, 'folders') and self.folders:
            self.ghostfold_sh = os.path.join(self.folders["GhostFold"], "ghostfold.sh")
            self.csv_to_fasta_py = os.path.join(self.folders["HelpScripts"], "pipe_csv_to_fasta.py")
            self.ghostfold_confidence_py = os.path.join(self.folders["HelpScripts"], "pipe_ghostfold_confidence.py")

    def _extract_job_name(self) -> str:
        """Extract job name from output folder structure."""
        folder_parts = self.output_folder.split(os.sep)
        for i, part in enumerate(folder_parts):
            if "_GhostFold" in part:
                if i == 0:
                    raise ValueError(f"Invalid output folder structure: {self.output_folder}")
                return folder_parts[i-1]

        raise ValueError(f"Could not extract job name from output folder: {self.output_folder}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences from various sources."""
        self.folders = pipeline_folders
        self._setup_file_paths()

        if self.input_is_tool_output:
            tool_output: ToolOutput = self.input_sequences

            source_sequences = []
            for seq_type in ["sequences", "queries_csv", "fa_files"]:
                seq_files = tool_output.get_output_files(seq_type)
                if seq_files:
                    source_sequences = seq_files
                    break

            if not source_sequences:
                raise ValueError(f"No sequence outputs found from {tool_output.tool_type}")

            self.input_sources["sequences"] = source_sequences[0]
            self.dependencies.append(tool_output.config)

        elif isinstance(self.input_sequences, list):
            if self.input_sequences:
                self.input_sources["sequences"] = self.input_sequences[0]
            else:
                raise ValueError("Empty sequence list provided")

        elif isinstance(self.input_sequences, str):
            if self.input_sequences.endswith('.csv'):
                csv_source = os.path.join(pipeline_folders.get("data", os.getcwd()), self.input_sequences)
                if os.path.exists(csv_source):
                    self.input_sources["sequences"] = csv_source
                else:
                    raise ValueError(f"CSV file not found: {csv_source}")
            elif self.input_sequences.endswith('.fasta') or self.input_sequences.endswith('.fa'):
                fasta_source = os.path.join(pipeline_folders.get("data", os.getcwd()), self.input_sequences)
                if os.path.exists(fasta_source):
                    self.input_sources["sequences"] = fasta_source
                else:
                    raise ValueError(f"FASTA file not found: {fasta_source}")
            else:
                if not self.name:
                    raise ValueError("name parameter required when providing direct sequence")
                self.input_sources["direct_sequence"] = self.input_sequences
        else:
            raise ValueError(f"Unsupported input type: {type(self.input_sequences)}")

    def get_config_display(self) -> List[str]:
        """Get GhostFold configuration display lines."""
        config_lines = super().get_config_display()

        if self.input_is_tool_output:
            config_lines.append(f"INPUT: {self.input_sequences.tool_type} output")
        else:
            config_lines.append(f"INPUT: {self.input_sequences}")

        mode = "MSA-only" if self.msa_only else "Full prediction"
        config_lines.append(f"MODE: {mode}")

        if not self.msa_only:
            config_lines.extend([
                f"NUM RECYCLE: {self.num_recycle}",
                f"NUM MODELS: {self.num_models}",
                f"NUM SEEDS: {self.num_seeds}"
            ])

        if self.subsample:
            config_lines.append("SUBSAMPLE: enabled")

        if self.mask_msa is not None:
            config_lines.append(f"MSA MASKING: {self.mask_msa*100:.1f}%")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate GhostFold execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# GhostFold execution script\n"
        script_content += "# Generated by BioPipelines\n\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_prepare_sequences()
        script_content += self.generate_script_run_ghostfold()

        if not self.msa_only:
            script_content += self.generate_script_extract_best_rank()
            script_content += self.generate_script_extract_confidence()

        script_content += self.generate_script_create_msas_table()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_prepare_sequences(self) -> str:
        """Generate the sequence preparation part of the script."""
        if "direct_sequence" in self.input_sources:
            sequence = self.input_sources["direct_sequence"]
            job_base = self.name or "sequence"

            return f"""echo "Creating FASTA from direct sequence"
cat > {self.queries_fasta} << EOF
>{job_base}
{sequence}
EOF

# Also create CSV for tracking
cat > {self.queries_csv} << EOF
id,sequence
{job_base},{sequence}
EOF

"""
        elif "sequences" in self.input_sources:
            source_file = self.input_sources["sequences"]

            if source_file.endswith(".fasta") or source_file.endswith(".fa"):
                return f"""echo "Using FASTA from: {source_file}"
cp "{source_file}" "{self.queries_fasta}"

# Convert FASTA to CSV for tracking
python {self.csv_to_fasta_py} --reverse "{self.queries_fasta}" "{self.queries_csv}"

"""
            else:
                # CSV input - need to convert to FASTA for GhostFold
                return f"""echo "Converting CSV to FASTA: {source_file}"
cp "{source_file}" "{self.queries_csv}"
python {self.csv_to_fasta_py} "{source_file}" "{self.queries_fasta}"

"""
        else:
            raise ValueError("No sequence input configured")

    def generate_script_run_ghostfold(self) -> str:
        """Generate the GhostFold execution part of the script."""
        job_base = self.name or self._extract_job_name()

        # Build GhostFold options
        gf_mode = "--msa-only" if self.msa_only else ""
        gf_options = ""

        if self.subsample:
            gf_options += " --subsample"

        if self.mask_msa is not None:
            gf_options += f" --mask_msa {self.mask_msa}"

        # GhostFold expects the project_name to be the output folder
        # and outputs go into project_name/predictions/<seq_id>/
        return f"""echo "Running GhostFold"
echo "Mode: {'MSA-only' if self.msa_only else 'Full prediction'}"
mkdir -p "{self.msa_folder}"
mkdir -p "{self.folding_folder}"

# Run GhostFold
cd "$(dirname {self.ghostfold_sh})"
bash {self.ghostfold_sh} \\
    --project_name "{self.output_folder}" \\
    --fasta_file "{self.queries_fasta}" \\
    {gf_mode} {gf_options}

cd - > /dev/null

# Move MSA files to MSAs folder
echo "Organizing MSA files"
find "{self.output_folder}" -name "*.a3m" -type f | while read msa_file; do
    if [ "$(dirname "$msa_file")" != "{self.msa_folder}" ]; then
        cp "$msa_file" "{self.msa_folder}/"
    fi
done

"""

    def generate_script_extract_best_rank(self) -> str:
        """Generate script section to extract best structures from ColabFold output."""
        return f"""echo "Extracting best structures"
# GhostFold/ColabFold creates files like:
# - <project>/predictions/<seq_id>/<seq_id>_unrelaxed_rank_001_*.pdb
# We want to copy the best ones to main directory as: <seq_id>.pdb

# Find and copy rank_001 structures
find "{self.output_folder}/predictions" -name "*_rank_001_*.pdb" -type f 2>/dev/null | while read pdb_file; do
    basename_file=$(basename "$pdb_file")

    # Extract sequence ID (everything before _relaxed_rank or _unrelaxed_rank)
    if [[ "$basename_file" == *"_relaxed_rank_"* ]]; then
        seq_id=$(echo "$basename_file" | sed 's/_relaxed_rank_001_.*//')
        # Prefer relaxed over unrelaxed
        cp "$pdb_file" "{self.output_folder}/$seq_id.pdb"
        echo "Extracted: $basename_file -> $seq_id.pdb"
    elif [[ "$basename_file" == *"_unrelaxed_rank_"* ]]; then
        seq_id=$(echo "$basename_file" | sed 's/_unrelaxed_rank_001_.*//')
        # Only copy unrelaxed if relaxed doesn't exist
        if [ ! -f "{self.output_folder}/$seq_id.pdb" ]; then
            cp "$pdb_file" "{self.output_folder}/$seq_id.pdb"
            echo "Extracted: $basename_file -> $seq_id.pdb"
        fi
    fi
done

# Also check Folding folder (if GhostFold puts them there)
shopt -s nullglob
for file in "{self.folding_folder}"/*_rank_001_*.pdb; do
    if [ -f "$file" ]; then
        basename_file=$(basename "$file")
        if [[ "$basename_file" == *"_relaxed_rank_"* ]]; then
            base=$(echo "$basename_file" | sed 's/_relaxed_rank_001_.*/.pdb/')
            cp "$file" "{self.output_folder}/$base"
        elif [[ "$basename_file" == *"_unrelaxed_rank_"* ]]; then
            base=$(echo "$basename_file" | sed 's/_unrelaxed_rank_001_.*/.pdb/')
            if [ ! -f "{self.output_folder}/$base" ]; then
                cp "$file" "{self.output_folder}/$base"
            fi
        fi
    fi
done

"""

    def generate_script_extract_confidence(self) -> str:
        """Generate script section to extract confidence metrics."""
        confidence_csv = os.path.join(self.output_folder, "confidence.csv")

        return f"""echo "Extracting confidence metrics"
python {self.ghostfold_confidence_py} "{self.output_folder}" "{confidence_csv}"

"""

    def generate_script_create_msas_table(self) -> str:
        """Generate script section to create MSAs CSV table."""
        msa_csv = os.path.join(self.output_folder, "msas.csv")

        return f"""echo "Creating MSAs table"
# Create MSAs CSV with id, sequence_id, sequence, msa_file columns
echo "id,sequence_id,sequence,msa_file" > "{msa_csv}"

# Read sequences from queries CSV and match with MSA files
if [ -f "{self.queries_csv}" ]; then
    tail -n +2 "{self.queries_csv}" | while IFS=, read -r seq_id sequence rest; do
        msa_file="{self.msa_folder}/$seq_id.a3m"
        if [ -f "$msa_file" ]; then
            echo "$seq_id,$seq_id,$sequence,$msa_file" >> "{msa_csv}"
        fi
    done
fi

echo "MSAs table created at {msa_csv}"

"""

    def _predict_structure_outputs(self) -> List[str]:
        """Predict structure files that GhostFold will generate."""
        structure_files = []
        sequence_ids = self._get_sequence_ids()

        if sequence_ids:
            for seq_id in sequence_ids:
                pdb_path = os.path.join(self.output_folder, f"{seq_id}.pdb")
                structure_files.append(pdb_path)

        if not structure_files and not self.msa_only:
            raise ValueError("Could not determine sequence IDs for structure prediction")

        return structure_files

    def _predict_msa_outputs(self) -> List[str]:
        """Predict MSA files that GhostFold will generate."""
        msa_files = []
        sequence_ids = self._get_sequence_ids()

        if sequence_ids:
            for seq_id in sequence_ids:
                msa_path = os.path.join(self.msa_folder, f"{seq_id}.a3m")
                msa_files.append(msa_path)

        return msa_files

    def _get_sequence_ids(self) -> List[str]:
        """Get sequence IDs from input sources."""
        sequence_ids = []

        # Try standardized input first
        if hasattr(self, 'standardized_input') and self.standardized_input:
            if hasattr(self.standardized_input, 'sequence_ids'):
                sequence_ids = self.standardized_input.sequence_ids
            elif '_data' in dir(self.standardized_input) and 'sequence_ids' in self.standardized_input._data:
                sequence_ids = self.standardized_input._data['sequence_ids']

        # Try ToolOutput
        if not sequence_ids and self.input_is_tool_output:
            tool_output: ToolOutput = self.input_sequences
            if hasattr(tool_output.config, '_predict_sequence_ids'):
                sequence_ids = tool_output.config._predict_sequence_ids()

        # Try dependencies
        if not sequence_ids and self.dependencies:
            for dep in self.dependencies:
                if hasattr(dep, '_predict_sequence_ids'):
                    sequence_ids = dep._predict_sequence_ids()
                    break

        # Handle direct sequence
        if not sequence_ids and "direct_sequence" in self.input_sources and self.name:
            sequence_ids = [self.name]

        return sequence_ids

    def get_output_files(self) -> Dict[str, List[str]]:
        """Get expected output files after GhostFold execution."""
        if not hasattr(self, 'queries_csv') or self.queries_csv is None:
            self._setup_file_paths()

        sequence_ids = self._get_sequence_ids()
        msa_files = self._predict_msa_outputs()
        msa_csv = os.path.join(self.output_folder, "msas.csv")

        # Build tables dict
        tables = {
            "msas": TableInfo(
                name="msas",
                path=msa_csv,
                columns=["id", "sequence_id", "sequence", "msa_file"],
                description="Synthetic MSA files generated by GhostFold",
                count=len(msa_files)
            )
        }

        # Base output
        output = {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [self.queries_csv],
            "sequence_ids": sequence_ids,
            "msas": msa_files,
            "msa_ids": sequence_ids,
            "tables": tables,
            "output_folder": self.output_folder,
            "msa_folder": self.msa_folder,
            "queries_fasta": self.queries_fasta,
            "queries_csv": [self.queries_csv]
        }

        # Add structure outputs if not msa_only mode
        if not self.msa_only:
            structure_files = self._predict_structure_outputs()
            structure_ids = [os.path.splitext(os.path.basename(f))[0] for f in structure_files]
            confidence_csv = os.path.join(self.output_folder, "confidence.csv")

            output["structures"] = structure_files
            output["structure_ids"] = structure_ids
            output["pdbs"] = structure_files  # Legacy alias

            # Add confidence table
            tables["confidence"] = TableInfo(
                name="confidence",
                path=confidence_csv,
                columns=["id", "structure", "plddt", "max_pae", "ptm"],
                description="GhostFold/ColabFold confidence metrics",
                count=len(structure_files)
            )

            # Add structures table
            tables["structures"] = TableInfo(
                name="structures",
                path=self.queries_csv,
                columns=["id", "sequence"],
                description="Input sequences used for GhostFold prediction",
                count=len(sequence_ids)
            )

        return output

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including GhostFold-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "ghostfold_params": {
                "name": self.name,
                "msa_only": self.msa_only,
                "num_recycle": self.num_recycle,
                "num_models": self.num_models,
                "num_seeds": self.num_seeds,
                "subsample": self.subsample,
                "mask_msa": self.mask_msa,
                "input_type": "tool_output" if self.input_is_tool_output else "direct"
            }
        })
        return base_dict
