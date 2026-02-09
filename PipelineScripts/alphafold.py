# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
AlphaFold2/ColabFold configuration for protein structure prediction.

Handles sequence folding with AlphaFold2 using ColabFold implementation,
integrating with upstream sequence generation tools and providing comprehensive
ranking and analysis capabilities.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream, create_map_table
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream, create_map_table


class AlphaFold(BaseConfig):
    """
    Configuration for AlphaFold2/ColabFold structure prediction.

    Predicts protein structures from amino acid sequences.

    Example:
        # Using Sequence tool
        proteins = Sequence(["MKTVRQ...", "AETGFT..."], ids=["p1", "p2"])
        af = AlphaFold(proteins=proteins)

        # Using output from another tool
        af = AlphaFold(proteins=mpnn_output)
    """

    TOOL_NAME = "AlphaFold"

    # Lazy path descriptors
    queries_csv = Path(lambda self: os.path.join(self.output_folder, f"{self.pipeline_name}_queries.csv"))
    queries_fasta = Path(lambda self: os.path.join(self.output_folder, f"{self.pipeline_name}_queries.fasta"))
    confidence_csv = Path(lambda self: os.path.join(self.output_folder, "confidence.csv"))
    folding_folder = Path(lambda self: os.path.join(self.output_folder, "Folding"))
    msas_folder = Path(lambda self: os.path.join(self.output_folder, "MSAs"))
    msa_csv = Path(lambda self: os.path.join(self.output_folder, "msas.csv"))
    structures_map = Path(lambda self: os.path.join(self.output_folder, "structures_map.csv"))
    missing_csv = Path(lambda self: os.path.join(self.output_folder, "missing.csv"))

    # Helper script paths
    colabfold_batch = Path(lambda self: os.path.join(self.folders["AlphaFold"], "colabfold-conda/bin/colabfold_batch"))
    fa_to_csv_fasta_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_fa_to_csv_fasta.py"))
    alphafold_confidence_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_alphafold_confidence.py"))
    alphafold_msas_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_alphafold_msas.py"))
    propagate_missing_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_propagate_missing.py"))

    def __init__(self,
                 proteins: Union[DataStream, StandardizedOutput],
                 num_relax: int = 0,
                 num_recycle: int = 3,
                 rand_seed: int = 0,
                 **kwargs):
        """
        Initialize AlphaFold configuration.

        Args:
            proteins: Input protein sequences as DataStream or StandardizedOutput.
                      Use Sequence("MKTVRQ...") to create from raw sequence strings.
            num_relax: Number of best models to relax with AMBER
            num_recycle: Number of recycling iterations (default 3)
            rand_seed: Random seed for reproducible results (0 = random)
        """
        # Store original input for upstream missing table lookup
        self.proteins = proteins

        # Resolve input to DataStream
        if isinstance(proteins, StandardizedOutput):
            self.sequences_stream: DataStream = proteins.streams.sequences
        elif isinstance(proteins, DataStream):
            self.sequences_stream = proteins
        else:
            raise ValueError(f"proteins must be DataStream or StandardizedOutput, got {type(proteins)}")

        # Store AlphaFold-specific parameters
        self.num_relax = num_relax
        self.num_recycle = num_recycle
        self.rand_seed = rand_seed

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate AlphaFold-specific parameters."""
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("proteins parameter is required and must not be empty")

        if self.num_relax < 0:
            raise ValueError("num_relax cannot be negative")

        if self.num_recycle < 1:
            raise ValueError("num_recycle must be at least 1")

        if self.rand_seed < 0:
            raise ValueError("rand_seed cannot be negative")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get AlphaFold configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"INPUT PROTEINS: {len(self.sequences_stream)} sequences",
            f"NUM RELAX: {self.num_relax}",
            f"NUM RECYCLE: {self.num_recycle}"
        ])

        if self.rand_seed > 0:
            config_lines.append(f"RAND SEED: {self.rand_seed}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate AlphaFold execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# AlphaFold execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_prepare_sequences()
        script_content += self._generate_script_run_alphafold()
        script_content += self._generate_script_extract_best_rank()
        script_content += self._generate_script_extract_confidence()
        script_content += self._generate_script_create_msas_table()
        script_content += self._generate_missing_table_propagation()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_prepare_sequences(self) -> str:
        """Generate the sequence preparation part of the script."""
        # Get source file from input DataStream
        source_file = self.sequences_stream.map_table

        # Check if we need to convert from FASTA folder to CSV
        # Only treat as ProteinMPNN seqs folder if it ends with /seqs or is a directory
        if source_file.endswith("seqs") and not source_file.endswith(".csv"):
            # ProteinMPNN output - convert .fa files to CSV
            return f"""echo "Converting ProteinMPNN .fa files to queries CSV"
python {self.fa_to_csv_fasta_py} {source_file} {self.queries_csv} {self.queries_fasta}

"""
        else:
            # Direct CSV or FASTA file - copy with error checking
            return f"""echo "Using sequences from: {source_file}"
if [ -f "{source_file}" ]; then
    cp "{source_file}" "{self.queries_csv}"
    echo "Successfully copied sequences file"
else
    echo "ERROR: Sequence file not found: {source_file}"
    echo "This usually means the previous step failed to generate the expected output"
    exit 1
fi

"""

    def _generate_script_run_alphafold(self) -> str:
        """Generate the AlphaFold execution part of the script."""
        # Build AlphaFold options
        af_options = ""
        if self.num_relax > 0:
            af_options += f" --amber --use-gpu-relax --num-relax {self.num_relax}"
        if self.num_recycle != 3:
            af_options += f" --num-recycle {self.num_recycle}"
        if self.rand_seed != 0:
            af_options += f" --random-seed {self.rand_seed}"

        return f"""echo "Running AlphaFold2/ColabFold"
echo "Options: {af_options}"
echo "Output folder: {self.output_folder}"

# Create Folding subfolder for raw ColabFold outputs
mkdir -p "{self.folding_folder}"

# Run ColabFold batch - output to Folding subfolder
{self.colabfold_batch} {self.queries_csv} "{self.folding_folder}" {af_options}

"""

    def _generate_script_extract_best_rank(self) -> str:
        """Generate script to extract best rank structures."""
        msa_section = f"""
# Create MSAs subfolder
mkdir -p "{self.msas_folder}"

# Copy MSA files to MSAs subfolder
echo "Extracting MSA files"
for msa_file in *.a3m; do
    if [ -f "$msa_file" ]; then
        cp "$msa_file" "{self.msas_folder}/"
        echo "Copied MSA: $msa_file"
    fi
done
"""

        return f"""echo "Extracting best structures from Folding subfolder"
# AlphaFold creates files like:
# - sequenceid_unrelaxed_rank_001_alphafold2_ptm_model_N_seed_SSS.pdb
# - sequenceid_relaxed_rank_001_alphafold2_ptm_model_N_seed_SSS.pdb
# - sequenceid.a3m (MSA file, not generated in single_sequence mode)
# We want to copy the best ones to main directory as: sequenceid.pdb

cd "{self.folding_folder}"

# Handle relaxed format (preferred if both exist)
for file in *_relaxed_rank_001_*.pdb; do
    if [ -f "$file" ]; then
        # Extract sequence ID (everything before _relaxed_rank_001)
        base=$(echo "$file" | sed 's/_relaxed_rank_001_.*/.pdb/')
        cp "$file" "{self.output_folder}/$base"
        echo "Extracted: $file -> $base"
    fi
done

# Handle unrelaxed format (only if relaxed doesn't exist)
for file in *_unrelaxed_rank_001_*.pdb; do
    if [ -f "$file" ]; then
        # Extract sequence ID (everything before _unrelaxed_rank_001)
        base=$(echo "$file" | sed 's/_unrelaxed_rank_001_.*/.pdb/')
        # Only copy if the relaxed version doesn't already exist
        if [ ! -f "{self.output_folder}/$base" ]; then
            cp "$file" "{self.output_folder}/$base"
            echo "Extracted: $file -> $base"
        else
            echo "Skipped $file (relaxed version already exists as $base)"
        fi
    fi
done
{msa_section}
cd - > /dev/null

"""

    def _generate_script_extract_confidence(self) -> str:
        """Generate script section to extract confidence metrics from JSON files."""
        return f"""echo "Extracting confidence metrics from JSON files"
python {self.alphafold_confidence_py} "{self.folding_folder}" "{self.output_folder}" "{self.confidence_csv}"

"""

    def _generate_script_create_msas_table(self) -> str:
        """Generate script section to create MSAs CSV table."""
        return f"""echo "Creating MSAs table"
python {self.alphafold_msas_py} "{self.msas_folder}" "{self.queries_csv}" "{self.msa_csv}"

"""

    def _generate_missing_table_propagation(self) -> str:
        """Generate script section to propagate missing.csv from upstream tools."""
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.proteins,
            self.sequences_stream
        )

        if not upstream_missing_path:
            return ""

        upstream_folder = os.path.dirname(upstream_missing_path)

        return f"""
# Propagate missing table from upstream tools
echo "Checking for upstream missing sequences..."
if [ -f "{upstream_missing_path}" ]; then
    echo "Found upstream missing.csv - propagating to current tool"
    python {self.propagate_missing_py} \\
        --upstream-folders "{upstream_folder}" \\
        --output-folder "{self.output_folder}"
else
    echo "No upstream missing.csv found - creating empty missing.csv"
    echo "id,removed_by,cause" > "{self.missing_csv}"
fi

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after AlphaFold execution."""
        # Use sequence IDs from input to predict structure files
        sequence_ids = self.sequences_stream.ids

        # Generate structure file paths
        structure_files = []
        for seq_id in sequence_ids:
            pdb_path = os.path.join(self.output_folder, f"{seq_id}.pdb")
            structure_files.append(pdb_path)

        # Create map_table for structures
        create_map_table(self.structures_map, sequence_ids, files=structure_files)

        structures = DataStream(
            name="structures",
            ids=sequence_ids,
            files=structure_files,
            map_table=self.structures_map,
            format="pdb"
        )

        # MSA files (AlphaFold generates .a3m files)
        msa_ids = []
        msa_files = []
        for seq_id in sequence_ids:
            msa_file = os.path.join(self.msas_folder, f"{seq_id}.a3m")
            msa_ids.append(seq_id)
            msa_files.append(msa_file)

        msas = DataStream(
            name="msas",
            ids=msa_ids,
            files=msa_files,
            map_table=self.msa_csv,
            format="a3m"
        )

        # Tables
        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_map,
                columns=["id", "file"],
                description="AlphaFold predicted structures",
                count=len(sequence_ids)
            ),
            "confidence": TableInfo(
                name="confidence",
                path=self.confidence_csv,
                columns=["id", "structure", "plddt", "max_pae", "ptm"],
                description="AlphaFold confidence metrics extracted from best rank models",
                count=len(sequence_ids)
            ),
            "msas": TableInfo(
                name="msas",
                path=self.msa_csv,
                columns=["id", "sequence_id", "sequence", "msa_file"],
                description="MSA files for sequence recycling between predictions",
                count=len(msa_files)
            )
        }

        # Check for upstream missing table
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.proteins,
            self.sequences_stream
        )
        if upstream_missing_path:
            tables["missing"] = TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "cause"],
                description="IDs removed by upstream tools with removal reason",
                count="variable"
            )

        return {
            "structures": structures,
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "msas": msas,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including AlphaFold-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "af_params": {
                "num_relax": self.num_relax,
                "num_recycle": self.num_recycle,
                "rand_seed": self.rand_seed
            }
        })
        return base_dict
