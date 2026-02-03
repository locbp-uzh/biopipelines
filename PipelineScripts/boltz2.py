"""
Boltz2 configuration for protein-ligand complex prediction.

Handles apo and holo structure prediction with MSA caching,
ligand binding affinity calculation, and comprehensive analysis.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream, create_map_table
    from .combinatorics import generate_combinatorics_config, get_mode, predict_output_ids, Bundle, Each
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream, create_map_table
    from combinatorics import generate_combinatorics_config, get_mode, predict_output_ids, Bundle, Each


def _unwrap_combinatorics(value):
    """
    Unwrap Bundle/Each wrappers to get the underlying source.

    Returns the first source inside the wrapper, or the value itself if not wrapped.
    """
    if isinstance(value, (Bundle, Each)):
        if value.sources:
            return _unwrap_combinatorics(value.sources[0])
    return value


class Boltz2(BaseConfig):
    """
    Boltz2 configuration for protein-ligand complex prediction.

    Predicts both apo (protein-only) and holo (protein-ligand) structures
    with automatic MSA management and comprehensive scoring.
    """

    TOOL_NAME = "Boltz2"

    # Path descriptors - lazy evaluation after output_folder is set
    apo_config_folder = Path(lambda self: os.path.join(self.output_folder, "ApoConfig"))
    bound_config_folder = Path(lambda self: os.path.join(self.output_folder, "BoundConfig"))
    library_folder = Path(lambda self: os.path.join(self.output_folder, "Library"))
    apo_prediction_folder = Path(lambda self: os.path.join(self.output_folder, "ApoPredictions"))
    bound_prediction_folder = Path(lambda self: os.path.join(self.output_folder, "BoundPredictions"))
    msa_cache_folder = Path(lambda self: os.path.join(self.output_folder, "MSAs"))
    config_files_dir = Path(lambda self: os.path.join(self.output_folder, "config_files"))

    # Configuration files
    base_config_file = Path(lambda self: os.path.join(self.output_folder, "base_config.txt"))
    bound_config_file = Path(lambda self: os.path.join(self.bound_config_folder, "bound_config.yaml"))
    queries_csv = Path(lambda self: os.path.join(self.output_folder, f"{self.get_effective_job_name() or 'prediction'}_queries.csv"))
    queries_fasta = Path(lambda self: os.path.join(self.output_folder, f"{self.get_effective_job_name() or 'prediction'}_queries.fasta"))
    expanded_library_csv = Path(lambda self: os.path.join(self.output_folder, "expanded_smiles_library.csv"))
    fasta_files_list_file = Path(lambda self: os.path.join(self.output_folder, ".input_fasta_files.txt"))
    sequence_ids_file = Path(lambda self: os.path.join(self.output_folder, "sequence_ids.csv"))

    # Output files
    library_scores_csv = Path(lambda self: os.path.join(self.library_folder, "library_scores.csv"))
    config_txt_file = Path(lambda self: os.path.join(self.output_folder, f"{self.get_effective_job_name() or 'prediction'}_config.txt"))
    results_zip = Path(lambda self: os.path.join(self.output_folder, f"{self.get_effective_job_name() or 'prediction'}.zip"))
    confidence_csv = Path(lambda self: os.path.join(self.output_folder, "confidence_scores.csv"))
    affinity_csv = Path(lambda self: os.path.join(self.output_folder, "affinity_scores.csv"))
    sequences_csv = Path(lambda self: os.path.join(self.output_folder, "sequences.csv"))
    msas_csv = Path(lambda self: os.path.join(self.output_folder, "msas.csv"))
    missing_csv = Path(lambda self: os.path.join(self.output_folder, "missing.csv"))
    ligands_csv = Path(lambda self: os.path.join(self.output_folder, "ligands.csv"))

    # Helper script paths
    boltz_config_unified_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_boltz_config_unified.py"))
    boltz_postprocessing_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_boltz_postprocessing.py"))
    boltz_msa_copy_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_boltz_msa_copy.py"))
    boltz_completion_check_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_boltz_completion_check.py"))
    propagate_missing_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_propagate_missing.py"))

    def __init__(self,
                 # Primary input parameters
                 config: Optional[str] = None,
                 proteins: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ligands: Optional[Union[str, DataStream, StandardizedOutput]] = None,
                 msas: Optional[StandardizedOutput] = None,
                 # Core prediction parameters
                 affinity: bool = True,
                 output_format: str = "pdb",
                 msa_server: str = "public",
                 # Advanced prediction parameters
                 recycling_steps: Optional[int] = None,
                 diffusion_samples: Optional[int] = None,
                 use_potentials: bool = False,
                 # Template parameters
                 template: Optional[str] = None,
                 template_chain_ids: Optional[List[str]] = None,
                 template_force: bool = True,
                 template_threshold: float = 5.0,
                 # Pocket constraint parameters
                 pocket_residues: Optional[List[int]] = None,
                 pocket_max_distance: float = 7.0,
                 pocket_force: bool = True,
                 # Glycosylation parameters
                 glycosylation: Optional[Dict[str, List[int]]] = None,
                 # Covalent linkage parameters
                 covalent_linkage: Optional[Dict[str, Any]] = None,
                 **kwargs):
        """
        Initialize Boltz2 configuration.

        Args:
            config: Direct YAML configuration string
            proteins: Protein sequences as DataStream or StandardizedOutput
            ligands: Single SMILES string, DataStream, or StandardizedOutput with compounds
            msas: MSA files from previous Boltz2 run (StandardizedOutput with msas table)
            affinity: Whether to calculate binding affinity
            output_format: Output format ("pdb" or "mmcif")
            msa_server: MSA generation ("public" or "local")
            recycling_steps: Number of recycling steps
            diffusion_samples: Number of diffusion samples
            use_potentials: Enable potentials for improved structure prediction
            template: Path to PDB template file
            template_chain_ids: List of chain IDs to apply template to
            template_force: Whether to force template usage
            template_threshold: RMSD threshold for template matching
            pocket_residues: List of residue positions defining binding pocket
            pocket_max_distance: Maximum distance for pocket constraint
            pocket_force: Whether to force pocket constraint
            glycosylation: Dict mapping chain IDs to Asn positions for N-glycosylation
            covalent_linkage: Dict specifying covalent attachment
            **kwargs: Additional parameters
        """
        self.config = config
        self.msas = msas

        # Combinatorics modes
        self._proteins_mode = get_mode(proteins)
        self._ligands_mode = get_mode(ligands)

        # Unwrap Bundle/Each to get the underlying source
        unwrapped_proteins = _unwrap_combinatorics(proteins)
        unwrapped_ligands = _unwrap_combinatorics(ligands)

        # Store raw inputs for combinatorics config generation
        self.proteins = proteins
        self.ligands = ligands

        # Resolve protein input to DataStream
        self.proteins_stream: Optional[DataStream] = None
        if unwrapped_proteins is not None:
            if isinstance(unwrapped_proteins, StandardizedOutput):
                self.proteins_stream = unwrapped_proteins.streams.sequences
            elif isinstance(unwrapped_proteins, DataStream):
                self.proteins_stream = unwrapped_proteins
            else:
                raise ValueError(f"proteins must be DataStream or StandardizedOutput, got {type(unwrapped_proteins)}")

        # Resolve ligand input
        self.ligands_stream: Optional[DataStream] = None
        self.ligands_smiles: Optional[str] = None
        if unwrapped_ligands is not None:
            if isinstance(unwrapped_ligands, StandardizedOutput):
                self.ligands_stream = unwrapped_ligands.streams.compounds
            elif isinstance(unwrapped_ligands, DataStream):
                self.ligands_stream = unwrapped_ligands
            elif isinstance(unwrapped_ligands, str):
                # Direct SMILES string
                self.ligands_smiles = unwrapped_ligands
            else:
                raise ValueError(f"ligands must be str (SMILES), DataStream, or StandardizedOutput, got {type(unwrapped_ligands)}")

        # Override affinity to False if no ligands
        if self.ligands is None and self.ligands_smiles is None and affinity:
            print("Warning: No ligands detected, setting affinity=False")
            self.affinity = False
        else:
            self.affinity = affinity

        self.output_format = output_format
        self.msa_server = msa_server
        self.recycling_steps = recycling_steps
        self.diffusion_samples = diffusion_samples
        self.use_potentials = use_potentials

        # Template parameters
        self.template = template
        self.template_chain_ids = template_chain_ids
        self.template_force = template_force
        self.template_threshold = template_threshold

        # Pocket constraint parameters
        self.pocket_residues = pocket_residues
        self.pocket_max_distance = pocket_max_distance
        self.pocket_force = pocket_force

        # Glycosylation and covalent linkage
        self.glycosylation = glycosylation
        self.covalent_linkage = covalent_linkage

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate Boltz2-specific parameters."""
        # Must have some form of input
        has_input = any([
            self.config is not None,
            self.proteins_stream is not None
        ])
        if not has_input:
            raise ValueError("Either config or proteins parameter is required")

        # Cannot specify both config and proteins
        if self.config is not None and self.proteins_stream is not None:
            raise ValueError("Cannot specify both config and proteins")

        # Validate enum values
        if self.output_format not in ["pdb", "mmcif"]:
            raise ValueError("output_format must be 'pdb' or 'mmcif'")

        if self.msa_server not in ["public", "local"]:
            raise ValueError("msa_server must be 'public' or 'local'")

        if self.recycling_steps is not None and (not isinstance(self.recycling_steps, int) or self.recycling_steps < 1):
            raise ValueError("recycling_steps must be a positive integer")

        if self.diffusion_samples is not None and (not isinstance(self.diffusion_samples, int) or self.diffusion_samples < 1):
            raise ValueError("diffusion_samples must be a positive integer")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files."""
        self.folders = pipeline_folders

    def _write_combinatorics_config(self) -> str:
        """Write combinatorics config file at pipeline time."""
        config_path = os.path.join(self.output_folder, "combinatorics_config.json")
        generate_combinatorics_config(
            config_path,
            sequences=self.proteins,
            compounds=self.ligands
        )
        return config_path

    def _generate_extra_config_params(self) -> str:
        """Generate extra CLI parameters for pipe_boltz_config_unified.py."""
        extra_params = []

        if self.template:
            extra_params.append(f'--template "{self.template}"')
            if self.template_chain_ids:
                extra_params.append(f'--template-chains "{",".join(self.template_chain_ids)}"')
            if self.template_force:
                extra_params.append('--template-force')
            extra_params.append(f'--template-threshold {self.template_threshold}')

        if self.pocket_residues:
            extra_params.append(f'--pocket-residues "{self.pocket_residues}"')
            extra_params.append(f'--pocket-max-distance {self.pocket_max_distance}')
            if self.pocket_force:
                extra_params.append('--pocket-force')

        if self.glycosylation:
            glyco_json = json.dumps(self.glycosylation)
            extra_params.append(f"--glycosylation '{glyco_json}'")

        if self.covalent_linkage:
            covalent_json = json.dumps(self.covalent_linkage)
            extra_params.append(f"--covalent-linkage '{covalent_json}'")

        return " ".join(extra_params)

    def _get_msa_table_flag(self) -> str:
        """Get the MSA table flag for the config generator."""
        if not self.msas:
            return ""

        if hasattr(self.msas, 'tables'):
            if hasattr(self.msas.tables, '_tables') and 'msas' in self.msas.tables._tables:
                msa_table_path = self.msas.tables._tables['msas'].path
                return f'--msa-table "{msa_table_path}"'
            elif hasattr(self.msas.tables, 'msas'):
                msa_table = self.msas.tables.msas
                if hasattr(msa_table, 'path'):
                    return f'--msa-table "{msa_table.path}"'

        return ""

    def generate_script(self, script_path: str) -> str:
        """Generate bash script for Boltz2 execution."""
        boltz_cache_folder = self.folders["BoltzCache"]
        msa_option = "" if self.msa_server == "local" else " --use_msa_server"

        script_content = "#!/bin/bash\n"
        script_content += "# Boltz2 execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        # Create basic folder structure
        script_content += f"""# Create output folders
mkdir -p {os.path.join(self.output_folder, "predictions")}
mkdir -p {self.msa_cache_folder}

"""

        # Handle MSA recycling if provided
        if self.msas:
            script_content += self._generate_msa_recycling_section()

        # Write combinatorics config
        combinatorics_config_path = self._write_combinatorics_config()

        effective_job_name = self.get_effective_job_name() or "prediction"
        config_file_path = os.path.join(self.output_folder, f"{effective_job_name}.yaml")

        if self.config:
            # Direct YAML configuration
            script_content += f"""
echo "Using direct YAML configuration"
cat > {config_file_path} << 'EOF'
{self.config}
EOF

"""
            uses_unified_config = False
        else:
            # Use unified config generation
            script_content += self._generate_unified_config_section(combinatorics_config_path, effective_job_name)
            uses_unified_config = True

        # Build Boltz2 options
        boltz_options = f"--cache {boltz_cache_folder} --out_dir {self.output_folder}{msa_option} --output_format {self.output_format}"

        if self.recycling_steps is not None:
            boltz_options += f" --recycling_steps {self.recycling_steps}"

        if self.diffusion_samples is not None:
            boltz_options += f" --diffusion_samples {self.diffusion_samples}"

        if self.use_potentials:
            boltz_options += " --use_potentials"

        # Run Boltz2 prediction
        if uses_unified_config:
            script_content += f"""
echo "Running Boltz2 prediction on individual config files"
for config_file in {self.config_files_dir}/*.yaml; do
    if [ -f "$config_file" ]; then
        echo "Processing config: $config_file"
        boltz predict "$config_file" {boltz_options}
    fi
done

"""
        else:
            script_content += f"""
echo "Running Boltz2 prediction"
boltz predict {config_file_path} {boltz_options}

"""

        # Post-process results
        script_content += self._generate_postprocess_section()

        # Propagate missing table
        script_content += self._generate_missing_table_propagation()

        # Completion check
        script_content += self._generate_completion_check_section()

        return script_content

    def _generate_unified_config_section(self, combinatorics_config_path: str, config_name: str) -> str:
        """Generate script section for unified config generation."""
        script = f"""
echo "Generating Boltz2 configurations using unified config generator"
mkdir -p {self.config_files_dir}

"""
        # Generate ligands CSV if using direct SMILES string
        if self.ligands_smiles:
            ligand_id = config_name if config_name != "prediction" else "ligand"
            script += f"""# Create ligands CSV from direct SMILES string
cat > {self.ligands_csv} << 'EOF'
id,format,smiles,ccd
{ligand_id},smiles,{self.ligands_smiles},
EOF

"""

        # Build command for unified config generator
        cmd_parts = [
            f'python {self.boltz_config_unified_py}',
            f'--combinatorics-config "{combinatorics_config_path}"',
            f'--output-dir "{self.config_files_dir}"',
        ]

        msa_table_flag = self._get_msa_table_flag()
        if msa_table_flag:
            cmd_parts.append(msa_table_flag)

        if self.affinity:
            cmd_parts.append('--affinity')

        extra_params = self._generate_extra_config_params()
        if extra_params:
            cmd_parts.append(extra_params)

        script += '# Generate config files\n'
        script += ' \\\n    '.join(cmd_parts) + '\n\n'

        return script

    def _generate_msa_recycling_section(self) -> str:
        """Generate script section for MSA recycling."""
        if not self.msas:
            return ""

        # Get MSA table path from msas input
        msa_table_path = None
        if hasattr(self.msas, 'tables'):
            if hasattr(self.msas.tables, '_tables') and 'msas' in self.msas.tables._tables:
                msa_table_path = self.msas.tables._tables['msas'].path
            elif hasattr(self.msas.tables, 'msas'):
                msa_table = self.msas.tables.msas
                if hasattr(msa_table, 'path'):
                    msa_table_path = msa_table.path

        if msa_table_path:
            return f"""
echo "Recycling MSAs from previous prediction"
python {self.boltz_msa_copy_py} \\
    --msa-table "{msa_table_path}" \\
    --output-folder "{self.msa_cache_folder}"

"""
        return ""

    def _generate_postprocess_section(self) -> str:
        """Generate script section for post-processing."""
        return f"""
echo "Post-processing Boltz2 results"
python {self.boltz_postprocessing_py} {self.output_folder} {self.output_folder} {self.sequence_ids_file}

if [ $? -ne 0 ]; then
    echo "Error: Post-processing failed"
    exit 1
fi

echo "Post-processing completed"

"""

    def _generate_missing_table_propagation(self) -> str:
        """Generate script section to propagate missing.csv from upstream tools."""
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.proteins,
            self.proteins_stream
        )

        if not upstream_missing_path:
            return ""

        upstream_folder = os.path.dirname(upstream_missing_path)
        structure_ext = ".pdb" if self.output_format == "pdb" else ".cif"
        msa_ext = ".csv" if self.msa_server == "public" else ".a3m"

        return f"""
# Propagate missing table from upstream tools
echo "Checking for upstream missing sequences..."
if [ -f "{upstream_missing_path}" ]; then
    echo "Found upstream missing.csv - propagating to current tool"
    python {self.propagate_missing_py} \\
        --upstream-folders "{upstream_folder}" \\
        --output-folder "{self.output_folder}" \\
        --structure-ext "{structure_ext}" \\
        --msa-ext "{msa_ext}"
else
    echo "No upstream missing.csv found"
fi

"""

    def _generate_completion_check_section(self) -> str:
        """Generate completion check section using helper script."""
        # Check if there's an upstream missing table
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.proteins,
            self.proteins_stream
        )

        if upstream_missing_path:
            # Use helper script for completion check with missing sequence filtering
            return f"""
# Completion check with missing sequence filtering
python {self.boltz_completion_check_py} \\
    --output-folder "{self.output_folder}" \\
    --sequence-ids-file "{self.sequence_ids_file}" \\
    --missing-file "{upstream_missing_path}" \\
    --output-format "{self.output_format}"

if [ $? -ne 0 ]; then
    echo "Boltz2 failed - some outputs missing"
    exit 1
fi

echo "Boltz2 completed successfully"
""" + self.generate_completion_check_footer()
        else:
            # Standard completion check
            return self.generate_completion_check_footer()

    def _predict_sequence_ids(self) -> List[str]:
        """Predict sequence IDs from input sources using combinatorics module."""
        return predict_output_ids(
            bundled_name="bundled_complex",
            sequences=self.proteins,
            compounds=self.ligands
        )

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after Boltz2 execution."""
        predicted_ids = self._predict_sequence_ids()

        # Structure files
        structure_ext = ".pdb" if self.output_format == "pdb" else ".cif"
        structure_files = [os.path.join(self.output_folder, f"{seq_id}{structure_ext}") for seq_id in predicted_ids]

        # Create structures DataStream
        structures_map = os.path.join(self.output_folder, "structures_map.csv")
        create_map_table(structures_map, predicted_ids, files=structure_files)

        structures = DataStream(
            name="structures",
            ids=predicted_ids,
            files=structure_files,
            map_table=structures_map,
            format=self.output_format
        )

        # MSA files
        msa_ext = ".csv" if self.msa_server == "public" else ".a3m"
        msa_files = [os.path.join(self.msa_cache_folder, f"{seq_id}{msa_ext}") for seq_id in predicted_ids]

        msas = DataStream(
            name="msas",
            ids=predicted_ids,
            files=msa_files,
            map_table=self.msas_csv,
            format="csv" if self.msa_server == "public" else "a3m"
        )

        # Sequences DataStream (from input)
        sequence_ids = list(self.proteins_stream.ids) if self.proteins_stream else predicted_ids
        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[],  # Sequences are value-based, stored in map_table
            map_table=self.sequences_csv,
            format="sequence"
        )

        # Compounds DataStream
        if self.ligands_stream:
            compounds = self.ligands_stream
        elif self.ligands_smiles:
            # Single SMILES - create compound entry
            effective_job_name = self.get_effective_job_name() or "ligand"
            compound_ids = [effective_job_name]
            compounds = DataStream(
                name="compounds",
                ids=compound_ids,
                files=[],
                map_table=self.ligands_csv,
                format="smiles"
            )
        else:
            compounds = DataStream.empty("compounds", "smiles")

        # Tables
        tables = {
            "structures": TableInfo(
                name="structures",
                path=structures_map,
                columns=["id", "file"],
                description="Boltz2 predicted structures",
                count=len(predicted_ids)
            ),
            "confidence": TableInfo(
                name="confidence",
                path=self.confidence_csv,
                columns=["id", "input_file", "confidence_score", "ptm", "iptm", "complex_plddt", "complex_iplddt"],
                description="Boltz2 confidence scores",
                count=len(predicted_ids)
            ),
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence"],
                description="Input protein sequences",
                count=len(sequence_ids)
            ),
            "msas": TableInfo(
                name="msas",
                path=self.msas_csv,
                columns=["id", "sequence_id", "sequence", "msa_file"],
                description="MSA files for recycling",
                count=len(predicted_ids)
            )
        }

        if self.affinity:
            tables["affinity"] = TableInfo(
                name="affinity",
                path=self.affinity_csv,
                columns=["id", "input_file", "affinity_pred_value", "affinity_probability_binary"],
                description="Boltz2 affinity predictions",
                count=len(predicted_ids)
            )

        if self.ligands_stream or self.ligands_smiles:
            tables["compounds"] = TableInfo(
                name="compounds",
                path=self.ligands_csv if self.ligands_smiles else (self.ligands_stream.map_table if self.ligands_stream else ""),
                columns=["id", "format", "smiles", "ccd"],
                description="Ligand compounds",
                count=len(compounds.ids) if compounds else 0
            )

        # Check for upstream missing table
        upstream_missing_path = self._get_upstream_missing_table_path(
            self.proteins,
            self.proteins_stream
        )
        if upstream_missing_path:
            tables["missing"] = TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "structure", "msa"],
                description="Sequences filtered out by upstream tools",
                count="variable"
            )

        return {
            "structures": structures,
            "sequences": sequences,
            "compounds": compounds,
            "msas": msas,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def get_config_display(self) -> List[str]:
        """Get configuration display lines for pipeline output."""
        config_lines = super().get_config_display()

        if self.ligands_smiles:
            config_lines.append(f"Ligands: {self.ligands_smiles}")
        elif self.ligands_stream:
            config_lines.append(f"Ligands: {len(self.ligands_stream)} compounds")

        config_lines.extend([
            f"Output format: {self.output_format}",
            f"MSA server: {self.msa_server}",
            f"Affinity calculation: {self.affinity}"
        ])

        if self.recycling_steps is not None:
            config_lines.append(f"Recycling steps: {self.recycling_steps}")

        if self.diffusion_samples is not None:
            config_lines.append(f"Diffusion samples: {self.diffusion_samples}")

        if self.use_potentials:
            config_lines.append(f"Use potentials: {self.use_potentials}")

        if self.template:
            config_lines.append(f"Template: {os.path.basename(self.template)}")
            if self.template_chain_ids:
                config_lines.append(f"Template chains: {', '.join(self.template_chain_ids)}")

        if self.pocket_residues:
            config_lines.append(f"Pocket residues: {self.pocket_residues}")

        if self.glycosylation:
            config_lines.append(f"Glycosylation: {self.glycosylation}")

        if self.covalent_linkage:
            config_lines.append(f"Covalent linkage: {self.covalent_linkage}")

        return config_lines

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "boltz2_params": {
                "config": self.config,
                "proteins": str(self.proteins) if self.proteins else None,
                "ligands": str(self.ligands) if self.ligands else self.ligands_smiles,
                "msas": str(self.msas) if self.msas else None,
                "affinity": self.affinity,
                "output_format": self.output_format,
                "msa_server": self.msa_server,
                "recycling_steps": self.recycling_steps,
                "diffusion_samples": self.diffusion_samples,
                "use_potentials": self.use_potentials,
                "template": self.template,
                "template_chain_ids": self.template_chain_ids,
                "template_force": self.template_force,
                "template_threshold": self.template_threshold,
                "pocket_residues": self.pocket_residues,
                "pocket_max_distance": self.pocket_max_distance,
                "pocket_force": self.pocket_force,
                "glycosylation": self.glycosylation,
                "covalent_linkage": self.covalent_linkage
            }
        })
        return base_dict
