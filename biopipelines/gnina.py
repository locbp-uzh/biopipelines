# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
GNINA docking and scoring configuration.

GNINA combines AutoDock Vina search with a convolutional neural network (CNN)
for pose scoring, giving more accurate binding pose prediction than Vina alone.
Supports single-ligand docking, compound-library screening, and multi-conformer
docking with statistical analysis across independent runs.

The binding box can be defined in three ways (checked in order):
  1. Explicit center + size parameters
  2. autobox_ligand — path or DataStream of a reference ligand
  3. Auto-detected from crystal ligand HETATM records already in the input PDB

Reference:
    McNutt et al. (2021) GNINA 1.0: molecular docking with deep learning.
    J Cheminform 13, 43. https://github.com/gnina/gnina
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream, create_map_table
    from .config_manager import ConfigManager
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream, create_map_table
    from config_manager import ConfigManager


class Gnina(BaseConfig):
    """
    GNINA molecular docking with CNN-based pose scoring.

    Performs protein-ligand docking using the GNINA binary (AutoDock Vina + CNN
    scoring). Supports multi-conformer docking with statistical analysis across
    independent runs for robust ranking.

    Main workflows:
        - Single-ligand docking: dock one compound against one or more proteins.
        - Library screening: dock a CompoundLibrary against a target.
        - Conformer docking: generate RDKit conformers per ligand, dock each
          independently, and rank by combined Vina score + conformer strain energy.

    Outputs:
        - docking_results: all accepted poses with Vina and CNN scores.
        - conformer_ranking: per-conformer aggregated statistics including
          pose consistency (fraction of runs that converge to the same pose)
          and pseudo_binding_energy (best Vina + conformer strain energy).
        - best_poses/: combined protein+ligand PDB files for each best pose.

    Usage:
        with Pipeline(...):
            Resources(gpu="A100", memory="16GB", time="8:00:00")
            protein = PDB("9RTM")
            ligand = Ligand(smiles="CN(C)c1ccc2...", ids="TMR")
            docking = Gnina(structures=protein, compounds=ligand)
    """

    TOOL_NAME = "Gnina"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Install the GNINA binary on an HPC system.

        Downloads the pre-built GNINA binary from the official GitHub releases
        and creates a versioned symlink. After installation, configure the
        ``gnina`` section in config.yaml with the CUDA modules and
        ``ld_library_path`` appropriate for your HPC environment:

            gnina:
              modules:
                - "cuda/<version>"
                - "cudnn/<version>"
              ld_library_path: "/path/to/cuda/lib"

        To find available modules and library paths on your system::

            module avail cuda && module avail cudnn
            module show cuda/<version>
        """
        repo_dir = folders.get("Gnina", "")
        binary = os.path.join(repo_dir, "gnina.1.3.2")
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -f "{binary}" ]; then
    echo "GNINA already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        return f"""echo "=== Installing GNINA ==="
{skip}mkdir -p {repo_dir}
cd {repo_dir}
wget https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2
chmod +x gnina.1.3.2
ln -sf gnina.1.3.2 gnina

echo "=== GNINA installation complete ==="
"""

    # Lazy path descriptors
    gnina_binary = Path(lambda self: os.path.join(self.folders["Gnina"], "gnina"))
    helper_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_gnina.py"))
    docking_results_csv = Path(lambda self: os.path.join(self.output_folder, "docking_results.csv"))
    conformer_ranking_csv = Path(lambda self: os.path.join(self.output_folder, "conformer_ranking.csv"))
    config_json = Path(lambda self: os.path.join(self.output_folder, "gnina_config.json"))
    structures_json = Path(lambda self: os.path.join(self.output_folder, "structures_ds.json"))
    compounds_json = Path(lambda self: os.path.join(self.output_folder, "compounds_ds.json"))
    structures_map = Path(lambda self: os.path.join(self.output_folder, "structures_map.csv"))
    missing_csv = Path(lambda self: os.path.join(self.output_folder, "missing.csv"))
    propagate_missing_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_propagate_missing.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 compounds: Union[DataStream, StandardizedOutput],
                 autobox_ligand: Union[DataStream, StandardizedOutput, str, None] = None,
                 center: Optional[str] = None,
                 size: Union[float, str, None] = None,
                 autobox_add: float = 4.0,
                 exhaustiveness: int = 32,
                 num_modes: int = 9,
                 num_runs: int = 5,
                 seed: int = 42,
                 cnn_scoring: str = "rescore",
                 generate_conformers: bool = False,
                 num_conformers: int = 50,
                 energy_window: float = 2.0,
                 conformer_rmsd: float = 1.0,
                 conformer_energies: Optional[tuple] = None,
                 cnn_score_threshold: float = 0.5,
                 rmsd_threshold: float = 2.0,
                 protonate: bool = True,
                 pH: float = 7.4,
                 **kwargs):
        """
        Initialize GNINA docking configuration.

        Args:
            structures: Input protein structures as DataStream or StandardizedOutput.
            compounds: Input ligands as DataStream or StandardizedOutput. Accepts
                output from Ligand(), CompoundLibrary(), or any tool that produces
                a compounds stream (e.g. BoltzGen).
            autobox_ligand: Reference ligand for automatic binding box definition.
                Accepts DataStream, StandardizedOutput, or a file path string.
                If None and no center+size given, auto-detects from crystal ligand
                HETATM records in the input PDB.
            center: Explicit box center as "x,y,z" string (e.g. "10.5,-2.3,8.1").
                Required together with size when not using autobox.
            size: Box dimensions in Angstroms. Either a single float (cubic box)
                or "x,y,z" string for an asymmetric box.
            autobox_add: Padding added around the autobox ligand in Angstroms
                (default 4.0). Only used when autobox_ligand is set or
                auto-detected.
            exhaustiveness: Search exhaustiveness — higher values explore more
                conformational space but take longer (default 32).
            num_modes: Number of docked poses output per GNINA run (default 9).
            num_runs: Number of independent docking runs per conformer with
                incremented seeds. More runs improve pose consistency statistics
                (default 5).
            seed: Base random seed. Each run uses seed + run_index (default 42).
            cnn_scoring: CNN scoring mode. Options: "rescore" (re-score Vina poses
                with CNN, default), "refinement" (optimize with CNN), "none"
                (Vina only), "rescore_only" (CNN score only, no Vina).
            generate_conformers: If True, generate RDKit ETKDGv3 conformers for
                each ligand, filter by energy window, cluster by RMSD, and dock
                each representative independently (default False).
            num_conformers: Number of conformers to generate before energy
                filtering and RMSD clustering (default 50).
            energy_window: Maximum relative MMFF94 energy (kcal/mol above the
                minimum energy conformer) to retain after generation (default 2.0).
            conformer_rmsd: Heavy-atom RMSD cutoff (Å) for Butina clustering of
                conformers. Only the lowest-energy representative from each cluster
                is docked. Set to 0 to disable clustering (default 1.0).
            conformer_energies: Pre-computed conformer energies as a
                (TableInfo, "column_name") tuple. Used when generate_conformers
                is False but the input SDF contains multiple pre-computed
                conformers with known energies.
            cnn_score_threshold: Minimum CNNscore to accept a pose (0–1).
                Poses below this threshold are discarded (default 0.5).
            rmsd_threshold: RMSD cutoff (Å) for pose consistency clustering
                across runs. Larger clusters indicate more reproducible poses
                (default 2.0).
            protonate: If True, add hydrogens to the protein with OpenBabel at
                the specified pH before docking (default True). Requires obabel
                in PATH.
            pH: Protonation pH for OpenBabel (default 7.4).
        """
        # Keep original input for upstream missing table detection
        self.structures_input = structures

        # Resolve structures input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # Resolve compounds input to DataStream
        if isinstance(compounds, StandardizedOutput):
            self.compounds_stream: DataStream = compounds.streams.compounds
        elif isinstance(compounds, DataStream):
            self.compounds_stream = compounds
        else:
            raise ValueError(f"compounds must be DataStream or StandardizedOutput, got {type(compounds)}")

        # Resolve autobox_ligand
        self.autobox_ligand_path = None
        if autobox_ligand is not None:
            if isinstance(autobox_ligand, StandardizedOutput):
                stream = autobox_ligand.streams.structures
                if len(stream) > 0:
                    _, self.autobox_ligand_path = stream[0]
            elif isinstance(autobox_ligand, DataStream):
                if len(autobox_ligand) > 0:
                    _, self.autobox_ligand_path = autobox_ligand[0]
            elif isinstance(autobox_ligand, str):
                self.autobox_ligand_path = autobox_ligand
            else:
                raise ValueError(f"autobox_ligand must be DataStream, StandardizedOutput, or str, got {type(autobox_ligand)}")

        # Box parameters
        self.center = center
        self.size = size
        self.autobox_add = autobox_add

        # Docking parameters
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.num_runs = num_runs
        self.seed = seed
        self.cnn_scoring = cnn_scoring

        # Conformer parameters
        self.generate_conformers = generate_conformers
        self.num_conformers = num_conformers
        self.energy_window = energy_window
        self.conformer_rmsd = conformer_rmsd
        self.conformer_energies = conformer_energies

        # Filtering & analysis
        self.cnn_score_threshold = cnn_score_threshold
        self.rmsd_threshold = rmsd_threshold

        # Protonation
        self.protonate = protonate
        self.pH = pH

        # Resolve conformer_energies table reference if provided
        self.conformer_energies_ref = None
        if conformer_energies is not None:
            self.validate_table_reference(conformer_energies)
            self.conformer_energies_ref = self.resolve_table_reference(conformer_energies)

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate GNINA parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        if not self.compounds_stream or len(self.compounds_stream) == 0:
            raise ValueError("compounds parameter is required and must not be empty")

        has_explicit_box = self.center is not None and self.size is not None
        has_autobox = self.autobox_ligand_path is not None
        if not has_explicit_box and not has_autobox:
            print("  Note: No explicit box defined. Will auto-detect from crystal "
                  "ligand if available in the input PDB structure.")

        if self.exhaustiveness <= 0:
            raise ValueError("exhaustiveness must be positive")
        if self.num_modes <= 0:
            raise ValueError("num_modes must be positive")
        if self.num_runs <= 0:
            raise ValueError("num_runs must be positive")
        if not 0 <= self.cnn_score_threshold <= 1:
            raise ValueError("cnn_score_threshold must be between 0 and 1")
        if self.rmsd_threshold <= 0:
            raise ValueError("rmsd_threshold must be positive")
        if self.energy_window <= 0:
            raise ValueError("energy_window must be positive")
        if self.num_conformers <= 0:
            raise ValueError("num_conformers must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and folder paths."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get GNINA configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"STRUCTURES: {len(self.structures_stream)} proteins",
            f"COMPOUNDS: {len(self.compounds_stream)} ligands",
            f"EXHAUSTIVENESS: {self.exhaustiveness}",
            f"NUM RUNS: {self.num_runs}",
            f"NUM MODES: {self.num_modes}",
            f"SEED: {self.seed}",
            f"CNN SCORING: {self.cnn_scoring}",
            f"CNN SCORE THRESHOLD: {self.cnn_score_threshold}",
        ])
        if self.center is not None:
            config_lines.append(f"BOX CENTER: {self.center}")
            config_lines.append(f"BOX SIZE: {self.size}")
        if self.autobox_ligand_path is not None:
            config_lines.append(f"AUTOBOX LIGAND: {self.autobox_ligand_path}")
            config_lines.append(f"AUTOBOX ADD: {self.autobox_add}")
        if self.generate_conformers:
            config_lines.append(f"GENERATE CONFORMERS: {self.num_conformers}")
            config_lines.append(f"ENERGY WINDOW: {self.energy_window} kcal/mol")
            config_lines.append(f"CONFORMER RMSD CUTOFF: {self.conformer_rmsd} Å")
        if self.protonate:
            config_lines.append(f"PROTONATE: pH {self.pH}")
        return config_lines

    def _write_config_json(self):
        """Write configuration and DataStreams to JSON files at pipeline time."""
        os.makedirs(self.output_folder, exist_ok=True)

        with open(self.structures_json, 'w') as f:
            json.dump(self.structures_stream.to_dict(), f, indent=2)

        with open(self.compounds_json, 'w') as f:
            json.dump(self.compounds_stream.to_dict(), f, indent=2)

        box_config = {"autobox_add": self.autobox_add}
        if self.center is not None and self.size is not None:
            box_config["center"] = self.center
            box_config["size"] = self.size
        if self.autobox_ligand_path is not None:
            box_config["autobox_ligand"] = self.autobox_ligand_path

        upstream_missing = self._get_upstream_missing_table_path(self.structures_input)

        config = {
            "output_folder": self.output_folder,
            "gnina_binary": self.gnina_binary,
            "structures_json": self.structures_json,
            "compounds_json": self.compounds_json,
            "docking_results_csv": self.docking_results_csv,
            "conformer_ranking_csv": self.conformer_ranking_csv,
            "missing_csv": upstream_missing,
            "box": box_config,
            "exhaustiveness": self.exhaustiveness,
            "num_modes": self.num_modes,
            "num_runs": self.num_runs,
            "seed": self.seed,
            "cnn_scoring": self.cnn_scoring,
            "generate_conformers": self.generate_conformers,
            "num_conformers": self.num_conformers,
            "energy_window": self.energy_window,
            "conformer_rmsd": self.conformer_rmsd,
            "conformer_energies_ref": self.conformer_energies_ref,
            "cnn_score_threshold": self.cnn_score_threshold,
            "rmsd_threshold": self.rmsd_threshold,
            "protonate": self.protonate,
            "pH": self.pH,
        }

        with open(self.config_json, 'w') as f:
            json.dump(config, f, indent=2)

    def _get_gnina_modules(self) -> List[str]:
        """Get CUDA modules for GNINA from the gnina section of config.yaml."""
        config_manager = ConfigManager()
        gnina_config = config_manager._config.get('gnina', {})
        return gnina_config.get('modules', [])

    def _get_gnina_ld_library_path(self) -> str:
        """Get LD_LIBRARY_PATH for GNINA from the gnina section of config.yaml."""
        config_manager = ConfigManager()
        gnina_config = config_manager._config.get('gnina', {})
        return gnina_config.get('ld_library_path', '')

    def generate_script(self, script_path: str) -> str:
        """Generate GNINA execution script."""
        self._write_config_json()

        script_content = "#!/bin/bash\n"
        script_content += "# GNINA docking execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        modules = self._get_gnina_modules()
        if modules:
            script_content += "# Load CUDA modules for GNINA\n"
            script_content += "module load " + " ".join(modules) + "\n"

        ld_path = self._get_gnina_ld_library_path()
        if ld_path:
            script_content += f'export LD_LIBRARY_PATH="{ld_path}:$LD_LIBRARY_PATH"\n'

        script_content += "\n"
        script_content += self._generate_script_run_gnina()
        script_content += self._generate_missing_table_propagation()
        script_content += self.generate_completion_check_footer()
        return script_content

    def _generate_missing_table_propagation(self) -> str:
        """Generate script section to propagate missing.csv from upstream tools."""
        upstream_missing_path = self._get_upstream_missing_table_path(self.structures_input)

        if not upstream_missing_path:
            return ""

        upstream_folder = os.path.dirname(upstream_missing_path)

        return f"""
# Propagate missing table from upstream tools
echo "Checking for upstream missing structures..."
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

    def _generate_script_run_gnina(self) -> str:
        """Generate the GNINA execution part of the script."""
        return f"""echo "Starting GNINA docking"
echo "Structures: {len(self.structures_stream)} proteins"
echo "Compounds: {len(self.compounds_stream)} ligands"
echo "Output folder: {self.output_folder}"
echo "Exhaustiveness: {self.exhaustiveness}"
echo "Runs per conformer: {self.num_runs}"

python {self.helper_py} {self.config_json}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after GNINA execution."""
        upstream_missing_path = self._get_upstream_missing_table_path(self.structures_input)

        tables = {
            "docking_results": TableInfo(
                name="docking_results",
                path=self.docking_results_csv,
                columns=["id", "protein_id", "ligand_id", "conformer_id", "run",
                         "pose", "vina_score", "cnn_score", "cnn_affinity"],
                description="All accepted docked poses with Vina and CNN scores",
                count=0
            ),
            "conformer_ranking": TableInfo(
                name="conformer_ranking",
                path=self.conformer_ranking_csv,
                columns=["id", "protein_id", "ligand_id", "conformer_id",
                         "best_vina", "mean_vina", "std_vina", "best_cnn_score",
                         "pose_consistency", "conformer_energy",
                         "pseudo_binding_energy", "best_pose_file"],
                description="Aggregated conformer ranking with statistical analysis",
                count=0
            ),
        }

        if upstream_missing_path:
            tables["missing"] = TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "cause"],
                description="IDs removed by upstream tools with removal reason",
                count="variable"
            )

        best_poses_dir = os.path.join(self.output_folder, "best_poses")
        protein_ids = self.structures_stream.ids
        ligand_ids = self.compounds_stream.ids

        # Always predict conformer 0 as the representative best pose per
        # (protein, ligand) pair. When generate_conformers=True the actual
        # number of clustered conformers is only known at execution time, so
        # we do not attempt to enumerate them here.
        structure_ids = [
            f"{prot}_{lig}_conf0"
            for prot in protein_ids
            for lig in ligand_ids
        ]
        structure_files = [
            os.path.join(best_poses_dir, f"{sid}_best.pdb")
            for sid in structure_ids
        ]
        create_map_table(self.structures_map, structure_ids, files=structure_files)
        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=structure_files,
            map_table=self.structures_map,
            format="pdb",
        )

        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "gnina_params": {
                "center": self.center,
                "size": self.size,
                "autobox_ligand": self.autobox_ligand_path,
                "autobox_add": self.autobox_add,
                "exhaustiveness": self.exhaustiveness,
                "num_modes": self.num_modes,
                "num_runs": self.num_runs,
                "seed": self.seed,
                "cnn_scoring": self.cnn_scoring,
                "generate_conformers": self.generate_conformers,
                "num_conformers": self.num_conformers,
                "energy_window": self.energy_window,
                "conformer_rmsd": self.conformer_rmsd,
                "cnn_score_threshold": self.cnn_score_threshold,
                "rmsd_threshold": self.rmsd_threshold,
                "protonate": self.protonate,
                "pH": self.pH,
            }
        })
        return base_dict
