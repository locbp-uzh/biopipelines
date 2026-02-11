# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
GNINA docking and scoring configuration.

Combines AutoDock Vina scoring with CNN-based pose scoring for molecular docking.
Supports multi-conformer docking with statistical analysis for robust ranking.
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .config_manager import ConfigManager
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from config_manager import ConfigManager


class Gnina(BaseConfig):
    """
    GNINA molecular docking with CNN-based pose scoring.

    Performs protein-ligand docking using the GNINA binary (AutoDock Vina + CNN scoring).
    Supports multi-conformer docking with statistical analysis for robust ranking.
    """

    TOOL_NAME = "Gnina"

    # Lazy path descriptors
    gnina_binary = Path(lambda self: os.path.join(self.folders["Gnina"], "gnina"))
    helper_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_gnina.py"))
    docking_results_csv = Path(lambda self: os.path.join(self.output_folder, "docking_results.csv"))
    conformer_ranking_csv = Path(lambda self: os.path.join(self.output_folder, "conformer_ranking.csv"))
    config_json = Path(lambda self: os.path.join(self.output_folder, "gnina_config.json"))
    structures_json = Path(lambda self: os.path.join(self.output_folder, "structures_ds.json"))
    compounds_json = Path(lambda self: os.path.join(self.output_folder, "compounds_ds.json"))

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
                 conformer_energies: Optional[tuple] = None,
                 cnn_score_threshold: float = 0.5,
                 rmsd_threshold: float = 2.0,
                 protonate: bool = True,
                 pH: float = 7.4,
                 **kwargs):
        """
        Initialize GNINA docking configuration.

        Args:
            structures: Input protein structures as DataStream or StandardizedOutput
            compounds: Input ligands as DataStream or StandardizedOutput
            autobox_ligand: Reference ligand for automatic box definition
            center: Box center coordinates as "x,y,z" string
            size: Box dimensions as uniform float or "x,y,z" string
            autobox_add: Padding around autobox in Angstroms (default: 4.0)
            exhaustiveness: Search exhaustiveness (default: 32)
            num_modes: Poses per GNINA run (default: 9)
            num_runs: Independent runs per conformer (default: 5)
            seed: Random seed (default: 42)
            cnn_scoring: CNN scoring mode (default: "rescore")
            generate_conformers: Generate conformers with RDKit (default: False)
            num_conformers: Conformers to generate before filtering (default: 50)
            energy_window: Max kcal/mol above minimum to keep (default: 2.0)
            conformer_energies: Pre-computed energies as (table, "column") tuple
            cnn_score_threshold: Min CNNscore to keep a pose (default: 0.5)
            rmsd_threshold: RMSD cutoff for pose clustering in Angstroms (default: 2.0)
            protonate: Add hydrogens to protein (default: True)
            pH: Protonation pH (default: 7.4)
        """
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
                # Extract first file from structures stream
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

        # Must provide either center+size or autobox_ligand
        has_explicit_box = self.center is not None and self.size is not None
        has_autobox = self.autobox_ligand_path is not None
        if not has_explicit_box and not has_autobox:
            raise ValueError(
                "Must provide either center+size for explicit box "
                "or autobox_ligand for automatic box definition"
            )

        if self.exhaustiveness <= 0:
            raise ValueError("exhaustiveness must be positive")

        if self.num_modes <= 0:
            raise ValueError("num_modes must be positive")

        if self.num_runs <= 0:
            raise ValueError("num_runs must be positive")

        if self.cnn_score_threshold < 0 or self.cnn_score_threshold > 1:
            raise ValueError("cnn_score_threshold must be between 0 and 1")

        if self.rmsd_threshold <= 0:
            raise ValueError("rmsd_threshold must be positive")

        if self.energy_window <= 0:
            raise ValueError("energy_window must be positive")

        if self.num_conformers <= 0:
            raise ValueError("num_conformers must be positive")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and dependencies."""
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

        if self.protonate:
            config_lines.append(f"PROTONATE: pH {self.pH}")

        return config_lines

    def _write_config_json(self):
        """Write configuration and DataStreams to JSON files at pipeline time."""
        os.makedirs(self.output_folder, exist_ok=True)

        # Write structures DataStream
        with open(self.structures_json, 'w') as f:
            json.dump(self.structures_stream.to_dict(), f, indent=2)

        # Write compounds DataStream
        with open(self.compounds_json, 'w') as f:
            json.dump(self.compounds_stream.to_dict(), f, indent=2)

        # Build box config
        box_config = {}
        if self.center is not None and self.size is not None:
            box_config["center"] = self.center
            box_config["size"] = self.size
        if self.autobox_ligand_path is not None:
            box_config["autobox_ligand"] = self.autobox_ligand_path
            box_config["autobox_add"] = self.autobox_add

        # Write full config
        config = {
            "output_folder": self.output_folder,
            "gnina_binary": self.gnina_binary,
            "structures_json": self.structures_json,
            "compounds_json": self.compounds_json,
            "docking_results_csv": self.docking_results_csv,
            "conformer_ranking_csv": self.conformer_ranking_csv,
            "box": box_config,
            "exhaustiveness": self.exhaustiveness,
            "num_modes": self.num_modes,
            "num_runs": self.num_runs,
            "seed": self.seed,
            "cnn_scoring": self.cnn_scoring,
            "generate_conformers": self.generate_conformers,
            "num_conformers": self.num_conformers,
            "energy_window": self.energy_window,
            "conformer_energies_ref": self.conformer_energies_ref,
            "cnn_score_threshold": self.cnn_score_threshold,
            "rmsd_threshold": self.rmsd_threshold,
            "protonate": self.protonate,
            "pH": self.pH,
        }

        with open(self.config_json, 'w') as f:
            json.dump(config, f, indent=2)

    def _get_gnina_modules(self) -> List[str]:
        """Get CUDA modules for GNINA from config.yaml."""
        config_manager = ConfigManager()
        gnina_config = config_manager._config.get('gnina', {})
        return gnina_config.get('modules', [])

    def _get_gnina_ld_library_path(self) -> str:
        """Get LD_LIBRARY_PATH for GNINA from config.yaml."""
        config_manager = ConfigManager()
        gnina_config = config_manager._config.get('gnina', {})
        return gnina_config.get('ld_library_path', '')

    def generate_script(self, script_path: str) -> str:
        """Generate GNINA execution script."""
        # Write config and DataStreams to JSON at pipeline time
        self._write_config_json()

        script_content = "#!/bin/bash\n"
        script_content += "# GNINA docking execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        # Load CUDA modules from config
        modules = self._get_gnina_modules()
        if modules:
            script_content += "# Load CUDA modules for GNINA\n"
            script_content += "module load " + " ".join(modules) + "\n"

        # Set LD_LIBRARY_PATH from config
        ld_path = self._get_gnina_ld_library_path()
        if ld_path:
            script_content += f'export LD_LIBRARY_PATH="{ld_path}:$LD_LIBRARY_PATH"\n'

        script_content += "\n"
        script_content += self._generate_script_run_gnina()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_gnina(self) -> str:
        """Generate the GNINA execution part of the script."""
        script = f"""echo "Starting GNINA docking"
echo "Structures: {len(self.structures_stream)} proteins"
echo "Compounds: {len(self.compounds_stream)} ligands"
echo "Output folder: {self.output_folder}"
echo "Exhaustiveness: {self.exhaustiveness}"
echo "Runs per conformer: {self.num_runs}"

# Run GNINA docking pipeline
python {self.helper_py} {self.config_json}

"""
        return script

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after GNINA execution."""
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

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration with all GNINA parameters."""
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
                "cnn_score_threshold": self.cnn_score_threshold,
                "rmsd_threshold": self.rmsd_threshold,
                "protonate": self.protonate,
                "pH": self.pH,
            }
        })
        return base_dict
