# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
GNINA docking and scoring configuration.

GNINA combines AutoDock Vina search with a convolutional neural network (CNN)
for pose scoring, giving more accurate binding pose prediction than Vina alone.
Supports single-ligand docking, compound-library screening, and multi-conformer
docking with statistical analysis across independent runs.

Three modes (``mode=``):
  - "docking" (default): full Vina search + CNN rescore. The ligand 3-D pose is
    built from chemistry (``compounds``); the binding box is defined by explicit
    center+size, an autobox_ligand, or auto-detected crystal HETATM records.
  - "score": no search — score the ligand exactly where it already sits in the
    input complex (``gnina --score_only``). Each complex is paired with the one
    compound it carries (single-compound stream, else structures-map provenance);
    that ligand's HETATM records are extracted and bond-order-templated against
    its SMILES before scoring.
  - "minimize": local energy minimization of the in-pocket pose, then score
    (``gnina --minimize``). Same input shape as "score"; the refined pose is
    emitted. In both no-search modes the box is autoboxed on the scored ligand.

The docking binding box can be defined in three ways (checked in order):
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
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .config_manager import ConfigManager
    from .biopipelines_io import Resolve
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from config_manager import ConfigManager
    from biopipelines_io import Resolve


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
        - docking_summary: per-conformer aggregated statistics including
          mean scores across runs, pose consistency (fraction of runs that
          converge to the same pose), and pseudo_binding_energy (best Vina +
          conformer strain energy).
        - best_poses/: combined protein+ligand PDB files for each best pose.

    Usage:
        with Pipeline(...):
            Resources(gpu="A100", memory="16GB", time="8:00:00")
            protein = PDB("9RTM", convert="pdb")
            ligand = Ligand(smiles="CN(C)c1ccc2...", ids="TMR")
            docking = Gnina(structures=protein, compounds=ligand)
    """

    TOOL_NAME = "Gnina"
    TOOL_VERSION = "2.2"
    # gnina is a binary whose flag surface the wrapper types only part of.
    FORWARD_UNKNOWN_KWARGS = "argparse"
    ENV_NAME = "gnina"
    BACKEND = "gnina"

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
        env = cls._install_env(env_manager)
        try:
            from .config_manager import ConfigManager
        except ImportError:
            from config_manager import ConfigManager
        scheduler = ConfigManager().get_scheduler()

        repo_dir = folders.get("Gnina", "")
        binary = os.path.join(repo_dir, "gnina.1.3.2")

        # The GNINA wrapper's pipe script needs RDKit + OpenBabel (`obabel`) +
        # numpy/pandas. On the cluster these come from the shared `biopipelines`
        # env (config maps Gnina there). Colab has no `biopipelines` micromamba
        # env, so we create a dedicated `gnina` env there (config.colab.yaml maps
        # Gnina -> gnina). Only emitted on Colab to leave the cluster path intact.
        # On Colab a dedicated `gnina` env supplies the pipe script's
        # RDKit/OpenBabel deps (config.colab.yaml maps Gnina -> gnina); the
        # cluster reuses the shared `biopipelines` env and needs no env here.
        # Guard the env creation independently of the binary download so each
        # half is (re)done only when *its own* artefact is missing — on Colab
        # the binary may persist on Drive while the env is ephemeral, or vice
        # versa, and we must repair only the missing one.
        env_block = ""
        env_check = "true"  # cluster: nothing to verify beyond the binary.
        if scheduler == "colab":
            biopipelines = folders.get("biopipelines", "")
            install = cls._env_install_block(env, env_manager, biopipelines)
            env_check = cls._env_exists_check(env, env_manager)
            env_block = f"""# Create the {env} env (skip if it already exists).
if ! {env_check}; then
    {install}
else
    echo "{env} environment already exists, skipping creation."
fi
"""

        skip = "" if force_reinstall else f"""# Check if already installed
if [ -f "{binary}" ] && {env_check}; then
    echo "GNINA already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        return f"""echo "=== Installing GNINA ==="
{skip}{env_block}mkdir -p "{repo_dir}"
cd "{repo_dir}"
# Download the prebuilt binary only if absent (it may persist on Drive while
# the env was lost); -nc keeps a present file untouched.
if [ ! -f "{binary}" ]; then
    wget -nc https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2
    chmod +x gnina.1.3.2
fi
ln -sf gnina.1.3.2 gnina

# Verify installation
if [ -x "{binary}" ] && {env_check}; then
    touch "$INSTALL_SUCCESS"
    echo "=== GNINA installation complete ==="
else
    echo "ERROR: GNINA verification failed (binary missing/not executable or {env} env absent)"
    exit 1
fi
"""

    # Lazy path descriptors — canonical sub-layout.
    #   configuration/  — gnina_config.json + input DataStream JSONs.
    #   execution/      — prepared_proteins/, conformers/, docking/ (scratch).
    #   structures/     — best-pose PDBs + structures_map.csv.
    #   tables/         — docking_results, docking_summary, missing.
    gnina_binary = Path(lambda self: os.path.join(self.folders["Gnina"], "gnina"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_gnina.py"))
    docking_results_csv = Path(lambda self: self.table_path("docking_results"))
    docking_summary_csv = Path(lambda self: self.table_path("docking_summary"))
    scores_csv = Path(lambda self: self.table_path("scores"))
    config_json = Path(lambda self: self.configuration_path("gnina_config.json"))
    structures_json = Path(lambda self: self.configuration_path("structures_ds.json"))
    compounds_json = Path(lambda self: self.configuration_path("compounds_ds.json"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    autobox_ds_json = Path(lambda self: self.configuration_path("autobox_ligand_ds.json"))

    MODES = ("docking", "score", "minimize")

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 compounds: Union[DataStream, StandardizedOutput],
                 mode: str = "docking",
                 scoring: str = "vina",
                 autobox_ligand: Union[DataStream, StandardizedOutput, str, None] = None,
                 center: Optional[str] = None,
                 size: Union[float, str, None] = None,
                 autobox_add: float = 4.0,
                 exhaustiveness: int = 8,
                 num_modes: int = 9,
                 num_runs: int = 1,
                 seed: int = 42,
                 cnn_scoring: str = "rescore",
                 generate_conformers: bool = False,
                 num_conformers: int = 50,
                 energy_window: float = 2.0,
                 conformer_rmsd: float = 1.0,
                 conformer_energies: Optional[tuple] = None,
                 cnn_score_threshold: float = 0.0,
                 rmsd_threshold: float = 2.0,
                 protonate: bool = True,
                 pH: float = 7.4,
                 dock_timeout: Optional[int] = 1800,
                 **kwargs):
        """
        Initialize GNINA docking configuration.

        Args:
            structures: Input protein structures as DataStream or StandardizedOutput.
                In mode="score"/"minimize" these are complexes that already carry
                the ligand as HETATM records (e.g. a Boltz2 co-fold or a crystal
                complex).
            compounds: Input ligands as DataStream or StandardizedOutput. Accepts
                output from Ligand(), CompoundLibrary(), or any tool that produces
                a compounds stream (e.g. BoltzGen). In mode="score"/"minimize"
                each complex is scored against the one compound it carries: a
                single-compound stream broadcasts to every complex, otherwise the
                complex→compound link is read from the structures-map provenance
                (compounds.id). The paired compound's code locates the HETATM
                residue and its SMILES templates that pose's bond orders, so each
                compound must carry SMILES; a code-only Ligand is rejected.
            mode: One of "docking" (default, full Vina search + CNN rescore),
                "score" (no search; score the in-pocket pose as-is via
                --score_only), or "minimize" (local minimization of the in-pocket
                pose then score, via --minimize). The docking-only parameters
                (exhaustiveness, num_modes, num_runs, generate_conformers and the
                conformer knobs, rmsd_threshold, explicit center/size/
                autobox_ligand) must stay at their defaults when mode is not
                "docking".
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
                conformational space but take longer (default 8, matching the
                upstream GNINA CLI default).
            num_modes: Number of docked poses output per GNINA run (default 9,
                matching the upstream GNINA CLI default).
            num_runs: Number of independent docking runs per conformer with
                incremented seeds. More runs improve pose consistency statistics
                (default 1).
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
                Poses below this threshold are discarded (default 0.0, keep
                everything). GNINA itself has no such filter — it sorts poses by
                CNNscore and returns them all — so a non-zero value is a choice
                this wrapper adds, not upstream behaviour. CNNscore answers "is
                this pose near the crystal pose", which is uncorrelated with
                binding strength (measured rho +0.025 over 578 de novo designs)
                and is trained on natural complexes, so designed pockets and
                synthetic ligands score low across the board. Set it only to
                screen for pose confidence; leave it at 0 when every input needs
                a score.
            rmsd_threshold: RMSD cutoff (Å) for pose consistency clustering
                across runs. Larger clusters indicate more reproducible poses
                (default 2.0).
            protonate: If True, add hydrogens to the protein with OpenBabel at
                the specified pH before docking (default True). Requires obabel
                in PATH.
            pH: Protonation pH for OpenBabel (default 7.4).
            dock_timeout: Per-run wall-clock cap in seconds for a single GNINA
                docking invocation (default 1800 = 30 min). A flexible ligand
                (e.g. an aminoglycoside) at high exhaustiveness can occasionally
                produce a Vina search that never terminates and would otherwise
                hang the whole job; a run that exceeds this is killed, warned,
                and treated as a failed run (the other runs/poses still count).
                Set generously above a normal dock — scale it up if you raise
                exhaustiveness. None or 0 disables the cap.

        Output:
            mode="docking":
                Streams: structures (.pdb) — best-pose complexes.
                Tables:
                    docking_results: id | structures.id | compounds.id | conformer_id | run | pose | vina_score | cnn_score | cnn_affinity
                    docking_summary: id | structures.id | compounds.id | conformer_id | best_vina | mean_vina | std_vina | best_cnn_score | mean_cnn_affinity | std_cnn_affinity | pose_consistency | conformer_energy | pseudo_binding_energy | best_pose_file
                    missing: id | removed_by | kind | cause
            mode="score"/"minimize":
                Streams: structures (.pdb) — the scored complex ("score" re-emits
                    the input pose unchanged; "minimize" emits the refined pose).
                Tables:
                    scores: id | structures.id | compounds.id | vina_affinity | cnn_score | cnn_affinity | cnn_vs | cnn_affinity_variance
                    missing: id | removed_by | kind | cause
        """
        if "backend" in kwargs:
            raise ValueError(
                "backend was removed; use the Gnina class for GNINA and the "
                "Vina class for AutoDock Vina.")

        self.mode = mode
        self.backend = self.BACKEND
        self.scoring = scoring

        # Keep original input for upstream missing table detection
        self.structures_input = structures

        # Resolve structures input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        # Keep original input for upstream missing table detection
        self.compounds_input = compounds

        # Resolve compounds input to DataStream
        if isinstance(compounds, StandardizedOutput):
            self.compounds_stream: DataStream = compounds.streams.compounds
        elif isinstance(compounds, DataStream):
            self.compounds_stream = compounds
        else:
            raise ValueError(f"compounds must be DataStream or StandardizedOutput, got {type(compounds)}")

        # Resolve autobox_ligand — store stream for runtime resolution
        self.autobox_ligand_stream = None
        self.autobox_ligand_id = None
        self.autobox_ligand_path = None  # Only set for string paths (no runtime resolution needed)
        if autobox_ligand is not None:
            # Dispatch on type before asking anything about the stream: hanging the
            # str case and the type error off a length test made an empty DataStream
            # report "must be DataStream ... got DataStream".
            if isinstance(autobox_ligand, StandardizedOutput):
                self.autobox_ligand_stream = autobox_ligand.streams.structures
            elif isinstance(autobox_ligand, DataStream):
                self.autobox_ligand_stream = autobox_ligand
            elif isinstance(autobox_ligand, str):
                self.autobox_ligand_path = autobox_ligand
            else:
                raise ValueError(f"autobox_ligand must be DataStream, StandardizedOutput, or str, got {type(autobox_ligand)}")
            # Only a single-item stream is a global reference; a per-structure
            # stream is resolved by id at runtime.
            if self.autobox_ligand_stream is not None and len(self.autobox_ligand_stream) == 1:
                self.autobox_ligand_id = self.autobox_ligand_stream.ids[0]

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
        self.dock_timeout = dock_timeout
        self.rmsd_threshold = rmsd_threshold

        # Protonation
        self.protonate = protonate
        self.pH = pH

        # Store conformer_energies table reference if provided
        self.conformer_energies_ref = None
        if conformer_energies is not None:
            self.conformer_energies_ref = conformer_energies

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate GNINA parameters."""
        if self.mode not in self.MODES:
            raise ValueError(f"mode must be one of {self.MODES}, got {self.mode!r}")

        if self.backend == "vina":
            self._validate_vina_backend()
        elif self.scoring != "vina":
            raise ValueError(
                "scoring only applies to Vina; Gnina selects its scoring "
                "function with cnn_scoring instead."
            )

        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        if not self.structures_stream.has_only_formats("pdb"):
            raise ValueError(
                f"GNINA requires PDB format structures, got '{self.structures_stream.format}'. "
                f"Use convert='pdb' in PDB() or RCSB() to ensure PDB format output."
            )

        if not self.compounds_stream or len(self.compounds_stream) == 0:
            raise ValueError("compounds parameter is required and must not be empty")

        if self.mode != "docking":
            self._validate_no_search_mode()
        else:
            has_explicit_box = self.center is not None and self.size is not None
            has_autobox = self.autobox_ligand_path is not None or self.autobox_ligand_stream is not None
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

    VINA_SCORING = ("vina", "vinardo", "ad4")

    def _validate_vina_backend(self):
        """Reject GNINA-only parameters on the Vina tool.

        Vina has no CNN, so these knobs cannot be honored. Raising keeps a run
        that set them from looking like it applied them.
        """
        if self.scoring not in self.VINA_SCORING:
            raise ValueError(
                f"scoring must be one of {self.VINA_SCORING}, got {self.scoring!r}")

        # Defaults come from Gnina's own signature: hard-coding them here meant a
        # change to one (85d3134 moved cnn_score_threshold to 0.0) turned the
        # inherited default into an "offender" and Vina stopped constructing at all.
        import inspect
        inherited = inspect.signature(Gnina.__init__).parameters
        cnn_only = ("cnn_scoring", "cnn_score_threshold")
        offenders = [name for name in cnn_only
                     if getattr(self, name) != inherited[name].default]
        if offenders:
            raise ValueError(
                f"Vina has no CNN scoring; remove "
                f"{', '.join(sorted(offenders))} (these only apply to Gnina)."
            )

    def _validate_no_search_mode(self):
        """Reject docking-only parameters set non-default in score/minimize mode.

        These knobs drive the Vina search / conformer pipeline, which the
        no-search modes skip entirely. Raising (rather than silently ignoring)
        keeps a misconfigured score run from looking like it honored them.
        """
        docking_only = {
            "center": (self.center, None),
            "size": (self.size, None),
            "autobox_ligand": (self.autobox_ligand_stream or self.autobox_ligand_path, None),
            "exhaustiveness": (self.exhaustiveness, 8),
            "num_modes": (self.num_modes, 9),
            "num_runs": (self.num_runs, 1),
            "generate_conformers": (self.generate_conformers, False),
            "num_conformers": (self.num_conformers, 50),
            "energy_window": (self.energy_window, 2.0),
            "conformer_rmsd": (self.conformer_rmsd, 1.0),
            "conformer_energies": (self.conformer_energies, None),
            "rmsd_threshold": (self.rmsd_threshold, 2.0),
        }
        offenders = [name for name, (value, default) in docking_only.items()
                     if value != default]
        if offenders:
            raise ValueError(
                f"mode={self.mode!r} ignores docking-only parameters; remove "
                f"{', '.join(sorted(offenders))} (these only apply to mode='docking')."
            )

    @property
    def docking_binary(self) -> str:
        """The binary to invoke.

        Vina is installed by conda into the tool env, so it resolves off PATH and
        has no repo folder; GNINA is a standalone download under folders["Gnina"].
        """
        if self.backend == "vina":
            return "vina"
        return self.gnina_binary

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files and folder paths."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get GNINA configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"BACKEND: {self.backend}",
            f"MODE: {self.mode}",
            f"STRUCTURES: {len(self.structures_stream)} proteins",
            f"COMPOUNDS: {len(self.compounds_stream)} ligands",
            f"SEED: {self.seed}",
        ])
        config_lines.append(
            f"SCORING: {self.scoring}" if self.backend == "vina"
            else f"CNN SCORING: {self.cnn_scoring}"
        )
        if self.mode == "docking":
            config_lines.extend([
                f"EXHAUSTIVENESS: {self.exhaustiveness}",
                f"NUM RUNS: {self.num_runs}",
                f"NUM MODES: {self.num_modes}",
            ])
            if self.backend != "vina":
                config_lines.append(f"CNN SCORE THRESHOLD: {self.cnn_score_threshold}")
        else:
            config_lines.append(f"AUTOBOX ADD: {self.autobox_add}")
        if self.center is not None:
            config_lines.append(f"BOX CENTER: {self.center}")
            config_lines.append(f"BOX SIZE: {self.size}")
        if self.autobox_ligand_stream is not None:
            config_lines.append(f"AUTOBOX LIGAND: {self.autobox_ligand_id}")
            config_lines.append(f"AUTOBOX ADD: {self.autobox_add}")
        elif self.autobox_ligand_path is not None:
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
        """Write configuration and DataStreams to JSON files at pipeline time.

        Scratch folders (prepared_proteins/, conformers/, docking/) all go
        under execution/; best-pose PDBs land in the structures/ stream
        folder via the explicit best_poses_dir key.
        """
        self.structures_stream.save_json(self.structures_json)
        self.compounds_stream.save_json(self.compounds_json)

        box_config = {"autobox_add": self.autobox_add}
        if self.center is not None and self.size is not None:
            box_config["center"] = self.center
            box_config["size"] = self.size
        if self.autobox_ligand_stream is not None:
            # Serialize autobox ligand DataStream for runtime resolution
            self.autobox_ligand_stream.save_json(self.autobox_ds_json)
            if len(self.autobox_ligand_stream) == 1:
                # One shared reference: the generated script rewrites this sentinel
                # with the resolved path before the runtime reads it.
                box_config["autobox_ligand"] = "__RESOLVE_AUTOBOX_LIGAND__"
            else:
                # Per-structure references are looked up by receptor id at runtime,
                # so there is no single path to name here -- leaving the sentinel
                # would make every receptor try to open it literally.
                box_config["autobox_ligand_ds"] = self.autobox_ds_json
        elif self.autobox_ligand_path is not None:
            box_config["autobox_ligand"] = self.autobox_ligand_path

        upstream_missing = self._collect_upstream_missing_paths(self.structures_input, self.compounds_input)

        config = {
            "mode": self.mode,
            "backend": self.backend,
            "scoring": self.scoring,
            "output_folder": self.execution_folder,
            "best_poses_dir": self.stream_folder("structures"),
            # Vina comes from its conda env's PATH; GNINA is a standalone
            # download under its own repo folder.
            "gnina_binary": self.docking_binary,
            "structures_json": self.structures_json,
            "compounds_json": self.compounds_json,
            "docking_results_csv": self.docking_results_csv,
            "docking_summary_csv": self.docking_summary_csv,
            "scores_csv": self.scores_csv,
            "structures_map": self.structures_map,
            "upstream_missing": upstream_missing,
            # Score/minimize always emit the missing table to record local
            # failures (no HETATM, SMILES mismatch); docking only on upstream.
            "missing_csv": self.missing_csv if (upstream_missing or self.mode != "docking") else None,
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
            "dock_timeout": self.dock_timeout,
            "rmsd_threshold": self.rmsd_threshold,
            "protonate": self.protonate,
            "pH": self.pH,
            # Argv tokens, not a shell string: the pipe script appends them to the list it hands to
            # subprocess.run with no shell, so nothing re-parses them.
            "extra_args": self.extra_args_tokens(),
        }

        with open(self.config_json, 'w') as f:
            json.dump(config, f, indent=2)

    def _get_gnina_overrides(self) -> Dict[str, Any]:
        """Read the per-tool overrides for GNINA from the active config.

        Lives at ``tool_overrides.gnina`` in the YAML — sites that need
        to pin CUDA modules or set LD_LIBRARY_PATH for GNINA fill in
        this block. Returns an empty dict when not configured (other
        sites don't need either knob).
        """
        config_manager = ConfigManager()
        overrides = config_manager._config.get('tool_overrides', {}) or {}
        # Keyed on 'gnina' because that is what sites configure, but the settings
        # are CUDA modules and library paths for the GPU binary. Vina is CPU-only
        # conda-forge software in its own env, so loading GNINA's CUDA stack into
        # its step is at best noise and at worst a library conflict.
        if self.BACKEND != "gnina":
            return {}
        return overrides.get('gnina', {}) or {}

    def _get_gnina_modules(self) -> List[str]:
        """CUDA modules to load before invoking gnina (tool_overrides.gnina.modules)."""
        modules = self._get_gnina_overrides().get('modules', []) or []
        for i, mod in enumerate(modules):
            _validate_freeform_string(f"tool_overrides.gnina.modules[{i}]", mod)
        return modules

    def _get_gnina_ld_library_path(self) -> str:
        """LD_LIBRARY_PATH to prepend (tool_overrides.gnina.ld_library_path)."""
        ld_path = self._get_gnina_overrides().get('ld_library_path', '') or ''
        _validate_freeform_string("tool_overrides.gnina.ld_library_path", ld_path)
        return ld_path

    def generate_script(self, script_path: str) -> str:
        """Generate GNINA execution script."""
        self._write_config_json()

        script_content = "#!/bin/bash\n"
        script_content += "# GNINA docking execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.warn_container_unsupported()

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
        """Propagate the union of upstream `missing` manifests into tables/missing.csv."""
        return self.generate_missing_propagation(
            self.structures_input, self.compounds_input, missing_csv=self.missing_csv
        )

    def _generate_script_resolve_autobox(self) -> str:
        """Generate bash snippet to resolve autobox ligand at runtime."""
        if self.autobox_ligand_stream is None:
            return ""

        # A stream of one is a single shared reference and can be resolved here.
        # A stream with one entry per structure must be resolved per receptor at
        # runtime -- resolving index 0 and reusing that path silently boxed every
        # receptor on the first structure's ligand.
        if len(self.autobox_ligand_stream) != 1:
            return f"""# Per-structure autobox: the runtime resolves one reference per receptor id
jq --arg ds "{self.autobox_ds_json}" '.box.autobox_ligand_ds = $ds' "{self.config_json}" > "{self.config_json}.tmp" && mv "{self.config_json}.tmp" "{self.config_json}"

"""

        return f"""# Resolve autobox ligand file at runtime (valid_set skips filtered ids)
AUTOBOX_LIGAND_ID={Resolve.stream_ids(self.autobox_ds_json, index=0, valid_set=True)}
AUTOBOX_LIGAND={Resolve.stream_item(self.autobox_ds_json, '$AUTOBOX_LIGAND_ID')}
echo "Resolved autobox ligand: $AUTOBOX_LIGAND"
jq --arg path "$AUTOBOX_LIGAND" '.box.autobox_ligand = $path' "{self.config_json}" > "{self.config_json}.tmp" && mv "{self.config_json}.tmp" "{self.config_json}"

"""

    def _generate_script_run_gnina(self) -> str:
        """Generate the GNINA execution part of the script."""
        return f"""{self._generate_script_resolve_autobox()}echo "Starting GNINA docking"
echo "Structures: {len(self.structures_stream)} proteins"
echo "Compounds: {len(self.compounds_stream)} ligands"
echo "Output folder: {self.output_folder}"
echo "Exhaustiveness: {self.exhaustiveness}"
echo "Runs per conformer: {self.num_runs}"

python {self.helper_py} {self.config_json}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after GNINA execution."""
        # Vina has no CNN, so those columns are absent rather than blank.
        is_vina = self.backend == "vina"
        cnn_pose_cols = [] if is_vina else ["cnn_score", "cnn_affinity"]
        cnn_summary_cols = [] if is_vina else [
            "best_cnn_score", "mean_cnn_affinity", "std_cnn_affinity"]
        cnn_score_cols = [] if is_vina else [
            "cnn_score", "cnn_affinity", "cnn_vs", "cnn_affinity_variance"]
        engine = "Vina" if is_vina else "GNINA"

        if self.mode == "docking":
            tables = {
                "docking_results": TableInfo(
                    name="docking_results",
                    path=self.docking_results_csv,
                    columns=["id", "structures.id", "compounds.id", "conformer_id",
                             "run", "pose", "vina_score"] + cnn_pose_cols,
                    description=f"All accepted docked poses with {engine} scores"
                ),
                "docking_summary": TableInfo(
                    name="docking_summary",
                    path=self.docking_summary_csv,
                    columns=["id", "structures.id", "compounds.id", "conformer_id",
                             "best_vina", "mean_vina", "std_vina"] + cnn_summary_cols
                            + ["pose_consistency", "conformer_energy",
                               "pseudo_binding_energy", "best_pose_file"],
                    description="Per-conformer aggregated docking statistics across runs"
                ),
            }
        else:
            tables = {
                "scores": TableInfo(
                    name="scores",
                    path=self.scores_csv,
                    columns=["id", "structures.id", "compounds.id",
                             "vina_affinity"] + cnn_score_cols,
                    description=f"{engine} {self.mode} scores for the in-pocket ligand pose"
                ),
            }

        if self.mode != "docking" or self._collect_upstream_missing_paths(
                self.structures_input, self.compounds_input):
            tables["missing"] = self.missing_table_info(self.missing_csv)

        best_poses_dir = self.stream_folder("structures")
        protein_ids = list(self.structures_stream.ids)

        if self.mode == "docking":
            # Docking poses every ligand into every protein -> prot × lig.
            ligand_ids = list(self.compounds_stream.ids)
            structure_ids = [
                f"{prot}_{lig}"
                for prot in protein_ids
                for lig in ligand_ids
            ]
        else:
            # Score/minimize scores the one ligand each complex already carries,
            # so the output keys on the complex's own id (one row per complex).
            structure_ids = protein_ids
        suffix = "_best.pdb" if self.mode == "docking" else f"_{self.mode}.pdb"
        structure_files = [os.path.join(best_poses_dir, f"<id>{suffix}")]

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
                "mode": self.mode,
                "backend": self.backend,
                "scoring": self.scoring,
                "center": self.center,
                "size": self.size,
                "autobox_ligand": self.autobox_ligand_id or self.autobox_ligand_path,
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


class Vina(Gnina):
    """AutoDock Vina docking — the same pipeline as Gnina without CNN scoring.

    Vina runs on CPU and has aarch64 builds, so it is the docking route on
    architectures where GNINA's prebuilt x86-64 binary cannot run. It keeps the
    box definition, conformer generation, multi-run pose-consistency analysis and
    the three modes; the CNN parameters are absent because Vina has no CNN, and
    the score tables carry no ``cnn_*`` columns.

    Vina reads PDBQT only, so receptors and ligands are converted with OpenBabel
    on the way in and the poses converted back to SDF on the way out. Receptor
    preparation via OpenBabel is the pragmatic route rather than the
    Meeko/ADFR-preferred one — adequate for ranking, but note it if you need
    publication-grade absolute affinities.

    Usage:
        with Pipeline(...):
            Resources(cpus=32, memory="16GB", time="4:00:00")
            protein = PDB("9RTM", convert="pdb")
            ligand = Ligand(smiles="CN(C)c1ccc2...", ids="TMR")
            docking = Vina(structures=protein, compounds=ligand)
    """

    TOOL_NAME = "Vina"
    TOOL_VERSION = "2.2"
    ENV_NAME = "vina"
    BACKEND = "vina"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Install AutoDock Vina from conda-forge.

        Unlike GNINA there is no binary to download: conda-forge ships `vina`
        for every platform including linux-aarch64, alongside the RDKit and
        OpenBabel the pipe script needs.
        """
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check(env, env_manager)
        install = cls._env_install_block(env, env_manager, biopipelines)

        run = cls._env_run(env, env_manager)

        # The skip path tests the binary, not just the env: a venv shim can exist
        # while the symlinks below are missing.
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check} && {run}vina --version >/dev/null 2>&1; then
    echo "Vina already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        return f"""echo "=== Installing AutoDock Vina ==="
{skip}{install}
{cls._link_conda_binaries(env, env_manager, ("vina", "vina_split", "obabel"))}
# Verify installation: the binary must run and OpenBabel must be present for the
# PDBQT conversions the wrapper depends on.
if {run}vina --version >/dev/null 2>&1 && {run}obabel -V >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== Vina installation complete ==="
else
    echo "ERROR: Vina verification failed (vina or obabel not runnable in the env)"
    exit 1
fi
"""

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 compounds: Union[DataStream, StandardizedOutput],
                 mode: str = "docking",
                 scoring: str = "vina",
                 autobox_ligand: Union[DataStream, StandardizedOutput, str, None] = None,
                 center: Optional[str] = None,
                 size: Union[float, str, None] = None,
                 autobox_add: float = 4.0,
                 exhaustiveness: int = 8,
                 num_modes: int = 9,
                 num_runs: int = 1,
                 seed: int = 42,
                 generate_conformers: bool = False,
                 num_conformers: int = 50,
                 energy_window: float = 2.0,
                 conformer_rmsd: float = 1.0,
                 conformer_energies: Optional[tuple] = None,
                 rmsd_threshold: float = 2.0,
                 protonate: bool = True,
                 pH: float = 7.4,
                 dock_timeout: Optional[int] = 1800,
                 **kwargs):
        """AutoDock Vina docking. See Gnina for the shared argument details.

        Args:
            scoring: Vina scoring function — "vina" (default), "vinardo", or
                "ad4". Replaces Gnina's cnn_scoring, which has no counterpart.

        Output:
            mode="docking":
                Streams: structures (.pdb) — best-pose complexes.
                Tables:
                    docking_results: id | structures.id | compounds.id | conformer_id | run | pose | vina_score
                    docking_summary: id | structures.id | compounds.id | conformer_id | best_vina | mean_vina | std_vina | pose_consistency | conformer_energy | pseudo_binding_energy | best_pose_file
                    missing: id | removed_by | kind | cause
            mode="score"/"minimize":
                Streams: structures (.pdb) — the scored complex.
                Tables:
                    scores: id | structures.id | compounds.id | vina_affinity
                    missing: id | removed_by | kind | cause
        """
        super().__init__(
            structures=structures,
            compounds=compounds,
            mode=mode,
            scoring=scoring,
            autobox_ligand=autobox_ligand,
            center=center,
            size=size,
            autobox_add=autobox_add,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes,
            num_runs=num_runs,
            seed=seed,
            generate_conformers=generate_conformers,
            num_conformers=num_conformers,
            energy_window=energy_window,
            conformer_rmsd=conformer_rmsd,
            conformer_energies=conformer_energies,
            rmsd_threshold=rmsd_threshold,
            protonate=protonate,
            pH=pH,
            dock_timeout=dock_timeout,
            **kwargs,
        )
