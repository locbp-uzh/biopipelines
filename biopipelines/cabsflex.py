# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
CABS-Flex tool for fast simulation of protein structure flexibility.

CABS-Flex performs coarse-grained Monte Carlo simulations to model protein
flexibility, producing an ensemble of conformational models and per-residue
RMSF (Root Mean Square Fluctuation) profiles.

Reference: Kurcinski et al. (2019) Bioinformatics, 35(4):694-695
Repository: https://bitbucket.org/lcbio/cabsflex
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, IndexedTableContainer
    from .file_paths import Path
    from .datastream import DataStream, create_map_table
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, IndexedTableContainer
    from file_paths import Path
    from datastream import DataStream, create_map_table


class CABSflex(BaseConfig):
    """
    CABS-Flex: fast simulation of protein structure flexibility.

    Runs coarse-grained Monte Carlo simulations on input protein structures
    to generate conformational ensembles and flexibility profiles (RMSF).
    Optionally rebuilds final models to all-atom representation using MODELLER.

    WARNING: MODELLER requires a license key. Obtain one (free for academics)
    at https://salilab.org/modeller/registration.html and set the environment
    variable KEY_MODELLER before running.
    """

    TOOL_NAME = "CABSflex"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        if env_manager == "pip":
            skip = "" if force_reinstall else """# Check if already installed
if python -c "import CABS" 2>/dev/null; then
    echo "CABS-Flex already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
            return f"""echo "=== Installing CABS-Flex (pip) ==="
{skip}echo "WARNING: MODELLER requires a license key."
echo "Get one at https://salilab.org/modeller/registration.html"
echo "Then set: export KEY_MODELLER=your_key"
pip install -q modeller cabs dssp

echo "=== CABS-Flex installation complete ==="
"""
        skip = "" if force_reinstall else """# Check if already installed
if conda list -n CABSflex cabs 2>/dev/null | grep -q cabs; then
    echo "CABS-Flex already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        return f"""echo "=== Installing CABS-Flex ==="
{skip}echo "WARNING: MODELLER requires a license key."
echo "Get one at https://salilab.org/modeller/registration.html"
echo "Then set: export KEY_MODELLER=your_key"
{env_manager} create -y -n CABSflex python=2.7
{env_manager} install -y -n CABSflex -c salilab modeller dssp
{env_manager} install -y -n CABSflex -c lcbio cabs

echo "=== CABS-Flex installation complete ==="
"""

    # Lazy path descriptors
    structures_ds_json = Path(lambda self: os.path.join(self.output_folder, "structures.json"))
    structures_map = Path(lambda self: os.path.join(self.output_folder, "structures_map.csv"))
    rmsf_all_csv = Path(lambda self: os.path.join(self.output_folder, "rmsf_all.csv"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_cabsflex.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 num_models: int = 10,
                 mc_cycles: int = 50,
                 mc_steps: int = 50,
                 mc_annealing: int = 20,
                 temperature: Optional[str] = None,
                 flexibility: Optional[str] = None,
                 filtering_count: int = 1000,
                 aa_rebuild: bool = False,
                 restraints: str = "ss2",
                 restraints_gap: int = 3,
                 restraints_min: float = 3.8,
                 restraints_max: float = 8.0,
                 restraints_reduce: Optional[float] = None,
                 weighted_fit: Optional[str] = None,
                 pdb_output: str = "M",
                 max_parallel: int = 1,
                 **kwargs):
        """
        Initialize CABS-Flex configuration.

        Args:
            structures: Input protein structures as DataStream or StandardizedOutput
            num_models: Number of cluster medoids / final models (default: 10)
            mc_cycles: Monte Carlo cycles; total snapshots = mc_annealing × mc_cycles (default: 50)
            mc_steps: Monte Carlo cycles between trajectory frames (default: 50)
            mc_annealing: Temperature annealing cycles (default: 20)
            temperature: Temperature range as "TINIT TFINAL" (default: "1.4 1.4")
            flexibility: Residue flexibility: float (0=flexible, 1=stiff), 'bf', 'bfi', 'bfg', or filename
            filtering_count: Number of low-energy models for clustering (default: 1000)
            aa_rebuild: Rebuild to all-atom with MODELLER (default: False). Requires MODELLER license.
            restraints: Restraint mode: 'all', 'ss1', 'ss2' (default: 'ss2')
            restraints_gap: Min gap along chain for restraints (default: 3)
            restraints_min: Min distance in Angstroms for restraints (default: 3.8)
            restraints_max: Max distance in Angstroms for restraints (default: 8.0)
            restraints_reduce: Reduce restraints by factor in [0,1] to speed up computation. None = no reduction.
            weighted_fit: Fit method: 'gauss', 'flex', 'ss', 'off', or filename. None = CABSflex default.
            pdb_output: Which structures to save: 'A' (all), 'R' (replicas), 'F' (filtered),
                        'C' (clusters), 'M' (models), 'S' (starting), 'N' (none). Combinable, e.g. 'RM'. (default: 'M')
            max_parallel: Max structures to run concurrently (default: 1 = sequential)
            **kwargs: Additional parameters

        Output:
            Streams:
                structures: PDB ensemble models (num_models per input structure)
                images: SVG plots (RMSF, RMSD, energy) per input structure
            Tables:
                rmsf_all: id | chain | resi | rmsf  (merged, all structures)
                rmsf.<input_id>: id | chain | resi | rmsf  (per input structure)
                structures: id | file | structures.id
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        self.num_models = num_models
        self.mc_cycles = mc_cycles
        self.mc_steps = mc_steps
        self.mc_annealing = mc_annealing
        self.temperature = temperature
        self.flexibility = flexibility
        self.filtering_count = filtering_count
        self.aa_rebuild = aa_rebuild
        self.restraints = restraints
        self.restraints_gap = restraints_gap
        self.restraints_min = restraints_min
        self.restraints_max = restraints_max
        self.restraints_reduce = restraints_reduce
        self.weighted_fit = weighted_fit
        self.pdb_output = pdb_output
        self.max_parallel = max_parallel

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate CABS-Flex parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if self.num_models < 1:
            raise ValueError("num_models must be >= 1")
        if self.mc_cycles < 1:
            raise ValueError("mc_cycles must be >= 1")
        if self.mc_steps < 1:
            raise ValueError("mc_steps must be >= 1")
        if self.mc_annealing < 1:
            raise ValueError("mc_annealing must be >= 1")
        if self.restraints not in ("all", "ss1", "ss2"):
            raise ValueError(f"restraints must be 'all', 'ss1', or 'ss2', got: {self.restraints}")
        if self.restraints_reduce is not None and not (0 <= self.restraints_reduce <= 1):
            raise ValueError(f"restraints_reduce must be in [0, 1], got: {self.restraints_reduce}")
        if self.max_parallel < 1:
            raise ValueError(f"max_parallel must be >= 1, got: {self.max_parallel}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"INPUT STRUCTURES: {len(self.structures_stream)} files",
            f"NUM MODELS: {self.num_models}",
            f"MC CYCLES: {self.mc_cycles}",
            f"MC ANNEALING: {self.mc_annealing}",
            f"AA REBUILD: {self.aa_rebuild}",
            f"RESTRAINTS: {self.restraints}",
            f"RESTRAINTS REDUCE: {self.restraints_reduce}",
            f"WEIGHTED FIT: {self.weighted_fit}",
            f"MAX PARALLEL: {self.max_parallel}"
        ])
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate CABS-Flex execution script."""
        import json

        # Serialize structures DataStream to JSON
        with open(self.structures_ds_json, 'w') as f:
            json.dump(self.structures_stream.to_dict(), f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# CABS-Flex execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_cabsflex()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_cabsflex(self) -> str:
        """Generate the CABS-Flex execution part of the script.

        CABSflex runs under Python 2.7, so execution is done in pure bash.
        Post-processing (file copying, RMSF parsing) runs under biopipelines
        (Python 3) via the helper script.
        """
        # Only pass flags that differ from CABSflex defaults
        flags = []
        if self.mc_cycles != 50:
            flags.append(f"-y {self.mc_cycles}")
        if self.mc_steps != 50:
            flags.append(f"-s {self.mc_steps}")
        if self.mc_annealing != 20:
            flags.append(f"-a {self.mc_annealing}")
        if self.num_models != 10:
            flags.append(f"-k {self.num_models}")
        if self.filtering_count != 1000:
            flags.append(f"-n {self.filtering_count}")
        if self.restraints != "ss2" or self.restraints_gap != 3 or self.restraints_min != 3.8 or self.restraints_max != 8.0:
            flags.append(f"-g {self.restraints} {self.restraints_gap} {self.restraints_min} {self.restraints_max}")
        if self.temperature:
            flags.append(f"-t {self.temperature}")
        if self.flexibility is not None:
            flags.append(f"-f {self.flexibility}")
        if self.weighted_fit is not None:
            flags.append(f"--weighted-fit {self.weighted_fit}")
        if self.restraints_reduce is not None:
            flags.append(f"--protein-restraints-reduce {self.restraints_reduce}")
        if self.aa_rebuild:
            flags.append("-A")
        if self.pdb_output != "A":
            flags.append(f"-o {self.pdb_output}")

        flags_str = " ".join(flags)

        # Generate commands for each structure (runs under Py2.7 CABSflex env)
        n_structures = len(self.structures_stream.ids)
        run_parallel = self.max_parallel > 1 and n_structures > 1

        cabsflex_cmds = ""
        if run_parallel:
            cabsflex_cmds += f'\necho "Running {n_structures} structures in parallel (max {self.max_parallel} concurrent)"\n'
            cabsflex_cmds += "PIDS=()\nFAILED=0\n"

        for sid, sfile in zip(self.structures_stream.ids, self.structures_stream.files):
            work_dir = os.path.join(self.output_folder, sid)
            cmd = f'CABSflex -i "{sfile}" -w "{work_dir}"'
            if flags_str:
                cmd += f" {flags_str}"

            if run_parallel:
                log_file = os.path.join(self.output_folder, f"{sid}.log")
                cabsflex_cmds += f"""
# Wait if we've reached the concurrency limit
while [ "${{#PIDS[@]}}" -ge {self.max_parallel} ]; do
    # Wait for any one job to finish, then remove completed PIDs
    NEW_PIDS=()
    for PID in "${{PIDS[@]}}"; do
        if kill -0 "$PID" 2>/dev/null; then
            NEW_PIDS+=("$PID")
        else
            wait "$PID" || FAILED=$((FAILED + 1))
        fi
    done
    PIDS=("${{NEW_PIDS[@]}}")
    if [ "${{#PIDS[@]}}" -ge {self.max_parallel} ]; then
        sleep 1
    fi
done
mkdir -p "{work_dir}"
(
    echo "=== Processing {sid} ==="
    {cmd}
) > "{log_file}" 2>&1 &
PIDS+=($!)
"""
            else:
                cabsflex_cmds += f"""
echo "=== Processing {sid} ==="
mkdir -p "{work_dir}"
{cmd}
if [ $? -ne 0 ]; then
    echo "Error: CABS-Flex failed for {sid}"
    exit 1
fi
"""

        if run_parallel:
            cabsflex_cmds += """
# Wait for remaining jobs
for PID in "${PIDS[@]}"; do
    wait "$PID" || FAILED=$((FAILED + 1))
done

if [ "$FAILED" -ne 0 ]; then
    echo "Error: $FAILED CABS-Flex job(s) failed"
    exit 1
fi
"""

        # Post-processing runs under biopipelines (Python 3)
        postproc_env = self.activate_environment(name="biopipelines")

        return f"""echo "Running CABS-Flex flexibility simulation"
echo "Input structures: {len(self.structures_stream)} files"
echo "Models per structure: {self.num_models}"

# --- Run CABSflex for each structure (Python 2.7 env) ---
{cabsflex_cmds}

echo "=== CABSflex runs complete, starting post-processing ==="

# --- Switch to Python 3 for post-processing ---
{postproc_env}

python "{self.helper_script}" \\
    --structures "{self.structures_ds_json}" \\
    --output_dir "{self.output_folder}" \\
    --rmsf_all_csv "{self.rmsf_all_csv}" \\
    --structures_map "{self.structures_map}" \\
    --num_models {self.num_models}

if [ $? -eq 0 ]; then
    echo "CABS-Flex completed successfully"
else
    echo "Error: CABS-Flex post-processing failed"
    exit 1
fi

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after CABS-Flex execution."""
        input_ids = self.structures_stream.ids

        # --- Output structures: num_models per input ---
        structure_ids = []
        structure_files = []
        provenance = {"structures": []}
        for input_id in input_ids:
            for model_idx in range(1, self.num_models + 1):
                output_id = f"{input_id}_{model_idx}"
                output_file = os.path.join(self.output_folder, f"{output_id}.pdb")
                structure_ids.append(output_id)
                structure_files.append(output_file)
                provenance["structures"].append(input_id)

        create_map_table(self.structures_map, structure_ids, files=structure_files, provenance=provenance)

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=structure_files,
            map_table=self.structures_map,
            format="pdb"
        )

        # --- Output images: SVG plots per input ---
        image_ids = []
        image_files = []
        svg_names = ["RMSF_seq.svg", "E_RMSD_A_total.svg", "RMSD_frame_A_replica_0.svg"]
        for input_id in input_ids:
            for svg_name in svg_names:
                img_id = f"{input_id}_{svg_name.replace('.svg', '')}"
                img_file = os.path.join(self.output_folder, f"{input_id}_{svg_name}")
                image_ids.append(img_id)
                image_files.append(img_file)

        images = DataStream(
            name="images",
            ids=image_ids,
            files=image_files,
            format="svg"
        )

        # --- Output tables ---
        rmsf_columns = ["id", "chain", "resi", "rmsf"]

        # Per-ID RMSF tables: res.tables.rmsf["<id>"]
        rmsf_per_id = IndexedTableContainer(
            name="rmsf",
            columns=rmsf_columns,
            description="Per-residue RMSF"
        )
        for input_id in input_ids:
            rmsf_path = os.path.join(self.output_folder, f"{input_id}_RMSF.csv")
            rmsf_per_id.add(input_id, path=rmsf_path)

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_map,
                columns=["id", "file", "structures.id"],
                description="CABS-Flex ensemble model structures",
                count=len(structure_ids)
            ),
            "rmsf_all": TableInfo(
                name="rmsf_all",
                path=self.rmsf_all_csv,
                columns=rmsf_columns,
                description="Per-residue RMSF from all input structures (merged)",
                count="variable"
            ),
            "rmsf": rmsf_per_id
        }

        return {
            "structures": structures,
            "images": images,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "cabsflex_params": {
                "num_models": self.num_models,
                "mc_cycles": self.mc_cycles,
                "mc_steps": self.mc_steps,
                "mc_annealing": self.mc_annealing,
                "temperature": self.temperature,
                "flexibility": self.flexibility,
                "filtering_count": self.filtering_count,
                "aa_rebuild": self.aa_rebuild,
                "restraints": self.restraints,
                "restraints_reduce": self.restraints_reduce,
                "weighted_fit": self.weighted_fit,
                "pdb_output": self.pdb_output,
                "max_parallel": self.max_parallel
            }
        })
        return base_dict
