# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Aggrescan3D (A3D): structure-based prediction of protein aggregation propensity.

A3D scores each residue of a folded protein for aggregation propensity using an
experimentally-derived aggregation scale combined with structural context
(solvent exposure). Negative scores indicate soluble residues, positive scores
indicate aggregation-prone residues; the per-structure average is a global
solubility/aggregation readout. A3D complements sequence-based solubility
predictors by acting on 3-D coordinates.

This wrapper runs A3D in static mode only: per-residue and global aggregation
scoring of an input structure. The FoldX-backed repair / mutation-enhancement
modes (-f, -m, -am) and the CABS-flex-backed dynamic mode (-d) are out of scope.

Reference:
    Paper: https://academic.oup.com/nar/article/47/W1/W300/5485072
    Standalone protocols: https://doi.org/10.1101/2020.09.09.276915
    Repository: https://bitbucket.org/lcbio/aggrescan3d
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import Resolve
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import Resolve


class Aggrescan3D(BaseConfig):
    """
    Aggrescan3D (A3D): structure-based aggregation propensity scoring.

    Runs A3D in static mode on each input structure, producing a per-residue
    aggregation profile and a per-structure global score.

    Inputs:
        structures: protein structures (StandardizedOutput or DataStream).

    Parameters:
        chains:   optional chain selection passed to A3D's -C (e.g. "A" or "AB").
                  None analyses all chains.
        distance: distance value (Angstrom) used in the A3D score calculation
                  (A3D -D, default 10).
        max_parallel: max structures to run concurrently (default 1 = sequential).

    Outputs:
        Streams:
            aggregation: per-residue A3D scores, one CSV per input structure
                         (resi-csv format; columns id, chain, resi, score).
                         Thresholdable by Selection to pick aggregation-prone
                         residues.
            structures:  A3D output PDB per input (B-factor field replaced with
                         the A3D score), format pdb.
            images:      per-chain A3D score plots (PNG).
        Tables:
            scores:          id | structures.id | avg_score | min_score |
                             max_score | n_residues | n_aggregation_prone
                             (per-structure global solubility/aggregation summary)
            aggregation_all: id | chain | resi | score (merged, all structures)
    """

    TOOL_NAME = "Aggrescan3D"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_manager} run -n Aggrescan3D python -c "import aggrescan" >/dev/null 2>&1; then
    echo "Aggrescan3D already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("Aggrescan3D", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("Aggrescan3D", env_manager, biopipelines)
        return f"""echo "=== Installing Aggrescan3D ==="
{skip}{remove_block}
{env_block}

# Verify installation
if {env_manager} run -n Aggrescan3D python -c "import aggrescan" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== Aggrescan3D installation complete ==="
else
    echo "ERROR: Aggrescan3D verification failed (cannot import aggrescan)"
    exit 1
fi
"""

    # Lazy path descriptors — canonical layout.
    #   _configuration/ — structures input DataStream JSON.
    #   _execution/     — per-structure A3D work dirs (scratch: A3D.csv, etc.).
    #   structures/     — A3D output PDBs + structures_map.
    #   images/         — per-chain score PNGs.
    #   aggregation/    — per-residue CSVs + aggregation_map.
    #   tables/         — scores (per-structure summary) + aggregation_all (merged).
    structures_ds_json = Path(lambda self: self.configuration_path("structures.json"))
    worklist_tsv = Path(lambda self: self.configuration_path("worklist.tsv"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    aggregation_map = Path(lambda self: self.stream_map_path("aggregation"))
    scores_csv = Path(lambda self: self.table_path("scores"))
    aggregation_all_csv = Path(lambda self: self.table_path("aggregation_all"))
    helper_script = Path(lambda self: self.pipe_script_path("pipe_aggrescan3d.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 chains: Optional[str] = None,
                 distance: float = 10.0,
                 max_parallel: int = 1,
                 **kwargs):
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(
                f"structures must be DataStream or StandardizedOutput, got {type(structures).__name__}"
            )

        self.chains = chains
        self.distance = distance
        self.max_parallel = max_parallel

        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if self.distance <= 0:
            raise ValueError(f"distance must be > 0, got: {self.distance}")
        if self.max_parallel < 1:
            raise ValueError(f"max_parallel must be >= 1, got: {self.max_parallel}")
        _validate_freeform_string("chains", self.chains)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.extend([
            f"INPUT STRUCTURES: {len(self.structures_stream)} files",
            f"CHAINS: {self.chains if self.chains else 'all'}",
            f"DISTANCE: {self.distance}",
            f"MAX PARALLEL: {self.max_parallel}",
        ])
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_ds_json)

        script_content = "#!/bin/bash\n"
        script_content += "# Aggrescan3D execution script\n"
        script_content += self.generate_completion_check_header()
        # Start in biopipelines (Python 3): the ID/worklist resolution and the
        # post-processing both need it. The A3D binary runs under its own
        # Python-2.7 env in the middle (see _generate_run_aggrescan); resolving
        # stream IDs there would crash, because resolve_stream_ids.py imports
        # the Python-3 biopipelines package.
        script_content += self.activate_environment(name="biopipelines")
        script_content += self._generate_run_aggrescan()
        script_content += self.generate_completion_check_footer()
        return script_content

    def _generate_run_aggrescan(self) -> str:
        """Run A3D per structure under the Aggrescan3D (Python 2.7) env, then
        post-process under biopipelines (Python 3)."""
        flags = ["-v 4"]
        if self.chains:
            flags.append(f"-C {self.chains}")
        if self.distance != 10.0:
            flags.append(f"-D {self.distance}")
        flags_str = " ".join(flags)

        n_structures = len(self.structures_stream)
        run_parallel = self.max_parallel > 1 and n_structures > 1

        container_prefix = self.container_prefix()
        execution_folder = self.execution_folder
        a3d_env = self.activate_environment()  # Aggrescan3D (Python 2.7) env
        postproc_env = self.activate_environment(name="biopipelines")

        if run_parallel:
            run_cmds = f"""
echo "Running structures in parallel (max {self.max_parallel} concurrent)"
PIDS=()
FAILED=0
while IFS=$'\\t' read -r struct_id STRUCT_FILE; do
    [ -z "$struct_id" ] && continue
    WORK_DIR="{execution_folder}/$struct_id"
    LOG_FILE="{execution_folder}/$struct_id.log"

    while [ "${{#PIDS[@]}}" -ge {self.max_parallel} ]; do
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
    (
        echo "=== Processing $struct_id ==="
        {container_prefix}aggrescan -i "$STRUCT_FILE" -w "$WORK_DIR" {flags_str}
    ) > "$LOG_FILE" 2>&1 &
    PIDS+=($!)
done < "{self.worklist_tsv}"

for PID in "${{PIDS[@]}}"; do
    wait "$PID" || FAILED=$((FAILED + 1))
done

if [ "$FAILED" -ne 0 ]; then
    echo "Error: $FAILED Aggrescan3D job(s) failed"
    exit 1
fi
"""
        else:
            run_cmds = f"""
while IFS=$'\\t' read -r struct_id STRUCT_FILE; do
    [ -z "$struct_id" ] && continue
    WORK_DIR="{execution_folder}/$struct_id"

    echo "=== Processing $struct_id ==="
    {container_prefix}aggrescan -i "$STRUCT_FILE" -w "$WORK_DIR" {flags_str}
    if [ $? -ne 0 ]; then
        echo "Error: Aggrescan3D failed for $struct_id"
        exit 1
    fi
done < "{self.worklist_tsv}"
"""

        return f"""echo "Running Aggrescan3D static-mode scoring"
echo "Input structures: {n_structures} files"

# Resolve the id->file worklist under biopipelines (Python 3) BEFORE switching
# to A3D's Python-2.7 env — resolve_stream_ids.py imports the Python-3 package
# and would crash under py2.7. The bash loop below then reads the TSV with no
# Python in the py2.7 section.
python "{self.helper_script}" \\
    --emit-worklist \\
    --structures "{self.structures_ds_json}" \\
    --worklist "{self.worklist_tsv}"
if [ $? -ne 0 ]; then
    echo "Error: failed to resolve Aggrescan3D worklist"
    exit 1
fi

# A3D's Python-2.7 matplotlib (per-chain plots) needs a headless backend. On
# Colab the host exports an inline backend the py2.7 env cannot import; Agg is
# correct on any headless runtime (Colab and HPC alike).
export MPLBACKEND=Agg

# --- Switch to the Aggrescan3D (Python 2.7) env and run A3D per structure ---
{a3d_env}
{run_cmds}

echo "=== Aggrescan3D runs complete, starting post-processing ==="

# --- Switch back to Python 3 for post-processing ---
{postproc_env}

python "{self.helper_script}" \\
    --structures "{self.structures_ds_json}" \\
    --work_root "{execution_folder}" \\
    --structures_dir "{self.stream_folder('structures')}" \\
    --images_dir "{self.stream_folder('images')}" \\
    --aggregation_dir "{self.stream_folder('aggregation')}" \\
    --structures_map "{self.structures_map}" \\
    --aggregation_map "{self.aggregation_map}" \\
    --scores_csv "{self.scores_csv}" \\
    --aggregation_all_csv "{self.aggregation_all_csv}"

if [ $? -eq 0 ]; then
    echo "Aggrescan3D completed successfully"
else
    echo "Error: Aggrescan3D post-processing failed"
    exit 1
fi

"""

    def get_output_files(self) -> Dict[str, Any]:
        input_ids = list(self.structures_stream.ids)

        # Output structures: A3D output.pdb, one per input (ids preserved).
        structures = DataStream(
            name="structures",
            ids=input_ids,
            files=[self.stream_path("structures", "<id>.pdb")],
            map_table=self.structures_map,
            format="pdb",
        )

        # Per-residue aggregation profile: one resi-csv per input structure.
        aggregation_files = [self.stream_path("aggregation", "<id>_A3D.csv")]
        aggregation_stream = DataStream(
            name="aggregation",
            ids=input_ids,
            files=aggregation_files,
            map_table=self.aggregation_map,
            format="resi-csv",
        )

        # Per-chain score plots (PNG). A3D names them <chainID>.png; the pipe
        # script prefixes the structure id so ids stay unique across inputs.
        images = DataStream(
            name="images",
            ids=input_ids,
            files=[self.stream_path("images", "<id>.png")],
            format="png",
        )

        tables = {
            "scores": TableInfo(
                name="scores",
                path=self.scores_csv,
                columns=["id", "structures.id", "avg_score", "min_score",
                         "max_score", "n_residues", "n_aggregation_prone"],
                description="Per-structure A3D global aggregation/solubility summary "
                            "(more negative avg_score = more soluble)",
            ),
            "aggregation_all": TableInfo(
                name="aggregation_all",
                path=self.aggregation_all_csv,
                columns=["id", "chain", "resi", "score"],
                description="Per-residue A3D aggregation scores from all input structures (merged)",
            ),
        }

        return {
            "structures": structures,
            "images": images,
            "aggregation": aggregation_stream,
            "tables": tables,
            "output_folder": self.output_folder,
        }
