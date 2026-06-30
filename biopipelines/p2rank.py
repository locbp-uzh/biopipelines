# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""P2Rank tool: template-free ligand-binding-site prediction.

For each input structure, runs the `prank predict` binary (bioconda) and
reports the predicted pockets (ranked by ligandability) together with their
per-residue scores. P2Rank scores and clusters points on the protein's
solvent-accessible surface; no ligand or template is required.

Reference: Krivak & Hoksza (2018) P2Rank. https://github.com/rdk/p2rank
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


class P2Rank(BaseConfig):
    """
    P2Rank: predict ligand-binding pockets on protein structures.

    Inputs:
        structures: input PDB/CIF structures.
        config: P2Rank model preset (passed as `prank predict -c <config>`).
                One of "default", "alphafold" (tuned for AlphaFold models,
                uses pLDDT in the B-factor column), or "conservation"
                (requires precomputed conservation scores). Default "default".
        threads: number of CPU threads P2Rank uses (`prank predict -threads`).
                 Default 1.
        visualizations: emit P2Rank's PyMOL/visualization artefacts
                 (`prank predict -visualizations`). Default False (faster).

    Outputs:
        Streams:
            residues:  per-residue resi-csv (one <id>.csv per input) with columns
                       id | chain | resi | resn | pocket_idx | score | probability.
                       Consumable by the Selection tool to retrieve pocket residues,
                       e.g. ``Selection.add(p2rank.streams.residues, include="probability>0.5")``.
        Tables:
            pockets:   id | pocket_idx | rank | score | probability |
                       n_residues | residues | center_x | center_y | center_z
                       (one row per predicted pocket; `residues` is a
                       chain-aware selection string e.g. "A12+A45-47")
            residues:  id | chain | resi | resn | pocket_idx | score | probability
                       (one row per protein residue; pocket_idx is the pocket
                       a residue belongs to, or 0 if unassigned)
            summary:   id | n_pockets | top_score | top_probability | top_residues
                       (top_residues = residues of the rank-1 pocket)
            missing:   id | cause
    """

    TOOL_NAME = "P2Rank"
    TOOL_VERSION = "1.1"

    _CONFIGS = ("default", "alphafold", "conservation")

    # P2Rank is distributed as a JVM release tarball (the bioconda package was
    # removed). We download it into the P2Rank repository folder and symlink its
    # `prank` launcher onto the env's PATH so the pipe script can call it bare.
    P2RANK_VERSION = "2.5.1"
    P2RANK_URL = (
        "https://github.com/rdk/p2rank/releases/download/"
        f"{P2RANK_VERSION}/p2rank_{P2RANK_VERSION}.tar.gz"
    )

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("P2Rank", "")
        prank_dir = f"{repo_dir}/p2rank_{cls.P2RANK_VERSION}"
        env_check = cls._env_exists_check("p2rank", env_manager)
        # Skip-check is file-based (no `prank --help`): booting the JVM just to
        # probe install state is slow and, on memory-constrained runtimes like
        # Colab, can crash the kernel. The jar + launcher existing is sufficient.
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check} && [ -f "{prank_dir}/bin/p2rank.jar" ] && [ -x "$(dirname "$({env_manager} run -n p2rank which java)")/prank" ]; then
    echo "P2Rank already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("p2rank", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("p2rank", env_manager, biopipelines)
        return f"""echo "=== Installing P2Rank ==="
{skip}{remove_block}
{env_block}

# Download + unpack the P2Rank release tarball (bioconda no longer ships it).
mkdir -p "{repo_dir}"
if [ ! -x "{prank_dir}/prank" ]; then
    echo "Downloading P2Rank {cls.P2RANK_VERSION}"
    curl -fsSL -o "{repo_dir}/p2rank.tar.gz" "{cls.P2RANK_URL}"
    tar xzf "{repo_dir}/p2rank.tar.gz" -C "{repo_dir}"
    rm -f "{repo_dir}/p2rank.tar.gz"
    chmod +x "{prank_dir}/prank"
fi

# Put `prank` on the env's PATH (env bin = parent of its `java`) so the pipe
# script's bare `prank` resolves after activation. We CANNOT symlink the
# launcher: it derives its install dir from `dirname ${{BASH_SOURCE[0]}}`, which
# through a symlink resolves to the env bin (no p2rank.jar there) and crashes
# with ClassNotFoundException. Write a thin wrapper that execs the real
# launcher by absolute path so its BASH_SOURCE points at the real install.
ENV_BIN="$(dirname "$({env_manager} run -n p2rank which java)")"
# Constrain the JVM via JAVA_OPTS (the prank launcher appends its own -Xmx2048m
# AFTER $JAVA_OPTS, and the last -Xmx wins, so we keep 2g and add the rest).
# On Colab the JVM logs "Cgroup memory controller path seems to have moved...
# detected limits won't be accurate" and then sizes heap/GC threads against the
# WRONG (host-sized) limits, ballooning memory until the kernel watchdog kills
# the session. -XX:-UseContainerSupport stops it trusting the broken cgroup,
# -XX:MaxRAM caps the ceiling it reasons from, SerialGC + ActiveProcessorCount=2
# stop the parallel GC spawning a host-sized thread pool.
cat > "$ENV_BIN/prank" <<WRAP
#!/usr/bin/env bash
export JAVA_OPTS="-Xmx2048m -XX:-UseContainerSupport -XX:MaxRAM=3g -XX:+UseSerialGC -XX:ActiveProcessorCount=2 \\$JAVA_OPTS"
exec "{prank_dir}/prank" "\\$@"
WRAP
chmod +x "$ENV_BIN/prank"

# Lightweight verification: the jar is unpacked and the launcher is on the env
# PATH. We deliberately do NOT boot the JVM (`prank --help`) here — that adds a
# slow cold start to every install and some P2Rank builds exit non-zero on
# --help, which would fail the check spuriously. The first real run exercises
# the JVM.
if [ -f "{prank_dir}/bin/p2rank.jar" ] && [ -x "$ENV_BIN/prank" ]; then
    touch "$INSTALL_SUCCESS"
    echo "=== P2Rank installation complete ==="
else
    echo "ERROR: P2Rank verification failed (p2rank.jar or prank launcher missing)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    pockets_csv = Path(lambda self: self.table_path("pockets"))
    residues_csv = Path(lambda self: self.table_path("residues"))
    residues_map = Path(lambda self: self.stream_map_path("residues"))
    summary_csv = Path(lambda self: self.table_path("summary"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_p2rank.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 config: str = "default",
                 threads: int = 1,
                 visualizations: bool = False,
                 **kwargs):
        self.structures = structures
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        self.config = config
        self.threads = int(threads)
        self.visualizations = bool(visualizations)
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if self.config not in self._CONFIGS:
            raise ValueError(f"config must be one of {self._CONFIGS}, got '{self.config}'")
        if self.threads < 1:
            raise ValueError("threads must be at least 1")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"STRUCTURES: {len(self.structures_stream)} files")
        lines.append(f"CONFIG: {self.config}")
        lines.append(f"THREADS: {self.threads}")
        lines.append(f"VISUALIZATIONS: {self.visualizations}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        upstream_missing_path = self._get_upstream_missing_table_path(self.structures)
        upstream_missing_flag = f' \\\n    --upstream-missing "{upstream_missing_path}"' if upstream_missing_path else ""
        script = "#!/bin/bash\n"
        script += "# P2Rank binding-site prediction script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Running P2Rank on {len(self.structures_stream)} structure(s) (config={self.config})"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --config "{self.config}" \\
    --threads {self.threads} \\
    --visualizations {1 if self.visualizations else 0} \\
    --scratch-dir "{self.extras_path()}" \\
    --pockets-csv "{self.pockets_csv}" \\
    --residues-csv "{self.residues_csv}" \\
    --residues-dir "{self.stream_folder('residues')}" \\
    --residues-map-csv "{self.residues_map}" \\
    --summary-csv "{self.summary_csv}" \\
    --container-prefix "{self.container_prefix()}" \\
    --missing-csv "{self.missing_csv}"{upstream_missing_flag}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        residues_stream = DataStream(
            name="residues",
            ids=self.structures_stream.ids,
            files=[self.stream_path("residues", "<id>.csv")],
            map_table=self.residues_map,
            format="resi-csv",
        )
        tables = {
            "pockets": TableInfo(
                name="pockets",
                path=self.pockets_csv,
                columns=["id", "pocket_idx", "rank", "score", "probability",
                         "n_residues", "residues", "center_x", "center_y", "center_z"],
                description="P2Rank predicted pockets (one row per pocket)",
            ),
            "residues": TableInfo(
                name="residues",
                path=self.residues_csv,
                columns=["id", "chain", "resi", "resn", "pocket_idx", "score", "probability"],
                description="P2Rank per-residue ligandability scores",
            ),
            "summary": TableInfo(
                name="summary",
                path=self.summary_csv,
                columns=["id", "n_pockets", "top_score", "top_probability", "top_residues"],
                description="P2Rank per-structure pocket summary",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="IDs removed (upstream or local failure) with removal reason",
            ),
        }
        return {
            "residues": residues_stream,
            "tables": tables,
            "output_folder": self.output_folder,
        }
