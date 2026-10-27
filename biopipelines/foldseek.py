# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Foldseek tool: structural clustering at TM-score thresholds.

Runs `foldseek easy-cluster` (one protein chain per structure) or
`foldseek easy-multimercluster` (whole complexes) once per TM-score threshold,
and reports for every structure the cluster it falls in and the number of
clusters the whole set splits into. The block mirrors the Foldseek block of the
RosettaCommons Chaperone framework.

Reference: van Kempen et al., Nature Biotechnology, doi:10.1038/s41587-023-01773-0. https://github.com/steineggerlab/foldseek
"""

import os
from typing import Dict, List, Any, Optional, Union

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


DEFAULT_TM_THRESHOLDS = (0.5, 0.6, 0.7, 0.8)
# Foldseek --alignment-type: 1 is TM-align, 2 is the 3Di+AA local alignment.
ALIGNMENT_TYPES = (1, 2)
DEFAULT_MULTIMER_TM_THRESHOLD = 0.65
DEFAULT_INTERFACE_LDDT_THRESHOLD = 0.65


def threshold_tag(tm: float) -> str:
    """Column suffix of a TM-score threshold: 0.5 -> 'tm50', 0.65 -> 'tm65'."""
    return f"tm{int(round(tm * 100)):02d}"


class Foldseek(BaseConfig):
    """
    Foldseek: cluster structures at TM-score thresholds.

    Answers "how many distinct folds are in this set, and which designs share
    one" by Foldseek's structural alignment, once per threshold in
    `tm_thresholds`. Cluster counts describe diversity under this protocol; they
    are not a per-design quality score.

    Not StructureCluster. StructureCluster is a base-env triage tool: a
    superposition-free distance-matrix proxy for TM-score, ranked on metrics you
    hand it, built for thousands of models. Foldseek runs the real aligner
    (TM-align by default), at several fixed TM-score thresholds, optionally on
    whole complexes, and ranks nothing. Reach for StructureCluster to choose
    which families to look at; reach for Foldseek to report diversity at
    published TM-score cutoffs.

    TM-SCORE THRESHOLD (read this before comparing counts across studies)
    Foldseek keeps an alignment when its TM-score exceeds the threshold
    (`--tmscore-threshold`). Foldseek's default normalizes that TM-score by the
    aligned region, and no coverage requirement is set, so a short well-aligned
    segment can join two structures. Counts therefore tend to understate
    diversity.

    MULTIMER MODE
    With `multimer=True` the structures are clustered as complexes with
    `foldseek easy-multimercluster`. `tm_thresholds` then sets
    `--chain-tm-threshold` (every chain must align above it), while the whole
    complex must also exceed `multimer_tm_threshold` and the interface must
    exceed `interface_lddt_threshold`. Chain names must not contain underscores.
    Without it every structure is expected to hold ONE protein chain; a file
    with several chains is recorded as a failure instead of being assigned to
    one of its chains by guess.

    Inputs:
        structures: PDB or mmCIF structures.
        tm_thresholds: TM-score thresholds, each in (0, 1) with at most two
            decimals (default [0.5, 0.6, 0.7, 0.8]). One Foldseek run each.
        multimer: cluster complexes with easy-multimercluster instead of single
            chains with easy-cluster (default False).
        alignment_type: Foldseek `--alignment-type`: 1 TM-align, 2 3Di+AA local
            alignment followed by the TM-score filter (default 1).
        min_seq_id: Foldseek `--min-seq-id`, the minimum sequence identity to
            cluster (default 0.0, no requirement). Single-chain mode only.
        multimer_tm_threshold: Foldseek `--multimer-tm-threshold` (default 0.65
            when `multimer=True`). Multimer mode only.
        interface_lddt_threshold: Foldseek `--interface-lddt-threshold` (default
            0.65 when `multimer=True`). Multimer mode only.
        threads: CPU threads Foldseek uses (`--threads`, default 1). Match the
            CPUs requested with Resources.

    Outputs:
        Tables:
            tm_clusters: id | cluster_<t> | representative_<t> | n_clusters_<t>
                for every threshold, where <t> is the percent form of the
                threshold (tm50, tm60, tm70, tm80 by default). `cluster_<t>` is
                the cluster number (1 is the largest), `representative_<t>` the
                id of the cluster's representative structure, `n_clusters_<t>`
                the number of clusters the whole set splits into (the same value
                on every row). One row per structure that entered; a structure
                Foldseek did not report carries NaN.
            missing: id | removed_by | kind | cause

    Usage::

        folds = ESMFold2(proteins=seqs)
        diversity = Foldseek(structures=folds, threads=8)

        # one structure per cluster at TM 0.5
        reps = Panda(tables=diversity.tables.tm_clusters,
                     operations=[Panda.filter("id == representative_tm50")],
                     pool=folds)

    See Also:
        StructureCluster: distance-matrix triage clustering with per-family
            ranking, no external binary.
    """

    TOOL_NAME = "Foldseek"
    TOOL_VERSION = "1.0"
    ENV_NAME = "foldseek"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check(env, env_manager)
        pull = cls._container_pull_block(folders, force_reinstall)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "Foldseek already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        return f"""echo "=== Installing Foldseek ==="
{pull}{skip}{remove_block}
{env_block}

if {cls._env_run(env, env_manager)}foldseek version >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== Foldseek installation complete ==="
else
    echo "ERROR: Foldseek verification failed (foldseek binary not on PATH)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    tm_clusters_csv = Path(lambda self: self.table_path("tm_clusters"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_foldseek.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 tm_thresholds: Optional[List[float]] = None,
                 multimer: bool = False,
                 alignment_type: int = 1,
                 min_seq_id: float = 0.0,
                 multimer_tm_threshold: Optional[float] = None,
                 interface_lddt_threshold: Optional[float] = None,
                 threads: int = 1,
                 **kwargs):
        self.structures = structures  # raw handle kept for missing-propagation
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        self.tm_thresholds = sorted(float(t) for t in (
            DEFAULT_TM_THRESHOLDS if tm_thresholds is None else tm_thresholds))
        self.multimer = bool(multimer)
        self.alignment_type = int(alignment_type)
        self.min_seq_id = float(min_seq_id)
        self.multimer_tm_threshold = multimer_tm_threshold
        self.interface_lddt_threshold = interface_lddt_threshold
        if self.multimer:
            if multimer_tm_threshold is None:
                self.multimer_tm_threshold = DEFAULT_MULTIMER_TM_THRESHOLD
            if interface_lddt_threshold is None:
                self.interface_lddt_threshold = DEFAULT_INTERFACE_LDDT_THRESHOLD
        self.threads = int(threads)
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.tm_thresholds:
            raise ValueError("tm_thresholds must hold at least one TM-score threshold")
        for t in self.tm_thresholds:
            if not 0.0 < t < 1.0 or abs(t * 100 - round(t * 100)) > 1e-6:
                raise ValueError(
                    f"tm_thresholds entries must lie in (0, 1) with at most two decimals, got {t}")
        if len({threshold_tag(t) for t in self.tm_thresholds}) != len(self.tm_thresholds):
            raise ValueError(f"tm_thresholds must be distinct, got {self.tm_thresholds}")
        if self.alignment_type not in ALIGNMENT_TYPES:
            raise ValueError(
                f"alignment_type must be one of {ALIGNMENT_TYPES} (1 TM-align, 2 3Di+AA), "
                f"got {self.alignment_type}")
        if not 0.0 <= self.min_seq_id <= 1.0:
            raise ValueError(f"min_seq_id must lie in [0, 1], got {self.min_seq_id}")
        if self.threads < 1:
            raise ValueError("threads must be at least 1")
        if self.multimer:
            if self.min_seq_id > 0.0:
                raise ValueError("min_seq_id applies to single-chain clustering only (multimer=False)")
            for name in ("multimer_tm_threshold", "interface_lddt_threshold"):
                if not 0.0 <= getattr(self, name) <= 1.0:
                    raise ValueError(f"{name} must lie in [0, 1], got {getattr(self, name)}")
        elif self.multimer_tm_threshold is not None or self.interface_lddt_threshold is not None:
            raise ValueError(
                "multimer_tm_threshold and interface_lddt_threshold apply to multimer=True only")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"TM THRESHOLDS: {', '.join(str(t) for t in self.tm_thresholds)}")
        lines.append(f"CLUSTERING: {'complexes (easy-multimercluster)' if self.multimer else 'single chains (easy-cluster)'}")
        lines.append(f"ALIGNMENT TYPE: {self.alignment_type}")
        lines.append(f"THREADS: {self.threads}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self._write_config_yaml()
        script = "#!/bin/bash\n"
        script += "# Foldseek structural clustering script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Running Foldseek on {len(self.structures_stream)} structure(s)"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-json "{self.structures_json}" \\
    --scratch-dir "{self.execution_path()}" \\
    --tm-clusters-csv "{self.tm_clusters_csv}" \\
    --local-missing-csv "{self.local_missing_csv}" \\
    --container-prefix "{self.container_prefix()}"
"""
        script += self.generate_missing_propagation(
            self.structures, local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        """Dump the parameters to configuration/config.yaml for the pipe to read."""
        import yaml
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump({
                "tm_thresholds": list(self.tm_thresholds),
                "multimer": self.multimer,
                "alignment_type": self.alignment_type,
                "min_seq_id": self.min_seq_id,
                "multimer_tm_threshold": self.multimer_tm_threshold,
                "interface_lddt_threshold": self.interface_lddt_threshold,
                "threads": self.threads,
            }, f, sort_keys=False)

    def get_output_files(self) -> Dict[str, Any]:
        columns = ["id"]
        for t in self.tm_thresholds:
            tag = threshold_tag(t)
            columns += [f"cluster_{tag}", f"representative_{tag}", f"n_clusters_{tag}"]
        tables = {
            "tm_clusters": TableInfo(
                name="tm_clusters",
                path=self.tm_clusters_csv,
                columns=columns,
                description="Foldseek cluster number, representative and cluster count per TM-score threshold",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
