# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
StructureCluster: group predicted structures by fold and rank the groups.

Answers the triage question at the top of a large design campaign -- "what
topologies did I actually get, and which of them fold well?" -- without asking
the user to look at ten thousand models. It clusters a structures stream by
backbone shape, then ranks the clusters on whatever per-id metrics are handed
to it (mean pLDDT, mean ipTM, mean buried ligand surface, ...).

Superposition-free by construction. The similarity is a distance-matrix RMSD
(dRMSD) between length-normalised C-alpha traces, so no alignment is ever
computed. That is what makes 10,000 structures tractable in the base env with
numpy alone; it is also the tool's main approximation, and the docstring says so
where a user will see it.
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


# Cluster labels are assigned AFTER ranking, so cluster_001 is the best cluster
# under `rank_by`. Zero-padded to keep lexical and numeric order identical --
# a plain str(i) would sort cluster_10 before cluster_2 in every spreadsheet.
_LABEL_FMT = "cluster_{:03d}"

_VALID_SS = ("helix_frac", "strand_frac", "coil_frac")


class StructureCluster(BaseConfig):
    """
    StructureCluster: fold-similarity clustering of a structures stream, with
    per-cluster ranking on user-supplied metrics.

    Intended for the first look at a large de novo campaign: how many distinct
    topologies came out, how big each family is, and which families fold with
    high confidence. It is a triage instrument, not a structural-biology
    measurement -- see SIMILARITY below.

    SIMILARITY (read this before setting `threshold`)
    -------------------------------------------------
    Each structure is reduced to a length-normalised C-alpha trace: the chain's
    C-alpha coordinates are resampled to `n_points` evenly spaced positions
    along the residue index, and the full pairwise distance matrix of those
    points is flattened. That vector is translation- and rotation-invariant, so
    two structures are compared by

        dRMSD = sqrt(mean((D_a - D_b)**2))            [Angstrom]

    with no superposition. dRMSD is converted to a bounded score with the
    TM-score normalisation,

        similarity = 1 / (1 + (dRMSD / d0)**2),
        d0 = max(0.5, 1.24 * (L - 15)**(1/3) - 1.8),  L = mean residue count

    so `similarity` lands in (0, 1] and the familiar 0.5 fold-identity cutoff is
    a sensible starting `threshold`.

    This is NOT TM-score and this tool is NOT TM-align. dRMSD on a resampled
    trace is a cheaper, stricter proxy: it is sensitive to overall size and to
    internal register in a way an optimal superposition is not, and it cannot
    recognise two folds related by a large rigid-body domain motion. Treat 0.5
    as a starting point to be tuned against the run's own distribution (the
    `assignments` table carries `similarity_to_representative` for exactly
    that), not as a value transferred from the TM-score literature.

    CLUSTERING
    ----------
    Sphere exclusion (leader / Taylor-Butina). Structures are visited in
    ranking order -- best first, by `rank_by` when metrics are supplied, by id
    otherwise -- and each unassigned structure becomes a representative that
    absorbs every remaining unassigned structure within `threshold` of it. So a
    cluster's representative is its best-scoring member, which is the model you
    want to open. Cost is O(N * n_clusters), not O(N^2), which is why 10,000
    inputs are routine.

    REPRESENTATIVE vs MEDOID -- they are different members and both are
    reported. The representative is the cluster's best-scoring model under
    `rank_by` (sphere exclusion visits in ranking order), so it is the one to
    take forward; it can sit at the edge of its own cluster. The medoid
    maximises mean similarity to the other members, so it is the most typical
    fold of the family and the one to look at when asking what this topology
    actually is. For a singleton they coincide.

    Note the consequence of the greedy rule: a structure joins the FIRST
    representative that captures it, not necessarily its nearest one. Sphere
    exclusion trades that for the linear cost and for representatives that are
    the best rather than the most central members.

    TOPOLOGY DESCRIPTORS
    --------------------
    Secondary-structure fractions come from C-alpha geometry alone (P-SEA-style
    i->i+2 / i->i+3 / i->i+4 distance criteria), because the whole point is to
    stay in the base env and run on ten thousand models. They are an estimate,
    good enough to separate an all-alpha bundle from an alpha/beta sandwich, and
    they are NOT DSSP. Run the `DSSP` tool on the cluster representatives when a
    real hydrogen-bond assignment matters.

    The distance bands are narrow, so the fractions are sensitive to backbone
    noise: on crystallographic input they are close to DSSP (86% helix for
    bacteriorhodopsin, 68% strand for a beta-sandwich design), but ~1 A of
    per-atom coordinate noise pulls them toward coil by 20-30 percentage points.
    Read them as a comparison BETWEEN models from the same predictor, not as
    absolute secondary-structure content.

    `relative_contact_order` is the classic fold-locality descriptor: the mean
    sequence separation of C-alpha contacts, divided by chain length. Low values
    mean local, helical topologies; high values mean the fold is stitched
    together by long-range (typically beta) contacts. It separates topologies
    that happen to share a secondary-structure composition, and it is far less
    noise-sensitive than the secondary-structure fractions.

    Inputs:
        structures: Predicted structures (PDB or mmCIF; ESMFold2's mmCIF is
            read directly). The protein chain analysed is `chain`, or the
            longest protein chain when `chain` is unset.
        metrics: Per-id table(s) carrying the columns named in `rank_by`. A
            TableInfo (``fold.tables.confidence``), a StandardizedOutput, a
            path, or a list of any of those; multiple tables are joined on
            ``id``. Optional -- without it, clusters are ranked by size.
        rank_by: Metric columns to average per cluster and rank the clusters on,
            in priority order (default ``["plddt", "iptm"]``). A column absent
            from every metrics table raises at config time only if it can be
            checked there; otherwise the run records NaN for it rather than
            failing, and the cluster ranking falls back to the next column.
        ascending: Sort direction for `rank_by`, a bool or a per-column list
            (default False -- higher is better, which is right for pLDDT, ipTM
            and buried surface alike).
        threshold: Similarity cutoff in (0, 1] for joining a cluster
            (default 0.5). Higher means tighter, more numerous clusters.
        n_points: C-alpha trace resampling length (default 64). Raising it
            sharpens discrimination between similar folds and costs
            O(n_points^2) memory per structure.
        chain: Chain id to analyse. Unset (default) takes the longest protein
            chain, which is the designed chain in every de novo pipeline.
        min_residues: Structures shorter than this are dropped from the
            clustering with ``kind="filter"`` (default 20). They still get a
            row in `assignments`, with a null cluster.
        max_structures: Safety ceiling on inputs (default 0, meaning no limit).
            Set it to fail fast rather than discover an OOM after a long run.

    Outputs:
        Tables:
            assignments: id | cluster | cluster_rank | is_representative |
                is_medoid | similarity_to_representative | n_residues | chain |
                radius_of_gyration | helix_frac | strand_frac | coil_frac |
                relative_contact_order | <each rank_by column>
            clusters: cluster | cluster_rank | size | fraction |
                representative | medoid | mean_<m> | median_<m> | min_<m> | max_<m>
                (for each rank_by column m) | mean_n_residues |
                mean_helix_frac | mean_strand_frac | mean_coil_frac |
                mean_relative_contact_order | mean_radius_of_gyration
            missing: id | removed_by | kind | cause

    Both tables carry one row per entity that entered -- a structure that could
    not be parsed keeps its `assignments` row with NaN values, so a downstream
    merge never silently shortens.

    Usage::

        folds = ESMFold2(proteins=seqs, ligands=dye)
        burial = SASA(structures=folds, ligand=dye, mode="ligand")

        families = StructureCluster(
            structures=folds,
            metrics=[folds.tables.confidence, burial.tables.sasa],
            rank_by=["plddt", "iptm", "delta_sasa"],
            threshold=0.5,
        )

        # the representative of each cluster, as an actual structures stream
        reps = Panda(tables=families.tables.assignments,
                     operations=[Panda.filter("is_representative == True"),
                                 Panda.sort("cluster_rank")],
                     pool=folds)

    See Also:
        EnsembleAnalysis: per-residue RMSF within one ensemble (a different
            question -- flexibility, not fold identity).
        DSSP: real hydrogen-bond secondary structure for the representatives.
        Panda: merging `assignments` back into a scoring table, or selecting
            representatives as a stream via ``pool=``.
    """

    TOOL_NAME = "StructureCluster"
    TOOL_VERSION = "1.1"

    # No _install_script: runs in the biopipelines base env (numpy + pandas +
    # the shared pdb_parser, which already reads mmCIF). Nothing to install.

    # ---------------------------------------------------------------- paths
    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    metrics_json = Path(lambda self: self.configuration_path("metrics.json"))
    assignments_csv = Path(lambda self: self.table_path("assignments"))
    clusters_csv = Path(lambda self: self.table_path("clusters"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_structure_cluster.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 metrics: Optional[Union[TableInfo, StandardizedOutput, str, List[Any]]] = None,
                 rank_by: Optional[List[str]] = None,
                 ascending: Union[bool, List[bool]] = False,
                 threshold: float = 0.5,
                 n_points: int = 64,
                 chain: str = "",
                 min_residues: int = 20,
                 max_structures: int = 0,
                 **kwargs):
        self.structures = structures  # raw handle kept for missing-propagation
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(
                f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        if metrics is None:
            self.metrics_input: List[Any] = []
        elif isinstance(metrics, list):
            self.metrics_input = list(metrics)
        else:
            self.metrics_input = [metrics]

        self.rank_by = list(rank_by) if rank_by else ["plddt", "iptm"]
        self.ascending = ascending
        self.threshold = threshold
        self.n_points = n_points
        self.chain = chain
        self.min_residues = min_residues
        self.max_structures = max_structures
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not 0.0 < self.threshold <= 1.0:
            raise ValueError(
                f"threshold must be a similarity in (0, 1], got {self.threshold}. "
                "It is compared against 1/(1+(dRMSD/d0)^2), not against an Angstrom "
                "distance -- 0.5 mirrors the TM-score fold-identity convention.")
        if self.n_points < 8:
            raise ValueError(
                f"n_points must be >= 8 to describe a fold, got {self.n_points}")
        if self.min_residues < 4:
            raise ValueError(f"min_residues must be >= 4, got {self.min_residues}")
        if self.max_structures < 0:
            raise ValueError(
                f"max_structures must be >= 0 (0 disables the ceiling), got {self.max_structures}")
        if not self.rank_by:
            raise ValueError(
                "rank_by must name at least one metric column, or be left at its default. "
                "To rank purely by cluster size, omit `metrics` instead.")
        if isinstance(self.ascending, list) and len(self.ascending) != len(self.rank_by):
            raise ValueError(
                f"ascending has {len(self.ascending)} entries but rank_by has "
                f"{len(self.rank_by)}; pass a single bool or one per column.")
        if self.metrics_input and not self.rank_by:
            raise ValueError("metrics were supplied but rank_by is empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders
        # Resolved here rather than in __init__: a TableInfo's path is only
        # final once the producing tool has its execution order.
        self.metrics_paths = [self._resolve_table_path(t) for t in self.metrics_input]

    def _resolve_table_path(self, table_input: Any) -> str:
        """TableInfo / TableReference / StandardizedOutput / str -> csv path.

        Mirrors Panda._resolve_table_path; the TableInfo branch must come first
        because TableInfo.__getattr__ answers ANY attribute with a
        TableReference, so the duck-typed check below would swallow it.
        """
        if isinstance(table_input, TableInfo):
            return table_input.info.path
        if hasattr(table_input, "path") and hasattr(table_input, "column"):
            return table_input.path
        if isinstance(table_input, str):
            return table_input
        tables_attr = getattr(table_input, "tables", None)
        if tables_attr is not None and hasattr(tables_attr, "_tables"):
            for name, info in tables_attr._tables.items():
                if name not in ("missing", "result"):
                    return info.info.path
            raise ValueError(
                f"metrics input {table_input!r} exposes no usable table "
                "(only 'missing'/'result'). Pass the table explicitly, e.g. "
                "fold.tables.confidence.")
        raise ValueError(
            f"metrics entries must be a TableInfo, StandardizedOutput or path string, "
            f"got {type(table_input).__name__}")

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"SIMILARITY THRESHOLD: {self.threshold} (dRMSD-based proxy, not TM-align)")
        lines.append(f"TRACE POINTS: {self.n_points}")
        lines.append(f"RANK BY: {', '.join(self.rank_by)}")
        lines.append(f"CHAIN: {self.chain or 'longest protein chain'}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self._write_config_yaml()

        script = "#!/bin/bash\n"
        script += "# StructureCluster: fold-similarity clustering + per-cluster ranking\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Clustering structures by fold similarity"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-json "{self.structures_json}" \\
    --metrics-json "{self.metrics_json}" \\
    --assignments-csv "{self.assignments_csv}" \\
    --clusters-csv "{self.clusters_csv}" \\
    --local-missing-csv "{self.local_missing_csv}"
"""
        # Structures dropped upstream must stay excused rather than counting as
        # this tool's failures; `metrics` is a table, not a filterable stream,
        # so only `structures` is propagated.
        script += self.generate_missing_propagation(
            self.structures, local_missing=self.local_missing_csv,
            missing_csv=self.missing_csv,
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        """Params to configuration/config.yaml; only paths travel as flags."""
        import json
        import yaml

        ascending = self.ascending
        if isinstance(ascending, bool):
            ascending = [ascending] * len(self.rank_by)

        with open(self.config_yaml, "w") as f:
            yaml.safe_dump({
                "threshold": float(self.threshold),
                "n_points": int(self.n_points),
                "chain": self.chain,
                "min_residues": int(self.min_residues),
                "max_structures": int(self.max_structures),
                "rank_by": list(self.rank_by),
                "ascending": [bool(a) for a in ascending],
            }, f, sort_keys=False)

        # Written even when empty, so the pipe's contract is one shape.
        with open(self.metrics_json, "w") as f:
            json.dump({"tables": list(getattr(self, "metrics_paths", []))}, f, indent=2)

    def get_output_files(self) -> Dict[str, Any]:
        # Metric columns ride along in `assignments` so a single table answers
        # "which cluster is this design in, and how good was it".
        assignment_cols = [
            "id", "cluster", "cluster_rank", "is_representative", "is_medoid",
            "similarity_to_representative", "n_residues", "chain",
            "radius_of_gyration", "helix_frac", "strand_frac", "coil_frac",
            "relative_contact_order",
        ] + list(self.rank_by)

        cluster_cols = ["cluster", "cluster_rank", "size", "fraction",
                        "representative", "medoid"]
        for m in self.rank_by:
            cluster_cols += [f"mean_{m}", f"median_{m}", f"min_{m}", f"max_{m}"]
        cluster_cols += [
            "mean_n_residues", "mean_helix_frac", "mean_strand_frac",
            "mean_coil_frac", "mean_relative_contact_order",
            "mean_radius_of_gyration",
        ]

        tables = {
            "assignments": TableInfo(
                name="assignments",
                path=self.assignments_csv,
                columns=assignment_cols,
                description="Per-structure cluster membership and fold descriptors",
            ),
            "clusters": TableInfo(
                name="clusters",
                path=self.clusters_csv,
                columns=cluster_cols,
                description="Per-cluster size, representative and mean metrics, best first",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
