# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""EnsembleAnalysis: per-residue RMSF + ensemble-level metrics from a
conformer ensemble.

Given a structural ensemble, the tool superposes all conformers (least-squares
on a chosen atom set) and reports per-residue fluctuation. It fills the gap left
by CABSflex (the only other RMSF producer), which is coupled to its own
coarse-grained sampling: EnsembleAnalysis works on *any* ensemble, regardless of
where the conformers came from.

The ensemble partition is always defined explicitly by a required ``groups=``
stream (the Consensus convention): its ids ARE the output ensemble ids, known at
config time. Each structures-stream id is matched to a group and its model(s)
pooled into that group's ensemble. The two input shapes compose freely under
this rule:

  - Multi-model file: one structures-stream id whose file holds several MODEL
    records (NMR ensembles, PLACER dumps). Map it to a group of its own (a
    one-id groups stream) to treat the file as a single ensemble.
  - Conformer set: a structures stream of single-model files (e.g. BioEmu's
    per-conformer PDBs `<seq>_1..N`); pass ``groups=`` the upstream sampler's
    input (e.g. the folded sequence) so the pooled ensemble id is the parent.

Outputs a `resi-csv` rmsf stream (`id, chain, resi, rmsf, rmsd_mean`) matching
CABSflex's schema, so `Selection` can threshold on it unchanged
(`Selection.add(ens.streams.rmsf, include="rmsf<=1.0")`), plus two standalone
tables: per-ensemble summary (mean/max RMSF, radius-of-gyration spread) and
per-frame metrics (RMSD-to-reference, Rg).

Runs in the `biopipelines` env (numpy + the shared pdb_parser); no external tool,
no dedicated env.
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


# Atom sets the superposition + RMSF may be computed over. CA is the standard
# RMSF convention (matches CABSflex's Cα profile); backbone widens to N/CA/C/O.
_SELECTIONS = ("CA", "backbone")
_REFERENCES = ("mean", "first")


class EnsembleAnalysis(BaseConfig):
    """
    EnsembleAnalysis: per-residue RMSF and ensemble metrics from a conformer
    ensemble.

    Inputs:
        structures: the ensemble, as a structures stream or StandardizedOutput.
                    Each id's file may hold one model or several MODEL records;
                    all of a group's members' models are pooled into the
                    group's ensemble.
        groups:     required stream whose ids define the ensemble partition
                    (Consensus convention). Input ids are matched to these group
                    ids and pooled per group; the output ensemble ids are the
                    group ids. Pass the upstream sampler's input (e.g. the
                    folded sequence) so the collapsed ensemble id is known up
                    front; for a single multi-model file, pass a one-id stream.
        selection:  atom set for superposition and RMSF — "CA" (default) or
                    "backbone" (N, CA, C, O).
        reference:  superposition reference — "mean" (default; iterative
                    align-to-mean) or "first" (align all frames to the first).

    Outputs:
        Streams:
            rmsf:  per-residue resi-csv, one CSV per ensemble id.
                   Columns: id, chain, resi, rmsf, rmsd_mean.
        Tables:
            residues: id | chain | resi | rmsf | rmsd_mean  (merged over all
                      ensembles; the per-residue rmsf stream, concatenated).
            ensemble: id | n_frames | n_residues | mean_rmsf | max_rmsf |
                      rg_mean | rg_std  (one row per ensemble)
            frames:   id | frame | rmsd_to_ref | rg  (one row per conformer)
    """

    TOOL_NAME = "EnsembleAnalysis"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """No dedicated env: runs in the biopipelines env (numpy + pdb_parser)."""
        return """echo "=== EnsembleAnalysis ==="
echo "Uses the biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== EnsembleAnalysis ready ==="
"""

    # Canonical layout.
    #   _configuration/ — input structures DataStream JSON.
    #   rmsf/           — per-ensemble resi-csv + rmsf_map.
    #   tables/         — ensemble (summary) + frames (per-conformer).
    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    groups_json = Path(lambda self: self.configuration_path(".groups.json"))
    rmsf_map = Path(lambda self: self.stream_map_path("rmsf"))
    residues_csv = Path(lambda self: self.table_path("residues"))
    ensemble_csv = Path(lambda self: self.table_path("ensemble"))
    frames_csv = Path(lambda self: self.table_path("frames"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_ensemble_analysis.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 groups: Union[DataStream, StandardizedOutput],
                 selection: str = "CA",
                 reference: str = "mean",
                 **kwargs):
        # Keep original inputs for upstream missing-table detection.
        self.structures_input = structures
        self.groups_input = groups

        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(
                f"structures must be DataStream or StandardizedOutput, got {type(structures)}"
            )

        # Explicit grouping (Consensus convention): the groups stream's ids ARE
        # the output ensemble ids, known at config time, so the declared ids
        # always match what the pipe script writes. Each input id is matched to
        # a group and its model(s) pooled into that group's ensemble.
        self.groups_stream: DataStream = self._resolve_group_stream(groups)

        self.selection = selection
        self.reference = reference
        super().__init__(**kwargs)

    @staticmethod
    def _resolve_group_stream(obj) -> DataStream:
        if isinstance(obj, DataStream):
            return obj
        if isinstance(obj, StandardizedOutput):
            candidates = [ds for _, ds in obj.streams.items()
                          if ds is not None and len(ds) > 0]
            if not candidates:
                raise ValueError("groups: no non-empty stream found in output")
            # The group ids define the partition; streams sharing one id set
            # yield the same grouping, so the pick is harmless. Differing id
            # sets would change the result, so refuse to guess.
            if len({tuple(ds.ids) for ds in candidates}) > 1:
                raise ValueError(
                    f"groups: output streams disagree on ids "
                    f"({[(c.name, len(c)) for c in candidates]}); pass an explicit "
                    f"stream (e.g. groups=tool.streams.structures)"
                )
            return candidates[0]
        raise ValueError(f"groups must be a DataStream or StandardizedOutput, got {type(obj)}")

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.groups_stream or len(self.groups_stream) == 0:
            raise ValueError("groups is required and must not be empty")
        if self.selection not in _SELECTIONS:
            raise ValueError(
                f"selection must be one of {_SELECTIONS}, got: {self.selection!r}"
            )
        if self.reference not in _REFERENCES:
            raise ValueError(
                f"reference must be one of {_REFERENCES}, got: {self.reference!r}"
            )

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"INPUT STRUCTURES: {len(self.structures_stream)}")
        lines.append(f"GROUPS: {len(self.groups_stream)} ids")
        lines.append(f"SELECTION: {self.selection}")
        lines.append(f"REFERENCE: {self.reference}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self.groups_stream.save_json(self.groups_json)
        groups_arg = f' \\\n    --groups-json "{self.groups_json}"'

        script = "#!/bin/bash\n"
        script += "# EnsembleAnalysis: per-residue RMSF + ensemble metrics\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Analyzing ensemble flexibility for {len(self.structures_stream)} input(s)"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --selection "{self.selection}" \\
    --reference "{self.reference}" \\
    --rmsf-dir "{self.stream_folder('rmsf')}" \\
    --rmsf-map "{self.rmsf_map}" \\
    --residues-csv "{self.residues_csv}" \\
    --ensemble-csv "{self.ensemble_csv}" \\
    --frames-csv "{self.frames_csv}"{groups_arg}
"""
        # A group whose members were all dropped upstream produces no <id>.csv;
        # propagate the upstream manifest so the completion check excuses it.
        script += self.generate_missing_propagation(
            self.structures_input, self.groups_input, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        # One ensemble per group; the ensemble ids ARE the groups stream ids
        # (Consensus convention), known at config time, so the declared ids
        # match what the pipe script writes.
        rmsf_stream = DataStream(
            name="rmsf",
            ids=self.groups_stream.ids,
            files=[self.stream_path("rmsf", "<id>.csv")],
            map_table=self.rmsf_map,
            format="resi-csv",
        )
        tables = {
            "residues": TableInfo(
                name="residues",
                path=self.residues_csv,
                columns=["id", "chain", "resi", "rmsf", "rmsd_mean"],
                description="Per-residue RMSF + RMSD-to-mean, merged over all ensembles",
            ),
            "ensemble": TableInfo(
                name="ensemble",
                path=self.ensemble_csv,
                columns=["id", "n_frames", "n_residues", "mean_rmsf",
                         "max_rmsf", "rg_mean", "rg_std"],
                description="Per-ensemble flexibility summary",
            ),
            "frames": TableInfo(
                name="frames",
                path=self.frames_csv,
                columns=["id", "frame", "rmsd_to_ref", "rg"],
                description="Per-conformer RMSD-to-reference and radius of gyration",
            ),
        }

        # Excuse groups whose members were all dropped upstream: declare a
        # `missing` table when any input axis carries one, so the completion
        # check doesn't fail the never-written <group>.csv.
        if self._collect_upstream_missing_paths(self.structures_input, self.groups_input):
            tables["missing"] = self.missing_table_info(self.missing_csv)

        return {
            "rmsf": rmsf_stream,
            "tables": tables,
            "output_folder": self.output_folder,
        }
