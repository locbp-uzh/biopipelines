# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""HBDesigner tool: design hydrogen-bonding networks onto a protein backbone.

HBDesigner places highly-connected hydrogen-bonding networks that satisfy the
requested design constraints onto an input backbone. For each input structure it
samples candidate networks, packs and scores them with PyRosetta, and retains the
top-ranked designs.

This wrapper runs HBDesigner once per input structure (a runtime loop over the
input ids), so an upstream backbone ensemble propagates fully into design. Each
run keeps `top_k` ranked designs, surfaced as `<parent>_<1..top_k>`.

The network is steered by three optional constraints, each a plain string
(broadcast to every input) or a per-input table column reference (resolved by id
match at runtime):
  * guide_res:   residues to center the design around (PDB chain/resnum, e.g. "A12,B13").
  * guide_seq:   amino-acid types to constrain the network ("X" leaves a position free).
  * anchor_res:  residues that must appear in every designed network.

Reference:
    Repo: https://github.com/Kuhlman-Lab/HBDesigner
    Docs: https://rosettacommons.github.io/HBDesigner/
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import Resolve, TableReference
    from .combinatorics import generate_multiplied_ids_pattern
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import Resolve, TableReference
    from combinatorics import generate_multiplied_ids_pattern


def _normalize_constraint_arg(name, value):
    """Normalize a per-input constraint argument to a (kind, token) pair.

    Each constraint (`guide_res` / `guide_seq` / `anchor_res`) accepts either a
    plain string — broadcast to every input structure — or a per-input column
    reference (a ``TableReference`` from ``tool.tables.X.col``, or a
    ``(TableInfo, "col")`` tuple), resolved by id match at runtime. Mirrors
    RFdiffusion's contigs convention.

    Returns ``("literal", str)`` or ``("table", "TABLE_REFERENCE:path:col")``.
    Only a literal is validated by ``_validate_freeform_string``; a table
    reference reaches bash through a resolved subshell, not raw interpolation.
    """
    if value is None:
        return ("literal", "")
    if isinstance(value, TableReference):
        return ("table", str(value))
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], str):
        return ("table", str(TableReference(value[0].info.path, value[1])))
    if isinstance(value, str):
        return ("literal", value)
    raise ValueError(
        f"{name} must be a string or a (table, column) reference, got {type(value)}"
    )


class HBDesigner(BaseConfig):
    """
    HBDesigner: design hydrogen-bonding networks onto an input backbone.

    Inputs:
        structures:  backbone PDBs (StandardizedOutput or DataStream). HBDesigner
                     runs once per input structure.
        n_res:       hydrogen-bonding network size, 2-6 residues (default 2).
        n_samples:   unique samples generated before packing/scoring (default 100).
        top_k:       designs retained per input after ranking, surfaced as
                     ``<parent>_<1..top_k>`` (default 5). Use top_k=1 for the
                     single best design per input.
        design_model: "design_020" (moderate noise, default) or "design_002"
                     (low noise).
        guide_res / guide_seq / anchor_res: optional design constraints, each a
                     string (broadcast) or a (table, column) reference (per-input).
        seed:        random seed for reproducibility (default None).

    Outputs:
        Streams:
            structures (.pdb) — the designed HBNet backbones.
            sequences (csv)   — one-letter sequence per designed backbone
                                (HBDesigner emits no FASTA; extracted from the PDB).
        Tables:
            sequences: id | structures.id | sequence
            summary: id | structures.id | Scaffold | Output_PDB | Rank |
                     HB_Score_full | HB_Score_hb | Avg_Burial | saturation |
                     buried_heavy_unsats | buried_unsat_Hpol | network
                     (one row per retained design; the metric columns mirror
                     HBDesigner's per-run stats CSV)
            missing: id | removed_by | kind | cause (ranks not produced when a
                     run yields fewer valid designs than top_k)
    """

    TOOL_NAME = "HBDesigner"
    TOOL_VERSION = "1.1"

    DESIGN_MODELS = ("design_002", "design_020")

    # ------------------------------------------------------------------
    # Install
    # ------------------------------------------------------------------

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False,
                        device="gpu", **kwargs):
        """Clone Kuhlman-Lab/HBDesigner, create the ``hbdesigner`` env from the
        device-matching vendored spec (``environments/hbdesigner.<device>.yaml``
        + ``hbdesigner.<device>.pip.txt``, vendored from upstream env_gpu.yaml /
        env_cpu.yaml), then ``pip install -e .`` the HBDesigner package from the
        clone (editable: run_hbdesigner resolves its model_weights/ relative to
        the package location, which must be the repo). Verification: the
        run_hbdesigner entry point resolves inside the env.

        device: "gpu" (default, torch+cu128) or "cpu" (torch+cpu). The env name
        is "hbdesigner" either way; device only selects which spec is installed.

        HBDesigner depends on PyRosetta, which is free for academic use but
        requires a paid commercial license; the user accepts that license by
        installing.
        """
        if device not in ("gpu", "cpu"):
            raise ValueError(f"device must be 'gpu' or 'cpu', got {device!r}")
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("HBDesigner", "")
        parent_dir = os.path.dirname(repo_dir)
        env_yaml = f"{biopipelines}/environments/hbdesigner.{device}.yaml"
        env_pip = f"{biopipelines}/environments/hbdesigner.{device}.pip.txt"

        env_check = cls._env_exists_check("hbdesigner", env_manager)
        repo_check = f'[ -d "{repo_dir}/.git" ]'
        entry_check = (
            f'{env_manager} run -n hbdesigner '
            f'python -c "import shutil,sys; sys.exit(0 if shutil.which(\'run_hbdesigner\') else 1)"'
        )
        # Env name is shared across devices, so the skip must also match the
        # recorded device marker or a cpu<->gpu switch would reuse the wrong env.
        device_marker = f"{repo_dir}/.hbdes_device"
        device_check = f'[ "$(cat "{device_marker}" 2>/dev/null)" = "{device}" ]'

        skip = "" if force_reinstall else f"""# Check if already installed (same device)
if {repo_check} && {env_check} && {device_check} && {entry_check} >/dev/null 2>&1; then
    echo "HBDesigner ({device}) already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        # Unconditional (not just force_reinstall): skip also declines on a
        # device switch, and env create fails on an existing env.
        remove_block = cls._env_remove_block("hbdesigner", env_manager)

        return f"""echo "=== Installing HBDesigner ({device}) ==="
{skip}{remove_block}
mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}/.git" ]; then
    git clone https://github.com/Kuhlman-Lab/HBDesigner.git "{repo_dir}"
fi

# Create the hbdesigner env from the device-matching vendored spec + pip layer.
{env_manager} env create -f "{env_yaml}" -y
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create hbdesigner environment."
    exit 1
fi
{env_manager} run -n hbdesigner pip install -r "{env_pip}"
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to install hbdesigner pip layer."
    exit 1
fi

# Editable: run_hbdesigner resolves model_weights/ relative to the package, which must be the repo.
cd "{repo_dir}"
{env_manager} run -n hbdesigner pip install -e .

# Verify: the run_hbdesigner entry point resolves inside the env.
if {entry_check} >/dev/null 2>&1; then
    echo "{device}" > "{device_marker}"
    touch "$INSTALL_SUCCESS"
    echo "=== HBDesigner installation complete ==="
else
    echo "ERROR: HBDesigner verification failed (run_hbdesigner not found in env)"
    exit 1
fi
"""

    # ------------------------------------------------------------------
    # Lazy path descriptors
    # ------------------------------------------------------------------
    #   summary_table       — TableInfo CSV of per-design metrics.
    #   structures_json     — serialized input DataStream (config-time artifact).
    #   constraints_json    — per-input resolved {guide_res, guide_seq,
    #                          anchor_res}, written at runtime by the resolver.
    #   collect_py          — pipe script that runs run_hbdesigner per input,
    #                          collects designed PDBs into the structures stream,
    #                          and builds the summary table.
    #   constraints_py      — resolver (one pass over all ids -> JSON).
    #   constraints_reader_py — per-id reader the bash loop calls per input.
    #   update_map_py       — shared script that writes structures_map.csv from
    #                          the actual output files.
    summary_table = Path(lambda self: self.table_path("summary"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    sequences_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))
    structures_json = Path(lambda self: self.configuration_path("input_structures.json"))
    constraints_json = Path(lambda self: self.configuration_path("constraints.json"))
    constraints_args_json = Path(lambda self: self.configuration_path("constraints_args.json"))
    collect_py = Path(lambda self: self.pipe_script_path("pipe_hbdesigner.py"))
    constraints_py = Path(lambda self: self.pipe_script_path("pipe_hbdesigner_constraints.py"))
    constraints_reader_py = Path(lambda self: self.pipe_script_path("resolve_hbdesigner_constraints.py"))
    update_map_py = Path(lambda self: self.pipe_script_path("pipe_update_structures_map.py"))

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 n_res: int = 2,
                 n_samples: int = 100,
                 top_k: int = 5,
                 design_model: str = "design_020",
                 guide_res: Optional[Union[str, TableReference, "tuple"]] = None,
                 guide_seq: Optional[Union[str, TableReference, "tuple"]] = None,
                 anchor_res: Optional[Union[str, TableReference, "tuple"]] = None,
                 seed: Optional[int] = None,
                 **kwargs):
        # Resolve input structures — iterated at execution time (never ids[0] at
        # config time, which breaks under lazy IDs).
        if isinstance(structures, StandardizedOutput):
            self.structures_stream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(
                f"structures must be DataStream or StandardizedOutput, got {type(structures)}"
            )

        self.n_res = n_res
        self.n_samples = n_samples
        self.top_k = top_k
        self.design_model = design_model
        self.seed = seed

        # Literal -> broadcast; table ref -> resolved per-id at runtime. Keep originals for display.
        self.guide_res = guide_res
        self.guide_seq = guide_seq
        self.anchor_res = anchor_res
        self._guide_res_arg = _normalize_constraint_arg("guide_res", guide_res)
        self._guide_seq_arg = _normalize_constraint_arg("guide_seq", guide_seq)
        self._anchor_res_arg = _normalize_constraint_arg("anchor_res", anchor_res)

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate HBDesigner-specific parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures cannot be empty")

        # run_hbdesigner reads PDB files specifically (--pdb); reject non-PDB or
        # mixed-format streams before execution rather than failing at runtime.
        if not self.structures_stream.has_only_formats("pdb"):
            raise ValueError(
                f"HBDesigner requires PDB structures; got format "
                f"{self.structures_stream.format!r}"
            )

        if not isinstance(self.n_res, int) or not (2 <= self.n_res <= 6):
            raise ValueError("n_res must be an integer in 2-6")
        if not isinstance(self.n_samples, int) or self.n_samples < 1:
            raise ValueError("n_samples must be a positive integer")
        if not isinstance(self.top_k, int) or self.top_k < 1:
            raise ValueError("top_k must be a positive integer")
        if self.design_model not in self.DESIGN_MODELS:
            raise ValueError(
                f"design_model must be one of {self.DESIGN_MODELS}, got {self.design_model!r}"
            )
        if self.seed is not None and not isinstance(self.seed, int):
            raise ValueError("seed must be an integer or None")

        # Only literal constraints are interpolated raw into bash; table
        # references reach the shell through a resolved subshell.
        for name, (kind, token) in (("guide_res", self._guide_res_arg),
                                     ("guide_seq", self._guide_seq_arg),
                                     ("anchor_res", self._anchor_res_arg)):
            if kind == "literal" and token:
                _validate_freeform_string(name, token)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get HBDesigner configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"STRUCTURES: {', '.join(self.structures_stream.ids)}",
            f"N_RES: {self.n_res}",
            f"N_SAMPLES: {self.n_samples}",
            f"TOP_K: {self.top_k}",
            f"DESIGN_MODEL: {self.design_model}",
        ])
        if self._guide_res_arg[1]:
            config_lines.append(f"GUIDE_RES: {self.guide_res}")
        if self._guide_seq_arg[1]:
            config_lines.append(f"GUIDE_SEQ: {self.guide_seq}")
        if self._anchor_res_arg[1]:
            config_lines.append(f"ANCHOR_RES: {self.anchor_res}")
        if self.seed is not None:
            config_lines.append(f"SEED: {self.seed}")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate HBDesigner execution script."""
        self.structures_stream.save_json(self.structures_json)

        script_content = "#!/bin/bash\n"
        script_content += "# HBDesigner execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_hbdesigner()
        script_content += self._generate_script_collect()
        script_content += self._generate_script_update_structures_map()
        script_content += self.generate_completion_check_footer()
        return script_content

    def _fixed_options(self) -> str:
        """run_hbdesigner flags shared by every input (no --pdb / --out_dir)."""
        opts = (
            f"--n_res {self.n_res}"
            f" --n_samples {self.n_samples}"
            f" --top_k {self.top_k}"
            f" --design_model {self.design_model}"
        )
        if self.seed is not None:
            opts += f" --seed {self.seed}"
        return opts

    def _generate_script_run_hbdesigner(self) -> str:
        """Run HBDesigner once per input structure, into a per-id out_dir.

        Per-input constraints (literal broadcast or table-column reference) are
        pre-resolved in one pass into constraints.json; the loop reads each id's
        values and appends only the non-empty ones as flags. Each run writes its
        designed PDBs and a CSV summary under <execution>/<id>/.
        """
        fixed = self._fixed_options()
        exec_root = self.execution_path()

        with open(self.constraints_args_json, "w") as f:
            json.dump({
                "structures_json": str(self.structures_json),
                "guide_res": self._guide_res_arg[1] or "-",
                "guide_seq": self._guide_seq_arg[1] or "-",
                "anchor_res": self._anchor_res_arg[1] or "-",
                "output_json": str(self.constraints_json),
            }, f, indent=2)

        return f"""echo "Resolving per-input constraints"
python {self.constraints_py} "{self.constraints_args_json}"

echo "Starting HBDesigner"
echo "Output folder: {self.output_folder}"

for STRUCT_ID in {Resolve.stream_ids(self.structures_json)}; do
    INPUT_PDB={Resolve.stream_item(self.structures_json, '$STRUCT_ID')}
    OUT_DIR="{exec_root}/$STRUCT_ID"
    mkdir -p "$OUT_DIR"
    # Clear any stale crash marker from a prior attempt so it reflects only this run.
    rm -f "$OUT_DIR/.hbdes_failed"
    OPTS=$(python {self.constraints_reader_py} "{self.constraints_json}" "$STRUCT_ID")
    GUIDE_RES=$(echo "$OPTS" | sed -n '1p')
    GUIDE_SEQ=$(echo "$OPTS" | sed -n '2p')
    ANCHOR_RES=$(echo "$OPTS" | sed -n '3p')
    HBDES_OPTIONS="--pdb \\"$INPUT_PDB\\" --out_dir \\"$OUT_DIR\\" {fixed}"
    if [ -n "$GUIDE_RES" ]; then
        HBDES_OPTIONS="$HBDES_OPTIONS --guide_res \\"$GUIDE_RES\\""
    fi
    if [ -n "$GUIDE_SEQ" ]; then
        HBDES_OPTIONS="$HBDES_OPTIONS --guide_seq \\"$GUIDE_SEQ\\""
    fi
    if [ -n "$ANCHOR_RES" ]; then
        HBDES_OPTIONS="$HBDES_OPTIONS --anchor_res \\"$ANCHOR_RES\\""
    fi
    echo "Running HBDesigner for $STRUCT_ID"
    if ! eval run_hbdesigner $HBDES_OPTIONS; then
        # Record a crash so the collector reports a real failure for this input
        # (kind="failure") rather than excusing its absent output as a filter.
        echo "ERROR: run_hbdesigner failed for $STRUCT_ID" >&2
        touch "$OUT_DIR/.hbdes_failed"
    fi
done

"""

    def _generate_script_collect(self) -> str:
        """Collect per-input designs into the structures + sequences streams and the summary table."""
        structures_dir = self.stream_folder("structures")
        exec_root = self.execution_path()
        return f"""echo "Collecting HBDesigner designs"
python {self.collect_py} \\
    --structures-json "{self.structures_json}" \\
    --exec-root "{exec_root}" \\
    --structures-dir "{structures_dir}" \\
    --summary-csv "{self.summary_table}" \\
    --sequences-csv "{self.sequences_csv}" \\
    --missing-csv "{self.missing_csv}" \\
    --top-k {self.top_k}

"""

    def _generate_script_update_structures_map(self) -> str:
        """Write structures_map.csv from the actual designed PDBs.

        Each id is ``<parent>_<n>``; recover the ``structures.id`` parent from
        the suffix.
        """
        structures_map = self.stream_map_path("structures")
        structures_dir = self.stream_folder("structures")
        prov_arg = ' --provenance-from-suffix "structures.id"'
        return f"""echo "Writing structures map from actual output files"
python {self.update_map_py} --structures-map "{structures_map}" --output-folder "{structures_dir}"{prov_arg}

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after HBDesigner execution."""
        # <parent>_<1..top_k> fan-out via the shared multiplier helper (compact, lazy-safe).
        suffix_pattern = f"<1..{self.top_k}>"
        structure_ids = generate_multiplied_ids_pattern(
            self.structures_stream.ids, suffix_pattern,
            input_stream_name="structures"
        )
        structures_dir = self.stream_folder("structures")
        file_template = [os.path.join(structures_dir, "<id>.pdb")]
        structures_map = self.stream_map_path("structures")

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=file_template,
            map_table=structures_map,
            format="pdb"
        )

        # Content-bearing stream (map_table == TableInfo path); the collect script
        # extracts the sequence from each PDB since HBDesigner emits no FASTA.
        sequences = DataStream(
            name="sequences",
            ids=structure_ids,
            files=[],
            map_table=self.sequences_csv,
            format="csv"
        )

        # id + structures.id are framework columns; the rest mirror HBDesigner's
        # per-run stats CSV (carried through verbatim by the collect script).
        tables = {
            "summary": TableInfo(
                name="summary",
                path=self.summary_table,
                columns=["id", "structures.id", "Scaffold", "Output_PDB", "Rank",
                         "HB_Score_full", "HB_Score_hb", "Avg_Burial", "saturation",
                         "buried_heavy_unsats", "buried_unsat_Hpol", "network"],
                description="HBDesigner per-design metrics"
            ),
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "structures.id", "sequence"],
                description="Sequences of the designed HBNet backbones"
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Ranks not produced (fewer valid designs than top_k)"
            )
        }

        return {
            "structures": structures,
            "sequences": sequences,
            "tables": tables,
            "output_folder": self.output_folder
        }
