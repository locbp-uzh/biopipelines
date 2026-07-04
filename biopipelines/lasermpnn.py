# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""LASErMPNN: ligand-conditioned inverse folding with all-atom sidechain packing.

LASErMPNN redesigns a protein sequence around a bound ligand and packs the
sidechains in a single pass, emitting full-atom complex PDBs (designed sequence
+ rotamers + ligand). The ligand is read directly from the input PDB's HETATM
records — there is no ligand-code argument.

This wrapper runs `python -m LASErMPNN.run_batch_inference` over all input
structures at once (fed as a path list so each input gets its own output
subdir), then collects the designs into the structures and sequences streams.
Each input yields `num_sequences` designs, surfaced as `<parent>_<1..N>`.

Fixed/redesigned positions are controlled through input-PDB B-factors with
`--fix_beta` (B=1 fixed, B=0 designed); when `fixed`/`redesigned` are given the
wrapper stamps B-factors onto a copy of each input PDB at runtime.

Reference:
    Repo: https://github.com/polizzilab/LASErMPNN
"""

import os
import json
from typing import Dict, List, Any, Optional, Union, Tuple

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import generate_multiplied_ids_pattern
    from .biopipelines_io import Resolve, TableReference
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_multiplied_ids_pattern
    from biopipelines_io import Resolve, TableReference


# Friendly model name -> weights filename shipped in the clone's model_weights/.
MODEL_WEIGHTS = {
    "default": "laser_weights_0p1A_nothing_heldout.pt",
    "ligandmpnn_split": "laser_weights_0p1A_noise_ligandmpnn_split.pt",
    "soluble": "soluble_weights_no_heldout_drop_clusters_optstep_65000.pt",
}


def _normalize_positions_arg(name, value):
    """Normalize a fixed/redesigned argument to a (kind, token) pair.

    Accepts a plain PyMOL selection string (broadcast to every input) or a
    per-input column reference (a ``TableReference`` from ``tool.tables.X.col``
    or a ``(TableInfo, "col")`` tuple), resolved by id match at runtime. Mirrors
    LigandMPNN's fixed/redesigned convention.

    Returns ``("literal", str)`` or ``("table", "TABLE_REFERENCE:path:col")``.
    """
    if value is None or value == "":
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


class LASErMPNN(BaseConfig):
    """LASErMPNN: ligand-conditioned inverse folding + all-atom packing."""

    TOOL_NAME = "LASErMPNN"
    TOOL_VERSION = "1.0"

    # ------------------------------------------------------------------
    # Install
    # ------------------------------------------------------------------

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False,
                        device="gpu", **kwargs):
        """Clone polizzilab/LASErMPNN and create the ``lasermpnn`` env from the
        device-matching vendored spec (``environments/lasermpnn.<device>.yaml`` +
        ``lasermpnn.<device>.pip.txt``, minimal conda + a pip layer with torch,
        the matching PyG extension wheels, prody, and pydssp).

        The repo ships its model weights in ``model_weights/`` and has **no**
        ``setup.py``; the notebook runs it as ``python -m LASErMPNN.run_batch_inference``
        from the clone's parent directory, so there is nothing to ``pip install -e``.
        Verification imports the package with the parent on ``PYTHONPATH``.

        device: "gpu" (default, torch+cu128) or "cpu" (torch+cpu). The env name
        is "lasermpnn" either way; device only selects which spec is installed.
        """
        if device not in ("gpu", "cpu"):
            raise ValueError(f"device must be 'gpu' or 'cpu', got {device!r}")
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("LASErMPNN", "")
        parent_dir = os.path.dirname(repo_dir)
        env_yaml = f"{biopipelines}/environments/lasermpnn.{device}.yaml"
        env_pip = f"{biopipelines}/environments/lasermpnn.{device}.pip.txt"

        env_check = cls._env_exists_check("lasermpnn", env_manager)
        repo_check = f'[ -d "{repo_dir}/.git" ]'
        weights_check = f'[ -f "{repo_dir}/model_weights/{MODEL_WEIGHTS["default"]}" ]'
        import_check = (
            f'PYTHONPATH="{parent_dir}" {env_manager} run -n lasermpnn '
            f'python -c "import LASErMPNN"'
        )
        # Env name is shared across devices, so the skip must also match the
        # recorded device marker or a cpu<->gpu switch would reuse the wrong env.
        device_marker = f"{repo_dir}/.laser_device"
        device_check = f'[ "$(cat "{device_marker}" 2>/dev/null)" = "{device}" ]'

        skip = "" if force_reinstall else f"""# Check if already installed (same device)
if {repo_check} && {weights_check} && {env_check} && {device_check} && {import_check} >/dev/null 2>&1; then
    echo "LASErMPNN ({device}) already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        # Unconditional (not just force_reinstall): the skip also declines on a
        # device switch, and env create fails on an existing env.
        remove_block = cls._env_remove_block("lasermpnn", env_manager)

        return f"""echo "=== Installing LASErMPNN ({device}) ==="
{skip}{remove_block}
mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}/.git" ]; then
    git clone https://github.com/polizzilab/LASErMPNN.git "{repo_dir}"
fi

# Create the lasermpnn env from the device-matching vendored spec + pip layer.
{env_manager} env create -f "{env_yaml}" -y
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create lasermpnn environment."
    exit 1
fi
{env_manager} run -n lasermpnn pip install -r "{env_pip}"
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to install lasermpnn pip layer."
    exit 1
fi

# No setup.py in the repo: the package is run as a module from the parent dir.
if {import_check} >/dev/null 2>&1; then
    echo "{device}" > "{device_marker}"
    touch "$INSTALL_SUCCESS"
    echo "=== LASErMPNN installation complete ==="
else
    echo "ERROR: LASErMPNN verification failed (cannot import LASErMPNN)"
    exit 1
fi
"""

    # ------------------------------------------------------------------
    # Lazy path descriptors
    # ------------------------------------------------------------------
    #   sequences_csv    — content-bearing sequences stream (map == table).
    #   summary/missing  — standalone TableInfo CSVs.
    #   structures_json  — serialized input DataStream (config-time artifact).
    #   positions_json   — per-input resolved {fixed, designed}, written at runtime.
    #   run_py           — pipe script: stamp B-factors, run inference, collect.
    #   positions_py     — resolver (one pass over all ids -> positions JSON).
    sequences_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    structures_json = Path(lambda self: self.configuration_path("input_structures.json"))
    positions_json = Path(lambda self: self.configuration_path("positions.json"))
    positions_args_json = Path(lambda self: self.configuration_path("positions_args.json"))
    run_py = Path(lambda self: self.pipe_script_path("pipe_lasermpnn.py"))
    positions_py = Path(lambda self: self.pipe_script_path("pipe_lasermpnn_bfactor.py"))
    repo_dir = Path(lambda self: self.folders["LASErMPNN"])

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 num_sequences: int = 1,
                 temperature: Optional[float] = None,
                 first_shell_temperature: Optional[float] = None,
                 chi_temperature: Optional[float] = None,
                 model: str = "default",
                 designs_per_batch: int = 30,
                 inputs_per_pass: int = 5,
                 disabled_residues: str = "X,C",
                 fixed: Union[str, Tuple['TableInfo', str], TableReference] = "",
                 redesigned: Union[str, Tuple['TableInfo', str], TableReference] = "",
                 chain: str = "A",
                 repack_only: bool = False,
                 repack_all: bool = False,
                 ignore_ligand: bool = False,
                 constrain_ala_gly: bool = False,
                 ala_budget: int = 4,
                 gly_budget: int = 0,
                 **kwargs):
        """Initialize LASErMPNN configuration.

        Args:
            structures: Input structures (StandardizedOutput or DataStream). PDBs
                        must carry the bound ligand as HETATM — LASErMPNN reads it
                        directly (no ligand-code argument).
            num_sequences: Designs to generate per input (upstream designs_per_input).
            temperature: Sequence sampling temperature (--sequence_temp). None (default)
                        uses the model's built-in temperature.
            first_shell_temperature: Temperature for binding-site residues
                        (--first_shell_sequence_temp). None (default) ties it to
                        the global temperature.
            chi_temperature: Rotamer (chi) sampling temperature (--chi_temp). None default.
            model: Weights to use — "default", "ligandmpnn_split", or "soluble".
            designs_per_batch: Designs per GPU batch (upstream default 30; tune for GPU memory).
            inputs_per_pass: Input files processed per GPU pass (upstream default 5).
            disabled_residues: Comma-separated one-letter codes to forbid
                        (upstream default "X,C" — omits Cys).
            fixed: Positions to hold at their input identity. PyMOL selection string
                        (broadcast) or a (table, column) reference (per-input). Stamped
                        as B-factor 1.0 with --fix_beta. Mutually exclusive with redesigned.
            redesigned: Positions to design; everything else on `chain` is held fixed.
                        Same accepted forms as `fixed`.
            chain: Default chain for chainless position input (default "A").
            repack_only: Keep the input sequence, only repack sidechains
                        (--repack_only_input_sequence).
            repack_all: Repack all residues, including fixed ones (--repack_all).
            ignore_ligand: Ignore the ligand during design (--ignore_ligand).
            constrain_ala_gly: Cap Ala/Gly over-sampling in exposed non-secondary-
                        structure regions (-c). Upstream flags this as generally helpful.
            ala_budget: Max Ala in the constrained region (upstream default 4).
            gly_budget: Max Gly in the constrained region (upstream default 0).

        Output:
            Streams:
                structures (.pdb) — full-atom designed complexes (seq + rotamers + ligand).
                sequences (csv)   — one row per design, extracted from designs.fasta.
            Tables:
                sequences: id | structures.id | sequence | score
                missing:   id | removed_by | kind | cause
        """
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(
                f"structures must be DataStream or StandardizedOutput, got {type(structures)}"
            )

        self.num_sequences = num_sequences
        self.temperature = temperature
        self.first_shell_temperature = first_shell_temperature
        self.chi_temperature = chi_temperature
        self.model = model
        self.designs_per_batch = designs_per_batch
        self.inputs_per_pass = inputs_per_pass
        self.disabled_residues = disabled_residues
        self.chain = chain
        self.repack_only = repack_only
        self.repack_all = repack_all
        self.ignore_ligand = ignore_ligand
        self.constrain_ala_gly = constrain_ala_gly
        self.ala_budget = ala_budget
        self.gly_budget = gly_budget

        self.fixed = fixed
        self.redesigned = redesigned
        self._fixed_arg = _normalize_positions_arg("fixed", fixed)
        self._redesigned_arg = _normalize_positions_arg("redesigned", redesigned)

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate LASErMPNN-specific parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures cannot be empty")

        # run_batch_inference reads PDB files (ligand from HETATM); reject non-PDB
        # or mixed-format streams before execution.
        if not self.structures_stream.has_only_formats("pdb"):
            raise ValueError(
                f"LASErMPNN requires PDB structures; got format "
                f"{self.structures_stream.format!r}"
            )

        if not isinstance(self.num_sequences, int) or self.num_sequences < 1:
            raise ValueError("num_sequences must be a positive integer")
        if not isinstance(self.designs_per_batch, int) or self.designs_per_batch < 1:
            raise ValueError("designs_per_batch must be a positive integer")
        if not isinstance(self.inputs_per_pass, int) or self.inputs_per_pass < 1:
            raise ValueError("inputs_per_pass must be a positive integer")

        for name, val in (("temperature", self.temperature),
                          ("first_shell_temperature", self.first_shell_temperature),
                          ("chi_temperature", self.chi_temperature)):
            if val is not None and (not isinstance(val, (int, float)) or val < 0):
                raise ValueError(f"{name} must be a non-negative number or None")

        if self.model not in MODEL_WEIGHTS:
            raise ValueError(f"model must be one of {sorted(MODEL_WEIGHTS)}, got {self.model!r}")

        if not isinstance(self.ala_budget, int) or self.ala_budget < 0:
            raise ValueError("ala_budget must be a non-negative integer")
        if not isinstance(self.gly_budget, int) or self.gly_budget < 0:
            raise ValueError("gly_budget must be a non-negative integer")

        if self._fixed_arg[1] and self._redesigned_arg[1]:
            raise ValueError("fixed and redesigned are mutually exclusive")

        # Only literal selections are interpolated raw into bash / stamped; table
        # references reach the resolver through a resolved subshell.
        _validate_freeform_string("disabled_residues", self.disabled_residues)
        _validate_freeform_string("chain", self.chain)
        if self._fixed_arg[0] == "literal" and self._fixed_arg[1]:
            _validate_freeform_string("fixed", self._fixed_arg[1])
        if self._redesigned_arg[0] == "literal" and self._redesigned_arg[1]:
            _validate_freeform_string("redesigned", self._redesigned_arg[1])

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input files."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get LASErMPNN configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"STRUCTURES: {', '.join(self.structures_stream.ids)}",
            f"NUM_SEQUENCES: {self.num_sequences}",
            f"MODEL: {self.model}",
            f"TEMPERATURE: {self.temperature if self.temperature is not None else 'model default'}",
            f"DISABLED_RESIDUES: {self.disabled_residues}",
        ])
        if self._fixed_arg[1]:
            config_lines.append(f"FIXED: {self.fixed}")
        if self._redesigned_arg[1]:
            config_lines.append(f"REDESIGNED: {self.redesigned}")
        if self.constrain_ala_gly:
            config_lines.append(f"ALA/GLY BUDGET: {self.ala_budget}/{self.gly_budget}")
        return config_lines

    # ------------------------------------------------------------------
    # Script generation
    # ------------------------------------------------------------------

    def generate_script(self, script_path: str) -> str:
        """Generate LASErMPNN execution script."""
        self.structures_stream.save_json(self.structures_json)

        has_positions = bool(self._fixed_arg[1] or self._redesigned_arg[1])

        script_content = "#!/bin/bash\n"
        script_content += "# LASErMPNN execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_device_detection()
        if has_positions:
            script_content += self._generate_resolve_positions()
        script_content += self._generate_run_and_collect(has_positions)
        script_content += self.generate_completion_check_footer()
        return script_content

    def _generate_device_detection(self) -> str:
        """Pick a torch device string at runtime (upstream defaults to cpu).

        run_batch_inference defaults --device to 'cpu', so without this it would
        run on CPU even on a GPU node. Probe torch in the activated env.
        """
        return """echo "Detecting compute device"
LASER_DEVICE=$(python -c "import torch; print('cuda:0' if torch.cuda.is_available() else 'cpu')")
echo "LASErMPNN device: $LASER_DEVICE"

"""

    def _generate_resolve_positions(self) -> str:
        """Resolve fixed/redesigned selections per input into positions.json."""
        with open(self.positions_args_json, "w") as f:
            json.dump({
                "structures_json": str(self.structures_json),
                "fixed": self._fixed_arg[1] or "-",
                "redesigned": self._redesigned_arg[1] or "-",
                "default_chain": self.chain,
                "output_json": str(self.positions_json),
            }, f, indent=2)

        return f"""echo "Resolving fixed/redesigned positions"
python {self.positions_py} resolve "{self.positions_args_json}"

"""

    def _run_options(self) -> str:
        """run_batch_inference flags shared by every input (no positionals/device).

        Values are unquoted: the pipe script shlex-splits this string into an
        argv list for subprocess.run (no shell), so embedding bash quotes here
        would pass literal quote characters through to the binary.
        """
        weights = os.path.join(self.repo_dir, "model_weights", MODEL_WEIGHTS[self.model])
        opts = (
            f"--model_weights_path {weights}"
            f" --designs_per_batch {self.designs_per_batch}"
            f" --inputs_processed_simultaneously {self.inputs_per_pass}"
            f" --disabled_residues {self.disabled_residues}"
            f" --output_fasta"
        )
        if self.temperature is not None:
            opts += f" --sequence_temp {self.temperature}"
        if self.first_shell_temperature is not None:
            opts += f" --first_shell_sequence_temp {self.first_shell_temperature}"
        if self.chi_temperature is not None:
            opts += f" --chi_temp {self.chi_temperature}"
        if self.repack_only:
            opts += " --repack_only_input_sequence"
        if self.repack_all:
            opts += " --repack_all"
        if self.ignore_ligand:
            opts += " --ignore_ligand"
        if self.constrain_ala_gly:
            opts += (
                " -c"
                f" --ala_budget {self.ala_budget}"
                f" --gly_budget {self.gly_budget}"
            )
        return opts

    def _generate_run_and_collect(self, has_positions: bool) -> str:
        """Stamp B-factors (if any), run inference, and collect designs.

        The whole thing runs in one pipe script: it prepares an input list (the
        original PDBs, or B-factor-stamped copies), invokes run_batch_inference
        as a module from the clone's parent dir (the repo has no setup.py), then
        renames the per-input ``design_{n}.pdb`` outputs into the structures
        stream and builds the sequences/missing tables from ``designs.fasta``.
        """
        parent_dir = os.path.dirname(self.repo_dir)
        exec_root = self.execution_path()
        structures_dir = self.stream_folder("structures")
        structures_map = self.stream_map_path("structures")
        positions_arg = f' --positions-json "{self.positions_json}"' if has_positions else ""
        fix_beta_flag = " --fix-beta" if has_positions else ""

        return f"""echo "Running LASErMPNN and collecting designs"
python {self.run_py} \\
    --structures-json "{self.structures_json}" \\
    --repo-parent "{parent_dir}" \\
    --exec-root "{exec_root}" \\
    --structures-dir "{structures_dir}" \\
    --structures-map "{structures_map}" \\
    --sequences-csv "{self.sequences_csv}" \\
    --missing-csv "{self.missing_csv}" \\
    --num-sequences {self.num_sequences} \\
    --device "$LASER_DEVICE" \\
    --run-options {self._quote(self._run_options())}{positions_arg}{fix_beta_flag}

"""

    @staticmethod
    def _quote(s: str) -> str:
        """Single-quote a bash argument, escaping embedded single quotes."""
        return "'" + s.replace("'", "'\\''") + "'"

    # ------------------------------------------------------------------
    # Output prediction
    # ------------------------------------------------------------------

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after LASErMPNN execution."""
        # <parent>_<1..num_sequences> fan-out via the shared multiplier helper.
        suffix_pattern = f"<1..{self.num_sequences}>"
        design_ids = generate_multiplied_ids_pattern(
            self.structures_stream.ids, suffix_pattern,
            input_stream_name="structures"
        )
        structures_dir = self.stream_folder("structures")
        file_template = [os.path.join(structures_dir, "<id>.pdb")]

        structures = DataStream(
            name="structures",
            ids=design_ids,
            files=file_template,
            map_table=self.stream_map_path("structures"),
            format="pdb"
        )

        # Content-bearing stream (map_table == TableInfo path); the run script
        # extracts each design's sequence from designs.fasta.
        sequences = DataStream(
            name="sequences",
            ids=design_ids,
            files=[],
            map_table=self.sequences_csv,
            format="csv"
        )

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "structures.id", "sequence", "score"],
                description="LASErMPNN designed sequences with log-probability score"
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Designs not produced (fewer valid designs than num_sequences)"
            )
        }

        return {
            "structures": structures,
            "sequences": sequences,
            "tables": tables,
            "output_folder": self.output_folder
        }
