# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

import os
import re
from typing import Dict, List, Any, Union, Optional

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import generate_multiplied_ids_pattern
    from .chain_rows import chain_row_ids, normalize_chains, validate_chains
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_multiplied_ids_pattern
    from chain_rows import chain_row_ids, normalize_chains, validate_chains


MODEL_TYPES = ("protein_mpnn", "ligand_mpnn")

# `<chain><residue number>[insertion code]`, the id spelling upstream parses.
RESIDUE_ID = re.compile(r"^[A-Za-z]+\d+[A-Za-z]*$")
CHAIN_ID = re.compile(r"^[A-Za-z]+$")
TOKEN = re.compile(r"^[A-Za-z0-9]{1,4}$")


def as_id_list(value) -> Optional[List[str]]:
    """A selection given as "A35,B40" or as a list, as a list; None stays None."""
    if value is None:
        return None
    if isinstance(value, str):
        return [part.strip() for part in value.split(",") if part.strip()]
    return list(value)


class FoundryMPNN(BaseConfig):
    """
    FoundryMPNN: ProteinMPNN / LigandMPNN sequence design from the RosettaCommons Foundry package.

    Foundry (https://github.com/RosettaCommons/foundry) ships a re-implementation of ProteinMPNN and
    LigandMPNN on AtomWorks and runs it with its `mpnn` command. It reads the original weights, which
    `FoundryMPNN.install()` downloads. It shares the `foundry` environment with `RFdiffusion3` and `RF3`.

    Which MPNN wrapper to use. This is the same model family as `ProteinMPNN` and `LigandMPNN`, not a new
    model, so the choice is about the code, not the biology:
      * Upstream's README asks that the original repositories be used for benchmarking and comparison until
        its re-implementation is validated, and calls its inference API unstable. `ProteinMPNN` and
        `LigandMPNN` wrap those original repositories, so they are the ones to cite or compare against.
      * Those two also offer what this one does not: per-position `fixed` / `redesigned` selections as PyMOL
        strings or per-structure table columns, a bound ligand named by a compounds stream with a
        `design_within` radius, side-chain packing, duplicate removal and a bias file.
      * This one takes the non-protein context of a CIF or PDB (ligands, ions, nucleic acids) from the file
        itself through AtomWorks, with no `ligand=` stream, designs scoped by `<chain><number>` ids as
        upstream spells them, and needs no environment beyond the one RFdiffusion3 and RF3 already use.
    Choose `model_type="ligand_mpnn"` when that non-protein context should condition the design and
    `"protein_mpnn"` for a bare backbone.

    Inputs:
        structures: backbones, PDB or mmCIF.
        model_type: "protein_mpnn" (default) or "ligand_mpnn".
        num_sequences: sequences sampled per batch (upstream `batch_size`).
        num_batches: batches per structure (upstream `number_of_batches`), so a structure yields
            `num_sequences * num_batches` designs.
        temperature: sampling temperature (upstream default 0.1, which is also the default here).
        seed: random seed, set per structure; unseeded when None.
        structure_noise: Gaussian noise in Angstrom added to the backbone before design (upstream default 0).
        designed_residues / fixed_residues: residue ids to redesign or to keep, as a list or "A35,B40,C52".
        designed_chains / fixed_chains: chain ids to redesign or to keep, as a list or "A,B".
            Upstream accepts one of the four per structure and refuses a mix, and so does this.
            Every structure receives the same selection, so a residue the structure lacks fails that
            structure (in `missing`), not the step.
        omit: residue tokens never sampled (3-letter codes, e.g. ["CYS"]). Replaces upstream's default
            ["UNK"], so keep "UNK" in the list when giving one.
        bias: {token: logit bias} applied at every designed position, e.g. {"ALA": -1.0}.
        chains: which chains of the backbone become `sequences` rows, with the same meaning as for
            `ProteinMPNN` and `LigandMPNN`. Upstream writes every chain of a design into one unseparated
            record, so the chains are told apart from the written structure. None (default) expects a
            single-chain backbone and keeps ids as `<structure>_<n>`; "A"/["A"] emits only that chain with
            the same ids; ["A","B"] emits `<structure>_<n>_<chain>`; "all" emits every chain with lazy
            `[_<?>]` ids. A multi-chain backbone under the default is a failure in `missing`, not a fused
            sequence. Only chains of standard amino acids count; nucleic-acid and ligand chains are context.
        checkpoint_path: an explicit weights file (upstream `checkpoint_path`), for weights other than the
            ones `FoundryMPNN.install()` fetched.
        is_legacy_weights: whether the checkpoint uses the original repositories' parameter ordering
            (default True, which is what the downloaded weights are and what upstream's README prescribes
            for them).

    Outputs:
        Streams:
            sequences: one row per designed chain, ids `<structure>_<n>` (n = 1..num_sequences*num_batches,
                       batch-major), with a `_<chain>` suffix when `chains` selects several.
            designs:   one row per design, before the chain split; the key chain rows are grouped by.
        Tables:
            sequences: id | design | structures.id | source_pdb | chain | sequence | sample | T | seed |
                       overall_confidence | ligand_confidence | seq_rec | ligand_seq_rec
                       `overall_confidence` and `ligand_confidence` are upstream's `confidence` and
                       `ligand_interface_confidence` (exp of the negative log-likelihood of the sampled
                       sequence) and `seq_rec` and `ligand_seq_rec` its `sequence_recovery` and
                       `ligand_interface_sequence_recovery`. A value a given release does not write is NaN;
                       the confidences are absent from the 0.2.0 release.
            designs:   id | structures.id | source_pdb | n_chains | overall_confidence | ligand_confidence |
                       seq_rec | ligand_seq_rec
            missing:   id | removed_by | kind | cause
    """

    TOOL_NAME = "FoundryMPNN"
    TOOL_VERSION = "1.1"
    # The foundry environment RFdiffusion3 already builds; the entry in each config's `environments:` names it.
    ENV_NAME = "foundry"

    MPNN_COMMAND = "mpnn"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        # Not delegated to RFdiffusion3: its install fetches only rfd3 weights and its "already installed" shortcut would skip these.
        env = cls._install_env(env_manager)
        ckpt_dir = folders.get("FoundryMPNN", "")
        if not ckpt_dir:
            raise ValueError(
                "FoundryMPNN: folders.repositories.FoundryMPNN is not set in the active config; "
                "it names the folder the MPNN weights are downloaded into")
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check(env, env_manager)
        run = cls._env_run(env, env_manager)
        have_weights = (f'ls "{ckpt_dir}"/proteinmpnn_*.pt >/dev/null 2>&1 '
                        f'&& ls "{ckpt_dir}"/ligandmpnn_*.pt >/dev/null 2>&1')
        skip = "" if force_reinstall else f"""# Check if already fully installed (env + both MPNN checkpoints)
if {env_check} && {have_weights}; then
    echo "FoundryMPNN already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        return f"""echo "=== Installing FoundryMPNN ({env}) ==="
{skip}{remove_block}
# The foundry env is shared with RFdiffusion3 and RF3: create it only when absent
if ! {env_check}; then
    {env_block}
else
    echo "{env} environment already exists, skipping creation."
fi

# Download the ProteinMPNN and LigandMPNN weights (skip if both are already there)
mkdir -p "{ckpt_dir}"
if ! ({have_weights}); then
    {run}foundry install proteinmpnn ligandmpnn --checkpoint-dir "{ckpt_dir}"
else
    echo "Checkpoint dir {ckpt_dir} already holds both MPNN checkpoints, skipping weight download."
fi

# Verify installation
if ({have_weights}) && {run}python -c "import mpnn.inference_engines.mpnn" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== FoundryMPNN installation complete ==="
else
    echo "ERROR: FoundryMPNN verification failed (checkpoints missing or mpnn not importable)"
    exit 1
fi
"""

    # Lazy path descriptors
    structures_json = Path(lambda self: self.configuration_path("input_structures.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    mpnn_out_folder = Path(lambda self: self.execution_path("mpnn"))
    checkpoint_dir = Path(lambda self: self.folders["FoundryMPNN"])

    sequences_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))
    designs_csv = Path(lambda self: self.stream_path("designs", "designs.csv"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))

    run_py = Path(lambda self: self.pipe_script_path("pipe_foundry_mpnn.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 model_type: str = "protein_mpnn",
                 num_sequences: int = 1,
                 num_batches: int = 1,
                 temperature: float = 0.1,
                 seed: Optional[int] = None,
                 structure_noise: Optional[float] = None,
                 designed_residues: Union[str, List[str], None] = None,
                 fixed_residues: Union[str, List[str], None] = None,
                 designed_chains: Union[str, List[str], None] = None,
                 fixed_chains: Union[str, List[str], None] = None,
                 omit: Optional[List[str]] = None,
                 bias: Optional[Dict[str, float]] = None,
                 chains: Union[str, List[str], None] = None,
                 checkpoint_path: Optional[str] = None,
                 is_legacy_weights: bool = True,
                 **kwargs):
        # Keep the original object: only a StandardizedOutput carries the upstream missing manifest.
        self.structures = structures
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        self.model_type = model_type
        self.num_sequences = num_sequences
        self.num_batches = num_batches
        self.temperature = temperature
        self.seed = seed
        self.structure_noise = structure_noise
        self.designed_residues = as_id_list(designed_residues)
        self.fixed_residues = as_id_list(fixed_residues)
        self.designed_chains = as_id_list(designed_chains)
        self.fixed_chains = as_id_list(fixed_chains)
        self.omit = list(omit) if omit is not None else None
        self.bias = dict(bias) if bias is not None else None
        self.chains = normalize_chains(chains)
        self.checkpoint_path = checkpoint_path
        self.is_legacy_weights = is_legacy_weights

        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")

        if self.model_type not in MODEL_TYPES:
            raise ValueError(f"model_type must be one of {MODEL_TYPES}, got {self.model_type!r}")

        for name in ("num_sequences", "num_batches"):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 1:
                raise ValueError(f"{name} must be a positive integer")

        if (not isinstance(self.temperature, (int, float)) or isinstance(self.temperature, bool)
                or self.temperature <= 0):
            raise ValueError("temperature must be a positive number")

        if self.seed is not None and (not isinstance(self.seed, int) or isinstance(self.seed, bool)):
            raise ValueError("seed must be an integer")

        if self.structure_noise is not None and (
                not isinstance(self.structure_noise, (int, float)) or isinstance(self.structure_noise, bool)
                or self.structure_noise < 0):
            raise ValueError("structure_noise must be a non-negative number of Angstrom")

        for name in ("designed_residues", "fixed_residues"):
            ids = getattr(self, name)
            if ids is not None and (not ids or any(not RESIDUE_ID.match(i) for i in ids)):
                raise ValueError(
                    f"{name} must be residue ids such as \"A35\" or \"B40A\" (chain, number, optional "
                    f"insertion code), got {ids!r}")
        for name in ("designed_chains", "fixed_chains"):
            ids = getattr(self, name)
            if ids is not None and (not ids or any(not CHAIN_ID.match(i) for i in ids)):
                raise ValueError(f"{name} must be chain ids such as \"A\", got {ids!r}")

        # Upstream's own scope rules, refused here so the failure is not found on the cluster.
        if self.designed_residues is not None and self.fixed_residues is not None:
            raise ValueError("pass designed_residues or fixed_residues, not both")
        if self.designed_chains is not None and self.fixed_chains is not None:
            raise ValueError("pass designed_chains or fixed_chains, not both")
        residue_scope = self.designed_residues is not None or self.fixed_residues is not None
        chain_scope = self.designed_chains is not None or self.fixed_chains is not None
        if residue_scope and chain_scope:
            raise ValueError("a residue selection and a chain selection cannot be mixed in one design")

        if self.omit is not None and (not isinstance(self.omit, list)
                                      or any(not isinstance(t, str) or not TOKEN.match(t) for t in self.omit)):
            raise ValueError("omit must be a list of residue tokens such as [\"UNK\", \"CYS\"]")
        if self.bias is not None:
            if not isinstance(self.bias, dict) or not self.bias:
                raise ValueError("bias must be a non-empty {token: value} dict")
            for token, value in self.bias.items():
                if not isinstance(token, str) or not TOKEN.match(token):
                    raise ValueError(f"bias key {token!r} is not a residue token such as \"ALA\"")
                if not isinstance(value, (int, float)) or isinstance(value, bool):
                    raise ValueError(f"bias[{token!r}] must be a number")

        validate_chains(self.chains)

        if not isinstance(self.is_legacy_weights, bool):
            raise ValueError("is_legacy_weights must be a bool")
        if self.checkpoint_path is not None:
            if not isinstance(self.checkpoint_path, str) or not self.checkpoint_path.strip():
                raise ValueError("checkpoint_path must be a non-empty path")
            _validate_freeform_string("checkpoint_path", self.checkpoint_path)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"MODEL: {self.model_type}")
        lines.append(f"NUM SEQUENCES: {self.num_sequences} x {self.num_batches} batch(es)")
        lines.append(f"TEMPERATURE: {self.temperature}")
        lines.append(f"CHAINS: {self.chains if self.chains else 'single (expected)'}")
        for name in ("designed_residues", "fixed_residues", "designed_chains", "fixed_chains"):
            if getattr(self, name) is not None:
                lines.append(f"{name.upper().replace('_', ' ')}: {','.join(getattr(self, name))}")
        if self.checkpoint_path:
            lines.append(f"CHECKPOINT: {self.checkpoint_path}")
        return lines

    def _config_dict(self) -> Dict[str, Any]:
        """Everything the pipe script needs besides paths; unset optional keys are left out."""
        config: Dict[str, Any] = {
            "model_type": self.model_type,
            "is_legacy_weights": self.is_legacy_weights,
            "num_sequences": self.num_sequences,
            "num_batches": self.num_batches,
            "temperature": self.temperature,
            "chains": self.chains,
        }
        optional = {
            "checkpoint_path": self.checkpoint_path,
            "seed": self.seed,
            "structure_noise": self.structure_noise,
            "designed_residues": self.designed_residues,
            "fixed_residues": self.fixed_residues,
            "designed_chains": self.designed_chains,
            "fixed_chains": self.fixed_chains,
            "omit": self.omit,
            "bias": self.bias,
        }
        config.update({key: value for key, value in optional.items() if value is not None})
        return config

    def _write_config_yaml(self):
        import yaml
        with open(self.config_yaml, "w") as handle:
            yaml.safe_dump(self._config_dict(), handle, sort_keys=False)

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self._write_config_yaml()

        script = "#!/bin/bash\n"
        script += "# FoundryMPNN script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += self._generate_run_section()
        script += self.generate_missing_propagation(
            self.structures, local_missing=self.local_missing_csv, missing_csv=self.missing_csv)
        script += self.generate_completion_check_footer()
        return script

    def _generate_run_section(self) -> str:
        # Architecture (b): the pipe stays on the host and prefixes only its `mpnn` call, once per structure.
        return f"""echo "Running FoundryMPNN ({self.model_type})"

# mpnn resolves its default weights through this search path
export FOUNDRY_CHECKPOINT_DIRS="{self.checkpoint_dir}"

python "{self.run_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-json "{self.structures_json}" \\
    --out-dir "{self.mpnn_out_folder}" \\
    --sequences-csv "{self.sequences_csv}" \\
    --designs-csv "{self.designs_csv}" \\
    --local-missing-csv "{self.local_missing_csv}" \\
    --container-prefix "{self.container_prefix()}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        # Designs per structure = num_sequences x num_batches, numbered batch-major.
        total = self.num_sequences * self.num_batches
        design_ids = generate_multiplied_ids_pattern(
            self.structures_stream.ids, f"<1..{total}>", input_stream_name="structures")
        sequence_ids = chain_row_ids(design_ids, self.chains)

        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[],
            map_table=self.sequences_csv,
            format="csv",
        )
        # The grouping key for the chain rows above: one row per design.
        designs = DataStream(
            name="designs",
            ids=design_ids,
            files=[],
            map_table=self.designs_csv,
            format="csv",
        )

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "design", "structures.id", "source_pdb", "chain", "sequence", "sample", "T",
                         "seed", "overall_confidence", "ligand_confidence", "seq_rec", "ligand_seq_rec"],
                description="FoundryMPNN sequence design results, one row per designed chain",
            ),
            "designs": TableInfo(
                name="designs",
                path=self.designs_csv,
                columns=["id", "structures.id", "source_pdb", "n_chains", "overall_confidence",
                         "ligand_confidence", "seq_rec", "ligand_seq_rec"],
                description="One row per design, before the chain split",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }

        return {
            "sequences": sequences,
            "designs": designs,
            "tables": tables,
            "output_folder": self.output_folder,
        }
