# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

import os
from typing import Dict, List, Any, Union, Optional

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import (generate_combinatorics_config, generate_multiplied_ids_pattern,
                                predict_output_ids_with_provenance)
    from .datastream_resolver import resolve_input_to_datastream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import (generate_combinatorics_config, generate_multiplied_ids_pattern,
                               predict_output_ids_with_provenance)
    from datastream_resolver import resolve_input_to_datastream


class RF3(BaseConfig):
    """
    RF3: RoseTTAFold3 all-atom co-folding of proteins, nucleic acids and ligands.

    The structure-prediction network of RosettaCommons Foundry (https://github.com/RosettaCommons/foundry),
    run through its `rf3 fold` command. It takes the same axes as `Boltz2` and `OpenFold3`, so the three are
    interchangeable in a pipeline, but it is not a drop-in replacement in what it predicts:

      * there is no affinity head, so no affinity table;
      * covalent bonds, templates, cyclic chains and ground-truth conformers exist upstream but are not
        wrapped here (a covalent modification needs the `bonds` list of an upstream JSON or a prepared CIF).

    It shares the `foundry` environment with `RFdiffusion3`, so one environment serves design and prediction.
    Its weights are its own: `RF3.install()` downloads them into the `RF3` repositories folder.

    Inputs:
        proteins: protein sequences.
        ssDNA / dsDNA / ssRNA / dsRNA: nucleic-acid sequences. A double-stranded axis becomes two chains, the
            second the reverse complement. The chain type is passed to upstream explicitly, because a string
            such as "ACGT" is also a valid protein sequence.
        ligands: compounds, by SMILES or CCD code. Upstream's own guidance is that a CCD code folds better
            than a SMILES string.
        msas: precomputed alignments for the protein axis, one per protein id (`.a3m`, `.a3m.gz` or `.fasta`).
            Upstream cannot compute an MSA itself, so without this the proteins fold single-sequence.
            With it, every protein chain of an input needs an alignment: an input with a chain that has none is not predicted and lands in `missing` as a `failure`.
            Upstream says it pairs alignments across chains on the fly and does not accept pre-paired ones.
        n_recycles: trunk recycles (upstream default 10). A single recycle is refused while early stopping is on, because upstream spends the first recycle on the early-stopping check and has none left to fold with.
        diffusion_batch_size: structures sampled per input from one forward pass of the trunk (upstream default 5).
            Always passed explicitly, so the declared ids of `top_only=False` are exact.
        num_steps: diffusion sampling steps. Upstream's config default is 50; its README states 200 and says 50
            performs similarly and more than twice as fast. Left to upstream unless given.
        seed: model seed. Left to upstream when None.
        early_stopping_plddt_threshold: an input whose mean all-atom pLDDT after the first recycle falls below
            this (0 to 1) is abandoned: no structure is written, and the id lands in `missing` as a `filter`.
            Upstream's default is 0.5 and so is this one, because the cost of a doomed prediction is the point
            of the option. 0 turns it off. Without an MSA the upstream examples trigger it.
        top_only: keep only the best sample of each input as `<id>` (default; best by upstream's ranking score,
            0.8 * ipTM + 0.2 * pTM, minus a penalty for clashing chains). False surfaces every sample as
            `<id>_1..N`, N being `diffusion_batch_size`.
        checkpoint_path: an explicit checkpoint file (upstream `ckpt_path`), for a model other than the one
            `RF3.install()` fetched.
        Any other keyword is forwarded to `rf3 fold` as a Hydra `key=value` override (for example
        `template_selection`, `cyclic_chains`, `residue_renaming_dict`). Only the parameters above are
        validated or reflected in the declared outputs.

    Outputs:
        Streams:
            structures: one `<id>.cif` per predicted complex. The B-factor column carries per-atom pLDDT on
                        a 0 to 1 scale (upstream's `annotate_b_factor_with_plddt`, on in its config).
            pae:        one `<id>.json` per structure, with the same ids as `structures`. Upstream's full
                        confidence file for that sample (`*_confidences.json`), copied unchanged. A JSON object
                        with keys `atom_chain_ids`, `atom_plddts` (N_atom floats, 0 to 1), `pae`
                        (N_token x N_token nested lists, in Angstrom, 0 to 32), `token_chain_ids`, `token_res_ids`.
                        Read the matrix as `numpy.array(json.load(f)["pae"])`. Tokens run over all chains,
                        one per standard residue and one per atom for ligands and modified residues.
                        `token_res_ids` is a running token index in upstream's source, not a residue number.
        Tables:
            confidence: id | plddt | ptm | iptm | pae | pde | has_clash | ranking_score | early_stopped | seed | sample
                        one row per id that entered, NaN where RF3 reported nothing. `plddt`, `pae` and `pde`
                        are upstream's `overall_*` means; `ptm`, `iptm`, `has_clash` and `ranking_score` are
                        its own. An input stopped early has `early_stopped` True and NaN everywhere else.
                        Upstream reports an `iptm` of 0 for a single chain.
            missing:    id | removed_by | kind | cause
    """

    TOOL_NAME = "RF3"
    TOOL_VERSION = "1.4"
    # `rf3 fold` is a Hydra entry point: an untyped kwarg is rendered as one more `key=value` override.
    FORWARD_UNKNOWN_KWARGS = "hydra"
    # The foundry environment RFdiffusion3 already builds; the entry in each config's `environments:` names it.
    ENV_NAME = "foundry"

    PREDICT_COMMAND = "rf3 fold"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        # Not delegated to RFdiffusion3: its install fetches only rfd3 weights and its "already installed" shortcut would skip these.
        env = cls._install_env(env_manager)
        ckpt_dir = folders.get("RF3", "")
        if not ckpt_dir:
            raise ValueError(
                "RF3: folders.repositories.RF3 is not set in the active config; "
                "it names the folder the RF3 weights are downloaded into")
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check(env, env_manager)
        run = cls._env_run(env, env_manager)
        have_weights = f'ls "{ckpt_dir}"/*.ckpt >/dev/null 2>&1'
        skip = "" if force_reinstall else f"""# Check if already fully installed (env + RF3 checkpoint)
if {env_check} && {have_weights}; then
    echo "RF3 already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        return f"""echo "=== Installing RF3 ({env}) ==="
{skip}{remove_block}
# The foundry env is shared with RFdiffusion3: create it only when absent
if ! {env_check}; then
    {env_block}
else
    echo "{env} environment already exists, skipping creation."
fi

# Download the RF3 weights (skip if a checkpoint is already there)
mkdir -p "{ckpt_dir}"
if ! {have_weights}; then
    {run}foundry install rf3 --checkpoint-dir "{ckpt_dir}"
else
    echo "Checkpoint dir {ckpt_dir} already holds an RF3 checkpoint, skipping weight download."
fi

# Verify installation
if {have_weights} && {run}python -c "import rf3.inference_engines.rf3" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== RF3 installation complete ==="
else
    echo "ERROR: RF3 verification failed (checkpoint missing or rf3 not importable)"
    exit 1
fi
"""

    # Configuration-time artifacts
    combinatorics_config_file = Path(lambda self: self.configuration_path("combinatorics_config.json"))
    msas_json = Path(lambda self: self.configuration_path("msas.json"))
    # Written at EXECUTION time: one JSON list of every example, so `rf3 fold` loads the model once.
    queries_json = Path(lambda self: self.configuration_path("queries.json"))
    skipped_csv = Path(lambda self: self.configuration_path("skipped_queries.csv"))

    prediction_folder = Path(lambda self: self.execution_folder)
    checkpoint_dir = Path(lambda self: self.folders["RF3"])

    structures_map_csv = Path(lambda self: self.stream_map_path("structures"))
    pae_folder = Path(lambda self: self.stream_folder("pae"))
    pae_map_csv = Path(lambda self: self.stream_map_path("pae"))
    confidence_csv = Path(lambda self: self.table_path("confidence"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))

    config_py = Path(lambda self: self.pipe_script_path("pipe_rf3_config.py"))
    postprocess_py = Path(lambda self: self.pipe_script_path("pipe_rf3_postprocess.py"))

    def __init__(self,
                 proteins: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ssDNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 dsDNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ssRNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 dsRNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ligands: Optional[Union[DataStream, StandardizedOutput]] = None,
                 msas: Optional[Union[DataStream, StandardizedOutput]] = None,
                 n_recycles: Optional[int] = None,
                 diffusion_batch_size: int = 5,
                 num_steps: Optional[int] = None,
                 seed: Optional[int] = None,
                 early_stopping_plddt_threshold: float = 0.5,
                 top_only: bool = True,
                 checkpoint_path: Optional[str] = None,
                 **kwargs):
        # Raw handles kept for missing-propagation; the resolved streams for the config.
        self.proteins = proteins
        self.ssDNA = ssDNA
        self.dsDNA = dsDNA
        self.ssRNA = ssRNA
        self.dsRNA = dsRNA
        self.ligands = ligands
        self.msas = msas

        self.proteins_stream = self._axis_stream(proteins, "sequences")
        self.ssDNA_stream = self._axis_stream(ssDNA, "sequences")
        self.dsDNA_stream = self._axis_stream(dsDNA, "sequences")
        self.ssRNA_stream = self._axis_stream(ssRNA, "sequences")
        self.dsRNA_stream = self._axis_stream(dsRNA, "sequences")
        self.ligands_stream = self._axis_stream(ligands, "compounds")
        self.msas_stream = self._axis_stream(msas, "msas")

        self.n_recycles = n_recycles
        self.diffusion_batch_size = diffusion_batch_size
        self.num_steps = num_steps
        self.seed = seed
        self.early_stopping_plddt_threshold = early_stopping_plddt_threshold
        self.top_only = top_only
        self.checkpoint_path = checkpoint_path

        super().__init__(**kwargs)

    @staticmethod
    def _axis_stream(value, fallback_stream: str):
        if value is None:
            return None
        return resolve_input_to_datastream(value, fallback_stream=fallback_stream)

    def validate_params(self):
        polymers = [self.proteins, self.ssDNA, self.dsDNA, self.ssRNA, self.dsRNA]
        if not any(axis is not None for axis in polymers) and self.ligands is None:
            raise ValueError(
                "RF3 needs at least one input axis: proteins, ssDNA, dsDNA, ssRNA, dsRNA or ligands")

        for name in ("n_recycles", "diffusion_batch_size", "num_steps"):
            value = getattr(self, name)
            if value is None and name != "diffusion_batch_size":
                continue
            if not isinstance(value, int) or isinstance(value, bool) or value < 1:
                raise ValueError(f"{name} must be a positive integer")

        if self.seed is not None and (not isinstance(self.seed, int) or isinstance(self.seed, bool)):
            raise ValueError("seed must be an integer")

        threshold = self.early_stopping_plddt_threshold
        if (not isinstance(threshold, (int, float)) or isinstance(threshold, bool)
                or not 0 <= threshold <= 1):
            raise ValueError(
                "early_stopping_plddt_threshold must be a number from 0 to 1 (0 turns early stopping off)")

        if self.n_recycles == 1 and threshold > 0:
            raise ValueError(
                "n_recycles=1 needs early_stopping_plddt_threshold=0: upstream takes the first recycle for the "
                "early-stopping check and has none left to fold with, so every input fails")

        if not isinstance(self.top_only, bool):
            raise ValueError("top_only must be a bool")

        if self.checkpoint_path is not None:
            if not isinstance(self.checkpoint_path, str) or not self.checkpoint_path.strip():
                raise ValueError("checkpoint_path must be a non-empty path")
            _validate_freeform_string("checkpoint_path", self.checkpoint_path)

        if self.msas is not None:
            if self.proteins is None:
                raise ValueError("msas supplies alignments for proteins, so a proteins axis is required")
            fmt = getattr(getattr(self, "msas_stream", None), "format", None)
            if fmt is not None and str(fmt).lower() == "csv":
                raise ValueError("RF3 reads MSAs as a3m or fasta; this msas stream is csv. Use "
                                 "MMseqs2(output_format=\"a3m\") or MSA(..., convert=\"a3m\")")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        axes = [name for name in ("proteins", "ssDNA", "dsDNA", "ssRNA", "dsRNA", "ligands")
                if getattr(self, name) is not None]
        lines.append(f"AXES: {', '.join(axes)}")
        lines.append(f"MSA: {'precomputed' if self.msas is not None else 'none (single sequence)'}")
        lines.append(f"DIFFUSION SAMPLES: {self.diffusion_batch_size}")
        if self.n_recycles is not None:
            lines.append(f"RECYCLES: {self.n_recycles}")
        if self.num_steps is not None:
            lines.append(f"DIFFUSION STEPS: {self.num_steps}")
        if self.seed is not None:
            lines.append(f"SEED: {self.seed}")
        lines.append(f"EARLY STOPPING pLDDT: {self.early_stopping_plddt_threshold}")
        if self.checkpoint_path:
            lines.append(f"CHECKPOINT: {self.checkpoint_path}")
        return lines

    def _axis_kwargs(self) -> Dict:
        kwargs = {}
        if self.proteins is not None:
            kwargs["proteins"] = (self.proteins, "sequences", "protein")
        if self.ssDNA is not None:
            kwargs["ssDNA"] = (self.ssDNA, "sequences", "ssdna")
        if self.dsDNA is not None:
            kwargs["dsDNA"] = (self.dsDNA, "sequences", "dsdna")
        if self.ssRNA is not None:
            kwargs["ssRNA"] = (self.ssRNA, "sequences", "ssrna")
        if self.dsRNA is not None:
            kwargs["dsRNA"] = (self.dsRNA, "sequences", "dsrna")
        if self.ligands is not None:
            kwargs["ligands"] = (self.ligands, "compounds", "ligand")
        return kwargs

    def _missing_input_sources(self):
        return tuple(axis for axis in
                     (self.proteins, self.ssDNA, self.dsDNA, self.ssRNA, self.dsRNA, self.ligands)
                     if axis is not None)

    def generate_script(self, script_path: str) -> str:
        generate_combinatorics_config(self.combinatorics_config_file, **self._axis_kwargs())

        script = "#!/bin/bash\n"
        script += "# RF3 script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += self._generate_config_section()
        script += self._generate_predict_section()
        script += self._generate_postprocess_section()
        script += self.generate_missing_propagation(
            *self._missing_input_sources(),
            local_missing=self.local_missing_csv,
            missing_csv=self.missing_csv,
        )
        script += self.generate_completion_check_footer()
        return script

    def _generate_config_section(self) -> str:
        """Build the inputs JSON at execution time, when the sequences actually exist.

        The wrapper cannot write it: the sequences live in upstream streams whose content is
        produced by the step before this one.
        """
        msa_flag = ""
        if self.msas_stream:
            self.msas_stream.save_json(self.msas_json)
            msa_flag = f' \\\n    --msas-json "{self.msas_json}"'
        return f"""echo "Building RF3 inputs JSON"
python "{self.config_py}" \\
    --combinatorics-config "{self.combinatorics_config_file}" \\
    --queries-json "{self.queries_json}" \\
    --skipped-csv "{self.skipped_csv}"{msa_flag}

"""

    def _predict_overrides(self) -> List[str]:
        """The Hydra `key=value` overrides for `rf3 fold`, in the order they are emitted."""
        overrides = [f'inputs="{self.queries_json}"',
                     f'out_dir="{self.prediction_folder}"',
                     f"diffusion_batch_size={self.diffusion_batch_size}",
                     f"early_stopping_plddt_threshold={self.early_stopping_plddt_threshold}"]
        if self.n_recycles is not None:
            overrides.append(f"n_recycles={self.n_recycles}")
        if self.num_steps is not None:
            overrides.append(f"num_steps={self.num_steps}")
        if self.seed is not None:
            overrides.append(f"seed={self.seed}")
        if self.checkpoint_path:
            overrides.append(f'ckpt_path="{self.checkpoint_path}"')
        overrides += self.extra_args_bash_tokens()
        return overrides

    def _generate_predict_section(self) -> str:
        joined = " \\\n    ".join(self._predict_overrides())
        return self.extra_args_echo() + f"""echo "Running RF3 prediction"

# rf3 resolves its weights through this search path
export FOUNDRY_CHECKPOINT_DIRS="{self.checkpoint_dir}"
if ! ls "${{FOUNDRY_CHECKPOINT_DIRS}}"/*.ckpt >/dev/null 2>&1; then
    echo "ERROR: no RF3 checkpoint under ${{FOUNDRY_CHECKPOINT_DIRS}}; run RF3.install() first"
    exit 1
fi

{self.container_prefix()}{self.PREDICT_COMMAND} \\
    {joined}

"""

    def _generate_postprocess_section(self) -> str:
        return f"""echo "Collecting RF3 predictions"
python "{self.postprocess_py}" \\
    --prediction-folder "{self.prediction_folder}" \\
    --queries-json "{self.queries_json}" \\
    --skipped-csv "{self.skipped_csv}" \\
    --structures-dir "{self.stream_folder('structures')}" \\
    --structures-map-csv "{self.structures_map_csv}" \\
    --confidence-csv "{self.confidence_csv}" \\
    --pae-dir "{self.pae_folder}" \\
    --pae-map-csv "{self.pae_map_csv}" \\
    --local-missing-csv "{self.local_missing_csv}" \\
    --top-only "{str(self.top_only)}" \\
    --diffusion-batch-size {self.diffusion_batch_size} \\
    --early-stopping-plddt-threshold {self.early_stopping_plddt_threshold}

"""

    def get_output_files(self) -> Dict[str, Any]:
        predicted_ids, _provenance = predict_output_ids_with_provenance(**self._axis_kwargs())

        # With top_only off every sample is surfaced as <id>_1..N, and the declared ids must say so.
        if self.top_only:
            structure_ids = predicted_ids
        else:
            structure_ids = generate_multiplied_ids_pattern(
                predicted_ids, f"<1..{self.diffusion_batch_size}>")

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=[self.stream_path("structures", "<id>.cif")],
            map_table=self.structures_map_csv,
            format="cif",
        )

        pae = DataStream(
            name="pae",
            ids=structure_ids,
            files=[self.stream_path("pae", "<id>.json")],
            map_table=self.pae_map_csv,
            format="json",
        )

        tables = {
            "confidence": TableInfo(
                name="confidence",
                path=self.confidence_csv,
                columns=["id", "plddt", "ptm", "iptm", "pae", "pde", "has_clash",
                         "ranking_score", "early_stopped", "seed", "sample"],
                description="Per-complex confidence scores from RF3",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }

        return {
            "structures": structures,
            "pae": pae,
            "tables": tables,
            "output_folder": self.output_folder,
        }
