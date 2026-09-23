# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

import os
from typing import Dict, List, Any, Union, Optional

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import (generate_combinatorics_config, get_mode,
                                predict_output_ids_with_provenance)
    from .datastream_resolver import resolve_input_to_datastream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import (generate_combinatorics_config, get_mode,
                               predict_output_ids_with_provenance)
    from datastream_resolver import resolve_input_to_datastream


OUTPUT_FORMATS = ("cif", "pdb")


class OpenFold3(BaseConfig):
    """
    OpenFold3: open-source co-folding of proteins, nucleic acids and ligands.

    The OpenFold Consortium / AlQuraishi lab reimplementation of AlphaFold3, under Apache 2.0 —
    weights included, no access request. It takes the same axes as `Boltz2` so the two are
    interchangeable in a pipeline, but it is **not** a drop-in replacement in what it predicts:

      * there is no affinity head, so no `affinity` parameter and no affinity table;
      * ligands are non-covalent, so there is no `covalent_linkage` equivalent.

    For a covalent probe or a binding-affinity readout, `Boltz2` remains the tool. OpenFold3 is
    for the structure itself, and for cross-checking a prediction against a second model.

    Inputs:
        proteins: protein sequences.
        ssDNA / dsDNA / ssRNA / dsRNA: nucleic-acid sequences. A double-stranded axis becomes
            two chains, the second the reverse complement.
        ligands: compounds, by SMILES or CCD code.
        msas: precomputed alignments (`.a3m`, `.sto` or `.npz`). Supplying them turns the MSA
            server off, since a run cannot take both.
        use_msa_server: generate alignments through the ColabFold server (default: on, unless
            `msas` is given). True together with `msas` is refused.
        num_diffusion_samples: structures sampled per query (upstream default 5).
        num_model_seeds: random seeds per query (upstream default 1).
        seeds: the explicit seed values, when the run has to be reproducible by seed.
        output_format: "cif" (default) or "pdb". A `.pdb` carries per-atom pLDDT in its B-factor.
        top_only: keep only the best sample as `<id>` (default). False surfaces every sample as
            `<id>_1..N`.
        inference_ckpt_name: a checkpoint from the upstream list (default "openfold3_p2_v1").
        inference_ckpt_path: an explicit checkpoint file — the way to run a privately fine-tuned
            model, such as a federated checkpoint, against the same pipeline.
        template: a CIF template applied to every protein chain.
        template_chain_ids: which chains the template applies to.
        low_mem: add the upstream `low_mem` preset, for a GPU that cannot hold the default.
        devices: GPUs to distribute over (upstream default 1).
        runner_yaml: a user-supplied runner YAML, merged under what the wrapper writes. The
            escape hatch for upstream settings this wrapper does not name.

    Outputs:
        Streams:
            structures: one `<id>.cif` (or `.pdb`) per predicted complex.
        Tables:
            confidence: id | plddt | ptm | iptm | gpde | has_clash | disorder |
                        ranking_score | seed | sample
            missing:    id | removed_by | kind | cause
    """

    TOOL_NAME = "OpenFold3"
    TOOL_VERSION = "1.4"
    ENV_NAME = "openfold3"

    # Upstream spells its flags with underscores. The auto-generated option table in the docs
    # renders them with dashes, but both worked examples -- the README quick start and the
    # precomputed-MSA how-to -- use underscores, so that is what is emitted.
    PREDICT_COMMAND = "run_openfold predict"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        run = cls._env_run(env, env_manager)
        # The weights are NOT fetched on first prediction: run_openfold validates the checkpoint
        # up front and refuses ("cowardly refusing to perform inference"), so a pipeline that
        # installed cleanly still dies at the first query. `setup_openfold` is the downloader,
        # and it runs on every install -- including the already-installed path, or an env made
        # before this step existed never acquires its weights.
        weights = f"""
echo "Ensuring OpenFold3 parameters"
if ! {run}setup_openfold --non-interactive; then
    echo "ERROR: setup_openfold could not download the model parameters"
    exit 1
fi
"""
        # Without this every install recreated the env and re-resolved openfold3.
        skip = "" if force_reinstall else f"""if {cls._env_exists_check(env, env_manager)} && {run}python -c "import openfold3" >/dev/null 2>&1; then
    echo "OpenFold3 already installed, skipping the env. Use force_reinstall=True to rebuild it."
{weights}
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        return f"""echo "=== Installing OpenFold3 ==="
{skip}{remove_block}
{env_block}
{weights}
if {run}python -c "import openfold3" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== OpenFold3 installation complete ==="
else
    echo "ERROR: OpenFold3 verification failed"
    exit 1
fi
"""

    # Configuration-time artifacts
    combinatorics_config_file = Path(lambda self: self.configuration_path("combinatorics_config.json"))
    runner_yaml_file = Path(lambda self: self.configuration_path("runner.yaml"))
    # Written at EXECUTION time: one JSON holding every query, since run_openfold folds a whole
    # batch in one process rather than one file at a time.
    queries_json = Path(lambda self: self.configuration_path("queries.json"))

    prediction_folder = Path(lambda self: self.execution_folder)

    structures_map_csv = Path(lambda self: self.stream_map_path("structures"))
    confidence_csv = Path(lambda self: self.table_path("confidence"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))

    config_py = Path(lambda self: self.pipe_script_path("pipe_openfold3_config.py"))
    postprocess_py = Path(lambda self: self.pipe_script_path("pipe_openfold3_postprocess.py"))

    def __init__(self,
                 proteins: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ssDNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 dsDNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ssRNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 dsRNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ligands: Optional[Union[DataStream, StandardizedOutput]] = None,
                 msas: Optional[Union[DataStream, StandardizedOutput]] = None,
                 use_msa_server: Optional[bool] = None,
                 num_diffusion_samples: Optional[int] = None,
                 num_model_seeds: Optional[int] = None,
                 seeds: Optional[List[int]] = None,
                 output_format: str = "cif",
                 top_only: bool = True,
                 inference_ckpt_name: Optional[str] = None,
                 inference_ckpt_path: Optional[str] = None,
                 template: Optional[str] = None,
                 template_chain_ids: Optional[List[str]] = None,
                 low_mem: bool = False,
                 devices: Optional[int] = None,
                 runner_yaml: Optional[str] = None,
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

        # Unset means: the server unless precomputed MSAs were given, as the docstring promises.
        self.use_msa_server = (msas is None) if use_msa_server is None else use_msa_server
        self.num_diffusion_samples = num_diffusion_samples
        self.num_model_seeds = num_model_seeds
        self.seeds = seeds
        self.output_format = output_format
        self.top_only = top_only
        self.inference_ckpt_name = inference_ckpt_name
        self.inference_ckpt_path = inference_ckpt_path
        self.template = template
        self.template_chain_ids = template_chain_ids
        self.low_mem = low_mem
        self.devices = devices
        self.runner_yaml = runner_yaml

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
                "OpenFold3 needs at least one input axis: proteins, ssDNA, dsDNA, ssRNA, "
                "dsRNA or ligands")

        if self.output_format not in OUTPUT_FORMATS:
            raise ValueError(
                f"output_format must be one of {OUTPUT_FORMATS}, got {self.output_format!r}")

        for name in ("num_diffusion_samples", "num_model_seeds", "devices"):
            value = getattr(self, name)
            if value is None:
                continue
            if not isinstance(value, int) or isinstance(value, bool) or value < 1:
                raise ValueError(f"{name} must be a positive integer")

        if self.seeds is not None:
            if not isinstance(self.seeds, (list, tuple)) or not self.seeds:
                raise ValueError("seeds must be a non-empty list of integers")
            if any(not isinstance(s, int) or isinstance(s, bool) for s in self.seeds):
                raise ValueError("seeds must contain integers only")
            # Upstream takes both, and disagreeing values are not reconcilable: the count says
            # how many seeds to draw, the list says which. Refuse rather than pick one.
            if self.num_model_seeds is not None and self.num_model_seeds != len(self.seeds):
                raise ValueError(
                    f"num_model_seeds ({self.num_model_seeds}) contradicts the {len(self.seeds)} "
                    f"value(s) in seeds; pass one or the other")

        if self.inference_ckpt_name and self.inference_ckpt_path:
            raise ValueError(
                "pass inference_ckpt_name or inference_ckpt_path, not both")

        if self.msas is not None and self.use_msa_server:
            raise ValueError(
                "msas supplies precomputed alignments, so use_msa_server must be False; "
                "OpenFold3 cannot take both a precomputed MSA and a server-generated one")

        if self.template_chain_ids and not self.template:
            raise ValueError("template_chain_ids has no effect without template")

        if not self.top_only and self.runner_yaml and self.num_diffusion_samples is None:
            # The sample count then lives in a YAML read only at run time, and the declared ids would guess it.
            raise ValueError("top_only=False with runner_yaml needs num_diffusion_samples, so the "
                             "per-sample ids can be declared")

        if self.msas is not None:
            fmt = getattr(getattr(self, "msas_stream", None), "format", None)
            if fmt is not None and str(fmt).lower() == "csv":
                raise ValueError("OpenFold3 reads MSAs as a3m; this msas stream is csv. Use "
                                 "MMseqs2(output_format=\"a3m\") or MSA(..., convert=\"a3m\")")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        axes = [name for name in ("proteins", "ssDNA", "dsDNA", "ssRNA", "dsRNA", "ligands")
                if getattr(self, name) is not None]
        lines.append(f"AXES: {', '.join(axes)}")
        lines.append(f"OUTPUT FORMAT: {self.output_format}")
        lines.append(f"MSA: {'ColabFold server' if self.use_msa_server else 'precomputed'}")
        if self.num_diffusion_samples is not None:
            lines.append(f"DIFFUSION SAMPLES: {self.num_diffusion_samples}")
        if self.seeds is not None:
            lines.append(f"SEEDS: {', '.join(str(s) for s in self.seeds)}")
        elif self.num_model_seeds is not None:
            lines.append(f"MODEL SEEDS: {self.num_model_seeds}")
        if self.inference_ckpt_path:
            lines.append(f"CHECKPOINT: {self.inference_ckpt_path}")
        elif self.inference_ckpt_name:
            lines.append(f"CHECKPOINT: {self.inference_ckpt_name}")
        if self.low_mem:
            lines.append("MEMORY PRESET: low_mem")
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

    def _write_runner_yaml(self):
        """Settings upstream reads from a runner YAML rather than the command line.

        Kept as a file rather than more flags because it is also the run's readable record of
        what the model was asked to do, and because `runner_yaml` lets a user add keys this
        wrapper does not name without the wrapper having to grow a parameter for each.
        """
        import yaml

        runner: Dict[str, Any] = {}
        presets = ["predict"] + (["low_mem"] if self.low_mem else [])
        runner["model_update"] = {"presets": presets}
        runner["output_writer_settings"] = {"structure_format": self.output_format}
        if self.seeds is not None:
            runner["experiment_settings"] = {"seeds": list(self.seeds)}
        if self.devices is not None:
            runner["pl_trainer_args"] = {"devices": self.devices}

        if self.runner_yaml:
            with open(self.runner_yaml) as handle:
                supplied = yaml.safe_load(handle) or {}
            # The user's file wins: it is the escape hatch, and a wrapper default silently
            # overriding what someone wrote explicitly is the opposite of an escape hatch.
            for key, value in supplied.items():
                if isinstance(value, dict) and isinstance(runner.get(key), dict):
                    runner[key].update(value)
                else:
                    runner[key] = value

        with open(self.runner_yaml_file, "w") as handle:
            yaml.safe_dump(runner, handle, sort_keys=False)

    def generate_script(self, script_path: str) -> str:
        generate_combinatorics_config(self.combinatorics_config_file, **self._axis_kwargs())
        self._write_runner_yaml()

        script = "#!/bin/bash\n"
        script += "# OpenFold3 script\n"
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
        """Build the query JSON at execution time, when the sequences actually exist.

        The wrapper cannot write it: the sequences live in upstream streams whose content is
        produced by the step before this one.
        """
        msa_flag = f' \\\n    --msas-json "{self.configuration_path("msas.json")}"' if self.msas_stream else ""
        if self.msas_stream:
            self.msas_stream.save_json(self.configuration_path("msas.json"))
        template_flag = ""
        if self.template:
            template_flag = f' \\\n    --template "{self.template}"'
            if self.template_chain_ids:
                template_flag += f' \\\n    --template-chains "{",".join(self.template_chain_ids)}"'
        return f"""echo "Building OpenFold3 query JSON"
python "{self.config_py}" \\
    --combinatorics-config "{self.combinatorics_config_file}" \\
    --queries-json "{self.queries_json}"{msa_flag}{template_flag}

"""

    def _generate_predict_section(self) -> str:
        options = [f'--query_json="{self.queries_json}"',
                   f'--output_dir="{self.prediction_folder}"',
                   f'--runner_yaml="{self.runner_yaml_file}"',
                   f'--use_msa_server={"True" if self.use_msa_server else "False"}']
        if self.num_diffusion_samples is not None:
            options.append(f"--num_diffusion_samples={self.num_diffusion_samples}")
        if self.num_model_seeds is not None:
            options.append(f"--num_model_seeds={self.num_model_seeds}")
        if self.inference_ckpt_path:
            options.append(f'--inference_ckpt_path="{self.inference_ckpt_path}"')
        elif self.inference_ckpt_name:
            options.append(f'--inference_ckpt_name="{self.inference_ckpt_name}"')
        forwarded = self.extra_args_bash()
        if forwarded:
            options.append(forwarded)
        cp = self.container_prefix()
        joined = " ".join(options)
        return self.extra_args_echo() + f"""echo "Running OpenFold3 prediction"
{cp}{self.PREDICT_COMMAND} {joined}

"""

    def _generate_postprocess_section(self) -> str:
        return f"""echo "Collecting OpenFold3 predictions"
python "{self.postprocess_py}" \\
    --prediction-folder "{self.prediction_folder}" \\
    --combinatorics-config "{self.combinatorics_config_file}" \\
    --queries-json "{self.queries_json}" \\
    --structures-dir "{self.stream_folder('structures')}" \\
    --structures-map-csv "{self.structures_map_csv}" \\
    --confidence-csv "{self.confidence_csv}" \\
    --local-missing-csv "{self.local_missing_csv}" \\
    --output-format "{self.output_format}" \\
    --top-only "{str(self.top_only)}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        predicted_ids, provenance = predict_output_ids_with_provenance(**self._axis_kwargs())

        extension = ".pdb" if self.output_format == "pdb" else ".cif"
        # With top_only off every sample is surfaced as <id>_1..N, so the declared ids have to
        # say so or the completion check demands files that were never meant to exist.
        if self.top_only:
            structure_ids = predicted_ids
        else:
            # The postprocess numbers every sample of every seed, so the count is their product.
            seeds = len(self.seeds) if self.seeds else (self.num_model_seeds or 1)
            samples = (self.num_diffusion_samples or 5) * seeds
            structure_ids = [f"{i}_<1..{samples}>" for i in predicted_ids]

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=[self.stream_path("structures", f"<id>{extension}")],
            map_table=self.structures_map_csv,
            format=self.output_format,
        )

        tables = {
            "confidence": TableInfo(
                name="confidence",
                path=self.confidence_csv,
                columns=["id", "plddt", "ptm", "iptm", "gpde", "has_clash",
                         "disorder", "ranking_score", "seed", "sample"],
                description="Aggregated per-complex confidence scores from OpenFold3",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }

        return {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder,
        }
