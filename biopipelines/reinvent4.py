# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""REINVENT4: generative design of small molecules by sampling a pretrained prior.

Wraps the `sampling` run mode of REINVENT 4 for four generators: de novo design (`reinvent`), scaffold decoration (`libinvent`), linker design (`linkinvent`) and molecule optimization (`mol2mol`). The sampled molecules become a first-class `compounds` stream.

Not implemented: transfer learning, staged learning (reinforcement learning with scoring components), the `scoring` run mode, the PepInvent generator, the transformer variants of LibInvent and LinkInvent, and the optional sampling filters (`filter.smarts`, `target_smiles_path`, `isomeric_smiles`, `randomize_smiles`, `tb_logdir`).

Reference:
    Loeffler et al. (2024) Reinvent 4: Modern AI-driven generative molecule design. J. Cheminform. 16, 20.
    https://github.com/MolecularAI/REINVENT4
"""

import json
import os
import re
from typing import Dict, List, Any, Union, Optional

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import generate_multiplied_ids_pattern
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_multiplied_ids_pattern


GENERATORS = ("reinvent", "libinvent", "linkinvent", "mol2mol")
CONDITIONAL_GENERATORS = ("libinvent", "linkinvent", "mol2mol")
SAMPLE_STRATEGIES = ("multinomial", "beamsearch")

# Zenodo record of the public priors: version 1 (2025-06-11), the one upstream's docs link.
PRIORS_RECORD = 15641297
PRIORS_URL = f"https://zenodo.org/api/records/{PRIORS_RECORD}/files"
DEFAULT_PROCESSOR = "cu126"

# md5 of each file as the Zenodo record reports it.
PRIOR_MD5 = {
    "reinvent.prior": "f268eb072f4fca69ca9434768d3cd461",
    "libinvent.prior": "fea30c293909216fe183bb0d1b7ca432",
    "linkinvent.prior": "65af13963aa2bf80fc941c0400288c8e",
    "mol2mol_similarity.prior": "451893c2d1818d9164bc73ea5a5bd7ef",
    "mol2mol_medium_similarity.prior": "78b75e478bfd7dcccee8ece010bc0b07",
    "mol2mol_high_similarity.prior": "35d80728589cfc69e2c5cc02a7fcab71",
    "mol2mol_scaffold.prior": "601d09cba7f57b1c717be3afbc550b8f",
    "mol2mol_scaffold_generic.prior": "d256fc90224e60bd84238697ada4a660",
    "mol2mol_mmp.prior": "d9051a0309e6a58a0e5262b1ccb39bdf",
    "pubchem_ecfp4_with_count_with_rank_reinvent4_dict_voc.prior": "97a16081c3b5897b9e900d2f448e0bf6",
}

GENERATOR_PRIOR = {
    "reinvent": "reinvent.prior",
    "libinvent": "libinvent.prior",
    "linkinvent": "linkinvent.prior",
}

# The Mol2Mol priors differ in the kind of similarity their training pairs have.
SIMILARITY_PRIOR = {
    "similarity": "mol2mol_similarity.prior",
    "medium_similarity": "mol2mol_medium_similarity.prior",
    "high_similarity": "mol2mol_high_similarity.prior",
    "scaffold": "mol2mol_scaffold.prior",
    "scaffold_generic": "mol2mol_scaffold_generic.prior",
    "mmp": "mol2mol_mmp.prior",
    "pubchem": "pubchem_ecfp4_with_count_with_rank_reinvent4_dict_voc.prior",
}
DEFAULT_SIMILARITY = "medium_similarity"

# Header of the CSV `reinvent` writes per generator, mapped to the column this tool emits.
UPSTREAM_COLUMNS = {
    "reinvent": {"SMILES": "smiles", "NLL": "nll"},
    "libinvent": {"SMILES": "smiles", "NLL": "nll", "Scaffold": "scaffold", "R-groups": "r_groups"},
    "linkinvent": {"SMILES": "smiles", "NLL": "nll", "Warheads": "warheads", "Linker": "linker"},
    "mol2mol": {"SMILES": "smiles", "NLL": "nll", "Input_SMILES": "input_smiles", "Tanimoto": "tanimoto"},
}

# Canonical compounds schema; generators with input compounds add a provenance column.
COMPOUND_COLUMNS = ["id", "format", "code", "smiles", "ccd"]
PROVENANCE_COLUMN = "compounds.id"
GENERATED_CODE = "LIG"
INPUT_COLUMNS = ["id", "smiles", "num_compounds"]

_PROCESSOR_RE = re.compile(r"^[A-Za-z0-9._-]+$")


def prior_file(generator: str, similarity: Optional[str] = None) -> str:
    """File name of the prior a generator (and, for mol2mol, a similarity type) samples from."""
    if generator == "mol2mol":
        return SIMILARITY_PRIOR[similarity or DEFAULT_SIMILARITY]
    return GENERATOR_PRIOR[generator]


# Resumable, checksum-verified download of one prior from the pinned Zenodo record.
_FETCH_PRIOR = """fetch_prior() {
    local name="$1" md5="$2" target="@PRIOR_DIR@/$1"
    if [ -s "$target" ] && [ "$(md5sum "$target" | cut -d' ' -f1)" = "$md5" ]; then
        echo "Prior $name already present"
        return 0
    fi
    echo "Downloading prior $name"
    curl -fL --retry 3 -C - -o "$target.part" "@PRIORS_URL@/$name/content" || return 1
    if [ "$(md5sum "$target.part" | cut -d' ' -f1)" != "$md5" ]; then
        echo "ERROR: checksum mismatch for $name"
        rm -f "$target.part"
        return 1
    fi
    mv "$target.part" "$target"
}"""


def _is_int(value) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


class REINVENT4(BaseConfig):
    """
    REINVENT4: sample molecules from a pretrained REINVENT 4 prior.

    Four generators, picked with `generator`:

        reinvent    de novo design; no input. Builds molecules atom by atom with an RNN.
        libinvent   scaffold decoration; `compounds` holds scaffolds, each with its attachment points marked `[*:0]`, `[*:1]`, ... The R-groups at those points are generated.
        linkinvent  linker design; `compounds` holds warhead pairs, each warhead with one attachment point `*`, the two joined by `|` (as `*c1ccc(F)cc1|*c1ccncc1`). The linker between them is generated.
        mol2mol     molecule optimization; `compounds` holds the molecules to optimize. Analogues within the similarity the prior was trained on are generated.

    The generators that start from input compounds run REINVENT once per input compound, so every output keeps the id of the compound it came from and one input that REINVENT rejects fails alone. The prior is loaded again for each input.

    Inputs:
        generator:      "reinvent", "libinvent", "linkinvent" or "mol2mol" (default "reinvent").
        compounds:      the input compounds, read from the `smiles` column of the compounds map_table. Required for libinvent, linkinvent and mol2mol; not accepted by reinvent.
        num_compounds:  molecules sampled per input compound, or in total for reinvent (default 100; upstream `num_smiles`). Duplicates and invalid SMILES are removed by REINVENT, so fewer may come back. For mol2mol with beam search it is the beam size.
        similarity:     mol2mol only: which prior, that is which kind of similarity the generated analogues have to the input (default "medium_similarity"), as the REINVENT 4 paper describes the priors:
                        "similarity": Tanimoto >= 0.5;
                        "medium_similarity": 0.5 <= Tanimoto < 0.7;
                        "high_similarity": Tanimoto >= 0.7;
                        "scaffold": same Murcko scaffold;
                        "scaffold_generic": same unlabelled Murcko scaffold;
                        "mmp": matched molecular pairs;
                        "pubchem": trained on PubChem pairs with Tanimoto >= 0.5 on ECFP4 count fingerprints, with a ranking loss that ties the NLL to the similarity.
                        The first six are trained on ChEMBL.
        sample_strategy: mol2mol only: "multinomial" (stochastic, upstream default) or "beamsearch" (deterministic).
        temperature:    mol2mol only: temperature of multinomial sampling, > 0 (upstream default 1.0); below 1 is more deterministic, above 1 more random.
        seed:           random seed passed to REINVENT (`--seed`), so a run can be repeated; an integer >= 1, since upstream treats 0 as unset. Unset leaves it random.

    Outputs:
        Streams:
            compounds:  value-based compounds stream of the generated molecules.
                        ids are <input id>_<n> for libinvent, linkinvent and mol2mol and <pipeline name>_<n> for reinvent, n counting from 1 in the order REINVENT wrote them.
        Tables:
            compounds:  id | format | code | smiles | ccd | compounds.id
                        (the stream's map_table; format is "smiles", code "LIG", ccd empty; compounds.id, the input compound's id, is absent for reinvent)
            samples:    id | compounds.id | smiles | nll | <generator columns>
                        (one row per generated molecule; compounds.id is absent for reinvent. nll is the negative log-likelihood of the molecule under the prior, rounded to two decimals by REINVENT: lower means the prior finds it more probable.
                        libinvent adds scaffold and r_groups (the scaffold as given to the model and the R-groups it generated), linkinvent adds warheads and linker, mol2mol adds input_smiles (the input as given to the model) and tanimoto (REINVENT's Tanimoto similarity of the analogue to its input).)
            inputs:     id | smiles | num_compounds
                        (libinvent, linkinvent and mol2mol only; one row per input compound; num_compounds is how many molecules it yielded, empty if its run failed)
            missing:    id | removed_by | kind | cause
                        (an input that yielded no molecule is kind="filter"; an input with no SMILES, or one REINVENT could not run, is kind="failure")

    Install:
        REINVENT4.install()                        # PyTorch build for CUDA 12.6
        REINVENT4.install(processor="cpu")         # any tag upstream's install.py takes: cu126, cu128, cpu, rocm6.4, xpu, mac
        Clones MolecularAI/REINVENT4 and installs it with its own install.py without the optional scoring dependencies, then downloads the priors from Zenodo record 15641297 and checks each md5.

    Usage:
        with Pipeline(...):
            Resources(gpu="any", memory="16GB", time="2:00:00")
            scaffolds = CompoundLibrary({"s1": "[*:0]Cc2ccc1cncc(C[*:1])c1c2"})
            decorated = REINVENT4(generator="libinvent", compounds=scaffolds, num_compounds=200)

            lowest = Panda(tables=decorated.tables.samples,
                           operations=[Panda.sort("nll"), Panda.head(50)],
                           pool=decorated)
    """

    TOOL_NAME = "REINVENT4"
    TOOL_VERSION = "1.2"
    ENV_NAME = "REINVENT4"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False,
                        processor=DEFAULT_PROCESSOR, **kwargs):
        """Create the env, install REINVENT 4 from a clone with upstream's install.py, fetch the priors.

        install.py runs pip without checking its exit status, so the install is verified by `reinvent --version`, which imports torch and rdkit, and by importing pandas, which the pipe script needs. An existing env is kept, so a failed download is retried by calling install() again; force_reinstall rebuilds the env and re-runs install.py on the existing clone.
        """
        if not isinstance(processor, str) or not _PROCESSOR_RE.match(processor):
            raise ValueError(f"processor must be a PyTorch build tag such as 'cu126' or 'cpu', got {processor!r}")
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("REINVENT4", "")
        prior_dir = f"{repo_dir}/priors"
        env_check = cls._env_exists_check(env, env_manager)
        env_run = cls._env_run(env, env_manager)
        cli_check = f'{env_run}reinvent --version >/dev/null 2>&1 && {env_run}python -c "import pandas" >/dev/null 2>&1'
        priors_check = " && ".join(f'[ -s "{prior_dir}/{name}" ]' for name in PRIOR_MD5)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check} \\
   && {cli_check} \\
   && {priors_check}; then
    echo "REINVENT4 already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        fetch_function = _FETCH_PRIOR.replace("@PRIOR_DIR@", prior_dir).replace("@PRIORS_URL@", PRIORS_URL)
        fetch_calls = "\n".join(f'fetch_prior "{name}" "{md5}" || exit 1' for name, md5 in PRIOR_MD5.items())
        return f"""echo "=== Installing REINVENT4 ==="
{skip}{remove_block}
mkdir -p "$(dirname "{repo_dir}")"
if [ ! -f "{repo_dir}/install.py" ]; then
    git clone --depth 1 https://github.com/MolecularAI/REINVENT4.git "{repo_dir}" || exit 1
fi

if {env_check}; then
    echo "{env} environment already present, not recreating it"
else
{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create {env} environment."
    exit 1
fi
fi

cd "{repo_dir}" || exit 1
{env_run}python install.py {processor} -d none

mkdir -p "{prior_dir}"
{fetch_function}
{fetch_calls}

# Verify: reinvent starts (imports torch and rdkit) AND every prior is in place.
if {cli_check} \\
   && {priors_check}; then
    touch "$INSTALL_SUCCESS"
    echo "=== REINVENT4 installation complete ==="
else
    echo "ERROR: REINVENT4 verification failed (reinvent does not start, pandas does not import, or a prior is missing). Retry with force_reinstall=True to rebuild the env."
    exit 1
fi
"""

    config_json = Path(lambda self: self.configuration_path("reinvent4_config.json"))
    compounds_json = Path(lambda self: self.configuration_path("compounds_ds.json"))
    sampling_dir = Path(lambda self: self.execution_path("sampling"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    compounds_csv = Path(lambda self: self.stream_path("compounds", "compounds.csv"))
    samples_csv = Path(lambda self: self.table_path("samples"))
    inputs_csv = Path(lambda self: self.table_path("inputs"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_reinvent4.py"))

    def __init__(self,
                 generator: str = "reinvent",
                 compounds: Optional[Union[DataStream, StandardizedOutput]] = None,
                 num_compounds: int = 100,
                 similarity: Optional[str] = None,
                 sample_strategy: Optional[str] = None,
                 temperature: Optional[float] = None,
                 seed: Optional[int] = None,
                 **kwargs):
        self.generator = generator
        self.compounds = compounds  # raw handle kept for missing-propagation
        if compounds is None:
            self.compounds_stream: Optional[DataStream] = None
        elif isinstance(compounds, StandardizedOutput):
            self.compounds_stream = compounds.streams.compounds
        elif isinstance(compounds, DataStream):
            self.compounds_stream = compounds
        else:
            raise ValueError(
                f"compounds must be DataStream or StandardizedOutput, got {type(compounds)}"
            )
        self.num_compounds = num_compounds
        self.similarity = similarity
        self.sample_strategy = sample_strategy
        self.temperature = temperature
        self.seed = seed
        super().__init__(**kwargs)

    @property
    def conditional(self) -> bool:
        return self.generator in CONDITIONAL_GENERATORS

    def validate_params(self):
        if self.generator not in GENERATORS:
            raise ValueError(f"generator must be one of {list(GENERATORS)}, got {self.generator!r}")

        if self.conditional:
            if self.compounds_stream is None or len(self.compounds_stream) == 0:
                raise ValueError(
                    f"generator={self.generator!r} needs compounds (the molecules to sample around), "
                    f"and it must not be empty")
        elif self.compounds_stream is not None:
            raise ValueError(
                "generator='reinvent' designs de novo and takes no compounds; "
                "use libinvent, linkinvent or mol2mol to sample around given molecules")

        if not _is_int(self.num_compounds) or self.num_compounds < 1:
            raise ValueError(f"num_compounds must be an integer >= 1, got {self.num_compounds!r}")
        if self.seed is not None and (not _is_int(self.seed) or self.seed < 1):
            raise ValueError(f"seed must be an integer >= 1 (upstream treats 0 as unset), got {self.seed!r}")

        for name, value in (("similarity", self.similarity),
                            ("sample_strategy", self.sample_strategy),
                            ("temperature", self.temperature)):
            if value is not None and self.generator != "mol2mol":
                raise ValueError(f"{name} only applies to generator='mol2mol', got generator={self.generator!r}")
        if self.similarity is not None and self.similarity not in SIMILARITY_PRIOR:
            raise ValueError(f"similarity must be one of {list(SIMILARITY_PRIOR)}, got {self.similarity!r}")
        if self.sample_strategy is not None and self.sample_strategy not in SAMPLE_STRATEGIES:
            raise ValueError(
                f"sample_strategy must be one of {list(SAMPLE_STRATEGIES)}, got {self.sample_strategy!r}")
        if self.temperature is not None and (
                isinstance(self.temperature, bool) or not isinstance(self.temperature, (int, float))
                or self.temperature <= 0):
            raise ValueError(f"temperature must be a number > 0, got {self.temperature!r}")
        if self.sample_strategy == "beamsearch" and self.temperature is not None:
            raise ValueError("temperature applies to multinomial sampling only; drop it for sample_strategy='beamsearch'")

    @staticmethod
    def compound_columns(generator: str) -> List[str]:
        """Columns of the compounds stream's map_table."""
        if generator in CONDITIONAL_GENERATORS:
            return COMPOUND_COLUMNS + [PROVENANCE_COLUMN]
        return list(COMPOUND_COLUMNS)

    @staticmethod
    def sample_columns(generator: str) -> List[str]:
        """Columns of the samples table: id, input provenance, smiles, nll, then what upstream reports for the generator."""
        columns = ["id", "smiles", "nll"]
        if generator in CONDITIONAL_GENERATORS:
            columns = ["id", PROVENANCE_COLUMN, "smiles", "nll"]
        if generator == "libinvent":
            return columns + ["scaffold", "r_groups"]
        if generator == "linkinvent":
            return columns + ["warheads", "linker"]
        if generator == "mol2mol":
            return columns + ["input_smiles", "tanimoto"]
        return columns

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"GENERATOR: {self.generator}")
        if self.conditional:
            lines.append(f"INPUTS: {len(self.compounds_stream)} compound(s)")
            lines.append(f"COMPOUNDS PER INPUT: {self.num_compounds}")
        else:
            lines.append(f"COMPOUNDS: {self.num_compounds}")
        if self.generator == "mol2mol":
            lines.append(f"SIMILARITY: {self.similarity or DEFAULT_SIMILARITY}")
            if self.sample_strategy is not None:
                lines.append(f"SAMPLE STRATEGY: {self.sample_strategy}")
            if self.temperature is not None:
                lines.append(f"TEMPERATURE: {self.temperature}")
        if self.seed is not None:
            lines.append(f"SEED: {self.seed}")
        return lines

    def _pipe_config(self) -> Dict[str, Any]:
        """The JSON the pipe script reads."""
        repository = self.folders.get("REINVENT4")
        if not repository:
            raise ValueError(
                "folders.repositories.REINVENT4 is not set in the active config; REINVENT4.install() "
                "clones the repository there and downloads the priors into its priors/ folder")
        return {
            "generator": self.generator,
            "model_file": os.path.join(repository, "priors", prior_file(self.generator, self.similarity)),
            "num_compounds": self.num_compounds,
            "sample_strategy": self.sample_strategy,
            "temperature": self.temperature,
            "seed": self.seed,
            "compounds_json": self.compounds_json if self.conditional else None,
            "id_prefix": None if self.conditional else self.pipeline_name,
            "sampling_dir": self.sampling_dir,
            "compounds_csv": self.compounds_csv,
            "samples_csv": self.samples_csv,
            "inputs_csv": self.inputs_csv if self.conditional else None,
            "local_missing_csv": self.local_missing_csv,
        }

    def generate_script(self, script_path: str) -> str:
        os.makedirs(os.path.dirname(self.config_json), exist_ok=True)
        if self.conditional:
            self.compounds_stream.save_json(self.compounds_json)
        with open(self.config_json, "w", encoding="utf-8") as f:
            json.dump(self._pipe_config(), f, indent=2)

        script = "#!/bin/bash\n"
        script += "# REINVENT4 sampling script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        # The helper stays on the host; only the reinvent subprocesses enter a container.
        script += f"""echo "Sampling with REINVENT4 generator {self.generator}"
python "{self.helper_py}" \\
    --config "{self.config_json}" \\
    --container-prefix "{self.container_prefix()}"
"""
        script += self.generate_missing_propagation(
            self.compounds, local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        # How many molecules survive REINVENT's validity and duplicate filter is only known at run time.
        if self.conditional:
            ids = generate_multiplied_ids_pattern(
                self.compounds_stream.ids, "[<#>]", input_stream_name="compounds")
        else:
            ids = [f"{self.pipeline_name}_[<#>]"]

        compounds_stream = DataStream(
            name="compounds",
            ids=ids,
            files=[],
            map_table=self.compounds_csv,
            format="csv",
        )

        tables = {
            "compounds": TableInfo(
                name="compounds",
                path=self.compounds_csv,
                columns=self.compound_columns(self.generator),
                description="Generated molecules (the compounds stream's content table)",
            ),
            "samples": TableInfo(
                name="samples",
                path=self.samples_csv,
                columns=self.sample_columns(self.generator),
                description="Negative log-likelihood and generator-specific columns per generated molecule",
            ),
        }
        if self.conditional:
            tables["inputs"] = TableInfo(
                name="inputs",
                path=self.inputs_csv,
                columns=list(INPUT_COLUMNS),
                description="Number of molecules each input compound yielded",
            )
        tables["missing"] = self.missing_table_info(self.missing_csv)

        return {
            "compounds": compounds_stream,
            "tables": tables,
            "output_folder": self.output_folder,
        }
