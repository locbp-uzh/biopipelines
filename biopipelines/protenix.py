# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

import json
import os
from typing import Dict, List, Any, Union, Optional

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import (generate_combinatorics_config,
                                predict_output_ids_with_provenance)
    from .datastream_resolver import resolve_input_to_datastream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import (generate_combinatorics_config,
                               predict_output_ids_with_provenance)
    from datastream_resolver import resolve_input_to_datastream


# Checkpoints `protenix pred -n` accepts that need only the sequences; the constraint and ESM ones need more.
MODEL_NAMES = (
    "protenix-v2",
    "protenix_base_default_v1.0.0",
    "protenix_base_20250630_v1.0.0",
    "protenix_base_default_v0.5.0",
    "protenix_mini_default_v0.5.0",
    "protenix_tiny_default_v0.5.0",
)
DEFAULT_MODEL = "protenix-v2"
# Exit status of pipe_protenix_weights.py for a checkpoint that could not be downloaded.
CHECKPOINT_UNAVAILABLE = 3

# Upstream's defaults, to declare how many structures a query writes.
UPSTREAM_SEEDS = 1
UPSTREAM_SAMPLES = 5

COVALENT_BOND_KEYS = ("entity1", "copy1", "position1", "atom1",
                      "entity2", "copy2", "position2", "atom2")
COVALENT_BOND_REQUIRED = ("entity1", "position1", "atom1", "entity2", "position2", "atom2")


class Protenix(BaseConfig):
    """
    Protenix: ByteDance's trainable reproduction of AlphaFold3, co-folding proteins, nucleic acids and ligands.

    Defaults to `protenix-v2` (about 464 M parameters).
    The upstream README states that the weights of `protenix-v2` are proprietary, are not under an open-source license and may not be redistributed or transferred to a third party without the rights holder's written consent; the code is Apache-2.0.
    Upstream's download URL for that checkpoint refused the request (HTTP 403) when checked on 2026-10-28, where the other checkpoints downloaded.
    The install then still succeeds, with a warning naming where the file goes; any other `model_name` is unaffected.
    It takes the same axes as `Boltz2` and `OpenFold3`, so a prediction can be cross-checked against another model by changing the tool name.
    There is no affinity head, and covalent bonds are explicit atom pairs (`covalent_bonds`), not a selection.
    Contact and pocket constraints exist upstream only for the older constraint checkpoint, and templates need kalign plus per-chain search results or a sequence database; neither is wired here.

    Inputs:
        proteins: protein sequences. Letters other than the 20 standard residues and X fail that query.
        ssDNA / dsDNA / ssRNA / dsRNA: nucleic-acid sequences in upper case, A/C/G/T/N for DNA and A/C/G/U/N for RNA. A double-stranded axis becomes two chains, the second the reverse complement.
        ligands: compounds, by SMILES or CCD code. A single-atom ion is a ligand with its CCD code (`MG`), which upstream builds the same way as its `ion` entity.
        msas: precomputed alignments (`.a3m`) for the protein chains, from `MMseqs2(output_format="a3m")` or `MSA`.
            Every protein chain of a query needs one; a query with a chain that has none is not predicted and lands in `missing`.
            They are passed as unpaired alignments; upstream's ColabFold note says Protenix pairs chains through species information that ColabFold-style alignments lack.
        use_msa_server: search alignments through upstream's MSA service at protenix-server.com, to which the sequences are sent (default: on, unless `msas` is given).
            False without `msas` folds every chain single-sequence.
            True together with `msas` is refused, since upstream then searches every chain of the query and overwrites the supplied alignments.
        model_name: checkpoint, one of MODEL_NAMES (default "protenix-v2").
            Cycles and diffusion steps are set per model by upstream (10 and 200 for the base models and v2, 4 and 5 for mini and tiny).
            `protenix-v2` refuses a complex of more than 2560 tokens, and that query lands in `missing`.
        num_diffusion_samples: structures sampled per seed (upstream default 5).
        seeds: random seeds (upstream default 101). The structures of a query are the seeds times the samples.
        top_only: keep only the structure with the best ranking score as `<id>` (default).
            False surfaces every structure as `<id>_1..N`, ordered by seed and then by upstream's sample index.
        covalent_bonds: covalent bonds between entities, as upstream's list of dicts with `entity1`, `position1`, `atom1`, `entity2`, `position2`, `atom2` and optionally `copy1` and `copy2` (both or neither).
            Entities are numbered from 1 in the order the query lists them: proteins, ssDNA, dsDNA (strand, then reverse complement), ssRNA, dsRNA, ligands, each axis in the order of its records.
            The same bonds apply to every query, so an entity number must mean the same thing in all of them.
        include_pae: write the PAE matrix of every kept structure to the `pae` stream (default True).
            It makes upstream dump a full-data JSON per sample, which is large; the wrapper reads the matrix out of it and deletes the JSONs.
        no_kernels: run the triangle operations with plain PyTorch instead of the cuEquivariance kernels, for a GPU those do not support.
        The layer norm always runs in plain PyTorch (`LAYERNORM_TYPE=torch`), because upstream's default fused layer norm is compiled at the first import and needs nvcc and CUDA_HOME, which a node with only the pip CUDA wheels lacks; export `LAYERNORM_TYPE=fast_layernorm` beforehand to compile it.

    Outputs:
        Streams:
            structures: one `<id>.cif` per predicted complex.
            pae:        one `<id>.npz` per structure, with the same ids as `structures`.
                        The array key is `pae`, a float32 N_token x N_token matrix in Angstrom, the layout Boltz2 writes.
                        Rows and columns are Protenix's tokens over all chains: one per standard residue, one per atom for ligands and modified residues.
                        Absent when `include_pae=False`.
        Tables:
            confidence: id | plddt | ptm | iptm | gpde | has_clash | disorder | ranking_score | seed | sample
                        plddt is 0 to 100, ptm and iptm 0 to 1.
                        `seed` and `sample` identify the structure the row describes; `sample` is upstream's index within the seed, 0 for the highest ranking score of that seed.
            missing:    id | removed_by | kind | cause
    """

    TOOL_NAME = "Protenix"
    TOOL_VERSION = "1.3"
    ENV_NAME = "protenix"

    PREDICT_COMMAND = "protenix pred"

    @staticmethod
    def _import_free_check(run):
        """Shell test that the package and its entry point are present without importing protenix.

        Importing protenix JIT-compiles a CUDA extension and fails on a node with no CUDA toolkit, so even `protenix --version` cannot verify an install made on a CPU node.
        """
        metadata = "import importlib.metadata as m; m.version('protenix')"
        entry_point = "import shutil, sys; sys.exit(0 if shutil.which('protenix') else 1)"
        return f'{run}python -c "{metadata}" >/dev/null 2>&1 && {run}python -c "{entry_point}" >/dev/null 2>&1'

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        run = cls._env_run(env, env_manager)
        cache = folders.get("ProtenixCache", "")
        root = f'export PROTENIX_ROOT_DIR="{cache}"\nmkdir -p "$PROTENIX_ROOT_DIR"\n' if cache else ""
        fetch_py = os.path.join(folders.get("pipe_scripts") or os.path.join(biopipelines, "pipe_scripts"),
                                "pipe_protenix_weights.py")
        # Upstream downloads these on the first prediction, which races between jobs and needs the network on a compute node.
        # Status 3 is a checkpoint upstream's URL refused: the env and the CCD files are fine, so it warns rather than failing the install.
        weights = f"""
echo "Ensuring the Protenix checkpoint and the CCD files"
{root}{run}python "{fetch_py}" --model-name "{DEFAULT_MODEL}"
FETCH_RC=$?
if [ "$FETCH_RC" -eq {CHECKPOINT_UNAVAILABLE} ]; then
    echo "WARNING: Protenix is installed without the {DEFAULT_MODEL} checkpoint; see above for where it goes"
elif [ "$FETCH_RC" -ne 0 ]; then
    echo "ERROR: the CCD files could not be downloaded"
    exit 1
fi
"""
        skip = "" if force_reinstall else f"""if {cls._env_exists_check(env, env_manager)} && {cls._import_free_check(run)}; then
    echo "Protenix already installed, skipping the env. Use force_reinstall=True to rebuild it."
{weights}
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        return f"""echo "=== Installing Protenix ==="
{skip}{remove_block}
{env_block}
{weights}
if {cls._import_free_check(run)}; then
    touch "$INSTALL_SUCCESS"
    echo "=== Protenix installation complete ==="
else
    echo "ERROR: Protenix verification failed"
    exit 1
fi
"""

    # Configuration-time artifacts
    combinatorics_config_file = Path(lambda self: self.configuration_path("combinatorics_config.json"))
    covalent_bonds_file = Path(lambda self: self.configuration_path("covalent_bonds.json"))
    msas_json_file = Path(lambda self: self.configuration_path("msas.json"))
    # Written at execution time, one JSON per query: one model load, and a failure stays in its own file.
    queries_folder = Path(lambda self: self.configuration_path("queries"))
    query_ids_csv = Path(lambda self: self.configuration_path("query_ids.csv"))
    skipped_csv = Path(lambda self: self.configuration_path("skipped_queries.csv"))

    prediction_folder = Path(lambda self: self.execution_folder)

    structures_map_csv = Path(lambda self: self.stream_map_path("structures"))
    pae_folder = Path(lambda self: self.stream_folder("pae"))
    pae_map_csv = Path(lambda self: self.stream_map_path("pae"))
    confidence_csv = Path(lambda self: self.table_path("confidence"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))

    config_py = Path(lambda self: self.pipe_script_path("pipe_protenix_config.py"))
    postprocess_py = Path(lambda self: self.pipe_script_path("pipe_protenix_postprocess.py"))

    def __init__(self,
                 proteins: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ssDNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 dsDNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ssRNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 dsRNA: Optional[Union[DataStream, StandardizedOutput]] = None,
                 ligands: Optional[Union[DataStream, StandardizedOutput]] = None,
                 msas: Optional[Union[DataStream, StandardizedOutput]] = None,
                 use_msa_server: Optional[bool] = None,
                 model_name: str = DEFAULT_MODEL,
                 num_diffusion_samples: Optional[int] = None,
                 seeds: Optional[List[int]] = None,
                 top_only: bool = True,
                 covalent_bonds: Optional[List[Dict[str, Any]]] = None,
                 include_pae: bool = True,
                 no_kernels: bool = False,
                 **kwargs):
        # Raw handles for missing-propagation, resolved streams for the config.
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

        # Unset means the server, unless precomputed MSAs were given.
        self.use_msa_server = (msas is None) if use_msa_server is None else use_msa_server
        self.model_name = model_name
        self.num_diffusion_samples = num_diffusion_samples
        self.seeds = seeds
        self.top_only = top_only
        self.covalent_bonds = covalent_bonds
        self.include_pae = include_pae
        self.no_kernels = no_kernels

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
                "Protenix needs at least one input axis: proteins, ssDNA, dsDNA, ssRNA, "
                "dsRNA or ligands")

        if self.model_name not in MODEL_NAMES:
            raise ValueError(
                f"model_name must be one of {MODEL_NAMES}, got {self.model_name!r}")

        value = self.num_diffusion_samples
        if value is not None and (not isinstance(value, int) or isinstance(value, bool) or value < 1):
            raise ValueError("num_diffusion_samples must be a positive integer")

        if self.seeds is not None:
            if not isinstance(self.seeds, (list, tuple)) or not self.seeds:
                raise ValueError("seeds must be a non-empty list of integers")
            if any(not isinstance(s, int) or isinstance(s, bool) for s in self.seeds):
                raise ValueError("seeds must contain integers only")

        for name in ("top_only", "include_pae", "no_kernels"):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be True or False")

        if self.msas is not None and self.use_msa_server:
            raise ValueError(
                "msas supplies precomputed alignments, so use_msa_server must be False; "
                "Protenix would search every protein chain of a query through the server "
                "and overwrite the supplied alignments")

        if self.msas is not None:
            fmt = getattr(getattr(self, "msas_stream", None), "format", None)
            if fmt is not None and str(fmt).lower() != "a3m":
                raise ValueError(f"Protenix reads MSAs as a3m; this msas stream is {fmt}. Use "
                                 "MMseqs2(output_format=\"a3m\") or MSA(..., convert=\"a3m\")")

        if self.covalent_bonds is not None:
            self._validate_covalent_bonds(self.covalent_bonds)

    @staticmethod
    def _validate_covalent_bonds(bonds):
        if not isinstance(bonds, (list, tuple)) or not bonds:
            raise ValueError("covalent_bonds must be a non-empty list of dicts")
        for index, bond in enumerate(bonds):
            label = f"covalent_bonds[{index}]"
            if not isinstance(bond, dict):
                raise ValueError(f"{label} must be a dict, got {type(bond).__name__}")
            unknown = sorted(set(bond) - set(COVALENT_BOND_KEYS))
            if unknown:
                raise ValueError(f"{label}: unknown key(s) {unknown}; expected {list(COVALENT_BOND_KEYS)}")
            missing = [key for key in COVALENT_BOND_REQUIRED if bond.get(key) in (None, "")]
            if missing:
                raise ValueError(f"{label}: missing key(s) {missing}")
            if (bond.get("copy1") is None) != (bond.get("copy2") is None):
                raise ValueError(f"{label}: copy1 and copy2 must be given together or not at all")
            for key in ("entity1", "entity2", "position1", "position2", "copy1", "copy2"):
                if bond.get(key) is None:
                    continue
                try:
                    number = int(bond[key])
                except (TypeError, ValueError):
                    raise ValueError(f"{label}: {key} must be a positive integer, got {bond[key]!r}")
                if isinstance(bond[key], bool) or number < 1:
                    raise ValueError(f"{label}: {key} counts from 1, got {bond[key]!r}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        axes = [name for name in ("proteins", "ssDNA", "dsDNA", "ssRNA", "dsRNA", "ligands")
                if getattr(self, name) is not None]
        lines.append(f"AXES: {', '.join(axes)}")
        lines.append(f"MODEL: {self.model_name}")
        if self.msas is not None:
            lines.append("MSA: precomputed")
        elif self.use_msa_server:
            lines.append("MSA: upstream service")
        else:
            lines.append("MSA: none (single sequence)")
        if self.num_diffusion_samples is not None:
            lines.append(f"DIFFUSION SAMPLES: {self.num_diffusion_samples}")
        if self.seeds is not None:
            lines.append(f"SEEDS: {', '.join(str(s) for s in self.seeds)}")
        if self.covalent_bonds:
            lines.append(f"COVALENT BONDS: {len(self.covalent_bonds)}")
        lines.append(f"PAE: {'written' if self.include_pae else 'not written'}")
        if self.no_kernels:
            lines.append("KERNELS: plain PyTorch")
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

    def _structures_per_query(self) -> int:
        seeds = len(self.seeds) if self.seeds else UPSTREAM_SEEDS
        return seeds * (self.num_diffusion_samples or UPSTREAM_SAMPLES)

    def _cache_root(self) -> Optional[str]:
        return (getattr(self, "folders", None) or {}).get("ProtenixCache") or None

    def generate_script(self, script_path: str) -> str:
        generate_combinatorics_config(self.combinatorics_config_file, **self._axis_kwargs())
        if self.covalent_bonds:
            with open(self.covalent_bonds_file, "w") as handle:
                json.dump([dict(bond) for bond in self.covalent_bonds], handle, indent=2)

        script = "#!/bin/bash\n"
        script += "# Protenix script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        if self._cache_root():
            script += f'export PROTENIX_ROOT_DIR="{self._cache_root()}"\n\n'
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
        """Build the per-query JSONs at execution time, when the sequences actually exist.

        The wrapper cannot write them: the sequences live in upstream streams whose content is
        produced by the step before this one.
        """
        flags = ""
        if self.msas_stream:
            self.msas_stream.save_json(self.msas_json_file)
            flags += f' \\\n    --msas-json "{self.msas_json_file}"'
        if self.covalent_bonds:
            flags += f' \\\n    --covalent-bonds-json "{self.covalent_bonds_file}"'
        return f"""echo "Building Protenix query JSONs"
python "{self.config_py}" \\
    --combinatorics-config "{self.combinatorics_config_file}" \\
    --queries-dir "{self.queries_folder}" \\
    --query-ids-csv "{self.query_ids_csv}" \\
    --skipped-csv "{self.skipped_csv}"{flags}

"""

    def _generate_predict_section(self) -> str:
        use_msa = self.msas is not None or self.use_msa_server
        options = [f'-i "{self.queries_folder}"',
                   f'-o "{self.prediction_folder}"',
                   f"-n {self.model_name}",
                   "--use_default_params true",
                   f"--use_msa {'true' if use_msa else 'false'}"]
        if self.seeds is not None:
            options.append(f"-s {','.join(str(s) for s in self.seeds)}")
        if self.num_diffusion_samples is not None:
            options.append(f"-e {self.num_diffusion_samples}")
        if self.include_pae:
            options.append("--need_atom_confidence true")
        if self.no_kernels:
            options.append("--trimul_kernel torch --triatt_kernel torch")
        cp = self.container_prefix()
        return f"""echo "Running Protenix prediction"
# Upstream's default fused layer norm JIT-compiles with nvcc at the first import; the PyTorch one needs no CUDA toolkit
export LAYERNORM_TYPE="${{LAYERNORM_TYPE:-torch}}"
{cp}{self.PREDICT_COMMAND} {" ".join(options)}

"""

    def _generate_postprocess_section(self) -> str:
        pae_flags = ""
        if self.include_pae:
            pae_flags = f'    --pae-dir "{self.pae_folder}" \\\n    --pae-map-csv "{self.pae_map_csv}" \\\n'
        return f"""echo "Collecting Protenix predictions"
python "{self.postprocess_py}" \\
    --prediction-folder "{self.prediction_folder}" \\
    --query-ids-csv "{self.query_ids_csv}" \\
    --skipped-csv "{self.skipped_csv}" \\
    --structures-dir "{self.stream_folder('structures')}" \\
    --structures-map-csv "{self.structures_map_csv}" \\
    --confidence-csv "{self.confidence_csv}" \\
{pae_flags}    --local-missing-csv "{self.local_missing_csv}" \\
    --top-only "{str(self.top_only)}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        predicted_ids, provenance = predict_output_ids_with_provenance(**self._axis_kwargs())

        # With top_only off every structure is <id>_1..N, and declaring the bare id would fail the completion check.
        if self.top_only:
            structure_ids = predicted_ids
        else:
            structure_ids = [f"{i}_<1..{self._structures_per_query()}>" for i in predicted_ids]

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=[self.stream_path("structures", "<id>.cif")],
            map_table=self.structures_map_csv,
            format="cif",
        )

        pae = DataStream.empty("pae", "npz")
        if self.include_pae:
            pae = DataStream(
                name="pae",
                ids=structure_ids,
                files=[self.stream_path("pae", "<id>.npz")],
                map_table=self.pae_map_csv,
                format="npz",
            )

        tables = {
            "confidence": TableInfo(
                name="confidence",
                path=self.confidence_csv,
                columns=["id", "plddt", "ptm", "iptm", "gpde", "has_clash",
                         "disorder", "ranking_score", "seed", "sample"],
                description="Per-complex confidence scores from Protenix",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }

        return {
            "structures": structures,
            "pae": pae,
            "tables": tables,
            "output_folder": self.output_folder,
        }
