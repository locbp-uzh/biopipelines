# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""BioEmu2 tool: ligand-capable equilibrium-ensemble emulation.

bioemu2 is the ligand-capable successor of bioemu. Each run samples
``num_samples`` conformers of a system that holds protein chains and
(optionally) ligands, plus optional contact and covalent-bond constraints.
The wrapper writes one system YAML per run, calls ``bioemu2
sample-equilibrium``, splits the sampled trajectory into per-conformer PDBs
and keeps the compact trajectory and the physicality report.

Reference: https://github.com/microsoft/bioemu2
"""

import os
from typing import Any, Dict, List, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import (generate_combinatorics_config, generate_multiplied_ids_pattern,
                                predict_output_ids_with_provenance)
    from .datastream_resolver import resolve_input_to_datastream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import (generate_combinatorics_config, generate_multiplied_ids_pattern,
                               predict_output_ids_with_provenance)
    from datastream_resolver import resolve_input_to_datastream


# The values `bioemu2 sample-equilibrium --sampler` and `--physicality-filter` accept.
SAMPLERS = ("edm", "edm_gd", "edm_fk", "edm_fkc")
PHYSICALITY_MODES = ("high_quality", "high_quality_per_residue", "no_violations", "no_filter")
MSA_FORMATS = ("a3m", "csv")
COVALENT_LINKAGE_KEYS = ("chain", "position", "protein_atom", "ligand_atom")
CONTACT_TOKEN_KEYS = ("asym_id", "res_id", "atom_name")
# bioemu2 declares requires-python >=3.12; the env spec, the install check and the tests all read this.
MIN_PYTHON = (3, 12)
PYTHON_CHECK = f"import sys; sys.exit(0 if sys.version_info >= {MIN_PYTHON} else 1)"


class BioEmu2(BaseConfig):
    """
    BioEmu2: ligand-capable equilibrium-ensemble emulation.

    Each combination of a sequence and a ligand is one run, and one system for bioemu2.
    The protein chains are A, B, ... in the order the sequences are given, and the ligands take the next chain letters.
    A ligand is passed to bioemu2 by CCD code when its compounds row has one and by SMILES otherwise, and bioemu2 places it from scratch: no 3D template is seeded.
    A compounds row with neither (a `Ligand(codes=...)` that only names a residue) fails its run.

    Inputs:
        sequences: protein sequences, as a sequences stream or a StandardizedOutput. Each row is one chain.
                   `Each` (the default) makes one run per sequence; `Bundle` puts several sequences in one system.
        ligands: optional compounds stream or StandardizedOutput. Without it the ensemble is apo.
                   `Each` makes one run per ligand and sequence; `Bundle` puts several ligands in one system.
        num_samples: conformers per run (default 100).
        batch_size: samples denoised together per forward pass (default 16). It caps GPU memory at long sequence times large batch.
        seed: RNG seed (default None: every run draws its own).
        sampler: `edm` (default, unsteered), `edm_gd` (gradient guidance), `edm_fk` (Feynman-Kac resampling) or `edm_fkc` (Feynman-Kac corrector).
                   The steered samplers use upstream's default physicality potentials.
        contacts: optional list of upstream contact dicts, `{"token1": {"asym_id": "A", "res_id": 3}, "token2": {"asym_id": "B", "atom_name": "C1"}, "max_distance": 6.0}`.
                   A token is a chain with a 1-based residue index, or a chain with a ligand atom name. They go into the system YAML unchanged.
        covalent_linkage: optional dict `{"chain": "A", "position": 3, "protein_atom": "SG", "ligand_atom": "C1"}` bonding a protein atom to an atom of the first ligand.
                   It becomes an upstream `bonds` entry, with the residue name read from the sequence. The ligand must be given by CCD code, because the bond names the ligand's residue and atoms.
                   Requires `ligands`.
        msas: optional precomputed alignments (msas stream, `.a3m` or Boltz-style `.csv`), matched to a sequence by its id (the `sequences.id` column, else the row id).
                   The alignment is written into bioemu2's ColabFold cache for this step, so no search is made for that sequence. A sequence with no alignment is searched through ColabFold as usual.
                   Only a system with one protein chain can use it. Formerly `precomputed_msa`.
        physicality_filter: `high_quality` (default), `high_quality_per_residue`, `no_violations` or `no_filter`.
                   Frames that fail stay in the trajectory and the structures; the mode only decides what the physicality table reports.

    Install:
        BioEmu2.install()   # Python 3.12 env, bioemu2[cuda], and the bioemu2-preview checkpoint from Hugging Face

    Outputs:
        Streams:
            structures:   one `<run>_<1..num_samples>.pdb` per sampled conformer.
            trajectories: one `<run>.xtc` per run, with its `<run>_topology.cif` named in the map's `topology` column.
        Tables:
            summary:     id | sequences.id | ligands.id | n_samples | n_residues
                         One row per run, with n_samples and n_residues empty when the run failed. n_residues counts every residue of the sampled system, ligands included.
            physicality: id | run.id | frame | status | reason
                         One row per conformer, with `id` the structure id and `frame` its zero-based index in the trajectory.
                         `status` is `passed`, `rejected` or `not_checked`, and `reason` the first failed check (empty when none).
            missing:     id | removed_by | kind | cause
    """

    TOOL_NAME = "BioEmu2"
    TOOL_VERSION = "1.2"
    ENV_NAME = "bioemu2"

    PARAMETER_ALIASES = {"precomputed_msa": "msas"}
    DEPRECATED_ALIASES = ("precomputed_msa",)

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        run = cls._env_run(env, env_manager)
        env_check = cls._env_exists_check(env, env_manager)
        pipe_scripts = folders.get("pipe_scripts") or f"{biopipelines}/pipe_scripts"
        weights_py = f"{pipe_scripts}/pipe_bioemu2_weights.py"
        requirements = f"{biopipelines}/environments/{env}.pip.txt"

        python_ok = f'{run}python -c "{PYTHON_CHECK}" >/dev/null 2>&1'
        installed = (f'{python_ok} && {run}python -c "import bioemu2, mdtraj, torch" >/dev/null 2>&1'
                     f' && {run}bioemu2 sample-equilibrium --help >/dev/null 2>&1')
        weights_cached = f'{run}python "{weights_py}" --check-only >/dev/null 2>&1'

        skip = "" if force_reinstall else f"""# Check if already installed (env and checkpoint)
if {installed} && {weights_cached}; then
    echo "BioEmu2 already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        stale_block = f"""# An env built from an earlier spec ran an older Python, which bioemu2 does not install on
if {env_check} && ! {python_ok}; then
    echo "The {env} environment runs Python older than {MIN_PYTHON[0]}.{MIN_PYTHON[1]}, which bioemu2 does not support; rebuilding it."
{cls._env_remove_block(env, env_manager)}
fi
"""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        return f"""echo "=== Installing BioEmu2 ==="
{skip}{remove_block}
{stale_block}
if {env_check}; then
    echo "{env} environment present, installing its requirements into it"
    {run}pip install -r "{requirements}"
else
{env_block}
fi

echo "Fetching the BioEmu2 checkpoint from Hugging Face"
if ! {run}python "{weights_py}"; then
    echo "WARNING: the BioEmu2 checkpoint could not be fetched; the first sampling run downloads it instead"
fi

# Verify installation
if {installed}; then
    touch "$INSTALL_SUCCESS"
    echo "=== BioEmu2 installation complete ==="
else
    echo "ERROR: BioEmu2 verification failed (needs Python {MIN_PYTHON[0]}.{MIN_PYTHON[1]} or newer, importable bioemu2, mdtraj and torch, and a working bioemu2 CLI)"
    exit 1
fi
"""

    combinatorics_config_file = Path(lambda self: self.configuration_path("combinatorics_config.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    msas_json = Path(lambda self: self.configuration_path("msas.json"))
    runs_folder = Path(lambda self: self.execution_path("runs"))
    msa_cache_folder = Path(lambda self: self.execution_path("msa_cache"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    trajectories_map = Path(lambda self: self.stream_map_path("trajectories"))
    summary_csv = Path(lambda self: self.table_path("summary"))
    physicality_csv = Path(lambda self: self.table_path("physicality"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_bioemu2.py"))

    def __init__(self,
                 sequences: Union[DataStream, StandardizedOutput],
                 ligands: Optional[Union[DataStream, StandardizedOutput]] = None,
                 num_samples: int = 100,
                 batch_size: int = 16,
                 seed: Optional[int] = None,
                 sampler: str = "edm",
                 contacts: Optional[List[Dict[str, Any]]] = None,
                 covalent_linkage: Optional[Dict[str, Any]] = None,
                 msas: Optional[Union[DataStream, StandardizedOutput]] = None,
                 physicality_filter: str = "high_quality",
                 **kwargs):
        # Raw handles for missing-propagation and the combinatorics config, resolved streams for validation.
        self.sequences = sequences
        self.ligands = ligands
        self.msas = msas
        self.sequences_stream = resolve_input_to_datastream(sequences, fallback_stream="sequences")
        self.ligands_stream = resolve_input_to_datastream(ligands, fallback_stream="compounds")
        self.msas_stream = resolve_input_to_datastream(msas, fallback_stream="msas")

        self.num_samples = num_samples
        self.batch_size = batch_size
        self.seed = seed
        self.sampler = sampler
        self.contacts = contacts
        self.covalent_linkage = covalent_linkage
        self.physicality_filter = physicality_filter
        super().__init__(**kwargs)

    def validate_params(self):
        if self.sequences_stream is None or not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")
        if self.ligands is not None and (not self.ligands_stream or len(self.ligands_stream) == 0):
            raise ValueError("ligands must not be empty; pass None for an apo ensemble")

        for name, minimum in (("num_samples", 1), ("batch_size", 1)):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}, got {value!r}")
        if self.seed is not None and (not isinstance(self.seed, int) or isinstance(self.seed, bool)):
            raise ValueError(f"seed must be an integer or None, got {self.seed!r}")
        if self.sampler not in SAMPLERS:
            raise ValueError(f"sampler must be one of {SAMPLERS}, got {self.sampler!r}")
        if self.physicality_filter not in PHYSICALITY_MODES:
            raise ValueError(
                f"physicality_filter must be one of {PHYSICALITY_MODES}, got {self.physicality_filter!r}")

        if self.msas is not None:
            fmt = getattr(self.msas_stream, "format", None)
            if fmt is not None and str(fmt).lower() not in MSA_FORMATS:
                raise ValueError(f"msas must be {' or '.join(MSA_FORMATS)} alignments; this stream is {fmt}. "
                                 "Use MMseqs2(output_format=\"a3m\") or MSA(..., convert=\"a3m\")")

        if self.contacts is not None:
            self._validate_contacts(self.contacts)
        if self.covalent_linkage is not None:
            self._validate_covalent_linkage(self.covalent_linkage)
            if self.ligands is None:
                raise ValueError("covalent_linkage bonds a protein atom to a ligand atom and requires ligands")

    @staticmethod
    def _validate_contacts(contacts):
        if not isinstance(contacts, (list, tuple)) or not contacts:
            raise ValueError("contacts must be a non-empty list of dicts")
        for index, contact in enumerate(contacts):
            label = f"contacts[{index}]"
            if not isinstance(contact, dict):
                raise ValueError(f"{label} must be a dict, got {type(contact).__name__}")
            for key in ("token1", "token2"):
                token = contact.get(key)
                if not isinstance(token, dict):
                    raise ValueError(f"{label} needs {key} as a dict with asym_id and res_id or atom_name")
                unknown = sorted(set(token) - set(CONTACT_TOKEN_KEYS))
                if unknown:
                    raise ValueError(f"{label}.{key}: unknown key(s) {unknown}; expected {list(CONTACT_TOKEN_KEYS)}")
                if not token.get("asym_id") or not (token.get("res_id") is not None or token.get("atom_name")):
                    raise ValueError(f"{label}.{key} needs an asym_id and a res_id or an atom_name")
            distance = contact.get("max_distance")
            if distance is not None and (isinstance(distance, bool) or not isinstance(distance, (int, float)) or distance <= 0):
                raise ValueError(f"{label}.max_distance must be a positive number of Angstrom, got {distance!r}")

    @staticmethod
    def _validate_covalent_linkage(link):
        if not isinstance(link, dict):
            raise ValueError(f"covalent_linkage must be a dict, got {type(link).__name__}")
        missing = [key for key in COVALENT_LINKAGE_KEYS if link.get(key) in (None, "")]
        if missing:
            raise ValueError(f"covalent_linkage missing required keys: {missing}")
        unknown = sorted(set(link) - set(COVALENT_LINKAGE_KEYS))
        if unknown:
            raise ValueError(f"covalent_linkage: unknown key(s) {unknown}; expected {list(COVALENT_LINKAGE_KEYS)}")
        position = link["position"]
        if not isinstance(position, int) or isinstance(position, bool) or position < 1:
            raise ValueError(f"covalent_linkage.position counts residues from 1, got {position!r}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append("LIGANDS: none (apo)" if self.ligands is None else "LIGANDS: given")
        lines.append(f"NUM SAMPLES PER RUN: {self.num_samples}")
        lines.append(f"BATCH SIZE: {self.batch_size}")
        lines.append(f"SAMPLER: {self.sampler}")
        lines.append(f"PHYSICALITY FILTER: {self.physicality_filter}")
        lines.append("MSA: precomputed" if self.msas is not None else "MSA: ColabFold (upstream)")
        if self.contacts:
            lines.append(f"CONTACTS: {len(self.contacts)}")
        if self.covalent_linkage is not None:
            link = self.covalent_linkage
            lines.append(f"COVALENT LINKAGE: chain {link['chain']} residue {link['position']} "
                         f"{link['protein_atom']} to ligand {link['ligand_atom']}")
        return lines

    def _axis_kwargs(self) -> Dict[str, Any]:
        axes = {"sequences": (self.sequences, "sequences", "protein")}
        if self.ligands is not None:
            axes["ligands"] = (self.ligands, "compounds", "ligand")
        return axes

    def _input_sources(self):
        return tuple(axis for axis in (self.sequences, self.ligands) if axis is not None)

    def _write_config_yaml(self):
        """Dump the run parameters for the pipe to read."""
        import yaml
        config = {
            "num_samples": self.num_samples,
            "batch_size": self.batch_size,
            "seed": self.seed,
            "sampler": self.sampler,
            "physicality_filter": self.physicality_filter,
            "contacts": [dict(contact) for contact in self.contacts] if self.contacts else None,
            "covalent_linkage": dict(self.covalent_linkage) if self.covalent_linkage else None,
        }
        with open(self.config_yaml, "w") as handle:
            yaml.safe_dump(config, handle, sort_keys=False)

    def generate_script(self, script_path: str) -> str:
        generate_combinatorics_config(self.combinatorics_config_file, **self._axis_kwargs())
        self._write_config_yaml()
        msas_flag = ""
        if self.msas_stream is not None:
            self.msas_stream.save_json(self.msas_json)
            msas_flag = f' \\\n    --msas-json "{self.msas_json}" \\\n    --msa-cache-dir "{self.msa_cache_folder}"'

        script = "#!/bin/bash\n"
        script += "# BioEmu2 ligand-capable equilibrium-ensemble sampling\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        # The helper imports biopipelines and stays on the host; only the bioemu2 CLI it calls takes the prefix.
        script += f"""echo "BioEmu2: sampling {self.num_samples} conformers per run"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --combinatorics-config "{self.combinatorics_config_file}" \\
    --runs-dir "{self.runs_folder}" \\
    --structures-dir "{self.stream_folder('structures')}" \\
    --structures-map-csv "{self.structures_map}" \\
    --trajectories-dir "{self.stream_folder('trajectories')}" \\
    --trajectories-map-csv "{self.trajectories_map}" \\
    --summary-csv "{self.summary_csv}" \\
    --physicality-csv "{self.physicality_csv}" \\
    --local-missing-csv "{self.local_missing_csv}" \\
    --container-prefix "{self.container_prefix()}"{msas_flag}
"""
        script += self.generate_missing_propagation(
            *self._input_sources(), local_missing=self.local_missing_csv, missing_csv=self.missing_csv)
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        run_ids, _ = predict_output_ids_with_provenance(**self._axis_kwargs())
        structure_ids = generate_multiplied_ids_pattern(run_ids, f"<1..{self.num_samples}>")
        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=[self.stream_path("structures", "<id>.pdb")],
            map_table=self.structures_map,
            format="pdb",
        )
        trajectories = DataStream(
            name="trajectories",
            ids=run_ids,
            files=[self.stream_path("trajectories", "<id>.xtc")],
            map_table=self.trajectories_map,
            format="xtc",
        )
        tables = {
            "summary": TableInfo(
                name="summary",
                path=self.summary_csv,
                columns=["id", "sequences.id", "ligands.id", "n_samples", "n_residues"],
                description="BioEmu2 per-run ensemble summary",
            ),
            "physicality": TableInfo(
                name="physicality",
                path=self.physicality_csv,
                columns=["id", "run.id", "frame", "status", "reason"],
                description="Per-conformer physicality report from bioemu2",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "structures": structures,
            "trajectories": trajectories,
            "tables": tables,
            "output_folder": self.output_folder,
        }
