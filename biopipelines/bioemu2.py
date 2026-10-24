# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""BioEmu2 tool: ligand-capable equilibrium-ensemble emulation.

bioemu2 is the ligand-capable successor of bioemu. Each run samples
``num_samples`` conformers of a system described by a system YAML with
protein and (optionally) ligand entities, plus optional contact / covalent
constraints. The wrapper writes one YAML per input entity set, calls
``bioemu2 sample-equilibrium``, splits the sampled trajectory into
per-conformer PDBs and keeps the compact trajectory + physicality report.

Reference: Lewis et al. (2024) BioEmu, and the ligand-aware extension at
https://github.com/microsoft/bioemu2
"""

import os
from typing import Any, Dict, List, Optional, Union

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


SAMPLERS = ("edm", "edm_gd", "edm_fk", "edm_fkc")
PHYSICALITY_MODES = ("high_quality", "medium_quality", "low_quality", "no_filter")


class BioEmu2(BaseConfig):
    """
    BioEmu2: ligand-capable equilibrium-ensemble emulation.

    Each (sequence, ligand set, constraints) combination is one BioEmu2 run.
    Ligands are passed by SMILES and placed by bioemu2 itself from scratch —
    no 3D template is seeded. A precomputed MSA may be supplied per sequence
    to skip bioemu2's own ColabFold cold-start.

    Inputs:
        sequences: amino-acid sequences (sequences stream or StandardizedOutput).
                   Each id's sequence is read from the stream's ``sequence``
                   column.
        ligands:   optional ligand stream. For each (sequence, ligand) pair the
                   wrapper builds one system with the protein as chain A and
                   each ligand entity after it. Pass ``None`` for a plain apo
                   ensemble. SMILES are read from the ligand stream.
        num_samples: conformers per run (default 100).
        batch_size: samples denoised together per forward pass (default 16).
                   Caps GPU memory at long sequence × large batch.
        seed:      optional RNG seed. If ``None``, each run is reseeded.
        sampler:   one of ``edm`` (default), ``edm_gd``, ``edm_fk``, ``edm_fkc``.
                   The three steered variants need a ``potentials`` YAML.
        contacts:  optional list of contact-constraint dicts (bioemu2 YAML
                   schema); each dict is forwarded verbatim into the system
                   YAML's ``contacts`` list.
        covalent_linkage: optional covalent-bond dict (bioemu2 YAML schema);
                   forwarded verbatim into the system YAML.
        precomputed_msa: optional path (or stream) to a precomputed .a3m file.
                   When provided, the wrapper pre-populates bioemu2's
                   per-sequence cache so no ColabFold query is made.
        physicality_filter: ``high_quality`` (default), ``medium_quality``,
                   ``low_quality`` or ``no_filter``. Checked frames stay in the
                   trajectory; the mode only affects the physicality report.

    Install:
        BioEmu2.install()

    Outputs:
        Streams:
            structures:  per-conformer PDB, ids ``<run>_<1..num_samples>``.
            trajectories: one ``<run>.xtc`` + ``<run>_topology.cif`` per run.
        Tables:
            summary:  id | sequence.id | ligand.id | n_samples | n_residues
            physicality: id | frame | invalid_bond | clash | chirality | chain_break
            missing:  id | removed_by | kind | cause
    """

    TOOL_NAME = "BioEmu2"
    TOOL_VERSION = "1.0"
    ENV_NAME = "bioemu2"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Install the bioemu2 env. bioemu2 is a pure-pip package; no conda-only
        dependencies, no AF2 params to wire (bioemu2 uses ColabFold MMseqs2
        with its own on-disk cache under ``~/.cache/bioemu2``)."""
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        import_check = (
            f'{cls._env_run(env, env_manager)}'
            f'python -c "import bioemu2" >/dev/null 2>&1'
        )
        cli_check = (
            f'{cls._env_run(env, env_manager)}'
            f'bioemu2 --help >/dev/null 2>&1'
        )

        skip = "" if force_reinstall else f"""# Check if already installed
if {import_check} && {cli_check}; then
    echo "BioEmu2 already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        env_install_block = f"""if ! {import_check}; then
{env_block}
else
    echo "{env} environment already present, skipping creation."
fi"""
        return f"""echo "=== Installing BioEmu2 ==="
{skip}{remove_block}
{env_install_block}

# Verify installation
if {import_check} && {cli_check}; then
    touch "$INSTALL_SUCCESS"
    echo "=== BioEmu2 installation complete ==="
else
    echo "ERROR: BioEmu2 verification failed (cannot import bioemu2 or run its CLI)"
    exit 1
fi
"""

    sequences_json = Path(lambda self: self.configuration_path("sequences.json"))
    ligands_json = Path(lambda self: self.configuration_path("ligands.json"))
    runs_json = Path(lambda self: self.configuration_path("runs.json"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    trajectories_map = Path(lambda self: self.stream_map_path("trajectories"))
    summary_csv = Path(lambda self: self.table_path("summary"))
    physicality_csv = Path(lambda self: self.table_path("physicality"))
    missing_csv = Path(lambda self: self.table_path("missing"))
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
                 precomputed_msa: Optional[Union[DataStream, StandardizedOutput]] = None,
                 physicality_filter: str = "high_quality",
                 **kwargs):
        self.sequences = sequences
        if isinstance(sequences, StandardizedOutput):
            self.sequences_stream: DataStream = sequences.streams.sequences
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
        else:
            raise ValueError(
                f"sequences must be DataStream or StandardizedOutput, got {type(sequences).__name__}"
            )

        self.ligands = ligands
        if ligands is None:
            self.ligands_stream: Optional[DataStream] = None
        elif isinstance(ligands, StandardizedOutput):
            # Prefer the compounds stream (SMILES-bearing) when present,
            # otherwise fall back to a stream named "ligands".
            try:
                self.ligands_stream = ligands.streams.compounds
            except AttributeError:
                self.ligands_stream = ligands.streams.ligands
        elif isinstance(ligands, DataStream):
            self.ligands_stream = ligands
        else:
            raise ValueError(
                f"ligands must be DataStream, StandardizedOutput or None, got {type(ligands).__name__}"
            )

        self.precomputed_msa = precomputed_msa
        if precomputed_msa is None:
            self.precomputed_msa_stream: Optional[DataStream] = None
        elif isinstance(precomputed_msa, StandardizedOutput):
            self.precomputed_msa_stream = precomputed_msa.streams.msas
        elif isinstance(precomputed_msa, DataStream):
            self.precomputed_msa_stream = precomputed_msa
        else:
            raise ValueError(
                f"precomputed_msa must be DataStream, StandardizedOutput or None, "
                f"got {type(precomputed_msa).__name__}"
            )

        self.num_samples = int(num_samples)
        self.batch_size = int(batch_size)
        self.seed = int(seed) if seed is not None else None
        self.sampler = str(sampler)
        self.contacts = list(contacts) if contacts else []
        self.covalent_linkage = dict(covalent_linkage) if covalent_linkage else None
        self.physicality_filter = str(physicality_filter)
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")
        if self.num_samples < 1:
            raise ValueError("num_samples must be >= 1")
        if self.batch_size < 1:
            raise ValueError("batch_size must be >= 1")
        if self.sampler not in SAMPLERS:
            raise ValueError(f"sampler must be one of {SAMPLERS}, got {self.sampler!r}")
        if self.physicality_filter not in PHYSICALITY_MODES:
            raise ValueError(
                f"physicality_filter must be one of {PHYSICALITY_MODES}, "
                f"got {self.physicality_filter!r}"
            )
        if self.covalent_linkage is not None:
            required = {"chain", "position", "protein_atom", "ligand_atom"}
            missing = required - set(self.covalent_linkage)
            if missing:
                raise ValueError(
                    f"covalent_linkage missing required keys: {sorted(missing)}"
                )
            if self.ligands_stream is None or len(self.ligands_stream) == 0:
                raise ValueError("covalent_linkage requires a ligands stream")
        for i, c in enumerate(self.contacts):
            if not isinstance(c, dict) or "token1" not in c or "token2" not in c:
                raise ValueError(
                    f"contacts[{i}] must be a dict with token1 and token2 keys"
                )

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"SEQUENCES: {len(self.sequences_stream)}")
        if self.ligands_stream is not None:
            lines.append(f"LIGANDS: {len(self.ligands_stream)}")
            lines.append(f"RUNS (seq × lig): {self._n_runs()}")
        else:
            lines.append("LIGANDS: none (apo)")
            lines.append(f"RUNS: {len(self.sequences_stream)}")
        lines.append(f"NUM SAMPLES PER RUN: {self.num_samples}")
        lines.append(f"BATCH SIZE: {self.batch_size}")
        lines.append(f"SAMPLER: {self.sampler}")
        if self.contacts:
            lines.append(f"CONTACTS: {len(self.contacts)}")
        if self.covalent_linkage is not None:
            lines.append(
                "COVALENT LINKAGE: "
                f"chain {self.covalent_linkage['chain']} pos "
                f"{self.covalent_linkage['position']}.{self.covalent_linkage['protein_atom']} — "
                f"ligand {self.covalent_linkage['ligand_atom']}"
            )
        return lines

    def _n_runs(self) -> int:
        if self.ligands_stream is None:
            return len(self.sequences_stream)
        return len(self.sequences_stream) * len(self.ligands_stream)

    def _run_ids(self) -> List[str]:
        if self.ligands_stream is None:
            return list(self.sequences_stream.ids)
        return [f"{s}+{l}" for s in self.sequences_stream.ids
                for l in self.ligands_stream.ids]

    def generate_script(self, script_path: str) -> str:
        self.sequences_stream.save_json(self.sequences_json)
        ligands_arg = ""
        if self.ligands_stream is not None:
            self.ligands_stream.save_json(self.ligands_json)
            ligands_arg = f' --ligands-json "{self.ligands_json}"'
        msa_arg = ""
        if self.precomputed_msa_stream is not None:
            self.precomputed_msa_stream.save_json(self.runs_json)
            msa_arg = f' --precomputed-msa-json "{self.runs_json}"'
        seed_arg = f" --seed {self.seed}" if self.seed is not None else ""
        import json as _json
        contacts_arg = ""
        if self.contacts:
            contacts_arg = f" --contacts-json '{_json.dumps(self.contacts)}'"
        covalent_arg = ""
        if self.covalent_linkage is not None:
            covalent_arg = (
                f" --covalent-linkage-json '{_json.dumps(self.covalent_linkage)}'"
            )
        upstream_missing_path = self._get_upstream_missing_table_path(self.sequences)
        upstream_missing_flag = (
            f' --upstream-missing "{upstream_missing_path}"' if upstream_missing_path else ""
        )

        script = "#!/bin/bash\n"
        script += "# BioEmu2 ligand-capable equilibrium-ensemble sampling\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += self.warn_container_unsupported()
        script += f"""echo "BioEmu2: sampling {self._n_runs()} run(s), {self.num_samples} conformers each"
python "{self.helper_py}" \\
    --sequences-json "{self.sequences_json}"{ligands_arg} \\
    --num-samples {self.num_samples} \\
    --batch-size {self.batch_size} \\
    --sampler "{self.sampler}" \\
    --physicality-filter "{self.physicality_filter}" \\
    --structures-dir "{self.stream_folder('structures')}" \\
    --structures-map "{self.structures_map}" \\
    --trajectories-dir "{self.stream_folder('trajectories')}" \\
    --trajectories-map "{self.trajectories_map}" \\
    --scratch-dir "{self.extras_path()}" \\
    --summary-csv "{self.summary_csv}" \\
    --physicality-csv "{self.physicality_csv}" \\
    --missing-csv "{self.missing_csv}"{seed_arg}{contacts_arg}{covalent_arg}{msa_arg}{upstream_missing_flag}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        run_ids = self._run_ids()
        suffix_pattern = f"<1..{self.num_samples}>"
        structure_ids = generate_multiplied_ids_pattern(
            run_ids, suffix_pattern, input_stream_name="runs",
        )
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
                columns=["id", "sequence.id", "ligand.id", "n_samples", "n_residues"],
                description="BioEmu2 per-run ensemble summary",
            ),
            "physicality": TableInfo(
                name="physicality",
                path=self.physicality_csv,
                columns=["id", "frame", "invalid_bond", "clash",
                         "chirality", "chain_break"],
                description="Per-frame physicality checks reported by bioemu2",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="IDs removed (upstream or local failure) with removal reason",
            ),
        }
        return {
            "structures": structures,
            "trajectories": trajectories,
            "tables": tables,
            "output_folder": self.output_folder,
        }
