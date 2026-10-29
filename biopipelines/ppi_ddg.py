# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""PPIDdG tool: protein-protein binding energy with PyRosetta.

Relaxes each complex with FastRelax (optionally constrained to the input coordinates) and reports ddG = score(complex) - score(binder chain translated 10 000 A away), with the beta_nov16 score function and no repacking of the separated state. The numbers are Rosetta energy units, not measured binding free energies.

PyRosetta is free for academic and non-profit use and needs a paid commercial license from the University of Washington; installing it accepts the PyRosetta license (https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.PyRosetta.md).

Reference:
    Mirrors the `ppi_ddg` block of Chaperone, https://github.com/RosettaCommons/chaperone (MIT).
    PyRosetta: Chaudhury, Lyskov & Gray (2010) Bioinformatics 26, 689.
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .rosetta_scoring import DEFAULT_SEED, pyrosetta_install_script, validate_chain, validate_positive_int
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from rosetta_scoring import DEFAULT_SEED, pyrosetta_install_script, validate_chain, validate_positive_int


class PPIDdG(BaseConfig):
    """
    PPIDdG: protein-protein binding energy by FastRelax + ddG (PyRosetta).

    For each input complex, `replicates` independent FastRelax runs are performed on a copy of the structure. After each, ddG = score(relaxed complex) - score(`binder_chain` translated 10 000 A away) is computed with beta_nov16 and no repacking. The reported value is the mean over replicates; more negative means stronger predicted binding. Rosetta energy units, not kcal/mol and not a measured affinity; compare only across runs with identical chains, relaxation and score function.

    PyRosetta is free for academic use and needs a paid commercial license; installing it accepts the Rosetta license. Intended for complexes of the 20 standard amino acids: a ligand or other hetero residue makes Rosetta fail on that structure.

    Inputs:
        structures: PDB complexes (StandardizedOutput or DataStream); mmCIF is not read.
        binder_chain: chain translated away to give the unbound state (default "A").
        target_chain: the partner that stays in place (default "B"). It must be present in the structure; as in Chaperone's block, every chain other than `binder_chain` stays in place, so extra chains count as part of the target.
        relax_rounds: FastRelax repeats per replicate (default 1; roughly 1-3 min per structure on a CPU core, 3-5 gives lower, more converged energies).
        replicates: independent relax + ddG cycles; mean and standard deviation are reported (default 1, which gives a standard deviation of 0). Each replicate adds one full relax run.
        constrain_to_start_coords: harmonic coordinate constraints to the input during FastRelax, keeping the backbone near the supplied geometry (default True).
        constraint_weight: weight of those constraints; higher stays closer to the input (default 1.0).
        seed: base random seed (default 42). Every replicate of every structure is seeded from it with a distinct, fixed value, so a run repeats exactly; PyRosetta otherwise seeds itself from the system and a rerun gives other numbers. A different seed gives a different relaxation, with ddG differences of a few energy units at one replicate.

    Outputs:
        Streams: (none)
        Tables:
            ppi_ddg: id | ddg | ddg_std | n_replicates | relax_score
                     ddg is the mean over replicates, ddg_std the population standard deviation, relax_score the total beta_nov16 score of the relaxed complex from the last replicate. A structure that fails keeps its row with NaN values.
            missing: id | removed_by | kind | cause
    """

    TOOL_NAME = "PPIDdG"
    TOOL_VERSION = "1.1"
    ENV_NAME = "PyRosetta"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, env_name=None, **kwargs):
        return pyrosetta_install_script(cls, folders, env_manager, force_reinstall, env_name)

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    ddg_csv = Path(lambda self: self.table_path("ppi_ddg"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_ppi_ddg.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 binder_chain: str = "A",
                 target_chain: str = "B",
                 relax_rounds: int = 1,
                 replicates: int = 1,
                 constrain_to_start_coords: bool = True,
                 constraint_weight: float = 1.0,
                 seed: int = DEFAULT_SEED,
                 **kwargs):
        self.structures = structures
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        self.binder_chain = binder_chain
        self.target_chain = target_chain
        self.relax_rounds = relax_rounds
        self.replicates = replicates
        self.constrain_to_start_coords = constrain_to_start_coords
        self.constraint_weight = constraint_weight
        self.seed = seed
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.structures_stream.has_only_formats("pdb"):
            raise ValueError(
                f"PPIDdG reads PDB files; got format {self.structures_stream.format!r}")
        validate_chain("binder_chain", self.binder_chain)
        validate_chain("target_chain", self.target_chain)
        if self.binder_chain == self.target_chain:
            raise ValueError("binder_chain and target_chain must name different chains")
        validate_positive_int("relax_rounds", self.relax_rounds)
        validate_positive_int("replicates", self.replicates)
        if not isinstance(self.constrain_to_start_coords, bool):
            raise ValueError(
                f"constrain_to_start_coords must be a bool, got {self.constrain_to_start_coords!r}")
        if isinstance(self.constraint_weight, bool) or not isinstance(self.constraint_weight, (int, float)) \
                or self.constraint_weight < 0:
            raise ValueError(f"constraint_weight must be a number >= 0, got {self.constraint_weight!r}")
        validate_positive_int("seed", self.seed)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"CHAINS: binder {self.binder_chain}, target {self.target_chain}")
        lines.append(f"RELAX: {self.relax_rounds} round(s) x {self.replicates} replicate(s)")
        if self.constrain_to_start_coords:
            lines.append(f"COORDINATE CONSTRAINTS: weight {self.constraint_weight}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self._write_config_yaml()
        script = "#!/bin/bash\n"
        script += "# PPIDdG script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += self.warn_container_unsupported()
        script += f"""echo "Running PPIDdG on {len(self.structures_stream)} structure(s)"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-json "{self.structures_json}" \\
    --table-csv "{self.ddg_csv}" \\
    --local-missing-csv "{self.local_missing_csv}"
"""
        script += self.generate_missing_propagation(
            self.structures, local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        import yaml
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump({
                "binder_chain": self.binder_chain,
                "target_chain": self.target_chain,
                "relax_rounds": self.relax_rounds,
                "replicates": self.replicates,
                "constrain_to_start_coords": self.constrain_to_start_coords,
                "constraint_weight": float(self.constraint_weight),
                "seed": self.seed,
            }, f, sort_keys=False)

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "ppi_ddg": TableInfo(
                name="ppi_ddg",
                path=self.ddg_csv,
                columns=["id", "ddg", "ddg_std", "n_replicates", "relax_score"],
                description="Rosetta ddG of the complex after FastRelax, mean and spread over replicates",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {"tables": tables, "output_folder": self.output_folder}
