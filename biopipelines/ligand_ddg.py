# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""LigandDdG tool: protein-ligand binding energy with PyRosetta.

Two scoring schemes, both ported from the Chaperone blocks of RosettaCommons: "hard" (`ligand_ddg`, a minimization protocol read through Rosetta's Ddg filter with no repacking) and "soft" (`ligand_ddg_soft`, interface repacked with the hard beta_nov16 score and scored with the soft-repulsive beta_nov16_soft, trimmed mean over replicates). Values are Rosetta energy units, not measured binding free energies, and the two schemes are not numerically comparable.

Rosetta needs a .params file for a non-standard ligand. It is generated from the ligand in each structure with Rosetta's molfile_to_params.py (fetched at install time, since the PyRosetta wheel does not ship it), using the SMILES of the compounds stream for bond orders; a pre-made params file can be supplied instead.

PyRosetta is free for academic and non-profit use and needs a paid commercial license from the University of Washington; installing it accepts the PyRosetta license (https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.PyRosetta.md).

Reference:
    Mirrors the `ligand_ddg` and `ligand_ddg_soft` blocks of Chaperone, https://github.com/RosettaCommons/chaperone (MIT).
    PyRosetta: Chaudhury, Lyskov & Gray (2010) Bioinformatics 26, 689.
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .input_standardization import resolve_basic_input
    from .ligand import Ligand
    from .rosetta_scoring import pyrosetta_install_script, validate_chain, validate_positive_int
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from input_standardization import resolve_basic_input
    from ligand import Ligand
    from rosetta_scoring import pyrosetta_install_script, validate_chain, validate_positive_int

SCORINGS = ("soft", "hard")
DEFAULT_SOFT_REPLICATES = 5


class LigandDdG(BaseConfig):
    """
    LigandDdG: protein-ligand binding energy (PyRosetta).

    ddG is the Rosetta score of the complex minus the score with the protein chain translated 10 000 A away from the ligand; more negative means stronger predicted binding. Both schemes are Rosetta energy units, not kcal/mol and not measured affinities, and a "hard" value must not be compared with a "soft" one.

    PyRosetta is free for academic use and needs a paid commercial license; installing it accepts the PyRosetta license.

    Inputs:
        structures: PDB complexes (StandardizedOutput or DataStream); mmCIF is not read. Only `binder_chain` and the ligand are kept; other chains are dropped. Everything else on the binder chain stays as it is, hetero residues and waters included, except hydrogens, which are discarded for Rosetta to rebuild.
        ligand: compounds stream naming the ligand (`Ligand(smiles=..., codes="LIG")` or any compounds-producing tool); a bare string is shorthand for `Ligand(codes=...)`. The residue code is read from the stream's `code` column and must be the residue name of the ligand in the structures (at most 3 characters in a PDB). Unless `params_file` is given, the stream must carry a single SMILES, which supplies the bond orders and charges; a code-only ligand cannot be parameterized.
        scoring: "soft" (default) or "hard".
            soft: the interface within 10 A of the ligand is repacked with beta_nov16, and the complex is scored with beta_nov16_soft against the separated state; repeated `replicates` times on the same pose, then the best and worst values are dropped when there are at least three and the rest averaged. Forgiving of the small clashes of predicted structures.
            hard: backbone-constrained soft-then-hard minimization of the interface followed by Rosetta's Ddg filter (beta_nov16, no repacking) across the ligand's jump. A single run; `replicates` does not apply.
        binder_chain: the protein chain, translated away to give the unbound state (default "A").
        ligand_chain: the chain carrying the ligand (default "B"). Exactly one residue named like the ligand code must sit on it, and it must be a chain of its own, as predictors write it: the ligand of a crystal file usually shares the protein chain and has to be moved to its own first. A ligand with no atom within 6 A of the binder chain (the ligand of another copy in a crystal, say) is a failure for that structure, since no pair of atoms is in range and the ddG would read as exactly zero.
        replicates: pack + score cycles for scoring="soft" (default 5). Not accepted with scoring="hard".
        ph: pH for Dimorphite-DL protonation of a ligand the structure carries no hydrogens for (default 7.4). A ligand that already has hydrogens is used as it is. Ignored with `params_file`.
        params_file: path to a pre-made Rosetta .params file. The ligand is then taken from the structure unchanged, so its atom names and hydrogens must already match the params; no SMILES is needed and no params are generated.

    Outputs:
        Streams: (none)
        Tables:
            ligand_ddg: scoring="hard": id | ddg_norepack
                        scoring="soft": id | ddg_soft | ddg_soft_std | n_replicates_used
                        ddg_soft_std is the population standard deviation of the retained replicates and n_replicates_used their count. A structure that fails keeps its row with NaN values.
            missing: id | removed_by | kind | cause
    """

    TOOL_NAME = "LigandDdG"
    TOOL_VERSION = "1.1"
    ENV_NAME = "PyRosetta"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, env_name=None, **kwargs):
        return pyrosetta_install_script(cls, folders, env_manager, force_reinstall, env_name)

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    ligand_json = Path(lambda self: self.configuration_path("ligand.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    ddg_csv = Path(lambda self: self.table_path("ligand_ddg"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_ligand_ddg.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: Union[str, DataStream, StandardizedOutput],
                 scoring: str = "soft",
                 binder_chain: str = "A",
                 ligand_chain: str = "B",
                 replicates: Optional[int] = None,
                 ph: float = 7.4,
                 params_file: Optional[str] = None,
                 **kwargs):
        self.structures = structures
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        self.ligand_input = ligand
        self.ligand_stream: DataStream = resolve_basic_input(
            ligand, Ligand, "compounds", "codes", allow_none=False)
        self.scoring = scoring
        self.binder_chain = binder_chain
        self.ligand_chain = ligand_chain
        self.replicates = replicates
        self.ph = ph
        self.params_file = params_file
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.structures_stream.has_only_formats("pdb"):
            raise ValueError(
                f"LigandDdG reads PDB files; got format {self.structures_stream.format!r}")
        if not self.ligand_stream or len(self.ligand_stream) == 0:
            raise ValueError(
                "ligand (a compounds stream, e.g. Ligand(smiles=..., codes=...)) is required and must not be empty")
        if self.scoring not in SCORINGS:
            raise ValueError(f"scoring must be one of {SCORINGS}, got {self.scoring!r}")
        validate_chain("binder_chain", self.binder_chain)
        validate_chain("ligand_chain", self.ligand_chain)
        if self.binder_chain == self.ligand_chain:
            raise ValueError("binder_chain and ligand_chain must name different chains")
        if self.scoring == "hard" and self.replicates is not None:
            raise ValueError("replicates applies to scoring='soft' only; scoring='hard' runs once")
        if self.replicates is not None:
            validate_positive_int("replicates", self.replicates)
        if isinstance(self.ph, bool) or not isinstance(self.ph, (int, float)) or not 0 <= self.ph <= 14:
            raise ValueError(f"ph must be a number between 0 and 14, got {self.ph!r}")
        if self.params_file is not None and (not isinstance(self.params_file, str) or not self.params_file.strip()):
            raise ValueError("params_file must be a path to a Rosetta .params file")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"SCORING: {self.scoring}")
        lines.append(f"CHAINS: protein {self.binder_chain}, ligand {self.ligand_chain}")
        if self.scoring == "soft":
            lines.append(f"REPLICATES: {self._resolved_replicates()}")
        if self.params_file:
            lines.append(f"PARAMS: {self.params_file}")
        else:
            lines.append(f"LIGAND PARAMS: generated from the compounds stream SMILES at pH {self.ph}")
        return lines

    def _resolved_replicates(self) -> int:
        return DEFAULT_SOFT_REPLICATES if self.replicates is None else self.replicates

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self.ligand_stream.save_json(self.ligand_json)
        self._write_config_yaml()
        script = "#!/bin/bash\n"
        script += "# LigandDdG script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += self.warn_container_unsupported()
        script += f"""echo "Running LigandDdG ({self.scoring}) on {len(self.structures_stream)} structure(s)"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-json "{self.structures_json}" \\
    --ligand-json "{self.ligand_json}" \\
    --work-dir "{self.execution_path()}" \\
    --table-csv "{self.ddg_csv}" \\
    --local-missing-csv "{self.local_missing_csv}"
"""
        script += self.generate_missing_propagation(
            self.structures, self.ligand_input,
            local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        import yaml
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump({
                "scoring": self.scoring,
                "binder_chain": self.binder_chain,
                "ligand_chain": self.ligand_chain,
                "replicates": self._resolved_replicates() if self.scoring == "soft" else None,
                "ph": float(self.ph),
                "params_file": self.params_file or "",
            }, f, sort_keys=False)

    def get_output_files(self) -> Dict[str, Any]:
        if self.scoring == "hard":
            columns = ["id", "ddg_norepack"]
        else:
            columns = ["id", "ddg_soft", "ddg_soft_std", "n_replicates_used"]
        tables = {
            "ligand_ddg": TableInfo(
                name="ligand_ddg",
                path=self.ddg_csv,
                columns=columns,
                description=f"Rosetta protein-ligand ddG ({self.scoring} scoring)",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {"tables": tables, "output_folder": self.output_folder}
