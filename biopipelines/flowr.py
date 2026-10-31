# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
FLOWR generates 3D ligands inside a protein pocket by flow matching.

Wraps `flowr.gen.generate_from_pdb` of https://github.com/jule-c/flowr (Cremer et al., Nature Computational Science 2026; arXiv:2504.10564).
"""

import os
from typing import Dict, List, Any, Optional, Union

# Conditioning tasks the released checkpoint was trained for; pipe_flowr.py maps each to one upstream inpainting flag.
CONDITIONING = ("interaction", "functional_group", "scaffold", "linker", "substructure")

# The no-explicit-hydrogen checkpoint of Zenodo record 15737419; the record is a fixed version, and the md5 is the one Zenodo lists.
CHECKPOINT_NAME = "flowr_noHs.ckpt"
CHECKPOINT_URL = "https://zenodo.org/api/records/15737419/files/flowr_noHs.ckpt/content"
CHECKPOINT_MD5 = "e0d49765629c4b9a6847c4dcdbf678ad"

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import generate_multiplied_ids_pattern
    from .input_standardization import resolve_basic_input
    from .ligand import Ligand
    from .sele_utils import sele_to_chains_and_list
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_multiplied_ids_pattern
    from input_standardization import resolve_basic_input
    from ligand import Ligand
    from sele_utils import sele_to_chains_and_list


class FLOWR(BaseConfig):
    """
    FLOWR: pocket-conditioned generation of 3D ligands by flow matching.

    Each input structure is one pocket. FLOWR samples `num_molecules` new ligands in it and the step returns them as a `compounds` stream (SMILES) and as a `structures` stream (one SDF pose per molecule, in the coordinate frame of the input structure). The molecules are not scored by FLOWR: upstream writes no per-molecule score, so `pockets` only counts what survived its filters.

    The pocket is given in one of two ways, and exactly one is required:
        ligand:           the residue of a ligand bound in each structure, as `Ligand("ZIT")` or the bare code `"ZIT"` (the first copy when a structure holds several); its `smiles`, when the stream has one, is the bond-order template of the extracted reference. FLOWR cuts the pocket from the protein residues within `pocket_cutoff` of that ligand, and the reference sets the size of the molecules unless `sample_sizes` is on. Required for `conditioning`.
        pocket_residues:  a chain-aware selection such as `"A10-25+A41+B7"` (a bare chain letter takes the whole chain). Only those residues are the pocket, and `num_heavy_atoms` fixes the molecule size. No reference ligand, so no `conditioning`.

    Inputs:
        structures:         protein or complex structures (PDB).
        ligand:             see above.
        pocket_residues:    see above.
        num_molecules:      molecules requested per pocket (default 100). Upstream keeps sampling until it has this many valid molecules or reaches its round limit, then drops molecules above 0.9 Tanimoto similarity (ECFP4) to an earlier one, so fewer can come back.
        integration_steps:  Euler steps of the flow (default 100, upstream's default).
        conditioning:       None (de novo, default) or one of "interaction", "functional_group", "scaffold", "linker", "substructure". Atoms of the reference ligand stay fixed at the reference coordinates and the rest of the molecule is generated: "interaction" fixes the atoms that interact with the pocket (ProLIF, so the complex should carry hydrogens, see `protonate_pocket`), "functional_group" the functional-group atoms, "scaffold" the Murcko scaffold, "linker" everything except the linker, which is regenerated, and "substructure" the atoms matching `substructure`. When the reference has no atom of the requested kind, nothing is fixed. Only one task runs at a time.
        substructure:       SMILES of the fragment to keep fixed, required with `conditioning="substructure"`. It is matched against the reference ligand as upstream does; a query that does not match makes FLOWR fail for that pocket.
        sample_sizes:       draw each molecule's atom count uniformly within 10 % of the reference size (or of `num_heavy_atoms`), clipped to 5-80 heavy atoms. Default None: on for de novo generation from a reference ligand, as in upstream's example script, and off otherwise. Upstream supports `conditioning` only at the reference size, so sample_sizes=True with it is refused.
        num_heavy_atoms:    heavy atoms of the generated molecules; required with `pocket_residues` and not allowed with `ligand`.
        pocket_cutoff:      distance (A) from the reference ligand within which protein residues enter the pocket (default 6.0, upstream's default).
        protonate_pocket:   add hydrogens to the protein with hydride before the pocket is built (default False). Needed for "interaction" when the structure has none. Worked with "interaction" on streptavidin (1STP); on factor Xa (2P16) upstream's ProLIF step failed on the hydrogenated pocket (RDKit: explicit valence of a hydrogen is 4).
        protonate_ligands:  add hydrogens to the generated molecules and relax them in the pocket (upstream's `--protonate_generated_ligands`, default False). The poses are then the protonated ones; the SMILES never carry explicit hydrogens. Failed on both structures tried (2P16, 1STP): upstream's pocket reading raises an RDKit KekulizeException inside the protonation, swallows it and returns no molecule, so the pocket fails.
        batch_cost:         upstream's batch size knob (default 100, the value of its example script). Upstream recommends a GPU with at least 40 GB; lower it on a smaller card.
        seed:               random seed (default 42, upstream's default).

    Resources: GPU. The environment ships torch 2.5.1 with CUDA 11.8, as upstream's does.

    Not implemented: the checkpoint trained with explicit hydrogens, flexible (apo) pockets and `pocket_flex`, interaction recovery and the other evaluation scripts, and the sampler options `coord_noise_std`, `categorical_strategy`, `ode_sampling_strategy` and `corrector_iters`.

    Outputs:
        Streams:
            compounds:  one row per molecule, ids `<structure id>_<k>` with k the 1-based record of upstream's SDF (value-based csv).
            structures: one `<id>.sdf` pose per molecule, same ids.
        Tables:
            compounds:  id | structures.id | sample | format | code | ccd | smiles
                        `format` is "smiles", `code` is "LIG" and `ccd` is empty. `smiles` is the isomeric canonical SMILES RDKit reads off the pose.
            pockets:    id | pocket | n_requested | n_generated
                        one row per input structure. `pocket` is "ligand:<code>" or "residues:<selection>"; n_generated is NaN when FLOWR failed on the structure.
            missing:    id | removed_by | kind | cause
                        a pocket FLOWR failed on is a "failure" row under the structure id; a molecule id below `num_molecules` that upstream's filters did not return is a "filter" row.

    Usage::

        with Pipeline(project="Examples", job="FLOWR-demo"):
            Resources(gpu="A100", memory="32GB", time="2:00:00")
            FLOWR.install()
            complex_ = PDB("2p16", convert="pdb")
            gen = FLOWR(structures=complex_, ligand=Ligand("GG2", smiles="..."), num_molecules=50)
    """

    TOOL_NAME = "FLOWR"
    TOOL_VERSION = "1.3"
    ENV_NAME = "FLOWR"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("FLOWR", "")
        parent_dir = os.path.dirname(repo_dir)
        checkpoint = os.path.join(repo_dir, "checkpoints", CHECKPOINT_NAME)
        import_check = f'cd "{repo_dir}" && PYTHONPATH="{repo_dir}" {cls._env_run(env, env_manager)}python -c "import pandas, yaml, rdkit, flowr.gen.generate_from_pdb"'

        env_check = cls._env_exists_check(env, env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -d "{repo_dir}/flowr" ] && [ -s "{checkpoint}" ] && {env_check} \\
   && ( {import_check} ) >/dev/null 2>&1; then
    echo "FLOWR already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        return f"""echo "=== Installing FLOWR ==="
{skip}mkdir -p "{parent_dir}"
cd "{parent_dir}"
if [ ! -d "{repo_dir}/flowr" ]; then
    git clone https://github.com/jule-c/flowr.git "{repo_dir}"
fi

{remove_block}
# Create the env only when absent: a missing checkpoint must not re-create the CUDA env
if ! {env_check}; then
{env_block}
fi

mkdir -p "{repo_dir}/checkpoints"
if [ ! -s "{checkpoint}" ]; then
    curl -fsSL -o "{checkpoint}.part" "{CHECKPOINT_URL}" || {{ echo "ERROR: download of {CHECKPOINT_URL} failed"; exit 1; }}
    if [ "$(md5sum "{checkpoint}.part" | cut -d' ' -f1)" != "{CHECKPOINT_MD5}" ]; then
        echo "ERROR: the downloaded checkpoint does not match the md5 Zenodo lists for it"
        rm -f "{checkpoint}.part"
        exit 1
    fi
    mv "{checkpoint}.part" "{checkpoint}"
fi

# Importing upstream's entry module pulls in every dependency the run needs; the rest are what the pipe script imports.
if [ -s "{checkpoint}" ] && ( {import_check} ); then
    touch "$INSTALL_SUCCESS"
    echo "=== FLOWR installation complete ==="
else
    echo "ERROR: FLOWR verification failed (checkpoint missing, or pandas, yaml, rdkit or flowr.gen.generate_from_pdb does not import)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    ligand_json = Path(lambda self: self.configuration_path("ligand.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    runs_folder = Path(lambda self: self.execution_path("runs"))
    compounds_csv = Path(lambda self: self.stream_path("compounds", "compounds.csv"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    pockets_csv = Path(lambda self: self.table_path("pockets"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    checkpoint = Path(lambda self: os.path.join(self.folders["FLOWR"], "checkpoints", CHECKPOINT_NAME))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_flowr.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: Union[str, DataStream, StandardizedOutput, None] = None,
                 pocket_residues: Optional[str] = None,
                 num_molecules: int = 100,
                 integration_steps: int = 100,
                 conditioning: Optional[str] = None,
                 substructure: Optional[str] = None,
                 sample_sizes: Optional[bool] = None,
                 num_heavy_atoms: Optional[int] = None,
                 pocket_cutoff: float = 6.0,
                 protonate_pocket: bool = False,
                 protonate_ligands: bool = False,
                 batch_cost: int = 100,
                 seed: int = 42,
                 **kwargs):
        self.structures = structures
        self.ligand = ligand
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        self.ligand_stream: Optional[DataStream] = resolve_basic_input(
            ligand, Ligand, "compounds", "codes", allow_none=True)
        self.pocket_residues = pocket_residues
        self.num_molecules = num_molecules
        self.integration_steps = integration_steps
        self.conditioning = conditioning
        self.substructure = substructure
        self.sample_sizes_requested = sample_sizes
        self.sample_sizes = (conditioning is None and self.ligand_stream is not None) if sample_sizes is None else sample_sizes
        self.num_heavy_atoms = num_heavy_atoms
        self.pocket_cutoff = pocket_cutoff
        self.protonate_pocket = protonate_pocket
        self.protonate_ligands = protonate_ligands
        self.batch_cost = batch_cost
        self.seed = seed
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.structures_stream.has_only_formats("pdb"):
            raise ValueError(
                f"FLOWR requires PDB-format structures, got '{self.structures_stream.format}'. "
                "Use convert='pdb' upstream (e.g. PDB(..., convert='pdb')).")
        if (self.ligand_stream is None) == (self.pocket_residues is None):
            raise ValueError("give exactly one of ligand (a reference ligand bound in the structure) "
                             "and pocket_residues (a residue selection)")
        if self.ligand_stream is not None and len(self.ligand_stream) == 0:
            raise ValueError("ligand must not be an empty stream")
        for name in ("num_molecules", "integration_steps", "batch_cost"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be an integer of at least 1, got {value!r}")
        if isinstance(self.seed, bool) or not isinstance(self.seed, int) or self.seed < 0:
            raise ValueError(f"seed must be a non-negative integer, got {self.seed!r}")
        if self.conditioning is not None and self.conditioning not in CONDITIONING:
            raise ValueError(f"conditioning must be None or one of {CONDITIONING}, got {self.conditioning!r}")
        if self.conditioning is not None and self.ligand_stream is None:
            raise ValueError("conditioning fixes atoms of a reference ligand, so it needs ligand, not pocket_residues")
        if self.conditioning == "substructure":
            if not isinstance(self.substructure, str) or not self.substructure.strip():
                raise ValueError('conditioning="substructure" needs substructure, the SMILES of the fragment to keep')
        elif self.substructure is not None:
            raise ValueError('substructure is only used with conditioning="substructure"')
        if self.conditioning is not None and self.sample_sizes_requested:
            raise ValueError("sample_sizes=True is not supported with conditioning: upstream keeps the reference size")
        if self.ligand_stream is not None:
            if self.num_heavy_atoms is not None:
                raise ValueError("num_heavy_atoms sets the size when there is no reference ligand; "
                                 "with ligand the reference sets it (and sample_sizes varies it)")
            if not isinstance(self.pocket_cutoff, (int, float)) or isinstance(self.pocket_cutoff, bool) \
                    or self.pocket_cutoff <= 0:
                raise ValueError(f"pocket_cutoff must be a positive number of A, got {self.pocket_cutoff!r}")
        else:
            if self.pocket_cutoff != 6.0:
                raise ValueError("pocket_cutoff is the distance from the reference ligand, so it needs ligand")
            if isinstance(self.num_heavy_atoms, bool) or not isinstance(self.num_heavy_atoms, int) \
                    or not 5 <= self.num_heavy_atoms <= 80:
                raise ValueError("pocket_residues needs num_heavy_atoms, an integer from 5 to 80 "
                                 f"(upstream's size range), got {self.num_heavy_atoms!r}")
            if not isinstance(self.pocket_residues, str):
                raise ValueError(f"pocket_residues must be a selection string, got {type(self.pocket_residues)}")
            _, residues = sele_to_chains_and_list(self.pocket_residues)
            if any(chain == "" for chain, _ in residues):
                raise ValueError(f"pocket_residues must name the chain of every residue, e.g. 'A10-25+B7', "
                                 f"got {self.pocket_residues!r}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        pocket = ("reference ligand (code resolved from the compounds stream at runtime)"
                  if self.ligand_stream is not None else f"residues {self.pocket_residues}")
        lines.extend([
            f"POCKET: {pocket}",
            f"MOLECULES PER POCKET: {self.num_molecules}",
            f"INTEGRATION STEPS: {self.integration_steps}",
            f"CONDITIONING: {self.conditioning or 'none (de novo)'}",
        ])
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        ligand_flag = ""
        if self.ligand_stream is not None:
            self.ligand_stream.save_json(self.ligand_json)
            ligand_flag = f' \\\n    --ligand-json "{self.ligand_json}"'
        self._write_config_yaml()
        script = "#!/bin/bash\n"
        script += "# FLOWR execution script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Generating {self.num_molecules} molecule(s) per pocket with FLOWR"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-json "{self.structures_json}"{ligand_flag} \\
    --flowr-repo "{self.folders["FLOWR"]}" \\
    --checkpoint "{self.checkpoint}" \\
    --runs-folder "{self.runs_folder}" \\
    --compounds-csv "{self.compounds_csv}" \\
    --structures-folder "{self.stream_folder("structures")}" \\
    --structures-map-csv "{self.structures_map}" \\
    --pockets-csv "{self.pockets_csv}" \\
    --local-missing-csv "{self.local_missing_csv}" \\
    --container-prefix "{self.container_prefix()}"
"""
        script += self.generate_missing_propagation(
            self.structures, self.ligand,
            local_missing=self.local_missing_csv, missing_csv=self.missing_csv,
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        import yaml
        config = {
            "num_molecules": self.num_molecules,
            "integration_steps": self.integration_steps,
            "conditioning": self.conditioning,
            "substructure": self.substructure,
            "sample_sizes": self.sample_sizes,
            "num_heavy_atoms": self.num_heavy_atoms,
            "pocket_cutoff": self.pocket_cutoff,
            "protonate_pocket": self.protonate_pocket,
            "protonate_ligands": self.protonate_ligands,
            "batch_cost": self.batch_cost,
            "seed": self.seed,
            "pocket_residues": self.pocket_residues,
        }
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump(config, f, sort_keys=False)

    def get_output_files(self) -> Dict[str, Any]:
        molecule_ids = generate_multiplied_ids_pattern(
            self.structures_stream.ids, f"<1..{self.num_molecules}>", input_stream_name="structures")
        compounds = DataStream(
            name="compounds",
            ids=molecule_ids,
            files=[],
            map_table=self.compounds_csv,
            format="csv",
        )
        structures = DataStream(
            name="structures",
            ids=molecule_ids,
            files=[self.stream_path("structures", "<id>.sdf")],
            map_table=self.structures_map,
            format="sdf",
        )
        tables = {
            "compounds": TableInfo(
                name="compounds",
                path=self.compounds_csv,
                columns=["id", "structures.id", "sample", "format", "code", "ccd", "smiles"],
                description="Generated molecules, one row per molecule",
            ),
            "pockets": TableInfo(
                name="pockets",
                path=self.pockets_csv,
                columns=["id", "pocket", "n_requested", "n_generated"],
                description="Molecules requested and returned for each input pocket",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "compounds": compounds,
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder,
        }
