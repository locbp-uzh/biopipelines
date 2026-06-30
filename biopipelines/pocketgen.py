# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PocketGen tool for ligand-aware full-atom pocket co-design.

PocketGen takes a (scaffold protein, target ligand) pair and jointly designs
the pocket residues' sequence and side-chain atomic positions so the resulting
pocket accommodates the ligand. One designed pocket per pair.

Reference:
    Zhang et al. (2024) Efficient generation of protein pockets with PocketGen.
    Nature Machine Intelligence. https://github.com/zaixizhang/PocketGen
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import predict_output_ids_with_provenance, generate_multiplied_ids_pattern
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import predict_output_ids_with_provenance, generate_multiplied_ids_pattern


# Google Drive file ID for the pretrained checkpoint (README link:
# https://drive.google.com/file/d/1cuvdiu3bXyni71A2hoeZSWT1NOsNfeD_/view).
POCKETGEN_CKPT_GDRIVE_ID = "1cuvdiu3bXyni71A2hoeZSWT1NOsNfeD_"


class PocketGen(BaseConfig):
    """
    PocketGen: full-atom pocket co-design conditioned on a target ligand.

    For each (scaffold, ligand) pair, predicts a redesigned pocket. Output IDs
    are the cartesian product <scaffold>+<ligand>; one sequence and one
    relaxed PDB per pair.

    Inputs:
        structures:   scaffold protein backbones with the ligand bound as
                      HETATM in the binding site. The bound coordinates
                      define the "edit pocket" PocketGen designs.
        ligand:       compounds stream carrying the canonical SMILES and the
                      3-letter `code` of the bound ligand. The stage script
                      uses ProLIF's idiom: extract the HETATM block matching
                      the `code` out of the scaffold PDB (keeps docked
                      coordinates) and apply the SMILES as a bond-order template
                      via AssignBondOrdersFromTemplate before writing the SDF
                      PocketGen consumes. Both `code` and `smiles` are read from
                      the stream at runtime (`code` must resolve to a single
                      value). Typically supplied as ``Ligand("CODE", smiles="...")``.

    PocketGen's preprocessing selects pocket residues by proximity to the
    ligand atoms in 3D space, so the staged SDF must have the bound
    coordinates. Pre-docked SDFs (e.g. DiffDock output as the ligand stream)
    are passed through verbatim.

    Outputs:
        Streams:
            structures:  relaxed PDB per pair (designed pocket grafted into the scaffold)
            sequences:   designed pocket sequence per pair
                         (content-bearing map_table: id | structures.id | compounds.id | sequence)
        Tables:
            sequences:   same columns as the stream's map_table
            missing:     id | removed_by | cause

    Usage::

        with Pipeline(project="Examples", job="PocketGen-demo"):
            Resources(gpu="A100", memory="32GB", time="2:00:00")
            scaffold = PDB("2p16", convert="pdb")     # ligand bound as HETATM
            lig = Ligand("GG2", smiles="...")          # CCD code + canonical SMILES
            pg = PocketGen(structures=scaffold, ligand=lig)
    """

    TOOL_NAME = "PocketGen"
    TOOL_VERSION = "1.1"

    # ------------------------------------------------------------------
    # Install
    # ------------------------------------------------------------------

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Clone zaixizhang/PocketGen, create its conda env, download the
        pretrained checkpoint from Google Drive.

        We do NOT patch upstream sources — generation runs through our own
        driver (pipe_pocketgen_driver.py) that imports utils.*/models.*
        directly. This keeps upstream pristine and avoids breakage on every
        upstream cosmetic change.

        Verification: env exists + torch + rdkit + esm import + checkpoint
        file present.
        """
        biopipelines = folders.get("biopipelines", "")
        repo_dir = folders.get("PocketGen", "")
        parent_dir = os.path.dirname(repo_dir)
        ckpt_path = os.path.join(repo_dir, "checkpoints", "pocketgen.pt")

        env_check = cls._env_exists_check("pocketgen", env_manager)
        repo_check = f'[ -d "{repo_dir}/utils" ] && [ -d "{repo_dir}/models" ]'
        ckpt_check = f'[ -f "{ckpt_path}" ]'
        skip = "" if force_reinstall else f"""# Check if already installed
if {repo_check} && {env_check} && {ckpt_check} \\
   && {env_manager} run -n pocketgen python -c "import torch, rdkit, esm" >/dev/null 2>&1; then
    echo "PocketGen already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        clone_block = f"""mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/zaixizhang/PocketGen.git "{repo_dir}"
fi"""

        remove_block = cls._env_remove_block("pocketgen", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("pocketgen", env_manager, biopipelines)

        ckpt_block = f"""# Fetch pretrained checkpoint (~hundreds of MB) from Google Drive.
mkdir -p "$(dirname {ckpt_path})"
if [ ! -f "{ckpt_path}" ]; then
    {env_manager} run -n pocketgen gdown --id {POCKETGEN_CKPT_GDRIVE_ID} -O "{ckpt_path}"
fi"""

        return f"""echo "=== Installing PocketGen ==="
{skip}{clone_block}

{remove_block}
{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create pocketgen environment."
    exit 1
fi

{ckpt_block}

# Verify installation
if {repo_check} && {ckpt_check} \\
   && {env_manager} run -n pocketgen python -c "import torch, rdkit, esm" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== PocketGen installation complete ==="
else
    echo "ERROR: PocketGen verification failed (repo, checkpoint, or imports missing)"
    exit 1
fi
"""

    # ------------------------------------------------------------------
    # Lazy path descriptors
    # ------------------------------------------------------------------

    structures_json = Path(lambda self: self.configuration_path(".input_structures.json"))
    ligand_json = Path(lambda self: self.configuration_path(".input_ligand.json"))
    staging_folder = Path(lambda self: self.configuration_path("staging"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    sequences_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    stage_py = Path(lambda self: self.pipe_script_path("pipe_pocketgen_stage.py"))
    driver_py = Path(lambda self: self.pipe_script_path("pipe_pocketgen_driver.py"))
    postprocess_py = Path(lambda self: self.pipe_script_path("pipe_pocketgen_postprocess.py"))
    pg_default_config = Path(lambda self: os.path.join(self.folders["PocketGen"], "configs", "train_model.yml"))
    pg_checkpoint = Path(lambda self: os.path.join(self.folders["PocketGen"], "checkpoints", "pocketgen.pt"))

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(
        self,
        structures: Union[DataStream, StandardizedOutput],
        ligand: Union[DataStream, StandardizedOutput],
        **kwargs,
    ):
        # Keep original inputs for predict_output_ids_with_provenance — the
        # helper plucks axis IDs by stream name from the StandardizedOutput,
        # so a bare DataStream is the wrong shape.
        self.structures_input = structures
        self.ligand_input = ligand

        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(
                f"structures must be DataStream or StandardizedOutput, got {type(structures).__name__}"
            )

        if isinstance(ligand, StandardizedOutput):
            self.ligand_stream: DataStream = ligand.streams.compounds
        elif isinstance(ligand, DataStream):
            self.ligand_stream = ligand
        else:
            raise ValueError(
                f"ligand must be DataStream or StandardizedOutput, got {type(ligand).__name__}"
            )

        super().__init__(**kwargs)

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.structures_stream.has_only_formats("pdb"):
            raise ValueError(
                f"PocketGen requires PDB-format structures, got '{self.structures_stream.format}'. "
                "Use convert='pdb' upstream."
            )
        if not self.ligand_stream or len(self.ligand_stream) == 0:
            raise ValueError("ligand parameter is required and must not be empty")

    # ------------------------------------------------------------------
    # Configure inputs
    # ------------------------------------------------------------------

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    # ------------------------------------------------------------------
    # Config display
    # ------------------------------------------------------------------

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append("LIGAND CODE: (resolved from ligand stream at runtime)")
        return lines

    # ------------------------------------------------------------------
    # Script generation
    # ------------------------------------------------------------------

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self.ligand_stream.save_json(self.ligand_json)

        script = "#!/bin/bash\n"
        script += "# PocketGen execution script\n"
        script += self.generate_completion_check_header()

        # Phase 1: stage pairs into the layout PocketGen wants (biopipelines env).
        script += "# --- Stage (scaffold, ligand) pairs (biopipelines env) ---\n"
        script += self.activate_environment(name="biopipelines")
        script += f"""python "{self.stage_py}" \\
    --structures-json "{self.structures_json}" \\
    --ligand-json "{self.ligand_json}" \\
    --staging-folder "{self.staging_folder}" \\
    --extras-dir "{self.extras_path()}"

if [ $? -ne 0 ]; then
    echo "Error: PocketGen staging failed"
    exit 1
fi

"""

        # Phase 2: run PocketGen via our own driver (pocketgen env). The
        # driver lives in pipe_scripts/ and imports utils.*/models.* from the
        # cloned repo at runtime — no upstream patching.
        script += "# --- PocketGen inference (pocketgen env, our driver) ---\n"
        script += self.activate_environment()  # pocketgen
        script += f"""{self.container_prefix()}python "{self.driver_py}" \\
    --pocketgen-repo "{self.folders["PocketGen"]}" \\
    --config "{self.pg_default_config}" \\
    --checkpoint "{self.pg_checkpoint}" \\
    --target "{self.staging_folder}"

if [ $? -ne 0 ]; then
    echo "Error: PocketGen inference failed"
    exit 1
fi

"""

        # Phase 3: post-process (biopipelines env).
        script += "# --- Post-processing (biopipelines env) ---\n"
        script += self.activate_environment(name="biopipelines")
        script += f"""python "{self.postprocess_py}" \\
    --staging-folder "{self.staging_folder}" \\
    --structures-folder "{self.stream_folder("structures")}" \\
    --sequences-folder "{self.stream_folder("sequences")}" \\
    --structures-map "{self.structures_map}" \\
    --sequences-csv "{self.sequences_csv}" \\
    --missing-csv "{self.missing_csv}"

if [ $? -ne 0 ]; then
    echo "Error: PocketGen post-processing failed"
    exit 1
fi

"""
        script += self.generate_completion_check_footer()
        return script

    # ------------------------------------------------------------------
    # Output prediction
    # ------------------------------------------------------------------

    # Upstream PocketGen hardcodes 8 samples per (scaffold, ligand) pair —
    # see generate_new.py: `datalist = [data for _ in range(8)]`.
    NUM_SAMPLES_PER_PAIR = 8

    def get_output_files(self) -> Dict[str, Any]:
        pair_ids, _ = predict_output_ids_with_provenance(
            structures=(self.structures_input, "structures"),
            compounds=(self.ligand_input, "compounds"),
        )

        # Multiply pair IDs by sample index 0..7 (compact pattern, lazy-safe).
        sample_pattern = f"<0..{self.NUM_SAMPLES_PER_PAIR - 1}>"
        sample_ids = generate_multiplied_ids_pattern(
            pair_ids, sample_pattern, input_stream_name="structures",
        )

        structures = DataStream(
            name="structures",
            ids=sample_ids,
            files=[self.stream_path("structures", "<id>.pdb")],
            map_table=self.structures_map,
            format="pdb",
        )

        # Content-bearing sequences stream: the CSV holds the lineage
        # (id -> scaffold/ligand provenance + sample index) and the designed
        # sequence.
        sequences = DataStream(
            name="sequences",
            ids=sample_ids,
            files=[],
            map_table=self.sequences_csv,
            format="csv",
        )

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "structures.id", "compounds.id", "sample", "sequence"],
                description="PocketGen designed pocket sequences (8 samples per pair)",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Pairs PocketGen could not process",
            ),
        }

        return {
            "structures": structures,
            "sequences": sequences,
            "tables": tables,
            "output_folder": self.output_folder,
        }
