# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PLACER tool for ligand-pose and sidechain/apo prediction in protein pockets.

PLACER (Protein-Ligand Atomic Conformational Ensemble Resolver) is a graph
neural network operating at the atomic level. From a partially corrupted input
structure it stochastically regenerates atomic coordinates, producing a
conformational ensemble scored per sample.

This wrapper exposes PLACER's two task modes, selected by the inputs given:

  * ligand-pose mode (``ligand`` provided) — the structure has the ligand bound
    as HETATM; PLACER corrupts and resamples that ligand's pose (plus pocket
    sidechains), keeping any other ligands fixed. Returns the full score set
    including the ligand-accuracy terms (prmsd, rmsd, kabsch).
  * sidechain/apo mode (``ligand`` omitted, ``target_res`` provided) — a protein
    residue is the crop center and PLACER resamples the surrounding sidechains.
    No ligand is predicted, so only the lDDT-style confidences (plddt,
    plddt_pde, fape) are returned.

Reference:
    Modeling protein-small molecule conformational ensembles with PLACER.
    PNAS 122(45):e2427161122 (2025). https://doi.org/10.1073/pnas.2427161122
    https://github.com/baker-laboratory/PLACER
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


class PLACER(BaseConfig):
    """
    PLACER: ligand-pose and sidechain/apo prediction in a protein pocket.

    Mode is inferred from the inputs:
      * ``ligand`` given  -> ligand-pose mode. Output IDs are the cartesian
        product ``<structure>+<ligand>`` multiplied by a sample index
        (``<structure>+<ligand>_<1..nsamples>``).
      * ``ligand`` omitted + ``target_res`` given -> sidechain/apo mode. Output
        IDs are ``<structure>_<1..nsamples>`` (no ligand axis).
    One predicted PDB per sample either way.

    Inputs:
        structures:   protein structures (PDB/mmCIF). PLACER parses HETATM
                      coordinates directly — no SDF staging needed. (mmCIF is
                      parsed correctly only for RCSB-sourced files; PDB is the
                      safe default.)
        ligand:       (ligand-pose mode) compounds stream carrying the 3-letter
                      residue ``code`` of the bound ligand. The driver resolves
                      the single ``code`` at runtime and passes it to PLACER's
                      ``predict_ligand`` selector (by ``name3``: every copy of
                      that residue is predicted, the rest fixed). Typically
                      ``Ligand("HEM")`` or ``Ligand("HEM", smiles="...")``.

    Parameters:
        target_res (str):  (sidechain/apo mode) protein residue used as the crop
                           center, ``chain-resno`` e.g. ``"A-149"``. Required
                           when ``ligand`` is omitted; rejected when it is given.
        exclude_sm (bool): drop ALL small molecules from the prediction (true
                           apo). Default False.
        nsamples (int):    number of samples per input (PLACER default 10;
                           upstream recommends 50-100).
        rerank (str):      optional metric to rank output models/CSV by, one of
                           {"prmsd", "plddt", "plddt_pde"}. None (default) keeps
                           PLACER's native order. (prmsd is unavailable in apo
                           mode — rank by plddt / plddt_pde there.)

    Outputs:
        Streams:
            structures:  predicted PDB per sample (prmsd in the B-factor column)
            compounds:   (ligand mode only) chemistry passthrough of the input
                         ligand (code/SMILES unchanged — pose is refined, not
                         identity)
        Tables:
            scores:      ligand mode -> id | structures.id | compounds.id |
                         sample | prmsd | plddt | plddt_pde | fape | rmsd | kabsch
                         apo mode    -> id | structures.id | sample | plddt |
                         plddt_pde | fape
            missing:     id | removed_by | cause (inputs PLACER could not process)

    Usage::

        # ligand-pose mode
        with Pipeline(project="Examples", job="PLACER-dock"):
            Resources(gpu="A100", memory="32GB", time="2:00:00")
            holo = PDB("4dtz")                     # ligand bound as HETATM
            lig = Ligand("LDP")                    # residue code present in 4dtz
            placer = PLACER(structures=holo, ligand=lig, nsamples=50, rerank="prmsd")

        # sidechain/apo mode
        with Pipeline(project="Examples", job="PLACER-sidechain"):
            Resources(gpu="A100", memory="32GB", time="2:00:00")
            apo = PDB("dnHEM1_apo")
            placer = PLACER(structures=apo, target_res="A-149", nsamples=50, rerank="plddt")
    """

    TOOL_NAME = "PLACER"
    TOOL_VERSION = "1.1"

    RERANK_CHOICES = ("prmsd", "plddt", "plddt_pde")

    # ------------------------------------------------------------------
    # Install
    # ------------------------------------------------------------------

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Clone baker-laboratory/PLACER and create its conda env.

        Model weights ship inside the repo (``weights/PLACER_model_1.pt``) — no
        separate download step. Generation runs through our own driver
        (pipe_placer_driver.py) that imports the upstream ``PLACER`` package
        from the clone at runtime; upstream is left pristine.

        Verification: env exists + repo present + ``import PLACER`` succeeds +
        the default weights file is on disk.
        """
        repo_dir = folders.get("PLACER", "")
        parent_dir = os.path.dirname(repo_dir)
        biopipelines = folders.get("biopipelines", "")
        weights_path = os.path.join(repo_dir, "weights", "PLACER_model_1.pt")

        env_check = cls._env_exists_check("placer", env_manager)
        repo_check = f'[ -f "{repo_dir}/PLACER.py" ]'
        weights_check = f'[ -f "{weights_path}" ]'
        # Verify the heavy binary deps are actually present, not just that the
        # env exists. A corrupt conda package cache can leave `mamba env create`
        # half-finished (torch/dgl never linked) while the run still reports
        # success — importing them explicitly makes that fail loudly. The
        # imports are echoed (no >/dev/null) so a failure is visible in the log.
        import_check = (
            f'{env_manager} run -n placer python -c '
            f'"import torch, dgl, e3nn; import sys; sys.path.insert(0, \'{repo_dir}\'); import PLACER; '
            f'print(\'PLACER verify OK; torch\', torch.__version__)"'
        )

        skip = "" if force_reinstall else f"""# Check if already installed
if {repo_check} && {env_check} && {weights_check} && {import_check}; then
    echo "PLACER already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        clone_block = f"""mkdir -p {parent_dir}
cd {parent_dir}
if [ ! -d "{repo_dir}" ]; then
    git clone https://github.com/baker-laboratory/PLACER.git "{repo_dir}"
fi"""

        remove_block = cls._env_remove_block("placer", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("placer", env_manager, biopipelines)

        return f"""echo "=== Installing PLACER ==="
{skip}{clone_block}

{remove_block}
{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create placer environment."
    exit 1
fi

# Verify installation
if {repo_check} && {weights_check} && {import_check}; then
    touch "$INSTALL_SUCCESS"
    echo "=== PLACER installation complete ==="
else
    echo "ERROR: PLACER verification failed (repo, weights, or 'import PLACER' missing)"
    exit 1
fi
"""

    # ------------------------------------------------------------------
    # Lazy path descriptors
    # ------------------------------------------------------------------

    structures_json = Path(lambda self: self.configuration_path(".input_structures.json"))
    ligand_json = Path(lambda self: self.configuration_path(".input_ligand.json"))
    bonds_json = Path(lambda self: self.configuration_path(".input_bonds.json"))
    runs_folder = Path(lambda self: self.execution_path("runs"))
    structures_map = Path(lambda self: self.stream_map_path("structures"))
    compounds_map = Path(lambda self: self.stream_map_path("compounds"))
    scores_csv = Path(lambda self: self.table_path("scores"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    driver_py = Path(lambda self: self.pipe_script_path("pipe_placer_driver.py"))
    postprocess_py = Path(lambda self: self.pipe_script_path("pipe_placer_postprocess.py"))

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(
        self,
        structures: Union[DataStream, StandardizedOutput],
        ligand: Optional[Union[DataStream, StandardizedOutput]] = None,
        target_res: Optional[str] = None,
        exclude_sm: bool = False,
        nsamples: int = 10,
        rerank: Optional[str] = None,
        bonds: Optional[Union[tuple, List[tuple]]] = None,
        **kwargs,
    ):
        # Keep original inputs for predict_output_ids_with_provenance — the
        # helper plucks axis IDs by stream name from the StandardizedOutput.
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

        # Mode is inferred from the inputs (see validate_params):
        #   ligand given            -> "ligand"   (predict_ligand by name3)
        #   ligand omitted + target -> "sidechain" (apo / sidechain repacking)
        if ligand is None:
            self.ligand_stream: Optional[DataStream] = None
            self.mode = "sidechain"
        else:
            if isinstance(ligand, StandardizedOutput):
                self.ligand_stream = ligand.streams.compounds
            elif isinstance(ligand, DataStream):
                self.ligand_stream = ligand
            else:
                raise ValueError(
                    f"ligand must be DataStream or StandardizedOutput, got {type(ligand).__name__}"
                )
            self.mode = "ligand"

        self.target_res = target_res
        self.exclude_sm = exclude_sm
        self.nsamples = nsamples
        self.rerank = rerank

        # Optional covalent bonds to enforce during resampling. Each is
        # (atom1, atom2, length) with atoms in '<residue>.<atom>' syntax, e.g.
        # ("A145.SG", "LIG.C12", 1.8) to tether a dye to a catalytic cysteine.
        # Normalize to a list of tuples.
        if bonds is None:
            self.bonds = []
        elif isinstance(bonds, tuple) and bonds and not isinstance(bonds[0], (tuple, list)):
            self.bonds = [bonds]          # a single (a1, a2, len) tuple
        else:
            self.bonds = list(bonds)
        self.list_unmatched_bonds = []

        super().__init__(**kwargs)

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        # PLACER reads PDB and RCSB mmCIF directly; the "pdb|cif" no-conversion
        # default (from PDB/RCSB with convert=None) is fine — the driver
        # dispatches on each file's extension at runtime.
        if not self.structures_stream.has_only_formats("pdb", "cif"):
            raise ValueError(
                f"PLACER requires PDB or mmCIF structures, got '{self.structures_stream.format}'."
            )

        # Mode-specific input requirements.
        if self.mode == "ligand":
            if not self.ligand_stream or len(self.ligand_stream) == 0:
                raise ValueError("ligand parameter is required and must not be empty")
            if self.target_res is not None:
                raise ValueError(
                    "target_res is only valid in sidechain/apo mode (ligand omitted); "
                    "in ligand-pose mode the bound ligand defines the crop center"
                )
        else:  # sidechain / apo
            if self.target_res is None:
                raise ValueError(
                    "PLACER needs either a ligand (ligand-pose mode) or target_res "
                    "(sidechain/apo mode). Provide one."
                )
            _validate_freeform_string("target_res", self.target_res)

        if not isinstance(self.exclude_sm, bool):
            raise ValueError(f"exclude_sm must be a bool, got {type(self.exclude_sm).__name__}")
        if not isinstance(self.nsamples, int) or self.nsamples <= 0:
            raise ValueError(f"nsamples must be a positive integer, got {self.nsamples!r}")
        if self.rerank is not None and self.rerank not in self.RERANK_CHOICES:
            raise ValueError(
                f"rerank must be one of {self.RERANK_CHOICES} or None, got {self.rerank!r}"
            )
        if self.bonds:
            if self.mode != "ligand":
                raise ValueError("bonds= is only valid in ligand-pose mode (a ligand to tether)")
            for b in self.bonds:
                if not (isinstance(b, (tuple, list)) and len(b) == 3):
                    raise ValueError(
                        f"each bond must be (atom1, atom2, length), got {b!r}"
                    )
                a1, a2, blen = b
                for a in (a1, a2):
                    if not isinstance(a, str) or "." not in a:
                        raise ValueError(
                            f"bond atoms must be '<residue>.<atom>' strings, got {a!r}"
                        )
                if not isinstance(blen, (int, float)):
                    raise ValueError(f"bond length must be numeric, got {blen!r}")

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
        lines.append(f"MODE: {self.mode}")
        lines.append(f"NSAMPLES: {self.nsamples}")
        lines.append(f"RERANK: {self.rerank if self.rerank is not None else '(none)'}")
        lines.append(f"EXCLUDE_SM: {self.exclude_sm}")
        if self.mode == "ligand":
            lines.append("LIGAND CODE: (resolved from ligand stream at runtime)")
        else:
            lines.append(f"TARGET_RES: {self.target_res}")
        return lines

    # ------------------------------------------------------------------
    # Script generation
    # ------------------------------------------------------------------

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        if self.mode == "ligand":
            self.ligand_stream.save_json(self.ligand_json)

        rerank_arg = f' --rerank "{self.rerank}"' if self.rerank is not None else ""
        exclude_sm_arg = " --exclude-sm" if self.exclude_sm else ""
        if self.mode == "ligand":
            mode_args = f' --mode ligand --ligand-json "{self.ligand_json}"'
        else:
            mode_args = f' --mode sidechain --target-res "{self.target_res}"'

        bonds_arg = ""
        if self.bonds:
            import json as _json
            with open(self.bonds_json, "w") as _bf:
                _json.dump([[a1, a2, blen] for (a1, a2, blen) in self.bonds], _bf)
            bonds_arg = f' --bonds-json "{self.bonds_json}"'

        script = "#!/bin/bash\n"
        script += "# PLACER execution script\n"
        script += self.generate_completion_check_header()
        script += self.warn_container_unsupported()

        # Phase 1: run PLACER per input via our own driver (placer env). The
        # driver imports the upstream PLACER package from the clone at runtime
        # and writes one {input}_model.pdb + {input}.csv per input.
        script += "# --- PLACER inference (placer env, our driver) ---\n"
        script += self.activate_environment()  # placer
        script += f"""python "{self.driver_py}" \\
    --placer-repo "{self.folders["PLACER"]}" \\
    --structures-json "{self.structures_json}" \\
    --runs-folder "{self.runs_folder}" \\
    --nsamples {self.nsamples}{rerank_arg}{exclude_sm_arg}{mode_args}{bonds_arg}

if [ $? -ne 0 ]; then
    echo "Error: PLACER inference failed"
    exit 1
fi

"""

        # Phase 2: collect ranked PDBs + scores into streams (biopipelines env).
        # The compounds passthrough is only written in ligand mode.
        compounds_args = (
            f' --ligand-json "{self.ligand_json}" --compounds-map "{self.compounds_map}"'
            if self.mode == "ligand" else ""
        )
        script += "# --- Post-processing (biopipelines env) ---\n"
        script += self.activate_environment(name="biopipelines")
        script += f"""python "{self.postprocess_py}" \\
    --mode {self.mode} \\
    --runs-folder "{self.runs_folder}" \\
    --structures-folder "{self.stream_folder("structures")}" \\
    --structures-map "{self.structures_map}" \\
    --scores-csv "{self.scores_csv}" \\
    --missing-csv "{self.missing_csv}"{compounds_args}

if [ $? -ne 0 ]; then
    echo "Error: PLACER post-processing failed"
    exit 1
fi

"""
        script += self.generate_completion_check_footer()
        return script

    # ------------------------------------------------------------------
    # Output prediction
    # ------------------------------------------------------------------

    def get_output_files(self) -> Dict[str, Any]:
        sample_pattern = f"<1..{self.nsamples}>"

        if self.mode == "ligand":
            # Two-axis output: <structure>+<ligand>, multiplied by sample index.
            base_ids, _ = predict_output_ids_with_provenance(
                structures=(self.structures_input, "structures"),
                compounds=(self.ligand_input, "compounds"),
            )
            score_columns = [
                "id", "structures.id", "compounds.id", "sample",
                "prmsd", "plddt", "plddt_pde", "fape", "rmsd", "kabsch",
            ]
        else:
            # Sidechain/apo: single axis (structures), no ligand prediction, so
            # no ligand-accuracy metrics (rmsd/kabsch/prmsd are not returned).
            base_ids = self.structures_stream.ids
            score_columns = [
                "id", "structures.id", "sample", "plddt", "plddt_pde", "fape",
            ]

        sample_ids = generate_multiplied_ids_pattern(
            base_ids, sample_pattern, input_stream_name="structures",
        )

        structures = DataStream(
            name="structures",
            ids=sample_ids,
            files=[self.stream_path("structures", "<id>.pdb")],
            map_table=self.structures_map,
            format="pdb",
        )

        tables = {
            "scores": TableInfo(
                name="scores",
                path=self.scores_csv,
                columns=score_columns,
                description="PLACER per-sample scores",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Inputs PLACER could not process",
            ),
        }

        out: Dict[str, Any] = {
            "structures": structures,
            "tables": tables,
            "output_folder": self.output_folder,
        }

        if self.mode == "ligand":
            # Chemistry passthrough: PLACER refines the pose but does not rename
            # the ligand, so ids/code/SMILES carry through unchanged.
            out["compounds"] = DataStream(
                name="compounds",
                ids=self.ligand_stream.ids,
                files=[],
                map_table=self.compounds_map,
                format="csv",
            )

        return out
