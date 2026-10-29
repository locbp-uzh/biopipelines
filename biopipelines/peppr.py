# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Peppr scores predicted biomolecular complexes against reference structures.

Wraps the peppr package (https://github.com/aivant/peppr), which matches the atoms of a prediction to its reference, including chain permutations and symmetric ligands, and then computes the CASP/CAPRI family of metrics.
"""

import os
from typing import Dict, List, Any, Optional, Union

# The metric keys the tool exposes, in table order; pipe_peppr.py builds one peppr metric for each.
PEPPR_METRICS = ("rmsd", "tm_score", "lddt", "lddt_pli", "fnat", "irmsd", "lrmsd", "dockq")

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream


class Peppr(BaseConfig):
    """
    Peppr: score predicted structures against reference structures with the peppr package.

    Atoms of each prediction are matched to its reference (chain permutations and symmetric ligands are resolved by peppr), then each requested metric is computed. A metric peppr cannot compute for a system, such as an interface metric on a monomer, is NaN. Nothing is filtered on a value; apply Panda.filter downstream.

    Inputs:
        structures: predicted structures (PDB or mmCIF).
        references: reference structures (PDB or mmCIF). A stream with one id scores every prediction; a stream with several is matched to the predictions by id (the same id, or a parent or child of it), so a refold finds the design it derives from. A sibling id is not a match. A prediction with no matching reference is not evaluated and is recorded in `missing`.
        metrics: keys to compute, a list drawn from the table below (default: all eight).

    Every molecule of a structure, polymer chain or ligand, must have its own chain id, and bonds must be known. biotite reads them: polymer bonds come from the Chemical Component Dictionary, and ligand bonds from CONECT records in a PDB file or from chem_comp_bond (else the dictionary) in an mmCIF file. A ligand in a PDB file without CONECT records has no bonds, so its id cannot be matched and lands in `missing`; give such structures as mmCIF. Hydrogens, solvent and monoatomic ions are dropped.

    Metrics (key: peppr class, what it reports):
        rmsd:     MonomerRMSD, backbone RMSD (Calpha, and C3' for nucleic acids) of each polymer chain after superposing it alone, averaged over chains weighted by atom count. Angstrom.
        tm_score: MonomerTMScore, TM-score of each protein chain, averaged the same way.
        lddt:     MonomerLDDTScore, lDDT of the contacts inside each polymer chain, averaged the same way. Contacts between chains are not included.
        lddt_pli: LDDTPLIScore, CASP lDDT-PLI, the lDDT of the polymer-ligand interactions. NaN when no ligand touches the polymer in the reference.
        fnat:     ContactFraction, fraction of the reference's contacts between two polymer chains that the prediction reproduces, averaged over chain pairs.
        irmsd:    InterfaceRMSD, interface RMSD between two polymer chains, averaged over pairs. Angstrom.
        lrmsd:    LigandRMSD, backbone RMSD of the shorter chain of each polymer pair in contact after superposing the longer one, averaged over pairs. Angstrom.
        dockq:    DockQScore, DockQ of each chain pair, polymer-ligand pairs included, averaged over pairs. 0 to 1, higher is better.
        The pair metrics (fnat, irmsd, lrmsd, dockq) are NaN for a structure with a single chain.

    Outputs:
        Streams: (none)
        Tables:
            peppr:   id | <one column per requested metric>
            missing: id | removed_by | kind | cause
    """

    TOOL_NAME = "Peppr"
    TOOL_VERSION = "1.1"
    ENV_NAME = "peppr"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        env = cls._install_env(env_manager)
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check(env, env_manager)
        env_run = cls._env_run(env, env_manager)
        imports = f'{env_run}python -c "import peppr, biotite, pandas, yaml"'
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check} && {imports} >/dev/null 2>&1; then
    echo "Peppr already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
        env_block = cls._env_install_block(env, env_manager, biopipelines)
        return f"""echo "=== Installing Peppr ==="
{skip}{remove_block}
if ! {env_check}; then
{env_block}
else
    echo "{env} environment already exists, skipping creation."
fi

# Importing peppr is the check that matters: it fails outright on a biotite it does not support.
if {imports} >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== Peppr installation complete ==="
else
    echo "ERROR: Peppr verification failed (import of peppr, biotite, pandas or yaml). Retry with force_reinstall=True to rebuild the env."
    {imports} 2>&1 | tail -5
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    references_json = Path(lambda self: self.configuration_path("references.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    peppr_csv = Path(lambda self: self.table_path("peppr"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_peppr.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 references: Union[DataStream, StandardizedOutput],
                 metrics: Optional[List[str]] = None,
                 **kwargs):
        self.structures = structures
        self.references = references
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        if isinstance(references, StandardizedOutput):
            self.references_stream: DataStream = references.streams.structures
        elif isinstance(references, DataStream):
            self.references_stream = references
        else:
            raise ValueError(f"references must be DataStream or StandardizedOutput, got {type(references)}")
        self.metrics = list(PEPPR_METRICS) if metrics is None else metrics
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.references_stream or len(self.references_stream) == 0:
            raise ValueError("references parameter is required and must not be empty")
        if not isinstance(self.metrics, (list, tuple)) or not self.metrics:
            raise ValueError(
                f"metrics must be a non-empty list drawn from {list(PEPPR_METRICS)}, got {self.metrics!r}")
        unknown = [m for m in self.metrics if m not in PEPPR_METRICS]
        if unknown:
            raise ValueError(f"unknown metric(s) {unknown}; supported: {list(PEPPR_METRICS)}")
        if len(set(self.metrics)) != len(self.metrics):
            raise ValueError(f"metrics lists a metric more than once: {list(self.metrics)}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        pairing = ("one reference for every structure" if len(self.references_stream) == 1
                   else "references matched to structures by id")
        lines.extend([
            f"STRUCTURES: {len(self.structures_stream)} files",
            f"REFERENCES: {len(self.references_stream)} files ({pairing})",
            f"METRICS: {', '.join(self.metrics)}",
        ])
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self.references_stream.save_json(self.references_json)
        self._write_config_yaml()
        script = "#!/bin/bash\n"
        script += "# Peppr scoring script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += self.warn_container_unsupported()
        script += f"""echo "Scoring {len(self.structures_stream)} structure(s) against {len(self.references_stream)} reference(s) with peppr"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-json "{self.structures_json}" \\
    --references-json "{self.references_json}" \\
    --peppr-csv "{self.peppr_csv}" \\
    --local-missing-csv "{self.local_missing_csv}"
"""
        script += self.generate_missing_propagation(
            self.structures, self.references,
            local_missing=self.local_missing_csv, missing_csv=self.missing_csv,
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        import yaml
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump({"metrics": list(self.metrics)}, f, sort_keys=False)

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "peppr": TableInfo(
                name="peppr",
                path=self.peppr_csv,
                columns=["id", *self.metrics],
                description="peppr metrics of each predicted structure against its reference",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
