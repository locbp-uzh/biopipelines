# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""IPSAE tool: scores predicted interfaces with ipSAE by running the upstream ipsae.py.

ipSAE condenses the predicted aligned error (PAE) matrix of a predicted complex into one confidence value per chain pair, counting only the residue pairs whose PAE is below a cutoff. The script also reports ipTM recomputed from the PAE, pDockQ, pDockQ2 and LIS. The formulas are upstream's: this wrapper stages each structure and PAE file, runs the unmodified script per id and parses what it writes.

Reference: Dunbrack (2025) Rēs ipSAE loquuntur: What's wrong with AlphaFold's ipTM score and how to fix it. bioRxiv, https://www.biorxiv.org/content/10.1101/2025.02.10.637595. Code: https://github.com/DunbrackLab/IPSAE
"""

import os
from typing import Dict, List, Any, Optional, Union

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

IPSAE_SCRIPT_URL = "https://raw.githubusercontent.com/DunbrackLab/IPSAE/main/ipsae.py"

STRUCTURE_FORMATS = ("pdb", "cif")
PAE_FORMATS = ("npz", "json")

# The score columns upstream prints per chain pair, lowercased.
SCORE_COLUMNS = ["ipsae", "ipsae_d0chn", "ipsae_d0dom", "iptm_af", "iptm_d0chn", "pdockq", "pdockq2", "lis"]
DETAIL_COLUMNS = ["n0res", "n0chn", "n0dom", "d0res", "d0chn", "d0dom", "nres1", "nres2", "dist1", "dist2"]
PAIR_COLUMNS = ["id", "chain1", "chain2", "type", "pae_cutoff", "dist_cutoff"] + SCORE_COLUMNS + DETAIL_COLUMNS
SUMMARY_COLUMNS = ["id", "binder_chain", "target_chain"] + SCORE_COLUMNS


class IPSAE(BaseConfig):
    """
    IPSAE: ipSAE interface scores for predicted complexes.

    Runs Dunbrack's `ipsae.py` (https://github.com/DunbrackLab/IPSAE) once per structure. Upstream reads a PAE matrix plus the structure and prints, for each chain pair, one row per direction (`type` = asym: the PAE of chain 2's residues with the structure aligned on chain 1's residues) and one row with their maximum (`type` = max).

    Inputs:
        structures: predicted complexes, PDB or mmCIF.
        pae: per-id PAE matrices, from a predictor that exposes a `pae` stream (`predictor.streams.pae`) or a DataStream of PAE files. Format `npz` (Boltz layout, key `pae`) or `json` (ColabFold or AlphaFold3 layout, key `pae` or `predicted_aligned_error`). A StandardizedOutput without a `pae` stream is refused.
        binder_chain: the chain whose interface is summarized (default "A").
        target_chains: chain id or list of chain ids scored against `binder_chain` (default None: every other chain).
        pae_cutoff: residue pairs with PAE below this value (A) enter the score (default 10.0).
        dist_cutoff: distance cutoff (A) upstream uses to count interface residues (`dist1`, `dist2`; default 10.0).
    The defaults for `binder_chain`, `target_chains` and `pae_cutoff` are those of the `ipsae` block of RosettaCommons Chaperone (https://github.com/RosettaCommons/chaperone); `dist_cutoff` is the 10 A used in upstream's AlphaFold3 and Boltz examples.

    Outputs:
        Tables:
            ipsae:        id | binder_chain | target_chain | ipsae | ipsae_d0chn | ipsae_d0dom | iptm_af | iptm_d0chn | pdockq | pdockq2 | lis
                          exactly one row per id that entered. The row is the `max`-type result for `binder_chain` against the target chain with the highest ipsae; `target_chain` names it. Every column is NaN when the id failed.
            ipsae_pairs:  id | chain1 | chain2 | type | pae_cutoff | dist_cutoff | ipsae | ipsae_d0chn | ipsae_d0dom | iptm_af | iptm_d0chn | pdockq | pdockq2 | lis | n0res | n0chn | n0dom | d0res | d0chn | d0dom | nres1 | nres2 | dist1 | dist2
                          every row upstream printed: both directions of each chain pair (`asym`) and their maximum (`max`), for all chains of the structure. An id that failed has no rows.
            missing:      id | removed_by | kind | cause

    `pdockq` and `pdockq2` need pLDDT and `iptm_af` needs the predictor's own ipTM. A PAE file that does not carry them (every `npz`, whose Boltz-layout siblings are not staged; an AlphaFold3 json without `atom_plddts`; a ColabFold json without `iptm`) leaves the column NaN rather than the zero upstream would compute from an absent input.
    """

    TOOL_NAME = "IPSAE"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        repo_dir = folders.get("IPSAE", "")
        if not repo_dir:
            raise ValueError(
                "IPSAE: folders.repositories.IPSAE is not set in the active config; "
                "it names the folder ipsae.py is downloaded into")
        script = f"{repo_dir}/ipsae.py"
        skip = "" if force_reinstall else f"""# Check if already installed
if [ -s "{script}" ]; then
    echo "IPSAE already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        return f"""echo "=== Installing IPSAE ==="
{skip}mkdir -p "{repo_dir}"
rm -f "{script}"
curl -fsSL -o "{script}.part" "{IPSAE_SCRIPT_URL}" || {{ echo "ERROR: download of {IPSAE_SCRIPT_URL} failed"; exit 1; }}
mv "{script}.part" "{script}"

# Verify the download parses as Python and that numpy, its only dependency, imports
if python -c "import ast, sys, numpy; ast.parse(open(sys.argv[1]).read())" "{script}" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== IPSAE installation complete ==="
else
    echo "ERROR: IPSAE verification failed (ipsae.py does not parse or numpy is missing)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    pae_json = Path(lambda self: self.configuration_path("pae.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    runs_folder = Path(lambda self: self.execution_path("runs"))
    ipsae_csv = Path(lambda self: self.table_path("ipsae"))
    pairs_csv = Path(lambda self: self.table_path("ipsae_pairs"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    ipsae_script = Path(lambda self: os.path.join(self.folders["IPSAE"], "ipsae.py"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_ipsae.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 pae: Union[DataStream, StandardizedOutput],
                 binder_chain: str = "A",
                 target_chains: Optional[Union[str, List[str]]] = None,
                 pae_cutoff: float = 10.0,
                 dist_cutoff: float = 10.0,
                 **kwargs):
        self.structures = structures
        self.pae = pae
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        if isinstance(pae, StandardizedOutput):
            self.pae_stream: DataStream = pae.streams.pae
            if self.pae_stream is None:
                raise ValueError(
                    f"pae: this tool output has no 'pae' stream (streams: {list(pae.streams.keys())}). "
                    f"Pass a predictor that exposes its PAE matrices as `.streams.pae`, "
                    f"or a DataStream of PAE files.")
        elif isinstance(pae, DataStream):
            self.pae_stream = pae
        else:
            raise ValueError(f"pae must be DataStream or StandardizedOutput, got {type(pae)}")
        self.binder_chain = binder_chain
        self.target_chains = [target_chains] if isinstance(target_chains, str) else target_chains
        self.pae_cutoff = pae_cutoff
        self.dist_cutoff = dist_cutoff
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if not self.structures_stream.has_only_formats(*STRUCTURE_FORMATS):
            raise ValueError(
                f"structures must be PDB or mmCIF, got format {self.structures_stream.format!r}")
        if not self.pae_stream or len(self.pae_stream) == 0:
            raise ValueError("pae parameter is required and must not be empty")
        if not self.pae_stream.has_only_formats(*PAE_FORMATS):
            raise ValueError(
                f"pae must be npz (Boltz layout) or json (ColabFold / AlphaFold3 layout), "
                f"got format {self.pae_stream.format!r}")
        self._check_chain_id(self.binder_chain, "binder_chain")
        if self.target_chains is not None:
            if not isinstance(self.target_chains, (list, tuple)) or not self.target_chains:
                raise ValueError("target_chains must be a chain id or a non-empty list of chain ids")
            for chain in self.target_chains:
                self._check_chain_id(chain, "target_chains")
            if self.binder_chain in self.target_chains:
                raise ValueError(f"target_chains must not contain binder_chain {self.binder_chain!r}")
            if len(set(self.target_chains)) != len(self.target_chains):
                raise ValueError(f"target_chains lists a chain more than once: {list(self.target_chains)}")
        for name in ("pae_cutoff", "dist_cutoff"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
                raise ValueError(f"{name} must be a positive number (A), got {value!r}")

    @staticmethod
    def _check_chain_id(chain, parameter):
        if not isinstance(chain, str) or not chain or any(c.isspace() for c in chain):
            raise ValueError(f"{parameter} entries must be non-empty chain ids without whitespace, got {chain!r}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"BINDER CHAIN: {self.binder_chain}")
        lines.append(f"TARGET CHAINS: {', '.join(self.target_chains) if self.target_chains else 'all other chains'}")
        lines.append(f"PAE CUTOFF: {self.pae_cutoff} A")
        lines.append(f"DIST CUTOFF: {self.dist_cutoff} A")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        self.pae_stream.save_json(self.pae_json)
        self._write_config_yaml()
        script = "#!/bin/bash\n"
        script += "# IPSAE script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Running IPSAE on {len(self.structures_stream)} structure(s)"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --structures-json "{self.structures_json}" \\
    --pae-json "{self.pae_json}" \\
    --ipsae-script "{self.ipsae_script}" \\
    --runs-dir "{self.runs_folder}" \\
    --ipsae-csv "{self.ipsae_csv}" \\
    --pairs-csv "{self.pairs_csv}" \\
    --local-missing-csv "{self.local_missing_csv}"
"""
        script += self.generate_missing_propagation(
            self.structures, self.pae, local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        import yaml
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump({
                "binder_chain": self.binder_chain,
                "target_chains": list(self.target_chains) if self.target_chains else None,
                "pae_cutoff": float(self.pae_cutoff),
                "dist_cutoff": float(self.dist_cutoff),
            }, f, sort_keys=False)

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "ipsae": TableInfo(
                name="ipsae",
                path=self.ipsae_csv,
                columns=SUMMARY_COLUMNS,
                description="ipSAE of the binder chain against its best-scoring target chain, one row per id",
            ),
            "ipsae_pairs": TableInfo(
                name="ipsae_pairs",
                path=self.pairs_csv,
                columns=PAIR_COLUMNS,
                description="Every chain-pair row upstream ipsae.py printed, asymmetric and max, per id",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
