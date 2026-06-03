# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""PLIP tool: Protein-Ligand Interaction Profiler.

Reads protein-ligand complex PDBs and reports detected non-covalent interactions
(hydrogen bonds, hydrophobic contacts, pi-stacking, salt bridges, halogen bonds,
water bridges, metal complexes).
"""

import os
from typing import Dict, List, Any, Union, Optional

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .input_standardization import resolve_basic_input
    from .ligand import Ligand
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from input_standardization import resolve_basic_input
    from ligand import Ligand


_MODES = ("ligand", "peptide", "intra")


class PLIP(BaseConfig):
    """
    PLIP: detect non-covalent interactions in protein complexes.

    Three analysis modes, selected by `mode`:
        "ligand" (default): protein-ligand interactions. Each complex must
                    contain at least one ligand in HETATM records.
        "peptide": protein-peptide / protein-protein interactions. The
                    chain(s) named in `chains` (deposited as ATOM records) are
                    treated as the peptide/partner "ligand"; PLIP profiles
                    their contacts with the rest of the structure (PLIP's
                    --peptides / --inter).
        "intra":  intra-chain interactions within the single chain named in
                    `chains` (PLIP's --intra).

    Inputs:
        structures: PDB complex structures.
        ligand: Optional compounds stream (Ligand(code="MK1"), or a multi-row
                    stream for several codes) restricting profiling to those
                    residue codes. The code(s) are read from the stream's
                    `code` column at runtime. Omit to profile all non-standard
                    residues. Only valid in mode="ligand".
        mode: "ligand" | "peptide" | "intra" (default "ligand").
        chains: Chain ID(s) to treat as the peptide/partner (mode="peptide")
                    or the single chain to profile (mode="intra"). A string or
                    list of strings. Required when mode != "ligand"; "intra"
                    takes exactly one chain.
        generate_pse: If True, also produce one PyMOL session (.pse) per
                    binding site detected (default: True). Requires PyMOL
                    in the env (shipped via pymol-open-source).

    Outputs:
        Streams:
            sessions: one <id>.pse per complex (the first binding site's
                      PyMOL session; multi-site complexes list all session
                      files in the summary table's `session_files` column).
                      Only present when generate_pse=True.
        Tables:
            interactions: id | ligand | interaction_type | residue | chain | resnum | distance | details
            summary:      id | n_hbonds | n_hydrophobic | n_pi_stacking | n_salt_bridges | n_halogen | n_water_bridges | session_files
    """

    TOOL_NAME = "PLIP"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        image = folders.get(f"container:{cls.TOOL_NAME}", "")
        if image:
            # Container mode: nothing to install. The image is expected to be
            # provided by the site; we just verify it exists and is executable.
            return f"""echo "=== PLIP (container mode) ==="
if [ ! -f "{image}" ]; then
    echo "ERROR: PLIP container image not found at {image}"
    exit 1
fi
touch "$INSTALL_SUCCESS"
echo "=== PLIP container ready ==="
"""

        env_check = cls._env_exists_check("plip", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "PLIP already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("plip", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("plip", env_manager, biopipelines)
        return f"""echo "=== Installing PLIP ==="
{skip}{remove_block}
{env_block}

if {env_manager} run -n plip python -c "import plip" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== PLIP installation complete ==="
else
    echo "ERROR: PLIP verification failed (cannot import plip)"
    exit 1
fi
"""

    structures_json = Path(lambda self: self.configuration_path("structures.json"))
    ligand_json = Path(lambda self: self.configuration_path("input_ligand.json"))
    sessions_map = Path(lambda self: self.stream_map_path("sessions"))
    interactions_csv = Path(lambda self: self.table_path("interactions"))
    summary_csv = Path(lambda self: self.table_path("summary"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_plip.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 ligand: Optional[Union[str, DataStream, StandardizedOutput]] = None,
                 mode: str = "ligand",
                 chains: Optional[Union[str, List[str]]] = None,
                 generate_pse: bool = True,
                 **kwargs):
        self.structures = structures
        self.ligand = ligand
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")
        # Optional ligand filter — a compounds stream; codes resolved at runtime.
        # A bare string is shorthand for an internal Ligand(code=...).
        self.ligand_stream: Optional[DataStream] = resolve_basic_input(
            ligand, Ligand, "compounds", "code")
        self.mode = mode
        self.chains: List[str] = (
            [] if chains is None else [chains] if isinstance(chains, str) else list(chains)
        )
        self.generate_pse = bool(generate_pse)
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if self.mode not in _MODES:
            raise ValueError(f"mode must be one of {_MODES}, got {self.mode!r}")
        if self.mode == "ligand":
            if self.chains:
                raise ValueError("chains is only valid when mode is 'peptide' or 'intra'")
        else:
            if not self.chains:
                raise ValueError(f"chains is required when mode is {self.mode!r}")
            if self.mode == "intra" and len(self.chains) != 1:
                raise ValueError("mode='intra' takes exactly one chain")
            if self.ligand_stream is not None:
                raise ValueError("ligand filter is only valid in mode='ligand'")
            for c in self.chains:
                _validate_freeform_string("chains", c)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"STRUCTURES: {len(self.structures_stream)} complexes")
        lines.append(f"MODE: {self.mode}")
        if self.chains:
            lines.append(f"CHAINS: {' '.join(self.chains)}")
        if self.ligand_stream:
            lines.append("LIGAND FILTER: (codes resolved from compounds stream at runtime)")
        lines.append(f"GENERATE PSE: {self.generate_pse}")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.structures_stream.save_json(self.structures_json)
        script = "#!/bin/bash\n"
        script += "# PLIP interaction profiling script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        # The ligand-code filter is read from the compounds stream inside the
        # pipe script; omitted entirely when no ligand stream is given.
        filter_arg = ""
        if self.ligand_stream is not None:
            self.ligand_stream.save_json(self.ligand_json)
            filter_arg = f' --ligand-json "{self.ligand_json}"'
        mode_arg = f' --mode {self.mode}'
        chains_arg = f' --chains {" ".join(self.chains)}' if self.chains else ""
        pse_arg = " --generate-pse" if self.generate_pse else ""
        pse_args = (
            f' --sessions-dir "{self.stream_folder("sessions")}"'
            f' --sessions-map-csv "{self.sessions_map}"'
            if self.generate_pse else ""
        )
        # When a container image is configured for PLIP, run the CLI via the
        # container and parse report.xml; otherwise call the Python API from
        # the active conda env. PLIP's image runscript already invokes
        # plipcmd.py as the entrypoint, so swap `exec` -> `run` and drop the
        # GPU flag (PLIP is CPU-only and --nv fails when nvidia-container-cli
        # is missing on CPU partitions).
        prefix_arg = ""
        if self.uses_container():
            prefix = (self.container_prefix()
                      .replace("exec", "run", 1)
                      .replace(" --nv ", " ").strip())
            prefix_arg = f' --container-prefix "{prefix}"'
        upstream_missing_flag = self.upstream_missing_flag(self.structures, self.ligand)
        script += f"""echo "Running PLIP on {len(self.structures_stream)} complex(es)"
python "{self.helper_py}" \\
    --structures-json "{self.structures_json}" \\
    --interactions-csv "{self.interactions_csv}" \\
    --summary-csv "{self.summary_csv}" \\
    --missing-csv "{self.missing_csv}"{mode_arg}{chains_arg}{filter_arg}{pse_arg}{pse_args}{prefix_arg}{upstream_missing_flag}
"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "interactions": TableInfo(
                name="interactions",
                path=self.interactions_csv,
                columns=["id", "ligand", "interaction_type", "residue", "chain",
                         "resnum", "distance", "details"],
                description="PLIP detected non-covalent interactions",
            ),
            "summary": TableInfo(
                name="summary",
                path=self.summary_csv,
                columns=["id", "n_hbonds", "n_hydrophobic", "n_pi_stacking",
                         "n_salt_bridges", "n_halogen", "n_water_bridges", "session_files"],
                description="PLIP per-complex interaction counts",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="IDs removed (upstream or local failure) with removal reason",
            ),
        }
        result: Dict[str, Any] = {
            "tables": tables,
            "output_folder": self.output_folder,
        }
        if self.generate_pse:
            result["sessions"] = DataStream(
                name="sessions",
                ids=self.structures_stream.ids,
                files=[self.stream_path("sessions", "<id>.pse")],
                map_table=self.sessions_map,
                format="pse",
            )
        return result
