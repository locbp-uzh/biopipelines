# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Ligand tool for fetching small molecule ligands from RCSB PDB, PubChem, or SMILES strings.

Downloads SDF files and converts them to PDB format with proper atom numbering.
Supports lookup by CCD code (RCSB), compound name, CID, or CAS number (PubChem).
Also supports direct SMILES input for custom molecules.
Fetches ligands with priority-based lookup: local_folder -> Ligands/ -> RCSB/PubChem download.
"""

import os
import re
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


class Ligand(BaseConfig):
    """
    Pipeline tool for fetching small molecule ligands from RCSB PDB, PubChem, or SMILES strings.

    Downloads SDF files and converts to PDB format with proper atom numbering.
    Implements priority-based lookup: checks local_folder (if provided), then
    Ligands/ folder, then downloads from RCSB or PubChem based on lookup type.
    Also supports direct SMILES input for custom molecules.
    """

    TOOL_NAME = "Ligand"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Ligand ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== Ligand ready ==="
"""

    # Lazy path descriptors
    compounds_csv = Path(lambda self: os.path.join(self.output_folder, "compounds.csv"))
    failed_csv = Path(lambda self: os.path.join(self.output_folder, "failed_downloads.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "fetch_config.json"))
    ligand_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_ligand.py"))

    def __init__(self,
                 lookup: Optional[Union[str, List[str]]] = None,
                 ids: Optional[Union[str, List[str]]] = None,
                 codes: Optional[Union[str, List[str]]] = None,
                 source: Optional[str] = None,
                 local_folder: Optional[str] = None,
                 output_format: str = "pdb",
                 smiles: Optional[Union[str, List[str]]] = None,
                 **kwargs):
        """
        Initialize Ligand tool.

        Args:
            lookup: Lookup value(s) for fetching. Can be:
                    - RCSB CCD codes: "ATP", "GDP", "HEM"
                    - PubChem CID: "2244"
                    - PubChem CAS: "50-78-2"
                    - PubChem name: "aspirin", "caffeine"
                    Can be None if using smiles parameter instead.
            ids: Output identifier(s) for filenames (e.g., "my_ligand" -> my_ligand.pdb).
                 If not provided, defaults to lookup values (for lookup) or "smiles_1", "smiles_2", etc. (for smiles).
            codes: 3-letter PDB residue code(s) to use in the PDB file (e.g., "LIG").
                   If not provided, defaults to lookup[:3].upper() (for lookup) or "LIG" (for smiles).
            source: Force source ("rcsb" or "pubchem"). If None, auto-detects:
                    - 1-3 uppercase alphanumeric -> RCSB (CCD)
                    - Purely numeric -> PubChem (CID)
                    - XX-XX-X format -> PubChem (CAS)
                    - Otherwise -> PubChem (name)
                    Ignored when using smiles parameter.
            local_folder: Custom local folder to check first (before Ligands/). Default: None
            output_format: Output format - "pdb" or "cif". CIF includes explicit bond orders
                          which is recommended for tools like RFdiffusion3. Default: "pdb"
            smiles: SMILES string(s) for direct molecule input. Bypasses lookup entirely.
                    Molecules are converted to 3D structures using RDKit.
            **kwargs: Additional parameters

        Fetch Priority (for lookup):
            For each ligand, searches in order:
            1. local_folder (if parameter provided)
            2. ./Ligands/ folder in repository
            3. Download from RCSB or PubChem (saved to both Ligands/ and output folder)

        Examples:
            # Simplest: just lookup (ids and codes derived automatically)
            lig = Ligand("ATP")  # ids=["ATP"], codes=["ATP"]

            # Multiple ligands with just lookup
            lig = Ligand(["ATP", "GDP"])  # ids=["ATP", "GDP"], codes=["ATP", "GDP"]

            # With custom id
            lig = Ligand("ATP", ids="atp_ligand")

            # PubChem ligand by compound name
            lig = Ligand("aspirin", ids="aspirin_lig", codes="ASP")

            # PubChem by CID
            lig = Ligand("2157", ids="caffeine", codes="CAF")

            # PubChem by CAS number
            lig = Ligand("15687-27-1", ids="ibuprofen", codes="IBU")

            # Multiple ligands with explicit ids and codes
            lig = Ligand(
                ["ATP", "aspirin"],
                ids=["lig1", "lig2"],
                codes=["ATP", "LIG"]
            )

            # Force source
            lig = Ligand("ATP", source="pubchem")

            # Direct SMILES input (single)
            lig = Ligand(smiles="CCO")  # ids=["smiles_1"], codes=["LIG"]

            # Direct SMILES input (multiple)
            lig = Ligand(smiles=["CCO", "CC(=O)O"])  # ids=["smiles_1", "smiles_2"], codes=["LIG", "LIG"]

            # SMILES with custom ids and codes
            lig = Ligand(smiles="CCO", ids="ethanol", codes="ETH")
        """
        # Handle smiles input
        if smiles is not None:
            if isinstance(smiles, str):
                self.smiles_values = [smiles]
            else:
                self.smiles_values = list(smiles)
        else:
            self.smiles_values = []

        # Handle lookup
        if lookup is not None:
            if isinstance(lookup, str):
                self.lookup_values = [lookup]
            else:
                self.lookup_values = list(lookup)
        else:
            self.lookup_values = []

        # Validate: must have at least one of lookup or smiles
        if not self.lookup_values and not self.smiles_values:
            raise ValueError("Must provide at least one of 'lookup' or 'smiles'")

        # Total number of ligands
        total_count = len(self.lookup_values) + len(self.smiles_values)

        # Handle ids - default based on input type
        if ids is not None:
            if isinstance(ids, str):
                self.custom_ids = [ids]
            else:
                self.custom_ids = list(ids)
        else:
            # Default ids: lookup values for lookup, "smiles_N" for smiles
            self.custom_ids = self.lookup_values.copy()
            for i in range(len(self.smiles_values)):
                self.custom_ids.append(f"smiles_{i + 1}")

        # Handle codes - default based on input type
        if codes is not None:
            if isinstance(codes, str):
                self.residue_codes = [codes.upper()]
            else:
                self.residue_codes = [c.upper() for c in codes]
        else:
            # Default codes: lookup[:3].upper() for lookup, "LIG" for smiles
            self.residue_codes = [lv.upper()[:3] for lv in self.lookup_values]
            for _ in range(len(self.smiles_values)):
                self.residue_codes.append("LIG")

        # Validate lengths
        if len(self.custom_ids) != total_count:
            raise ValueError(f"Length mismatch: ids has {len(self.custom_ids)} items but total ligands is {total_count}")
        if len(self.residue_codes) != total_count:
            raise ValueError(f"Length mismatch: codes has {len(self.residue_codes)} items but total ligands is {total_count}")

        # Validate source
        if source is not None and source not in ["rcsb", "pubchem"]:
            raise ValueError(f"Invalid source: {source}. Must be 'rcsb', 'pubchem', or None")
        self.source = source

        self.local_folder = local_folder

        # Validate and store output format
        if output_format not in ["pdb", "cif"]:
            raise ValueError(f"Invalid output_format: {output_format}. Must be 'pdb' or 'cif'")
        self.output_format = output_format

        # Warn about non-standard residue codes (standard PDB is 1-3 alphanumeric)
        for code in self.residue_codes:
            if not code:
                raise ValueError("Residue code cannot be empty")
            if len(code) > 3 or not code.isalnum():
                print(f"  Warning: Non-standard residue code '{code}'. Standard PDB uses 1-3 alphanumeric characters.")

        # Initialize base class
        super().__init__(**kwargs)

    def _detect_lookup_type(self, lookup: str) -> str:
        """
        Detect the type of lookup value.

        Returns: "ccd" (RCSB), "cid" (PubChem), "cas" (PubChem), or "name" (PubChem)
        """
        # CCD codes: 1-3 uppercase alphanumeric characters
        if re.match(r'^[A-Z0-9]{1,3}$', lookup.upper()) and not lookup.isdigit():
            return "ccd"

        # PubChem CID: purely numeric
        if lookup.isdigit():
            return "cid"

        # CAS number: XX-XX-X format (digits-digits-digit)
        if re.match(r'^\d+-\d+-\d$', lookup):
            return "cas"

        # Default: compound name
        return "name"

    def validate_params(self):
        """Validate Ligand parameters."""
        if not self.custom_ids:
            raise ValueError("ids cannot be empty")

        if not self.lookup_values and not self.smiles_values:
            raise ValueError("Must have at least one of lookup or smiles")

        if not self.residue_codes:
            raise ValueError("codes cannot be empty")

        total_count = len(self.lookup_values) + len(self.smiles_values)
        if len(self.custom_ids) != total_count:
            raise ValueError(f"ids length ({len(self.custom_ids)}) must match total ligands ({total_count})")

        if len(self.residue_codes) != total_count:
            raise ValueError(f"codes length ({len(self.residue_codes)}) must match total ligands ({total_count})")

        for code in self.residue_codes:
            if not code:
                raise ValueError("Residue code cannot be empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters and check for local files."""
        self.folders = pipeline_folders

        # Check which files exist locally and which will need to be downloaded
        repo_ligands_folder = pipeline_folders.get('Ligands', '')
        self.found_locally = []
        self.needs_download = []

        # Process lookup values
        for lookup, custom_id in zip(self.lookup_values, self.custom_ids[:len(self.lookup_values)]):
            found = False
            local_path = None

            # Check local_folder if specified
            if self.local_folder:
                local_path = os.path.join(self.local_folder, f"{lookup}.pdb")
                if os.path.exists(local_path):
                    self.found_locally.append((lookup, local_path))
                    found = True
                    print(f"  Found ligand {lookup} locally: {local_path}")

            # Check Ligands/ folder
            if not found and repo_ligands_folder:
                local_path = os.path.join(repo_ligands_folder, f"{lookup}.pdb")
                if os.path.exists(local_path):
                    self.found_locally.append((lookup, local_path))
                    found = True
                    print(f"  Found ligand {lookup} locally: {local_path}")

            if not found:
                # Determine source for display
                lookup_type = self._detect_lookup_type(lookup)
                effective_source = self.source if self.source else ("rcsb" if lookup_type == "ccd" else "pubchem")

                self.needs_download.append(lookup)
                print(f"  Ligand {lookup} not found locally, will download from {effective_source}")

        # Process SMILES values (these are always "generated" at runtime)
        smiles_start_idx = len(self.lookup_values)
        for i, smiles in enumerate(self.smiles_values):
            custom_id = self.custom_ids[smiles_start_idx + i]
            smiles_preview = smiles[:30] + "..." if len(smiles) > 30 else smiles
            print(f"  SMILES ligand {custom_id}: {smiles_preview} (will be generated with RDKit)")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"IDS: {', '.join(self.custom_ids)} ({len(self.custom_ids)} ligands)",
            f"CODES: {', '.join(self.residue_codes)}",
            f"FORMAT: {self.output_format.upper()}",
        ])

        if self.lookup_values:
            config_lines.append(f"LOOKUP: {', '.join(self.lookup_values)}")
            config_lines.append(f"SOURCE: {self.source if self.source else 'auto-detect'}")
            config_lines.append(f"LOCAL_FOLDER: {self.local_folder if self.local_folder else 'None (uses Ligands/)'}")

        if self.smiles_values:
            smiles_preview = [s[:20] + "..." if len(s) > 20 else s for s in self.smiles_values]
            config_lines.append(f"SMILES: {', '.join(smiles_preview)} ({len(self.smiles_values)} molecules)")

        # Add status of files found/not found
        if hasattr(self, 'found_locally') and self.found_locally:
            config_lines.append(f"FOUND_LOCALLY: {', '.join([lookup for lookup, _ in self.found_locally])}")
        if hasattr(self, 'needs_download') and self.needs_download:
            config_lines.append(f"NEEDS_DOWNLOAD: {', '.join(self.needs_download)}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate script to fetch ligands from local folders or RCSB/PubChem."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# Ligand execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_ligand()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_ligand(self) -> str:
        """Generate the ligand fetching part of the script."""
        import json

        repo_ligands_folder = self.folders['Ligands']

        config_data = {
            "custom_ids": self.custom_ids,
            "residue_codes": self.residue_codes,
            "lookup_values": self.lookup_values,
            "smiles_values": self.smiles_values,
            "source": self.source,
            "local_folder": self.local_folder,
            "output_format": self.output_format,
            "repo_ligands_folder": repo_ligands_folder,
            "output_folder": self.output_folder,
            "compounds_table": self.compounds_csv,
            "failed_table": self.failed_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        total_count = len(self.lookup_values) + len(self.smiles_values)
        lookup_info = f"Lookup: {', '.join(self.lookup_values)}" if self.lookup_values else ""
        smiles_info = f"SMILES: {len(self.smiles_values)} molecule(s)" if self.smiles_values else ""

        return f"""echo "Processing {total_count} ligands"
echo "Custom IDs: {', '.join(self.custom_ids)}"
echo "Residue codes: {', '.join(self.residue_codes)}"
{f'echo "{lookup_info}"' if lookup_info else ''}
{f'echo "{smiles_info}"' if smiles_info else ''}
echo "Output folder: {self.output_folder}"

python "{self.ligand_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after ligand fetching."""
        # Generate structure file paths using custom IDs
        ext = self.output_format  # "pdb" or "cif"
        structure_files = [os.path.join(self.output_folder, f"{custom_id}.{ext}")
                          for custom_id in self.custom_ids]

        total_count = len(self.lookup_values) + len(self.smiles_values)

        tables = {
            "compounds": TableInfo(
                name="compounds",
                path=self.compounds_csv,
                columns=["id", "format", "code", "lookup", "source", "ccd", "cid", "cas", "smiles", "name", "formula", "file_path"],
                description="Successfully fetched/generated ligand files with metadata",
                count=total_count
            ),
            "failed": TableInfo(
                name="failed",
                path=self.failed_csv,
                columns=["lookup", "error_message", "source", "attempted_path"],
                description="Failed ligand fetches with error details",
                count="variable"
            )
        }

        # Create DataStreams
        structures = DataStream(
            name="structures",
            ids=self.custom_ids.copy(),
            files=structure_files,
            format=self.output_format
        )

        compounds = DataStream(
            name="compounds",
            ids=self.custom_ids.copy(),
            files=[self.compounds_csv],
            map_table=self.compounds_csv,
            format="csv"
        )

        return {
            "structures": structures,
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": compounds,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "custom_ids": self.custom_ids,
                "residue_codes": self.residue_codes,
                "lookup_values": self.lookup_values,
                "smiles_values": self.smiles_values,
                "source": self.source,
                "local_folder": self.local_folder,
                "output_format": self.output_format
            }
        })
        return base_dict
