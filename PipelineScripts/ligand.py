"""
Ligand tool for fetching small molecule ligands from RCSB PDB or PubChem.

Downloads SDF files and converts them to PDB format with proper atom numbering.
Supports lookup by CCD code (RCSB), compound name, CID, or CAS number (PubChem).
Fetches ligands with priority-based lookup: local_folder -> Ligands/ -> RCSB/PubChem download.
"""

import os
import re
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class Ligand(BaseConfig):
    """
    Pipeline tool for fetching small molecule ligands from RCSB PDB or PubChem.

    Downloads SDF files and converts to PDB format with proper atom numbering.
    Implements priority-based lookup: checks local_folder (if provided), then
    Ligands/ folder, then downloads from RCSB or PubChem based on lookup type.
    """

    # Tool identification
    TOOL_NAME = "Ligand"
    

    def __init__(self,
                 lookup: Union[str, List[str]],
                 ids: Optional[Union[str, List[str]]] = None,
                 codes: Optional[Union[str, List[str]]] = None,
                 source: Optional[str] = None,
                 local_folder: Optional[str] = None,
                 output_format: str = "pdb",
                 **kwargs):
        """
        Initialize Ligand tool.

        Args:
            lookup: Lookup value(s) for fetching. Can be:
                    - RCSB CCD codes: "ATP", "GDP", "HEM"
                    - PubChem CID: "2244"
                    - PubChem CAS: "50-78-2"
                    - PubChem name: "aspirin", "caffeine"
            ids: Output identifier(s) for filenames (e.g., "my_ligand" -> my_ligand.pdb).
                 If not provided, defaults to lookup values.
            codes: 3-letter PDB residue code(s) to use in the PDB file (e.g., "LIG").
                   If not provided, defaults to lookup[:3].upper().
            source: Force source ("rcsb" or "pubchem"). If None, auto-detects:
                    - 1-3 uppercase alphanumeric -> RCSB (CCD)
                    - Purely numeric -> PubChem (CID)
                    - XX-XX-X format -> PubChem (CAS)
                    - Otherwise -> PubChem (name)
            local_folder: Custom local folder to check first (before Ligands/). Default: None
            output_format: Output format - "pdb" or "cif". CIF includes explicit bond orders
                          which is recommended for tools like RFdiffusion3. Default: "pdb"
            **kwargs: Additional parameters

        Fetch Priority:
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
        """
        # Handle lookup (required)
        if isinstance(lookup, str):
            self.lookup_values = [lookup]
        else:
            self.lookup_values = list(lookup)

        # Handle ids - default to lookup if not provided
        if ids is not None:
            if isinstance(ids, str):
                self.custom_ids = [ids]
            else:
                self.custom_ids = list(ids)
        else:
            # Default ids to lookup values
            self.custom_ids = self.lookup_values.copy()

        # Handle codes - default to lookup[:3].upper() if not provided
        if codes is not None:
            if isinstance(codes, str):
                self.residue_codes = [codes.upper()]
            else:
                self.residue_codes = [c.upper() for c in codes]
        else:
            # Default codes to lookup values (truncated to 3 chars, uppercased)
            self.residue_codes = [lv.upper()[:3] for lv in self.lookup_values]

        # Validate lengths
        if len(self.custom_ids) != len(self.lookup_values):
            raise ValueError(f"Length mismatch: ids has {len(self.custom_ids)} items but lookup has {len(self.lookup_values)} items")
        if len(self.custom_ids) != len(self.residue_codes):
            raise ValueError(f"Length mismatch: ids has {len(self.custom_ids)} items but codes has {len(self.residue_codes)} items")

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

        if not self.lookup_values:
            raise ValueError("lookup cannot be empty")

        if not self.residue_codes:
            raise ValueError("codes cannot be empty")

        if len(self.custom_ids) != len(self.lookup_values):
            raise ValueError("ids and lookup must have same length")

        if len(self.custom_ids) != len(self.residue_codes):
            raise ValueError("ids and codes must have same length")

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

        for lookup, custom_id in zip(self.lookup_values, self.custom_ids):
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

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"IDS: {', '.join(self.custom_ids)} ({len(self.custom_ids)} ligands)",
            f"LOOKUP: {', '.join(self.lookup_values)}",
            f"CODES: {', '.join(self.residue_codes)}",
            f"SOURCE: {self.source if self.source else 'auto-detect'}",
            f"FORMAT: {self.output_format.upper()}",
            f"LOCAL_FOLDER: {self.local_folder if self.local_folder else 'None (uses Ligands/)'}"
        ])

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
        compounds_table = os.path.join(self.output_folder, "compounds.csv")
        failed_table = os.path.join(self.output_folder, "failed_downloads.csv")
        repo_ligands_folder = self.folders['Ligands']

        config_file = os.path.join(self.output_folder, "fetch_config.json")
        config_data = {
            "custom_ids": self.custom_ids,
            "residue_codes": self.residue_codes,
            "lookup_values": self.lookup_values,
            "source": self.source,
            "local_folder": self.local_folder,
            "output_format": self.output_format,
            "repo_ligands_folder": repo_ligands_folder,
            "output_folder": self.output_folder,
            "compounds_table": compounds_table,
            "failed_table": failed_table
        }

        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Fetching {len(self.lookup_values)} ligands"
echo "Lookup values: {', '.join(self.lookup_values)}"
echo "Custom IDs: {', '.join(self.custom_ids)}"
echo "Residue codes: {', '.join(self.residue_codes)}"
echo "Source: {self.source if self.source else 'auto-detect'}"
echo "Priority: {'local_folder -> ' if self.local_folder else ''}Ligands/ -> download"
echo "Output folder: {self.output_folder}"

python "{os.path.join(self.folders['HelpScripts'], 'pipe_ligand.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Successfully fetched ligands"
    echo "Compounds table: {compounds_table}"
    if [ -f "{failed_table}" ]; then
        echo "Failed downloads logged: {failed_table}"
    fi
else
    echo "Error: Failed to fetch ligands"
    exit 1
fi

"""

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after ligand fetching.

        Returns:
            Dictionary with output file paths and table information
        """
        # Generate structure file paths using custom IDs
        ext = self.output_format  # "pdb" or "cif"
        structure_files = []
        for custom_id in self.custom_ids:
            structure_files.append(os.path.join(self.output_folder, f"{custom_id}.{ext}"))

        # Structure IDs are the custom IDs
        structure_ids = self.custom_ids.copy()

        # Output tables
        compounds_csv = os.path.join(self.output_folder, "compounds.csv")
        failed_csv = os.path.join(self.output_folder, "failed_downloads.csv")

        # Define tables that will be created
        tables = {
            "compounds": TableInfo(
                name="compounds",
                path=compounds_csv,
                columns=["id", "format", "code", "lookup", "source", "ccd", "cid", "cas", "smiles", "name", "formula", "file_path"],
                description="Successfully fetched ligand files with metadata",
                count=len(self.lookup_values)
            ),
            "failed": TableInfo(
                name="failed",
                path=failed_csv,
                columns=["lookup", "error_message", "source", "attempted_path"],
                description="Failed ligand fetches with error details",
                count="variable"
            )
        }

        return {
            "structures": structure_files,
            "structure_ids": structure_ids,
            "compounds": [compounds_csv],
            "compound_ids": structure_ids,  # Same as structure_ids
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
                "source": self.source,
                "local_folder": self.local_folder,
                "output_format": self.output_format
            }
        })
        return base_dict
