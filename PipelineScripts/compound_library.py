"""
CompoundLibrary configuration for generating and processing ligand libraries.

Handles dictionary-based SMILES library generation, CSV output generation,
and optional covalent ligand CCD/PKL file preparation for Boltz2.
"""

import os
import json
import csv
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


class CompoundLibrary(BaseConfig):
    """
    CompoundLibrary configuration for dictionary-based SMILES library generation.

    Expands dictionaries with substitution keys into complete compound libraries,
    generates CSV files with standardized format, and optionally prepares
    covalent ligand files for Boltz2.
    """

    TOOL_NAME = "CompoundLibrary"

    # Lazy path descriptors
    compounds_csv = Path(lambda self: os.path.join(self.output_folder, "compounds.csv"))
    compound_properties_csv = Path(lambda self: os.path.join(self.output_folder, "compound_properties.csv"))
    summary_file = Path(lambda self: os.path.join(self.output_folder, "summary.txt"))
    library_dict_json = Path(lambda self: os.path.join(self.output_folder, "library_dict.json"))
    covalent_folder = Path(lambda self: os.path.join(self.output_folder, "covalent_library") if self.covalent else None)
    covalent_compounds_csv = Path(lambda self: os.path.join(self.output_folder, "covalent_library", "compounds.csv") if self.covalent else None)
    compound_expansion_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_compound_library.py"))
    smiles_properties_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_smiles_properties.py"))
    covalent_generation_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_compound_library.py"))

    def __init__(self,
                 library: Union[str, Dict[str, Union[str, List[str]]]],
                 primary_key: Optional[str] = None,
                 covalent: bool = False,
                 validate_smiles: bool = True,
                 conformer_method: str = "UFF",
                 **kwargs):
        """
        Initialize CompoundLibrary configuration.

        Args:
            library: Dictionary with expansion keys or path to existing CSV library
            primary_key: Root key for expansion when library is a dictionary
            covalent: Generate CCD/PKL files for covalent ligand binding (calls runtime script)
            validate_smiles: Validate SMILES strings during expansion
            conformer_method: Method for conformer generation ("UFF", "OpenFF", "DFT")
            **kwargs: Additional parameters
        """
        # Store CompoundLibrary-specific parameters
        self.library = library
        self.primary_key = primary_key
        self.covalent = covalent
        self.validate_smiles = validate_smiles
        self.conformer_method = conformer_method

        # Track library source type
        self.library_dict = None
        self.library_csv = None
        self.expanded_compounds = []
        self.compound_ids = []

        # Initialize base class
        super().__init__(**kwargs)

        # For dictionary libraries, do expansion immediately to get compound IDs
        if isinstance(self.library, dict):
            self.library_dict = self.library
            self._expand_library()

    def validate_params(self):
        """Validate CompoundLibrary-specific parameters."""
        if not self.library:
            raise ValueError("library parameter is required")

        # Validate library format
        if isinstance(self.library, dict):
            # Dictionary-based library
            if self.primary_key and self.primary_key not in self.library:
                raise ValueError(f"primary_key '{self.primary_key}' not found in library dictionary")
        elif isinstance(self.library, str):
            # CSV file path
            if not self.library.endswith('.csv'):
                raise ValueError("library file must have .csv extension")
        else:
            raise ValueError("library must be a dictionary or CSV file path")

        # Validate conformer method for covalent ligands
        if self.covalent:
            valid_methods = ["UFF", "OpenFF", "DFT"]
            if self.conformer_method not in valid_methods:
                raise ValueError(f"conformer_method must be one of: {valid_methods}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input library sources."""
        self.folders = pipeline_folders

        if isinstance(self.library, str):
            # CSV file - check if it exists
            if os.path.exists(self.library):
                self.library_csv = self.library
            else:
                # Try in project directory
                project_path = os.path.join(pipeline_folders["biopipelines"], self.library)
                if os.path.exists(project_path):
                    self.library_csv = project_path
                else:
                    raise ValueError(f"Library CSV file not found: {self.library}")
        elif isinstance(self.library, dict):
            # Dictionary-based library (expansion already done in __init__)
            if not self.library_dict:
                self.library_dict = self.library
                self._expand_library()

    def _expand_library(self):
        """Expand dictionary library into individual compounds and generate IDs."""
        if not self.library_dict:
            return

        # Use the same expansion logic as boltz_compound_library.py but with <key> format
        library = self.library_dict.copy()
        library_keys = list(library.keys())
        primary_lib_key = self.primary_key

        # Find the primary key in the base configuration
        if primary_lib_key and primary_lib_key not in library:
            raise ValueError(f"Primary key '{primary_lib_key}' not found in library")

        # Expand each base compound from the primary key
        if primary_lib_key:
            final_compounds = []
            primary_value = library[primary_lib_key]

            # Handle both string and list values for primary key
            if isinstance(primary_value, str):
                # Single SMILES string
                final_compounds.append({'smiles': primary_value, 'branching': {}})
            elif isinstance(primary_value, list):
                # List of SMILES strings
                for base in primary_value:
                    final_compounds.append({'smiles': base, 'branching': {}})
            else:
                raise ValueError(f"Primary key '{primary_lib_key}' must be a string or list of strings")

            no_new_branching = False
            while not no_new_branching:
                no_new_branching = True
                updated_compounds = []

                # Process every compound in current list
                for compound in final_compounds:
                    key_found = False
                    # Check for every library key in the current compound's SMILES
                    for key in library_keys:
                        # Use <key> format instead of *key*
                        key_pattern = f"<{key}>"
                        if key_pattern in compound['smiles']:
                            key_found = True
                            no_new_branching = False
                            # For every possible substitution for the key, create a new compound
                            key_options = library[key]
                            # Handle both string and list values for expansion keys
                            if isinstance(key_options, str):
                                key_options = [key_options]

                            for option in key_options:
                                new_smiles = compound['smiles'].replace(key_pattern, option, 1)
                                new_branching = compound['branching'].copy()
                                new_branching[key] = option
                                updated_compounds.append({'smiles': new_smiles, 'branching': new_branching})
                            # Process one key per compound per iteration
                            break

                    if not key_found:
                        updated_compounds.append(compound)

                final_compounds = updated_compounds

            # Generate compound IDs based on the pattern from boltz_compound_library.py
            num_compounds = len(final_compounds)
            characters = 4
            if num_compounds > 9: characters = 3
            if num_compounds > 99: characters = 2
            if num_compounds > 999: characters = 1
            if num_compounds > 99999: characters = 0

            compound_ids = []
            for u_l_n in range(num_compounds):
                u_l_n_str = str(u_l_n)
                n0 = 5 - characters - len(u_l_n_str)
                zeros_str = '0' * n0
                compound_name = primary_lib_key if num_compounds == 1 else primary_lib_key[:characters] + zeros_str + u_l_n_str
                compound_ids.append(compound_name)

            # Store expanded compounds
            self.expanded_compounds = final_compounds
            self.compound_ids = compound_ids
        else:
            # There is no expansion, every key corresponds to an item
            compound_ids = []
            final_compounds = []
            for name, smiles in library.items():
                compound_ids.append(name)
                final_compounds.append({'smiles': smiles, 'branching': {}})
            self.expanded_compounds = final_compounds
            self.compound_ids = compound_ids

    def generate_script(self, script_path: str) -> str:
        """
        Generate bash script for CompoundLibrary processing.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        script_content = "#!/bin/bash\n"
        script_content += "# CompoundLibrary processing script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += "echo \"Processing compound library\"\n"

        # Create output directories
        script_content += f"""
# Create output directories
mkdir -p "{self.output_folder}"
"""
        if self.covalent:
            script_content += f'mkdir -p "{self.covalent_folder}"\n'

        if self.library_csv:
            # Process existing CSV file
            script_content += f"""
echo "Loading compound library from CSV: {os.path.basename(self.library_csv)}"
cp "{self.library_csv}" "{self.compounds_csv}"
"""
        elif self.library_dict:
            # Generate from dictionary
            # Write JSON files at pipeline time (not SLURM time)
            os.makedirs(self.output_folder, exist_ok=True)
            with open(self.library_dict_json, 'w') as f:
                json.dump(self.library_dict, f, indent=2)

            compounds_data_json = os.path.join(self.output_folder, "compounds_data.json")
            with open(compounds_data_json, 'w') as f:
                json.dump({'expanded_compounds': self.expanded_compounds, 'compound_ids': self.compound_ids}, f, indent=2)

            script_content += f"""
echo "Generating compound library from dictionary ({len(self.expanded_compounds)} compounds)"

# Generate CSV with expanded compounds
python3 -c "
import csv
import json

# Load compound data from JSON file
with open('{compounds_data_json}', 'r') as f:
    data = json.load(f)
expanded_compounds = data['expanded_compounds']
compound_ids = data['compound_ids']

# Write CSV file with standardized format
with open('{self.compounds_csv}', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)

    # Write header with standardized columns
    header = ['id', 'format', 'smiles', 'ccd']

    # Add branching columns
    all_branch_keys = set()
    for comp in expanded_compounds:
        all_branch_keys.update(comp['branching'].keys())
    all_branch_keys = sorted(list(all_branch_keys))
    header.extend(all_branch_keys)

    writer.writerow(header)

    # Write compound data
    for i, compound_data in enumerate(expanded_compounds):
        compound_id = compound_ids[i]
        row = [
            compound_id,
            'smiles',  # format
            compound_data['smiles'],
            ''  # ccd (empty for non-covalent)
        ]

        # Add branching information
        for key in all_branch_keys:
            row.append(compound_data['branching'].get(key, ''))

        writer.writerow(row)

print(f'Generated compound library: {{len(expanded_compounds)}} compounds')
"
"""

        # Generate covalent ligand files if requested
        if self.covalent:
            script_content += f"""
echo "Generating covalent ligand CCD/PKL files"
# Note: This requires Boltz2 cache folder to be available
# The script will generate CCD/PKL files based on the compounds CSV
if [ -n "$BOLTZ_CACHE_FOLDER" ] && [ -d "$BOLTZ_CACHE_FOLDER" ]; then
    python3 "{self.covalent_generation_py}" \\
        "$BOLTZ_CACHE_FOLDER" \\
        "base_config.txt" \\
        "{self.compounds_csv}" \\
        "{self.covalent_folder}" \\
        "{self.covalent_folder}" \\
        "{self.covalent_compounds_csv}" \\
        "{self.conformer_method}"
else
    echo "Warning: BOLTZ_CACHE_FOLDER not set or not found. Skipping covalent file generation."
    echo "To generate covalent files, ensure Boltz2 environment is available and BOLTZ_CACHE_FOLDER is set."
fi
"""

        # Generate summary
        library_type = "Dictionary" if self.library_dict else "CSV file"
        primary_key_str = self.primary_key if self.primary_key else "None"
        covalent_str = str(self.covalent)
        conformer_method_str = self.conformer_method
        compounds_csv_basename = os.path.basename(self.compounds_csv)
        is_covalent = self.covalent

        script_content += f"""
echo "Generating library summary"
python3 -c "
import pandas as pd
import os

# Read library file
df = pd.read_csv('{self.compounds_csv}')
compound_count = len(df)

# Write summary
with open('{self.summary_file}', 'w') as f:
    f.write('Compound Library Summary\\n')
    f.write('========================\\n')
    f.write(f'Library type: {library_type}\\n')
    if '{primary_key_str}' != 'None':
        f.write(f'Primary key: {primary_key_str}\\n')
    f.write(f'Total compounds: {{compound_count}}\\n')
    f.write(f'Covalent ligands: {covalent_str}\\n')
    f.write(f'Conformer method: {conformer_method_str}\\n')
    f.write(f'Output file: {compounds_csv_basename}\\n')
    if {is_covalent}:
        f.write(f'Covalent library folder: covalent_library/\\n')

print(f'Library processed: {{compound_count}} compounds')
print(f'Output: {self.compounds_csv}')
"
"""

        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """
        Get expected output files after CompoundLibrary processing.

        Returns:
            Dictionary with DataStream objects and tables
        """
        # Generate predicted compound IDs if not already done
        if not self.compound_ids and self.library_dict:
            if self.covalent:
                # For covalent: use library keys as-is (they're already expanded)
                self.compound_ids = list(self.library_dict.keys())
            else:
                # For simple library: use library keys as-is
                self.compound_ids = list(self.library_dict.keys())

        # Build tables with rich metadata
        tables = {}
        columns = ["id", "format", "smiles", "ccd"]
        if self.library_dict:
            # Add branching columns
            all_branch_keys = set()
            for comp in self.expanded_compounds:
                all_branch_keys.update(comp['branching'].keys())
            columns.extend(sorted(list(all_branch_keys)))

        tables["compounds"] = TableInfo(
            name="compounds",
            path=self.compounds_csv,
            columns=columns,
            description="Generated compound library with SMILES and metadata",
            count=len(self.expanded_compounds) if self.expanded_compounds else 0
        )

        if self.covalent and self.covalent_compounds_csv:
            tables["covalent_compounds"] = TableInfo(
                name="covalent_compounds",
                path=self.covalent_compounds_csv,
                columns=["id", "format", "smiles", "ccd"],
                description="Covalent compound library with CCD identifiers"
            )

        # Create compounds DataStream
        compounds = DataStream(
            name="compounds",
            ids=self.compound_ids,
            files=[],  # Value-based format - data is in map_table, not individual files
            map_table=self.compounds_csv,
            format="csv"
        )

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": DataStream.empty("sequences", "fasta"),
            "compounds": compounds,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def get_config_display(self) -> List[str]:
        """Get configuration display lines for pipeline output."""
        config_lines = super().get_config_display()

        if isinstance(self.library, dict):
            config_lines.append(f"LIBRARY: Dictionary ({len(self.library)} keys)")
            if self.primary_key:
                config_lines.append(f"PRIMARY KEY: {self.primary_key}")
        else:
            config_lines.append(f"LIBRARY: {os.path.basename(self.library)}")

        config_lines.extend([
            f"COVALENT LIGANDS: {self.covalent}",
            f"VALIDATE SMILES: {self.validate_smiles}"
        ])

        if self.covalent:
            config_lines.append(f"CONFORMER METHOD: {self.conformer_method}")

        return config_lines

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including CompoundLibrary-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "compound_library_params": {
                "library": self.library if isinstance(self.library, str) else "<dictionary>",
                "primary_key": self.primary_key,
                "covalent": self.covalent,
                "validate_smiles": self.validate_smiles,
                "conformer_method": self.conformer_method,
                "num_compounds": len(self.expanded_compounds) if self.expanded_compounds else 0
            }
        })
        return base_dict
