# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
CompoundLibrary configuration for generating and processing ligand libraries.

Handles dictionary-based SMILES library generation, CSV output generation,
and optional covalent ligand CCD/PKL file preparation for Boltz2.
"""

import os
import json
import csv
import itertools
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

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== CompoundLibrary ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== CompoundLibrary ready ==="
"""

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

        Output:
            Streams: compounds (.csv)
            Tables:
                compounds: id | format | smiles | ccd | ...branching_keys
                covalent_compounds: id | format | smiles | ccd (if covalent=True)
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
        self.library_cdxml = None
        self.expanded_compounds = []
        self.compound_ids = []

        # Initialize base class
        super().__init__(**kwargs)

        # For dictionary libraries, do expansion immediately to get compound IDs
        if isinstance(self.library, dict):
            self.library_dict = self.library
            self._expand_library()
        elif isinstance(self.library, str) and self.library.endswith('.cdxml'):
            self.library_cdxml = self.library
            self._expand_cdxml()

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
            if not (self.library.endswith('.csv') or self.library.endswith('.cdxml')):
                raise ValueError("library file must have .csv or .cdxml extension")
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

        if isinstance(self.library, str) and self.library.endswith('.cdxml'):
            # CDXML file - resolve path
            if os.path.exists(self.library):
                self.library_cdxml = self.library
            else:
                project_path = os.path.join(pipeline_folders["biopipelines"], self.library)
                if os.path.exists(project_path):
                    self.library_cdxml = project_path
                else:
                    raise ValueError(f"CDXML file not found: {self.library}")
            # Expand if not already done (file resolved via project path)
            if not self.expanded_compounds:
                self._expand_cdxml()
        elif isinstance(self.library, str):
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
                                # Strip <> from option for display (e.g. "<o-hydroxyphenyl>" -> "o-hydroxyphenyl")
                                clean_option = option.strip("<>") if option.startswith("<") and option.endswith(">") else option
                                new_branching[key] = clean_option
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

    def _expand_cdxml(self):
        """Expand CDXML file with R-group labels into enumerated compounds."""
        try:
            from rdkit import Chem
            from rdkit.Chem import rdmolops
        except ImportError:
            raise ImportError(
                "RDKit is required for CDXML R-group enumeration. "
                "Install with: conda install -c conda-forge rdkit"
            )

        cdxml_path = self.library_cdxml
        if not cdxml_path or not os.path.exists(cdxml_path):
            raise ValueError(f"CDXML file not found: {cdxml_path}")

        # Parse CDXML file
        mols = Chem.MolsFromCDXMLFile(cdxml_path)
        if not mols:
            raise ValueError(f"No molecules found in CDXML file: {cdxml_path}")

        # Filter out None entries (failed parses)
        valid_mols = [m for m in mols if m is not None]
        if not valid_mols:
            raise ValueError(f"All molecules failed to parse from CDXML file: {cdxml_path}")
        if len(valid_mols) < 2:
            raise ValueError(
                "CDXML file must contain at least 2 molecules (1 core + 1 fragment). "
                f"Found {len(valid_mols)} valid molecule(s)."
            )

        # Find R-group dummy atoms on each molecule
        # Dummy atom: atomic num == 0, atom map num > 0
        def get_rgroup_positions(mol):
            positions = []
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 0 and atom.GetAtomMapNum() > 0:
                    positions.append(atom.GetAtomMapNum())
            return positions

        # Identify core: molecule with the most R-group dummy atoms
        mol_rgroups = [(mol, get_rgroup_positions(mol)) for mol in valid_mols]
        if not any(positions for _, positions in mol_rgroups):
            raise ValueError(
                "No R-group labels (R1, R2, ...) found in the CDXML file. "
                "Draw R-group labels on the core scaffold and fragments in ChemDraw."
            )

        # Core = molecule with the most R-group positions
        core_idx = max(range(len(mol_rgroups)), key=lambda i: len(mol_rgroups[i][1]))
        core_mol, core_positions = mol_rgroups[core_idx]

        if not core_positions:
            raise ValueError("No R-group labels found on any molecule in the CDXML file.")

        core_position_set = set(core_positions)

        # Group fragments by R-group position
        # Each non-core molecule must have exactly one R-group dummy atom
        fragments_by_position = {}  # position_num -> list of (mol, smiles)
        for i, (mol, positions) in enumerate(mol_rgroups):
            if i == core_idx:
                continue
            if len(positions) == 0:
                # Molecule without R-group labels — skip with warning
                continue
            if len(positions) > 1:
                raise ValueError(
                    f"R-group fragment (molecule {i+1}) has {len(positions)} R-group labels "
                    f"(R{', R'.join(str(p) for p in positions)}). "
                    "Each fragment must have exactly one R-group label."
                )
            pos = positions[0]
            if pos not in core_position_set:
                raise ValueError(
                    f"Fragment has R-group label R{pos}, but the core scaffold only has "
                    f"positions: {', '.join(f'R{p}' for p in sorted(core_position_set))}."
                )
            if pos not in fragments_by_position:
                fragments_by_position[pos] = []
            fragments_by_position[pos].append(mol)

        # Validate: every core position must have at least one fragment
        missing_positions = core_position_set - set(fragments_by_position.keys())
        if missing_positions:
            raise ValueError(
                f"No fragments found for core position(s): "
                f"{', '.join(f'R{p}' for p in sorted(missing_positions))}. "
                "Draw at least one fragment with each R-group label."
            )

        # Sort positions for deterministic enumeration
        sorted_positions = sorted(fragments_by_position.keys())
        position_labels = [f"R{p}" for p in sorted_positions]
        fragment_groups = [fragments_by_position[p] for p in sorted_positions]

        # Get SMILES for each fragment (for branching metadata)
        def fragment_smiles(mol):
            """Get SMILES for a fragment, removing the dummy atom for display."""
            try:
                # Create an editable copy and remove dummy atoms for display SMILES
                rwmol = Chem.RWMol(mol)
                dummy_indices = [a.GetIdx() for a in rwmol.GetAtoms()
                                 if a.GetAtomicNum() == 0]
                for idx in sorted(dummy_indices, reverse=True):
                    rwmol.RemoveAtom(idx)
                try:
                    Chem.SanitizeMol(rwmol)
                    return Chem.MolToSmiles(rwmol)
                except Exception:
                    return Chem.MolToSmiles(mol)
            except Exception:
                return "?"

        # Enumerate all combinations
        final_compounds = []
        failed_count = 0
        for combo in itertools.product(*fragment_groups):
            assembled = self._assemble_molecule(core_mol, combo)
            if assembled is None:
                failed_count += 1
                continue
            smiles = Chem.MolToSmiles(assembled)
            branching = {}
            for label, frag in zip(position_labels, combo):
                branching[label] = fragment_smiles(frag)
            final_compounds.append({'smiles': smiles, 'branching': branching})

        if not final_compounds:
            raise ValueError(
                f"All {failed_count} R-group combinations failed to assemble. "
                "Check that R-group labels are correctly placed on attachment bonds."
            )

        if failed_count > 0:
            import warnings
            warnings.warn(
                f"{failed_count} R-group combination(s) failed to assemble and were skipped."
            )

        # Generate compound IDs (same logic as dict mode)
        base_name = os.path.splitext(os.path.basename(self.library_cdxml))[0]
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
            compound_name = base_name if num_compounds == 1 else base_name[:characters] + zeros_str + u_l_n_str
            compound_ids.append(compound_name)

        self.expanded_compounds = final_compounds
        self.compound_ids = compound_ids

    @staticmethod
    def _assemble_molecule(core, fragments):
        """
        Assemble a molecule from core scaffold and R-group fragments using molzip.

        Args:
            core: RDKit Mol - core scaffold with R-group dummy atoms
            fragments: tuple of RDKit Mol - one fragment per R-group position

        Returns:
            Assembled RDKit Mol, or None on failure
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import rdmolops

            # Combine core + all fragments into one Mol
            combined = Chem.RWMol(core)
            for frag in fragments:
                combined = Chem.CombineMols(combined, frag)

            # Use molzip to pair dummy atoms by atom map number
            params = Chem.rdmolops.MolzipParams()
            params.label = Chem.rdmolops.MolzipLabel.AtomMapNumber
            assembled = Chem.molzip(combined, params)

            # Sanitize
            Chem.SanitizeMol(assembled)
            return assembled
        except Exception:
            return None

    def _generate_csv_script(self, compounds_data_json: str, source_label: str) -> str:
        """Generate inline Python script to write compounds CSV from JSON data."""
        return f"""
echo "Generating compound library from {source_label} ({len(self.expanded_compounds)} compounds)"

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
        elif self.library_dict or self.library_cdxml:
            # Generate from dictionary or CDXML (both use pre-expanded compounds)
            os.makedirs(self.output_folder, exist_ok=True)

            if self.library_dict:
                with open(self.library_dict_json, 'w') as f:
                    json.dump(self.library_dict, f, indent=2)
                source_label = "dictionary"
            else:
                source_label = f"CDXML ({os.path.basename(self.library_cdxml)})"

            compounds_data_json = os.path.join(self.output_folder, "compounds_data.json")
            with open(compounds_data_json, 'w') as f:
                json.dump({'expanded_compounds': self.expanded_compounds, 'compound_ids': self.compound_ids}, f, indent=2)

            script_content += self._generate_csv_script(compounds_data_json, source_label)

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
        library_type = "Dictionary" if self.library_dict else ("CDXML" if self.library_cdxml else "CSV file")
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
        if self.library_dict or self.library_cdxml:
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
        elif self.library_cdxml:
            config_lines.append(f"LIBRARY: CDXML ({os.path.basename(self.library_cdxml)})")
            # Show R-group positions
            if self.expanded_compounds:
                all_positions = set()
                for comp in self.expanded_compounds:
                    all_positions.update(comp['branching'].keys())
                if all_positions:
                    config_lines.append(f"R-GROUP POSITIONS: {', '.join(sorted(all_positions))}")
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
        if self.library_cdxml:
            library_type = "cdxml"
        elif self.library_dict:
            library_type = "dictionary"
        else:
            library_type = "csv"

        base_dict.update({
            "compound_library_params": {
                "library": self.library if isinstance(self.library, str) else "<dictionary>",
                "library_type": library_type,
                "primary_key": self.primary_key,
                "covalent": self.covalent,
                "validate_smiles": self.validate_smiles,
                "conformer_method": self.conformer_method,
                "num_compounds": len(self.expanded_compounds) if self.expanded_compounds else 0
            }
        })
        return base_dict
