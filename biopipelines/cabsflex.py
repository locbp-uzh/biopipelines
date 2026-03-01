# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
CABS-Flex tool for fast simulation of protein structure flexibility.

CABS-Flex performs coarse-grained Monte Carlo simulations to model protein
flexibility, producing an ensemble of conformational models and per-residue
RMSF (Root Mean Square Fluctuation) profiles.

Reference: Kurcinski et al. (2019) Bioinformatics, 35(4):694-695
Repository: https://bitbucket.org/lcbio/cabsflex
"""

import os
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo
    from .file_paths import Path
    from .datastream import DataStream, create_map_table
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo
    from file_paths import Path
    from datastream import DataStream, create_map_table


class CABSflex(BaseConfig):
    """
    CABS-Flex: fast simulation of protein structure flexibility.

    Runs coarse-grained Monte Carlo simulations on input protein structures
    to generate conformational ensembles and flexibility profiles (RMSF).
    Optionally rebuilds final models to all-atom representation using MODELLER.

    WARNING: MODELLER requires a license key. Obtain one (free for academics)
    at https://salilab.org/modeller/registration.html and set the environment
    variable KEY_MODELLER before running.
    """

    TOOL_NAME = "CABSflex"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        if env_manager == "pip":
            skip = "" if force_reinstall else """# Check if already installed
if python -c "import CABS" 2>/dev/null; then
    echo "CABS-Flex already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
            return f"""echo "=== Installing CABS-Flex (pip) ==="
{skip}echo "WARNING: MODELLER requires a license key."
echo "Get one at https://salilab.org/modeller/registration.html"
echo "Then set: export KEY_MODELLER=your_key"
pip install -q modeller cabs

echo "=== CABS-Flex installation complete ==="
"""
        skip = "" if force_reinstall else """# Check if already installed
if conda list -n CABSflex cabs 2>/dev/null | grep -q cabs; then
    echo "CABS-Flex already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        return f"""echo "=== Installing CABS-Flex ==="
{skip}echo "WARNING: MODELLER requires a license key."
echo "Get one at https://salilab.org/modeller/registration.html"
echo "Then set: export KEY_MODELLER=your_key"
{env_manager} create -y -n CABSflex python=2.7
{env_manager} install -y -n CABSflex -c salilab modeller
{env_manager} install -y -n CABSflex -c lcbio cabs

echo "=== CABS-Flex installation complete ==="
"""

    # Lazy path descriptors
    structures_ds_json = Path(lambda self: os.path.join(self.output_folder, "structures.json"))
    structures_map = Path(lambda self: os.path.join(self.output_folder, "structures_map.csv"))
    rmsf_all_csv = Path(lambda self: os.path.join(self.output_folder, "rmsf_all.csv"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_cabsflex.py"))

    def __init__(self,
                 structures: Union[DataStream, StandardizedOutput],
                 num_models: int = 10,
                 mc_cycles: int = 50,
                 mc_steps: int = 50,
                 mc_annealing: int = 20,
                 temperature: Optional[str] = None,
                 flexibility: Optional[str] = None,
                 filtering_count: int = 1000,
                 aa_rebuild: bool = True,
                 restraints: str = "ss2",
                 restraints_gap: int = 3,
                 restraints_min: float = 3.8,
                 restraints_max: float = 8.0,
                 weighted_fit: str = "gauss",
                 **kwargs):
        """
        Initialize CABS-Flex configuration.

        Args:
            structures: Input protein structures as DataStream or StandardizedOutput
            num_models: Number of cluster medoids / final models (default: 10)
            mc_cycles: Monte Carlo cycles between trajectory frames (default: 50)
            mc_steps: Monte Carlo steps (default: 50)
            mc_annealing: Temperature annealing cycles (default: 20)
            temperature: Temperature range as "TINIT TFINAL" (default: "1.4 1.4")
            flexibility: Residue flexibility: float, 'bf', 'bfi', 'bfg', or filename (default: None = 1.0)
            filtering_count: Number of low-energy models for clustering (default: 1000)
            aa_rebuild: Rebuild to all-atom with MODELLER (default: True)
            restraints: Restraint mode: 'all', 'ss1', 'ss2' (default: 'ss2')
            restraints_gap: Min gap along chain for restraints (default: 3)
            restraints_min: Min distance in Angstroms for restraints (default: 3.8)
            restraints_max: Max distance in Angstroms for restraints (default: 8.0)
            weighted_fit: Fit method: 'gauss', 'flex', 'ss', 'off', or filename (default: 'gauss')
            **kwargs: Additional parameters

        Output:
            Streams:
                structures: PDB ensemble models (num_models per input structure)
                images: SVG plots (RMSF, RMSD, energy) per input structure
            Tables:
                rmsf_all: id | chain | resi | rmsf  (merged, all structures)
                rmsf.<input_id>: id | chain | resi | rmsf  (per input structure)
                structures: id | file | structures.id
        """
        # Resolve input to DataStream
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(f"structures must be DataStream or StandardizedOutput, got {type(structures)}")

        self.num_models = num_models
        self.mc_cycles = mc_cycles
        self.mc_steps = mc_steps
        self.mc_annealing = mc_annealing
        self.temperature = temperature
        self.flexibility = flexibility
        self.filtering_count = filtering_count
        self.aa_rebuild = aa_rebuild
        self.restraints = restraints
        self.restraints_gap = restraints_gap
        self.restraints_min = restraints_min
        self.restraints_max = restraints_max
        self.weighted_fit = weighted_fit

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate CABS-Flex parameters."""
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if self.num_models < 1:
            raise ValueError("num_models must be >= 1")
        if self.mc_cycles < 1:
            raise ValueError("mc_cycles must be >= 1")
        if self.mc_steps < 1:
            raise ValueError("mc_steps must be >= 1")
        if self.mc_annealing < 1:
            raise ValueError("mc_annealing must be >= 1")
        if self.restraints not in ("all", "ss1", "ss2"):
            raise ValueError(f"restraints must be 'all', 'ss1', or 'ss2', got: {self.restraints}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"INPUT STRUCTURES: {len(self.structures_stream)} files",
            f"NUM MODELS: {self.num_models}",
            f"MC CYCLES: {self.mc_cycles}",
            f"MC ANNEALING: {self.mc_annealing}",
            f"AA REBUILD: {self.aa_rebuild}",
            f"RESTRAINTS: {self.restraints}",
            f"WEIGHTED FIT: {self.weighted_fit}"
        ])
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate CABS-Flex execution script."""
        import json

        # Serialize structures DataStream to JSON
        with open(self.structures_ds_json, 'w') as f:
            json.dump(self.structures_stream.to_dict(), f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# CABS-Flex execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_cabsflex()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_cabsflex(self) -> str:
        """Generate the CABS-Flex execution part of the script."""
        # Build CABSflex flags
        flags = []
        flags.append(f"-y {self.mc_cycles}")
        flags.append(f"-s {self.mc_steps}")
        flags.append(f"-a {self.mc_annealing}")
        flags.append(f"-k {self.num_models}")
        flags.append(f"-n {self.filtering_count}")
        flags.append(f"-g {self.restraints} {self.restraints_gap} {self.restraints_min} {self.restraints_max}")
        flags.append(f"--weighted-fit {self.weighted_fit}")

        if self.temperature:
            flags.append(f"-t {self.temperature}")
        if self.flexibility is not None:
            flags.append(f"-f {self.flexibility}")
        if self.aa_rebuild:
            flags.append("-A")

        # Output only models (M) to save disk space
        flags.append("-o M")

        flags_str = " ".join(flags)

        ds_json = self.structures_ds_json

        return f"""echo "Running CABS-Flex flexibility simulation"
echo "Input structures: {len(self.structures_stream)} files"
echo "Models per structure: {self.num_models}"
echo "Output: {self.output_folder}"

# Read structures JSON and iterate
STRUCTURES_JSON="{ds_json}"
OUTPUT_DIR="{self.output_folder}"

# Process each structure
python -c "
import json, sys
with open('$STRUCTURES_JSON') as f:
    ds = json.load(f)
ids = ds['ids']
files = ds['files']
if len(files) == 1 and len(ids) > 1:
    files = files * len(ids)
for i, sid in enumerate(ids):
    print(sid + '\\t' + files[i])
" | while IFS=$'\\t' read -r STRUCT_ID STRUCT_FILE; do

    echo "=== Processing $STRUCT_ID ==="

    WORK_DIR="$OUTPUT_DIR/${{STRUCT_ID}}"
    mkdir -p "$WORK_DIR"

    CABSflex -i "$STRUCT_FILE" -w "$WORK_DIR" {flags_str}

    if [ $? -ne 0 ]; then
        echo "Error: CABS-Flex failed for $STRUCT_ID"
        continue
    fi

    # Copy model PDBs to output folder with proper naming
    MODEL_IDX=0
    for MODEL_PDB in "$WORK_DIR/output_pdbs"/model_*.pdb; do
        if [ -f "$MODEL_PDB" ]; then
            MODEL_IDX=$((MODEL_IDX + 1))
            cp "$MODEL_PDB" "$OUTPUT_DIR/${{STRUCT_ID}}_${{MODEL_IDX}}.pdb"
        fi
    done
    echo "Copied $MODEL_IDX models for $STRUCT_ID"

    # Copy SVG plots to output folder
    for SVG_FILE in "$WORK_DIR/plots"/*.svg; do
        if [ -f "$SVG_FILE" ]; then
            SVG_NAME=$(basename "$SVG_FILE")
            cp "$SVG_FILE" "$OUTPUT_DIR/${{STRUCT_ID}}_${{SVG_NAME}}"
        fi
    done

    # Copy and reformat RMSF.csv with proper headers
    # Format: tab-separated, no header, columns like "A2\\t3.516"
    # where A2 = chain(A) + residue_index(2)
    if [ -f "$WORK_DIR/plots/RMSF.csv" ]; then
        echo "id,chain,resi,rmsf" > "$OUTPUT_DIR/${{STRUCT_ID}}_RMSF.csv"
        while IFS=$'\\t' read -r RESID RMSF_VAL; do
            # Parse chain letter and residue number from e.g. "A2"
            CHAIN=$(echo "$RESID" | head -c 1)
            RESI=$(echo "$RESID" | tail -c +2)
            echo "$STRUCT_ID,$CHAIN,$RESI,$RMSF_VAL" >> "$OUTPUT_DIR/${{STRUCT_ID}}_RMSF.csv"
        done < "$WORK_DIR/plots/RMSF.csv"
        echo "RMSF saved to $OUTPUT_DIR/${{STRUCT_ID}}_RMSF.csv"
    fi

done

echo "CABS-Flex simulation complete"

# Build merged RMSF and structures map
python "{self.helper_script}" \\
    --structures "{self.structures_ds_json}" \\
    --output_dir "{self.output_folder}" \\
    --rmsf_all_csv "{self.rmsf_all_csv}" \\
    --structures_map "{self.structures_map}" \\
    --num_models {self.num_models}

if [ $? -eq 0 ]; then
    echo "CABS-Flex analysis completed successfully"
else
    echo "Error: CABS-Flex post-processing failed"
    exit 1
fi

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after CABS-Flex execution."""
        input_ids = self.structures_stream.ids

        # --- Output structures: num_models per input ---
        structure_ids = []
        structure_files = []
        provenance = {"structures": []}
        for input_id in input_ids:
            for model_idx in range(1, self.num_models + 1):
                output_id = f"{input_id}_{model_idx}"
                output_file = os.path.join(self.output_folder, f"{output_id}.pdb")
                structure_ids.append(output_id)
                structure_files.append(output_file)
                provenance["structures"].append(input_id)

        create_map_table(self.structures_map, structure_ids, files=structure_files, provenance=provenance)

        structures = DataStream(
            name="structures",
            ids=structure_ids,
            files=structure_files,
            map_table=self.structures_map,
            format="pdb"
        )

        # --- Output images: SVG plots per input ---
        image_ids = []
        image_files = []
        svg_names = ["RMSF_seq.svg", "E_RMSD_A_total.svg", "RMSD_frame_A_replica_0.svg"]
        for input_id in input_ids:
            for svg_name in svg_names:
                img_id = f"{input_id}_{svg_name.replace('.svg', '')}"
                img_file = os.path.join(self.output_folder, f"{input_id}_{svg_name}")
                image_ids.append(img_id)
                image_files.append(img_file)

        images = DataStream(
            name="images",
            ids=image_ids,
            files=image_files,
            format="svg"
        )

        # --- Output tables ---
        rmsf_columns = ["id", "chain", "resi", "rmsf"]

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_map,
                columns=["id", "file", "structures.id"],
                description="CABS-Flex ensemble model structures",
                count=len(structure_ids)
            ),
            "rmsf_all": TableInfo(
                name="rmsf_all",
                path=self.rmsf_all_csv,
                columns=rmsf_columns,
                description="Per-residue RMSF from all input structures (merged)",
                count="variable"
            )
        }

        # Per-ID RMSF tables: res.tables.rmsf_<id>
        for input_id in input_ids:
            rmsf_path = os.path.join(self.output_folder, f"{input_id}_RMSF.csv")
            tables[f"rmsf_{input_id}"] = TableInfo(
                name=f"rmsf_{input_id}",
                path=rmsf_path,
                columns=rmsf_columns,
                description=f"Per-residue RMSF for {input_id}",
                count="variable"
            )

        return {
            "structures": structures,
            "images": images,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "cabsflex_params": {
                "num_models": self.num_models,
                "mc_cycles": self.mc_cycles,
                "mc_steps": self.mc_steps,
                "mc_annealing": self.mc_annealing,
                "temperature": self.temperature,
                "flexibility": self.flexibility,
                "filtering_count": self.filtering_count,
                "aa_rebuild": self.aa_rebuild,
                "restraints": self.restraints,
                "weighted_fit": self.weighted_fit
            }
        })
        return base_dict
