# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ADMET-AI property prediction tool.

Predicts 41 ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity)
properties and 8 physicochemical properties for molecules from SMILES strings.
Optionally compares predictions against DrugBank approved drugs (percentile ranks).
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


class ADMETAI(BaseConfig):
    """
    Pipeline tool for ADMET property prediction of compound libraries.

    Takes compounds (SMILES strings) and predicts 41 ADMET endpoints and
    8 physicochemical properties using the ADMET-AI Chemprop model.
    Optionally generates DrugBank percentile comparison.

    Commonly used for:
    - Drug-likeness evaluation of compound libraries
    - ADMET property profiling early in drug discovery
    - Comparing predictions against approved drugs
    - Filtering compounds by pharmacokinetic properties
    """

    TOOL_NAME = "ADMETAI"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        if env_manager == "pip":
            skip = "" if force_reinstall else """# Check if already installed
if python -c "from admet_ai import ADMETModel" 2>/dev/null; then
    echo "ADMET-AI already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
            return f"""echo "=== Installing ADMET-AI (pip) ==="
{skip}pip install admet-ai
echo "=== ADMET-AI installation complete ==="
"""
        env_file = os.path.join(folders.get("biopipelines", "."), "Environments", "admet_ai.yaml")
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_manager} env list 2>/dev/null | grep -q "AdmetAIEnv"; then
    echo "ADMET-AI already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        return f"""echo "=== Installing ADMET-AI ==="
{skip}{env_manager} env create -f "{env_file}" -y
echo "=== ADMET-AI installation complete ==="
"""

    # Lazy path descriptors
    predictions_csv = Path(lambda self: os.path.join(self.output_folder, "predictions.csv"))
    drugbank_csv = Path(lambda self: os.path.join(self.output_folder, "drugbank.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "admet_ai_config.json"))
    compounds_ds_json = Path(lambda self: os.path.join(self.output_folder, "compounds.json"))
    admet_ai_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_admet_ai.py"))

    def __init__(self,
                 compounds: Union[DataStream, StandardizedOutput],
                 drugbank: bool = True,
                 **kwargs):
        """
        Initialize ADMET-AI property prediction tool.

        Args:
            compounds: Input compounds as DataStream or StandardizedOutput.
                       Must have a map_table with a 'smiles' column.
            drugbank: If True (default), also generate DrugBank percentile comparison.
            **kwargs: Additional parameters

        Output:
            Streams: (none)
            Tables:
                predictions: id | smiles | (41 ADMET + 8 physicochemical property columns)
                drugbank: id | smiles | (percentile columns for each property, only if drugbank=True)
        """
        # Resolve input to DataStream
        if isinstance(compounds, StandardizedOutput):
            self.compounds_stream: DataStream = compounds.streams.compounds
        elif isinstance(compounds, DataStream):
            self.compounds_stream = compounds
        else:
            raise ValueError(f"compounds must be DataStream or StandardizedOutput, got {type(compounds)}")

        self.drugbank = drugbank

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate ADMET-AI parameters."""
        if not self.compounds_stream or len(self.compounds_stream) == 0:
            raise ValueError("compounds cannot be empty")

        if not getattr(self.compounds_stream, 'map_table', None):
            raise ValueError(
                "compounds DataStream must have a map_table with SMILES data. "
                "Use Ligand() or CompoundLibrary() to provide compounds."
            )

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input compounds."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"COMPOUNDS: {len(self.compounds_stream)} molecules",
            f"DRUGBANK COMPARISON: {'yes' if self.drugbank else 'no'}",
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate ADMET-AI execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# ADMET-AI execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_admet_ai()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_admet_ai(self) -> str:
        """Generate the ADMET-AI prediction part of the script."""
        import json

        # Serialize compounds DataStream to JSON for HelpScript to load
        self.compounds_stream.save_json(self.compounds_ds_json)

        config_data = {
            "compounds_json": self.compounds_ds_json,
            "predictions_csv": self.predictions_csv,
            "drugbank_csv": self.drugbank_csv,
            "drugbank": self.drugbank,
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Running ADMET-AI property prediction"
echo "Compounds: {len(self.compounds_stream)}"
echo "DrugBank comparison: {'yes' if self.drugbank else 'no'}"
echo "Output predictions: {self.predictions_csv}"

python "{self.admet_ai_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after ADMET-AI prediction."""
        tables = {
            "predictions": TableInfo(
                name="predictions",
                path=self.predictions_csv,
                columns=["id", "smiles"],  # + dynamic ADMET columns from ADMETModel
                description="ADMET property predictions (41 ADMET + 8 physicochemical properties)",
                count=len(self.compounds_stream) if self.compounds_stream else None
            ),
        }

        if self.drugbank:
            tables["drugbank"] = TableInfo(
                name="drugbank",
                path=self.drugbank_csv,
                columns=["id", "smiles"],  # + dynamic percentile columns
                description="DrugBank approved drug percentile comparison for each property",
                count=len(self.compounds_stream) if self.compounds_stream else None
            )

        return {
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "drugbank": self.drugbank,
            }
        })
        return base_dict
