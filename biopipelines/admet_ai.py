# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ADMET-AI predictions for compound libraries.

Runs Chemprop-RDKit ADMET endpoint predictions over a compounds stream and
writes one row per input compound with all upstream-reported properties.

Reference:
    Swanson et al. (2024) ADMET-AI: a machine learning ADMET platform for
    evaluation of large-scale chemical libraries. Bioinformatics 40, btae416.
    https://github.com/swansonk14/admet_ai
"""

import os
import json
from typing import Dict, List, Any, Union

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
    ADMET-AI predictions for a compounds stream.

    Reads SMILES from the compounds map_table, batches them through
    ``ADMETModel.predict()`` once, and writes a single CSV with one row per
    input compound. Columns are ``id``, ``smiles``, and every ADMET endpoint
    that the upstream model returns (~40+ regression/classification scores).

    Usage:
        with Pipeline(...):
            Resources(memory="8GB", time="1:00:00")
            library = CompoundLibrary("my_library.csv")
            admet = ADMETAI(compounds=library)
    """

    TOOL_NAME = "ADMETAI"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Install ADMET-AI in a dedicated conda environment."""
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("admet_ai", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check}; then
    echo "ADMET-AI already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("admet_ai", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("admet_ai", env_manager, biopipelines)
        return f"""echo "=== Installing ADMET-AI ==="
{skip}{remove_block}
{env_block}

# Verify installation — instantiating ADMETModel downloads the bundled
# Chemprop-RDKit weights, so a successful import + construction is the
# right verification (and pre-warms the model cache for the first run).
if {env_manager} run -n admet_ai python -c "from admet_ai import ADMETModel; ADMETModel()" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== ADMET-AI installation complete ==="
else
    echo "ERROR: ADMET-AI verification failed (cannot import admet_ai or load model)"
    exit 1
fi
"""

    # Lazy path descriptors — canonical sub-layout.
    admet_csv = Path(lambda self: self.table_path("admet"))
    config_json = Path(lambda self: self.configuration_path("admet_ai_config.json"))
    compounds_json = Path(lambda self: self.configuration_path("compounds_ds.json"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_admet_ai.py"))

    def __init__(self,
                 compounds: Union[DataStream, StandardizedOutput],
                 **kwargs):
        """
        Initialize ADMET-AI prediction tool.

        Args:
            compounds: Input compounds as DataStream or StandardizedOutput.
                SMILES are read from the ``smiles`` column of the compounds
                map_table.

        Output:
            Streams: (none)
            Tables:
                admet: id | smiles | <ADMET endpoint columns from upstream>
        """
        # Resolve compounds input to DataStream
        if isinstance(compounds, StandardizedOutput):
            self.compounds_stream: DataStream = compounds.streams.compounds
        elif isinstance(compounds, DataStream):
            self.compounds_stream = compounds
        else:
            raise ValueError(
                f"compounds must be DataStream or StandardizedOutput, got {type(compounds)}"
            )

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate ADMET-AI parameters."""
        if not self.compounds_stream or len(self.compounds_stream) == 0:
            raise ValueError("compounds parameter is required and must not be empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure folder paths."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.append(f"COMPOUNDS: {len(self.compounds_stream)} molecules")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate ADMET-AI execution script."""
        self.compounds_stream.save_json(self.compounds_json)

        config_data = {
            "compounds_json": self.compounds_json,
            "output_csv": self.admet_csv,
        }
        with open(self.config_json, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# ADMET-AI execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""echo "Running ADMET-AI predictions"
echo "Compounds: {len(self.compounds_stream)}"
echo "Output: {self.admet_csv}"

python "{self.helper_py}" --config "{self.config_json}"

"""
        script_content += self.generate_completion_check_footer()
        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after ADMET-AI execution.

        Column list is the static prefix (``id``, ``smiles``) — the rest are
        ADMET endpoints filled in by the upstream model at runtime. Listing
        them all here would couple the framework to a specific admet-ai
        version's property set.
        """
        tables = {
            "admet": TableInfo(
                name="admet",
                path=self.admet_csv,
                columns=["id", "smiles"],
                description="ADMET-AI endpoint predictions (one row per input compound)"
            )
        }

        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
