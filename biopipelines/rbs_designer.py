# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
RBSDesigner tool for designing synthetic ribosome binding sites.

Takes DNA sequences (e.g. from DNAEncoder) and designs RBS sequences using the
Salis thermodynamic model (Salis, Mirsky & Voigt, Nature Biotechnology 2009)
with simulated annealing. Outputs full genes ready for synthesis:
pre_sequence + rbs_sequence + cds_dna.
"""

import os
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


# TIR presets mapping
TIR_PRESETS = {
    "low": 100,
    "medium": 1000,
    "high": 10000,
    "maximum": 100000,
}


class RBSDesigner(BaseConfig):
    """
    Pipeline tool for designing synthetic ribosome binding sites (RBS).

    Uses the Salis thermodynamic model to predict translation initiation rates
    and designs RBS sequences via simulated annealing to match a target TIR.
    Requires ViennaRNA for RNA free energy calculations.

    Citation: Salis, Mirsky & Voigt, Nature Biotechnology 2009.
    """

    TOOL_NAME = "RBSDesigner"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        biopipelines = folders.get("biopipelines", "")
        if env_manager == "pip":
            skip = "" if force_reinstall else """# Check if already installed
if python -c "import RNA; import Bio" 2>/dev/null; then
    echo "RBSDesigner deps already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
            return f"""echo "=== Installing RBSDesigner (pip) ==="
{skip}pip install ViennaRNA biopython

echo "=== RBSDesigner installation complete ==="
"""
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_manager} env list 2>/dev/null | grep -q "rbs_designer"; then
    echo "RBSDesigner environment already installed, skipping. Use force_reinstall=True to reinstall."
    exit 0
fi
"""
        return f"""echo "=== Installing RBSDesigner ==="
{skip}echo "Creating rbs_designer environment with ViennaRNA from bioconda..."
{env_manager} env create -f {biopipelines}/Environments/rbs_designer.yaml -y
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create rbs_designer environment."
    echo "ViennaRNA requires the bioconda channel with flexible channel priority."
    echo "Check {biopipelines}/Environments/rbs_designer.yaml for details."
    exit 1
fi

echo "=== RBSDesigner installation complete ==="
"""

    # Lazy path descriptors
    rbs_csv = Path(lambda self: os.path.join(self.output_folder, "rbs.csv"))
    info_txt = Path(lambda self: os.path.join(self.output_folder, "rbs_info.txt"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "rbs_designer_config.json"))
    designer_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_rbs_designer.py"))

    def __init__(self,
                 sequences: Union[DataStream, StandardizedOutput],
                 tir: Union[str, int, float] = "medium",
                 pre_sequence: str = "",
                 **kwargs):
        """
        Initialize RBS designer tool.

        Args:
            sequences: Input DNA sequences as DataStream or StandardizedOutput
                       (typically from DNAEncoder).
            tir: Target translation initiation rate. Either a numeric value or
                 a preset string: "low" (100), "medium" (1000), "high" (10000),
                 "maximum" (100000).
            pre_sequence: Optional fixed 5'UTR DNA sequence to prepend before
                          the designed RBS in the full gene output.
            **kwargs: Additional parameters.

        Examples:
            # Design RBS for medium expression
            rbs = RBSDesigner(sequences=dna, tir="medium")

            # Design RBS for specific TIR with 5'UTR prefix
            rbs = RBSDesigner(sequences=dna, tir=5000, pre_sequence="AATTAA")
        """
        # Resolve input to DataStream
        if isinstance(sequences, StandardizedOutput):
            self.sequences_stream: DataStream = sequences.streams.sequences
        elif isinstance(sequences, DataStream):
            self.sequences_stream = sequences
        else:
            raise ValueError(f"sequences must be DataStream or StandardizedOutput, got {type(sequences)}")

        # Resolve TIR preset or numeric value
        if isinstance(tir, str):
            tir_lower = tir.lower()
            if tir_lower in TIR_PRESETS:
                self.tir = TIR_PRESETS[tir_lower]
                self.tir_preset = tir_lower
            else:
                raise ValueError(
                    f"Invalid TIR preset '{tir}'. Must be one of: {list(TIR_PRESETS.keys())} "
                    f"or a numeric value."
                )
        else:
            self.tir = float(tir)
            self.tir_preset = None

        self.pre_sequence = pre_sequence.upper()

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate RBSDesigner parameters."""
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")

        if self.tir <= 0:
            raise ValueError(f"TIR must be positive, got {self.tir}")

        # Validate pre_sequence is valid DNA
        if self.pre_sequence:
            valid_bases = set("ATCG")
            invalid = set(self.pre_sequence) - valid_bases
            if invalid:
                raise ValueError(
                    f"pre_sequence contains invalid characters: {invalid}. "
                    f"Must be DNA (A, T, C, G only)."
                )

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        tir_display = f"{self.tir:.0f}"
        if self.tir_preset:
            tir_display += f" ({self.tir_preset})"

        config_lines.extend([
            f"SEQUENCES: {len(self.sequences_stream)} sequences",
            f"TARGET TIR: {tir_display}",
        ])

        if self.pre_sequence:
            config_lines.append(f"PRE_SEQUENCE: {self.pre_sequence}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate script to perform RBS design."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# RBSDesigner execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_run_design()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_run_design(self) -> str:
        """Generate the RBS design part of the script."""
        import json

        config_data = {
            "sequences_csv": self.sequences_stream.map_table,
            "tir": self.tir,
            "pre_sequence": self.pre_sequence,
            "rbs_output": self.rbs_csv,
            "info_output": self.info_txt,
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Designing RBS sequences"
echo "Input sequences: {self.sequences_stream.map_table}"
echo "Target TIR: {self.tir}"
echo "Output folder: {self.output_folder}"

python "{self.designer_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after RBS design."""
        sequence_ids = self.sequences_stream.ids.copy()

        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[],
            map_table=self.rbs_csv,
            format="dna"
        )

        tables = {
            "rbs": TableInfo(
                name="rbs",
                path=self.rbs_csv,
                columns=[
                    "id", "dna_sequence", "rbs_sequence", "full_gene",
                    "dg_total", "tir_predicted", "target_tir", "target_dg",
                    "spacing", "dg_mrna_rrna", "dg_start", "dg_spacing",
                    "dg_mrna", "dg_standby",
                ],
                description="RBS design results with thermodynamic parameters",
                count=len(sequence_ids)
            )
        }

        return {
            "sequences": sequences,
            "tables": tables,
            "output_folder": self.output_folder,
            "info": self.info_txt,
        }
