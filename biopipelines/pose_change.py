# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
PoseChange - Measure ligand pose distance between reference and sample structures.

Calculates RMSD and distance metrics for ligand poses between a reference
holo structure (e.g., from X-ray crystallography) and designed holo structures.
Useful for validating whether design tools reproduce known binding poses.
"""

import os
from typing import Dict, List, Any, Union, Optional

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import Resolve
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import Resolve


class PoseChange(BaseConfig):
    """
    Measures ligand pose distance between reference and sample structures.

    Compares ligand conformations between a reference holo structure (e.g., XRC)
    and designed/predicted holo structures. Calculates RMSD, centroid distance,
    and orientation metrics for ligand poses.

    Commonly used for:
    - Validating binding pose prediction accuracy
    - Comparing designed structures to experimental references
    - Analyzing ligand pose consistency across designs
    """

    TOOL_NAME = "PoseChange"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        # Runs in ProteinEnv (PyMOL). Delegate the install so a notebook calling
        # PoseChange.install() doesn't have to know it depends on PyMOL.
        from .pymol import PyMOL
        return PyMOL._install_script(folders, env_manager=env_manager,
                                     force_reinstall=force_reinstall, **kwargs)

    # Lazy path descriptors
    analysis_csv = Path(lambda self: self.table_path("pose_analysis"))
    config_file = Path(lambda self: self.configuration_path("pose_config.json"))
    samples_ds_json = Path(lambda self: self.configuration_path("samples_structures.json"))
    reference_ds_json = Path(lambda self: self.configuration_path("reference_structure.json"))
    reference_ligand_json = Path(lambda self: self.configuration_path(".reference_ligand_compounds.json"))
    sample_ligand_json = Path(lambda self: self.configuration_path(".sample_ligand_compounds.json"))
    pose_change_py = Path(lambda self: self.pipe_script_path("pipe_pose_change.py"))

    def __init__(self,
                 reference_structure: Union[DataStream, StandardizedOutput],
                 sample_structures: Union[DataStream, StandardizedOutput],
                 reference_ligand: Union[DataStream, StandardizedOutput],
                 sample_ligand: Union[DataStream, StandardizedOutput, None] = None,
                 reference_alignment: Optional[str] = None,
                 target_alignment: Optional[str] = None,
                 **kwargs):
        """
        Initialize PoseChange tool.

        Args:
            reference_structure: Reference structure(s). One reference compares every
                sample against it; multiple references are matched to each sample by
                id/provenance at runtime (e.g. each refold vs the input pose it derives
                from).
            sample_structures: Sample structures to compare against the reference(s).
            reference_ligand: The ligand in the reference structure. Pass a
                ``Ligand`` / any compounds-producing tool output — the residue
                ``code`` is read from the compounds stream's map_table at runtime
                (Ligand Contract).
            sample_ligand: The ligand in the sample structures (default: same as
                reference_ligand). Pass a ``Ligand`` / compounds stream.
            reference_alignment: PyMOL selection for reference structure alignment
                                (default: "not resn {reference_ligand}" — everything except the ligand,
                                built at runtime once the ligand code is resolved)
            target_alignment: PyMOL selection for target structure alignment
                             (default: "not resn {sample_ligand}" — everything except the ligand,
                             built at runtime once the ligand code is resolved)
            **kwargs: Additional parameters passed to BaseConfig

        Output:
            Streams: (none)
            Tables:
                changes: id | target_structure | reference_structure | ligand_rmsd | centroid_distance | orientation_angle | orientation_axis | alignment_rmsd | num_ligand_atoms
        """
        # Resolve reference structure to DataStream
        if isinstance(reference_structure, StandardizedOutput):
            self.reference_stream: DataStream = reference_structure.streams.structures
        elif isinstance(reference_structure, DataStream):
            self.reference_stream = reference_structure
        else:
            raise ValueError(f"reference_structure must be DataStream or StandardizedOutput, got {type(reference_structure)}")

        # Resolve sample structures to DataStream
        if isinstance(sample_structures, StandardizedOutput):
            self.samples_stream: DataStream = sample_structures.streams.structures
        elif isinstance(sample_structures, DataStream):
            self.samples_stream = sample_structures
        else:
            raise ValueError(f"sample_structures must be DataStream or StandardizedOutput, got {type(sample_structures)}")

        # Ligand residue codes are always read from the compounds stream at
        # runtime. sample_ligand defaults to the reference ligand stream.
        self.reference_ligand_stream = self._resolve_ligand(reference_ligand, "reference_ligand")
        self.sample_ligand_stream = (self.reference_ligand_stream if sample_ligand is None
                                     else self._resolve_ligand(sample_ligand, "sample_ligand"))

        # Alignment selections embed the ligand code, unknown until runtime, so
        # an unspecified alignment is built as "not resn <code>" by the script.
        self.reference_alignment = reference_alignment
        self.target_alignment = target_alignment

        super().__init__(**kwargs)

    @staticmethod
    def _resolve_ligand(ligand, param_name):
        """Resolve a Ligand/compounds output to its compounds DataStream."""
        if isinstance(ligand, StandardizedOutput):
            return ligand.streams.compounds
        if isinstance(ligand, DataStream):
            return ligand
        raise ValueError(f"{param_name} must be a Ligand/compounds DataStream or "
                         f"StandardizedOutput, got {type(ligand)}")

    def validate_params(self):
        """Validate tool parameters."""
        if not self.reference_stream or len(self.reference_stream) == 0:
            raise ValueError("reference_structure cannot be empty")

        if not self.samples_stream or len(self.samples_stream) == 0:
            raise ValueError("sample_structures cannot be empty")

        if len(self.reference_ligand_stream) == 0:
            raise ValueError("reference_ligand compounds stream is empty")
        if len(self.sample_ligand_stream) == 0:
            raise ValueError("sample_ligand compounds stream is empty")

        # Explicit alignments, when given, must be non-empty. Unspecified ones
        # (None) are built at runtime from the resolved ligand code.
        if self.reference_alignment is not None and not self.reference_alignment:
            raise ValueError("reference_alignment cannot be empty")
        if self.target_alignment is not None and not self.target_alignment:
            raise ValueError("target_alignment cannot be empty")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        runtime = "(code from compounds stream at runtime)"
        ref_aln_display = self.reference_alignment if self.reference_alignment else "(not resn <ligand> at runtime)"
        tgt_aln_display = self.target_alignment if self.target_alignment else "(not resn <ligand> at runtime)"

        config_lines.extend([
            f"REFERENCE: {len(self.reference_stream)} structure",
            f"SAMPLE STRUCTURES: {len(self.samples_stream)}",
            f"REFERENCE LIGAND: {runtime}",
            f"SAMPLE LIGAND: {runtime}",
            f"REFERENCE ALIGNMENT: {ref_aln_display}",
            f"TARGET ALIGNMENT: {tgt_aln_display}",
            f"METRICS: RMSD, centroid, orientation"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate PoseChange execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# PoseChange execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_pose_change()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_pose_change(self) -> str:
        """Generate the pose distance analysis part of the script."""
        import json

        # Serialize DataStreams to JSON for runtime resolution
        self.samples_stream.save_json(self.samples_ds_json)
        self.reference_stream.save_json(self.reference_ds_json)

        reference_id = self.reference_stream.ids[0]
        # >1 reference: each sample is matched to its own reference by id/provenance
        # in the pipe script. ==1: the single reference path is resolved in bash below.
        multi_reference = len(self.reference_stream) > 1

        # Ligand codes and (default) alignments are resolved at runtime from the
        # compounds streams, so they go in as empty placeholders the script fills.
        config_data = {
            "reference_json": self.reference_ds_json,
            "reference_id": reference_id,
            "multi_reference": multi_reference,
            "reference_ligand": "",
            "samples_json": self.samples_ds_json,
            "ligand": "",
            "reference_alignment": self.reference_alignment or "",
            "target_alignment": self.target_alignment or "",
            "output_csv": self.analysis_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Resolve both ligand residue codes from their compounds streams.
        # Colons are stripped for the PyMOL residue match.
        self.reference_ligand_stream.save_json(self.reference_ligand_json)
        self.sample_ligand_stream.save_json(self.sample_ligand_json)
        ref_lig_id = self.reference_ligand_stream.ids[0]
        smp_lig_id = self.sample_ligand_stream.ids[0]
        resolve_blocks = (
            f'REF_LIG_RAW={Resolve.stream_item(self.reference_ligand_json, ref_lig_id, column="code")}\n'
            'REF_LIGAND="${REF_LIG_RAW//:/}"\n'
            f'SMP_LIG_RAW={Resolve.stream_item(self.sample_ligand_json, smp_lig_id, column="code")}\n'
            'SMP_LIGAND="${SMP_LIG_RAW//:/}"\n'
        )

        # Patch the resolved codes in, then build "not resn <code>" alignment
        # defaults for any alignment the caller left unset.
        patch_block = f"""# Patch config with resolved reference PDB path and ligand codes.
jq --arg path "$REFERENCE_PDB" '.reference_pdb = $path' "{self.config_file}" > "{self.config_file}.tmp" && mv "{self.config_file}.tmp" "{self.config_file}"
jq --arg v "$REF_LIGAND" '.reference_ligand = $v' "{self.config_file}" > "{self.config_file}.tmp" && mv "{self.config_file}.tmp" "{self.config_file}"
jq --arg v "$SMP_LIGAND" '.ligand = $v' "{self.config_file}" > "{self.config_file}.tmp" && mv "{self.config_file}.tmp" "{self.config_file}"

# Build "not resn <code>" alignment defaults where none was supplied.
if [ -z "$(jq -r '.reference_alignment' "{self.config_file}")" ]; then
    jq --arg v "not resn $REF_LIGAND" '.reference_alignment = $v' "{self.config_file}" > "{self.config_file}.tmp" && mv "{self.config_file}.tmp" "{self.config_file}"
fi
if [ -z "$(jq -r '.target_alignment' "{self.config_file}")" ]; then
    jq --arg v "not resn $SMP_LIGAND" '.target_alignment = $v' "{self.config_file}" > "{self.config_file}.tmp" && mv "{self.config_file}.tmp" "{self.config_file}"
fi"""

        return f"""REFERENCE_ID={Resolve.stream_ids(self.reference_ds_json, index=0)}
REFERENCE_PDB={Resolve.stream_item(self.reference_ds_json, '$REFERENCE_ID')}
{resolve_blocks}
echo "Running pose distance analysis"
echo "Reference: $REFERENCE_PDB"
echo "Reference ligand: $REF_LIGAND"
echo "Sample structures: {len(self.samples_stream)}"
echo "Sample ligand: $SMP_LIGAND"
echo "Output: {self.analysis_csv}"

{patch_block}

python "{self.pose_change_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files."""
        columns = [
            "id",
            "target_structure",
            "reference_structure",
            "ligand_rmsd",
            "centroid_distance",
            "orientation_angle",
            "orientation_axis",
            "alignment_rmsd",
            "num_ligand_atoms",
            "reference_alignment",
            "target_alignment"
        ]

        tables = {
            "changes": TableInfo(
                name="changes",
                path=self.analysis_csv,
                columns=columns,
                description="Ligand pose distance analysis comparing sample poses to reference"
            )
        }

        return {
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "reference_alignment": self.reference_alignment,
                "target_alignment": self.target_alignment
            }
        })
        return base_dict
