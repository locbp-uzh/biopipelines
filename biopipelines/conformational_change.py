# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
ConformationalChange analysis for measuring structural changes between reference and target structures.

Analyzes protein structures to quantify conformational changes via RMSD computed
by PyMOL's align, super, or cealign methods.
"""

import os
from typing import Dict, List, Any, Optional, Tuple, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .biopipelines_io import TableReference
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from biopipelines_io import TableReference


class ConformationalChange(BaseConfig):
    """
    Pipeline tool for analyzing conformational changes between reference and target structures.

    Takes reference structures and target structures as input, aligns them using specified method,
    and calculates multiple metrics to quantify conformational changes in selected regions.

    Generates CSV with comprehensive conformational change analysis.

    Commonly used for:
    - Conformational change analysis between different states
    - Binding-induced structural changes
    - Flexibility region identification
    - Comparison of predicted vs experimental structures
    """

    # Tool identification
    TOOL_NAME = "ConformationalChange"
    TOOL_VERSION = "1.2"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        # Runs in ProteinEnv (PyMOL). Delegate the install so a notebook calling
        # ConformationalChange.install() doesn't have to know it depends on PyMOL.
        from .pymol import PyMOL
        return PyMOL._install_script(folders, env_manager=env_manager,
                                     force_reinstall=force_reinstall, **kwargs)

    # Lazy path descriptors
    analysis_csv = Path(lambda self: self.table_path("conformational_change_analysis"))
    config_file = Path(lambda self: self.configuration_path("conformational_change_config.json"))
    reference_ds_json = Path(lambda self: self.configuration_path("reference_structures.json"))
    target_ds_json = Path(lambda self: self.configuration_path("target_structures.json"))
    analysis_py = Path(lambda self: self.pipe_script_path("pipe_conformational_change.py"))

    def __init__(self,
                 reference_structures: Union[DataStream, StandardizedOutput],
                 target_structures: Union[DataStream, StandardizedOutput],
                 selection: Optional[Union[str, Tuple['TableInfo', str]]] = None,
                 alignment: str = "align",
                 atoms: str = "all",
                 pairing: str = "sequence",
                 cycles: int = 5,
                 cutoff: float = 2.0,
                 frame: Optional[Union[str, Tuple['TableInfo', str]]] = None,
                 **kwargs):
        """
        Initialize conformational change analysis tool.

        Args:
            reference_structures: Reference structures as DataStream or StandardizedOutput
            target_structures: Target structures as DataStream or StandardizedOutput
            selection: Region specification. Accepts:
                      - None: Compare all atoms (whole structure RMSD)
                      - String: '10-20+30-40' (fixed residue ranges for all structures)
                      - Table column reference: (table, "column_name") for per-structure selections
            alignment: Alignment method - "align", "super", or "cealign" (default: "align").
                      Only used when pairing="sequence".
            atoms: Which atoms to use for alignment. Options:
                  - "all" (default): all atoms
                  - "CA": alpha-carbon only
                  - "backbone": backbone atoms (CA+C+N+O)
                  - Any '+'-separated atom names, e.g. "CA+CB"
            pairing: How atoms in the two structures are put into correspondence:
                  - "sequence" (default): PyMOL align/super/cealign. Pairs by SEQUENCE
                    similarity, so residues it cannot match are silently excluded from
                    the RMSD. Correct for comparing homologues.
                  - "ordered": cmd.fit(matchmaker=-1) — pairs the Nth atom of one
                    selection with the Nth of the other. Correct when both structures
                    share a numbering scheme and atom order but NOT their sequence,
                    e.g. a design and the refold of an inverse-folded sequence.
                  - "identifier": cmd.fit(matchmaker=0) — pairs on full atom identifiers
                    (chain/resi/resn/name). Note this includes the RESIDUE NAME, so
                    positions that were mutated are dropped.
            cycles: Outlier-rejection cycles (default 5, PyMOL's default). Each cycle
                  discards the worst-fitting atom pairs and refits, so the reported RMSD
                  describes only the atoms that survived. Set 0 to measure every atom.
            cutoff: Rejection threshold in Angstrom for those cycles (default 2.0).
            frame: Optional selection to superpose on before measuring `selection`.
                  With frame set, the fit is computed on `frame` and the RMSD is then
                  measured over `selection` in that frame WITHOUT refitting. Use it to
                  ask "given the structures are aligned on their fixed core, how far is
                  the designed part from where it was designed?" — a segment allowed its
                  own superposition can score well while sitting in the wrong place.
            **kwargs: Additional parameters

        Selection Syntax (string):
            - '10-20' → residues 10 to 20
            - '10-20+30-40' → residues 10-20 and 30-40
            - '145+147+150' → specific residues 145, 147, and 150

        Alignment Methods:
            - "align": PyMOL align (sequence-dependent, fast)
            - "super": PyMOL super (structure-based superposition)
            - "cealign": PyMOL cealign (combinatorial extension alignment)

        Output:
            Streams: (none)
            Tables:
                changes: id | reference_structure | target_structure | selection |
                         num_aligned_atoms | RMSD | RMSD_before | num_atoms_before |
                         num_residues_aligned | atoms_dropped_pct

        RMSD is measured over the atoms that survived refinement; RMSD_before and
        num_atoms_before are the same quantities before any outlier rejection. When
        those differ materially the reported RMSD is not describing the whole selection
        — atoms_dropped_pct makes that visible instead of silent.
        """
        # Resolve reference structures to DataStream
        if isinstance(reference_structures, StandardizedOutput):
            self.reference_stream: DataStream = reference_structures.streams.structures
        elif isinstance(reference_structures, DataStream):
            self.reference_stream = reference_structures
        else:
            raise ValueError(f"reference_structures must be DataStream or StandardizedOutput, got {type(reference_structures)}")

        # Resolve target structures to DataStream
        if isinstance(target_structures, StandardizedOutput):
            self.target_stream: DataStream = target_structures.streams.structures
        elif isinstance(target_structures, DataStream):
            self.target_stream = target_structures
        else:
            raise ValueError(f"target_structures must be DataStream or StandardizedOutput, got {type(target_structures)}")

        self.selection_spec = selection
        self.alignment_method = alignment
        self.atoms = atoms
        self.pairing = pairing
        self.cycles = cycles
        self.cutoff = cutoff
        self.frame_spec = frame

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate ConformationalChange parameters."""
        if not self.reference_stream or len(self.reference_stream) == 0:
            raise ValueError("reference_structures cannot be empty")

        if not self.target_stream or len(self.target_stream) == 0:
            raise ValueError("target_structures cannot be empty")

        if self.alignment_method not in ["align", "super", "cealign"]:
            raise ValueError(f"Alignment method must be 'align', 'super', or 'cealign', got: {self.alignment_method}")

        if self.pairing not in ("sequence", "ordered", "identifier"):
            raise ValueError(
                f"pairing must be 'sequence', 'ordered', or 'identifier', got: {self.pairing}")
        if self.pairing != "sequence" and self.alignment_method != "align":
            raise ValueError(
                f"alignment={self.alignment_method!r} only applies to pairing='sequence'; "
                f"pairing={self.pairing!r} uses cmd.fit")
        # bool is an int subclass, so cycles=True would otherwise pass as 1 cycle.
        if isinstance(self.cycles, bool) or not isinstance(self.cycles, int) or self.cycles < 0:
            raise ValueError(f"cycles must be a non-negative integer, got: {self.cycles!r}")
        if isinstance(self.cutoff, bool) or not isinstance(self.cutoff, (int, float)):
            raise ValueError(f"cutoff must be a number, got: {type(self.cutoff).__name__}")
        if self.cutoff <= 0:
            raise ValueError(f"cutoff must be positive, got: {self.cutoff!r}")
        if self.pairing == "sequence" and self.alignment_method == "cealign" and self.cycles != 5:
            raise ValueError("cealign has no outlier-rejection cycles; leave cycles at its default")
        if self.frame_spec is not None and self.selection_spec is None:
            raise ValueError("frame requires selection: there is nothing to measure inside the frame")

        if isinstance(self.selection_spec, str):
            _validate_freeform_string("selection", self.selection_spec)
        if isinstance(self.frame_spec, str):
            _validate_freeform_string("frame", self.frame_spec)
        _validate_freeform_string("atoms", self.atoms)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input structures."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        if self.selection_spec is None:
            selection_display = "All atoms (whole structure)"
        elif isinstance(self.selection_spec, TableReference):
            selection_display = f"Column reference: {self.selection_spec.column}"
        else:
            selection_display = self.selection_spec

        config_lines.extend([
            f"REFERENCE STRUCTURES: {len(self.reference_stream)} files",
            f"TARGET STRUCTURES: {len(self.target_stream)} files",
            f"SELECTION: {selection_display}",
            f"PAIRING: {self.pairing}",
            f"ALIGNMENT METHOD: {self.alignment_method if self.pairing == 'sequence' else 'cmd.fit'}",
            f"ATOMS: {self.atoms}",
            f"CYCLES: {self.cycles} (cutoff {self.cutoff} A)"
            + ("  [no outlier rejection]" if self.cycles == 0 else ""),
            f"FRAME: {self.frame_spec if self.frame_spec is not None else 'measured selection itself'}",
            f"METRICS: RMSD ({self.atoms} atoms)"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate conformational change analysis execution script."""
        script_content = "#!/bin/bash\n"
        script_content += "# ConformationalChange execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_analysis()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_analysis(self) -> str:
        """Generate the conformational change analysis part of the script."""
        import json

        # Serialize DataStreams to JSON for pipe_script to load
        self.reference_stream.save_json(self.reference_ds_json)
        self.target_stream.save_json(self.target_ds_json)

        # Handle selection input
        if self.selection_spec is None:
            selection_config = {"type": "all"}
        elif isinstance(self.selection_spec, TableReference):
            selection_config = {"type": "table_column", "table_path": self.selection_spec.path, "column_name": self.selection_spec.column}
        else:
            selection_config = {"type": "fixed", "value": self.selection_spec}

        if self.frame_spec is None:
            frame_config = None
        elif isinstance(self.frame_spec, TableReference):
            frame_config = {"type": "table_column", "table_path": self.frame_spec.path,
                            "column_name": self.frame_spec.column}
        else:
            frame_config = {"type": "fixed", "value": self.frame_spec}

        config_data = {
            "reference_structures_json": self.reference_ds_json,
            "target_structures_json": self.target_ds_json,
            "selection": selection_config,
            "frame": frame_config,
            "alignment_method": self.alignment_method,
            "atoms": self.atoms,
            "pairing": self.pairing,
            "cycles": self.cycles,
            "cutoff": self.cutoff,
            "output_csv": self.analysis_csv
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Running conformational change analysis"
echo "Reference structures: {len(self.reference_stream)}"
echo "Target structures: {len(self.target_stream)}"
echo "Selection: {self.selection_spec}"
echo "Alignment method: {self.alignment_method}"
echo "Output: {self.analysis_csv}"

python "{self.analysis_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after conformational change analysis."""
        tables = {
            "changes": TableInfo(
                name="changes",
                path=self.analysis_csv,
                columns=["id", "reference_structure", "target_structure", "selection",
                        "num_aligned_atoms", "RMSD", "RMSD_before", "num_atoms_before",
                        "num_residues_aligned", "atoms_dropped_pct"],
                description="Conformational change analysis between reference and target structures"
            )
        }

        return {
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        if isinstance(self.selection_spec, TableReference):
            selection_str = f"table_column:{self.selection_spec.column}"
        else:
            selection_str = str(self.selection_spec)
        base_dict.update({
            "tool_params": {
                "selection": selection_str,
                "alignment_method": self.alignment_method,
                "atoms": self.atoms
            }
        })
        return base_dict
