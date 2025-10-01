"""
StitchSequences tool for combining sequences from multiple sequence generation tools.

Takes sequences from multiple tools and combines them by stitching the sequences
together based on position specifications. Takes a base sequence and overlays
other sequences at specified positions to create combined sequences.
"""

import os
import pandas as pd
from typing import Dict, List, Any, Optional, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, DatasheetInfo


class StitchSequences(BaseConfig):
    """
    Pipeline tool for stitching sequences from multiple sequence generation tools.

    Combines sequences from different tools by overlaying/replacing specific positions.
    Takes a base sequence (first tool) and overlays other sequences at specified positions,
    creating combined sequences.
    """

    # Tool identification
    TOOL_NAME = "StitchSequences"
    DEFAULT_ENV = "ProteinEnv"
    DEFAULT_RESOURCES = {"gpu": "T4", "memory": "4GB", "time": "0:30:00"}

    def __init__(self,
                 sequences: List[Union[ToolOutput, StandardizedOutput]],
                 selections: Union[List[Union[str, ToolOutput]], str] = None,
                 **kwargs):
        """
        Initialize StitchSequences configuration.

        Args:
            sequences: List of sequence tools/outputs to combine. First is base, others overlay.
                      e.g., [tool1, tool2] - tool1 provides base, tool2 overlays at positions
            positions: Position specifications for each sequence (except first base sequence).
                      Can be:
                      - List of strings: ["", "10-20+30-40"] - empty string means base, string means overlay positions
                      - List of datasheet references: ["", distances.datasheets.analysis.within]
                      - Single string: "10-20+30-40" - applies to second sequence only
            **kwargs: Additional parameters

        Position Syntax:
            String format:
            - "" → base sequence (no overlay)
            - "10-20" → overlay at residues 10 to 20
            - "10-20+30-40" → overlay at residues 10-20 and 30-40
            - "145+147+150" → overlay at specific residues 145, 147, and 150

            Datasheet format:
            - datasheet_reference → use position values from datasheet column

        Examples:
            # Combine base sequence with overlay at specific positions
            sequences = pipeline.add(StitchSequences(
                sequences=[base_tool, overlay_tool],
                selections=["", "10-20+30-40"]  # overlay_tool at positions 10-20 and 30-40
            ))

            # Use datasheet to specify positions
            sequences = pipeline.add(StitchSequences(
                sequences=[base_tool, overlay_tool],
                selections=["", tool.datasheets.analysis.positions]
            ))

            # Single overlay (second sequence only)
            sequences = pipeline.add(StitchSequences(
                sequences=[base_tool, overlay_tool],
                selections="10-20+30-40"  # applies to overlay_tool
            ))
        """
        self.input_sequences = sequences
        self.position_specs = selections

        # Validate input
        if not sequences or len(sequences) < 2:
            raise ValueError("At least 2 sequence sources are required")

        # Normalize selections parameter
        if selections is None:
            # Default: first is base, others have empty positions (no overlay)
            self.position_specs = [""] + [""] * (len(sequences) - 1)
        elif isinstance(selections, str):
            # Single string: applies to second sequence only
            if len(sequences) != 2:
                raise ValueError("Single selection string only valid with exactly 2 sequence sources")
            self.position_specs = ["", selections]
        elif isinstance(selections, list):
            if len(selections) != len(sequences):
                raise ValueError(f"Selection specifications ({len(selections)}) must match sequence count ({len(sequences)})")
            self.position_specs = selections
        else:
            raise ValueError("Selections must be string, list, or None")

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        for seq_source in sequences:
            if hasattr(seq_source, 'config'):
                self.dependencies.append(seq_source.config)

        for pos_spec in self.position_specs:
            if hasattr(pos_spec, 'config'):
                self.dependencies.append(pos_spec.config)

    def validate_params(self):
        """Validate StitchSequences parameters."""
        if not self.input_sequences:
            raise ValueError("sequences parameter is required")

        if len(self.input_sequences) < 2:
            raise ValueError("At least 2 sequence sources required")

        # Validate all sequence sources are valid types
        for i, seq_source in enumerate(self.input_sequences):
            if not isinstance(seq_source, (ToolOutput, StandardizedOutput)):
                raise ValueError(f"Sequence source {i} must be ToolOutput or StandardizedOutput")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences from previous tools."""
        self.folders = pipeline_folders

        # Get sequence datasheets from each source - no assumptions or fallbacks
        self.sequence_sources = []

        for i, seq_source in enumerate(self.input_sequences):
            # Each source MUST have a sequences datasheet
            if not hasattr(seq_source, 'datasheets'):
                raise ValueError(f"Sequence source {i} must have datasheets")

            if not hasattr(seq_source.datasheets, 'sequences'):
                raise ValueError(f"Sequence source {i} must have datasheets.sequences")

            # Get the sequences datasheet path directly
            sequences_datasheet = seq_source.datasheets.sequences
            if hasattr(sequences_datasheet, 'path'):
                sequences_file = sequences_datasheet.path
            elif isinstance(sequences_datasheet, str):
                sequences_file = sequences_datasheet
            else:
                raise ValueError(f"Sequence source {i} datasheets.sequences must have path or be string path")

            source_info = {
                'tool': seq_source,
                'sequences_file': sequences_file,
                'tool_name': seq_source.__class__.__name__
            }

            self.sequence_sources.append(source_info)

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        config_lines.extend([
            f"SEQUENCE SOURCES: {len(self.input_sequences)}",
            f"BASE SEQUENCE: {self.input_sequences[0].__class__.__name__}",
        ])

        for i in range(1, len(self.input_sequences)):
            positions = self.position_specs[i]
            pos_display = positions if isinstance(positions, str) else f"Datasheet: {positions}"
            config_lines.append(f"OVERLAY {i}: {self.input_sequences[i].__class__.__name__} at {pos_display}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate sequence stitching script.

        Args:
            script_path: Path where script should be written

        Returns:
            Script content as string
        """
        # Output files
        stitched_sequences_csv = os.path.join(self.output_folder, "sequences.csv")

        # Create config file for the stitching process
        config_file = os.path.join(self.output_folder, "stitch_config.json")

        # Build configuration data
        config_data = {
            "sequence_sources": [],
            "position_specs": [],
            "output_csv": stitched_sequences_csv
        }

        # Add sequence source files
        for i, source in enumerate(self.sequence_sources):
            config_data["sequence_sources"].append({
                "index": i,
                "sequences_file": source['sequences_file'],
                "tool_name": source['tool'].__class__.__name__
            })

        # Add position specifications
        for i, pos_spec in enumerate(self.position_specs):
            if isinstance(pos_spec, str):
                config_data["position_specs"].append({
                    "index": i,
                    "type": "fixed",
                    "value": pos_spec
                })
            else:
                # Datasheet reference
                config_data["position_specs"].append({
                    "index": i,
                    "type": "datasheet",
                    "datasheet_path": getattr(pos_spec, 'datasheet_path', ''),
                    "column_name": getattr(pos_spec, 'column_name', 'positions')
                })

        # Write config file
        import json
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Generate script content
        script_content = f"""#!/bin/bash
# StitchSequences execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running sequence stitching"
echo "Sequence sources: {len(self.input_sequences)}"
echo "Base sequence: {self.input_sequences[0].__class__.__name__}"
"""

        # Add overlay echo statements for each overlay tool
        for i in range(1, len(self.input_sequences)):
            script_content += f'echo "Overlay {i}: {self.input_sequences[i].__class__.__name__}"\n'

        script_content += f"""echo "Output: {stitched_sequences_csv}"

# Run Python stitching script
python "{os.path.join(self.folders['HelpScripts'], 'pipe_stitch_sequences.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Sequence stitching completed successfully"
    echo "Results written to: {stitched_sequences_csv}"
else
    echo "Error: Sequence stitching failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """
        Get expected output files after sequence stitching.

        Returns:
            Dictionary with output file paths
        """
        sequences_csv = os.path.join(self.output_folder, "sequences.csv")

        # Get sequence IDs from base tool - no fallbacks
        base_tool = self.input_sequences[0]
        if not hasattr(base_tool, 'sequence_ids'):
            raise ValueError("Base sequence tool must have sequence_ids")
        base_sequence_ids = base_tool.sequence_ids

        datasheets = {
            "sequences": DatasheetInfo(
                name="sequences",
                path=sequences_csv,
                columns=["id", "sequence"],
                description="Stitched sequences combining multiple MPNN tools",
                count=len(base_sequence_ids)
            )
        }

        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [sequences_csv],
            "sequence_ids": base_sequence_ids,
            "datasheets": datasheets,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "sequence_count": len(self.input_sequences),
                "selection_specs": [str(pos) for pos in self.position_specs],
                "base_tool": self.input_sequences[0].__class__.__name__,
                "overlay_tools": [tool.__class__.__name__ for tool in self.input_sequences[1:]]
            }
        })
        return base_dict