"""
StitchSequences tool for combining sequences with segment substitutions.

Takes a template sequence and substitutes specific regions with alternative
sequences, generating all combinations (Cartesian product).
"""

import os
import json
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo


class StitchSequences(BaseConfig):
    """
    Pipeline tool for stitching sequences with segment substitutions.

    Takes a template sequence and replaces specific regions with alternative
    sequences, generating all Cartesian product combinations.

    Usage:
        StitchSequences(
            template=tool1,  # or raw sequence string
            substitutions={
                "11-19": tool2,      # ToolOutput with sequences
                "31-44": ["AAAA", "BBBB", "CCCC"]  # or list of raw sequences
            }
        )
    """

    TOOL_NAME = "StitchSequences"
    DEFAULT_ENV = None

    def __init__(self,
                 template: Union[str, ToolOutput, StandardizedOutput],
                 substitutions: Dict[str, Union[List[str], ToolOutput, StandardizedOutput]] = None,
                 id_map: Dict[str, str] = {"*": "*_<N>"},
                 **kwargs):
        """
        Initialize StitchSequences configuration.

        Args:
            template: Base sequence - can be:
                - Raw sequence string
                - ToolOutput/StandardizedOutput with sequences
            substitutions: Dictionary mapping position ranges to replacement options:
                - Keys: Position strings like "11-19" or "31-44+50-55"
                - Values: List of raw sequences OR ToolOutput with sequences
            id_map: ID mapping pattern for matching table IDs to sequence IDs
            **kwargs: Additional parameters

        Position Syntax:
            - "10-20" → positions 10 to 20 (inclusive, 1-indexed)
            - "10-20+30-40" → positions 10-20 and 30-40
            - "145+147+150" → specific positions 145, 147, and 150

        Examples:
            # With ToolOutput
            stitched = StitchSequences(
                template=proteinmpnn_output,
                substitutions={
                    "11-19": ligandmpnn_output,
                    "31-44": ["AAAA", "BBBB"]
                }
            )

            # With raw sequences
            stitched = StitchSequences(
                template="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIAAALEHHHHHH",
                substitutions={
                    "11-19": ["AAAAAAAA", "BBBBBBBB", "CCCCCCCC"],
                    "31-44": ["DDDDDDDDDDDDDD", "EEEEEEEEEEEEEE"]
                }
            )

            # Concatenation mode (no template, integer keys)
            stitched = StitchSequences(
                substitutions={
                    1: ["AAAA", "BBBB"],      # First segment options
                    2: ["CCCC", "DDDD", "EEEE"]  # Second segment options
                }
            )
            # Output: 2 × 3 = 6 concatenated sequences
        """
        # Handle concatenation mode: no template, integer keys
        if template is None and substitutions:
            # Check if all keys are integers
            if all(isinstance(k, int) for k in substitutions.keys()):
                # Auto-generate template and convert keys
                num_segments = len(substitutions)
                template = "A" * num_segments
                substitutions = {str(k): v for k, v in substitutions.items()}

        if template is None:
            raise ValueError("template parameter is required")

        self.template = template
        self.substitutions = substitutions or {}
        self.id_map = id_map

        # Validate substitution keys are position strings or table references
        for pos_key in self.substitutions.keys():
            if not isinstance(pos_key, (str, tuple)):
                raise ValueError(f"Substitution key must be position string or table reference, got {type(pos_key)}")

        # Initialize base class
        super().__init__(**kwargs)

        # Set up dependencies
        if isinstance(self.template, (ToolOutput, StandardizedOutput)):
            if hasattr(self.template, 'config'):
                self.dependencies.append(self.template.config)

        for options in self.substitutions.values():
            if isinstance(options, (ToolOutput, StandardizedOutput)):
                if hasattr(options, 'config'):
                    self.dependencies.append(options.config)

    def validate_params(self):
        """Validate StitchSequences parameters."""
        if self.template is None:
            raise ValueError("template is required")

        for pos_range, options in self.substitutions.items():
            self._parse_position_range(pos_range)
            if options is None:
                raise ValueError(f"Substitution options for '{pos_range}' cannot be None")

    def _parse_position_range(self, pos_range: str) -> List[int]:
        """Parse position range string into list of positions."""
        positions = []
        parts = pos_range.split('+')

        for part in parts:
            part = part.strip()
            if '-' in part:
                start, end = map(int, part.split('-'))
                positions.extend(range(start, end + 1))
            else:
                positions.append(int(part))

        return sorted(set(positions))

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure inputs from template and substitutions."""
        self.folders = pipeline_folders

        self.template_info = self._extract_sequence_info(self.template, "template")

        # Process substitutions - keys can be strings or table references
        self.substitution_infos = {}
        for pos_key, options in self.substitutions.items():
            # Convert key to a serializable format
            if isinstance(pos_key, str):
                key_info = {"type": "fixed", "positions": pos_key}
                key_str = pos_key
            elif isinstance(pos_key, tuple) and len(pos_key) == 2:
                # Table reference (TableInfo, column_name)
                table_info, column_name = pos_key
                table_path = table_info.path if hasattr(table_info, 'path') else str(table_info)
                key_info = {"type": "table", "table_path": table_path, "column": column_name}
                key_str = f"table:{column_name}"
            else:
                raise ValueError(f"Invalid position key type: {type(pos_key)}")

            self.substitution_infos[key_str] = {
                "position_key": key_info,
                "sequences": self._extract_sequence_info(options, f"substitution[{key_str}]")
            }

    def _extract_sequence_info(self, source, name: str) -> Dict[str, Any]:
        """Extract sequence information from various source types."""
        if isinstance(source, str):
            return {
                "type": "raw",
                "sequences": [source],
                "source_name": "raw_sequence"
            }

        elif isinstance(source, list):
            if not all(isinstance(s, str) for s in source):
                raise ValueError(f"{name}: list must contain only strings")
            return {
                "type": "raw_list",
                "sequences": source,
                "source_name": "raw_sequences"
            }

        elif isinstance(source, (ToolOutput, StandardizedOutput)):
            if not hasattr(source, 'tables'):
                raise ValueError(f"{name}: ToolOutput must have tables")

            if not hasattr(source.tables, 'sequences'):
                raise ValueError(f"{name}: ToolOutput must have tables.sequences")

            sequences_table = source.tables.sequences
            if hasattr(sequences_table, 'path'):
                sequences_file = sequences_table.path
            elif isinstance(sequences_table, str):
                sequences_file = sequences_table
            else:
                raise ValueError(f"{name}: tables.sequences must have path or be string")

            return {
                "type": "tool_output",
                "sequences_file": sequences_file,
                "source_name": source.__class__.__name__,
                "sequence_ids": getattr(source, 'sequence_ids', [])
            }

        else:
            raise ValueError(f"{name}: unsupported type {type(source)}")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        if isinstance(self.template, str):
            template_display = f"raw sequence ({len(self.template)} chars)"
        elif isinstance(self.template, (ToolOutput, StandardizedOutput)):
            template_display = self.template.__class__.__name__
        else:
            template_display = str(type(self.template))

        config_lines.append(f"TEMPLATE: {template_display}")
        config_lines.append(f"SUBSTITUTIONS: {len(self.substitutions)} regions")

        for pos_key, options in self.substitutions.items():
            # Format position key
            if isinstance(pos_key, str):
                pos_display = pos_key
            elif isinstance(pos_key, tuple):
                pos_display = f"table:{pos_key[1]}"
            else:
                pos_display = str(pos_key)

            # Format options
            if isinstance(options, list):
                options_display = f"{len(options)} raw sequences"
            elif isinstance(options, (ToolOutput, StandardizedOutput)):
                options_display = options.__class__.__name__
            else:
                options_display = str(type(options))
            config_lines.append(f"  {pos_display}: {options_display}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate sequence stitching script."""
        output_csv = os.path.join(self.output_folder, "sequences.csv")
        config_file = os.path.join(self.output_folder, "stitch_config.json")

        config_data = {
            "template": self.template_info,
            "substitutions": self.substitution_infos,
            "id_map": self.id_map,
            "output_csv": output_csv
        }

        os.makedirs(os.path.dirname(config_file), exist_ok=True)
        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = f"""#!/bin/bash
# StitchSequences execution script
# Generated by BioPipelines pipeline system

{self.generate_completion_check_header()}

echo "Running sequence stitching"
echo "Template: {self.template_info['source_name']}"
echo "Substitution regions: {len(self.substitutions)}"
"""

        for key_str, info in self.substitution_infos.items():
            script_content += f'echo "  {key_str}: {info["sequences"]["source_name"]}"\n'

        script_content += f"""echo "Output: {output_csv}"

python "{os.path.join(self.folders['HelpScripts'], 'pipe_stitch_sequences.py')}" \\
  --config "{config_file}"

if [ $? -eq 0 ]; then
    echo "Sequence stitching completed successfully"
    echo "Results written to: {output_csv}"
else
    echo "Error: Sequence stitching failed"
    exit 1
fi

{self.generate_completion_check_footer()}
"""

        return script_content

    def get_output_files(self) -> Dict[str, List[str]]:
        """Get expected output files after sequence stitching."""
        sequences_csv = os.path.join(self.output_folder, "sequences.csv")
        predicted_ids = self._predict_output_sequence_ids()

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=sequences_csv,
                columns=["id", "sequence"],
                description="Stitched sequences with segment substitutions",
                count=len(predicted_ids)
            )
        }

        return {
            "structures": [],
            "structure_ids": [],
            "compounds": [],
            "compound_ids": [],
            "sequences": [sequences_csv],
            "sequence_ids": predicted_ids,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def _predict_output_sequence_ids(self) -> List[str]:
        """Predict output sequence IDs based on Cartesian product of all options."""
        import re

        # Get template sequence IDs
        if isinstance(self.template, str):
            template_ids = ["seq"]
        elif isinstance(self.template, list):
            template_ids = [f"seq_{i+1}" for i in range(len(self.template))]
        elif isinstance(self.template, (ToolOutput, StandardizedOutput)):
            if hasattr(self.template, 'sequence_ids'):
                template_ids = self.template.sequence_ids
            else:
                template_ids = ["seq"]
        else:
            template_ids = ["seq"]

        # Get counts for each substitution region
        substitution_counts = []
        for options in self.substitutions.values():
            if isinstance(options, list):
                substitution_counts.append(len(options))
            elif isinstance(options, (ToolOutput, StandardizedOutput)):
                if hasattr(options, 'sequence_ids'):
                    seq_ids = options.sequence_ids
                    base_ids = set()
                    for seq_id in seq_ids:
                        match = re.match(r'^(.+)_\d+$', seq_id)
                        base_ids.add(match.group(1) if match else seq_id)
                    if base_ids:
                        avg_per_base = len(seq_ids) // len(base_ids)
                        substitution_counts.append(max(1, avg_per_base))
                    else:
                        substitution_counts.append(1)
                else:
                    substitution_counts.append(1)
            else:
                substitution_counts.append(1)

        # Calculate total combinations per template ID
        total_per_template = 1
        for count in substitution_counts:
            total_per_template *= count

        # Generate predicted IDs
        predicted_ids = []
        for template_id in template_ids:
            match = re.match(r'^(.+)_\d+$', template_id)
            base_id = match.group(1) if match else template_id

            for n in range(1, total_per_template + 1):
                predicted_ids.append(f"{base_id}_{n}")

        return predicted_ids

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()

        if isinstance(self.template, str):
            template_summary = f"raw_sequence({len(self.template)} chars)"
        elif isinstance(self.template, (ToolOutput, StandardizedOutput)):
            template_summary = self.template.__class__.__name__
        else:
            template_summary = str(type(self.template))

        substitutions_summary = {}
        for pos_key, options in self.substitutions.items():
            # Format key
            if isinstance(pos_key, str):
                key_str = pos_key
            elif isinstance(pos_key, tuple):
                key_str = f"table:{pos_key[1]}"
            else:
                key_str = str(pos_key)

            # Format value
            if isinstance(options, list):
                substitutions_summary[key_str] = f"{len(options)} raw sequences"
            elif isinstance(options, (ToolOutput, StandardizedOutput)):
                substitutions_summary[key_str] = options.__class__.__name__
            else:
                substitutions_summary[key_str] = str(type(options))

        base_dict.update({
            "tool_params": {
                "template": template_summary,
                "substitutions": substitutions_summary,
                "id_map": self.id_map
            }
        })
        return base_dict
