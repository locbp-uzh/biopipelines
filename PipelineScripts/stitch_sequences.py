"""
StitchSequences tool for combining sequences with segment substitutions.

Takes a template sequence and substitutes specific regions with alternative
sequences, generating all combinations (Cartesian product).
"""

import os
import json
from typing import Dict, List, Any, Union
import re

try:
    from .base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, ToolOutput, StandardizedOutput, TableInfo

# Import ID mapping utilities from HelpScripts
try:
    from HelpScripts.id_map_utils import get_mapped_ids
except ImportError:
    # Fallback for when running from different directory
    import sys
    help_scripts_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'HelpScripts')
    if help_scripts_path not in sys.path:
        sys.path.insert(0, help_scripts_path)
    from id_map_utils import get_mapped_ids


class StitchSequences(BaseConfig):
    """
    Pipeline tool for stitching sequences with segment substitutions and indels.

    Takes a template sequence and applies two types of modifications:
    - substitutions: Position-to-position replacement from equal-length sequences
    - indels: Segment replacement where each contiguous segment is replaced

    Usage:
        StitchSequences(
            template=tool1,  # or raw sequence string
            substitutions={
                "11-19+31-44": tool2,  # Copy residues at these positions from tool2 sequences
            },
            indels={
                "50-55": ["GGGG", "SSSS"]  # Replace segment 50-55 with these sequences
            }
        )
    """

    TOOL_NAME = "StitchSequences"
    

    def __init__(self,
                 template: Union[str, ToolOutput, StandardizedOutput],
                 substitutions: Dict[str, Union[List[str], ToolOutput, StandardizedOutput]] = None,
                 indels: Dict[str, Union[List[str], ToolOutput, StandardizedOutput]] = None,
                 id_map: Dict[str, str] = {"*": "*_<N>"},
                 **kwargs):
        """
        Initialize StitchSequences configuration.

        Args:
            template: Base sequence - can be:
                - Raw sequence string
                - ToolOutput/StandardizedOutput with sequences
            substitutions: Position-to-position substitutions from equal-length sequences.
                For each position in the selection, the residue at that position in the
                substitution sequence replaces the residue at that position in the template.
                - Keys: Position strings like "11-19" or "11-19+31-44"
                - Values: ToolOutput with sequences (must be same length as template)
            indels: Segment replacements where each contiguous segment is replaced.
                For example, "6-7+9-10": "GP" replaces both segments 6-7 and 9-10 with "GP".
                - Keys: Position strings like "50-55" or "6-7+9-10+17-18"
                - Values: List of raw sequences (each segment replaced with full sequence)
            id_map: ID mapping pattern for matching table IDs to sequence IDs
            **kwargs: Additional parameters

        Position Syntax:
            - "10-20" → positions 10 to 20 (inclusive, 1-indexed)
            - "10-20+30-40" → positions 10-20 and 30-40
            - "145+147+150" → specific positions 145, 147, and 150

        Processing Order:
            1. Substitutions are applied first (same-length, position-to-position)
            2. Indels are applied second (can change sequence length)

        Examples:
            # Position-to-position substitution from ToolOutput
            stitched = StitchSequences(
                template=pmpnn,  # e.g., 180 residue sequences
                substitutions={
                    "11-19+31-44": lmpnn  # Also 180 residue sequences
                }
            )
            # Residues at positions 11-19 and 31-44 are copied from lmpnn

            # Segment replacement with indels
            stitched = StitchSequences(
                template="MKTAYIAKQRQISFVKSHFS...",
                indels={
                    "11-15": ["AAAAA", "GGGGG"],  # Replace segment with 5-char options
                    "20-22": ["XX", "YYY", "ZZZZ"]  # Can change length
                }
            )

            # Combined: substitutions then indels
            stitched = StitchSequences(
                template=pmpnn,
                substitutions={
                    "6-12+19+21": lmpnn  # Position-to-position from lmpnn
                },
                indels={
                    "50-55": ["LINKER", "GGG"]  # Replace segment 50-55
                }
            )

            # Concatenation mode (no template, integer keys)
            stitched = StitchSequences(
                indels={
                    1: ["AAAA", "BBBB"],      # First segment options
                    2: ["CCCC", "DDDD", "EEEE"]  # Second segment options
                }
            )
            # Output: 2 × 3 = 6 concatenated sequences
        """
        # Handle concatenation mode: no template, integer keys in indels
        if template is None and indels:
            if all(isinstance(k, int) for k in indels.keys()):
                num_segments = len(indels)
                template = "A" * num_segments
                indels = {str(k): v for k, v in indels.items()}

        if template is None:
            raise ValueError("template parameter is required")

        self.template = template
        self.substitutions = substitutions or {}
        self.indels = indels or {}
        self.id_map = id_map

        # Validate substitution keys are position strings or table references
        for pos_key in self.substitutions.keys():
            if not isinstance(pos_key, (str, tuple)):
                raise ValueError(f"Substitution key must be position string or table reference, got {type(pos_key)}")

        # Validate indel keys are position strings or table references
        for pos_key in self.indels.keys():
            if not isinstance(pos_key, (str, tuple)):
                raise ValueError(f"Indel key must be position string or table reference, got {type(pos_key)}")

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

        for options in self.indels.values():
            if isinstance(options, (ToolOutput, StandardizedOutput)):
                if hasattr(options, 'config'):
                    self.dependencies.append(options.config)

    def validate_params(self):
        """Validate StitchSequences parameters."""
        if self.template is None:
            raise ValueError("template is required")

        for pos_range, options in self.substitutions.items():
            if isinstance(pos_range, str):
                self._parse_position_range(pos_range)
            if options is None:
                raise ValueError(f"Substitution options for '{pos_range}' cannot be None")

        for pos_range, options in self.indels.items():
            if isinstance(pos_range, str):
                self._parse_position_range(pos_range)
            if options is None:
                raise ValueError(f"Indel options for '{pos_range}' cannot be None")

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
        """Configure inputs from template, substitutions, and indels."""
        self.folders = pipeline_folders

        self.template_info = self._extract_sequence_info(self.template, "template")

        # Process substitutions - keys can be strings or table references
        self.substitution_infos = {}
        for pos_key, options in self.substitutions.items():
            if isinstance(pos_key, str):
                key_info = {"type": "fixed", "positions": pos_key}
                key_str = pos_key
            elif isinstance(pos_key, tuple) and len(pos_key) == 2:
                table_info, column_name = pos_key
                table_path = table_info.path if hasattr(table_info, 'path') else str(table_info)
                key_info = {"type": "table", "table_path": table_path, "column": column_name}
                key_str = f"table:{column_name}"
            else:
                raise ValueError(f"Invalid substitution position key type: {type(pos_key)}")

            self.substitution_infos[key_str] = {
                "position_key": key_info,
                "sequences": self._extract_sequence_info(options, f"substitution[{key_str}]")
            }

        # Process indels - keys can be strings or table references
        self.indel_infos = {}
        for pos_key, options in self.indels.items():
            if isinstance(pos_key, str):
                key_info = {"type": "fixed", "positions": pos_key}
                key_str = pos_key
            elif isinstance(pos_key, tuple) and len(pos_key) == 2:
                table_info, column_name = pos_key
                table_path = table_info.path if hasattr(table_info, 'path') else str(table_info)
                key_info = {"type": "table", "table_path": table_path, "column": column_name}
                key_str = f"table:{column_name}"
            else:
                raise ValueError(f"Invalid indel position key type: {type(pos_key)}")

            self.indel_infos[key_str] = {
                "position_key": key_info,
                "sequences": self._extract_sequence_info(options, f"indel[{key_str}]")
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

            result = {
                "type": "tool_output",
                "sequences_file": sequences_file,
                "source_name": source.__class__.__name__,
                "sequence_ids": getattr(source, 'sequence_ids', [])
            }

            # Include structure files if available (for PDB residue number mapping)
            if hasattr(source, 'structures') and source.structures:
                result["structure_files"] = source.structures
                result["structure_ids"] = getattr(source, 'structure_ids', [])

            return result

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

        if self.substitutions:
            config_lines.append(f"SUBSTITUTIONS: {len(self.substitutions)} regions (position-to-position)")
            for pos_key, options in self.substitutions.items():
                if isinstance(pos_key, str):
                    pos_display = pos_key
                elif isinstance(pos_key, tuple):
                    pos_display = f"table:{pos_key[1]}"
                else:
                    pos_display = str(pos_key)

                if isinstance(options, list):
                    options_display = f"{len(options)} raw sequences"
                elif isinstance(options, (ToolOutput, StandardizedOutput)):
                    options_display = options.__class__.__name__
                else:
                    options_display = str(type(options))
                config_lines.append(f"  {pos_display}: {options_display}")

        if self.indels:
            config_lines.append(f"INDELS: {len(self.indels)} regions (segment replacement)")
            for pos_key, options in self.indels.items():
                if isinstance(pos_key, str):
                    pos_display = pos_key
                elif isinstance(pos_key, tuple):
                    pos_display = f"table:{pos_key[1]}"
                else:
                    pos_display = str(pos_key)

                if isinstance(options, list):
                    options_display = f"{len(options)} raw sequences"
                elif isinstance(options, (ToolOutput, StandardizedOutput)):
                    options_display = options.__class__.__name__
                else:
                    options_display = str(type(options))
                config_lines.append(f"  {pos_display}: {options_display}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate StitchSequences execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# StitchSequences execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_stitch_sequences()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_stitch_sequences(self) -> str:
        """Generate the sequence stitching execution part of the script."""
        output_csv = os.path.join(self.output_folder, "sequences.csv")
        config_file = os.path.join(self.output_folder, "stitch_config.json")

        config_data = {
            "template": self.template_info,
            "substitutions": self.substitution_infos,
            "indels": self.indel_infos,
            "id_map": self.id_map,
            "output_csv": output_csv
        }

        with open(config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = f"""echo "Running sequence stitching"
echo "Template: {self.template_info['source_name']}"
echo "Substitution regions: {len(self.substitutions)}"
echo "Indel regions: {len(self.indels)}"
"""

        for key_str, info in self.substitution_infos.items():
            script_content += f'echo "  Substitution {key_str}: {info["sequences"]["source_name"]}"\n'

        for key_str, info in self.indel_infos.items():
            script_content += f'echo "  Indel {key_str}: {info["sequences"]["source_name"]}"\n'

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
        """Predict output sequence IDs based on matched template IDs."""
        import re
        from itertools import product

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

        # Check if all substitutions and indels are raw lists with fixed positions
        all_sub_positions_fixed = all(isinstance(k, str) for k in self.substitutions.keys())
        all_subs_raw = all(isinstance(v, list) for v in self.substitutions.values())
        all_indel_positions_fixed = all(isinstance(k, str) for k in self.indels.keys())
        all_indels_raw = all(isinstance(v, list) for v in self.indels.values())

        if all_sub_positions_fixed and all_subs_raw and all_indel_positions_fixed and all_indels_raw:
            # Raw operations with fixed positions - Cartesian product applies
            substitution_counts = [len(v) for v in self.substitutions.values()]
            indel_counts = [len(v) for v in self.indels.values()]
            all_counts = substitution_counts + indel_counts

            predicted_ids = []
            for template_id in template_ids:
                match = re.match(r'^(.+)_\d+$', template_id)
                base_id = match.group(1) if match else template_id

                if all_counts:
                    index_ranges = [range(1, count + 1) for count in all_counts]
                    for combo in product(*index_ranges):
                        suffix = "_".join(str(idx) for idx in combo)
                        predicted_ids.append(f"{base_id}_{suffix}")
                else:
                    predicted_ids.append(base_id)
            return predicted_ids
        else:
            # Matched sequences from tool outputs - match IDs using id_map pattern
            return self._predict_matched_sequence_ids(template_ids)

    def _predict_matched_sequence_ids(self, template_ids: List[str]) -> List[str]:
        """
        Predict output IDs when substitutions/indels come from tool outputs.

        Uses get_mapped_ids for flexible matching that supports:
        - Exact match (source == target)
        - Child match (target = source + suffix)
        - Parent match (source = target + suffix)
        - Sibling match (common ancestor)

        Mirrors the ID matching logic in pipe_stitch_sequences.py.
        """
        from itertools import product

        # Get sequence_ids from each substitution source
        sub_ids_list = []
        for options in self.substitutions.values():
            if isinstance(options, (ToolOutput, StandardizedOutput)):
                sub_ids_list.append(options.sequence_ids)

        # Get sequence_ids from each indel source
        indel_ids_list = []
        for options in self.indels.values():
            if isinstance(options, (ToolOutput, StandardizedOutput)):
                indel_ids_list.append(options.sequence_ids)

        predicted_ids = []
        for template_id in template_ids:
            # Find matching IDs from each substitution source using get_mapped_ids
            sub_matched = []
            for sub_ids in sub_ids_list:
                matches = get_mapped_ids([template_id], sub_ids, self.id_map, unique=False)
                matched = matches.get(template_id, [])
                sub_matched.append(matched)

            # Find matching IDs from each indel source using get_mapped_ids
            indel_matched = []
            for indel_ids in indel_ids_list:
                matches = get_mapped_ids([template_id], indel_ids, self.id_map, unique=False)
                matched = matches.get(template_id, [])
                indel_matched.append(matched)

            # Combine all matched ID lists
            all_matched = sub_matched + indel_matched

            # Check if any source has no matches
            if any(len(m) == 0 for m in all_matched):
                # Skip this template - no valid combinations possible
                continue

            if not all_matched:
                predicted_ids.append(template_id)
                continue

            # Generate cartesian product of indices (1-based)
            index_ranges = [range(len(m)) for m in all_matched]
            for index_combo in product(*index_ranges):
                suffix = "_".join(str(idx + 1) for idx in index_combo)
                predicted_ids.append(f"{template_id}_{suffix}")

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
            if isinstance(pos_key, str):
                key_str = pos_key
            elif isinstance(pos_key, tuple):
                key_str = f"table:{pos_key[1]}"
            else:
                key_str = str(pos_key)

            if isinstance(options, list):
                substitutions_summary[key_str] = f"{len(options)} raw sequences"
            elif isinstance(options, (ToolOutput, StandardizedOutput)):
                substitutions_summary[key_str] = options.__class__.__name__
            else:
                substitutions_summary[key_str] = str(type(options))

        indels_summary = {}
        for pos_key, options in self.indels.items():
            if isinstance(pos_key, str):
                key_str = pos_key
            elif isinstance(pos_key, tuple):
                key_str = f"table:{pos_key[1]}"
            else:
                key_str = str(pos_key)

            if isinstance(options, list):
                indels_summary[key_str] = f"{len(options)} raw sequences"
            elif isinstance(options, (ToolOutput, StandardizedOutput)):
                indels_summary[key_str] = options.__class__.__name__
            else:
                indels_summary[key_str] = str(type(options))

        base_dict.update({
            "tool_params": {
                "template": template_summary,
                "substitutions": substitutions_summary,
                "indels": indels_summary,
                "id_map": self.id_map
            }
        })
        return base_dict
