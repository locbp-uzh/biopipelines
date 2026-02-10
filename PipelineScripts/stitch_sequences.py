# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
StitchSequences tool for combining sequences with segment substitutions.

Takes a template sequence and substitutes specific regions with alternative
sequences, generating all combinations (Cartesian product).
"""

import os
import json
from typing import Dict, List, Any, Union
import re
from itertools import product

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

# Import ID mapping utilities from HelpScripts
try:
    from HelpScripts.id_map_utils import get_mapped_ids
except ImportError:
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

    # Lazy path descriptors
    sequences_csv = Path(lambda self: os.path.join(self.output_folder, "sequences.csv"))
    missing_csv = Path(lambda self: os.path.join(self.output_folder, "missing.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "stitch_config.json"))
    helper_script = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_stitch_sequences.py"))

    def __init__(self,
                 template: Union[str, DataStream, StandardizedOutput],
                 substitutions: Dict[str, Union[List[str], DataStream, StandardizedOutput]] = None,
                 indels: Dict[str, Union[List[str], DataStream, StandardizedOutput]] = None,
                 id_map: Dict[str, str] = {"*": "*_<N>"},
                 remove_duplicates: bool = True,
                 **kwargs):
        """
        Initialize StitchSequences configuration.

        Args:
            template: Base sequence - can be:
                - Raw sequence string
                - DataStream/StandardizedOutput with sequences
            substitutions: Position-to-position substitutions from equal-length sequences.
                For each position in the selection, the residue at that position in the
                substitution sequence replaces the residue at that position in the template.
                - Keys: Position strings like "11-19" or "11-19+31-44"
                - Values: DataStream/StandardizedOutput with sequences (must be same length as template)
            indels: Segment replacements where each contiguous segment is replaced.
                For example, "6-7+9-10": "GP" replaces both segments 6-7 and 9-10 with "GP".
                - Keys: Position strings like "50-55" or "6-7+9-10+17-18"
                - Values: List of raw sequences (each segment replaced with full sequence)
            id_map: ID mapping pattern for matching table IDs to sequence IDs
            remove_duplicates: Remove duplicate sequences from output (default True)
            **kwargs: Additional parameters

        Position Syntax:
            - "10-20" → positions 10 to 20 (inclusive, 1-indexed)
            - "10-20+30-40" → positions 10-20 and 30-40
            - "145+147+150" → specific positions 145, 147, and 150

        Processing Order:
            1. Substitutions are applied first (same-length, position-to-position)
            2. Indels are applied second (can change sequence length)
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
        self.remove_duplicates = remove_duplicates

        # Resolve template to DataStream if needed
        self.template_stream = None
        self.template_sequence = None
        if isinstance(template, str):
            self.template_sequence = template
        elif isinstance(template, StandardizedOutput):
            if template.streams.sequences and len(template.streams.sequences) > 0:
                self.template_stream = template.streams.sequences
            else:
                raise ValueError("StandardizedOutput has no sequences for template")
        elif isinstance(template, DataStream):
            self.template_stream = template
        else:
            raise ValueError(f"template must be str, DataStream, or StandardizedOutput, got {type(template)}")

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
                table_path = table_info.info.path if hasattr(table_info, 'info') else str(table_info)
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
                table_path = table_info.info.path if hasattr(table_info, 'info') else str(table_info)
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

        elif isinstance(source, DataStream):
            if source.map_table:
                return {
                    "type": "tool_output",
                    "sequences_file": source.map_table,
                    "source_name": "DataStream",
                    "sequence_ids": source.ids or []
                }
            elif source.files:
                return {
                    "type": "tool_output",
                    "sequences_file": source.files[0],
                    "source_name": "DataStream",
                    "sequence_ids": source.ids or []
                }
            else:
                raise ValueError(f"{name}: DataStream has no files or map_table")

        elif isinstance(source, StandardizedOutput):
            if not hasattr(source, 'tables'):
                raise ValueError(f"{name}: StandardizedOutput must have tables")

            if not hasattr(source.tables, 'sequences'):
                raise ValueError(f"{name}: StandardizedOutput must have tables.sequences")

            sequences_table = source.tables.sequences
            if hasattr(sequences_table, 'info'):
                sequences_file = sequences_table.info.path
            elif isinstance(sequences_table, str):
                sequences_file = sequences_table
            else:
                raise ValueError(f"{name}: tables.sequences must have path or be string")

            result = {
                "type": "tool_output",
                "sequences_file": sequences_file,
                "source_name": source.__class__.__name__,
                "sequence_ids": source.streams.sequences.ids if source.streams.sequences else []
            }

            # Include structure files if available (for PDB residue number mapping)
            if source.streams.structures and len(source.streams.structures) > 0:
                result["structure_files"] = source.streams.structures.files
                result["structure_ids"] = source.streams.structures.ids or []

            return result

        else:
            raise ValueError(f"{name}: unsupported type {type(source)}")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        if isinstance(self.template, str):
            template_display = f"raw sequence ({len(self.template)} chars)"
        elif isinstance(self.template, (DataStream, StandardizedOutput)):
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
                elif isinstance(options, (DataStream, StandardizedOutput)):
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
                elif isinstance(options, (DataStream, StandardizedOutput)):
                    options_display = options.__class__.__name__
                else:
                    options_display = str(type(options))
                config_lines.append(f"  {pos_display}: {options_display}")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate StitchSequences execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        step_tool_name = os.path.basename(self.output_folder)

        # Collect upstream missing tables from all input sources
        # (template, substitutions, indels)
        all_sources = []
        if self.template_stream:
            all_sources.append(self.template_stream)
        for options in self.substitutions.values():
            if options is not None:
                all_sources.append(options)
        for options in self.indels.values():
            if options is not None:
                all_sources.append(options)

        upstream_missing_paths = []
        seen_paths = set()
        for source in all_sources:
            path = self._get_upstream_missing_table_path(source)
            if path and path not in seen_paths:
                upstream_missing_paths.append(path)
                seen_paths.add(path)

        config_data = {
            "template": self.template_info,
            "substitutions": self.substitution_infos,
            "indels": self.indel_infos,
            "id_map": self.id_map,
            "remove_duplicates": self.remove_duplicates,
            "output_csv": self.sequences_csv,
            "missing_csv": self.missing_csv,
            "step_tool_name": step_tool_name,
            "upstream_missing_paths": upstream_missing_paths
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# StitchSequences execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()

        script_content += f"""echo "Running sequence stitching"
echo "Template: {self.template_info['source_name']}"
echo "Substitution regions: {len(self.substitutions)}"
echo "Indel regions: {len(self.indels)}"
"""

        for key_str, info in self.substitution_infos.items():
            script_content += f'echo "  Substitution {key_str}: {info["sequences"]["source_name"]}"\n'

        for key_str, info in self.indel_infos.items():
            script_content += f'echo "  Indel {key_str}: {info["sequences"]["source_name"]}"\n'

        script_content += f"""echo "Output: {self.sequences_csv}"

python "{self.helper_script}" \\
  --config "{self.config_file}"

if [ $? -eq 0 ]; then
    echo "Sequence stitching completed successfully"
    echo "Results written to: {self.sequences_csv}"
else
    echo "Error: Sequence stitching failed"
    exit 1
fi

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after sequence stitching."""
        predicted_ids = self._predict_output_sequence_ids()

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence"],
                description="Stitched sequences with segment substitutions",
                count=len(predicted_ids)
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "cause"],
                description="IDs removed (duplicates or upstream) with removal reason",
                count="variable"
            )
        }

        sequences = DataStream(
            name="sequences",
            ids=predicted_ids,
            files=[self.sequences_csv],
            map_table=self.sequences_csv,
            format="csv"
        )

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": sequences,
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

    def _predict_output_sequence_ids(self) -> List[str]:
        """Predict output sequence IDs based on matched template IDs."""
        # Get template sequence IDs
        if isinstance(self.template, str):
            template_ids = ["seq"]
        elif isinstance(self.template, list):
            template_ids = [f"seq_{i+1}" for i in range(len(self.template))]
        elif self.template_stream:
            template_ids = self.template_stream.ids or ["seq"]
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
        """Predict output IDs when substitutions/indels come from tool outputs."""
        # Get sequence_ids from each substitution source
        sub_ids_list = []
        for options in self.substitutions.values():
            if isinstance(options, StandardizedOutput):
                sub_ids_list.append(options.streams.sequences.ids if options.streams.sequences else [])
            elif isinstance(options, DataStream):
                sub_ids_list.append(options.ids or [])

        # Get sequence_ids from each indel source
        indel_ids_list = []
        for options in self.indels.values():
            if isinstance(options, StandardizedOutput):
                indel_ids_list.append(options.streams.sequences.ids if options.streams.sequences else [])
            elif isinstance(options, DataStream):
                indel_ids_list.append(options.ids or [])

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
        elif isinstance(self.template, (DataStream, StandardizedOutput)):
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
            elif isinstance(options, (DataStream, StandardizedOutput)):
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
            elif isinstance(options, (DataStream, StandardizedOutput)):
                indels_summary[key_str] = options.__class__.__name__
            else:
                indels_summary[key_str] = str(type(options))

        base_dict.update({
            "tool_params": {
                "template": template_summary,
                "substitutions": substitutions_summary,
                "indels": indels_summary,
                "id_map": self.id_map,
                "remove_duplicates": self.remove_duplicates
            }
        })
        return base_dict
