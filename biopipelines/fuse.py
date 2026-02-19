# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Fuse configuration for fusion sequence generation.

Creates fusion sequences by linking multiple sequences with flexible linkers
of varying lengths, generating all possible combinations for downstream
folding and analysis. Works with both protein and DNA sequences.

Each input slot can be a DataStream or StandardizedOutput with multiple IDs.
The output is the cartesian product of all slot options × linker length options.
"""

import os
import json
from typing import Dict, List, Any, Union
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


class Fuse(BaseConfig):
    """
    Configuration for fusion sequence generation.

    Generates all combinations of sequence fusions with variable linker lengths,
    producing meaningful sequence IDs and structured output for downstream
    folding tools like AlphaFold. Works with both protein and DNA sequences.

    Each slot in the sequences list can contain multiple IDs (DataStream with N sequences).
    The output is the cartesian product of all slots × linker lengths.

    Example:
        # 3 sequences from pmpnn × 2 from lmpnn × linker lengths 2-3
        # = 3 × 2 × 2 = 12 fusion combinations
        fuse = Fuse(
            sequences=[pmpnn, lmpnn],
            linker_lengths=["2-3"]
        )
    """

    TOOL_NAME = "Fuse"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Fuse ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== Fuse ready ==="
"""

    # Lazy path descriptors
    queries_csv = Path(lambda self: os.path.join(self.output_folder, f"{self._get_job_base()}_queries.csv"))
    queries_fasta = Path(lambda self: os.path.join(self.output_folder, f"{self._get_job_base()}_queries.fasta"))
    fuse_config_json = Path(lambda self: os.path.join(self.output_folder, "fuse_config.json"))
    fuse_queries_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_fuse_queries.py"))

    def __init__(self,
                 sequences: List[Union[DataStream, StandardizedOutput]],
                 name: str = "",
                 linker: str = "GGGGSGGGGSGGGGSGGGGS",
                 linker_lengths: List[str] = None,
                 **kwargs):
        """
        Initialize Fuse configuration.

        Args:
            sequences: List of DataStreams or StandardizedOutputs, one per fusion slot.
                      Each slot can contain multiple IDs — the output is the
                      cartesian product of all slot options × linker lengths.
            name: Job name for output files
            linker: Linker sequence to use (will be truncated to specified lengths)
            linker_lengths: List of length ranges for each junction (e.g., ["1-6", "1-6"])
            **kwargs: Additional parameters
        """
        # Resolve input slots
        self.input_slots = self._resolve_slots(sequences)

        # Validate we have enough slots for fusion
        if len(self.input_slots) < 2:
            raise ValueError("Fuse requires at least 2 sequence slots for fusion")

        # Store Fuse-specific parameters
        self.name = name
        self.linker = linker

        # Set default linker_lengths based on number of slots if not provided
        if linker_lengths is None:
            expected_junctions = len(self.input_slots) - 1
            self.linker_lengths = ["1-6"] * expected_junctions
        else:
            self.linker_lengths = linker_lengths

        # Validate linker_lengths matches slot count
        expected_junctions = len(self.input_slots) - 1
        if len(self.linker_lengths) != expected_junctions:
            raise ValueError(f"linker_lengths must have {expected_junctions} entries for {len(self.input_slots)} slots")

        # Initialize base class
        super().__init__(**kwargs)

    def _resolve_slots(self, sequences: List[Union[DataStream, StandardizedOutput]]) -> List[Dict[str, Any]]:
        """
        Resolve sequence inputs to per-slot info dicts.

        Args:
            sequences: List of DataStream or StandardizedOutput objects

        Returns:
            List of slot dicts: [{"ids": [...], "map_table": str, "files": [...]}]
        """
        if not isinstance(sequences, list):
            raise ValueError(f"sequences must be a list, got {type(sequences)}")

        slots = []
        for i, item in enumerate(sequences):
            slot = self._resolve_single_slot(item, i)
            slots.append(slot)
        return slots

    def _resolve_single_slot(self, item: Union[DataStream, StandardizedOutput], index: int) -> Dict[str, Any]:
        """Resolve a single slot to its info dict."""
        if isinstance(item, StandardizedOutput):
            if item.streams.sequences and len(item.streams.sequences) > 0:
                ds = item.streams.sequences
            elif item.streams.structures and len(item.streams.structures) > 0:
                ds = item.streams.structures
            else:
                raise ValueError(f"Slot {index}: StandardizedOutput has no sequences or structures")
            return self._datastream_to_slot(ds, index)

        elif isinstance(item, DataStream):
            return self._datastream_to_slot(item, index)

        else:
            raise ValueError(f"Slot {index}: expected DataStream or StandardizedOutput, got {type(item)}")

    @staticmethod
    def _datastream_to_slot(ds: DataStream, index: int) -> Dict[str, Any]:
        """Convert a DataStream to a slot dict."""
        if len(ds) == 0:
            raise ValueError(f"Slot {index}: DataStream is empty")

        slot = {
            "ids": list(ds.ids),
            "map_table": ds.map_table or "",
            "files": list(ds.files) if ds.files else []
        }

        if not slot["map_table"] and not slot["files"]:
            raise ValueError(f"Slot {index}: DataStream has no map_table and no files")

        return slot

    def _get_job_base(self) -> str:
        """Get job base name for file naming."""
        return self.name if self.name else "fuse"

    @staticmethod
    def _parse_length_spec(spec: str) -> List[int]:
        """Parse length specification like '1-6' or '3+5-7'."""
        lengths = []
        if '+' in spec:
            for part in spec.split('+'):
                lengths.extend(Fuse._parse_length_spec(part.strip()))
        elif '-' in spec and not spec.startswith('-'):
            if spec.count('-') == 1:
                start, end = map(int, spec.split('-'))
                lengths.extend(range(start, end + 1))
            else:
                lengths.append(int(spec))
        else:
            lengths.append(int(spec))
        return lengths

    def _predict_sequence_ids(self) -> tuple:
        """
        Predict sequence IDs and provenance from the cartesian product of
        all slot IDs × linker length combinations.

        ID format: {S1_id}_{J0_len}_{S2_id}_{J1_len}_{S3_id}
        where S = slot sequence, J = junction linker length.

        Returns:
            (sequence_ids, provenance) where provenance is:
            {"sequences_1": [...], "sequences_2": [...], ...}
        """
        # Get per-slot IDs
        slot_id_lists = [slot["ids"] for slot in self.input_slots]

        # Get per-junction linker length options as strings
        junction_length_lists = []
        for spec in self.linker_lengths:
            lengths = self._parse_length_spec(spec)
            junction_length_lists.append([str(l) for l in lengths])

        # Build cartesian product axes: [slot0_ids, junction0_lengths, slot1_ids, junction1_lengths, slot2_ids]
        axes = []
        for i, slot_ids in enumerate(slot_id_lists):
            axes.append(slot_ids)
            if i < len(junction_length_lists):
                axes.append(junction_length_lists[i])

        # Generate all combinations
        sequence_ids = []
        provenance = {f"sequences_{i+1}": [] for i in range(len(self.input_slots))}

        for combo in product(*axes):
            # combo alternates: (slot0_id, junc0_len, slot1_id, junc1_len, slot2_id, ...)
            seq_id = "_".join(str(part) for part in combo)
            sequence_ids.append(seq_id)

            # Track provenance — slot IDs are at even indices (0, 2, 4, ...)
            for slot_idx in range(len(self.input_slots)):
                combo_idx = slot_idx * 2  # Slots are at 0, 2, 4, ...
                provenance[f"sequences_{slot_idx+1}"].append(combo[combo_idx])

        return sequence_ids, provenance

    def validate_params(self):
        """Validate Fuse-specific parameters."""
        if not self.input_slots:
            raise ValueError("input_slots parameter is required")

        if len(self.input_slots) < 2:
            raise ValueError("Fuse requires at least 2 sequence slots")

        if not self.linker:
            raise ValueError("linker sequence is required")

        if not self.linker_lengths:
            raise ValueError("linker_lengths is required")

        # Validate linker_lengths format (basic check)
        for length_spec in self.linker_lengths:
            if not isinstance(length_spec, str) or not any(c in length_spec for c in '0123456789'):
                raise ValueError(f"Invalid linker_lengths format: {length_spec}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input sequences."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get Fuse configuration display lines."""
        config_lines = super().get_config_display()

        slot_summary = ", ".join(f"{len(s['ids'])} seqs" for s in self.input_slots)
        config_lines.extend([
            f"SLOTS: {len(self.input_slots)} ({slot_summary})",
            f"LINKER: {self.linker[:20]}{'...' if len(self.linker) > 20 else ''}",
            f"LINKER_LENGTHS: {', '.join(self.linker_lengths)}"
        ])

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """
        Generate Fuse execution script.

        Writes a JSON config file and calls pipe_fuse_queries.py with --config.
        """
        # Write JSON config for the pipe script
        config_data = {
            "name": self._get_job_base(),
            "linker": self.linker,
            "linker_lengths": self.linker_lengths,
            "slots": self.input_slots,
            "output_csv": self.queries_csv,
            "output_fasta": self.queries_fasta
        }

        os.makedirs(self.output_folder, exist_ok=True)
        with open(self.fuse_config_json, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# Fuse execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""echo "Generating fusion sequence combinations"
echo "Slots: {len(self.input_slots)}"
echo "Linker: {self.linker}"
echo "Linker lengths: {', '.join(self.linker_lengths)}"

# Call fuse_queries.py with JSON config
python {self.fuse_queries_py} --config "{self.fuse_config_json}"

# Check if files were created successfully
if [ ! -f "{self.queries_csv}" ]; then
    echo "ERROR: Failed to generate queries CSV file"
    exit 1
fi

if [ ! -f "{self.queries_fasta}" ]; then
    echo "ERROR: Failed to generate queries FASTA file"
    exit 1
fi

echo "Successfully generated fusion sequences:"
echo "CSV file: {self.queries_csv}"
echo "FASTA file: {self.queries_fasta}"

# Show summary
NUM_SEQUENCES=$(tail -n +2 "{self.queries_csv}" | wc -l)
echo "Generated $NUM_SEQUENCES fusion sequence combinations"

"""
        script_content += self.generate_completion_check_footer()

        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """
        Get expected output files after Fuse execution.

        Returns:
            Dictionary with DataStream objects and tables
        """
        # Predict sequence IDs and provenance
        sequence_ids, provenance = self._predict_sequence_ids()

        # Build position column names dynamically based on number of slots
        # For n slots: S1, L1, S2, L2, ..., S(n-1), L(n-1), Sn
        num_slots = len(self.input_slots)
        position_columns = []
        for i in range(1, num_slots + 1):
            position_columns.append(f"S{i}")
            if i < num_slots:
                position_columns.append(f"L{i}")

        # Build provenance column names
        provenance_columns = [f"sequences_{i+1}.id" for i in range(num_slots)]

        # Organize tables by content type
        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.queries_csv,
                columns=["id", "sequence", "lengths"] + provenance_columns + position_columns,
                description="Fusion sequences with provenance and sequence/linker positions",
                count=len(sequence_ids)
            )
        }

        # Create sequences DataStream
        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[self.queries_fasta],
            map_table=self.queries_csv,
            format="fasta"
        )

        return {
            "structures": DataStream.empty("structures", "pdb"),
            "sequences": sequences,
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration including Fuse-specific parameters."""
        base_dict = super().to_dict()
        base_dict.update({
            "fuse_params": {
                "name": self.name,
                "linker": self.linker,
                "linker_lengths": self.linker_lengths,
                "num_slots": len(self.input_slots),
                "slot_sizes": [len(s["ids"]) for s in self.input_slots]
            }
        })
        return base_dict
