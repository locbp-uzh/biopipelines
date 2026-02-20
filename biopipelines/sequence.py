# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Sequence tool for creating sequence DataStreams from raw sequence strings, CSV/FASTA
files, or RCSB PDB codes.

Provides a standardized way to pass sequences to other tools like Boltz2, AlphaFold, etc.
Supports automatic type detection (protein/DNA/RNA) or explicit type specification.
"""

import os
import re
import pandas as pd
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


class Sequence(BaseConfig):
    """
    Pipeline tool for creating sequence DataStreams from raw sequence strings, CSV/FASTA
    files, or RCSB PDB codes.

    Converts raw sequence strings, CSV/FASTA files, or sequences fetched from RCSB
    into a standardized format that can be passed to other tools like Boltz2,
    AlphaFold, ProteinMPNN, etc.
    """

    TOOL_NAME = "Sequence"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== Sequence ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== Sequence ready ==="
"""

    # Lazy path descriptors
    sequences_csv = Path(lambda self: os.path.join(self.output_folder, "sequences.csv"))
    sequences_fasta = Path(lambda self: os.path.join(self.output_folder, "sequences.fasta"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "sequence_config.json"))
    sequence_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_sequence.py"))

    def __init__(self,
                 seq: Union[str, List[str]],
                 type: str = "auto",
                 ids: Optional[Union[str, List[str]]] = None,
                 **kwargs):
        """
        Initialize Sequence tool.

        Args:
            seq: Sequence input. Can be:
                 - Single sequence string: "MKTVRQERLKSIVRILERSKEPVSGAQ"
                 - Multiple sequences: ["MKTVRQ...", "AETGFT..."]
                 - Path to CSV file with 'id' and 'sequence' columns
                   (absolute, relative, or filename inside the Sequences/ folder)
                 - Path to a FASTA file (.fasta or .fa)
                   (absolute, relative, or filename inside the Sequences/ folder)
                 - 4-character RCSB PDB code: "4EQ7"
                   Fetches all protein chains from the RCSB FASTA API.
            type: Sequence type - "auto", "protein", "dna", or "rna".
                  "auto" will detect based on sequence content:
                  - Contains only ACGT -> DNA
                  - Contains only ACGU -> RNA
                  - Otherwise -> protein
            ids: Output identifier(s). Ignored when loading from a file or PDB code.
                 If not provided for raw sequences, defaults to "seq_1", "seq_2", etc.
            **kwargs: Additional parameters

        Examples:
            # Single protein sequence
            seq = Sequence("MKTVRQERLKSIVRILERSKEPVSGAQ")

            # With custom id
            seq = Sequence("MKTVRQ...", ids="my_protein")

            # From RCSB PDB code (fetches longest chain per entity)
            seq = Sequence("4EQ7")

            # Load from CSV (absolute path, or filename inside Sequences/)
            seq = Sequence("my_proteins.csv")
            seq = Sequence("/absolute/path/to/sequences.csv")

            # Load from FASTA (absolute path, or filename inside Sequences/)
            seq = Sequence("my_proteins.fasta")

            # Explicit type
            seq = Sequence("ACGTACGT", type="dna")
        """
        # Track source
        self.from_csv = False
        self.from_fasta = False
        self.from_pdb_code = False
        self.from_excel = False
        self._excel_source_path = None
        self.source_csv_path = None
        # Pending filename: set when a bare filename is given that needs resolution
        # against the Sequences/ folder in configure_inputs
        self._pending_filename = None

        if isinstance(seq, str):
            # --- RCSB PDB code ---
            if self._is_pdb_code(seq):
                self._load_from_pdb_code(seq, ids)
            # --- CSV/FASTA path or filename (may need Sequences/ resolution) ---
            elif self._has_sequence_extension(seq):
                if os.path.isfile(seq):
                    # Exists as given — load immediately
                    if ids is not None:
                        print(f"  Warning: 'ids' parameter ignored when loading from file (using ids from file)")
                    lower = seq.lower()
                    if lower.endswith('.xlsx') or lower.endswith('.xls'):
                        self._load_from_excel(seq)
                    elif lower.endswith('.csv'):
                        self._load_from_csv(seq)
                    else:
                        self._load_from_fasta(seq)
                else:
                    # Defer to configure_inputs to try joining with Sequences/ folder
                    self._pending_filename = seq
                    self.sequences = []
                    self.custom_ids = []
                    self._pending_ids_param = ids
            # --- Raw sequence string ---
            else:
                self.sequences = [seq]
                self._apply_ids(ids, 1)
        else:
            # List of raw sequence strings
            self.sequences = list(seq)
            if not self.sequences:
                raise ValueError("Must provide at least one sequence")
            for i, s in enumerate(self.sequences):
                if not s or not s.strip():
                    raise ValueError(f"Sequence at index {i} is empty")
            self._apply_ids(ids, len(self.sequences))

        # Validate and store type
        valid_types = ["auto", "protein", "dna", "rna"]
        if type not in valid_types:
            raise ValueError(f"Invalid type: {type}. Must be one of: {valid_types}")
        self.seq_type = type

        # Detect types (may be empty if pending; re-detected in configure_inputs)
        self.detected_types = [
            self._detect_type(s) if self.seq_type == "auto" else self.seq_type
            for s in self.sequences
        ]

        # Initialize base class
        super().__init__(**kwargs)

    # ------------------------------------------------------------------
    # Helpers: source detection
    # ------------------------------------------------------------------

    def _is_pdb_code(self, s: str) -> bool:
        """Return True if s looks like a 4-character alphanumeric RCSB PDB code."""
        return len(s) == 4 and s.isalnum() and not os.path.exists(s)

    def _has_sequence_extension(self, s: str) -> bool:
        """Return True if s has a CSV, Excel, or FASTA file extension."""
        lower = s.lower()
        return (lower.endswith('.csv') or lower.endswith('.xlsx') or lower.endswith('.xls')
                or lower.endswith('.fasta') or lower.endswith('.fa'))

    def _apply_ids(self, ids, count: int):
        """Set self.custom_ids from the ids parameter or generate defaults."""
        if ids is not None:
            if isinstance(ids, str):
                self.custom_ids = [ids]
            else:
                self.custom_ids = list(ids)
            if len(self.custom_ids) != count:
                raise ValueError(
                    f"Length mismatch: ids has {len(self.custom_ids)} items "
                    f"but sequences has {count} items"
                )
        else:
            self.custom_ids = [f"seq_{i + 1}" for i in range(count)]

    # ------------------------------------------------------------------
    # Loaders
    # ------------------------------------------------------------------

    def _load_from_pdb_code(self, pdb_code: str, ids):
        """Fetch protein chain sequences from the RCSB FASTA API."""
        try:
            import requests
        except ImportError:
            raise ImportError("'requests' module is required to fetch sequences from RCSB")

        pdb_id = pdb_code.upper()
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"
        print(f"  Fetching sequences for {pdb_id} from RCSB FASTA: {url}")

        response = requests.get(url, timeout=15,
                                headers={"User-Agent": "BioPipelines/1.0"})
        response.raise_for_status()

        fasta_text = response.text
        chains = self._parse_fasta_text(fasta_text)

        if not chains:
            raise ValueError(f"No protein sequences found in RCSB FASTA for {pdb_id}")

        self.sequences = [seq for _, seq in chains]
        if ids is not None:
            if isinstance(ids, str):
                self.custom_ids = [ids]
            else:
                self.custom_ids = list(ids)
            if len(self.custom_ids) != len(self.sequences):
                raise ValueError(
                    f"ids has {len(self.custom_ids)} items but RCSB returned "
                    f"{len(self.sequences)} chain(s) for {pdb_id}"
                )
        else:
            self.custom_ids = [header for header, _ in chains]

        self.from_pdb_code = True
        print(f"  Loaded {len(self.sequences)} chain(s) from RCSB {pdb_id}: "
              f"{', '.join(self.custom_ids)}")

    def _parse_fasta_text(self, fasta_text: str) -> List[tuple]:
        """
        Parse FASTA text into (id, sequence) pairs.
        Uses the first token of each header line as the id.
        """
        chains = []
        current_id = None
        current_seq_parts = []

        for line in fasta_text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seq = "".join(current_seq_parts)
                    if seq:
                        chains.append((current_id, seq))
                # Use first token (e.g. "4EQ7_1|Chain A" -> "4EQ7_1")
                header = line[1:].split("|")[0].strip()
                current_id = header
                current_seq_parts = []
            else:
                current_seq_parts.append(line)

        if current_id is not None:
            seq = "".join(current_seq_parts)
            if seq:
                chains.append((current_id, seq))

        return chains

    def _load_from_csv(self, csv_path: str):
        """Load sequences from a CSV file with 'id' and 'sequence' columns."""
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"CSV file not found: {csv_path}")

        self.source_csv_path = os.path.abspath(csv_path)
        self.from_csv = True

        df = pd.read_csv(csv_path)

        if 'id' not in df.columns:
            raise ValueError(f"CSV file must have 'id' column. Found: {list(df.columns)}")
        if 'sequence' not in df.columns:
            raise ValueError(f"CSV file must have 'sequence' column. Found: {list(df.columns)}")

        self.custom_ids = df['id'].astype(str).tolist()
        self.sequences = df['sequence'].astype(str).tolist()

        for i, (seq_id, seq) in enumerate(zip(self.custom_ids, self.sequences)):
            if not seq or not seq.strip():
                raise ValueError(f"Empty sequence for id '{seq_id}' at row {i}")

        print(f"  Loaded {len(self.sequences)} sequences from CSV: {csv_path}")

    def _load_from_excel(self, excel_path: str):
        """Load sequences from an Excel file with 'id' and 'sequence' columns."""
        if not os.path.exists(excel_path):
            raise FileNotFoundError(f"Excel file not found: {excel_path}")

        df = pd.read_excel(excel_path)

        if 'id' not in df.columns:
            raise ValueError(f"Excel file must have 'id' column. Found: {list(df.columns)}")
        if 'sequence' not in df.columns:
            raise ValueError(f"Excel file must have 'sequence' column. Found: {list(df.columns)}")

        self.custom_ids = df['id'].astype(str).tolist()
        self.sequences = df['sequence'].astype(str).tolist()

        for i, (seq_id, seq) in enumerate(zip(self.custom_ids, self.sequences)):
            if not seq or not seq.strip():
                raise ValueError(f"Empty sequence for id '{seq_id}' at row {i}")

        self.from_excel = True
        self._excel_source_path = os.path.abspath(excel_path)
        # source_csv_path is set in configure_inputs once output_folder is available
        print(f"  Loaded {len(self.sequences)} sequences from Excel: {excel_path}")

    def _load_from_fasta(self, fasta_path: str):
        """Load sequences from a FASTA file."""
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        with open(fasta_path, "r") as f:
            chains = self._parse_fasta_text(f.read())

        if not chains:
            raise ValueError(f"No sequences found in FASTA file: {fasta_path}")

        self.custom_ids = [h for h, _ in chains]
        self.sequences = [s for _, s in chains]
        self.from_fasta = True
        self.source_csv_path = os.path.abspath(fasta_path)
        print(f"  Loaded {len(self.sequences)} sequences from FASTA: {fasta_path}")

    # ------------------------------------------------------------------
    # Type detection
    # ------------------------------------------------------------------

    def _detect_type(self, seq: str) -> str:
        """
        Detect sequence type based on content.

        Args:
            seq: Sequence string

        Returns:
            "protein", "dna", or "rna"
        """
        # Clean sequence (remove whitespace, convert to uppercase)
        clean_seq = re.sub(r'\s', '', seq.upper())

        # Check for DNA (only ACGT)
        if re.match(r'^[ACGT]+$', clean_seq):
            return "dna"

        # Check for RNA (only ACGU)
        if re.match(r'^[ACGU]+$', clean_seq):
            return "rna"

        # Default to protein
        return "protein"

    def validate_params(self):
        """Validate Sequence parameters."""
        # Pending-filename sequences are validated after configure_inputs resolves them
        if self._pending_filename:
            return

        if not self.sequences:
            raise ValueError("sequences cannot be empty")

        if not self.custom_ids:
            raise ValueError("ids cannot be empty")

        if len(self.custom_ids) != len(self.sequences):
            raise ValueError(f"ids length ({len(self.custom_ids)}) must match sequences length ({len(self.sequences)})")

        # Validate sequence characters based on detected type
        valid_chars = {
            "protein": set("ACDEFGHIKLMNPQRSTVWY"),
            "dna": set("ACGT"),
            "rna": set("ACGU"),
        }
        for seq_id, seq, seq_type in zip(self.custom_ids, self.sequences, self.detected_types):
            allowed = valid_chars[seq_type]
            clean = seq.upper().replace(" ", "")
            invalid = set(clean) - allowed
            if invalid:
                raise ValueError(
                    f"Sequence '{seq_id}' contains invalid characters for {seq_type}: "
                    f"{''.join(sorted(invalid))}. "
                    f"Allowed characters: {''.join(sorted(allowed))}"
                )

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure input parameters."""
        self.folders = pipeline_folders

        # Resolve a pending path against the Sequences/ folder
        if self._pending_filename:
            sequences_folder = pipeline_folders.get("Sequences", "")
            candidate = os.path.join(sequences_folder, self._pending_filename)
            if not os.path.isfile(candidate):
                raise FileNotFoundError(
                    f"Sequence file not found: '{self._pending_filename}'. "
                    f"Tried as-is and as '{candidate}'."
                )
            if self._pending_ids_param is not None:
                print(f"  Warning: 'ids' parameter ignored when loading from file (using ids from file)")
            lower = self._pending_filename.lower()
            if lower.endswith('.xlsx') or lower.endswith('.xls'):
                self._load_from_excel(candidate)
            elif lower.endswith('.csv'):
                self._load_from_csv(candidate)
            else:
                self._load_from_fasta(candidate)
            self._pending_filename = None
            # Re-detect types now that sequences are loaded
            self.detected_types = [
                self._detect_type(s) if self.seq_type == "auto" else self.seq_type
                for s in self.sequences
            ]

        # For Excel inputs: write the converted CSV to sequences_csv (the lazy path)
        # so that source_csv_path and map_table both point to a real CSV file.
        if self.from_excel:
            os.makedirs(self.output_folder, exist_ok=True)
            df = pd.DataFrame({"id": self.custom_ids, "sequence": self.sequences})
            df.to_csv(self.sequences_csv, index=False)
            self.source_csv_path = self.sequences_csv
            print(f"  Converted Excel to CSV: {self.sequences_csv}")

        # Display sequence info
        for seq_id, seq, seq_type in zip(self.custom_ids, self.sequences, self.detected_types):
            seq_preview = seq[:30] + "..." if len(seq) > 30 else seq
            print(f"  Sequence {seq_id}: {seq_preview} (type: {seq_type}, length: {len(seq)})")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        # Show source
        if self.from_pdb_code:
            config_lines.append(f"SOURCE: RCSB PDB")
        elif self.from_excel and self._excel_source_path:
            config_lines.append(f"SOURCE (Excel): {self._excel_source_path}")
        elif (self.from_csv or self.from_fasta) and self.source_csv_path:
            config_lines.append(f"SOURCE: {self.source_csv_path}")
        elif self._pending_filename:
            config_lines.append(f"SOURCE: {self._pending_filename} (Sequences/ folder)")

        config_lines.extend([
            f"IDS: {', '.join(self.custom_ids[:5])}{'...' if len(self.custom_ids) > 5 else ''} ({len(self.custom_ids)} sequences)",
            f"TYPE: {self.seq_type}",
        ])

        # Show detected types if auto
        if self.seq_type == "auto":
            type_counts = {}
            for t in self.detected_types:
                type_counts[t] = type_counts.get(t, 0) + 1
            type_summary = ", ".join([f"{t}: {c}" for t, c in type_counts.items()])
            config_lines.append(f"DETECTED: {type_summary}")

        # Show sequence previews
        for seq_id, seq in zip(self.custom_ids[:3], self.sequences[:3]):
            seq_preview = seq[:20] + "..." if len(seq) > 20 else seq
            config_lines.append(f"  {seq_id}: {seq_preview}")
        if len(self.sequences) > 3:
            config_lines.append(f"  ... and {len(self.sequences) - 3} more")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate script to create sequence files."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# Sequence creation script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self.generate_script_run_sequence()
        script_content += self.generate_completion_check_footer()

        return script_content

    def generate_script_run_sequence(self) -> str:
        """Generate the sequence file creation part of the script."""
        import json

        config_data = {
            "custom_ids": self.custom_ids,
            "sequences": self.sequences,
            "types": self.detected_types,
            "output_folder": self.output_folder,
            "sequences_csv": self.sequences_csv,
            "sequences_fasta": self.sequences_fasta
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        return f"""echo "Creating sequence files for {len(self.sequences)} sequences"
echo "IDs: {', '.join(self.custom_ids)}"
echo "Output folder: {self.output_folder}"

python "{self.sequence_py}" --config "{self.config_file}"

"""

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after sequence file creation."""
        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence", "type", "length"],
                description="Sequence data with metadata",
                count=len(self.sequences)
            )
        }

        # sequences stream: value-based, sequence data lives in CSV map_table
        sequences = DataStream(
            name="sequences",
            ids=self.custom_ids.copy(),
            files=[],
            map_table=self.sequences_csv,
            format="csv"
        )

        # fasta stream: file-based, single .fasta file shared across all IDs
        fasta = DataStream(
            name="fasta",
            ids=self.custom_ids.copy(),
            files=[self.sequences_fasta],
            format="fasta"
        )

        return {
            "sequences": sequences,
            "fasta": fasta,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "custom_ids": self.custom_ids,
                "sequences": self.sequences,
                "seq_type": self.seq_type,
                "detected_types": self.detected_types,
                "from_csv": self.from_csv,
                "from_fasta": self.from_fasta,
                "from_pdb_code": self.from_pdb_code,
                "from_excel": self.from_excel,
                "excel_source_path": self._excel_source_path,
                "source_csv_path": self.source_csv_path
            }
        })
        return base_dict
