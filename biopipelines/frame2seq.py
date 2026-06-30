# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Frame2Seq tool for fast structure-conditioned sequence design.

Frame2seq is a non-autoregressive masked-language inverse-folding model that
generates N sequences per input backbone in a single forward pass. It serves
the same role as ProteinMPNN (structure -> sequence) but is materially faster
and reports a slightly higher native-sequence recovery on CATH 4.2.

Reference:
    Paper: https://arxiv.org/abs/2312.02447
    Repo:  https://github.com/dakpinaroglu/Frame2seq
"""

import os
from typing import Dict, List, Any, Tuple, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
    from .combinatorics import generate_multiplied_ids_pattern
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream
    from combinatorics import generate_multiplied_ids_pattern


def _selection_arg(value) -> str:
    """Serialize a fixed/redesigned selection for the pipe script.

    A plain selection string passes through; a TableReference (``table.column``)
    stringifies to a ``TABLE_REFERENCE:<path>:<column>`` token (mirrors
    ProteinMPNN). Empty selections become ``"-"``.
    """
    if not value:
        return "-"
    return str(value)


class Frame2Seq(BaseConfig):
    """
    Frame2Seq: structure-conditioned masked-language inverse folding.

    Predicts `num_sequences` candidate sequences for each input backbone in a
    single non-autoregressive pass. Output IDs follow the ProteinMPNN-style
    multiplier convention: <structure_id>_<n> for n in 1..num_sequences.

    Inputs:
        structures:      backbone PDBs (StandardizedOutput or DataStream)
        num_sequences:   number of sequences to sample per structure (>=1)
        temperature:     sampling temperature (>0)
        chain:           chain to redesign (default "A")
        omit_aa:         single-letter codes to exclude from sampling
                         (e.g. "CM" to omit Cys and Met). Empty disables.
        fixed:           residues to keep at their input identity. Like
                         ProteinMPNN/LigandMPNN, accepts either a chain-aware
                         selection string (e.g. "A10-20+A30") or a
                         (TableInfo, column) table reference. Empty disables.
        redesigned:      residues to redesign; the complement (all other
                         residues of `chain`) is held fixed. Same accepted
                         forms as `fixed`. Mutually exclusive with `fixed`.

    Outputs:
        Streams:
            sequences:  content-bearing CSV with id | sequence | score |
                        recovery | structures.id  (map_table doubles as data)
            fasta:      shared multi-record FASTA (all designed sequences)
        Tables:
            sequences:  same columns as the stream's map_table
            missing:    id | removed_by | cause (per-structure failures)
    """

    TOOL_NAME = "Frame2Seq"
    TOOL_VERSION = "1.1"

    # ------------------------------------------------------------------
    # Install
    # ------------------------------------------------------------------

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        """Create the frame2seq env from environments/frame2seq.<variant>.yaml
        and pip-install the package on top. Verification imports frame2seq +
        torch (the heaviest transitive dep)."""
        biopipelines = folders.get("biopipelines", "")
        env_check = cls._env_exists_check("frame2seq", env_manager)
        skip = "" if force_reinstall else f"""# Check if already installed
if {env_check} \\
   && {env_manager} run -n frame2seq python -c "import frame2seq, torch" >/dev/null 2>&1; then
    echo "frame2seq environment already installed, skipping. Use force_reinstall=True to reinstall."
    touch "$INSTALL_SUCCESS"
    exit 0
fi
"""
        remove_block = cls._env_remove_block("frame2seq", env_manager) if force_reinstall else ""
        env_block = cls._env_install_block("frame2seq", env_manager, biopipelines)
        return f"""echo "=== Installing Frame2Seq ==="
{skip}{remove_block}
{env_block}
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create frame2seq environment."
    exit 1
fi

# Verify installation (frame2seq + torch importable)
if {env_manager} run -n frame2seq python -c "import frame2seq, torch" >/dev/null 2>&1; then
    touch "$INSTALL_SUCCESS"
    echo "=== Frame2Seq installation complete ==="
else
    echo "ERROR: Frame2Seq verification failed (cannot import frame2seq or torch)"
    exit 1
fi
"""

    # ------------------------------------------------------------------
    # Lazy path descriptors — canonical sub-layout.
    #   _configuration/  .input_structures.json
    #   _execution/      raw frame2seq dumps (seqs/seqs.fasta + per-pdb csvs)
    #   sequences/       content-bearing: sequences.csv IS map_table + data
    #   tables/          standalone TableInfo (missing.csv)
    # ------------------------------------------------------------------

    structures_json = Path(lambda self: self.configuration_path(".input_structures.json"))
    frame2seq_out_folder = Path(lambda self: self.execution_folder)
    sequences_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))
    sequences_fasta = Path(lambda self: self.stream_path("sequences", "sequences.fasta"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    design_py = Path(lambda self: self.pipe_script_path("pipe_frame2seq.py"))

    # ------------------------------------------------------------------
    # Constructor
    # ------------------------------------------------------------------

    def __init__(
        self,
        structures: Union[DataStream, StandardizedOutput],
        num_sequences: int = 1,
        temperature: float = 1.0,
        chain: str = "A",
        omit_aa: str = "",
        fixed: Union[str, Tuple['TableInfo', str]] = "",
        redesigned: Union[str, Tuple['TableInfo', str]] = "",
        **kwargs,
    ):
        if isinstance(structures, StandardizedOutput):
            self.structures_stream: DataStream = structures.streams.structures
        elif isinstance(structures, DataStream):
            self.structures_stream = structures
        else:
            raise ValueError(
                f"structures must be DataStream or StandardizedOutput, got {type(structures).__name__}"
            )

        self.num_sequences = num_sequences
        self.temperature = temperature
        self.chain = chain
        self.omit_aa = omit_aa
        self.fixed = fixed
        self.redesigned = redesigned
        super().__init__(**kwargs)

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def validate_params(self):
        if not self.structures_stream or len(self.structures_stream) == 0:
            raise ValueError("structures parameter is required and must not be empty")
        if self.num_sequences < 1:
            raise ValueError("num_sequences must be at least 1")
        if self.temperature <= 0:
            raise ValueError("temperature must be positive")

        _validate_freeform_string("chain", self.chain)
        _validate_freeform_string("omit_aa", self.omit_aa)
        if isinstance(self.fixed, str):
            _validate_freeform_string("fixed", self.fixed)
        if isinstance(self.redesigned, str):
            _validate_freeform_string("redesigned", self.redesigned)
        if self.fixed and self.redesigned:
            raise ValueError("Specify either `fixed` or `redesigned`, not both")

        # omit_aa is a string of single-letter codes; reject anything else early.
        if self.omit_aa:
            valid = set("ACDEFGHIKLMNPQRSTVWY")
            bad = [c for c in self.omit_aa if c not in valid]
            if bad:
                raise ValueError(
                    f"omit_aa must contain only standard amino-acid letters; got invalid: {bad}"
                )

    # ------------------------------------------------------------------
    # Configure inputs
    # ------------------------------------------------------------------

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    # ------------------------------------------------------------------
    # Config display
    # ------------------------------------------------------------------

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.extend([
            f"INPUT STRUCTURES: {len(self.structures_stream)}",
            f"NUM SEQUENCES PER STRUCTURE: {self.num_sequences}",
            f"TEMPERATURE: {self.temperature}",
            f"CHAIN: {self.chain}",
        ])
        if self.omit_aa:
            lines.append(f"OMIT AA: {self.omit_aa}")
        if self.fixed:
            lines.append(f"FIXED: {self.fixed}")
        if self.redesigned:
            lines.append(f"REDESIGNED: {self.redesigned}")
        return lines

    # ------------------------------------------------------------------
    # Script generation
    # ------------------------------------------------------------------

    def generate_script(self, script_path: str) -> str:
        # Serialize DataStream so the pipe script can iterate at runtime.
        self.structures_stream.save_json(self.structures_json)

        script = "#!/bin/bash\n"
        script += "# Frame2Seq design script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()  # frame2seq env
        script += self.warn_container_unsupported()

        optional_flags = ""
        if self.omit_aa:
            optional_flags += f' --omit-aa "{self.omit_aa}"'

        # fixed/redesigned accept a selection string or a (TableInfo, column)
        # table reference, serialized exactly like ProteinMPNN. The pipe script
        # resolves them per structure (and computes the complement for
        # redesigned) against the input chain.
        fixed_arg = _selection_arg(self.fixed)
        redesigned_arg = _selection_arg(self.redesigned)

        script += f"""echo "Running Frame2Seq (chain={self.chain}, num_sequences={self.num_sequences}, T={self.temperature})"

python "{self.design_py}" \\
    --ds-json "{self.structures_json}" \\
    --output-folder "{self.frame2seq_out_folder}" \\
    --sequences-csv "{self.sequences_csv}" \\
    --sequences-fasta "{self.sequences_fasta}" \\
    --missing-csv "{self.missing_csv}" \\
    --num-sequences {self.num_sequences} \\
    --temperature {self.temperature} \\
    --chain "{self.chain}" \\
    --fixed "{fixed_arg}" \\
    --redesigned "{redesigned_arg}"{optional_flags}

if [ $? -ne 0 ]; then
    echo "Error: Frame2Seq design failed"
    exit 1
fi

"""
        script += self.generate_completion_check_footer()
        return script

    # ------------------------------------------------------------------
    # Output prediction
    # ------------------------------------------------------------------

    def get_output_files(self) -> Dict[str, Any]:
        # Multiplier IDs: <structure_id>_<1..num_sequences>
        suffix_pattern = f"<1..{self.num_sequences}>"
        sequence_ids = generate_multiplied_ids_pattern(
            self.structures_stream.ids,
            suffix_pattern,
            input_stream_name="structures",
        )

        # Content-bearing sequences stream: the CSV both lists IDs and holds
        # their sequences. Following the ProteinMPNN convention.
        sequences = DataStream(
            name="sequences",
            ids=sequence_ids,
            files=[],
            map_table=self.sequences_csv,
            format="csv",
        )

        # Shared-file FASTA twin so downstream tools that want FASTA can take
        # it as a single artifact (slice-by-id via the registered fasta slicer).
        fasta = DataStream(
            name="fasta",
            ids=sequence_ids,
            files=self.sequences_fasta,
            map_table=self.sequences_csv,
            format="fasta",
            metadata={
                "sequences_per_file": self.num_sequences,
                "contains_original": False,
                "source": "frame2seq_postprocessed",
            },
        )

        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "structures.id", "sequence", "score", "recovery"],
                description="Frame2Seq designed sequences (content-bearing map_table)",
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Structures Frame2Seq could not process, with reason",
            ),
        }

        return {
            "sequences": sequences,
            "fasta": fasta,
            "tables": tables,
            "output_folder": self.output_folder,
        }
