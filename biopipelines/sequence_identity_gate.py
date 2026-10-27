# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
SequenceIdentityGate: global sequence identity of each design against natural references.

A novelty gate for design campaigns: a design that copies a natural protein can
score well on every structural metric because the predictor may have memorized
it too, so this tool reports how close each design sits to the nearest natural
sequence. It mirrors the sequence identity gate of the RosettaCommons Chaperone
framework, applied to protein sequences.
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


class SequenceIdentityGate(BaseConfig):
    """
    SequenceIdentityGate: maximum global identity of each design to a set of
    natural reference sequences.

    For every design the tool aligns it against every reference with a global
    (Needleman-Wunsch) alignment and keeps the highest identity and the
    reference that gave it. Identity is the number of identical aligned
    positions divided by the alignment length, gaps included, so a design that
    shares a long stretch with a longer natural protein is not scored as if it
    were that protein. The alignment uses BLOSUM62 with a gap-open score of -10
    and a gap-extend score of -0.5 (the EMBOSS Needle defaults); end gaps are
    penalized. When several alignments tie on score, Biopython's first one
    decides.

    The tool never filters. `passed` states whether `max_identity` is strictly
    below `threshold`, and every design keeps its row, so the gate is applied
    downstream, as an explicit step::

        gate = SequenceIdentityGate(sequences=designs, references=Sequence("naturals.fasta"))
        novel = Panda(tables=gate.tables.identity_gate,
                      operations=[Panda.filter("passed == True")],
                      pool=designs)

    The whole sequence is scored. A motif copied unchanged from a natural
    protein therefore counts toward the identity, which raises `max_identity`
    for every design that carries it; set `threshold` for that campaign rather
    than reusing a value from another one.

    Inputs:
        sequences: protein sequences to score (StandardizedOutput or DataStream
            with a `sequences` stream). One row is one polymer chain, so a
            multi-chain design is scored chain by chain.
        references: natural protein sequences to compare against, as a
            `sequences` stream; `Sequence("references.fasta")` turns a FASTA
            file into one. Required and non-empty. A reference that cannot be
            scored stops the step rather than silently weakening the gate.
        threshold: maximum identity, as a fraction in [0, 1], a design may have
            and still pass (default 0.4). Set it per campaign.

    Outputs:
        Tables:
            identity_gate: id | max_identity | nearest_reference | passed
                `max_identity` is a fraction in [0, 1], `nearest_reference` the
                id of the reference that gave it (the first one on a tie), and
                `passed` is `max_identity < threshold`. A design that could not
                be scored carries NaN in all three.
            missing: id | removed_by | kind | cause

    Not covered: the identity of the whole sequence only. A motif-excluded
    score would need the motif residues mapped from structure numbering onto
    sequence positions, which no stream in the framework carries.
    """

    TOOL_NAME = "SequenceIdentityGate"
    TOOL_VERSION = "1.0"

    # No _install_script: runs in the biopipelines base env (Biopython is a base dependency).

    sequences_json = Path(lambda self: self.configuration_path("sequences.json"))
    references_json = Path(lambda self: self.configuration_path("references.json"))
    config_yaml = Path(lambda self: self.configuration_path("config.yaml"))
    identity_gate_csv = Path(lambda self: self.table_path("identity_gate"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_sequence_identity_gate.py"))

    def __init__(self,
                 sequences: Union[DataStream, StandardizedOutput],
                 references: Union[DataStream, StandardizedOutput],
                 threshold: float = 0.4,
                 **kwargs):
        self.sequences = sequences  # raw handle kept for missing-propagation
        self.sequences_stream = self._sequences_stream(sequences, "sequences")
        self.references_stream = self._sequences_stream(references, "references")
        self.threshold = float(threshold)
        super().__init__(**kwargs)

    @staticmethod
    def _sequences_stream(source, argument: str) -> DataStream:
        if isinstance(source, StandardizedOutput):
            return source.streams.sequences
        if isinstance(source, DataStream):
            return source
        raise ValueError(f"{argument} must be DataStream or StandardizedOutput, got {type(source)}")

    def validate_params(self):
        if not self.sequences_stream or len(self.sequences_stream) == 0:
            raise ValueError("sequences parameter is required and must not be empty")
        if not self.references_stream or len(self.references_stream) == 0:
            raise ValueError("references parameter is required and must not be empty")
        if not 0.0 <= self.threshold <= 1.0:
            raise ValueError(f"threshold must lie in [0, 1] as a fraction, got {self.threshold}")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"REFERENCES: {len(self.references_stream)} sequence(s)")
        lines.append(f"THRESHOLD: {self.threshold} (maximum identity to pass)")
        return lines

    def generate_script(self, script_path: str) -> str:
        self.sequences_stream.save_json(self.sequences_json)
        self.references_stream.save_json(self.references_json)
        self._write_config_yaml()
        script = "#!/bin/bash\n"
        script += "# SequenceIdentityGate script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Scoring {len(self.sequences_stream)} sequence(s) against {len(self.references_stream)} reference(s)"
python "{self.helper_py}" \\
    --config-yaml "{self.config_yaml}" \\
    --sequences-json "{self.sequences_json}" \\
    --references-json "{self.references_json}" \\
    --identity-gate-csv "{self.identity_gate_csv}" \\
    --local-missing-csv "{self.local_missing_csv}"
"""
        script += self.generate_missing_propagation(
            self.sequences, local_missing=self.local_missing_csv, missing_csv=self.missing_csv
        )
        script += self.generate_completion_check_footer()
        return script

    def _write_config_yaml(self):
        """Dump the parameters to configuration/config.yaml for the pipe to read."""
        import yaml
        with open(self.config_yaml, "w") as f:
            yaml.safe_dump({"threshold": self.threshold}, f, sort_keys=False)

    def get_output_files(self) -> Dict[str, Any]:
        tables = {
            "identity_gate": TableInfo(
                name="identity_gate",
                path=self.identity_gate_csv,
                columns=["id", "max_identity", "nearest_reference", "passed"],
                description="Maximum global identity of each design to the natural references and the gate verdict",
            ),
            "missing": self.missing_table_info(self.missing_csv),
        }
        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
