# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""UniProt tool: fetch sequences and annotations from the UniProt REST API.

Takes a list of UniProt accessions, produces a `sequences` stream populated from
the canonical FASTA endpoint plus a tables/annotations.csv with GO terms,
domains, organism, length, etc. for downstream filtering or design workflows.
"""

import os
from typing import Dict, List, Any, Union

try:
    from .base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from .file_paths import Path
    from .datastream import DataStream
except ImportError:
    import sys
    sys.path.append(os.path.dirname(__file__))
    from base_config import BaseConfig, StandardizedOutput, TableInfo, _validate_freeform_string
    from file_paths import Path
    from datastream import DataStream


class UniProt(BaseConfig):
    """
    UniProt: retrieve sequences and annotations from the UniProt REST API.

    Inputs:
        accessions: list of UniProt accession strings (e.g. ["P12345", "Q9Y6K9"]).

    Outputs:
        Streams:
            sequences: id | sequence | type | length  (content-bearing CSV at sequences/sequences.csv)
        Tables:
            annotations: id | name | organism | taxonomy_id | length | go_terms | pfam | ec_number | reviewed
    """

    TOOL_NAME = "UniProt"
    TOOL_VERSION = "1.0"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== UniProt ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== UniProt ready ==="
"""

    # Content-bearing sequences stream: sequences.csv IS the content + map.
    sequences_csv = Path(lambda self: self.stream_path("sequences", "sequences.csv"))
    annotations_csv = Path(lambda self: self.table_path("annotations"))
    config_json = Path(lambda self: self.configuration_path("uniprot_config.json"))
    helper_script = Path(lambda self: self.pipe_script_path("pipe_uniprot.py"))

    def __init__(self,
                 accessions: Union[str, List[str]],
                 **kwargs):
        """
        Args:
            accessions: One UniProt accession or a list of them.
            **kwargs: Forwarded to BaseConfig.
        """
        if isinstance(accessions, str):
            accessions = [accessions]
        self.accessions: List[str] = list(accessions)
        super().__init__(**kwargs)

    def validate_params(self):
        if not self.accessions:
            raise ValueError("accessions is required and must not be empty")
        for i, acc in enumerate(self.accessions):
            if not isinstance(acc, str) or not acc:
                raise ValueError(f"accessions[{i}] must be a non-empty string")
            _validate_freeform_string(f"accessions[{i}]", acc)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        lines = super().get_config_display()
        lines.append(f"ACCESSIONS: {len(self.accessions)} ({', '.join(self.accessions[:5])}{'...' if len(self.accessions) > 5 else ''})")
        return lines

    def generate_script(self, script_path: str) -> str:
        import json
        os.makedirs(os.path.dirname(self.config_json), exist_ok=True)
        with open(self.config_json, "w", encoding="utf-8") as f:
            json.dump({
                "accessions": self.accessions,
            }, f, indent=2)

        script = "#!/bin/bash\n"
        script += "# UniProt fetch script\n"
        script += self.generate_completion_check_header()
        script += self.activate_environment()
        script += f"""echo "Fetching {len(self.accessions)} UniProt accession(s)"
python "{self.helper_script}" \\
    --config "{self.config_json}" \\
    --sequences-csv "{self.sequences_csv}" \\
    --annotations-csv "{self.annotations_csv}"

if [ $? -ne 0 ]; then
    echo "Error: UniProt fetch failed"
    exit 1
fi

"""
        script += self.generate_completion_check_footer()
        return script

    def get_output_files(self) -> Dict[str, Any]:
        # Content-bearing sequences stream: id pattern is the accession list itself.
        sequences = DataStream(
            name="sequences",
            ids=list(self.accessions),
            files=[],  # value-based stream
            map_table=self.sequences_csv,
            format="csv",
        )
        tables = {
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence", "type", "length"],
                description="UniProt sequences",
            ),
            "annotations": TableInfo(
                name="annotations",
                path=self.annotations_csv,
                columns=["id", "name", "organism", "taxonomy_id", "length",
                         "go_terms", "pfam", "ec_number", "reviewed"],
                description="UniProt annotations",
            ),
        }
        return {
            "sequences": sequences,
            "structures": DataStream.empty("structures", "pdb"),
            "compounds": DataStream.empty("compounds", "sdf"),
            "tables": tables,
            "output_folder": self.output_folder,
        }
