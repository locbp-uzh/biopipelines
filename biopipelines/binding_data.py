# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Measured protein-ligand binding affinities for a compounds stream.

Queries the ChEMBL and BindingDB web services for experimentally determined
affinity records (Ki, Kd, IC50, EC50) of each input compound against its
protein targets.

References:
    Zdrazil et al. (2024) The ChEMBL Database in 2023. Nucleic Acids Res 52, D1180.
    https://www.ebi.ac.uk/chembl/api/data
    Gilson et al. (2016) BindingDB in 2015. Nucleic Acids Res 44, D1045.
    https://bindingdb.org/rwd/bind/BindingDBRESTfulAPI.jsp
"""

import os
import json
from typing import Dict, List, Any, Optional, Union

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


CHEMBL_API_URL = "https://www.ebi.ac.uk/chembl/api/data"
BINDINGDB_API_URL = "https://bindingdb.org/rest"

SOURCES = ("chembl", "bindingdb")

# ChEMBL exposes both the raw and the standardized measurement; only the
# standardized one is unit-harmonised across assays, so it is what we read.
AFFINITY_TYPES = ("Ki", "Kd", "IC50", "EC50")

MATCH_TYPES = ("exact", "similarity", "substructure")


class BindingData(BaseConfig):
    """
    Experimental binding-affinity records for a compounds stream.

    For each input compound the tool resolves the molecule in each configured
    source and collects its measured affinities against protein targets. Every
    input id yields at least one row: a compound with no record in any source
    is reported with empty measurement fields rather than dropped, so the
    ``affinities`` table stays joinable against the input stream.

    Affinity values are reported in nM. ChEMBL records additionally carry
    ``pchembl_value`` (-log10 of the molar value), which is comparable across
    measurement types; BindingDB records carry the relation ('=', '>', '<')
    separately from the numeric value.

    Sources are queried independently and their rows concatenated, so the same
    measurement curated by both databases appears twice, distinguished by the
    ``source`` column. Deduplicate downstream with ``Panda`` if needed.

    Usage:
        with Pipeline(...):
            Resources(time="1:00:00")
            library = CompoundLibrary("my_library.csv")
            affinities = BindingData(compounds=library, max_affinity=1000)

            # Keep only sub-micromolar dissociation constants
            Panda(tables=affinities.tables.affinities,
                  operations=[Panda.filter("affinity_type == 'Kd'")])
    """

    TOOL_NAME = "BindingData"
    TOOL_VERSION = "1.2"

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== BindingData ==="
echo "Uses biopipelines environment (no additional installation needed)."
touch "$INSTALL_SUCCESS"
echo "=== BindingData ready ==="
"""

    # Lazy path descriptors — canonical sub-layout.
    affinities_csv = Path(lambda self: self.table_path("affinities"))
    targets_csv = Path(lambda self: self.table_path("targets"))
    missing_csv = Path(lambda self: self.table_path("missing"))
    local_missing_csv = Path(lambda self: self.execution_path("local_missing.csv"))
    config_json = Path(lambda self: self.configuration_path("binding_data_config.json"))
    compounds_json = Path(lambda self: self.configuration_path("compounds_ds.json"))
    helper_py = Path(lambda self: self.pipe_script_path("pipe_binding_data.py"))

    def __init__(self,
                 compounds: Union[DataStream, StandardizedOutput],
                 sources: Union[str, List[str]] = "chembl",
                 affinity_types: Optional[List[str]] = None,
                 match: str = "exact",
                 similarity: float = 0.85,
                 max_affinity: Optional[float] = None,
                 organism: Optional[str] = None,
                 max_records: int = 1000,
                 **kwargs):
        """
        Initialize the BindingData tool.

        Args:
            compounds: Input compounds as DataStream or StandardizedOutput.
                SMILES are read from the ``smiles`` column of the compounds
                map_table.
            sources: Which databases to query — "chembl" (default), "bindingdb",
                or a list of both. Rows carry a ``source`` column.
            affinity_types: Measurement types to keep, any of "Ki", "Kd",
                "IC50", "EC50". None (default) keeps all four.
            match: How a compound is matched to database molecules —
                "exact" (default, ChEMBL flexmatch / BindingDB similarity 1.0),
                "similarity" (Tanimoto above ``similarity``), or "substructure"
                (ChEMBL only).
            similarity: Tanimoto cutoff in [0, 1] for ``match="similarity"``
                (default 0.85). Ignored for the other match modes.
            max_affinity: Keep only records at or below this value in nM.
                None (default) keeps every record regardless of potency.
            organism: Restrict to targets of this source organism, matched
                against the record's organism field (e.g. "Homo sapiens").
                None (default) keeps every organism.
            max_records: Cap on the affinity records fetched per compound per
                source (default 1000).

        Output:
            Streams: (none)
            Tables:
                affinities: id | smiles | source | source_molecule_id | target_id |
                    target_name | organism | affinity_type | relation | affinity_nm |
                    pchembl_value | assay_id | assay_description | assay_type |
                    document_id | match_similarity
                targets: target_id | source | target_name | organism | n_compounds |
                    n_records | best_affinity_nm
                missing: id | removed_by | kind | cause
        """
        # Resolve compounds input to DataStream
        if isinstance(compounds, StandardizedOutput):
            self.compounds_stream: DataStream = compounds.streams.compounds
        elif isinstance(compounds, DataStream):
            self.compounds_stream = compounds
        else:
            raise ValueError(
                f"compounds must be DataStream or StandardizedOutput, got {type(compounds)}"
            )

        # Kept for missing-propagation, which walks the raw input handle.
        self.compounds = compounds

        self.sources = [sources] if isinstance(sources, str) else list(sources)
        self.sources = [s.lower() for s in self.sources]
        self.affinity_types = (list(AFFINITY_TYPES) if affinity_types is None
                               else list(affinity_types))
        self.match = match.lower()
        self.similarity = similarity
        self.max_affinity = max_affinity
        self.organism = organism
        self.max_records = max_records

        super().__init__(**kwargs)

    def validate_params(self):
        """Validate BindingData parameters."""
        if not self.compounds_stream or len(self.compounds_stream) == 0:
            raise ValueError("compounds parameter is required and must not be empty")

        if not self.sources:
            raise ValueError(f"sources must name at least one of {SOURCES}")
        for source in self.sources:
            if source not in SOURCES:
                raise ValueError(f"Invalid source: {source}. Must be one of {SOURCES}")

        if not self.affinity_types:
            raise ValueError(f"affinity_types must name at least one of {AFFINITY_TYPES}")
        for affinity_type in self.affinity_types:
            if affinity_type not in AFFINITY_TYPES:
                raise ValueError(
                    f"Invalid affinity_type: {affinity_type}. Must be one of {AFFINITY_TYPES}"
                )

        if self.match not in MATCH_TYPES:
            raise ValueError(f"Invalid match: {self.match}. Must be one of {MATCH_TYPES}")

        if self.match == "substructure" and "bindingdb" in self.sources:
            raise ValueError(
                "match='substructure' is supported by ChEMBL only; BindingDB's API "
                "offers exact and similarity matching. Drop 'bindingdb' from sources "
                "or use match='similarity'."
            )

        if not 0.0 <= self.similarity <= 1.0:
            raise ValueError(f"similarity must be between 0 and 1, got {self.similarity}")

        if self.max_affinity is not None and self.max_affinity <= 0:
            raise ValueError(f"max_affinity must be positive, got {self.max_affinity}")

        if self.max_records < 1:
            raise ValueError(f"max_records must be at least 1, got {self.max_records}")

        _validate_freeform_string("organism", self.organism)

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Configure folder paths."""
        self.folders = pipeline_folders

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()
        config_lines.extend([
            f"COMPOUNDS: {len(self.compounds_stream)} molecules",
            f"SOURCES: {', '.join(self.sources)}",
            f"MATCH: {self.match}" + (f" (Tanimoto >= {self.similarity})"
                                      if self.match == "similarity" else ""),
            f"AFFINITY_TYPES: {', '.join(self.affinity_types)}",
        ])
        if self.max_affinity is not None:
            config_lines.append(f"MAX_AFFINITY: {self.max_affinity} nM")
        if self.organism is not None:
            config_lines.append(f"ORGANISM: {self.organism}")
        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate BindingData execution script."""
        self.compounds_stream.save_json(self.compounds_json)

        config_data = {
            "compounds_json": str(self.compounds_json),
            "sources": self.sources,
            "affinity_types": self.affinity_types,
            "match": self.match,
            "similarity": self.similarity,
            "max_affinity": self.max_affinity,
            "organism": self.organism,
            "max_records": self.max_records,
            "affinities_csv": str(self.affinities_csv),
            "targets_csv": str(self.targets_csv),
            "local_missing_csv": str(self.local_missing_csv),
        }
        with open(self.config_json, 'w') as f:
            json.dump(config_data, f, indent=2)

        script_content = "#!/bin/bash\n"
        script_content += "# BindingData execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += f"""echo "Querying binding affinities"
echo "Compounds: {len(self.compounds_stream)}"
echo "Sources: {', '.join(self.sources)}"
echo "Output: {self.affinities_csv}"

python "{self.helper_py}" --config "{self.config_json}"

"""
        script_content += self.generate_missing_propagation(
            self.compounds,
            local_missing=self.local_missing_csv,
            missing_csv=self.missing_csv,
        )
        script_content += self.generate_completion_check_footer()
        return script_content

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output tables after the affinity lookup."""
        tables = {
            "affinities": TableInfo(
                name="affinities",
                path=self.affinities_csv,
                columns=["id", "smiles", "source", "source_molecule_id",
                         "target_id", "target_name", "organism",
                         "affinity_type", "relation", "affinity_nm",
                         "pchembl_value", "assay_id", "assay_description",
                         "assay_type", "document_id", "match_similarity",
                         "queried_utc", "chembl_release"],
                description="Measured binding affinities per compound-target pair"
            ),
            "targets": TableInfo(
                name="targets",
                path=self.targets_csv,
                columns=["target_id", "source", "target_name", "organism",
                         "n_compounds", "n_records", "best_affinity_nm",
                         "queried_utc", "chembl_release"],
                description="Protein targets aggregated over the returned records"
            ),
            "missing": TableInfo(
                name="missing",
                path=self.missing_csv,
                columns=["id", "removed_by", "kind", "cause"],
                description="Compounds with no affinity record, and upstream-filtered ids"
            ),
        }

        return {
            "tables": tables,
            "output_folder": self.output_folder,
        }
