# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
RCSB tool for searching the RCSB PDB Search API v2 and downloading matching structures.

Performs programmatic searches using text, sequence similarity, structure similarity,
chemical similarity, sequence motifs, and structure motifs. Downloads matching structures
with the same output format as the PDB tool.
"""

import os
import json
from enum import Enum
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


class RCSBQuery:
    """Represents a single RCSB search query node."""

    def __init__(self, query_dict: Dict[str, Any]):
        self.query_dict = query_dict


# ---------------------------------------------------------------------------
# Value enumerations (used as arguments to StructureAttribute methods)
# ---------------------------------------------------------------------------

class PolymerEntityType(Enum):
    """Polymer Entity Type — mirrors RCSB 'Polymer Entity Type' attribute values."""
    PROTEIN           = "Protein"
    DNA               = "DNA"
    RNA               = "RNA"
    NA_HYBRID         = "NA-hybrid"
    OTHER             = "Other"


class SymmetryType(Enum):
    """Symmetry Type — mirrors RCSB 'Symmetry Type' attribute values (Assembly Features)."""
    ASYMMETRIC        = "Asymmetric"
    HOMO_2_MER        = "Homo 2-mer"
    HOMO_3_MER        = "Homo 3-mer"
    HOMO_4_MER        = "Homo 4-mer"
    HOMO_5_MER        = "Homo 5-mer"
    HOMO_6_MER        = "Homo 6-mer"
    HOMO_7_MER        = "Homo 7-mer"
    PSEUDO_SYMMETRIC  = "Pseudo-Symmetric"
    HETEROMERIC       = "Heteromeric"


class ExperimentalMethod(Enum):
    """Experimental Method — mirrors RCSB 'Experimental Method' attribute values."""
    X_RAY_DIFFRACTION      = "X-RAY DIFFRACTION"
    ELECTRON_MICROSCOPY    = "ELECTRON MICROSCOPY"
    SOLUTION_NMR           = "SOLUTION NMR"
    SOLID_STATE_NMR        = "SOLID-STATE NMR"
    NEUTRON_DIFFRACTION    = "NEUTRON DIFFRACTION"
    ELECTRON_CRYSTALLOGRAPHY = "ELECTRON CRYSTALLOGRAPHY"
    FIBER_DIFFRACTION      = "FIBER DIFFRACTION"
    EPR                    = "EPR"
    FLUORESCENCE_TRANSFER  = "FLUORESCENCE TRANSFER"
    INFRARED_SPECTROSCOPY  = "INFRARED SPECTROSCOPY"


# ---------------------------------------------------------------------------
# Attribute descriptor — holds the API path, produces RCSBQuery via methods
# ---------------------------------------------------------------------------

class _Attr:
    """
    Represents a single searchable RCSB attribute.

    Methods produce RCSBQuery objects that can be passed directly to RCSB().
    Each method name mirrors the kind of comparison being performed, making
    intent explicit at the call site.

    Example:
        RCSB(
            StructureAttribute.NonpolymerFeatures.ComponentIdentifier.is_any("FMN", "FAD"),
            StructureAttribute.PolymerMolecularFeatures.PolymerEntityType.equals(PolymerEntityType.PROTEIN),
            StructureAttribute.AssemblyFeatures.SymmetryType.equals(SymmetryType.ASYMMETRIC),
            StructureAttribute.PolymerMolecularFeatures.PolymerEntitySequenceLength.between(80, 300),
        )
    """

    def __init__(self, path: str):
        self._path = path

    def _make(self, operator: str, value: Any, negation: bool = False) -> RCSBQuery:
        params: Dict[str, Any] = {
            "attribute": self._path,
            "operator": operator,
            "value": value,
        }
        if negation:
            params["negation"] = True
        return RCSBQuery({
            "type": "terminal",
            "service": "text",
            "parameters": params,
        })

    # --- comparison methods ---

    def equals(self, value: Any) -> RCSBQuery:
        """Exact equality. Pass an Enum member or a plain value."""
        v = value.value if isinstance(value, Enum) else value
        return self._make("exact_match", v)

    def is_any(self, *values) -> RCSBQuery:
        """Value is one of the supplied options (OR). Pass Enum members or plain strings."""
        vs = [v.value if isinstance(v, Enum) else v for v in values]
        return self._make("in", vs)

    def contains(self, text: str) -> RCSBQuery:
        """Value contains the given words (case-insensitive substring search)."""
        return self._make("contains_words", text)

    def contains_phrase(self, text: str) -> RCSBQuery:
        """Value contains the given exact phrase."""
        return self._make("contains_phrase", text)

    def between(self, low: Any, high: Any, inclusive: bool = True) -> RCSBQuery:
        """Numeric or date range. Both bounds included by default."""
        return self._make("range", {
            "from": low,
            "to": high,
            "include_lower": inclusive,
            "include_upper": inclusive,
        })

    def greater_than(self, value: Any, or_equal: bool = False) -> RCSBQuery:
        """Greater than (or equal) comparison."""
        op = "greater_or_equal" if or_equal else "greater"
        return self._make(op, value)

    def less_than(self, value: Any, or_equal: bool = False) -> RCSBQuery:
        """Less than (or equal) comparison."""
        op = "less_or_equal" if or_equal else "less"
        return self._make(op, value)

    def exists(self) -> RCSBQuery:
        """Match entries where this attribute exists (is not null)."""
        return self._make("exists", None)

    def not_equals(self, value: Any) -> RCSBQuery:
        """Negated exact equality."""
        v = value.value if isinstance(value, Enum) else value
        return self._make("exact_match", v, negation=True)

    def not_any(self, *values) -> RCSBQuery:
        """Value is none of the supplied options."""
        vs = [v.value if isinstance(v, Enum) else v for v in values]
        return self._make("in", vs, negation=True)


# ---------------------------------------------------------------------------
# StructureAttribute — namespace mirroring RCSB advanced search categories
# ---------------------------------------------------------------------------

class StructureAttribute:
    """
    Typed attribute descriptors for RCSB advanced search, organised by the
    same categories shown at rcsb.org/search/advanced.

    Each attribute exposes comparison methods that return an RCSBQuery:
        .equals(value)
        .is_any(*values)
        .contains(text)
        .contains_phrase(text)
        .between(low, high)
        .greater_than(value, or_equal=False)
        .less_than(value, or_equal=False)
        .exists()
        .not_equals(value)
        .not_any(*values)

    Usage:
        from biopipelines.rcsb import StructureAttribute, PolymerEntityType, SymmetryType

        RCSB(
            StructureAttribute.NonpolymerFeatures.ComponentIdentifier.is_any("FMN"),
            StructureAttribute.PolymerMolecularFeatures.PolymerEntityType.equals(PolymerEntityType.PROTEIN),
            StructureAttribute.AssemblyFeatures.SymmetryType.equals(SymmetryType.ASYMMETRIC),
            StructureAttribute.PolymerMolecularFeatures.PolymerEntitySequenceLength.between(80, 300),
            max_results=100,
        )
    """

    # --- ID(s) and Keywords ---
    class IDsAndKeywords:
        EntryId                     = _Attr("rcsb_entry_container_identifiers.entry_id")
        StructureKeywords           = _Attr("struct_keywords.pdbx_keywords")
        AdditionalKeywords          = _Attr("struct_keywords.text")

    # --- Structure Details ---
    class StructureDetails:
        StructureTitle              = _Attr("struct.title")
        DepositDate                 = _Attr("rcsb_accession_info.deposit_date")
        ReleaseDate                 = _Attr("rcsb_accession_info.initial_release_date")
        RevisionDate                = _Attr("rcsb_accession_info.major_revision")

    # --- Entry Features ---
    class EntryFeatures:
        MolecularWeight             = _Attr("rcsb_entry_info.molecular_weight")
        NumberOfProteinEntities     = _Attr("rcsb_entry_info.polymer_entity_count_protein")
        NumberOfRnaEntities         = _Attr("rcsb_entry_info.polymer_entity_count_RNA")
        NumberOfDnaEntities         = _Attr("rcsb_entry_info.polymer_entity_count_DNA")
        TotalPolymerResidues        = _Attr("rcsb_entry_info.deposited_polymer_monomer_count")
        EntryPolymerComposition     = _Attr("rcsb_entry_info.polymer_composition")
        EntryPolymerTypes           = _Attr("rcsb_entry_info.selected_polymer_entity_types")

    # --- Polymer Molecular Features ---
    class PolymerMolecularFeatures:
        PolymerEntityDescription    = _Attr("rcsb_polymer_entity.pdbx_description")
        PolymerEntityType           = _Attr("entity_poly.rcsb_entity_polymer_type")
        PolymerEntitySequenceLength = _Attr("entity_poly.rcsb_sample_sequence_length")
        MonomerComponentId          = _Attr("entity_poly.pdbx_seq_one_letter_code_can")
        ScientificName              = _Attr("rcsb_entity_source_organism.scientific_name")
        TaxonomyId                  = _Attr("rcsb_entity_source_organism.ncbi_taxonomy_id")
        GeneName                    = _Attr("rcsb_gene_name.value")
        EnzymeClassification        = _Attr("rcsb_polymer_entity.rcsb_ec_lineage.id")

    # --- Nonpolymer (Ligand) Molecular Features ---
    class NonpolymerFeatures:
        ComponentIdentifier         = _Attr("rcsb_nonpolymer_entity_container_identifiers.comp_id")
        FormulaWeight               = _Attr("rcsb_nonpolymer_entity.formula_weight")
        LigandQscore                = _Attr("rcsb_nonpolymer_entity_instance_container_identifiers.ligand_qscore")

    # --- Assembly Features ---
    class AssemblyFeatures:
        SymmetryType                = _Attr("rcsb_struct_symmetry.type")
        SymmetrySymbol              = _Attr("rcsb_struct_symmetry.symbol")
        OligomericState             = _Attr("pdbx_struct_assembly.oligomeric_details")

    # --- Methods ---
    class Methods:
        ExperimentalMethod          = _Attr("rcsb_entry_info.experimental_method")
        Resolution                  = _Attr("rcsb_entry_info.resolution_combined")


class RCSB(BaseConfig):
    """
    Pipeline tool for searching the RCSB PDB Search API and downloading matching structures.

    Performs a search at configuration time to identify matching PDB entries, then
    downloads the structures at execution time. Outputs the same streams as the PDB tool.

    Example:
        # Text search
        results = RCSB(
            RCSB.Text("insulin receptor"),
            max_results=10
        )

        # Attribute search with resolution filter
        results = RCSB(
            RCSB.Attribute("rcsb_entity_source_organism.scientific_name", "exact_match", "Homo sapiens"),
            RCSB.Attribute("rcsb_entry_info.resolution_combined", "less", 2.0),
            max_results=20
        )

        # Sequence similarity search
        results = RCSB(
            RCSB.Sequence("MKTVRQERLKSIVRILERSKEPVSGAQ", identity_cutoff=0.9),
            max_results=50
        )

        # Structure similarity
        results = RCSB(
            RCSB.Structure("4HHB", assembly_id=1),
            max_results=10
        )

        # Sequence motif (PROSITE format)
        results = RCSB(
            RCSB.SeqMotif("C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H", pattern_type="prosite"),
            max_results=100
        )

        # Chemical similarity (by SMILES)
        results = RCSB(
            RCSB.Chemical("c1ccc(cc1)C(=O)O", match_type="fingerprint-similarity"),
            max_results=10
        )

        # Structure motif
        results = RCSB(
            RCSB.StrucMotif([
                {"label_comp_id": "HIS", "label_asym_id": "A", "label_seq_id": 94},
                {"label_comp_id": "HIS", "label_asym_id": "A", "label_seq_id": 96},
                {"label_comp_id": "HIS", "label_asym_id": "A", "label_seq_id": 119},
            ], pdb_id="4HHB"),
            max_results=10
        )

        # Combined queries (AND by default)
        results = RCSB(
            RCSB.Text("kinase"),
            RCSB.Attribute("rcsb_entry_info.resolution_combined", "less", 2.5),
            max_results=50,
            sort="resolution"
        )

        # Use downstream like PDB output
        af = AlphaFold(proteins=results)
    """

    TOOL_NAME = "RCSB"

    SEARCH_API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
    GRAPHQL_URL = "https://data.rcsb.org/graphql"

    # Sort field mapping for convenience names
    SORT_FIELDS = {
        "score": "score",
        "resolution": "rcsb_entry_info.resolution_combined",
        "release_date": "rcsb_entry_info.release_date",
    }

    @classmethod
    def _install_script(cls, folders, env_manager="mamba", force_reinstall=False, **kwargs):
        return """echo "=== RCSB ==="
echo "Uses biopipelines environment (no additional installation needed)."
echo "=== RCSB ready ==="
"""

    # Lazy path descriptors
    structures_csv = Path(lambda self: os.path.join(self.output_folder, "structures.csv"))
    sequences_csv = Path(lambda self: os.path.join(self.output_folder, "sequences.csv"))
    compounds_csv = Path(lambda self: os.path.join(self.output_folder, "compounds.csv"))
    failed_csv = Path(lambda self: os.path.join(self.output_folder, "failed_downloads.csv"))
    search_results_csv = Path(lambda self: os.path.join(self.output_folder, "search_results.csv"))
    entry_info_csv = Path(lambda self: os.path.join(self.output_folder, "entry_info.csv"))
    config_file = Path(lambda self: os.path.join(self.output_folder, "fetch_config.json"))
    pdb_py = Path(lambda self: os.path.join(self.folders["HelpScripts"], "pipe_pdb.py"))

    # --- Static methods for creating query objects ---

    @staticmethod
    def Text(value: str) -> RCSBQuery:
        """
        Full-text search across all fields.

        Supports boolean operators: + (AND), | (OR), - (NOT), quoted phrases.

        Args:
            value: Search text (e.g., "insulin receptor", "kinase + human")

        Returns:
            RCSBQuery for full-text search
        """
        return RCSBQuery({
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": value
            }
        })

    @staticmethod
    def Attribute(attribute: str, operator: str, value: Any, negation: bool = False) -> RCSBQuery:
        """
        Search by a specific RCSB attribute.

        Args:
            attribute: RCSB attribute path (e.g., "rcsb_entry_info.resolution_combined",
                      "rcsb_entity_source_organism.scientific_name")
            operator: Comparison operator. One of:
                      "exact_match", "contains_words", "contains_phrase",
                      "greater", "less", "greater_or_equal", "less_or_equal",
                      "equals", "range", "exists", "in"
            value: Value to compare against. For "range", use dict with
                   "from", "to", "include_lower", "include_upper" keys.
                   For "in", use a list of values.
            negation: If True, invert the operator logic (default: False)

        Returns:
            RCSBQuery for attribute search

        Example:
            RCSB.Attribute("rcsb_entry_info.resolution_combined", "less", 2.0)
            RCSB.Attribute("rcsb_entity_source_organism.scientific_name", "exact_match", "Homo sapiens")
            RCSB.Attribute("rcsb_entry_info.experimental_method", "exact_match", "X-RAY DIFFRACTION")
            RCSB.Attribute("rcsb_entry_info.release_date", "range", {"from": "2023-01-01", "to": "2024-01-01", "include_lower": True, "include_upper": False})
        """
        params = {
            "attribute": attribute,
            "operator": operator,
            "value": value
        }
        if negation:
            params["negation"] = True

        return RCSBQuery({
            "type": "terminal",
            "service": "text",
            "parameters": params
        })

    @staticmethod
    def Sequence(sequence: str, identity_cutoff: float = 0.9,
                 evalue_cutoff: float = 0.1,
                 sequence_type: str = "protein") -> RCSBQuery:
        """
        BLAST-like sequence similarity search.

        Args:
            sequence: Query sequence (amino acid or nucleotide)
            identity_cutoff: Minimum sequence identity (0.0-1.0, default: 0.9)
            evalue_cutoff: Maximum E-value threshold (default: 0.1)
            sequence_type: "protein" (default), "dna", or "rna"

        Returns:
            RCSBQuery for sequence similarity search
        """
        target = {
            "protein": "pdb_protein_sequence",
            "dna": "pdb_dna_sequence",
            "rna": "pdb_rna_sequence"
        }.get(sequence_type, "pdb_protein_sequence")

        return RCSBQuery({
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": evalue_cutoff,
                "identity_cutoff": identity_cutoff,
                "target": target,
                "value": sequence
            }
        })

    @staticmethod
    def SeqMotif(pattern: str, pattern_type: str = "prosite",
                 sequence_type: str = "protein") -> RCSBQuery:
        """
        Sequence motif search using PROSITE, simple, or regex patterns.

        Args:
            pattern: Motif pattern string
            pattern_type: "prosite" (default), "simple", or "regex"
            sequence_type: "protein" (default), "dna", or "rna"

        Returns:
            RCSBQuery for sequence motif search

        Example:
            RCSB.SeqMotif("C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H", pattern_type="prosite")
        """
        target = {
            "protein": "pdb_protein_sequence",
            "dna": "pdb_dna_sequence",
            "rna": "pdb_rna_sequence"
        }.get(sequence_type, "pdb_protein_sequence")

        return RCSBQuery({
            "type": "terminal",
            "service": "seqmotif",
            "parameters": {
                "value": pattern,
                "pattern_type": pattern_type,
                "target": target
            }
        })

    @staticmethod
    def Structure(entry_id: str, assembly_id: int = 1) -> RCSBQuery:
        """
        3D structure similarity search.

        Args:
            entry_id: PDB ID of the query structure (e.g., "4HHB")
            assembly_id: Assembly ID to use (default: 1)

        Returns:
            RCSBQuery for structure similarity search
        """
        return RCSBQuery({
            "type": "terminal",
            "service": "structure",
            "parameters": {
                "value": {
                    "entry_id": entry_id.upper(),
                    "assembly_id": str(assembly_id)
                },
                "operator": "strict_shape_match"
            }
        })

    @staticmethod
    def StrucMotif(residues: List[Dict[str, Any]], pdb_id: str,
                   assembly_id: int = 1) -> RCSBQuery:
        """
        Structure motif search using geometric residue arrangements.

        Args:
            residues: List of residue descriptors, each a dict with keys:
                      "label_comp_id" (residue name), "label_asym_id" (chain),
                      "label_seq_id" (sequence position)
            pdb_id: PDB ID containing the motif
            assembly_id: Assembly ID (default: 1)

        Returns:
            RCSBQuery for structure motif search

        Example:
            RCSB.StrucMotif([
                {"label_comp_id": "HIS", "label_asym_id": "A", "label_seq_id": 94},
                {"label_comp_id": "HIS", "label_asym_id": "A", "label_seq_id": 96},
                {"label_comp_id": "HIS", "label_asym_id": "A", "label_seq_id": 119},
            ], pdb_id="4HHB")
        """
        return RCSBQuery({
            "type": "terminal",
            "service": "strucmotif",
            "parameters": {
                "value": {
                    "entry_id": pdb_id.upper(),
                    "assembly_id": str(assembly_id),
                    "residue_ids": residues
                }
            }
        })

    @staticmethod
    def Chemical(value: str, match_type: str = "graph-relaxed",
                 descriptor_type: str = "SMILES") -> RCSBQuery:
        """
        Chemical similarity search by SMILES or InChI descriptor.

        For searching by CCD code (e.g., "ATP"), use RCSB.Attribute instead:
            RCSB.Attribute("rcsb_chem_comp_container_identifiers.comp_id", "in", ["ATP"])

        Args:
            value: SMILES or InChI string (e.g., "c1ccc(cc1)C(=O)O")
            match_type: Matching method. Options:
                        "graph-exact", "graph-strict", "graph-relaxed",
                        "graph-relaxed-stereo", "fingerprint-similarity",
                        "sub-struct-graph-exact", "sub-struct-graph-strict",
                        "sub-struct-graph-relaxed", "sub-struct-graph-relaxed-stereo"
                        (default: "graph-relaxed")
            descriptor_type: Type of the value - "SMILES" (default) or "InChI"

        Returns:
            RCSBQuery for chemical similarity search

        Example:
            RCSB.Chemical("c1ccc(cc1)C(=O)O", match_type="fingerprint-similarity")
            RCSB.Chemical("c1ccccc1", match_type="graph-relaxed")
        """
        return RCSBQuery({
            "type": "terminal",
            "service": "chemical",
            "parameters": {
                "value": value,
                "type": "descriptor",
                "descriptor_type": descriptor_type,
                "match_type": match_type
            }
        })

    # --- Instance methods ---

    def __init__(self,
                 *queries: RCSBQuery,
                 max_results: int = 10,
                 return_type: str = "entry",
                 sort: str = "score",
                 format: str = "pdb",
                 ids: Optional[Union[str, List[str]]] = None,
                 remove_waters: bool = True,
                 chain: str = "longest",
                 logical_operator: str = "and",
                 **kwargs):
        """
        Initialize RCSB search tool.

        Args:
            *queries: One or more RCSBQuery objects (from RCSB.Text, RCSB.Attribute, etc.).
                      Multiple queries are combined with the logical_operator.
            max_results: Maximum number of PDB entries to return (default: 10)
            return_type: Result granularity - "entry" (default), "assembly",
                        "polymer_entity", "polymer_instance"
            sort: Sort field - "score" (default), "resolution", "release_date",
                  or any RCSB attribute path
            format: Output structure format - "pdb" (default) or "cif"
            ids: Optional custom IDs. If None, uses PDB IDs from search results.
            remove_waters: Whether to remove water molecules (default: True)
            chain: Which chain to extract sequence from - "longest" (default) or chain letter
            logical_operator: How to combine multiple queries - "and" (default) or "or"
            **kwargs: Additional parameters

        Output:
            Streams: structures (.pdb/.cif), sequences (.csv), compounds (.csv)
            Tables:
                structures: id | pdb_id | file_path | format | file_size | source
                sequences: id | sequence
                compounds: id | code | format | smiles | ccd
                search_results: id | pdb_id | score
                failed: pdb_id | error_message | source | attempted_path
        """
        # Validate queries
        if not queries:
            raise ValueError("At least one query is required (e.g., RCSB.Text('kinase'))")

        self.queries = list(queries)
        for i, q in enumerate(self.queries):
            if not isinstance(q, RCSBQuery):
                raise ValueError(
                    f"Argument {i + 1} is not an RCSB query. "
                    f"Use RCSB.Text(), RCSB.Attribute(), RCSB.Sequence(), etc."
                )

        self.max_results = max_results
        self.return_type = return_type
        self.format = format.lower()
        self.remove_waters = remove_waters
        self.chain = chain
        self.logical_operator = logical_operator.lower()
        self.custom_ids = None

        # Handle sort
        if sort in self.SORT_FIELDS:
            self.sort_field = self.SORT_FIELDS[sort]
        else:
            self.sort_field = sort

        # Handle custom IDs
        if ids is not None:
            if isinstance(ids, str):
                self.custom_ids = [ids]
            else:
                self.custom_ids = list(ids)

        # Validate parameters
        if self.return_type not in ["entry", "assembly", "polymer_entity", "polymer_instance"]:
            raise ValueError(f"Invalid return_type: {self.return_type}")

        if self.format not in ["pdb", "cif"]:
            raise ValueError(f"Invalid format: {self.format}. Must be 'pdb' or 'cif'")

        if self.logical_operator not in ["and", "or"]:
            raise ValueError(f"Invalid logical_operator: {logical_operator}. Must be 'and' or 'or'")

        # Will be populated during configure_inputs
        self.pdb_ids = []
        self.search_scores = []
        self.result_ids = []

        super().__init__(**kwargs)

    def _build_query_json(self) -> Dict[str, Any]:
        """Build the full RCSB Search API query JSON."""
        # Add node_id to each terminal query (required by the API)
        nodes = []
        for i, q in enumerate(self.queries):
            node = dict(q.query_dict)
            node["node_id"] = i
            nodes.append(node)

        # Build query node
        if len(nodes) == 1:
            query_node = nodes[0]
        else:
            query_node = {
                "type": "group",
                "logical_operator": self.logical_operator,
                "nodes": nodes
            }

        # Build request
        request = {
            "query": query_node,
            "return_type": self.return_type,
            "request_options": {
                "paginate": {
                    "start": 0,
                    "rows": self.max_results
                }
            }
        }

        # Add sorting
        if self.sort_field != "score":
            request["request_options"]["sort"] = [{
                "sort_by": "attribute",
                "attribute": self.sort_field,
                "direction": "asc"
            }]

        return request

    def _extract_pdb_id(self, identifier: str) -> str:
        """Extract the PDB ID from a search result identifier.

        Different return_types produce different identifier formats:
        - entry: "4HHB"
        - assembly: "4HHB-1"
        - polymer_entity: "4HHB_1"
        - polymer_instance: "4HHB.A"
        """
        if self.return_type == "entry":
            return identifier
        elif self.return_type == "assembly":
            return identifier.split("-")[0]
        elif self.return_type == "polymer_entity":
            return identifier.split("_")[0]
        elif self.return_type == "polymer_instance":
            return identifier.split(".")[0]
        return identifier

    def _perform_search(self) -> List[Dict[str, Any]]:
        """
        Execute the RCSB search at configuration time.

        Returns:
            List of result dicts with 'identifier' and 'score' keys

        Raises:
            ValueError: If search fails or returns no results
        """
        import requests

        query_json = self._build_query_json()

        try:
            response = requests.post(
                self.SEARCH_API_URL,
                json=query_json,
                headers={"Content-Type": "application/json"},
                timeout=30
            )
        except requests.exceptions.RequestException as e:
            raise ValueError(f"RCSB search request failed: {e}")

        if response.status_code == 204:
            raise ValueError("RCSB search returned no results")

        if response.status_code != 200:
            raise ValueError(
                f"RCSB search failed with status {response.status_code}: "
                f"{response.text[:200]}"
            )

        data = response.json()
        results = data.get("result_set", [])
        total = data.get("total_count", 0)

        if not results:
            raise ValueError("RCSB search returned no results")

        print(f"  RCSB search: {total} total matches, returning top {len(results)}")
        return results

    def _fetch_entry_metadata(self, pdb_ids: List[str]) -> List[Dict[str, Any]]:
        """
        Fetch entry metadata from RCSB Data API via GraphQL.

        Retrieves title, resolution, method, organism, molecular weight,
        citation info, and release date for each PDB entry.

        Args:
            pdb_ids: List of PDB IDs to fetch metadata for

        Returns:
            List of metadata dicts, one per PDB ID (in same order)
        """
        import requests

        # GraphQL query for batch metadata
        ids_str = ", ".join(f'"{pid}"' for pid in pdb_ids)
        query = f"""
        {{
          entries(entry_ids: [{ids_str}]) {{
            rcsb_id
            struct {{
              title
            }}
            rcsb_entry_info {{
              resolution_combined
              experimental_method
              molecular_weight
              polymer_entity_count_protein
              deposited_polymer_monomer_count
            }}
            rcsb_primary_citation {{
              title
              rcsb_journal_abbrev
              year
              rcsb_authors
            }}
            rcsb_accession_info {{
              deposit_date
              initial_release_date
            }}
            polymer_entities {{
              rcsb_polymer_entity {{
                pdbx_description
              }}
              rcsb_entity_source_organism {{
                scientific_name
                ncbi_taxonomy_id
              }}
            }}
          }}
        }}
        """

        try:
            resp = requests.post(self.GRAPHQL_URL, json={"query": query}, timeout=15)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            print(f"  Warning: Could not fetch entry metadata: {e}")
            return [{} for _ in pdb_ids]

        entries_data = (data.get("data") or {}).get("entries") or []

        # Build lookup by PDB ID
        entry_lookup = {}
        for entry in entries_data:
            if entry is None:
                continue
            entry_lookup[entry.get("rcsb_id", "")] = entry

        # Build metadata in same order as pdb_ids
        metadata_list = []
        for pdb_id in pdb_ids:
            entry = entry_lookup.get(pdb_id, {})
            if not entry:
                metadata_list.append({})
                continue

            info = entry.get("rcsb_entry_info") or {}
            struct = entry.get("struct") or {}
            cit = entry.get("rcsb_primary_citation") or {}
            acc = entry.get("rcsb_accession_info") or {}

            # Collect organisms and entity descriptions
            organisms = []
            descriptions = []
            for pe in (entry.get("polymer_entities") or []):
                desc = (pe.get("rcsb_polymer_entity") or {}).get("pdbx_description", "")
                if desc and desc not in descriptions:
                    descriptions.append(desc)
                for org in (pe.get("rcsb_entity_source_organism") or []):
                    name = org.get("scientific_name", "")
                    if name and name not in organisms:
                        organisms.append(name)

            resolution = info.get("resolution_combined", [None])
            resolution_val = resolution[0] if isinstance(resolution, list) and resolution else resolution

            authors = cit.get("rcsb_authors") or []
            release_date = (acc.get("initial_release_date") or "")[:10]  # YYYY-MM-DD
            deposit_date = (acc.get("deposit_date") or "")[:10]

            metadata_list.append({
                "title": struct.get("title", ""),
                "resolution": resolution_val,
                "method": info.get("experimental_method", ""),
                "molecular_weight_kda": info.get("molecular_weight", ""),
                "protein_entity_count": info.get("polymer_entity_count_protein", ""),
                "residue_count": info.get("deposited_polymer_monomer_count", ""),
                "organism": "; ".join(organisms),
                "entity_description": "; ".join(descriptions),
                "citation_title": cit.get("title", ""),
                "citation_journal": cit.get("rcsb_journal_abbrev", ""),
                "citation_year": cit.get("year", ""),
                "citation_authors": "; ".join(authors[:5]) + ("..." if len(authors) > 5 else ""),
                "release_date": release_date,
                "deposit_date": deposit_date,
            })

        return metadata_list

    def validate_params(self):
        """Validate RCSB parameters."""
        if not self.queries:
            raise ValueError("At least one query is required")

        if self.max_results < 1:
            raise ValueError("max_results must be at least 1")

        if self.format not in ["pdb", "cif"]:
            raise ValueError("format must be 'pdb' or 'cif'")

    def configure_inputs(self, pipeline_folders: Dict[str, str]):
        """Perform RCSB search and configure inputs for download."""
        self.folders = pipeline_folders

        # Perform search
        results = self._perform_search()

        # Extract PDB IDs and scores
        seen_pdb_ids = set()
        for result in results:
            identifier = result.get("identifier", "")
            score = result.get("score", 0.0)
            pdb_id = self._extract_pdb_id(identifier)

            # Deduplicate PDB IDs (different entities/chains may map to same entry)
            if pdb_id not in seen_pdb_ids:
                seen_pdb_ids.add(pdb_id)
                self.pdb_ids.append(pdb_id)
                self.search_scores.append(score)
                self.result_ids.append(identifier)

        # Apply custom IDs if provided
        if self.custom_ids is not None:
            if len(self.custom_ids) != len(self.pdb_ids):
                raise ValueError(
                    f"ids has {len(self.custom_ids)} items but search returned "
                    f"{len(self.pdb_ids)} unique PDB entries"
                )
            self.output_ids = self.custom_ids
        else:
            self.output_ids = [pdb_id.lower() for pdb_id in self.pdb_ids]

        # Fetch entry metadata from RCSB Data API
        self.entry_metadata = self._fetch_entry_metadata(self.pdb_ids)

        # Print summary with metadata
        for pdb_id, score, meta in zip(self.pdb_ids, self.search_scores, self.entry_metadata):
            title = meta.get("title", "")
            if len(title) > 60:
                title = title[:57] + "..."
            organism = meta.get("organism", "")
            resolution = meta.get("resolution", "")
            res_str = f", {resolution}A" if resolution else ""
            org_str = f", {organism}" if organism else ""
            print(f"  {pdb_id}: {title}{res_str}{org_str} (score: {score:.2f})")

    def get_config_display(self) -> List[str]:
        """Get configuration display lines."""
        config_lines = super().get_config_display()

        # Summarize queries
        query_summaries = []
        for q in self.queries:
            service = q.query_dict.get("service", "unknown")
            if service == "full_text":
                query_summaries.append(f"Text(\"{q.query_dict['parameters']['value']}\")")
            elif service == "text":
                attr = q.query_dict["parameters"].get("attribute", "")
                op = q.query_dict["parameters"].get("operator", "")
                val = q.query_dict["parameters"].get("value", "")
                query_summaries.append(f"Attribute({attr}, {op}, {val})")
            elif service == "sequence":
                seq = q.query_dict["parameters"].get("value", "")[:20]
                query_summaries.append(f"Sequence({seq}...)")
            elif service == "structure":
                entry = q.query_dict["parameters"].get("value", {}).get("entry_id", "")
                query_summaries.append(f"Structure({entry})")
            elif service == "seqmotif":
                pattern = q.query_dict["parameters"].get("value", "")[:30]
                query_summaries.append(f"SeqMotif({pattern}...)")
            elif service == "strucmotif":
                entry = q.query_dict["parameters"].get("value", {}).get("entry_id", "")
                query_summaries.append(f"StrucMotif({entry})")
            elif service == "chemical":
                val = q.query_dict["parameters"].get("value", "")[:20]
                query_summaries.append(f"Chemical({val})")
            else:
                query_summaries.append(f"{service}(...)")

        config_lines.extend([
            f"QUERIES: {f' {self.logical_operator.upper()} '.join(query_summaries)}",
            f"MAX_RESULTS: {self.max_results}",
            f"RETURN_TYPE: {self.return_type}",
            f"SORT: {self.sort_field}",
            f"FORMAT: {self.format.upper()}",
        ])

        if self.pdb_ids:
            config_lines.append(f"FOUND: {', '.join(self.pdb_ids)} ({len(self.pdb_ids)} entries)")

        return config_lines

    def generate_script(self, script_path: str) -> str:
        """Generate RCSB execution script."""
        os.makedirs(self.output_folder, exist_ok=True)

        script_content = "#!/bin/bash\n"
        script_content += "# RCSB search + download execution script\n"
        script_content += self.generate_completion_check_header()
        script_content += self.activate_environment()
        script_content += self._generate_script_body()
        script_content += self.generate_completion_check_footer()

        return script_content

    def _generate_script_body(self) -> str:
        """Generate the download execution part of the script."""
        repo_pdbs_folder = self.folders['PDBs']

        # Build config for pipe_pdb.py (same format it expects)
        config_data = {
            "pdb_ids": self.pdb_ids,
            "custom_ids": self.output_ids,
            "format": self.format,
            "local_folder": None,
            "repo_pdbs_folder": repo_pdbs_folder,
            "biological_assembly": False,
            "remove_waters": self.remove_waters,
            "chain": self.chain,
            "output_folder": self.output_folder,
            "structures_table": self.structures_csv,
            "sequences_table": self.sequences_csv,
            "failed_table": self.failed_csv,
            "compounds_table": self.compounds_csv,
            "operations": []
        }

        with open(self.config_file, 'w') as f:
            json.dump(config_data, f, indent=2)

        # Write search results CSV
        self._write_search_results_csv()

        return f"""echo "RCSB Search: downloading {len(self.pdb_ids)} structures"
echo "PDB IDs: {', '.join(self.pdb_ids)}"
echo "Output IDs: {', '.join(self.output_ids)}"
echo "Format: {self.format.upper()}"
echo "Output folder: {self.output_folder}"

python "{self.pdb_py}" --config "{self.config_file}"

"""

    def _write_search_results_csv(self):
        """Write search results CSV with scores and entry info CSV with metadata."""
        import csv
        os.makedirs(self.output_folder, exist_ok=True)

        # Write search_results.csv (scores)
        with open(self.search_results_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["id", "pdb_id", "result_id", "score"])
            for output_id, pdb_id, result_id, score in zip(
                self.output_ids, self.pdb_ids, self.result_ids, self.search_scores
            ):
                writer.writerow([output_id, pdb_id, result_id, score])

        # Write entry_info.csv (metadata from RCSB Data API)
        info_columns = [
            "id", "pdb_id", "title", "resolution", "method",
            "molecular_weight_kda", "organism", "entity_description",
            "protein_entity_count", "residue_count",
            "citation_title", "citation_journal", "citation_year", "citation_authors",
            "release_date", "deposit_date"
        ]
        with open(self.entry_info_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(info_columns)
            for output_id, pdb_id, meta in zip(
                self.output_ids, self.pdb_ids, self.entry_metadata
            ):
                writer.writerow([
                    output_id,
                    pdb_id,
                    meta.get("title", ""),
                    meta.get("resolution", ""),
                    meta.get("method", ""),
                    meta.get("molecular_weight_kda", ""),
                    meta.get("organism", ""),
                    meta.get("entity_description", ""),
                    meta.get("protein_entity_count", ""),
                    meta.get("residue_count", ""),
                    meta.get("citation_title", ""),
                    meta.get("citation_journal", ""),
                    meta.get("citation_year", ""),
                    meta.get("citation_authors", ""),
                    meta.get("release_date", ""),
                    meta.get("deposit_date", ""),
                ])

    def get_output_files(self) -> Dict[str, Any]:
        """Get expected output files after search and download."""
        extension = ".pdb" if self.format == "pdb" else ".cif"
        structure_files = [os.path.join(self.output_folder, f"{oid}{extension}")
                          for oid in self.output_ids]

        tables = {
            "structures": TableInfo(
                name="structures",
                path=self.structures_csv,
                columns=["id", "pdb_id", "file_path", "format", "file_size", "source"],
                description="Successfully fetched structure files from RCSB search",
                count=len(self.pdb_ids)
            ),
            "sequences": TableInfo(
                name="sequences",
                path=self.sequences_csv,
                columns=["id", "sequence"],
                description="Protein sequences extracted from structures",
                count=len(self.pdb_ids)
            ),
            "compounds": TableInfo(
                name="compounds",
                path=self.compounds_csv,
                columns=["id", "code", "format", "smiles", "ccd"],
                description="Ligands extracted from PDB structures",
                count="variable"
            ),
            "search_results": TableInfo(
                name="search_results",
                path=self.search_results_csv,
                columns=["id", "pdb_id", "result_id", "score"],
                description="RCSB search results with scores",
                count=len(self.pdb_ids)
            ),
            "entry_info": TableInfo(
                name="entry_info",
                path=self.entry_info_csv,
                columns=["id", "pdb_id", "title", "resolution", "method",
                         "molecular_weight_kda", "organism", "entity_description",
                         "protein_entity_count", "residue_count",
                         "citation_title", "citation_journal", "citation_year",
                         "citation_authors", "release_date", "deposit_date"],
                description="Entry metadata from RCSB (title, organism, citation, etc.)",
                count=len(self.pdb_ids)
            ),
            "failed": TableInfo(
                name="failed",
                path=self.failed_csv,
                columns=["pdb_id", "error_message", "source", "attempted_path"],
                description="Failed structure fetches with error details",
                count="variable"
            )
        }

        structures = DataStream(
            name="structures",
            ids=list(self.output_ids),
            files=structure_files,
            map_table=self.structures_csv,
            format=self.format
        )

        sequences = DataStream(
            name="sequences",
            ids=list(self.output_ids),
            files=[self.sequences_csv],
            map_table=self.sequences_csv,
            format="csv"
        )

        compounds = DataStream(
            name="compounds",
            ids=[],
            files=[],
            map_table=self.compounds_csv,
            format="csv"
        )

        return {
            "structures": structures,
            "sequences": sequences,
            "compounds": compounds,
            "tables": tables,
            "output_folder": self.output_folder
        }

    def to_dict(self) -> Dict[str, Any]:
        """Serialize configuration."""
        base_dict = super().to_dict()
        base_dict.update({
            "tool_params": {
                "queries": [q.query_dict for q in self.queries],
                "max_results": self.max_results,
                "return_type": self.return_type,
                "sort_field": self.sort_field,
                "format": self.format,
                "pdb_ids": self.pdb_ids,
                "output_ids": self.output_ids if hasattr(self, 'output_ids') else [],
                "search_scores": self.search_scores,
                "remove_waters": self.remove_waters,
                "chain": self.chain,
                "logical_operator": self.logical_operator
            }
        })
        return base_dict
