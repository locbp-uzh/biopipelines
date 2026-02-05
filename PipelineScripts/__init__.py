"""
PipelineScripts Module

Unified pipeline system for protein modeling workflows.
Provides tool configurations, pipeline orchestration, and seamless integration.
"""

from .pipeline import Pipeline, Bundle, Each
from .base_config import BaseConfig, ToolOutput
from .combinatorics import Bundle, Each, CombinatoricsConfig, generate_combinatorics_config, get_mode
from .datastream import DataStream, create_map_table
from .file_paths import Path
from .table_utils import get_table, get_table_path, list_tables, table_exists
from .rfdiffusion import RFdiffusion
from .rfdiffusion_allatom import RFdiffusionAllAtom, RFDAA_PrepareLigand
from .protein_mpnn import ProteinMPNN
from .alphafold import AlphaFold
from .ligand import Ligand
from .sequence import Sequence
from .ligand_mpnn import LigandMPNN
from .boltz2 import Boltz2
from .boltzgen import BoltzGen, BoltzGenMerge
from .compound_library import CompoundLibrary
from .extract_metrics import ExtractMetrics
from .folders import FolderManager
from .conformational_change import ConformationalChange
from .protein_ligand_contacts import ProteinLigandContacts
from .distance import Distance
from .angle import Angle
from .stitch_sequences import StitchSequences
from .split_chains import SplitChains
from .dynamic_bind import DynamicBind
from .dna_encoder import DNAEncoder
from .load import LoadOutput, LoadOutputs
from .sequence_metric_correlation import SequenceMetricCorrelation
from .pymol import PyMOL
from .fuse import Fuse
from .plot import Plot
from .panda import Panda
from .mmseqs2 import MMseqs2, MMseqs2Server
from .mmseqs2_lcf import MMseqs2LCF, MMseqs2ServerLCF
from .converters import *
from .entities import *

__all__ = [
    # Core pipeline components
    'Pipeline',
    'BaseConfig',
    'ToolOutput',
    'FolderManager',

    # DataStream I/O
    'DataStream',
    'create_map_table',
    'DataStreamResolver',
    'resolve_to_datastream',

    # File path descriptor
    'Path',

    # Table utilities
    'get_table',
    'get_table_path',
    'list_tables',
    'table_exists',

    # Combinatorics for input handling
    'Bundle',
    'Each',
    'CombinatoricsConfig',
    'generate_combinatorics_config',
    'get_mode',

    # Modeling tools
    'RFdiffusion',
    'RFdiffusionAllAtom',
    'RFDAA_PrepareLigand',
    'ProteinMPNN',
    'AlphaFold',
    'Ligand',
    'Sequence',
    'LigandMPNN',
    'Boltz2',
    'BoltzGen',
    'BoltzGenMerge',
    'CompoundLibrary',
    'ExtractMetrics',
    'ConformationalChange',
    'ProteinLigandContacts',
    'Distance',
    'Angle',
    'StitchSequences',
    'SplitChains',
    'DynamicBind',
    'DNAEncoder',
    'LoadOutput',
    'LoadOutputs',
    'SequenceMetricCorrelation',
    'PyMOL',
    'Fuse',
    'Plot',
    'Panda',
    'MMseqs2',
    'MMseqs2Server',
    'MMseqs2LCF',
    'MMseqs2ServerLCF',

    # Utility functions
    'pdb_to_jsonl',
    'fasta_to_csv',
    'validate_pdb',
    'validate_fasta',
    'combine_fasta_files'
]