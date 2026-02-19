# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
biopipelines Module

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
from .contacts import Contacts
from .distance import Distance
from .angle import Angle
from .stitch_sequences import StitchSequences
from .split_chains import SplitChains
from .dna_encoder import DNAEncoder
from .rbs_designer import RBSDesigner
from .load import Load, LoadMultiple
from .sequence_metric_correlation import SequenceMetricCorrelation
from .pymol import PyMOL
from .fuse import Fuse
from .plot import Plot
from .panda import Panda
from .mmseqs2 import MMseqs2, MMseqs2Server
from .remap import ReMap
from .gnina import Gnina
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
    'resolve_input_to_datastream',

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
    'Contacts',
    'Distance',
    'Angle',
    'StitchSequences',
    'SplitChains',
    'DNAEncoder',
    'RBSDesigner',
    'Load',
    'LoadMultiple',
    'SequenceMetricCorrelation',
    'PyMOL',
    'Fuse',
    'Plot',
    'Panda',
    'MMseqs2',
    'MMseqs2Server',
    'ReMap',
    'Gnina',

    # Utility functions
    'pdb_to_jsonl',
    'fasta_to_csv',
    'validate_pdb',
    'validate_fasta',
    'combine_fasta_files'
]