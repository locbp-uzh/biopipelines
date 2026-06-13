# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
biopipelines Module

Unified pipeline system for protein modeling workflows.
Provides tool configurations, pipeline orchestration, and seamless integration.
"""

__version__ = "1.2.0"

from .pipeline import Pipeline, Bundle, Each, Folder
from .base_config import BaseConfig, ToolOutput, IndexedTableContainer
from .combinatorics import Bundle, Each, CombinatoricsConfig, generate_combinatorics_config, get_mode
from .datastream import DataStream, create_map_table
from .file_paths import Path
from .table_utils import get_table, get_table_path, list_tables, table_exists, get_indexed_table
from .rfdiffusion import RFdiffusion
from .rfdiffusion_allatom import RFdiffusionAllAtom, RFDAA_PrepareLigand
from .rfdiffusion2 import RFdiffusion2
from .rfdiffusion3 import RFdiffusion3
from .hbdesigner import HBDesigner
from .protein_mpnn import ProteinMPNN
from .frame2seq import Frame2Seq
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
from .dna_encoder import DNAEncoder
from .rbs_designer import RBSDesigner
from .load import Load, LoadMultiple
from .sequence_metric_correlation import SequenceMetricCorrelation
from .pymol import PyMOL
from .fuse import Fuse
from .plot import Plot
from .panda import Panda
from .msa import MSA
from .mmseqs2 import MMseqs2, MMseqs2Server
from .remap import ReMap
from .mock import Mock
from .gnina import Gnina
from .diffdock import DiffDock
from .pocketgen import PocketGen
from .placer import PLACER
from .dynamicbind import DynamicBind
from .neuralplexer import NeuralPLexer
from .posebusters import PoseBusters
from .admet_ai import ADMETAI
from .esmfold import ESMFold
from .uniprot import UniProt
from .rdkit_descriptors import RDKit
from .openbabel import OpenBabel
from .plip import PLIP
from .prodigy import Prodigy
from .prolif import ProLIF
from .fpocket import FPocket
from .p2rank import P2Rank
from .af2bind import AF2BIND
from .dssp import DSSP
from .apbs import APBS
from .openmm import OpenMM
from .reduce import Reduce
from .xtb import XTB
from .rtmscore import RTMScore
from .gems import GEMS
from .bioemu import BioEmu
from .cabsflex import CABSflex
from .ensemble_analysis import EnsembleAnalysis
from .scripting import Scripting
from .mutagenesis import Mutagenesis
from .selection import Selection
from .distance_selector import DistanceSelector
from .consensus import Consensus
from .mutation_profiler import MutationProfiler
from .mutation_composer import MutationComposer
from .pool import Pool
from .pdb import PDB
from .sasa import SASA
from .pose_change import PoseChange
from .bayesian_adjuster import BayesianAdjuster
from .thermompnn import ThermoMPNN
from .vespag import VespaG
from .aggrescan3d import Aggrescan3D
from .plm_sol import PLM_Sol
from .converters import *
from .entities import *

__all__ = [
    # Core pipeline components
    'Pipeline',
    'Folder',
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
    'IndexedTableContainer',
    'get_table',
    'get_table_path',
    'get_indexed_table',
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
    'RFdiffusion2',
    'RFdiffusion3',
    'HBDesigner',
    'ProteinMPNN',
    'Frame2Seq',
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
    'DNAEncoder',
    'RBSDesigner',
    'Load',
    'LoadMultiple',
    'SequenceMetricCorrelation',
    'PyMOL',
    'Fuse',
    'Plot',
    'Panda',
    'MSA',
    'MMseqs2',
    'MMseqs2Server',
    'ReMap',
    'Mock',
    'Gnina',
    'DiffDock',
    'PocketGen',
    'PLACER',
    'DynamicBind',
    'NeuralPLexer',
    'PoseBusters',
    'PoseBustersTool',
    'ADMETAI',
    'ESMFold',
    'UniProt',
    'RDKit',
    'OpenBabel',
    'PLIP',
    'Prodigy',
    'ProLIF',
    'FPocket',
    'P2Rank',
    'AF2BIND',
    'DSSP',
    'APBS',
    'OpenMM',
    'Reduce',
    'XTB',
    'RTMScore',
    'GEMS',
    'BioEmu',
    'CABSflex',
    'EnsembleAnalysis',
    'Scripting',
    'Mutagenesis',
    'Selection',
    'DistanceSelector',
    'Consensus',
    'MutationProfiler',
    'MutationComposer',
    'Pool',
    'PDB',
    'SASA',
    'PoseChange',
    'BayesianAdjuster',
    'ThermoMPNN',
    'VespaG',
    'Aggrescan3D',
    'PLM_Sol',

    # Utility functions
    'pdb_to_jsonl',
    'fasta_to_csv',
    'validate_pdb',
    'validate_fasta',
    'combine_fasta_files'
]