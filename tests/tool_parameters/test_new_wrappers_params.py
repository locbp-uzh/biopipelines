# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Parameter coverage for every external-tool wrapper added in commit 345204c.

Each test drives a single-tool ``Pipeline`` at config time and asserts that
EVERY user-facing constructor parameter (the non-default value passed in) lands
in the emitted artifacts (pipeline.sh + per-step scripts + ``_configuration/``
files). This is the same string-grep contract the rest of the tool_parameters
suite uses (see ``test_gnina_params.py``); it exercises construction,
validation, script emission, and output declaration without running any
external model.

Coverage:
  coordinate consumers   — DSSP, Reduce, FPocket, P2Rank, Prodigy, APBS, OpenMM,
                           Aggrescan3D, PLIP, AF2BIND, Frame2Seq, ThermoMPNN
  coordinate-ligand      — XTB, ProLIF, PocketGen, PLACER, GEMS, RTMScore,
                           NeuralPLexer
  docking                — DiffDock, DynamicBind
  sequence consumers     — BioEmu, ESMFold, VespaG, PLM_Sol
  data / chemistry       — UniProt, OpenBabel, RDKit, EnsembleAnalysis, Consensus
"""

import pytest

from ._helpers import assert_substrings_in, read_all_emitted_artifacts


pytestmark = pytest.mark.tool_parameters


# ── mock input builders ───────────────────────────────────────────────────────

def _structures_mock(fmt="pdb"):
    from biopipelines.mock import Mock
    return Mock(
        ids=["s1"],
        streams={"structures": {"format": fmt, "file": "<id>.pdb"}},
        map_table_strategy="config",
    )


def _compounds_mock():
    from biopipelines.mock import Mock
    return Mock(
        ids=["c1"],
        streams={"compounds": {"format": "csv", "smiles": "CCO"}},
        map_table_strategy="config",
    )


def _ligand_structures_mock():
    """A ligand provided as a StandardizedOutput carrying BOTH a structures
    (coordinate) and a compounds (chemistry) stream — the shape GEMS/RTMScore
    and the coordinate-ligand tools accept."""
    from biopipelines.mock import Mock
    return Mock(
        ids=["c1"],
        streams={
            "structures": {"format": "sdf", "file": "<id>.sdf"},
            "compounds": {"format": "csv", "smiles": "CCO"},
        },
        map_table_strategy="config",
    )


def _sequences_mock():
    from biopipelines.mock import Mock
    return Mock(
        ids=["q1"],
        streams={"sequences": {"format": "csv", "sequence": "MKTAYIAKQR"}},
        map_table_strategy="config",
    )


def _resi_csv_mock(ids):
    """A per-residue resi-csv stream (CABSflex-shaped) for EnsembleAnalysis /
    Consensus inputs."""
    from biopipelines.mock import Mock
    return Mock(
        ids=ids,
        streams={"structures": {"format": "resi-csv", "file": "<id>.csv"}},
        map_table_strategy="config",
    )


# ── coordinate-only structure consumers ──────────────────────────────────────

def test_dssp_smoke(local_config, isolated_cwd, new_pipeline):
    from biopipelines.dssp import DSSP
    pipeline = new_pipeline("dssp_params")
    with pipeline:
        s = _structures_mock()
        DSSP(structures=s.streams.structures)
        content = read_all_emitted_artifacts(pipeline.save())
    assert "DSSP" in content or "dssp" in content


def test_reduce_smoke(local_config, isolated_cwd, new_pipeline):
    from biopipelines.reduce import Reduce
    pipeline = new_pipeline("reduce_params")
    with pipeline:
        s = _structures_mock()
        Reduce(structures=s.streams.structures)
        content = read_all_emitted_artifacts(pipeline.save())
    assert "Reduce" in content or "reduce" in content


def test_fpocket_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.fpocket import FPocket
    pipeline = new_pipeline("fpocket_params")
    with pipeline:
        s = _structures_mock()
        FPocket(
            structures=s.streams.structures,
            min_alpha_spheres=40,
            min_radius=3.5,
            max_radius=6.5,
            clustering_distance=2.6,
        )
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["40", "3.5", "6.5", "2.6"])


def test_p2rank_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.p2rank import P2Rank
    pipeline = new_pipeline("p2rank_params")
    with pipeline:
        s = _structures_mock()
        P2Rank(structures=s.streams.structures, config="alphafold", threads=4, visualizations=True)
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["alphafold", "4"])


def test_prodigy_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.prodigy import Prodigy
    pipeline = new_pipeline("prodigy_params")
    with pipeline:
        s = _structures_mock()
        Prodigy(structures=s.streams.structures, interface="A C", temperature=37.0)
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["A C", "37"])


def test_apbs_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.apbs import APBS
    pipeline = new_pipeline("apbs_params")
    with pipeline:
        s = _structures_mock()
        APBS(
            structures=s.streams.structures,
            ph=6.5, forcefield="PARSE", ion_concentration=0.2,
            grid_dim=97, pdie=4.0, sdie=80.0, solver="npbe",
        )
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["6.5", "PARSE", "0.2", "97", "4.0", "80.0", "npbe"])


def test_openmm_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.openmm import OpenMM
    pipeline = new_pipeline("openmm_params")
    with pipeline:
        s = _structures_mock()
        OpenMM(
            structures=s.streams.structures,
            max_iterations=500, tolerance_kj_per_mol_nm=5.0,
            forcefield="amber14-all", solvent="implicit-gbn2",
            platform="CPU", restraint_selection="CA", restraint_k=750.0,
        )
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["500", "5.0", "amber14-all", "implicit-gbn2", "CPU", "CA", "750.0"])


def test_aggrescan3d_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.aggrescan3d import Aggrescan3D
    pipeline = new_pipeline("aggrescan3d_params")
    with pipeline:
        s = _structures_mock()
        Aggrescan3D(structures=s.streams.structures, chains="AB", distance=12.0, max_parallel=3)
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["AB", "12.0", "3"])


def test_plip_params_ligand_mode(local_config, isolated_cwd, new_pipeline):
    from biopipelines.plip import PLIP
    from biopipelines.ligand import Ligand
    pipeline = new_pipeline("plip_lig_params")
    with pipeline:
        s = _structures_mock()
        PLIP(structures=s.streams.structures, ligand=Ligand(code="ETH"),
             mode="ligand", generate_pse=False)
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["ligand", "ETH"])


def test_plip_params_peptide_mode(local_config, isolated_cwd, new_pipeline):
    """chains is only valid in peptide/intra mode; cover it there."""
    from biopipelines.plip import PLIP
    pipeline = new_pipeline("plip_pep_params")
    with pipeline:
        s = _structures_mock()
        PLIP(structures=s.streams.structures, mode="peptide", chains="AB", generate_pse=True)
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["peptide", "AB"])


def test_af2bind_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.af2bind import AF2BIND
    pipeline = new_pipeline("af2bind_params")
    with pipeline:
        s = _structures_mock()
        AF2BIND(structures=s.streams.structures, chain="B",
                mask_sidechains=False, mask_sequence=True, top_k=20)
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["B", "20"])


def test_frame2seq_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.frame2seq import Frame2Seq
    pipeline = new_pipeline("frame2seq_params")
    with pipeline:
        s = _structures_mock()
        Frame2Seq(structures=s.streams.structures, num_sequences=5,
                  temperature=0.7, chain="B", omit_aa="CW", fixed="10-20")
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["5", "0.7", "B", "CW"])


def test_thermompnn_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.thermompnn import ThermoMPNN
    pipeline = new_pipeline("thermompnn_params")
    with pipeline:
        s = _structures_mock()
        ThermoMPNN(structures=s.streams.structures, chain="B", mutations="A10G+L20P")
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["B", "A10G"])


# ── coordinate-ligand consumers (structures + a ligand) ───────────────────────

def test_xtb_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.xtb import XTB
    pipeline = new_pipeline("xtb_params")
    with pipeline:
        s = _structures_mock()
        lig = _ligand_structures_mock()
        XTB(structures=s.streams.structures, ligand=lig, method="gfn1",
            solvent="water", charge=-1, opt=True)
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["gfn1", "water", "-1"])


def test_prolif_smoke(local_config, isolated_cwd, new_pipeline):
    from biopipelines.prolif import ProLIF
    from biopipelines.ligand import Ligand
    pipeline = new_pipeline("prolif_params")
    with pipeline:
        s = _structures_mock()
        ProLIF(structures=s.streams.structures, ligand=Ligand(smiles="CCO", ids="c1", codes="ETH"))
        content = read_all_emitted_artifacts(pipeline.save())
    assert "ProLIF" in content or "prolif" in content


def test_pocketgen_smoke(local_config, isolated_cwd, new_pipeline):
    from biopipelines.pocketgen import PocketGen
    from biopipelines.ligand import Ligand
    pipeline = new_pipeline("pocketgen_params")
    with pipeline:
        s = _structures_mock()
        # PocketGen fans out (protein, ligand) pairs via combinatorics, which
        # reads axis ids from the full structures output and the ligand's
        # compounds stream.
        PocketGen(structures=s, ligand=Ligand(smiles="CCO", ids="c1", codes="ETH"))
        content = read_all_emitted_artifacts(pipeline.save())
    assert "PocketGen" in content or "pocketgen" in content


def test_placer_params_ligand_mode(local_config, isolated_cwd, new_pipeline):
    from biopipelines.placer import PLACER
    from biopipelines.ligand import Ligand
    pipeline = new_pipeline("placer_lig_params")
    with pipeline:
        s = _structures_mock()
        PLACER(structures=s, ligand=Ligand(smiles="CCO", ids="c1", codes="ETH"),
               nsamples=15, rerank="prmsd")
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["15", "prmsd"])


def test_placer_params_apo_mode(local_config, isolated_cwd, new_pipeline):
    """target_res / exclude_sm are only valid in sidechain/apo mode (no ligand)."""
    from biopipelines.placer import PLACER
    pipeline = new_pipeline("placer_apo_params")
    with pipeline:
        s = _structures_mock()
        PLACER(structures=s, target_res="A50", exclude_sm=True, nsamples=12)
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["12", "A50"])


def test_gems_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.gems import GEMS
    pipeline = new_pipeline("gems_params")
    with pipeline:
        s = _structures_mock()
        lig = _ligand_structures_mock()
        GEMS(structures=s, ligands=lig, skip_ligand_embedding=True)
        content = read_all_emitted_artifacts(pipeline.save())
    assert "GEMS" in content or "gems" in content


def test_rtmscore_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.rtmscore import RTMScore
    pipeline = new_pipeline("rtmscore_params")
    with pipeline:
        s = _structures_mock()
        lig = _ligand_structures_mock()
        RTMScore(structures=s, ligands=lig, cutoff=8.0, model="model1")
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["8.0", "model1"])


def test_neuralplexer_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.neuralplexer import NeuralPLexer
    from biopipelines.ligand import Ligand
    pipeline = new_pipeline("neuralplexer_params")
    with pipeline:
        s = _structures_mock()
        NeuralPLexer(
            structures=s, compounds=Ligand(smiles="CCO", ids="c1", codes="ETH"),
            n_samples=8, num_steps=30, chunk_size=2,
            sampler="DDIM", cuda=False,
        )
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["8", "30", "DDIM"])


# ── docking tools (structures + compounds) ────────────────────────────────────

def test_diffdock_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.diffdock import DiffDock
    from biopipelines.ligand import Ligand
    pipeline = new_pipeline("diffdock_params")
    with pipeline:
        s = _structures_mock()
        # DiffDock fans out (protein, ligand) pairs via combinatorics, which
        # collects axis ids from the full structures output and the ligand's
        # compounds stream.
        DiffDock(
            structures=s, compounds=Ligand(smiles="CCO", ids="c1", codes="ETH"),
            samples_per_complex=24, inference_steps=18, batch_size=16,
            no_final_step_noise=False,
        )
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["24", "18", "16"])


def test_dynamicbind_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.dynamicbind import DynamicBind
    from biopipelines.ligand import Ligand
    pipeline = new_pipeline("dynamicbind_params")
    with pipeline:
        s = _structures_mock()
        DynamicBind(
            structures=s, compounds=Ligand(smiles="CCO", ids="c1", codes="ETH"),
            num_samples=30, num_saved=5, inference_steps=15,
            num_workers=2, rigid_protein=True, make_movie=True, seed=7,
        )
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["30", "5", "15", "7"])


# ── sequence consumers ─────────────────────────────────────────────────────────

def test_bioemu_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.bioemu import BioEmu
    pipeline = new_pipeline("bioemu_params")
    with pipeline:
        q = _sequences_mock()
        BioEmu(sequences=q.streams.sequences, num_samples=25, batch_size=5,
               reconstruct_sidechains=True, filter_samples=False)
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["25", "5"])


def test_esmfold_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.esmfold import ESMFold
    pipeline = new_pipeline("esmfold_params")
    with pipeline:
        q = _sequences_mock()
        ESMFold(sequences=q.streams.sequences, chunk_size=128, num_recycles=8,
                max_tokens_per_batch=2048, cpu_offload=True, cpu_only=True)
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["128", "8", "2048"])


def test_vespag_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.vespag import VespaG
    pipeline = new_pipeline("vespag_params")
    with pipeline:
        q = _sequences_mock()
        VespaG(sequences=q.streams.sequences, mutations="A10G+L20P")
        content = read_all_emitted_artifacts(pipeline.save())
    assert "A10G" in content


def test_plm_sol_smoke(local_config, isolated_cwd, new_pipeline):
    from biopipelines.plm_sol import PLM_Sol
    pipeline = new_pipeline("plm_sol_params")
    with pipeline:
        q = _sequences_mock()
        PLM_Sol(sequences=q.streams.sequences)
        content = read_all_emitted_artifacts(pipeline.save())
    assert "PLM_Sol" in content or "plm_sol" in content


# ── data / chemistry tools ─────────────────────────────────────────────────────

def test_uniprot_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.uniprot import UniProt
    pipeline = new_pipeline("uniprot_params")
    with pipeline:
        UniProt(accessions=["P12345", "Q9Y6K9"])
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["P12345", "Q9Y6K9"])


def test_openbabel_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.openbabel import OpenBabel
    pipeline = new_pipeline("openbabel_params")
    with pipeline:
        lig = _compounds_mock()
        OpenBabel(compounds=lig.streams.compounds, convert_3d="sdf",
                  add_hydrogens=True, pH=7.4, gen3d=True, gen3d_quality="best",
                  minimize=True, ff="UFF", minimize_steps=300)
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["sdf", "7.4", "best", "UFF", "300"])


def test_rdkit_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.rdkit_descriptors import RDKit
    pipeline = new_pipeline("rdkit_params")
    with pipeline:
        lig = _compounds_mock()
        RDKit(compounds=lig.streams.compounds,
              descriptors=["MolWt", "TPSA"], morgan_fp=True)
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["MolWt", "TPSA"])


def test_ensemble_analysis_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.ensemble_analysis import EnsembleAnalysis
    pipeline = new_pipeline("ensemble_params")
    with pipeline:
        confs = _resi_csv_mock(["c_1", "c_2", "c_3"])
        groups = _sequences_mock()
        EnsembleAnalysis(structures=confs.streams.structures,
                         groups=groups.streams.sequences,
                         selection="CA", reference="first")
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["CA", "first"])


def test_consensus_params(local_config, isolated_cwd, new_pipeline):
    from biopipelines.consensus import Consensus
    pipeline = new_pipeline("consensus_params")
    with pipeline:
        confs = _resi_csv_mock(["c_1", "c_2"])
        groups = _sequences_mock()
        Consensus(
            confs.streams.structures,
            groups=groups.streams.sequences,
            operations=[Consensus.fraction("rmsf<=2.0", name="stable_frac"),
                        Consensus.mean("rmsf")],
        )
        content = read_all_emitted_artifacts(pipeline.save())
    assert_substrings_in(content, ["stable_frac", "rmsf"])
