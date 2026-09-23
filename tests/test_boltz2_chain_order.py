# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Boltz2 chain lettering depends on whether an axis iterates.

generate_configs() letters chains in the order entries are appended, and it
appends ITERATION axes before static-only ones. So the same set of entities
gets different chain IDs depending only on whether one of them is wrapped in
Each() -- which matters to every downstream tool that addresses residues by
chain (LigandMPNN's `redesigned`, PLIP/ProLIF selections, RMSD chain args).

These tests pin the behaviour rather than assert it is desirable. Changing it
would silently re-letter chains in existing pipelines, so it is documented
here instead: if you feed Boltz2 output to a chain-addressed tool, make the
chain explicit and do not assume the first protein is A.
"""

import sys
from pathlib import Path
import types

import pytest

PIPE_SCRIPTS = Path(__file__).resolve().parent.parent / "pipe_scripts"
if str(PIPE_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(PIPE_SCRIPTS))


BINDER = "MEEKVEEIKK"
TARGET = "QVCTGTDMKL"
SMILES = "O=C(C)Oc1ccccc1C(=O)O"


@pytest.fixture
def args():
    return types.SimpleNamespace(
        affinity=True, single_sequence=False, template=None,
        template_chains=None, glycosylation=None, covalent_linkage=None,
        contacts=None, disulfide_bonds=None, metal_coord=None,
        pocket_residues=None, pocket_binder=None, pocket_max_distance=None,
    )


def _axes(protein_iterated):
    """Ligand axis declared first, as the Boltz2 wrapper emits it.

    protein_iterated=False -> Bundle(binder, target): nothing iterates
    protein_iterated=True  -> Bundle(Each(binder), target)
    """
    proteins = {"entity_type": "protein", "mode": "bundle", "static_first": False}
    if protein_iterated:
        proteins["iterated"] = [{"id": "binder", "sequence": BINDER}]
        proteins["static"] = [{"id": "target", "sequence": TARGET}]
    else:
        proteins["iterated"] = []
        proteins["static"] = [{"id": "binder", "sequence": BINDER},
                              {"id": "target", "sequence": TARGET}]
    return {
        "ligands": {"entity_type": "ligand", "mode": "bundle", "iterated": [],
                    "static": [{"id": "ligand", "smiles": SMILES}],
                    "static_first": False},
        "proteins": proteins,
    }


def _by_chain(config):
    """{chain_id: entity_kind} for one generated config."""
    return {entry[next(iter(entry))]["id"]: next(iter(entry))
            for entry in config["sequences"]}


def test_static_only_bundle_puts_ligand_first(args, record_case):
    """Bundle(binder, target) + ligand -> ligand takes chain A."""
    from pipe_boltz_config_unified import generate_configs

    configs = generate_configs(_axes(protein_iterated=False),
                               {"by_id": {}, "by_seq": {}}, args)
    (_cid, config), = configs
    by_chain = _by_chain(config)

    record_case(input="Bundle(binder, target), ligand",
                expected="A=ligand", actual=f"A={by_chain['A']}")
    assert by_chain == {"A": "ligand", "B": "protein", "C": "protein"}


def test_iterated_protein_axis_takes_chain_a(args, record_case):
    """Bundle(Each(binder), target) + ligand -> binder takes chain A.

    Same entities as the test above; only Each() differs.
    """
    from pipe_boltz_config_unified import generate_configs

    configs = generate_configs(_axes(protein_iterated=True),
                               {"by_id": {}, "by_seq": {}}, args)
    (_cid, config), = configs
    by_chain = _by_chain(config)

    record_case(input="Bundle(Each(binder), target), ligand",
                expected="A=protein", actual=f"A={by_chain['A']}")
    assert by_chain == {"A": "protein", "B": "protein", "C": "ligand"}


def test_affinity_binder_follows_the_ligand_chain(args):
    """The affinity property must name whichever chain the ligand landed on."""
    from pipe_boltz_config_unified import generate_configs

    for iterated, expected in ((False, "A"), (True, "C")):
        configs = generate_configs(_axes(protein_iterated=iterated),
                                   {"by_id": {}, "by_seq": {}}, args)
        (_cid, config), = configs
        assert config["properties"] == [{"affinity": {"binder": expected}}]
