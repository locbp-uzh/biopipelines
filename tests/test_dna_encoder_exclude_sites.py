"""Codon choice must keep restriction sites out of the encoded DNA.

A site landing inside a synthesised gene ruins the cloning strategy it was chosen
for, and the site can straddle codon boundaries — so filtering one codon at a
time is not enough.
"""

import importlib.util
import os
import random
import sys

import pytest

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HELPER = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_dna_encoder.py")

pytest.importorskip("pandas")
pytest.importorskip("Bio")


def _load():
    spec = importlib.util.spec_from_file_location("pipe_dna_encoder", HELPER)
    m = importlib.util.module_from_spec(spec)
    sys.modules["pipe_dna_encoder"] = m
    spec.loader.exec_module(m)
    return m


@pytest.fixture(scope="module")
def mod():
    return _load()


@pytest.fixture(scope="module")
def ec_table(mod):
    freqs = mod.parse_cocoputs_table(mod.COCOPUTS_TABLES["Escherichia coli"])
    return mod.create_aa_codon_frequency_table(freqs)


def _has(dna, site):
    from itertools import product
    return site in dna


def test_reverse_complement_is_excluded_too(mod):
    """A site on the other strand cuts just as well."""
    expanded = mod.expand_sites(["GGTCTC"])          # BsaI, not palindromic
    assert "GGTCTC" in expanded
    assert "GAGACC" in expanded, "reverse complement must also be forbidden"


def test_palindrome_expands_to_itself(mod):
    assert mod.expand_sites(["GAATTC"]) == ["GAATTC"]


def test_common_sites_absent_from_encoded_dna(mod, ec_table):
    """The encoder must not emit any excluded site, on either strand."""
    random.seed(0)
    protein = ("LHHPVFQQESFTRQVLWKLLKVVKFGEVISYSHLAALAGNPAATAAVKTALSGNPVPILIPCHRVV"
               "QGDLDVGGYEGGLAVKEWLLAHEGHRLGKPATKAEIDAEMKTASAEAKRFMDQVKAYLDDP")
    sites = ["GAATTC", "GGATCC", "AAGCTT", "CATATG", "CTCGAG", "GCGGCCGC"]
    forbidden = mod.expand_sites(sites)
    for _ in range(20):
        dna = mod.encode_sequence_excluding(protein, ec_table, forbidden)
        assert len(dna) == 3 * len(protein)
        for s in forbidden:
            assert s not in dna, f"encoder emitted {s}"


def test_translation_is_preserved(mod, ec_table):
    """Excluding sites must not change the protein."""
    from Bio.Seq import Seq
    random.seed(1)
    protein = "MEEFRRKLAAGGSSWWYYCCDDEEHHIIKKLLMMNNPPQQRRSSTTVVWWYY"
    forbidden = mod.expand_sites(["GAATTC", "GGATCC", "GGTCTC", "CGTCTC"])
    dna = mod.encode_sequence_excluding(protein, ec_table, forbidden)
    assert str(Seq(dna).translate()) == protein


def test_site_spanning_a_codon_boundary_is_caught(mod, ec_table):
    """The dangerous case: no single codon contains the site.

    EcoRI's GAATTC across E-F is GAA|TTC — two perfectly ordinary codons that
    only form the site together. A per-codon filter would miss it entirely.
    """
    random.seed(2)
    protein = "EFEFEFEFEFEFEFEF"
    forbidden = mod.expand_sites(["GAATTC"])
    for _ in range(30):
        dna = mod.encode_sequence_excluding(protein, ec_table, forbidden)
        assert "GAATTC" not in dna
        from Bio.Seq import Seq
        assert str(Seq(dna).translate()) == protein


def test_unavoidable_site_raises_rather_than_emitting_it(mod, ec_table):
    """Methionine and tryptophan have one codon each: ATG and TGG.

    MW can only ever be ATGTGG, so excluding it is impossible — the encoder must
    say so instead of returning DNA that contains the site.
    """
    forbidden = mod.expand_sites(["ATGTGG"])
    with pytest.raises(ValueError, match="without the excluded sites"):
        mod.encode_sequence_excluding("MWMWMW", ec_table, forbidden, max_restarts=5)


def test_no_sites_matches_unconstrained_behaviour(mod, ec_table):
    """With nothing to exclude the output is still a valid encoding."""
    from Bio.Seq import Seq
    random.seed(3)
    protein = "MKVLATTGEER"
    dna = mod.encode_sequence_excluding(protein, ec_table, [])
    assert str(Seq(dna).translate()) == protein


def test_iupac_sites_are_honoured(mod, ec_table):
    """Degenerate recognition sequences (e.g. BstYI RGATCY) must also be avoided."""
    random.seed(4)
    protein = "MDSRIDSRIDSRIDSRI"
    forbidden = mod.expand_sites(["RGATCY"])
    dna = mod.encode_sequence_excluding(protein, ec_table, forbidden)
    for a in "AG":
        for b in "CT":
            assert f"{a}GATC{b}" not in dna
