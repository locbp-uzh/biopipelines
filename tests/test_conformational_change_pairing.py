"""ConformationalChange must be able to measure every atom it was asked about.

PyMOL's alignment defaults answer "how similar are these structures where they
match?", which is the right question for homologues and the wrong one for "did
this sequence fold the way it was designed". Three defaults conspire:

* `align` pairs by SEQUENCE, so residues it cannot match are excluded. After two
  rounds of MPNN the designed segment is ~41% identical to the design it came
  from, by construction.
* `cycles=5` then discards the worst-fitting survivors, and its `cutoff=2.0`
  coincides exactly with the < 2 A gate this project applies.
* `cealign` reports the best common fragment path, capping the measurement at
  16-32 atoms regardless of how long the segment is.

Measured on design 60_Panda_1_104_1, whose designed segment is 11.43 A from its
design over all 200 backbone atoms:

    align   (defaults)  ->  1.89 A over 105 atoms   <- would pass a < 2 A gate
    cealign (defaults)  ->  3.75 A over  32 atoms
    fit matchmaker=-1   -> 11.43 A over 200 atoms
    and in the frame of the superposed core, 18.89 A

So the tool could report 1.89 A for a segment sitting 19 A from where it was
designed. These tests pin the parameters that make the honest measurement
reachable.

The runtime tests drive ``pipe_conformational_change.align_and_compute_rmsd`` against a stand-in PyMOL command object, so they assert on what the code does with cycles/cutoff/matchmaker rather than on the text of its source.
"""

import importlib.util
import inspect
import math
import os
import sys

import pytest


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNTIME_SCRIPT = os.path.join(REPO_ROOT, "pipe_scripts", "pipe_conformational_change.py")


def _sig():
    from biopipelines.conformational_change import ConformationalChange
    return inspect.signature(ConformationalChange.__init__).parameters


# ── config-time behaviour ─────────────────────────────────────────────────────

def _pair_of_structure_mocks():
    from biopipelines.mock import Mock
    ref = Mock(ids=["r1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
               map_table_strategy="config")
    tgt = Mock(ids=["t1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
               map_table_strategy="config")
    return ref, tgt


@pytest.fixture
def make_change(local_config, isolated_cwd, new_pipeline):
    """Construct a ConformationalChange inside a live pipeline.

    Parameter validation runs during registration, so an invalid kwarg raises here -- which is what lets these tests assert on behaviour instead of grepping the source of ``validate_params``.
    """
    from biopipelines.conformational_change import ConformationalChange

    pipeline = new_pipeline("cc_pairing")
    with pipeline:
        ref, tgt = _pair_of_structure_mocks()

        def _make(**kwargs):
            return ConformationalChange(
                reference_structures=ref.streams.structures,
                target_structures=tgt.streams.structures,
                **kwargs,
            )
        yield _make


def test_pairing_cycles_cutoff_and_frame_are_exposed(make_change):
    """The four knobs must be real, accepted constructor parameters."""
    p = _sig()
    for name in ("pairing", "cycles", "cutoff", "frame"):
        assert name in p, f"ConformationalChange must expose {name!r}"
    # And they must actually be accepted together, not merely present in the signature.
    make_change(pairing="ordered", cycles=0, cutoff=1.5, frame="A1-10", selection="A20-30")


def test_defaults_preserve_existing_behaviour():
    """Existing pipelines must not silently change their numbers."""
    p = _sig()
    assert p["pairing"].default == "sequence"
    assert p["cycles"].default == 5      # PyMOL's own default
    assert p["cutoff"].default == 2.0
    assert p["frame"].default is None


@pytest.mark.parametrize("pairing", ["sequence", "ordered", "identifier"])
def test_valid_pairings_accepted(pairing, make_change):
    """Each documented pairing mode must construct without raising."""
    make_change(pairing=pairing)


def test_invalid_pairing_is_rejected(make_change):
    with pytest.raises(ValueError, match="pairing must be"):
        make_change(pairing="bogus")


@pytest.mark.parametrize("cycles", [-1, True, 2.5])
def test_negative_cycles_rejected(cycles, make_change):
    """bool is an int subclass, so cycles=True must be rejected too."""
    with pytest.raises(ValueError, match="cycles must be a non-negative integer"):
        make_change(cycles=cycles)


def test_frame_requires_a_selection_to_measure(make_change):
    with pytest.raises(ValueError, match="frame requires selection"):
        make_change(frame="A1-10")


def test_cealign_rejects_a_non_default_cycles(make_change):
    """cealign has no outlier-rejection cycles, so asking for them is a mistake."""
    with pytest.raises(ValueError, match="cealign has no outlier-rejection cycles"):
        make_change(alignment="cealign", cycles=3)


def test_non_sequence_pairing_rejects_an_alignment_method(make_change):
    """ordered/identifier go through cmd.fit, where `alignment` means nothing."""
    with pytest.raises(ValueError, match="only applies to pairing='sequence'"):
        make_change(pairing="ordered", alignment="cealign")


def test_table_reports_pre_refinement_numbers(make_change):
    """Without these columns, a 47% inflated pass rate is invisible in the output.

    RMSD alone cannot tell you that refinement threw away half the segment; the
    before/after pair and the dropped percentage are what make it checkable.
    """
    out = make_change()
    columns = set(vars(out.tables.changes))
    for col in ("RMSD", "RMSD_before", "num_aligned_atoms", "num_atoms_before",
                "num_residues_aligned", "atoms_dropped_pct"):
        assert col in columns, f"changes table must carry {col}, got {sorted(columns)}"


# ── runtime behaviour, driven against a stand-in PyMOL ────────────────────────

class _FakeCmd:
    """Stand-in for ``pymol.cmd``, recording how it was called.

    ``counts``/``coords`` keys are matched as substrings of the selection string PyMOL would receive, longest key first, so a test can give a different answer for a frame selection than for the measured one.
    """

    def __init__(self, counts=None, coords=None, align_result=None, fit_results=None):
        self.counts = counts or {}
        self.coords = coords or {}
        self.align_result = align_result
        self.fit_results = list(fit_results or [])
        self.align_calls = []
        self.super_calls = []
        self.cealign_calls = []
        self.fit_calls = []
        self.rms_cur_calls = []

    def _lookup(self, table, selection):
        for key in sorted(table, key=len, reverse=True):
            if key in selection:
                return table[key]
        raise AssertionError(f"no fake entry matches selection {selection!r}")

    def count_atoms(self, selection):
        return self._lookup(self.counts, selection)

    def align(self, mobile, target, cycles=5, cutoff=2.0):
        self.align_calls.append({"mobile": mobile, "target": target,
                                 "cycles": cycles, "cutoff": cutoff})
        return self.align_result

    def super(self, mobile, target, cycles=5, cutoff=2.0):
        self.super_calls.append({"mobile": mobile, "target": target,
                                 "cycles": cycles, "cutoff": cutoff})
        return self.align_result

    def cealign(self, ref, target):
        self.cealign_calls.append({"ref": ref, "target": target})
        return {"RMSD": 3.75, "alignment_length": 32}

    def fit(self, mobile, target, matchmaker=0, cycles=5, cutoff=2.0):
        self.fit_calls.append({"mobile": mobile, "target": target,
                               "matchmaker": matchmaker, "cycles": cycles,
                               "cutoff": cutoff})
        return self.fit_results[len(self.fit_calls) - 1]

    def rms_cur(self, mobile, target, matchmaker=-1):
        self.rms_cur_calls.append({"mobile": mobile, "target": target,
                                   "matchmaker": matchmaker})
        a = self._lookup(self.coords, mobile)
        b = self._lookup(self.coords, target)
        assert len(a) == len(b), "fake coords must be paired"
        total = sum((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
                    for (x1, y1, z1), (x2, y2, z2) in zip(a, b))
        return math.sqrt(total / len(a))


@pytest.fixture
def runtime():
    """Load pipe_conformational_change.py against a stubbed ``pymol``."""
    import types

    saved = {name: sys.modules.get(name) for name in ("pymol", "pymol.cmd", "pipe_cc_under_test")}
    stub = types.ModuleType("pymol")
    stub.cmd = types.SimpleNamespace()
    stub.finish_launching = lambda *a, **k: None
    sys.modules["pymol"] = stub
    sys.modules["pymol.cmd"] = stub.cmd

    spec = importlib.util.spec_from_file_location("pipe_cc_under_test", RUNTIME_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["pipe_cc_under_test"] = module
    spec.loader.exec_module(module)

    yield module

    for name, value in saved.items():
        if value is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = value


def test_runtime_passes_cycles_to_pymol(runtime, record_case):
    """The bug was calling cmd.align with no cycles argument at all."""
    fake = _FakeCmd(
        counts={"tgt": 200, "tgt and name CA": 26, "ref": 200},
        # (rmsd, n_after, _, rmsd_before, n_before, _, residues)
        align_result=(1.89, 105, 0, 11.43, 200, 0, 26),
    )
    runtime.cmd = fake
    runtime.align_and_compute_rmsd("ref", "tgt", None, "align", cycles=3, cutoff=1.5)

    record_case(
        input="align_and_compute_rmsd(..., cycles=3, cutoff=1.5)",
        expected=[{"cycles": 3, "cutoff": 1.5}],
        actual=[{k: c[k] for k in ("cycles", "cutoff")} for c in fake.align_calls],
    )
    assert len(fake.align_calls) == 1
    assert fake.align_calls[0]["cycles"] == 3, "align must receive the configured cycles"
    assert fake.align_calls[0]["cutoff"] == 1.5, "align must receive the configured cutoff"


def test_runtime_reports_pre_refinement_numbers(runtime, record_case):
    """The real 60_Panda_1_104_1 numbers: 1.89 A survives after 47% was dropped."""
    fake = _FakeCmd(
        counts={"tgt": 200, "tgt and name CA": 26, "ref": 200},
        align_result=(1.89, 105, 0, 11.43, 200, 0, 26),
    )
    runtime.cmd = fake
    result = runtime.align_and_compute_rmsd("ref", "tgt", None, "align")

    expected = {
        "RMSD": 1.89, "num_aligned_atoms": 105,
        "RMSD_before": 11.43, "num_atoms_before": 200,
        "num_residues_aligned": 26, "atoms_dropped_pct": 47.5,
    }
    record_case(input="align -> 1.89 A over 105 of 200 atoms",
                expected=expected, actual=result)
    assert result == expected


def test_runtime_warns_when_atoms_are_dropped(runtime, capsys):
    """Dropping a large fraction of atoms must be announced, not silent."""
    fake = _FakeCmd(
        counts={"tgt": 200, "tgt and name CA": 26, "ref": 200},
        align_result=(1.89, 105, 0, 11.43, 200, 0, 26),
    )
    runtime.cmd = fake
    runtime.align_and_compute_rmsd("ref", "tgt", None, "align")
    out = capsys.readouterr().out
    assert "refinement dropped 48% of atoms" in out, out  # 47.5 at .0f
    assert "200 -> 105" in out, out


def test_runtime_stays_quiet_when_nothing_is_dropped(runtime, capsys):
    fake = _FakeCmd(
        counts={"tgt": 200, "tgt and name CA": 26, "ref": 200},
        align_result=(1.89, 200, 0, 1.89, 200, 0, 26),
    )
    runtime.cmd = fake
    result = runtime.align_and_compute_rmsd("ref", "tgt", None, "align")
    assert result["atoms_dropped_pct"] == 0.0
    assert "refinement dropped" not in capsys.readouterr().out


@pytest.mark.parametrize("pairing,expected_matchmaker", [("ordered", -1), ("identifier", 0)])
def test_ordered_and_identifier_pairing_use_fit_matchmaker(
    runtime, pairing, expected_matchmaker,
):
    """ordered pairs the Nth atom with the Nth (-1); identifier pairs on ids (0)."""
    fake = _FakeCmd(
        counts={"tgt": 200, "tgt and name CA": 26, "ref": 200},
        fit_results=[11.43, 11.43],
    )
    runtime.cmd = fake
    result = runtime.align_and_compute_rmsd("ref", "tgt", None, "align", pairing=pairing)

    assert not fake.align_calls, "non-sequence pairing must not go through cmd.align"
    assert fake.fit_calls, "ordered/identifier pairing must use cmd.fit"
    assert fake.fit_calls[0]["matchmaker"] == expected_matchmaker
    assert result["RMSD"] == 11.43


def test_ordered_pairing_checks_atom_counts_match(runtime):
    """matchmaker=-1 pairs the Nth atom with the Nth; mismatched counts are silent nonsense."""
    fake = _FakeCmd(
        counts={"tgt": 200, "tgt and name CA": 26, "ref": 180},
        fit_results=[11.43, 11.43],
    )
    runtime.cmd = fake
    with pytest.raises(ValueError, match="same atom count on both sides"):
        runtime.align_and_compute_rmsd("ref", "tgt", None, "align", pairing="ordered")
    assert not fake.fit_calls, "cmd.fit must not run on mismatched selections"


def test_frame_measures_without_refitting(runtime, record_case):
    """A segment given its own superposition can fit itself while sitting elsewhere.

    60_Panda_1_104_1's segment scores 11.43 A on its own best fit and 18.89 A once
    the two cores are aligned — the second is the number that matters, because the
    docking box is placed via that same core superposition.

    Here the measured selection sits at displacements of 0, 4, 3 and 0 A from the reference, so the honest RMSD is sqrt((0+16+9+0)/4) = 2.5 A. The align call superposes the FRAME; the measurement must then be rms_cur on the selection, with no refit of its own.
    """
    measured_ref = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
    measured_tgt = [(0.0, 0.0, 0.0), (5.0, 0.0, 0.0), (2.0, 3.0, 0.0), (3.0, 0.0, 0.0)]
    fake = _FakeCmd(
        counts={
            "resi 1-10": 120, "resi 1-10) and name CA": 15,
            "resi 20-30": 4, "resi 20-30) and name CA": 1,
        },
        coords={"tgt and (chain A and resi 20-30)": measured_tgt,
                "ref and (chain A and resi 20-30)": measured_ref},
        align_result=(0.42, 120, 0, 0.42, 120, 0, 15),
    )
    runtime.cmd = fake
    result = runtime.align_and_compute_rmsd(
        "ref", "tgt", "A20-30", "align", frame="A1-10",
    )

    expected_rmsd = math.sqrt((0 + 16 + 9 + 0) / 4)
    record_case(input="frame=A1-10, selection=A20-30, displacements 0/4/3/0 A",
                expected=round(expected_rmsd, 6), actual=round(result["RMSD"], 6))

    assert fake.align_calls, "the frame itself must still be superposed"
    assert "resi 1-10" in fake.align_calls[0]["mobile"], \
        "align must superpose the frame, not the measured selection"
    assert len(fake.rms_cur_calls) == 1, \
        "with a frame, the measurement must use rms_cur (no refit)"
    assert "resi 20-30" in fake.rms_cur_calls[0]["mobile"]
    assert result["RMSD"] == pytest.approx(expected_rmsd) == pytest.approx(2.5)
    # No refit means the pre/post numbers are the same measurement.
    assert result["RMSD_before"] == pytest.approx(expected_rmsd)
    assert result["num_aligned_atoms"] == result["num_atoms_before"] == 4


def test_frame_rejects_a_measured_selection_with_mismatched_counts(runtime):
    """rms_cur answers 0.000 on unequal counts rather than raising."""
    fake = _FakeCmd(
        counts={
            "resi 1-10": 120, "resi 1-10) and name CA": 15,
            "tgt and (chain A and resi 20-30)": 4,
            "ref and (chain A and resi 20-30)": 7,
        },
        align_result=(0.42, 120, 0, 0.42, 120, 0, 15),
    )
    runtime.cmd = fake
    with pytest.raises(ValueError, match="pairs them by position"):
        runtime.align_and_compute_rmsd("ref", "tgt", "A20-30", "align", frame="A1-10")
    assert not fake.rms_cur_calls, "rms_cur must not run on mismatched selections"


def test_cealign_caps_the_measurement(runtime):
    """cealign reports its best common fragment path, not the whole selection."""
    fake = _FakeCmd(counts={"tgt": 200, "tgt and name CA": 26, "ref": 200})
    runtime.cmd = fake
    result = runtime.align_and_compute_rmsd("ref", "tgt", None, "cealign")
    assert fake.cealign_calls and not fake.align_calls
    assert result["RMSD"] == 3.75
    assert result["num_aligned_atoms"] == 32, \
        "cealign's alignment_length is the real atom count, not the selection size"


def test_unknown_alignment_and_pairing_are_rejected(runtime):
    fake = _FakeCmd(counts={"tgt": 200, "tgt and name CA": 26, "ref": 200})
    runtime.cmd = fake
    with pytest.raises(ValueError, match="Unknown alignment method"):
        runtime.align_and_compute_rmsd("ref", "tgt", None, "nonsense")
    with pytest.raises(ValueError, match="Unknown pairing"):
        runtime.align_and_compute_rmsd("ref", "tgt", None, "align", pairing="nonsense")


# ── still a source grep, deliberately ─────────────────────────────────────────

def test_frame_column_ids_are_mapped_not_looked_up_raw():
    """A frame column lives on a design-level table; targets are per-sequence folds.

    cuts.csv is keyed by design id (60_Panda_1_100) while the folds are
    60_Panda_1_100_1 — a raw dict lookup misses every one. The script already
    resolves the `selection` column through get_mapped_ids; `frame` must go
    through the same path rather than a string-prefix fallback.

    This one stays a source grep: the id resolution happens inside ``analyze_all_conformational_changes``, which needs a full runtime config (datastream JSONs, per-id map tables, a written selection CSV) and walks every structure pair, so driving it behaviourally would test the harness more than the mapping.
    """
    src = open(RUNTIME_SCRIPT, encoding="utf-8").read()
    assert "target_to_frame_id = get_mapped_ids(" in src, \
        "frame column ids must be resolved with get_mapped_ids"
    assert "no frame selection for" in src, \
        "an unresolvable frame must be reported, not silently treated as no frame"


def test_gnina_keeps_every_pose_by_default():
    """GNINA has no CNN score filter; the wrapper must not invent one silently.

    The old default of 0.5 discarded 463 of 578 docked designs. It selects on an
    axis uncorrelated with binding (Spearman +0.025 against affinity here) and its
    CNN is trained on natural complexes, so de novo pockets with a synthetic dye
    score low across the board — max 0.766, median 0.360. Every one of the eight
    best binders (-11.0 to -10.1 kcal/mol) fell below the cut.
    """
    from biopipelines.gnina import Gnina
    p = inspect.signature(Gnina.__init__).parameters
    assert p["cnn_score_threshold"].default == 0.0, \
        "a non-zero default silently drops poses GNINA itself would return"
