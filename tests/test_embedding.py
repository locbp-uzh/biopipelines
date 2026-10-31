"""Embedding: a registry of embedding models, one vector per input, and cosine similarities between chosen pairs.

The numerical core and the step loop run against a stand-in embedder, so the registry contract, the id handling, the `.npy` shapes and the pair arithmetic are pinned here without torch. The adapter block builds the Progres adapter over a torch-only stand-in of the `progres` package. The last block runs the real `progres` package on real weights and is skipped unless EMBEDDING_TEST_WEIGHTS names a checkpoint (the Progres and CIRPIN checkpoints of the CIRPIN repository both work) and `progres`, `torch_scatter` and `torch_geometric` import.
"""

import ast
import importlib.util
import math
import os
import re
import sys
import types

import numpy as np
import pytest

pytest.importorskip("yaml")

ROOT = os.path.join(os.path.dirname(__file__), "..")
SCRIPT = os.path.join(ROOT, "pipe_scripts", "pipe_embedding.py")
WRAPPER = os.path.join(ROOT, "biopipelines", "embedding.py")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def pipe():
    return load("pipe_embedding", SCRIPT)


# ── the registry ─────────────────────────────────────────────────────────────

def test_the_registry_holds_the_progres_family():
    from biopipelines.embedding import EMBEDDING_MODELS

    assert tuple(EMBEDDING_MODELS) == ("progres", "cirpin")
    for name, entry in EMBEDDING_MODELS.items():
        assert entry.name == name
        assert (entry.input_kind, entry.dimension, entry.node_module) == ("structures", 128, "node_dec")
        assert entry.weights_env == "TRAINED_MODEL_FP" and entry.extra_args == {"fileformat": "guess"}
        assert set(entry.formats) == {"pdb", "cif", "mmcif"}


def test_every_entry_declares_a_complete_contract_and_an_adapter_that_exists(pipe):
    from biopipelines.embedding import EMBEDDING_MODELS

    for entry in EMBEDDING_MODELS.values():
        assert entry.input_kind in ("structures", "sequences", "compounds")
        assert entry.formats and all(isinstance(f, str) and f for f in entry.formats)
        assert entry.weights_env.isupper() and entry.dimension > 0 and isinstance(entry.dimension, int)
        assert isinstance(entry.upstream_default, str) and isinstance(entry.node_module, str)
        assert isinstance(entry.extra_args, dict)
        assert entry.adapter in pipe.ADAPTERS


def test_where_upstream_ships_weights_is_declared_for_every_model():
    from biopipelines.embedding import EMBEDDING_MODELS

    assert "Zenodo" in EMBEDDING_MODELS["progres"].upstream_default
    assert "CIRPIN_model_5k_cp_epoch301.pt" in EMBEDDING_MODELS["cirpin"].upstream_default


def test_an_entry_is_frozen():
    from biopipelines.embedding import EMBEDDING_MODELS

    with pytest.raises(Exception):
        EMBEDDING_MODELS["progres"].dimension = 64


def test_an_unknown_model_names_the_registered_ones():
    from biopipelines.embedding import registry_entry

    for bad in ("esm2", "", None, "Progres"):
        with pytest.raises(ValueError, match=r"model must be one of \('progres', 'cirpin'\)"):
            registry_entry(bad)


def _registered_name_literals(path):
    """(line, name) of every string constant equal to a registered model name, with the ones that are the `name=` of an EmbeddingModel call marked."""
    from biopipelines.embedding import EMBEDDING_MODELS

    tree = ast.parse(open(path, encoding="utf-8").read())
    parents = {child: parent for parent in ast.walk(tree) for child in ast.iter_child_nodes(parent)}
    found = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value in EMBEDDING_MODELS:
            parent = parents[node]
            declared = (isinstance(parent, ast.keyword) and parent.arg == "name"
                        and isinstance(parents[parent], ast.Call) and getattr(parents[parent].func, "id", "") == "EmbeddingModel")
            found.append((node.lineno, node.value, declared))
    return found


def test_no_model_is_named_in_the_code_outside_the_registry():
    wrapper = _registered_name_literals(WRAPPER)
    assert sorted(name for _, name, declared in wrapper if declared) == ["cirpin", "progres"]
    assert [x for x in wrapper if not x[2]] == []
    assert _registered_name_literals(SCRIPT) == []


# ── cosine similarity ────────────────────────────────────────────────────────

@pytest.mark.parametrize("a,b,expected", [
    ([1, 0, 0], [1, 0, 0], 1.0),
    ([1, 0, 0], [0, 1, 0], 0.0),
    ([1, 0, 0], [-1, 0, 0], -1.0),
    ([1, 0], [1, 1], math.sqrt(2) / 2),
    ([3, 4], [4, 3], 24 / 25),
    ([2, 0, 0], [5, 0, 0], 1.0),
])
def test_cosine_of_known_vectors(pipe, a, b, expected):
    assert pipe.cosine_similarity(a, b) == pytest.approx(expected)
    assert pipe.cosine_similarity(b, a) == pytest.approx(expected)


def test_cosine_does_not_depend_on_the_length_of_the_vectors(pipe):
    rng = np.random.default_rng(0)
    a, b = rng.normal(size=128), rng.normal(size=128)
    assert pipe.cosine_similarity(7.5 * a, 0.01 * b) == pytest.approx(pipe.cosine_similarity(a, b))
    assert pipe.cosine_similarity(a, b) == pytest.approx(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def test_cosine_of_a_vector_with_itself_is_not_above_one_through_rounding(pipe):
    v = np.random.default_rng(1).normal(size=2560).astype(np.float32)
    assert pipe.cosine_similarity(v, v) <= 1.0


@pytest.mark.parametrize("a,b", [([0, 0], [1, 1]), ([1, 1], [0, 0]), ([1, float("nan")], [1, 1]), ([1, float("inf")], [1, 1])])
def test_cosine_of_a_zero_or_non_finite_vector_is_refused(pipe, a, b):
    with pytest.raises(ValueError, match="zero or non-finite"):
        pipe.cosine_similarity(a, b)


@pytest.mark.parametrize("a,b", [([1, 2, 3], [1, 2]), ([[1, 2]], [[1, 2]]), (1.0, 2.0)])
def test_cosine_of_vectors_of_different_shape_is_refused(pipe, a, b):
    with pytest.raises(ValueError, match="one length"):
        pipe.cosine_similarity(a, b)


def test_an_embedding_must_have_the_registered_length_and_finite_values(pipe):
    assert pipe.check_embedding([0.5, 0.5], 2).dtype == np.float32
    with pytest.raises(ValueError, match=r"shape \(3,\).*\(2,\)"):
        pipe.check_embedding([1, 2, 3], 2)
    with pytest.raises(ValueError, match="non-finite"):
        pipe.check_embedding([1.0, float("nan")], 2)


# ── pairs ────────────────────────────────────────────────────────────────────

def _write_pairs(tmp_path, first, second, extra=None):
    import pandas as pd

    frame = pd.DataFrame({"query": first, "target": second})
    if extra is not None:
        frame["note"] = extra
    path = tmp_path / "pairs.csv"
    frame.to_csv(path, index=False)
    return str(path)


def test_pairs_are_read_in_first_seen_order_and_a_repeated_pair_is_one(pipe, tmp_path):
    path = _write_pairs(tmp_path, ["a", "a", "c", "b"], ["b", "b", "d", "a"])
    assert pipe.read_pairs(f"TABLE_REFERENCE:{path}:query", f"TABLE_REFERENCE:{path}:target") == [
        ("a", "b"), ("c", "d"), ("b", "a")]


def test_the_pair_reader_is_the_one_the_aligners_use(pipe):
    from biopipelines import pair_alignment

    assert pipe.read_pairs is pair_alignment.read_pairs and pipe.pair_id is pair_alignment.pair_id


def test_ids_that_look_numeric_or_like_na_are_read_as_the_text_they_are(pipe, tmp_path):
    """A plain read_csv turned `007` into `7.0` and `NA` into an empty cell, so the pair of two ids that exist failed under another name."""
    path = tmp_path / "pairs.csv"
    path.write_text("query,target\n007,12\nNA,x1\nnull,None\n")
    pairs = pipe.read_pairs(f"TABLE_REFERENCE:{path}:query", f"TABLE_REFERENCE:{path}:target")
    assert pairs == [("007", "12"), ("NA", "x1"), ("null", "None")]
    rows, failures = pipe.pair_rows(pairs, {"007": unit(1, 0), "12": unit(1, 1), "NA": unit(0, 1), "x1": unit(1, 0)},
                                    {"007", "12", "NA", "x1"}, {"007", "12", "NA", "x1"})
    assert [r["id"] for r in rows[:2]] == ["007+12", "NA+x1"]
    assert not any(math.isnan(r["cosine_similarity"]) for r in rows[:2])
    assert [f[0] for f in failures] == ["null+None"] and "'null' is not in the input stream" in failures[0][2]


def test_an_empty_cell_of_the_pair_table_is_an_empty_id(pipe, tmp_path):
    path = _write_pairs(tmp_path, ["a", None], [None, "b"])
    assert pipe.read_pairs(f"TABLE_REFERENCE:{path}:query", f"TABLE_REFERENCE:{path}:target") == [("a", ""), ("", "b")]


def test_two_pairs_that_would_share_a_similarities_id_are_refused(pipe, tmp_path):
    path = _write_pairs(tmp_path, ["a", "a+b"], ["b+c", "c"])
    with pytest.raises(ValueError, match=r"both have the id 'a\+b\+c'"):
        pipe.read_pairs(f"TABLE_REFERENCE:{path}:query", f"TABLE_REFERENCE:{path}:target")


def test_a_pair_table_without_the_named_column_lists_its_columns(pipe, tmp_path):
    path = _write_pairs(tmp_path, ["a"], ["b"])
    with pytest.raises(ValueError, match=r"Column 'absent' is not in the pair table.*\['query', 'target'\]"):
        pipe.read_pairs(f"TABLE_REFERENCE:{path}:query", f"TABLE_REFERENCE:{path}:absent")


def test_a_missing_pair_table_is_refused(pipe, tmp_path):
    ref = f"TABLE_REFERENCE:{tmp_path / 'absent.csv'}:query"
    with pytest.raises(FileNotFoundError):
        pipe.read_pairs(ref, ref)


def unit(*values):
    v = np.array(values, dtype=np.float32)
    return v / np.linalg.norm(v)


def test_every_pair_row_is_kept_and_a_member_without_an_embedding_gives_nan_and_its_reason(pipe):
    vectors = {"a": unit(1, 0), "b": unit(1, 1), "c": unit(0, 1)}
    declared = {"a", "b", "c", "failed", "filtered"}
    attempted = {"a", "b", "c", "failed"}
    pairs = [("a", "b"), ("a", "ghost"), ("failed", "c"), ("c", "filtered"), ("", "a"), ("b", "a")]
    rows, failures = pipe.pair_rows(pairs, vectors, declared, attempted)
    assert [r["id"] for r in rows] == ["a+b", "a+ghost", "failed+c", "c+filtered", "+a", "b+a"]
    assert [list(r) for r in rows] == [pipe.SIMILARITY_COLUMNS] * 6
    sims = [r["cosine_similarity"] for r in rows]
    assert sims[0] == pytest.approx(math.sqrt(2) / 2) and sims[5] == sims[0]
    assert all(math.isnan(s) for s in sims[1:5])
    causes = {pid: cause for pid, _kind, cause in failures}
    assert causes["a+ghost"] == "'ghost' is not in the input stream"
    assert causes["failed+c"] == "'failed' failed to embed"
    assert causes["c+filtered"] == "'filtered' was dropped upstream (no file)"
    assert causes["+a"] == "a pair has an empty id"
    assert len(failures) == 4


def test_a_member_with_no_partner_in_the_stream_is_unmatched_and_one_that_failed_to_embed_is_a_failure(pipe):
    vectors = {"a": unit(1, 0), "c": unit(0, 1)}
    pairs = [("a", "ghost"), ("failed", "c"), ("c", "filtered"), ("", "a")]
    _, failures = pipe.pair_rows(pairs, vectors, {"a", "c", "failed", "filtered"}, {"a", "c", "failed"})
    assert {pid: kind for pid, kind, _cause in failures} == {
        "a+ghost": "unmatched", "failed+c": "failure", "c+filtered": "unmatched", "+a": "unmatched"}


def test_a_pair_of_a_structure_with_itself_is_one(pipe):
    rows, failures = pipe.pair_rows([("a", "a")], {"a": unit(3, 4)}, {"a"}, {"a"})
    assert rows[0]["cosine_similarity"] == pytest.approx(1.0) and failures == []


# ── the step loop, against a stand-in embedder ───────────────────────────────

DIM = 4


class FakeEmbedder:
    """Vectors by file name; `bad` raises, `notfinite` returns a non-finite vector, `wide` the wrong length, `nonodes` no nodes."""

    created = []
    VECTORS = {"a": unit(1, 0, 0, 0), "b": unit(1, 1, 0, 0), "c": unit(0, 0, 1, 0)}

    def __init__(self, entry, weights, device, want_nodes):
        self.device = "cpu"
        FakeEmbedder.created.append((entry.name, weights, device, want_nodes))
        self.want_nodes = want_nodes

    def embed(self, path):
        name = os.path.splitext(os.path.basename(path))[0]
        if name == "bad":
            raise ValueError("no C-alpha atoms")
        if name == "notfinite":
            return np.array([1, float("nan"), 0, 0], dtype=np.float32), np.zeros((5, 3), dtype=np.float32)
        if name == "wide":
            return np.ones(DIM + 1, dtype=np.float32), np.zeros((5, 3), dtype=np.float32)
        nodes = None if name == "nonodes" else np.full((len(name) + 5, 3), 1.0, dtype=np.float32)
        return self.VECTORS.get(name, unit(1, 1, 1, 1)), nodes


@pytest.fixture
def fake(pipe, monkeypatch, tmp_path):
    from biopipelines.embedding import EMBEDDING_MODELS, EmbeddingModel

    FakeEmbedder.created = []
    monkeypatch.setitem(EMBEDDING_MODELS, "fake", EmbeddingModel(
        name="fake", input_kind="structures", formats=("pdb",), adapter="fake_adapter",
        weights_env="FAKE_MODEL_WEIGHTS", upstream_default="nowhere", dimension=DIM, node_module="nodes"))
    monkeypatch.setitem(pipe.ADAPTERS, "fake_adapter", FakeEmbedder)
    weights = tmp_path / "weights.pt"
    weights.write_bytes(b"weights")
    monkeypatch.setenv("FAKE_MODEL_WEIGHTS", str(weights))
    return FakeEmbedder


def write_inputs(tmp_path, ids, node_embeddings=False, pairs=None, device="cpu"):
    import json
    import pandas as pd

    files = []
    for sid in ids:
        path = tmp_path / "in" / f"{sid}.pdb"
        path.parent.mkdir(exist_ok=True)
        path.write_text("ATOM")
        files.append(str(path))
    map_csv = tmp_path / "structures_map.csv"
    pd.DataFrame({"id": ids, "file": files}).to_csv(map_csv, index=False)
    stream = tmp_path / "inputs.json"
    stream.write_text(json.dumps({"name": "structures", "ids": ids, "files": files, "map_table": str(map_csv), "format": "pdb"}))
    config = {"model": "fake", "node_embeddings": node_embeddings, "device": device, "first": None, "second": None}
    if pairs is not None:
        table = _write_pairs(tmp_path, [p[0] for p in pairs], [p[1] for p in pairs])
        config["first"], config["second"] = f"TABLE_REFERENCE:{table}:query", f"TABLE_REFERENCE:{table}:target"
    config_yaml = tmp_path / "config.yaml"
    import yaml
    config_yaml.write_text(yaml.safe_dump(config))
    out = tmp_path / "007_Embedding"
    argv = [
        "pipe_embedding.py", "--config-yaml", str(config_yaml), "--inputs-json", str(stream),
        "--embeddings-dir", str(out / "embeddings"), "--embeddings-map-csv", str(out / "embeddings" / "embeddings_map.csv"),
        "--node-embeddings-dir", str(out / "node_embeddings"),
        "--node-embeddings-map-csv", str(out / "node_embeddings" / "node_embeddings_map.csv"),
        "--similarities-csv", str(out / "tables" / "similarities.csv"),
        "--local-missing-csv", str(out / "_execution" / "local_missing.csv"),
    ]
    return argv, out


def run_main(pipe, monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", argv)
    pipe.main()


def test_one_npy_per_id_of_the_registered_length_and_a_map_row_for_each(pipe, fake, monkeypatch, tmp_path):
    import pandas as pd

    argv, out = write_inputs(tmp_path, ["a", "b", "c"])
    run_main(pipe, monkeypatch, argv)

    emb_map = pd.read_csv(out / "embeddings" / "embeddings_map.csv")
    assert list(emb_map.columns) == ["id", "file"] and list(emb_map["id"]) == ["a", "b", "c"]
    assert sorted(os.listdir(out / "embeddings")) == ["a.npy", "b.npy", "c.npy", "embeddings_map.csv"]
    for sid, path in zip(emb_map["id"], emb_map["file"]):
        vector = np.load(path)
        assert vector.shape == (DIM,) and vector.dtype == np.float32
        assert np.allclose(vector, FakeEmbedder.VECTORS[sid])
    assert not (out / "node_embeddings").exists() and not (out / "tables" / "similarities.csv").exists()
    assert len(pd.read_csv(out / "_execution" / "local_missing.csv")) == 0


def test_the_weights_come_from_the_variable_the_registry_names_and_the_model_gets_its_arguments(pipe, fake, monkeypatch, tmp_path):
    argv, _ = write_inputs(tmp_path, ["a"], node_embeddings=True, device="cuda")
    run_main(pipe, monkeypatch, argv)
    assert fake.created == [("fake", os.environ["FAKE_MODEL_WEIGHTS"], "cuda", True)]


def test_node_embeddings_are_one_npy_per_id_of_shape_residues_by_dimension(pipe, fake, monkeypatch, tmp_path):
    import pandas as pd

    argv, out = write_inputs(tmp_path, ["a", "bcd"], node_embeddings=True)
    run_main(pipe, monkeypatch, argv)

    node_map = pd.read_csv(out / "node_embeddings" / "node_embeddings_map.csv")
    assert list(node_map["id"]) == ["a", "bcd"]
    assert sorted(os.listdir(out / "node_embeddings")) == ["a.npy", "bcd.npy", "node_embeddings_map.csv"]
    shapes = {sid: np.load(path).shape for sid, path in zip(node_map["id"], node_map["file"])}
    assert shapes == {"a": (6, 3), "bcd": (8, 3)}
    assert np.load(node_map.loc[0, "file"]).dtype == np.float32


def test_a_failed_id_has_no_file_no_map_row_and_a_failure_row_while_the_others_go_on(pipe, fake, monkeypatch, tmp_path):
    import pandas as pd

    argv, out = write_inputs(tmp_path, ["a", "bad", "c"])
    run_main(pipe, monkeypatch, argv)

    assert list(pd.read_csv(out / "embeddings" / "embeddings_map.csv")["id"]) == ["a", "c"]
    assert sorted(os.listdir(out / "embeddings")) == ["a.npy", "c.npy", "embeddings_map.csv"]
    missing = pd.read_csv(out / "_execution" / "local_missing.csv")
    assert list(missing.columns) == ["id", "removed_by", "kind", "cause"]
    assert list(missing["id"]) == ["bad"] and missing.loc[0, "kind"] == "failure" and "no C-alpha" in missing.loc[0, "cause"]
    assert missing.loc[0, "removed_by"] == "007_Embedding"


@pytest.mark.parametrize("name,cause", [("notfinite", "non-finite"), ("wide", "registry declares")])
def test_a_vector_of_the_wrong_length_or_with_non_finite_values_fails_that_id(pipe, fake, monkeypatch, tmp_path, name, cause):
    import pandas as pd

    argv, out = write_inputs(tmp_path, ["a", name])
    run_main(pipe, monkeypatch, argv)
    missing = pd.read_csv(out / "_execution" / "local_missing.csv")
    assert list(missing["id"]) == [name] and cause in missing.loc[0, "cause"]
    assert not (out / "embeddings" / f"{name}.npy").exists()


def test_node_embeddings_wanted_but_not_returned_fail_the_id_and_write_nothing_for_it(pipe, fake, monkeypatch, tmp_path):
    import pandas as pd

    argv, out = write_inputs(tmp_path, ["a", "nonodes"], node_embeddings=True)
    run_main(pipe, monkeypatch, argv)
    assert list(pd.read_csv(out / "_execution" / "local_missing.csv")["id"]) == ["nonodes"]
    assert not (out / "embeddings" / "nonodes.npy").exists() and not (out / "node_embeddings" / "nonodes.npy").exists()
    assert list(pd.read_csv(out / "node_embeddings" / "node_embeddings_map.csv")["id"]) == ["a"]


def test_the_similarity_table_has_a_row_per_pair_and_failed_pairs_reach_the_missing_table(pipe, fake, monkeypatch, tmp_path):
    import pandas as pd

    pairs = [("a", "b"), ("a", "c"), ("b", "ghost"), ("bad", "a"), ("a", "b")]
    argv, out = write_inputs(tmp_path, ["a", "b", "c", "bad"], pairs=pairs)
    run_main(pipe, monkeypatch, argv)

    sims = pd.read_csv(out / "tables" / "similarities.csv")
    assert list(sims.columns) == ["id", "first", "second", "cosine_similarity"]
    assert list(sims["id"]) == ["a+b", "a+c", "b+ghost", "bad+a"]
    assert sims["id"].is_unique
    assert sims.loc[0, "cosine_similarity"] == pytest.approx(math.sqrt(2) / 2)
    assert sims.loc[1, "cosine_similarity"] == pytest.approx(0.0)
    assert sims["cosine_similarity"].isna().tolist() == [False, False, True, True]

    missing = pd.read_csv(out / "_execution" / "local_missing.csv")
    assert list(missing["id"]) == ["bad", "b+ghost", "bad+a"]
    assert dict(zip(missing["id"], missing["kind"])) == {"bad": "failure", "b+ghost": "unmatched", "bad+a": "failure"}
    assert "'ghost' is not in the input stream" in list(missing["cause"])


def test_pair_ids_that_look_numeric_survive_the_step(pipe, fake, monkeypatch, tmp_path):
    """Through the step, not just the reader: `007` and `NA` are ids of the input stream and their pair is compared."""
    import pandas as pd

    argv, out = write_inputs(tmp_path, ["007", "NA", "b"], pairs=[("007", "NA"), ("NA", "b")])
    run_main(pipe, monkeypatch, argv)
    sims = pd.read_csv(out / "tables" / "similarities.csv", dtype=str, keep_default_na=False)
    assert list(sims["id"]) == ["007+NA", "NA+b"] and list(sims["first"]) == ["007", "NA"]
    assert (sims["cosine_similarity"] != "").all()
    assert len(pd.read_csv(out / "_execution" / "local_missing.csv")) == 0


def test_a_pair_of_structures_that_both_exist_is_not_a_failure_even_if_another_id_failed(pipe, fake, monkeypatch, tmp_path):
    import pandas as pd

    argv, out = write_inputs(tmp_path, ["a", "bad", "b"], pairs=[("a", "b")])
    run_main(pipe, monkeypatch, argv)
    assert list(pd.read_csv(out / "_execution" / "local_missing.csv")["id"]) == ["bad"]


def test_a_bad_pair_table_stops_the_step_before_anything_is_embedded(pipe, fake, monkeypatch, tmp_path):
    argv, out = write_inputs(tmp_path, ["a"], pairs=[("a", "a")])
    import yaml
    config = yaml.safe_load(open(argv[2]))
    config["second"] = config["second"].rsplit(":", 1)[0] + ":absent"
    open(argv[2], "w").write(yaml.safe_dump(config))
    with pytest.raises(ValueError, match="not in the pair table"):
        run_main(pipe, monkeypatch, argv)
    assert fake.created == [] and not (out / "embeddings").exists()


def test_without_weights_the_step_stops_before_a_model_is_built_and_says_where_upstream_keeps_them(pipe, fake, monkeypatch, tmp_path):
    monkeypatch.delenv("FAKE_MODEL_WEIGHTS")
    argv, _ = write_inputs(tmp_path, ["a"])
    with pytest.raises(FileNotFoundError, match="FAKE_MODEL_WEIGHTS is not set.*Upstream: nowhere"):
        run_main(pipe, monkeypatch, argv)
    assert fake.created == []


def test_weights_that_are_not_a_file_are_refused(pipe, fake, monkeypatch, tmp_path):
    monkeypatch.setenv("FAKE_MODEL_WEIGHTS", str(tmp_path / "absent.pt"))
    argv, _ = write_inputs(tmp_path, ["a"])
    with pytest.raises(FileNotFoundError, match="is not a file"):
        run_main(pipe, monkeypatch, argv)


def test_every_id_failing_is_a_step_failure_but_the_tables_are_still_written(pipe, fake, monkeypatch, tmp_path):
    import pandas as pd

    argv, out = write_inputs(tmp_path, ["bad"], pairs=[("bad", "bad")])
    with pytest.raises(SystemExit) as caught:
        run_main(pipe, monkeypatch, argv)
    assert caught.value.code == 1
    assert list(pd.read_csv(out / "tables" / "similarities.csv")["id"]) == ["bad+bad"]
    assert list(pd.read_csv(out / "_execution" / "local_missing.csv")["id"]) == ["bad", "bad+bad"]


def test_a_stream_whose_files_are_all_absent_is_a_step_failure_with_the_tables_written(pipe, fake, monkeypatch, tmp_path):
    import pandas as pd

    argv, out = write_inputs(tmp_path, ["a", "b"], pairs=[("a", "b")])
    for name in ("a", "b"):
        (tmp_path / "in" / f"{name}.pdb").unlink()
    with pytest.raises(SystemExit) as caught:
        run_main(pipe, monkeypatch, argv)
    assert caught.value.code == 1
    assert fake.created and len(pd.read_csv(out / "embeddings" / "embeddings_map.csv")) == 0
    sims = pd.read_csv(out / "tables" / "similarities.csv")
    assert list(sims["id"]) == ["a+b"] and sims["cosine_similarity"].isna().all()
    assert "was dropped upstream" in pd.read_csv(out / "_execution" / "local_missing.csv").loc[0, "cause"]


def test_the_weights_digest_names_the_file_that_ran(pipe, tmp_path, capsys, monkeypatch, fake):
    import hashlib

    argv, _ = write_inputs(tmp_path, ["a"])
    run_main(pipe, monkeypatch, argv)
    expected = hashlib.sha256(b"weights").hexdigest()[:12]
    assert pipe.weights_digest(os.environ["FAKE_MODEL_WEIGHTS"]) == expected
    assert expected in capsys.readouterr().out


@pytest.mark.parametrize("requested,cuda,expected", [("auto", True, "cuda"), ("auto", False, "cpu"), ("cpu", True, "cpu"), ("cuda", True, "cuda")])
def test_the_device_is_resolved(pipe, requested, cuda, expected):
    assert pipe.resolve_device(requested, cuda) == expected


def test_cuda_without_a_gpu_fails_instead_of_running_on_cpu(pipe):
    with pytest.raises(RuntimeError, match="no CUDA device"):
        pipe.resolve_device("cuda", False)
    with pytest.raises(ValueError, match="device"):
        pipe.resolve_device("tpu", True)


# ── the Progres adapter, over a torch-only stand-in of the progres package ───

class _StandInModel:
    """Built lazily so torch is only needed when the adapter tests run."""

    @staticmethod
    def build(torch):
        class Model(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.node_enc = torch.nn.Linear(2, 3)
                self.node_dec = torch.nn.Linear(3, 3)
                self.graph_dec = torch.nn.Linear(3, 4)

            def forward(self, data):
                feats = self.node_dec(self.node_enc(data.x)).unsqueeze(0)
                graph = self.graph_dec(feats.squeeze(0).sum(dim=0, keepdim=True))
                return torch.nn.functional.normalize(graph, dim=1)

        return Model


@pytest.fixture
def stand_in_progres(monkeypatch):
    torch = pytest.importorskip("torch")
    calls = []

    class Graph:
        def __init__(self, x):
            self.x = x

        def to(self, device):
            calls.append(("to", device))
            return self

    def read_coords(path, **kwargs):
        calls.append(("read_coords", os.path.basename(path), kwargs))
        n = int(os.path.splitext(os.path.basename(path))[0].split("_")[-1])
        return [[float(i), 0.0, 0.0] for i in range(n)]

    def coords_to_graph(coords):
        return Graph(torch.arange(len(coords) * 2, dtype=torch.float32).reshape(len(coords), 2))

    module = types.ModuleType("progres.progres")
    module.Model, module.read_coords, module.coords_to_graph = _StandInModel.build(torch), read_coords, coords_to_graph
    package = types.ModuleType("progres")
    package.progres = module
    monkeypatch.setitem(sys.modules, "progres", package)
    monkeypatch.setitem(sys.modules, "progres.progres", module)
    return torch, module, calls


def save_checkpoint(torch, module, path, wrapped=True):
    state = module.Model().state_dict()
    torch.save({"epoch": 7, "model": state} if wrapped else state, path)
    return str(path)


def test_the_adapter_loads_a_checkpoint_and_returns_a_unit_vector_and_the_per_node_output(pipe, stand_in_progres, tmp_path):
    from biopipelines.embedding import EMBEDDING_MODELS

    torch, module, calls = stand_in_progres
    weights = save_checkpoint(torch, module, tmp_path / "w.pt")
    embedder = pipe.build_embedder(EMBEDDING_MODELS["cirpin"], weights, "cpu", True)
    vector, nodes = embedder.embed(str(tmp_path / "s_12.pdb"))
    assert vector.shape == (4,) and np.linalg.norm(vector) == pytest.approx(1.0, abs=1e-6)
    assert nodes.shape == (12, 3) and nodes.dtype == np.float32
    graph = ("read_coords", "s_12.pdb", {"fileformat": "guess"})
    assert graph in calls and ("to", "cpu") in calls


def test_the_per_node_output_is_the_one_the_graph_decoder_pools(pipe, stand_in_progres, tmp_path):
    from biopipelines.embedding import EMBEDDING_MODELS

    torch, module, _ = stand_in_progres
    embedder = pipe.build_embedder(EMBEDDING_MODELS["progres"], save_checkpoint(torch, module, tmp_path / "w.pt"), "cpu", True)
    vector, nodes = embedder.embed(str(tmp_path / "s_9.pdb"))
    with torch.no_grad():
        again = torch.nn.functional.normalize(embedder.model.graph_dec(torch.tensor(nodes).sum(0, keepdim=True)), dim=1)[0].numpy()
    assert np.allclose(again, vector, atol=1e-6)


def test_without_node_embeddings_no_hook_is_attached_and_none_are_returned(pipe, stand_in_progres, tmp_path):
    from biopipelines.embedding import EMBEDDING_MODELS

    torch, module, _ = stand_in_progres
    embedder = pipe.build_embedder(EMBEDDING_MODELS["progres"], save_checkpoint(torch, module, tmp_path / "w.pt"), "cpu", False)
    assert len(embedder.model.node_dec._forward_hooks) == 0
    assert embedder.embed(str(tmp_path / "s_9.pdb"))[1] is None


def test_the_same_embedder_gives_each_structure_its_own_nodes(pipe, stand_in_progres, tmp_path):
    from biopipelines.embedding import EMBEDDING_MODELS

    torch, module, _ = stand_in_progres
    embedder = pipe.build_embedder(EMBEDDING_MODELS["progres"], save_checkpoint(torch, module, tmp_path / "w.pt"), "cpu", True)
    assert embedder.embed(str(tmp_path / "s_5.pdb"))[1].shape == (5, 3)
    assert embedder.embed(str(tmp_path / "s_8.pdb"))[1].shape == (8, 3)


@pytest.mark.parametrize("residues", [0, 1, 2])
def test_fewer_than_three_residues_fail_with_a_message_that_says_what_to_check(pipe, stand_in_progres, tmp_path, residues):
    from biopipelines.embedding import EMBEDDING_MODELS

    torch, module, _ = stand_in_progres
    embedder = pipe.build_embedder(EMBEDDING_MODELS["progres"], save_checkpoint(torch, module, tmp_path / "w.pt"), "cpu", False)
    with pytest.raises(ValueError, match=rf"{residues} C-alpha atom\(s\) in the first chain.*at least 3"):
        embedder.embed(str(tmp_path / f"s_{residues}.pdb"))


def test_a_checkpoint_that_is_not_a_model_dict_is_refused(pipe, stand_in_progres, tmp_path):
    from biopipelines.embedding import EMBEDDING_MODELS

    torch, module, _ = stand_in_progres
    bare = save_checkpoint(torch, module, tmp_path / "bare.pt", wrapped=False)
    with pytest.raises(ValueError, match="not a Progres-family checkpoint"):
        pipe.build_embedder(EMBEDDING_MODELS["progres"], bare, "cpu", False)


def test_a_checkpoint_of_another_architecture_fails_at_load_not_silently(pipe, stand_in_progres, tmp_path):
    from biopipelines.embedding import EMBEDDING_MODELS

    torch, module, _ = stand_in_progres
    other = tmp_path / "other.pt"
    torch.save({"epoch": 1, "model": {"node_enc.weight": torch.zeros(5, 5)}}, other)
    with pytest.raises(RuntimeError, match="Missing key|size mismatch|Unexpected key"):
        pipe.build_embedder(EMBEDDING_MODELS["cirpin"], str(other), "cpu", False)


def test_the_checkpoint_is_read_with_weights_only(pipe, stand_in_progres, tmp_path):
    import pickle
    from biopipelines.embedding import EMBEDDING_MODELS

    torch, module, _ = stand_in_progres

    class Payload:
        def __reduce__(self):
            return (print, ("arbitrary code ran",))

    hostile = tmp_path / "hostile.pt"
    torch.save({"epoch": 1, "model": module.Model().state_dict(), "extra": Payload()}, hostile)
    with pytest.raises(pickle.UnpicklingError):
        pipe.build_embedder(EMBEDDING_MODELS["progres"], str(hostile), "cpu", False)


def test_a_node_module_the_network_lacks_is_an_error_when_nodes_are_wanted(pipe, stand_in_progres, tmp_path, monkeypatch):
    import dataclasses
    from biopipelines.embedding import EMBEDDING_MODELS

    torch, module, _ = stand_in_progres
    entry = dataclasses.replace(EMBEDDING_MODELS["progres"], node_module="no_such_module")
    with pytest.raises(AttributeError, match="no_such_module"):
        pipe.build_embedder(entry, save_checkpoint(torch, module, tmp_path / "w.pt"), "cpu", True)


# ── the wrapper ──────────────────────────────────────────────────────────────

def make_weights(tmp_path, name="weights.pt"):
    path = tmp_path / name
    path.write_bytes(b"weights")
    return str(path)


def pair_columns(tmp_path):
    from biopipelines.base_config import TableInfo

    table = TableInfo("pairs", _write_pairs(tmp_path, ["a"], ["b"]), columns=["query", "target"])
    return table.query, table.target


def build(new_pipeline, tmp_path, **overrides):
    from biopipelines.embedding import Embedding
    from biopipelines.mock import Mock

    args = {"model": "cirpin", "weights": make_weights(tmp_path)}
    args.update(overrides)
    pipeline = new_pipeline("embedding")
    with pipeline:
        source = Mock(ids=["s1", "s2"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        Embedding(source, **args)
        script_path = pipeline.save()
    return pipeline.tools[-1], script_path


def test_defaults(local_config, isolated_cwd, new_pipeline, tmp_path):
    tool, _ = build(new_pipeline, tmp_path)
    assert (tool.model, tool.node_embeddings, tool.device, tool.first, tool.second) == ("cirpin", False, "cpu", None, None)
    assert tool.inputs_stream.name == "structures" and list(tool.inputs_stream.ids) == ["s1", "s2"]


def test_the_model_has_no_default(local_config, isolated_cwd, new_pipeline, tmp_path):
    from biopipelines.embedding import Embedding
    from biopipelines.mock import Mock

    with new_pipeline("embedding_nomodel"):
        source = Mock(ids=["s1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        with pytest.raises(TypeError):
            Embedding(source, weights=make_weights(tmp_path))


@pytest.mark.parametrize("model", ["progres", "cirpin"])
def test_every_registered_model_is_accepted(local_config, isolated_cwd, new_pipeline, tmp_path, model):
    assert build(new_pipeline, tmp_path, model=model)[0].model == model


@pytest.mark.parametrize("overrides,message", [
    ({"model": "esm2"}, "model must be one of"),
    ({"model": None}, "model must be one of"),
    ({"weights": ""}, "weights must be the path"),
    ({"weights": "   "}, "weights must be the path"),
    ({"weights": None}, "weights must be the path"),
    ({"weights": "https://example.org/model.pt"}, "not a URL"),
    ({"weights": "s3://bucket/model.pt"}, "not a URL"),
    ({"node_embeddings": 1}, "node_embeddings must be True or False"),
    ({"device": "gpu"}, "device must be one of"),
])
def test_bad_parameters_fail_at_configuration(local_config, isolated_cwd, new_pipeline, tmp_path, overrides, message):
    with pytest.raises(ValueError, match=message):
        build(new_pipeline, tmp_path, **overrides)


def test_weights_that_do_not_exist_are_named_with_where_upstream_keeps_them(local_config, isolated_cwd, new_pipeline, tmp_path, monkeypatch):
    import dataclasses
    from biopipelines.embedding import EMBEDDING_MODELS

    with pytest.raises(ValueError, match="weights file not found: .*absent.pt.*CIRPIN_model_5k_cp_epoch301.pt"):
        build(new_pipeline, tmp_path, model="cirpin", weights=str(tmp_path / "absent.pt"))
    with pytest.raises(ValueError, match="Zenodo"):
        build(new_pipeline, tmp_path, model="progres", weights=str(tmp_path / "absent.pt"))
    monkeypatch.setitem(EMBEDDING_MODELS, "bare", dataclasses.replace(EMBEDDING_MODELS["cirpin"], name="bare", upstream_default=""))
    with pytest.raises(ValueError, match="Upstream names no default file"):
        build(new_pipeline, tmp_path, model="bare", weights=str(tmp_path / "absent.pt"))


def test_a_directory_is_not_a_weights_file(local_config, isolated_cwd, new_pipeline, tmp_path):
    with pytest.raises(ValueError, match="weights file not found"):
        build(new_pipeline, tmp_path, weights=str(tmp_path))


def test_a_relative_weights_path_becomes_absolute(local_config, isolated_cwd, new_pipeline, tmp_path):
    make_weights(tmp_path, "rel.pt")
    tool, _ = build(new_pipeline, tmp_path, weights="rel.pt")
    assert os.path.isabs(tool._weights_path()) and os.path.samefile(tool._weights_path(), tmp_path / "rel.pt")


@pytest.mark.parametrize("fmt", ["pdb", "cif", "mmcif", "pdb|cif"])
def test_pdb_and_mmcif_streams_are_accepted(local_config, isolated_cwd, new_pipeline, tmp_path, fmt):
    from biopipelines.embedding import Embedding
    from biopipelines.mock import Mock

    with new_pipeline("embedding_fmt"):
        source = Mock(ids=["s1"], streams={"structures": {"format": fmt, "file": "<id>.x"}})
        Embedding(source, model="progres", weights=make_weights(tmp_path))


@pytest.mark.parametrize("fmt", ["sdf", "csv", "pdb|sdf"])
def test_other_formats_are_refused_with_the_accepted_ones(local_config, isolated_cwd, new_pipeline, tmp_path, fmt):
    from biopipelines.embedding import Embedding
    from biopipelines.mock import Mock

    with new_pipeline("embedding_badfmt"):
        source = Mock(ids=["s1"], streams={"structures": {"format": fmt, "file": "<id>.x"}})
        with pytest.raises(ValueError, match=r"reads structures in \('pdb', 'cif', 'mmcif'\)"):
            Embedding(source, model="progres", weights=make_weights(tmp_path))


def test_a_bare_value_is_not_an_input(local_config, isolated_cwd, new_pipeline, tmp_path):
    from biopipelines.embedding import Embedding
    from biopipelines.mock import Mock

    with new_pipeline("embedding_bare"):
        Mock(ids=["s1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        with pytest.raises(ValueError, match="inputs must be DataStream or StandardizedOutput"):
            Embedding("1UBQ", model="progres", weights=make_weights(tmp_path))


def test_an_output_without_the_stream_the_model_reads_is_refused(local_config, isolated_cwd, new_pipeline, tmp_path):
    from biopipelines.embedding import Embedding
    from biopipelines.mock import Mock

    with new_pipeline("embedding_nostream"):
        source = Mock(ids=["s1"], streams={"sequences": {"format": "csv", "values": ["MKT"]}})
        with pytest.raises(ValueError, match="no 'structures' stream"):
            Embedding(source, model="progres", weights=make_weights(tmp_path))


def test_the_input_stream_is_the_one_the_registry_entry_names(local_config, isolated_cwd, new_pipeline, tmp_path, monkeypatch):
    from biopipelines.embedding import EMBEDDING_MODELS, Embedding, EmbeddingModel
    from biopipelines.mock import Mock

    monkeypatch.setitem(EMBEDDING_MODELS, "seqmodel", EmbeddingModel(
        name="seqmodel", input_kind="sequences", formats=("csv",), adapter="progres_graph",
        weights_env="SEQ_WEIGHTS", upstream_default="", dimension=320))
    pipeline = new_pipeline("embedding_sequences")
    with pipeline:
        source = Mock(ids=["s1", "s2"], streams={"sequences": {"format": "csv", "values": ["MKT", "ACD"]},
                                                  "structures": {"format": "pdb", "file": "<id>.pdb"}})
        Embedding(source, model="seqmodel", weights=make_weights(tmp_path))
        pipeline.save()
    tool = pipeline.tools[-1]
    assert tool.inputs_stream.name == "sequences"
    script = tool.generate_script("unused")
    assert "export SEQ_WEIGHTS=" in script and "TRAINED_MODEL_FP" not in script


def test_a_model_without_a_per_node_output_refuses_node_embeddings(local_config, isolated_cwd, new_pipeline, tmp_path, monkeypatch):
    from biopipelines.embedding import EMBEDDING_MODELS, EmbeddingModel

    monkeypatch.setitem(EMBEDDING_MODELS, "flat", EmbeddingModel(
        name="flat", input_kind="structures", formats=("pdb",), adapter="progres_graph",
        weights_env="FLAT_WEIGHTS", upstream_default="", dimension=8))
    assert build(new_pipeline, tmp_path, model="flat")[0].node_embeddings is False
    with pytest.raises(ValueError, match="no per-node output"):
        build(new_pipeline, tmp_path, model="flat", node_embeddings=True)


# ── pairs at configuration ───────────────────────────────────────────────────

def test_pairs_are_two_columns_of_one_table(local_config, isolated_cwd, new_pipeline, tmp_path):
    first, second = pair_columns(tmp_path)
    tool, _ = build(new_pipeline, tmp_path, first=first, second=second)
    assert (tool.first.column, tool.second.column) == ("query", "target") and tool.first.path == tool.second.path


def test_a_table_and_column_tuple_is_a_column_reference(local_config, isolated_cwd, new_pipeline, tmp_path):
    from biopipelines.base_config import TableInfo

    path = _write_pairs(tmp_path, ["a"], ["b"])
    table = TableInfo("pairs", path, columns=["query", "target"])
    tool, _ = build(new_pipeline, tmp_path, first=(table, "query"), second=(path, "target"))
    assert str(tool.first) == f"TABLE_REFERENCE:{path}:query" and str(tool.second) == f"TABLE_REFERENCE:{path}:target"


@pytest.mark.parametrize("which", ["first", "second"])
def test_one_member_column_alone_is_refused(local_config, isolated_cwd, new_pipeline, tmp_path, which):
    first, second = pair_columns(tmp_path)
    with pytest.raises(ValueError, match="given together or not at all"):
        build(new_pipeline, tmp_path, **{which: first})


def test_columns_of_two_tables_are_refused(local_config, isolated_cwd, new_pipeline, tmp_path):
    from biopipelines.base_config import TableInfo

    one = TableInfo("one", _write_pairs(tmp_path, ["a"], ["b"]), columns=["query"])
    other = tmp_path / "other.csv"
    other.write_text("query\na\n")
    two = TableInfo("two", str(other), columns=["query"])
    with pytest.raises(ValueError, match="same table"):
        build(new_pipeline, tmp_path, first=one.query, second=two.query)


def test_a_literal_string_is_not_a_column_reference(local_config, isolated_cwd, new_pipeline, tmp_path):
    first, _ = pair_columns(tmp_path)
    with pytest.raises(ValueError, match="must be a column of the pair table"):
        build(new_pipeline, tmp_path, first="query", second=first)


# ── the generated script ─────────────────────────────────────────────────────

def test_the_script_exports_the_weights_path_under_the_registered_variable_and_runs_the_helper(local_config, isolated_cwd, new_pipeline, tmp_path):
    import yaml

    first, second = pair_columns(tmp_path)
    weights = make_weights(tmp_path)
    tool, _ = build(new_pipeline, tmp_path, weights=weights, first=first, second=second, node_embeddings=True, device="auto")
    script = tool.generate_script("unused")
    assert f"export TRAINED_MODEL_FP={weights}\n" in script or f"export TRAINED_MODEL_FP='{weights}'\n" in script
    assert script.index("export TRAINED_MODEL_FP") < script.index("pipe_embedding.py")
    assert 'pipe_embedding.py"' in script and "pipe_propagate_missing.py" in script and tool.local_missing_csv in script
    assert yaml.safe_load(open(tool.config_yaml)) == {
        "model": "cirpin", "node_embeddings": True, "device": "auto",
        "first": str(first), "second": str(second)}
    assert os.path.exists(tool.inputs_json)


def test_a_weights_path_with_shell_characters_is_quoted_not_expanded(local_config, isolated_cwd, new_pipeline, tmp_path):
    tricky = tmp_path / "my weights $(touch pwned) `id` \"q\".pt"
    try:
        tricky.write_bytes(b"w")
    except OSError:
        pytest.skip("this file system does not allow such a name")
    tool, _ = build(new_pipeline, tmp_path, weights=str(tricky))
    line = next(l for l in tool.generate_script("unused").splitlines() if l.startswith("export TRAINED_MODEL_FP="))
    import shlex
    assert shlex.split(line)[1] == f"TRAINED_MODEL_FP={tool._weights_path()}"


def test_every_flag_the_wrapper_emits_is_declared_by_the_pipe(local_config, isolated_cwd, new_pipeline, tmp_path):
    tool, _ = build(new_pipeline, tmp_path)
    call = tool.generate_script("unused").split("pipe_embedding.py")[1].split("pipe_propagate")[0]
    emitted = set(re.findall(r"^\s+(--[a-z-]+)", call, re.M))
    declared = set(re.findall(r'add_argument\("(--[a-z-]+)"', open(SCRIPT).read()))
    assert emitted == declared and len(emitted) == 8


def test_a_configured_container_is_reported_not_silently_ignored(local_config, isolated_cwd, new_pipeline, tmp_path, monkeypatch):
    from biopipelines.config_manager import ConfigManager

    monkeypatch.setattr(ConfigManager, "get_container_executor", lambda self: "apptainer")
    tool, _ = build(new_pipeline, tmp_path)
    plain = tool.generate_script("unused")
    assert "does not support container execution" not in plain and "apptainer" not in plain

    tool.folders = {**tool.folders, "container:Embedding": "/img/embedding.sif", "__container_binds__": ["/w"]}
    contained = tool.generate_script("unused")
    assert "WARNING: a container is configured for Embedding" in contained and "running in environment mode" in contained
    assert "apptainer exec" not in contained


# ── declared outputs ─────────────────────────────────────────────────────────

def test_the_declared_outputs_are_two_streams_and_the_missing_table(local_config, isolated_cwd, new_pipeline, tmp_path):
    tool, _ = build(new_pipeline, tmp_path)
    outputs = tool.get_output_files()
    stream = outputs["embeddings"]
    assert (stream.name, stream.format, list(stream.ids)) == ("embeddings", "npy", ["s1", "s2"])
    assert stream.files[0].endswith("<id>.npy") and stream.map_table == tool.embeddings_map
    assert len(outputs["node_embeddings"]) == 0 and outputs["node_embeddings"].format == "npy"
    assert set(outputs["tables"]) == {"missing"}
    assert outputs["tables"]["missing"].info.columns == ["id", "removed_by", "kind", "cause"]


def test_node_embeddings_and_pairs_declare_their_own_outputs(local_config, isolated_cwd, new_pipeline, tmp_path):
    first, second = pair_columns(tmp_path)
    tool, _ = build(new_pipeline, tmp_path, first=first, second=second, node_embeddings=True)
    outputs = tool.get_output_files()
    nodes = outputs["node_embeddings"]
    assert (nodes.name, nodes.format, list(nodes.ids)) == ("node_embeddings", "npy", ["s1", "s2"])
    assert nodes.files[0].endswith("<id>.npy") and nodes.map_table == tool.node_embeddings_map
    assert set(outputs["tables"]) == {"missing", "similarities"}
    assert outputs["tables"]["similarities"].info.columns == ["id", "first", "second", "cosine_similarity"]


def test_the_stream_maps_and_tables_sit_where_the_pipe_is_told_to_write(local_config, isolated_cwd, new_pipeline, tmp_path):
    first, second = pair_columns(tmp_path)
    tool, _ = build(new_pipeline, tmp_path, first=first, second=second, node_embeddings=True)
    script = tool.generate_script("unused")
    for path in (tool.embeddings_map, tool.node_embeddings_map, tool.similarities_csv, tool.local_missing_csv):
        assert f'"{path}"' in script


def test_the_stream_names_are_registered():
    from biopipelines.contract_enforcement import KNOWN_STREAM_NAMES

    assert {"embeddings", "node_embeddings"} <= KNOWN_STREAM_NAMES


# ── installation and registration ────────────────────────────────────────────

def test_the_install_builds_the_env_and_verifies_the_imports_without_downloading_weights():
    from biopipelines.embedding import Embedding

    script = Embedding._install_script({"biopipelines": "/bp"}, "mamba")
    assert "import torch, torch_scatter, torch_geometric, einops, progres.progres" in script
    assert 'touch "$INSTALL_SUCCESS"' in script and "already installed" in script
    assert "zenodo" not in script.lower() and "wget" not in script and "curl" not in script


def test_force_reinstall_drops_the_already_installed_shortcut():
    from biopipelines.embedding import Embedding

    assert "already installed" not in Embedding._install_script({"biopipelines": "/bp"}, "mamba", force_reinstall=True)


def test_install_can_be_called_in_a_pipeline(local_config, isolated_cwd, new_pipeline):
    from biopipelines.embedding import Embedding
    from biopipelines.mock import Mock

    with new_pipeline("embedding_install"):
        Mock(ids=["s1"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
        Embedding.install()


def test_the_env_is_mapped_in_every_variant_that_runs_it_and_its_files_exist():
    import yaml

    for config in ("config.cluster.yaml", "config.colab.yaml", "config.container.yaml"):
        environments = yaml.safe_load(open(os.path.join(ROOT, config), encoding="utf-8"))["environments"]
        assert environments["Embedding"] == "Embedding", config
    assert "Embedding" not in yaml.safe_load(open(os.path.join(ROOT, "config.daint.yaml"), encoding="utf-8"))["environments"]
    assert os.path.isfile(os.path.join(ROOT, "environments", "Embedding.yaml"))
    pip = open(os.path.join(ROOT, "environments", "Embedding.pip.txt"), encoding="utf-8").read()
    assert re.search(r"^progres$", pip, re.M) and re.search(r"^torch-scatter==", pip, re.M)


# ── against the real progres package and real weights ────────────────────────

def ca_pdb(chains):
    """CA-only PDB text; chains is a list of (chain id, residue count), each a helix."""
    lines, k = [], 1
    for chain, n in chains:
        for i in range(n):
            x, y, z = 2.3 * math.cos(math.radians(100 * i)), 2.3 * math.sin(math.radians(100 * i)), 1.5 * i
            lines.append(f"ATOM  {k:5d}  CA  ALA {chain}{i + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C")
            k += 1
        lines.append("TER")
    return "\n".join(lines + ["END"]) + "\n"


@pytest.fixture(scope="module")
def real_weights():
    path = os.environ.get("EMBEDDING_TEST_WEIGHTS", "")
    if not path or not os.path.isfile(path):
        pytest.skip("set EMBEDDING_TEST_WEIGHTS to a Progres or CIRPIN checkpoint to run the real-model tests")
    for module in ("torch", "torch_scatter", "torch_geometric", "einops", "progres.progres"):
        pytest.importorskip(module)
    return path


@pytest.fixture(scope="module")
def real_embedder(pipe, real_weights):
    from biopipelines.embedding import EMBEDDING_MODELS

    return pipe.build_embedder(EMBEDDING_MODELS["cirpin"], real_weights, "cpu", True)


def test_the_real_checkpoint_loads_and_gives_128_numbers_of_unit_length(real_embedder, tmp_path):
    path = tmp_path / "helix.pdb"
    path.write_text(ca_pdb([("A", 40)]))
    vector, nodes = real_embedder.embed(str(path))
    assert vector.shape == (128,) and vector.dtype == np.float32
    assert np.linalg.norm(vector) == pytest.approx(1.0, abs=1e-5)
    assert nodes.shape == (40, 128) and nodes.dtype == np.float32


def test_the_real_per_node_output_is_what_the_graph_decoder_pools(real_embedder, tmp_path):
    torch = real_embedder.torch
    path = tmp_path / "helix.pdb"
    path.write_text(ca_pdb([("A", 25)]))
    vector, nodes = real_embedder.embed(str(path))
    with torch.no_grad():
        again = torch.nn.functional.normalize(real_embedder.model.graph_dec(torch.tensor(nodes).sum(0, keepdim=True)), dim=1)[0].numpy()
    assert np.allclose(again, vector, atol=1e-5)


def test_the_real_reader_embeds_the_first_chain_only(real_embedder, tmp_path):
    one, two = tmp_path / "one.pdb", tmp_path / "two.pdb"
    one.write_text(ca_pdb([("A", 30)]))
    two.write_text(ca_pdb([("A", 30), ("B", 20)]))
    assert np.allclose(real_embedder.embed(str(one))[0], real_embedder.embed(str(two))[0], atol=1e-6)


def test_the_real_model_refuses_two_residues_with_the_adapter_message(real_embedder, tmp_path):
    path = tmp_path / "tiny.pdb"
    path.write_text(ca_pdb([("A", 2)]))
    with pytest.raises(ValueError, match="at least 3"):
        real_embedder.embed(str(path))
