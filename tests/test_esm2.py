"""ESM2: one forward pass per masked residue gives the pseudo-log-likelihood, one unmasked pass gives the embedding.

The numerical core works on numpy arrays and takes the model as a callable, so the masking enumeration and the aggregation are pinned here against stand-in models. The last block runs the 8M checkpoint and are skipped unless a folder holding it is named in ESM2_TEST_WEIGHTS (ESM2.install(models=["esm2_t6_8M_UR50D"]) writes one; so does pipe_esm2_weights.py).
"""

import importlib.util
import json
import math
import os
import re
import sys

import numpy as np
import pytest

pytest.importorskip("yaml")

ROOT = os.path.join(os.path.dirname(__file__), "..")
SCRIPT = os.path.join(ROOT, "pipe_scripts", "pipe_esm2.py")
WEIGHTS_SCRIPT = os.path.join(ROOT, "pipe_scripts", "pipe_esm2_weights.py")

# vocab.txt of facebook/esm2_t6_8M_UR50D, identical for the five released sizes read here.
VOCAB = ["<cls>", "<pad>", "<eos>", "<unk>", "L", "A", "G", "V", "S", "E", "R", "T", "I", "D", "P", "K", "Q", "N",
         "F", "Y", "M", "H", "W", "C", "X", "B", "U", "Z", "O", ".", "-", "<null_1>", "<mask>"]
MASK_ID = VOCAB.index("<mask>")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def pipe():
    return load("pipe_esm2", SCRIPT)


@pytest.fixture(scope="module")
def weights_pipe():
    return load("pipe_esm2_weights", WEIGHTS_SCRIPT)


def tokens(sequence):
    """[<cls>, residues, <eos>] with the ids of the real vocabulary."""
    return np.array([0] + [VOCAB.index(c) for c in sequence] + [2])


# ── the alphabet and the sequence check ──────────────────────────────────────

def test_the_residue_alphabet_is_the_letters_of_the_esm2_vocabulary(pipe):
    assert set(pipe.RESIDUE_ALPHABET) == {v for v in VOCAB if len(v) == 1 and v.isalpha()}
    assert len(pipe.RESIDUE_ALPHABET) == 25


def test_a_sequence_over_the_alphabet_is_returned_unchanged(pipe):
    assert pipe.validate_sequence("MKTAYIAKQRXBUZO") == "MKTAYIAKQRXBUZO"


@pytest.mark.parametrize("bad,letters", [("MKT1", "1"), ("mkt", "kmt"), ("MK T", " "), ("AC:DE", ":"), ("AC-DE", "-"), ("MKT*", "*")])
def test_letters_outside_the_alphabet_are_named(pipe, bad, letters):
    with pytest.raises(ValueError, match="outside the ESM-2 alphabet") as caught:
        pipe.validate_sequence(bad)
    assert str(caught.value).endswith(": " + letters)


@pytest.mark.parametrize("bad", ["", None, float("nan"), 12])
def test_an_empty_or_non_text_sequence_is_refused(pipe, bad):
    with pytest.raises(ValueError, match="empty or non-text"):
        pipe.validate_sequence(bad)


# ── which positions are masked ───────────────────────────────────────────────

def test_every_residue_is_masked_exactly_once_and_nothing_else(pipe):
    ids = tokens("MKTAYIA")
    seen = []

    def logits_fn(batch):
        seen.append(batch.copy())
        return np.zeros(batch.shape + (len(VOCAB),))

    pipe.masked_log_likelihoods(ids, MASK_ID, logits_fn, batch_size=3)
    rows = np.vstack(seen)
    assert rows.shape == (7, 9)
    for k, row in enumerate(rows):
        masked = np.flatnonzero(row == MASK_ID)
        assert list(masked) == [k + 1]
        kept = np.delete(np.arange(9), k + 1)
        assert np.array_equal(row[kept], ids[kept])
    assert (rows[:, 0] == 0).all() and (rows[:, -1] == 2).all()


@pytest.mark.parametrize("length,batch_size,sizes", [(7, 3, [3, 3, 1]), (6, 3, [3, 3]), (4, 100, [4]), (5, 1, [1] * 5)])
def test_the_masked_copies_go_through_in_batches_of_at_most_batch_size(pipe, length, batch_size, sizes):
    seen = []

    def logits_fn(batch):
        seen.append(len(batch))
        return np.zeros(batch.shape + (len(VOCAB),))

    pipe.masked_log_likelihoods(tokens("A" * length), MASK_ID, logits_fn, batch_size)
    assert seen == sizes


def test_the_input_token_ids_are_not_modified(pipe):
    ids = tokens("MKTAY")
    before = ids.copy()
    pipe.masked_log_likelihoods(ids, MASK_ID, lambda b: np.zeros(b.shape + (len(VOCAB),)), 2)
    assert np.array_equal(ids, before)


def test_a_sequence_with_no_residue_is_refused(pipe):
    with pytest.raises(ValueError, match="at least one residue"):
        pipe.masked_log_likelihoods(np.array([0, 2]), MASK_ID, lambda b: b, 4)


def test_a_logits_array_of_the_wrong_shape_is_refused(pipe):
    with pytest.raises(ValueError, match="shape"):
        pipe.masked_log_likelihoods(tokens("MKT"), MASK_ID, lambda b: np.zeros((len(b), 2, len(VOCAB))), 4)


# ── the pseudo-log-likelihood against stand-in models ────────────────────────

def test_a_uniform_model_gives_minus_log_vocabulary_per_residue(pipe):
    ids = tokens("MKTAYIAKQR")
    per_residue = pipe.masked_log_likelihoods(ids, MASK_ID, lambda b: np.zeros(b.shape + (len(VOCAB),)), 4)
    assert np.allclose(per_residue, -math.log(len(VOCAB)))
    total, mean = pipe.pseudo_log_likelihood(per_residue)
    assert total == pytest.approx(-10 * math.log(len(VOCAB)))
    assert mean == pytest.approx(-math.log(len(VOCAB)))


def oracle(true_ids, confidence):
    """A model that, at a masked position only, puts `confidence` on the true residue and 0 on every other token."""
    def logits_fn(batch):
        out = np.zeros(batch.shape + (len(VOCAB),))
        for b, row in enumerate(batch):
            for t in np.flatnonzero(row == MASK_ID):
                out[b, t, true_ids[t]] = confidence
        return out
    return logits_fn


def test_the_log_probability_is_read_at_the_masked_position_of_each_row(pipe):
    ids = tokens("MKTAYIAKQR")
    expected = confidence = 3.0
    expected -= math.log(math.exp(confidence) + len(VOCAB) - 1)
    per_residue = pipe.masked_log_likelihoods(ids, MASK_ID, oracle(ids, confidence), 3)
    assert np.allclose(per_residue, expected)


def test_each_position_reads_its_own_distribution(pipe):
    ids = tokens("ACDEF")
    table = np.random.default_rng(0).normal(size=(len(ids), len(VOCAB)))

    def logits_fn(batch):
        return np.broadcast_to(table, batch.shape + (len(VOCAB),)).copy()

    got = pipe.masked_log_likelihoods(ids, MASK_ID, logits_fn, 2)
    log_p = table - np.log(np.exp(table).sum(axis=1, keepdims=True))
    assert np.allclose(got, [log_p[i, ids[i]] for i in range(1, 6)])


def test_the_result_does_not_depend_on_the_batch_size(pipe):
    ids = tokens("MKTAYIAKQRQISFVKSHFSRQ")
    table = np.random.default_rng(1).normal(size=(len(ids), len(VOCAB)))
    logits_fn = lambda batch: np.broadcast_to(table, batch.shape + (len(VOCAB),)).copy()  # noqa: E731
    reference = pipe.masked_log_likelihoods(ids, MASK_ID, logits_fn, 1)
    for size in (2, 5, 21, 22, 64):
        assert np.allclose(pipe.masked_log_likelihoods(ids, MASK_ID, logits_fn, size), reference)


def test_the_log_softmax_is_stable_for_large_logits(pipe):
    out = pipe.log_softmax(np.array([[1000.0, 1000.0, -1000.0]]))
    assert np.allclose(out, [[-math.log(2), -math.log(2), -2000.0 - math.log(2)]])
    assert np.isfinite(out).all()


def test_sum_and_mean_of_the_per_residue_values(pipe):
    assert pipe.pseudo_log_likelihood(np.array([-1.0, -2.0, -6.0])) == (-9.0, -3.0)


# ── the mean-pooled embedding ────────────────────────────────────────────────

def test_the_embedding_averages_the_residue_positions_only(pipe):
    hidden = np.zeros((6, 3))
    hidden[0] = 1000.0
    hidden[-1] = -1000.0
    hidden[1:-1] = np.arange(12).reshape(4, 3)
    vector = pipe.mean_embedding(hidden)
    assert np.allclose(vector, [4.5, 5.5, 6.5])
    assert vector.dtype == np.float32 and vector.shape == (3,)


@pytest.mark.parametrize("shape", [(2, 4), (5,), (4, 3, 2)])
def test_a_hidden_state_without_a_residue_row_is_refused(pipe, shape):
    with pytest.raises(ValueError, match="hidden state"):
        pipe.mean_embedding(np.zeros(shape))


# ── the device and the weights folder ────────────────────────────────────────

@pytest.mark.parametrize("requested,cuda,expected", [("auto", True, "cuda"), ("auto", False, "cpu"), ("cpu", True, "cpu"), ("cpu", False, "cpu"), ("cuda", True, "cuda")])
def test_the_device_is_resolved(pipe, requested, cuda, expected):
    assert pipe.resolve_device(requested, cuda) == expected


def test_cuda_without_a_gpu_fails_instead_of_running_on_cpu(pipe):
    with pytest.raises(RuntimeError, match="no CUDA device"):
        pipe.resolve_device("cuda", False)


def test_an_unknown_device_is_refused(pipe):
    with pytest.raises(ValueError, match="device"):
        pipe.resolve_device("tpu", True)


def make_weights(folder, with_vocab=True, weights="model.safetensors"):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "config.json").write_text("{}")
    if with_vocab:
        (folder / "vocab.txt").write_text("\n".join(VOCAB))
    if weights:
        (folder / weights).write_bytes(b"")
    return str(folder)


def test_a_complete_weights_folder_passes_the_check(pipe, tmp_path):
    pipe.check_weights_folder(make_weights(tmp_path / "m"))
    pipe.check_weights_folder(make_weights(tmp_path / "sharded", weights="pytorch_model-00001-of-00002.bin"))


@pytest.mark.parametrize("kwargs,missing", [({"with_vocab": False}, "vocab.txt"), ({"weights": None}, "weights file")])
def test_an_incomplete_weights_folder_names_what_is_missing_and_how_to_install(pipe, tmp_path, kwargs, missing):
    with pytest.raises(FileNotFoundError, match=missing) as caught:
        pipe.check_weights_folder(make_weights(tmp_path / "m", **kwargs))
    assert "ESM2.install()" in str(caught.value)


def test_a_missing_weights_folder_is_refused(pipe, tmp_path):
    with pytest.raises(FileNotFoundError, match="does not exist"):
        pipe.check_weights_folder(str(tmp_path / "absent"))


# ── what the weights script downloads ────────────────────────────────────────

FILES_8M = [".gitattributes", "README.md", "config.json", "model.safetensors", "pytorch_model.bin",
            "special_tokens_map.json", "tf_model.h5", "tokenizer_config.json", "vocab.txt"]
FILES_3B = [".gitattributes", "README.md", "config.json", "pytorch_model-00001-of-00002.bin",
            "pytorch_model-00002-of-00002.bin", "pytorch_model.bin.index.json", "special_tokens_map.json",
            "tf_model-00001-of-00002.h5", "tf_model-00002-of-00002.h5", "tf_model.h5.index.json",
            "tokenizer_config.json", "vocab.txt"]


def test_a_repository_with_safetensors_downloads_them_and_not_the_pickles(weights_pipe):
    assert weights_pipe.select_files(FILES_8M) == [
        "config.json", "vocab.txt", "special_tokens_map.json", "tokenizer_config.json", "model.safetensors"]


def test_a_repository_with_only_sharded_pickles_downloads_the_shards_and_their_index(weights_pipe):
    assert weights_pipe.select_files(FILES_3B) == [
        "config.json", "vocab.txt", "special_tokens_map.json", "tokenizer_config.json",
        "pytorch_model-00001-of-00002.bin", "pytorch_model-00002-of-00002.bin", "pytorch_model.bin.index.json"]


def test_sharded_safetensors_bring_their_index(weights_pipe):
    files = ["config.json", "vocab.txt", "model-00001-of-00002.safetensors", "model-00002-of-00002.safetensors",
             "model.safetensors.index.json", "pytorch_model.bin"]
    assert weights_pipe.select_files(files)[-3:] == [
        "model-00001-of-00002.safetensors", "model-00002-of-00002.safetensors", "model.safetensors.index.json"]


@pytest.mark.parametrize("files,message", [
    (["config.json", "model.safetensors"], "vocab.txt"),
    (["vocab.txt", "model.safetensors"], "config.json"),
    (["config.json", "vocab.txt", "tf_model.h5"], "neither safetensors nor"),
])
def test_a_repository_missing_what_the_run_needs_is_refused(weights_pipe, files, message):
    with pytest.raises(ValueError, match=message):
        weights_pipe.select_files(files)


def test_a_folder_that_is_already_pinned_is_not_downloaded_again(weights_pipe, tmp_path, capsys):
    folder = tmp_path / "esm2_t6_8M_UR50D"
    folder.mkdir()
    (folder / "REVISION").write_text("facebook/esm2_t6_8M_UR50D@abc123\n")
    weights_pipe.download("esm2_t6_8M_UR50D", str(tmp_path))
    printed = capsys.readouterr().out
    assert "already pinned" in printed and "abc123" in printed


def test_every_model_the_wrapper_offers_is_a_choice_of_the_weights_script(weights_pipe):
    from biopipelines.esm2 import ESM2_MODELS
    assert weights_pipe.ESM2_MODELS == ESM2_MODELS


# ── the step loop, against a stand-in model ──────────────────────────────────

class FakeModel:
    """Uniform logits over the vocabulary; hidden rows carry the token id, so the pooled vector is the mean token id."""

    mask_id = MASK_ID
    device = "cpu"
    created = []
    batches = []

    def __init__(self, folder, device):
        FakeModel.created.append((folder, device))

    def encode(self, sequence):
        return tokens(sequence)

    def logits(self, batch):
        FakeModel.batches.append(len(batch))
        return np.zeros(batch.shape + (len(VOCAB),))

    def hidden(self, token_ids):
        return np.tile(np.asarray(token_ids, dtype=float)[:, None], (1, 4))


@pytest.fixture
def fake_model(pipe, monkeypatch):
    FakeModel.created, FakeModel.batches = [], []
    monkeypatch.setattr(pipe, "Esm2Model", FakeModel)
    return FakeModel


def write_inputs(tmp_path, sequences, compute_pll=True, compute_embeddings=True, batch_size=4, device="auto", weights=True):
    import pandas as pd

    ids = list(sequences)
    map_csv = tmp_path / "sequences_map.csv"
    pd.DataFrame({"id": ids, "sequence": [sequences[i] for i in ids]}).to_csv(map_csv, index=False)
    stream = tmp_path / "sequences.json"
    stream.write_text(json.dumps({"name": "sequences", "ids": ids, "files": [], "map_table": str(map_csv), "format": "csv"}))
    config = tmp_path / "config.yaml"
    config.write_text(json.dumps({"model": "esm2_t6_8M_UR50D", "compute_pll": compute_pll,
                                  "compute_embeddings": compute_embeddings, "batch_size": batch_size, "device": device}))
    out = tmp_path / "out"
    return [
        "pipe_esm2.py", "--config-yaml", str(config), "--sequences-json", str(stream),
        "--weights-dir", make_weights(tmp_path / "weights") if weights else str(tmp_path / "absent"),
        "--embeddings-dir", str(out / "embeddings"), "--embeddings-map-csv", str(out / "embeddings" / "embeddings_map.csv"),
        "--scores-csv", str(out / "tables" / "scores.csv"), "--local-missing-csv", str(out / "_execution" / "local_missing.csv"),
    ], out


def run_main(pipe, monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", argv)
    pipe.main()


def test_one_scores_row_per_id_with_a_failed_id_as_nan_and_a_missing_entry(pipe, fake_model, monkeypatch, tmp_path):
    import pandas as pd

    argv, out = write_inputs(tmp_path, {"a": "MKTAYIA", "bad": "MK1", "c": "ACDE"})
    run_main(pipe, monkeypatch, argv)

    scores = pd.read_csv(out / "tables" / "scores.csv")
    assert list(scores.columns) == pipe.SCORE_COLUMNS == ["id", "length", "pll", "pll_per_residue"]
    assert list(scores["id"]) == ["a", "bad", "c"]
    assert list(scores["length"].astype("Int64")) == [7, pd.NA, 4]
    assert scores.loc[0, "pll"] == pytest.approx(-7 * math.log(len(VOCAB)))
    assert scores.loc[2, "pll_per_residue"] == pytest.approx(-math.log(len(VOCAB)))
    assert scores.loc[1, ["pll", "pll_per_residue"]].isna().all()

    missing = pd.read_csv(out / "_execution" / "local_missing.csv")
    assert list(missing.columns) == ["id", "removed_by", "kind", "cause"]
    assert list(missing["id"]) == ["bad"] and missing.loc[0, "kind"] == "failure" and "1" in missing.loc[0, "cause"]


def test_the_embedding_files_and_their_map_cover_only_the_ids_that_succeeded(pipe, fake_model, monkeypatch, tmp_path):
    import pandas as pd

    argv, out = write_inputs(tmp_path, {"a": "MKTAYIA", "bad": "MK1", "c": "ACDE"})
    run_main(pipe, monkeypatch, argv)

    embedding_map = pd.read_csv(out / "embeddings" / "embeddings_map.csv")
    assert list(embedding_map["id"]) == ["a", "c"]
    assert sorted(os.listdir(out / "embeddings")) == ["a.npy", "c.npy", "embeddings_map.csv"]
    vector = np.load(embedding_map.loc[1, "file"])
    assert vector.dtype == np.float32 and vector.shape == (4,)
    assert np.allclose(vector, np.mean([VOCAB.index(c) for c in "ACDE"]))


def test_batch_size_and_device_reach_the_model(pipe, fake_model, monkeypatch, tmp_path):
    argv, _ = write_inputs(tmp_path, {"a": "MKTAYIA"}, batch_size=3, device="cpu")
    run_main(pipe, monkeypatch, argv)
    assert fake_model.batches == [3, 3, 1]
    assert fake_model.created == [(str(tmp_path / "weights"), "cpu")]


def test_without_pll_there_is_no_scores_table_and_no_forward_pass_per_residue(pipe, fake_model, monkeypatch, tmp_path):
    argv, out = write_inputs(tmp_path, {"a": "MKTAYIA"}, compute_pll=False)
    run_main(pipe, monkeypatch, argv)
    assert fake_model.batches == []
    assert not (out / "tables" / "scores.csv").exists()
    assert (out / "embeddings" / "a.npy").exists()


def test_without_embeddings_there_is_no_stream_map_and_no_file(pipe, fake_model, monkeypatch, tmp_path):
    argv, out = write_inputs(tmp_path, {"a": "MKTAYIA"}, compute_embeddings=False)
    run_main(pipe, monkeypatch, argv)
    assert (out / "tables" / "scores.csv").exists()
    assert not (out / "embeddings" / "embeddings_map.csv").exists()
    assert os.listdir(out / "embeddings") == []


def test_every_id_failing_is_a_step_failure_but_the_tables_are_still_written(pipe, fake_model, monkeypatch, tmp_path):
    import pandas as pd

    argv, out = write_inputs(tmp_path, {"a": "M1", "b": ""})
    with pytest.raises(SystemExit) as caught:
        run_main(pipe, monkeypatch, argv)
    assert caught.value.code == 1
    assert list(pd.read_csv(out / "tables" / "scores.csv")["id"]) == ["a", "b"]
    assert list(pd.read_csv(out / "_execution" / "local_missing.csv")["id"]) == ["a", "b"]


def test_a_stream_with_no_id_is_a_step_failure(pipe, fake_model, monkeypatch, tmp_path):
    argv, _ = write_inputs(tmp_path, {})
    with pytest.raises(SystemExit) as caught:
        run_main(pipe, monkeypatch, argv)
    assert caught.value.code == 1


def test_a_weights_folder_that_is_not_installed_stops_before_any_model_is_built(pipe, fake_model, monkeypatch, tmp_path):
    argv, _ = write_inputs(tmp_path, {"a": "MKT"}, weights=False)
    with pytest.raises(FileNotFoundError, match="ESM2.install"):
        run_main(pipe, monkeypatch, argv)
    assert fake_model.created == []


def test_the_pinned_revision_is_printed_when_recorded(pipe, fake_model, monkeypatch, tmp_path, capsys):
    argv, _ = write_inputs(tmp_path, {"a": "MKT"})
    (tmp_path / "weights" / "REVISION").write_text("facebook/esm2_t6_8M_UR50D@c731040\n")
    run_main(pipe, monkeypatch, argv)
    assert "facebook/esm2_t6_8M_UR50D@c731040" in capsys.readouterr().out


# ── the wrapper ──────────────────────────────────────────────────────────────

def build(new_pipeline, **overrides):
    from biopipelines.esm2 import ESM2
    from biopipelines.mock import Mock

    pipeline = new_pipeline("esm2")
    with pipeline:
        source = Mock(ids=["s1", "s2"], streams={"sequences": {"format": "csv", "values": ["MKTAYIA", "ACDEFG"]}})
        ESM2(sequences=source, **overrides)
        script_path = pipeline.save()
    return pipeline.tools[-1], script_path


def test_defaults(local_config, isolated_cwd, new_pipeline):
    tool, _ = build(new_pipeline)
    assert tool.model == "esm2_t33_650M_UR50D"
    assert (tool.compute_pll, tool.compute_embeddings, tool.batch_size, tool.device) == (True, True, 16, "auto")
    assert tool.sequences_stream.name == "sequences" and list(tool.sequences_stream.ids) == ["s1", "s2"]


def test_a_bare_sequence_string_is_refused(local_config, isolated_cwd, new_pipeline):
    from biopipelines.esm2 import ESM2
    from biopipelines.mock import Mock

    with new_pipeline("esm2_bare"):
        Mock(ids=["s1"], streams={"sequences": {"format": "csv", "values": ["MKT"]}})
        with pytest.raises(ValueError, match="sequences must be"):
            ESM2(sequences="MKTAYIA")


def test_every_offered_model_is_accepted_and_the_15b_one_is_not(local_config, isolated_cwd, new_pipeline):
    from biopipelines.esm2 import ESM2_MODELS

    assert ESM2_MODELS == ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D", "esm2_t30_150M_UR50D",
                           "esm2_t33_650M_UR50D", "esm2_t36_3B_UR50D")
    for model in ESM2_MODELS:
        assert build(new_pipeline, model=model)[0].model == model
    with pytest.raises(ValueError, match="model must be one of"):
        build(new_pipeline, model="esm2_t48_15B_UR50D")


@pytest.mark.parametrize("overrides,message", [
    ({"model": "facebook/esm2_t33_650M_UR50D"}, "model must be one of"),
    ({"model": ""}, "model must be one of"),
    ({"compute_pll": 1}, "compute_pll must be True or False"),
    ({"compute_embeddings": "yes"}, "compute_embeddings must be True or False"),
    ({"compute_pll": False, "compute_embeddings": False}, "nothing to compute"),
    ({"batch_size": 0}, "batch_size"),
    ({"batch_size": -4}, "batch_size"),
    ({"batch_size": 2.5}, "batch_size"),
    ({"batch_size": True}, "batch_size"),
    ({"device": "gpu"}, "device must be one of"),
])
def test_bad_parameters_fail_at_configuration(local_config, isolated_cwd, new_pipeline, overrides, message):
    with pytest.raises(ValueError, match=message):
        build(new_pipeline, **overrides)


def test_the_generated_script_runs_the_helper_on_the_pinned_model_folder(local_config, isolated_cwd, new_pipeline):
    import yaml

    tool, _ = build(new_pipeline, model="esm2_t12_35M_UR50D", batch_size=8, device="cpu", compute_embeddings=False)
    script = tool.generate_script("unused")
    assert "pipe_esm2.py" in script
    folder = os.path.join(tool.folders["ESM2"], "esm2_t12_35M_UR50D")
    assert f'--weights-dir "{folder}"' in script
    assert yaml.safe_load(open(tool.config_yaml)) == {
        "model": "esm2_t12_35M_UR50D", "compute_pll": True, "compute_embeddings": False, "batch_size": 8, "device": "cpu"}
    assert os.path.exists(tool.sequences_json)
    assert "pipe_propagate_missing.py" in script and tool.local_missing_csv in script


def test_every_flag_the_wrapper_emits_is_declared_by_the_pipe(local_config, isolated_cwd, new_pipeline):
    tool, _ = build(new_pipeline)
    call = tool.generate_script("unused").split("pipe_esm2.py")[1].split("pipe_propagate")[0]
    emitted = set(re.findall(r"^\s+(--[a-z-]+)", call, re.M))
    declared = set(re.findall(r'add_argument\("(--[a-z-]+)"', open(SCRIPT).read()))
    assert emitted == declared


def test_no_weights_folder_in_the_config_is_named(local_config, isolated_cwd, new_pipeline):
    tool, _ = build(new_pipeline)
    tool.folders = {}
    with pytest.raises(ValueError, match="repositories"):
        tool.generate_script("unused")


def test_the_declared_outputs_are_the_embeddings_stream_and_two_tables(local_config, isolated_cwd, new_pipeline, pipe):
    tool, _ = build(new_pipeline)
    outputs = tool.get_output_files()
    stream = outputs["embeddings"]
    assert (stream.name, stream.format, list(stream.ids)) == ("embeddings", "npy", ["s1", "s2"])
    assert stream.files[0].endswith("<id>.npy") and stream.map_table == tool.embeddings_map
    assert set(outputs["tables"]) == {"scores", "missing"}
    assert outputs["tables"]["scores"].info.columns == pipe.SCORE_COLUMNS
    assert outputs["tables"]["missing"].info.columns == ["id", "removed_by", "kind", "cause"]


def test_without_pll_the_scores_table_is_not_declared(local_config, isolated_cwd, new_pipeline):
    outputs = build(new_pipeline, compute_pll=False)[0].get_output_files()
    assert set(outputs["tables"]) == {"missing"} and len(outputs["embeddings"]) == 2


def test_without_embeddings_the_stream_is_an_empty_placeholder(local_config, isolated_cwd, new_pipeline):
    outputs = build(new_pipeline, compute_embeddings=False)[0].get_output_files()
    assert len(outputs["embeddings"]) == 0 and "scores" in outputs["tables"]


def test_embeddings_is_a_registered_stream_name():
    from biopipelines.contract_enforcement import KNOWN_STREAM_NAMES
    assert "embeddings" in KNOWN_STREAM_NAMES


# ── installation ─────────────────────────────────────────────────────────────

def test_the_install_pins_the_default_model_and_verifies_env_and_weights():
    from biopipelines.esm2 import ESM2

    script = ESM2._install_script({"ESM2": "/data/esm2", "biopipelines": "/bp"}, "mamba")
    assert '/bp/pipe_scripts/pipe_esm2_weights.py" --dest "/data/esm2" --models esm2_t33_650M_UR50D' in script
    assert '[ -f "/data/esm2/esm2_t33_650M_UR50D/REVISION" ]' in script
    assert "import torch, transformers, huggingface_hub" in script and 'touch "$INSTALL_SUCCESS"' in script
    assert "already installed" in script


def test_several_models_are_downloaded_and_each_is_verified():
    from biopipelines.esm2 import ESM2

    script = ESM2._install_script({"ESM2": "/data/esm2", "biopipelines": "/bp"}, "mamba",
                                  models=["esm2_t6_8M_UR50D", "esm2_t36_3B_UR50D"])
    assert "--models esm2_t6_8M_UR50D esm2_t36_3B_UR50D" in script
    assert script.count("/REVISION\" ]") == 4


def test_a_single_model_name_is_accepted_for_models():
    from biopipelines.esm2 import ESM2

    assert "--models esm2_t6_8M_UR50D\n" in ESM2._install_script({"ESM2": "/d", "biopipelines": "/bp"}, "mamba", models="esm2_t6_8M_UR50D")


@pytest.mark.parametrize("models", [["esm2_t48_15B_UR50D"], [], ["esm2_t6_8M_UR50D", "nonsense"]])
def test_unknown_or_no_models_are_refused_at_install(models):
    from biopipelines.esm2 import ESM2

    with pytest.raises(ValueError, match="models must be"):
        ESM2._install_script({"ESM2": "/d", "biopipelines": "/bp"}, "mamba", models=models)


def test_force_reinstall_drops_the_already_installed_shortcut():
    from biopipelines.esm2 import ESM2

    assert "already installed" not in ESM2._install_script({"ESM2": "/d", "biopipelines": "/bp"}, "mamba", force_reinstall=True)


def test_the_install_names_the_missing_config_key():
    from biopipelines.esm2 import ESM2

    with pytest.raises(ValueError, match="repositories.ESM2"):
        ESM2._install_script({}, "mamba")


def test_install_with_a_models_choice_can_be_called_in_a_pipeline(local_config, isolated_cwd, new_pipeline):
    from biopipelines.esm2 import ESM2
    from biopipelines.mock import Mock

    with new_pipeline("esm2_install"):
        Mock(ids=["s1"], streams={"sequences": {"format": "csv", "values": ["MKT"]}})
        ESM2.install(models=["esm2_t6_8M_UR50D"])


# ── against the real 8M checkpoint ───────────────────────────────────────────

SEQUENCE = "MKTAYIAKQRQISFVKSHFSRQ"


@pytest.fixture(scope="module")
def real_folder():
    folder = os.environ.get("ESM2_TEST_WEIGHTS", "")
    if not folder or not os.path.isfile(os.path.join(folder, "config.json")):
        pytest.skip("set ESM2_TEST_WEIGHTS to a folder holding facebook/esm2_t6_8M_UR50D to run the real-model tests")
    pytest.importorskip("torch")
    pytest.importorskip("transformers")
    return folder


@pytest.fixture(scope="module")
def real_model(pipe, real_folder):
    return pipe.Esm2Model(real_folder, "cpu")


def test_the_real_vocabulary_is_the_one_the_alphabet_and_the_tests_assume(real_folder):
    assert [line.strip() for line in open(os.path.join(real_folder, "vocab.txt")) if line.strip()] == VOCAB


def test_the_real_tokenizer_gives_cls_residues_eos_and_the_mask_id(real_model):
    assert list(real_model.encode(SEQUENCE)) == list(tokens(SEQUENCE))
    assert real_model.mask_id == MASK_ID


def test_the_real_pll_equals_a_naive_loop_of_single_masked_forward_passes(pipe, real_model):
    import torch

    ids = real_model.encode(SEQUENCE)
    got = pipe.masked_log_likelihoods(ids, real_model.mask_id, real_model.logits, 5)
    naive = []
    with torch.no_grad():
        for i in range(1, len(ids) - 1):
            x = torch.tensor(ids)[None].clone()
            x[0, i] = real_model.mask_id
            log_p = torch.log_softmax(real_model.model(input_ids=x).logits[0, i].double(), dim=-1)
            naive.append(log_p[int(ids[i])].item())
    assert np.allclose(got, naive, atol=1e-4)
    assert (got < 0).all() and pipe.pseudo_log_likelihood(got)[0] == pytest.approx(sum(naive), abs=1e-3)
    assert np.allclose(pipe.masked_log_likelihoods(ids, real_model.mask_id, real_model.logits, 1), got, atol=1e-4)


def test_the_real_embedding_is_the_mean_of_the_encoders_last_hidden_state_over_the_residues(pipe, real_model):
    import torch

    ids = real_model.encode(SEQUENCE)
    with torch.no_grad():
        last = real_model.model.esm(input_ids=torch.tensor(ids)[None]).last_hidden_state[0].numpy()
    vector = pipe.mean_embedding(real_model.hidden(ids))
    assert vector.shape == (real_model.model.config.hidden_size,) == (320,)
    assert np.allclose(vector, last[1:-1].mean(axis=0), atol=1e-5)


def test_the_real_model_runs_the_whole_step(pipe, real_folder, monkeypatch, tmp_path):
    import pandas as pd

    argv, out = write_inputs(tmp_path, {"a": SEQUENCE, "bad": "MK1", "b": "ACDEFGHIK"}, device="cpu")
    argv[argv.index("--weights-dir") + 1] = real_folder
    run_main(pipe, monkeypatch, argv)

    scores = pd.read_csv(out / "tables" / "scores.csv")
    assert list(scores["id"]) == ["a", "bad", "b"] and list(scores["length"].astype("Int64")) == [22, pd.NA, 9]
    assert (scores.loc[[0, 2], "pll"] < 0).all() and pd.isna(scores.loc[1, "pll"])
    assert np.allclose(scores["pll"][[0, 2]] / scores["length"][[0, 2]], scores["pll_per_residue"][[0, 2]])
    assert np.load(out / "embeddings" / "a.npy").shape == (320,)
    assert list(pd.read_csv(out / "_execution" / "local_missing.csv")["id"]) == ["bad"]
