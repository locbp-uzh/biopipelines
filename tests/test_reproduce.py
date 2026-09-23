"""Reproduction has to be honest about what it cannot promise.

The manifest says what a run was configured to be. It does not say which environment executed it, because the environment lives on the node and the manifest is written before the job leaves the laptop. So the run itself records a digest of its environments, and the comparison has to notice that two runs with an identical plan can still have executed against different software.

That is the trap these tests exist for: the manifest hash is sealed at `save()`, so a comparison that short-circuits on an equal hash would report "no difference" for exactly the case the provenance work is meant to catch — the same pipeline, rerun after someone updated a conda environment.

The rest holds `bp_reproduce` to saying what it does not know. A probe that could not reach the target, a run whose script was never preserved, an environment that was never digested: each has to be named, because "no differences found" and "nothing was checked" look identical to whoever reads the report.
"""

import json
import os

import pytest

from biopipelines import manifest, reproduce


class FakeFS:
    """Enough of the `fs` protocol for these readers, over a real temp tree."""

    def __init__(self, root):
        self.root = str(root)

    def join(self, *parts):
        return "/".join(str(p).rstrip("/") for p in parts if str(p) != "")

    def listdir(self, path):
        if not os.path.isdir(path):
            return []
        return sorted((n, os.path.isdir(os.path.join(path, n))) for n in os.listdir(path))

    def read_text(self, path):
        with open(path, encoding="utf-8") as handle:
            return handle.read()

    def is_file(self, path):
        return os.path.isfile(path)


def make_run(tmp_path, name, *, digests=None, script=True, record=None):
    job = tmp_path / name
    (job / "RunTime").mkdir(parents=True)
    body = record or {"schema": 1, "hash": "sha256:same", "created": "2026-01-01T00:00:00",
                      "biopipelines": "1.4.1", "commit": "abc1234", "variant": "cluster",
                      "scheduler": "slurm", "project": "P", "job": name,
                      "tools": [{"order": 1, "tool": "ProteinMPNN", "tool_version": "2.6",
                                 "environments": ["proteinmpnn"], "container_image": "",
                                 "parameters": {"passed": {}, "resolved": {"sampling_temp": 0.1}}}]}
    (job / manifest.FILENAME).write_text(json.dumps(body), encoding="utf-8")
    if script:
        (job / "RunTime" / "campaign.py").write_text("# pipeline\n", encoding="utf-8")
    if digests:
        environments = job / "environments"
        environments.mkdir()
        for env, digest in digests.items():
            (environments / f"{env}.sha256").write_text(digest + "\n", encoding="utf-8")
    return str(job)


class TestTheEnvironmentIsPartOfTheIdentity:

    def test_digests_are_read_back(self, tmp_path):
        job = make_run(tmp_path, "a", digests={"proteinmpnn": "d" * 64})
        assert manifest.environments(job, fs=FakeFS(tmp_path)) == {"proteinmpnn": "d" * 64}

    def test_an_unavailable_digest_is_not_recorded_as_a_value(self, tmp_path):
        job = make_run(tmp_path, "b", digests={"proteinmpnn": "unavailable"})
        assert manifest.environments(job, fs=FakeFS(tmp_path)) == {}, (
            "'unavailable' is the node saying it had no sha256sum, not an environment identity")

    def test_read_merges_them_into_the_manifest(self, tmp_path):
        job = make_run(tmp_path, "c", digests={"proteinmpnn": "e" * 64})
        record = manifest.read(job, fs=FakeFS(tmp_path))
        assert record["environments_resolved"] == {"proteinmpnn": "e" * 64}

    def test_a_changed_environment_is_seen_although_the_hash_matches(self, tmp_path):
        """The regression this file exists for."""
        fs = FakeFS(tmp_path)
        before = manifest.read(make_run(tmp_path, "before", digests={"proteinmpnn": "1" * 64}), fs=fs)
        after = manifest.read(make_run(tmp_path, "after", digests={"proteinmpnn": "2" * 64}), fs=fs)
        assert before["hash"] == after["hash"], "the plans are identical; only the software moved"

        differences = manifest.compare(before, after)
        assert differences, (
            "an equal hash hid a changed environment — the hash is sealed at save(), before the "
            "environment that runs the job is even known")
        assert any("proteinmpnn" in line for line in differences)

    def test_an_unrecorded_environment_is_said_to_be_unknown(self, tmp_path):
        fs = FakeFS(tmp_path)
        known = manifest.read(make_run(tmp_path, "known", digests={"proteinmpnn": "3" * 64}), fs=fs)
        blank = manifest.read(make_run(tmp_path, "blank"), fs=fs)
        lines = manifest.environment_differences(known, blank)
        assert lines and "cannot be compared" in lines[0], (
            "silence would claim the environments match, which an absent record cannot support")

    def test_two_runs_with_no_digests_at_all_report_nothing(self, tmp_path):
        fs = FakeFS(tmp_path)
        one = manifest.read(make_run(tmp_path, "one"), fs=fs)
        two = manifest.read(make_run(tmp_path, "two"), fs=fs)
        assert manifest.environment_differences(one, two) == []


class TestThePlan:

    def test_it_finds_the_preserved_script(self, tmp_path):
        job = make_run(tmp_path, "withscript")
        assert reproduce.preserved_script(job, FakeFS(tmp_path)).endswith("campaign.py")

    def test_a_run_with_no_script_says_there_is_nothing_to_rerun(self, tmp_path):
        job = make_run(tmp_path, "noscript", script=False)
        found = reproduce.plan(job, FakeFS(tmp_path))
        assert found["script"] is None
        assert "NOT PRESERVED" in reproduce.summarize(found)

    def test_a_run_with_no_manifest_says_so_rather_than_guessing(self, tmp_path):
        bare = tmp_path / "old"
        bare.mkdir()
        found = reproduce.plan(str(bare), FakeFS(tmp_path))
        assert "predates the manifest" in found["error"]
        assert reproduce.summarize(found) == found["error"]

    def test_an_unprobed_target_is_not_reported_as_matching(self, tmp_path):
        found = reproduce.plan(make_run(tmp_path, "unprobed"), FakeFS(tmp_path))
        text = reproduce.summarize(found, host="daint")
        assert "was not probed" in text
        assert "matches the record" not in text, (
            "an unprobed target reported as matching is the failure mode this whole tool is for")

    def test_the_summary_ends_with_the_command_and_its_cost(self, tmp_path):
        found = reproduce.plan(make_run(tmp_path, "ready"), FakeFS(tmp_path))
        text = reproduce.summarize(found)
        assert "bp_submit" in text and "campaign.py" in text
        assert "SPENDS COMPUTE" in text


class TestDifferences:

    RECORD = {"biopipelines": "1.4.1", "commit": "abc1234", "variant": "cluster",
              "scheduler": "slurm", "env_manager": "mamba"}

    def test_a_field_the_target_did_not_report_is_not_a_mismatch(self):
        assert reproduce.differences(self.RECORD, {"biopipelines": "1.4.1"}) == []

    def test_a_moved_commit_is_named(self):
        lines = reproduce.differences(self.RECORD, {"commit": "def5678"})
        assert len(lines) == 1 and "abc1234" in lines[0] and "def5678" in lines[0]

    def test_a_null_from_the_probe_is_not_a_mismatch(self):
        assert reproduce.differences(self.RECORD, {"commit": None}) == [], (
            "a probe that could not read the commit did not find a different one")

    def test_a_second_cluster_reports_every_field_that_moved(self):
        lines = reproduce.differences(self.RECORD, {"variant": "daint", "scheduler": "slurm",
                                                    "commit": "abc1234", "biopipelines": "1.5.0"})
        assert len(lines) == 2
        assert any("config variant" in line for line in lines)
        assert any("BioPipelines version" in line for line in lines)


class TestTheScriptCapturesTheEnvironment:

    def build(self, job, isolated_cwd, debug=False):
        """A pipeline on the cluster variant, whose tools actually name conda environments.

        The `local` fixture config gives its tools no environments, so nothing is there to
        capture and the block is correctly empty — which would make this test pass on a site
        where the feature does nothing.
        """
        from biopipelines import Pipeline, Resources
        from biopipelines.mock import Mock
        from biopipelines.protein_mpnn import ProteinMPNN

        pipeline = Pipeline(project="TestSuite", job=job, description="d", on_the_fly=False,
                            local_output=True, config="cluster", debug=debug)
        with pipeline:
            Resources(partition="gpu", time="01:00:00")
            source = Mock(ids=["d0"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
            ProteinMPNN(structures=source.streams.structures, num_sequences=1)
            script = pipeline.save()
        return pipeline, script

    def test_the_capture_block_is_emitted_without_debug(self, isolated_cwd):
        pipeline, script = self.build("envcap", isolated_cwd)
        assert pipeline.debug is False
        assert pipeline._iter_pipeline_envs(), "this site must name environments or nothing is proven"

        body = open(script, encoding="utf-8").read()
        assert "BioPipelines environment capture" in body, (
            "the environment is recorded only under debug, which is off when it matters")
        assert ".sha256" in body, "the capture writes no digest, so nothing can be compared"

    def test_the_export_is_redacted_before_it_is_kept(self, isolated_cwd):
        """Found by running a real campaign: the export carried a GitLab token in plaintext.

        `pip freeze` renders an editable install as the URL it was cloned from, and the cluster
        clone carries credentials in that URL. The capture wrote it into the group-readable
        output root. Nothing unredacted may survive the block.
        """
        _pipeline, script = self.build("redact", isolated_cwd)
        body = open(script, encoding="utf-8").read()
        capture = body.split("BioPipelines environment capture")[1].split("end environment capture")[0]

        assert "bp_redact" in capture, "the export is written without passing through redaction"
        assert "glpat" in capture, "the token patterns are not in the emitted filter"
        assert "umask 077; mktemp" in capture and 'rm -f "$bp_raw"' in capture, (
            "the unredacted capture must be written to a private scratch file and removed")
        for line in capture.splitlines():
            if "bp_digest" in line and ".sha256" in line:
                assert "$bp_raw" not in line, "the raw, unredacted file must never be digested"

    def test_redaction_is_defined_in_the_debug_block_too(self, isolated_cwd):
        """The debug capture runs its own `pip freeze`, and may run when the other block did not."""
        _pipeline, script = self.build("redact_debug", isolated_cwd, debug=True)
        body = open(script, encoding="utf-8").read()
        debug = body.split("BioPipelines debug capture")[1].split("end debug capture")[0]
        assert "bp_redact() {" in debug, (
            "a pipeline naming no environment skips the other block entirely, leaving bp_redact "
            "undefined here")
        for line in debug.splitlines():
            if "pip freeze" in line:
                assert "bp_redact" in line, f"unredacted pip freeze: {line.strip()}"

    def test_a_failed_export_is_not_digested(self, isolated_cwd):
        """Found on s3it: the block ran, the export failed, and it hashed the error message.

        That digest is stable and shaped exactly like an environment identity, so two runs that
        both failed to export would compare as the same environment — the false match the whole
        comparison exists to prevent.
        """
        _pipeline, script = self.build("failed_export", isolated_cwd)
        body = open(script, encoding="utf-8").read()
        capture = body.split("BioPipelines environment capture")[1].split("end environment capture")[0]
        assert "rm -f" in capture, "a failed export must remove any digest, not leave a stale one"
        for line in capture.splitlines():
            if "bp_digest" in line and ".sha256" in line:
                assert line.startswith("  "), (
                    "the digest is unconditional; it must sit inside the success branch")

    def test_the_export_keeps_build_strings(self, isolated_cwd):
        """`--no-builds` is right for a portable env file and wrong for identity."""
        _pipeline, script = self.build("builds", isolated_cwd)
        body = open(script, encoding="utf-8").read()
        capture = body.split("BioPipelines environment capture")[1].split("end environment capture")[0]
        assert "env export" in capture and "--no-builds" not in capture, (
            "two different builds of the same versions would digest identically")

    def test_debug_capture_still_generates_on_a_conda_site(self, isolated_cwd):
        """`cls._env_run` was undefined, so `debug=True` raised NameError at save() on every
        conda/mamba site -- which is every cluster variant."""
        from biopipelines import Pipeline, Resources
        from biopipelines.mock import Mock

        pipeline = Pipeline(project="TestSuite", job="dbg", description="d", on_the_fly=False,
                            local_output=True, config="cluster", debug=True)
        with pipeline:
            Resources(partition="gpu", time="01:00:00")
            Mock(ids=["a"], streams={"s": {"format": "pdb", "file": "<id>.pdb"}})
            pipeline.save()


class TestRedactionCoversTheTokensInUse:
    """The first version matched GitLab's `glpat-` and missed GitHub's `ghp_`.

    The lab has repositories on both hosts. A pattern that covers one and not the other is worse
    than none, because the export looks redacted.
    """

    @staticmethod
    def redact(text):
        import re
        import subprocess
        from biopipelines.pipeline import REDACT_SED

        done = subprocess.run(["sed"] + REDACT_SED.strip()[4:].split(" ", 0) or ["sed"],
                              input=text, capture_output=True, text=True)
        return done.stdout

    def test_every_token_prefix_is_covered(self):
        from biopipelines.pipeline import _REDACT_PATTERNS
        import re

        samples = {
            "glpat": "https://oauth2:glpat-AAAABBBBCCCCDDDD@gitlab.uzh.ch/locbp/x.git",
            "ghp_": "https://ghp_AAAABBBBCCCCDDDDEEEEFFFFGGGG@github.com/locbp-uzh/x.git",
            "gho_": "https://gho_AAAABBBBCCCCDDDD@github.com/locbp-uzh/x.git",
            "github_pat_": "https://github_pat_AAAABBBBCCCC@github.com/locbp-uzh/x.git",
        }
        # Translate the sed patterns to Python for a dependency-free check of the same regexes.
        rules = []
        for pattern in _REDACT_PATTERNS:
            _, find, replace, _flags = pattern.split("#")
            rules.append((find.replace("[:space:]", r"\s"), replace))
        for label, text in samples.items():
            out = text
            for find, replace in rules:
                out = re.sub(find, replace.replace("\1", r"\1"), out)
            assert label.rstrip("_-") not in out or "REDACTED" in out, (
                f"{label} survived redaction: {out}")
            assert "AAAABBBB" not in out, f"{label} token body survived: {out}"

    def test_a_token_outside_a_url_is_redacted_too(self):
        """The URL rule catches a token in a clone URL. The prefix rule is what catches one in a
        pip line, an env dump or a comment, where there is no `://…@` around it."""
        import re
        from biopipelines.pipeline import _REDACT_PATTERNS

        rules = []
        for pattern in _REDACT_PATTERNS:
            _, find, replace, _flags = pattern.split("#")
            rules.append((find.replace("[:space:]", r"\s"), replace))
        for token in ("ghp_AAAABBBBCCCCDDDD", "gho_AAAABBBBCCCCDDDD", "glpat-AAAABBBBCCCCDDDD"):
            out = f"# left over from a manual install: {token}"
            for find, replace in rules:
                out = re.sub(find, replace.replace("\1", r""), out)
            assert "AAAABBBB" not in out, f"bare {token[:5]} survived: {out}"

    @staticmethod
    def python_redact(text):
        """The sed rules applied with Python's re, flags and back-references translated."""
        import re
        from biopipelines.pipeline import _REDACT_PATTERNS
        for pattern in _REDACT_PATTERNS:
            _, find, replace, flags = pattern.split("#")
            find = find.replace("[:space:]", r"\s")
            replace = re.sub(r"\\(\d)", r"\\g<\1>", replace)
            text = re.sub(find, replace, text, flags=re.I if "I" in flags else 0)
        return text

    @pytest.mark.parametrize("line, secret", [
        ("  HF_TOKEN: hf_AAAABBBBCCCCDDDDEEEE", "AAAABBBB"),
        ("  WANDB_API_KEY: 0123456789abcdef", "0123456789"),
        ("  AWS_SECRET_ACCESS_KEY: abcdEFGH1234", "abcdEFGH"),
        ("PASSWORD=hunter2", "hunter2"),
        ("OPENAI_API_KEY=sk-AAAABBBBCCCCDDDDEEEE", "AAAABBBB"),
        ("https://conda.anaconda.org/t/abcdefgh1234/mychannel", "abcdefgh1234"),
        ("https://user:pa@ss@host.org/simple", "ss@host"),
        ("https://__token__:pypi-AAAABBBBCCCCDDDDEEEE@upload.pypi.org", "AAAABBBB"),
    ])
    def test_the_forms_1_5_0_missed_are_redacted(self, line, secret):
        assert secret not in self.python_redact(line)

    @pytest.mark.parametrize("line", [
        "tokenizers==0.15.0",
        "  - tokenizers=0.15.0=py310h1_0",
        "git+https://github.com/org/repo.git@v1.0",
    ])
    def test_package_lines_survive_intact(self, line):
        """Redacting a version pin would change the environment digest for no secret at all."""
        assert self.python_redact(line) == line

    def test_the_unredacted_capture_never_lands_in_the_output_root(self):
        import inspect
        from biopipelines import pipeline
        source = inspect.getsource(pipeline.Pipeline)
        assert "umask 077; mktemp" in source and ".txt.raw" not in source


class TestTheCaptureReachesEveryRunShape:

    def test_a_multi_batch_run_still_captures_the_environment(self, isolated_cwd):
        """A mixed-resource campaign runs `pipeline_batch<N>.sh`, never `RunTime/pipeline.sh`.

        The capture block had one call site, in the single-script generator, so every pipeline
        with more than one `Resources(...)` -- the normal shape on a cluster -- wrote a manifest
        with no environment digests beside it.
        """
        from biopipelines import Pipeline, Resources
        from biopipelines.mock import Mock
        from biopipelines.protein_mpnn import ProteinMPNN

        pipeline = Pipeline(project="TestSuite", job="batched", description="d",
                            on_the_fly=False, local_output=True, config="cluster")
        with pipeline:
            Resources(gpu="A100", time="00:10:00")
            source = Mock(ids=["d0"], streams={"structures": {"format": "pdb", "file": "<id>.pdb"}})
            Resources(gpu="none", time="00:10:00")
            ProteinMPNN(structures=source.streams.structures, num_sequences=1)
            pipeline.save()

        assert len(pipeline.batch_resources) > 1, (
            "two Resources() calls must produce two batches or this proves nothing")
        first = pipeline._generate_batch_script(0, 0, len(pipeline.tools))
        assert "BioPipelines environment capture" in first, (
            "the batch scripts are what actually run; without the block here a mixed-resource "
            "campaign records no environment at all")
        assert ".sha256" in first


def test_a_conda_pin_named_like_a_secret_keeps_its_version():
    """`tiktoken=0.5.1=pypi_0` matched the key rule and lost its version, changing the env digest."""
    import shutil
    import subprocess
    from biopipelines.pipeline import REDACT_SED
    if not shutil.which("bash"):
        pytest.skip("needs bash and GNU sed")
    text = "  - tiktoken=0.5.1=pypi_0\n  HF_TOKEN: hf_AAAABBBBCCCCDDDDEEEE\n"
    out = subprocess.run(["bash", "-c", REDACT_SED.strip()], input=text,
                         capture_output=True, text=True).stdout
    assert "tiktoken=0.5.1=pypi_0" in out and "AAAABBBB" not in out
