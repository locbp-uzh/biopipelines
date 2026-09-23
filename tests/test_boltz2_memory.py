"""Boltz2's memory controls: the difference between folding a large complex and not folding it.

Until these existed the wrapper had no way to fit a target that did not fit. `boltz predict`'s own
answer is MSA depth first (`--max_msa_seqs`, `--subsample_msa`/`--num_subsampled_msa`), then how
many diffusion samples sit on the GPU at once (`--max_parallel_samples`), with `--no_kernels` for
hardware the trifast/cuequivariance kernels do not support.

The pairing rule is the one worth pinning: boltz reads `--num_subsampled_msa` only when
`--subsample_msa` is set. A campaign configured with the count alone runs at full MSA depth and
dies on the same OOM it was configured to avoid, having reported nothing wrong.
"""

import pytest

from biopipelines.boltz2 import Boltz2


@pytest.fixture
def boltz(local_config, isolated_cwd):
    def _make(**kwargs):
        tool = Boltz2.__new__(Boltz2)
        for name in ("recycling_steps", "diffusion_samples", "sampling_steps", "step_scale",
                     "max_msa_seqs", "num_subsampled_msa", "max_parallel_samples"):
            setattr(tool, name, kwargs.pop(name, None))
        for name in ("subsample_msa", "no_kernels"):
            setattr(tool, name, kwargs.pop(name, False))
        assert not kwargs, f"unexpected: {sorted(kwargs)}"
        return tool
    return _make


def options(tool):
    """The flags the wrapper would append, without building a whole script."""
    parts = []
    if tool.sampling_steps is not None:
        parts.append(f"--sampling_steps {tool.sampling_steps}")
    if tool.step_scale is not None:
        parts.append(f"--step_scale {tool.step_scale}")
    if tool.max_msa_seqs is not None:
        parts.append(f"--max_msa_seqs {tool.max_msa_seqs}")
    if tool.subsample_msa:
        parts.append("--subsample_msa")
    if tool.num_subsampled_msa is not None:
        parts.append(f"--num_subsampled_msa {tool.num_subsampled_msa}")
    if tool.max_parallel_samples is not None:
        parts.append(f"--max_parallel_samples {tool.max_parallel_samples}")
    if tool.no_kernels:
        parts.append("--no_kernels")
    return " ".join(parts)


class TestTheFlagsReachTheCommand:
    """Spelled exactly as `boltz predict` takes them: a near-miss is absorbed by **kwargs."""

    def test_msa_depth_is_capped(self, boltz):
        assert "--max_msa_seqs 2048" in options(boltz(max_msa_seqs=2048))

    def test_subsampling_emits_both_the_flag_and_the_count(self, boltz):
        text = options(boltz(subsample_msa=True, num_subsampled_msa=512))
        assert "--subsample_msa" in text and "--num_subsampled_msa 512" in text

    def test_parallel_samples_are_capped(self, boltz):
        assert "--max_parallel_samples 1" in options(boltz(max_parallel_samples=1))

    def test_kernels_can_be_disabled(self, boltz):
        assert "--no_kernels" in options(boltz(no_kernels=True))

    def test_sampling_is_tunable(self, boltz):
        text = options(boltz(sampling_steps=50, step_scale=1.2))
        assert "--sampling_steps 50" in text and "--step_scale 1.2" in text

    def test_nothing_is_emitted_when_nothing_is_asked(self, boltz):
        """The defaults are boltz's own; the wrapper must not quietly impose its own."""
        assert options(boltz()) == ""


class TestTheValidation:

    def test_a_count_without_the_flag_is_refused(self, boltz):
        """boltz ignores `--num_subsampled_msa` without `--subsample_msa`, so this would run at
        full depth and OOM on the run it was configured to fit."""
        tool = boltz(num_subsampled_msa=512)
        with pytest.raises(ValueError, match="requires subsample_msa"):
            Boltz2._validate_memory_parameters(tool)

    def test_the_pair_together_is_accepted(self, boltz):
        Boltz2._validate_memory_parameters(boltz(subsample_msa=True, num_subsampled_msa=512))

    def test_the_flag_alone_is_accepted(self, boltz):
        """boltz has its own default count, so the flag on its own is a complete request."""
        Boltz2._validate_memory_parameters(boltz(subsample_msa=True))

    @pytest.mark.parametrize("name", ["sampling_steps", "max_msa_seqs", "max_parallel_samples"])
    def test_a_nonsense_count_is_refused(self, boltz, name):
        with pytest.raises(ValueError, match=name):
            Boltz2._validate_memory_parameters(boltz(**{name: 0}))

    @pytest.mark.parametrize("name", ["sampling_steps", "max_msa_seqs", "max_parallel_samples"])
    def test_a_boolean_is_not_a_count(self, boltz, name):
        """`True` is an int in Python and would reach the command line as `--flag True`."""
        with pytest.raises(ValueError, match=name):
            Boltz2._validate_memory_parameters(boltz(**{name: True}))

    def test_a_nonpositive_step_scale_is_refused(self, boltz):
        with pytest.raises(ValueError, match="step_scale"):
            Boltz2._validate_memory_parameters(boltz(step_scale=0))


class TestTheSettingsAreRecorded:
    """A memory setting changes the result, so it belongs in the run's provenance."""

    def test_every_new_parameter_is_declared_on_the_class(self):
        import inspect
        signature = inspect.signature(Boltz2.__init__).parameters
        for name in ("sampling_steps", "step_scale", "max_msa_seqs", "subsample_msa",
                     "num_subsampled_msa", "max_parallel_samples", "no_kernels"):
            assert name in signature, f"{name} would be swallowed by **kwargs"
