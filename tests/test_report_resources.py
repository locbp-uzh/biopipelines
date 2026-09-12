# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""The run page reports the resources a step asked for, and `gpus` is a count that only means something once a `gpu` spec asks for GPUs at all. `gpus` defaults to 1 whether or not one was requested, so a page that prints it unconditionally tells the reader a CPU-only job took a GPU. The scheduler is the authority here: `gpu_directive` returns no directive when the spec is None, "none" or empty, whatever the count says, and these tests hold the report to the same rule."""

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from biopipelines.schedulers import SlurmBackend  # noqa: E402


def load_renderer():
    spec = importlib.util.spec_from_file_location("pipeline_report", ROOT / "biopipelines" / "renderers" / "pipeline_report.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["pipeline_report"] = module
    spec.loader.exec_module(module)
    return module


DECLINED = [None, "none", "", "None"]


def test_a_declined_gpu_is_not_reported_as_one_gpu():
    fmt = load_renderer()._fmt_resources
    for spec in DECLINED:
        text = fmt({"gpu": spec, "gpus": 1, "memory": "4GB", "time": "00:20:00"})
        assert "gpus=1" not in text, f"spec {spec!r} reported a GPU: {text}"
        assert "gpu=none" in text
        assert "memory=4GB" in text and "time=00:20:00" in text


def test_a_requested_gpu_still_reports_its_count():
    fmt = load_renderer()._fmt_resources
    text = fmt({"gpu": "A100", "gpus": 2, "memory": "32GB", "time": "12:00:00"})
    assert "gpu=A100" in text and "gpus=2" in text


def test_the_report_agrees_with_the_scheduler():
    """Whatever the report says about GPUs must match whether the backend emits a directive, since that is what actually reaches SLURM."""
    fmt = load_renderer()._fmt_resources
    backend = SlurmBackend()
    for spec in DECLINED + ["A100", "gpu", "80GB"]:
        directive, _ = backend.gpu_directive(spec, 1, single_gpu_type=False)
        reported = fmt({"gpu": spec, "gpus": 1, "memory": "4GB", "time": "01:00:00"})
        claims_gpu = "gpus=1" in reported
        assert claims_gpu == bool(directive.strip()), (
            f"spec {spec!r}: report claims_gpu={claims_gpu} but directive is {directive!r}"
        )
