# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""RFdiffusion under fold conditioning: the ids the step declares are the files its run writes.

Fold conditioning (`scaffold_dir`, `target_ss`, `target_adj`) is one run named after the pipeline step, driven by `scaffoldguided.*` alone: upstream takes the target of a fold-conditioned binder from `scaffoldguided.target_path`, never from `inference.input_pdb`. `get_output_files` nevertheless keyed the ids off a `pdb` input whenever one was given, so a one-structure fold-conditioned step declared `p1_<1..2>` while the run wrote `<step>_1.pdb` and `<step>_2.pdb`: every declared output was absent and the step was marked FAILED. A `pdb` input that fold conditioning would silently ignore is now refused, and `hotspot_res` is accepted there with a `target_path`, the way upstream pairs them.

Nothing here runs RFdiffusion; the declared ids are compared with the `inference.output_prefix` of the generated script.
"""

import pathlib
import re

import pytest

from biopipelines.rfdiffusion import RFdiffusion

FOLD = {
    "scaffold_dir": {"scaffold_dir": "/scaffolds"},
    "target_ss_and_adj": {"target_ss": "/t_ss.pt", "target_adj": "/t_adj.pt"},
    "binder_with_target": {"scaffold_dir": "/scaffolds", "target_path": "/target.pdb",
                           "target_ss": "/t_ss.pt", "target_adj": "/t_adj.pt"},
}


def _structures(ids=("p1",)):
    from biopipelines.mock import Mock
    return Mock(ids=list(ids), streams={"structures": {"format": "pdb", "file": "<id>.pdb"}},
                map_table_strategy="config")


def _step(new_pipeline, with_pdb=False, **kwargs):
    """The tool's output and its step script, built inside a real pipeline."""
    pipeline = new_pipeline("fold_conditioning")
    with pipeline:
        if with_pdb:
            kwargs["pdb"] = _structures()
        built = RFdiffusion(num_designs=2, **kwargs)
        script_path = pipeline.save()
    script = next(pathlib.Path(script_path).parent.glob("*_RFdiffusion.sh")).read_text(encoding="utf-8", errors="replace")
    return built, script


class TestDeclaredIdsAreTheFilesTheRunWrites:

    @pytest.mark.parametrize("mode", list(FOLD))
    def test_the_declared_ids_are_the_output_prefix_of_the_run(self, mode, local_config, isolated_cwd, new_pipeline):
        built, script = _step(new_pipeline, **FOLD[mode])
        prefixes = re.findall(r'"inference\.output_prefix=([^"]+)"', script)
        assert len(prefixes) == 1, prefixes
        stem = pathlib.PurePosixPath(prefixes[0].replace("\\", "/")).name
        assert list(built.streams.structures.ids) == [f"{stem}_<1..2>"]

    @pytest.mark.parametrize("mode", list(FOLD))
    def test_no_accepted_combination_with_a_pdb_declares_ids_the_run_does_not_write(self, mode, local_config,
                                                                                    isolated_cwd, new_pipeline):
        try:
            built, script = _step(new_pipeline, with_pdb=True, contigs="A1-5/10-20", **FOLD[mode])
        except ValueError:
            return
        stem = pathlib.PurePosixPath(re.findall(r'"inference\.output_prefix=([^"]+)"', script)[0].replace("\\", "/")).name
        assert list(built.streams.structures.ids) == [f"{stem}_<1..2>"]

    @pytest.mark.parametrize("mode", list(FOLD))
    def test_the_provenance_command_names_no_parent(self, mode, local_config, isolated_cwd, new_pipeline):
        _, script = _step(new_pipeline, **FOLD[mode])
        update = script[script.index("pipe_update_structures_map.py"):]
        assert "--parents-json" not in update[:update.index("\n")]


class TestAPdbInputIsRefusedWhereFoldConditioningWouldIgnoreIt:

    @pytest.mark.parametrize("mode", list(FOLD))
    def test_pdb_with_fold_conditioning_is_refused(self, mode, local_config, isolated_cwd, new_pipeline):
        pipeline = new_pipeline("fold_with_pdb")
        with pytest.raises(ValueError, match="pdb"):
            with pipeline:
                RFdiffusion(pdb=_structures(), contigs="A1-5/10-20", **FOLD[mode])

    def test_pdb_without_fold_conditioning_still_declares_one_run_per_input(self, local_config, isolated_cwd,
                                                                            new_pipeline):
        pipeline = new_pipeline("pdb_only")
        with pipeline:
            built = RFdiffusion(pdb=_structures(("p1", "q2")), contigs="A1-5/10-20", num_designs=2)
        assert list(built.streams.structures.ids) == ["p1_<1..2>", "q2_<1..2>"]


class TestHotspotsFollowTheTargetFoldConditioningNames:

    def test_hotspots_with_a_target_path_are_accepted(self, local_config, isolated_cwd, new_pipeline):
        _, script = _step(new_pipeline, hotspot_res=["A59", "A83"], **FOLD["binder_with_target"])
        assert '"ppi.hotspot_res=[A59,A83]"' in script
        assert "inference.input_pdb" not in script

    def test_hotspots_without_any_target_are_still_refused(self, local_config, isolated_cwd, new_pipeline):
        pipeline = new_pipeline("hotspots_alone")
        with pytest.raises(ValueError, match="hotspot_res"):
            with pipeline:
                RFdiffusion(hotspot_res=["A59"], **FOLD["scaffold_dir"])
