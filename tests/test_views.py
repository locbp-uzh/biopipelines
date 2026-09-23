"""Showing a step's results instead of describing them.

`llm/pipelines.md` has told agents for a while to render results rather than answer in prose,
and the tool layer could not honor it: `bp-visualize` is a console script, so an agent with only
the MCP tools handed the user a shell command and a path.

The split these tests pin is the one that decides what a chat can actually show. A step's images
are files, so they can come back inline; the page's structure views are 3Dmol.js running in a
browser and cannot. Claiming otherwise in either direction is the failure worth catching — an
agent that thinks it showed a structure will not offer the page.
"""

import pytest

from biopipelines import views


class FakeSsh:
    """A cluster with a file tree, for the parts of the render that look before they act."""

    def __init__(self, tree=None, output="", code=0):
        self.host, self.repo, self.timeout = "cluster", "~/biopipelines", 30
        self.tree = tree or {}
        self.output, self.code = output, code
        self.commands = []

    def run(self, command, timeout=None, stdin=None):
        self.commands.append(command)
        return self.code, self.output, ""

    def join(self, *parts):
        return "/".join(str(p).rstrip("/") for p in parts if str(p) != "")

    def is_dir(self, path):
        return str(path) in self.tree

    def listdir(self, path):
        return list(self.tree.get(str(path), []))


STEP = "/runs/job_001/003_Boltz2"


def tree(**extra):
    base = {
        "/runs/job_001/003_Boltz2": [("structures", True), ("tables", True),
                                     ("_configuration", True), ("_extras", True)],
        "/runs/job_001/003_Boltz2/structures": [("d1.pdb", False)],
        "/runs/job_001/003_Boltz2/tables": [("scores.csv", False)],
        "/runs/job_001/003_Boltz2/_configuration": [("secret.png", False)],
        "/runs/job_001/003_Boltz2/_extras": [("003_Boltz2_view.html", False)],
    }
    base.update(extra)
    return base


class TestWhichFilesCanBeShown:
    def test_a_steps_images_are_found(self):
        ssh = FakeSsh(tree(**{"/runs/job_001/003_Boltz2/_extras": [
            ("003_Boltz2_view.html", False), ("affinity.png", False)]}))
        names = [name for _path, name in views.images(ssh, "/runs/job_001", "003_Boltz2")]
        assert names == ["affinity.png"]

    def test_framework_folders_are_not_searched_for_images(self):
        """`_configuration` holds plumbing; an image there is not a result worth showing."""
        names = [name for _p, name in views.images(FakeSsh(tree()), "/runs/job_001", "003_Boltz2")]
        assert names == []

    def test_non_images_are_left_alone(self):
        ssh = FakeSsh(tree(**{"/runs/job_001/003_Boltz2/tables": [
            ("scores.csv", False), ("plot.PNG", False)]}))
        names = [name for _p, name in views.images(ssh, "/runs/job_001", "003_Boltz2")]
        assert names == ["plot.PNG"], "the suffix check must not be case-sensitive"

    def test_the_number_shown_is_capped(self):
        many = [(f"p{i}.png", False) for i in range(20)]
        ssh = FakeSsh(tree(**{"/runs/job_001/003_Boltz2/tables": many}))
        assert len(views.images(ssh, "/runs/job_001", "003_Boltz2")) == views.MAX_IMAGES


class TestRender:
    def test_it_runs_on_the_login_node_with_no_scheduler(self):
        ssh = FakeSsh(tree(), output="wrote /runs/job_001/003_Boltz2/_extras/v.html\n")
        result = views.render(ssh, "/runs/job_001", "003_Boltz2")
        assert result["ok"]
        assert "bp-visualize" in ssh.commands[0]
        assert "sbatch" not in ssh.commands[0] and "srun" not in ssh.commands[0]

    def test_ordering_and_the_cap_reach_the_command(self):
        """'Show me the best five' is answered by the tool's own flags, not by post-filtering."""
        ssh = FakeSsh(tree())
        views.render(ssh, "/runs/job_001", "003_Boltz2",
                     descending="affinity.pred", max_items=5)
        assert "--descending" in ssh.commands[0] and "affinity.pred" in ssh.commands[0]
        assert "--max-items 5" in ssh.commands[0]

    def test_ascending_is_ignored_when_descending_is_given(self):
        ssh = FakeSsh(tree())
        views.render(ssh, "/runs/job_001", "003_Boltz2", descending="a.b", ascending="c.d")
        assert "--ascending" not in ssh.commands[0]

    def test_a_missing_step_is_reported_before_anything_runs(self):
        ssh = FakeSsh({})
        result = views.render(ssh, "/runs/job_001", "003_Boltz2")
        assert not result["ok"] and ssh.commands == []

    def test_a_failing_render_carries_its_output(self):
        ssh = FakeSsh(tree(), output="Traceback...\n", code=1)
        result = views.render(ssh, "/runs/job_001", "003_Boltz2")
        assert not result["ok"]
        assert "Traceback" in views.summarize(result)

    @pytest.mark.parametrize("said,expected", [
        ("wrote /a/b/view.html", "/a/b/view.html"),
        ("Written to: /a/b/view.html\n", "/a/b/view.html"),
        ("nothing useful", "/runs/job_001/003_Boltz2/_extras/003_Boltz2_view.html"),
    ])
    def test_the_page_path_comes_from_the_output_or_the_default(self, said, expected):
        assert views.page_path(said, FakeSsh(), "/runs/job_001", "003_Boltz2") == expected


class TestSummary:
    def test_it_says_where_to_open_the_interactive_view(self):
        result = {"ok": True, "page": "/remote/v.html", "failures": []}
        text = views.summarize(result, page_local="C:/tmp/v.html", shown=["affinity.png"])
        assert "C:/tmp/v.html" in text and "interactive" in text
        assert "affinity.png" in text

    def test_a_step_with_no_images_says_so_rather_than_implying_a_picture(self):
        """An agent that believes it showed a structure will not offer the page."""
        text = views.summarize({"ok": True, "page": "/remote/v.html", "failures": []})
        assert "nothing to show inline" in text
        assert "only" in text and "browser" in text

    def test_an_image_too_large_to_send_is_named_not_dropped(self):
        text = views.summarize({"ok": True, "page": "/p.html", "failures": []},
                               shown=["a.png"], skipped=["huge.png"])
        assert "huge.png" in text and "in the page instead" in text
