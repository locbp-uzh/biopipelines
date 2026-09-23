"""PROJECT.md and HISTORY.md, maintained by an agent without destroying what a human wrote.

The whole module rests on one rule: never overwrite. PROJECT.md carries the scientific
judgment — the goal, the inputs, and the validated protocol — and an agent that silently
rewrites it destroys the reason the folder is worth keeping. These tests exist to make that
rule hold under the obvious ways it would be broken.
"""

import io

from biopipelines import project_docs


def test_scaffold_creates_the_two_documents(tmp_path):
    result = project_docs.scaffold(tmp_path / "SNAP33")
    assert result == {"project": "created", "history": "created"}
    assert (tmp_path / "SNAP33" / "PROJECT.md").exists()
    assert (tmp_path / "SNAP33" / "HISTORY.md").exists()


def test_the_readme_is_opt_in(tmp_path):
    assert "readme" not in project_docs.scaffold(tmp_path / "A")
    assert project_docs.scaffold(tmp_path / "B", readme=True)["readme"] == "created"


def test_scaffold_never_overwrites_existing_content(tmp_path):
    project = tmp_path / "SNAP33"
    project.mkdir()
    written = "# SNAP33\n\nMy own goal section, hand-written.\n"
    io.open(project / "PROJECT.md", "w", encoding="utf-8").write(written)

    result = project_docs.scaffold(project)
    assert result == {"project": "kept", "history": "created"}
    assert io.open(project / "PROJECT.md", encoding="utf-8").read() == written, (
        "an agent must never clobber the document that carries the science")


def test_scaffold_on_an_established_project_fills_only_the_gap(tmp_path):
    """Running it on a project that predates the convention should add the missing half."""
    project = tmp_path / "TNT-Snap"
    project.mkdir()
    io.open(project / "PROJECT.md", "w", encoding="utf-8").write("# TNT-Snap\n")
    assert project_docs.scaffold(project)["history"] == "created"


def test_the_project_name_defaults_to_the_folder(tmp_path):
    project_docs.scaffold(tmp_path / "SNAP33")
    assert "# SNAP33" in io.open(tmp_path / "SNAP33" / "PROJECT.md", encoding="utf-8").read()


def test_history_appends_in_order_and_keeps_earlier_entries(tmp_path):
    project = tmp_path / "SNAP33"
    project_docs.scaffold(project)
    project_docs.append_history(project, "First run", "500 designs.", date="2026-06-01")
    project_docs.append_history(project, "Second run", date="2026-06-02")

    text = io.open(project / "HISTORY.md", encoding="utf-8").read()
    assert text.index("First run") < text.index("Second run"), "newest goes at the bottom"
    assert "500 designs." in text


def test_history_creates_the_file_when_it_is_missing(tmp_path):
    project_docs.append_history(tmp_path / "New", "Started", date="2026-06-01")
    assert (tmp_path / "New" / "HISTORY.md").exists()
    assert (tmp_path / "New" / "PROJECT.md").exists(), "scaffolding fills both"


def test_history_entries_are_readable_back(tmp_path):
    project = tmp_path / "SNAP33"
    project_docs.append_history(project, "Run 006", date="2026-06-01")
    project_docs.append_history(project, "Reorganized the tree", date="2026-09-14")
    assert project_docs.history_entries(project) == [
        {"date": "2026-06-01", "title": "Run 006"},
        {"date": "2026-09-14", "title": "Reorganized the tree"}]


def test_job_folders_are_found_newest_first(tmp_path):
    project = tmp_path / "SNAP33"
    for name in ("Binder_design_006", "Binder_design_007", "inputs", "scripts"):
        (project / name).mkdir(parents=True)
    # Names, not paths: the same call has to answer for a folder on a cluster, where a local
    # `Path` means nothing.
    assert project_docs.job_dirs(project) == ["Binder_design_007", "Binder_design_006"], (
        "only <name>_NNN folders are jobs; inputs/ and scripts/ are not")


def test_survey_reports_documents_jobs_and_entries(tmp_path):
    project = tmp_path / "SNAP33"
    (project / "Binder_design_006").mkdir(parents=True)
    project_docs.scaffold(project)
    project_docs.append_history(project, "Run 006", date="2026-06-01")

    state = project_docs.survey(project)
    assert state["documents"] == {"project": True, "history": True, "readme": False}
    assert state["jobs"] == ["Binder_design_006"]
    assert state["entries"] == [{"date": "2026-06-01", "title": "Run 006"}]


def test_survey_of_a_folder_that_is_not_a_project_is_empty_not_an_error(tmp_path):
    state = project_docs.survey(tmp_path / "nothing-here")
    assert state["jobs"] == [] and state["entries"] == []
    assert state["documents"] == {"project": False, "history": False, "readme": False}


def test_a_missing_folder_is_not_an_empty_project(tmp_path):
    """The distinction that matters: "not there" and "there but empty" need different answers.

    `bp_project` had no `host`, so a cluster project folder was read through local pathlib and
    reported as a project with no documents and no jobs — a description of an empty project,
    for a project holding two runs on another machine.
    """
    state = project_docs.survey(tmp_path / "NeverCreated")
    assert state["exists"] is False
    assert state["jobs"] == [] and state["entries"] == []
    assert not any(state["documents"].values())


def test_an_existing_but_bare_project_says_it_exists(tmp_path):
    bare = tmp_path / "Bare"
    bare.mkdir()
    state = project_docs.survey(bare)
    assert state["exists"] is True and not any(state["documents"].values())
