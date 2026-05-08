"""Unit tests for biopipelines.folders — placeholder resolution.

FolderManager resolves `<biopipelines>`, `<cwd>`, `<scratch>` etc. against
`folders.base` entries in the config. A regression here breaks every real
pipeline run, so the minimal guarantees are worth locking down.
"""


def test_folder_manager_resolves_placeholders(local_config, isolated_cwd, record_case):
    """FolderManager expands <biopipelines> and <cwd> into absolute paths."""
    from biopipelines.folders import FolderManager

    fm = FolderManager(project="TestSuite", job="folders", local_output=True)
    folders = fm.get_folders()

    record_case(
        input="FolderManager(local_output=True)",
        expected=("biopipelines / pipe_scripts keys present, no <placeholders>", True),
        actual=("keys",
                "biopipelines" in folders and "pipe_scripts" in folders
                and "<" not in folders["pipe_scripts"]),
    )
    assert "biopipelines" in folders
    assert "pipe_scripts" in folders
    # Placeholders must be fully resolved (no "<biopipelines>" left over).
    for name, path in folders.items():
        assert "<" not in path, (
            f"FolderManager left an unresolved placeholder in {name}: {path}"
        )


def test_folder_manager_biopipelines_output_contains_project_job(
    local_config, isolated_cwd, record_case,
):
    """biopipelines_output resolves to a real directory under cwd when
    local_output=True (doesn't need to include project/job — those are layered
    on by Pipeline itself — but must be rooted inside cwd)."""
    import os
    from biopipelines.folders import FolderManager

    fm = FolderManager(project="TestSuite", job="folders2", local_output=True)
    folders = fm.get_folders()
    out = folders["biopipelines_output"]

    record_case(
        input="biopipelines_output under cwd",
        expected=("absolute path rooted in cwd", True),
        actual=("is abs and under cwd",
                os.path.isabs(out) and str(isolated_cwd) in out),
    )
    assert os.path.isabs(out)
    assert str(isolated_cwd) in out
