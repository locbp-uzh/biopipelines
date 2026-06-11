"""Config overlay mechanism: committed base + gitignored .config.<v>.yaml.

Covers the pure merge/diff helpers and the end-to-end load path that deep-merges
a user overlay onto the committed base, so adding keys to the repo config never
forces a user to redefine their local edits.
"""

import textwrap

import pytest

from biopipelines.config_manager import ConfigManager, _deep_merge
from biopipelines.config_editor import _overlay_diff


# ── _deep_merge ───────────────────────────────────────────────────────────────

def test_deep_merge_recurses_into_nested_dicts():
    base = {"a": {"x": 1, "y": 2}, "b": 3}
    overlay = {"a": {"y": 20, "z": 30}, "c": 4}
    assert _deep_merge(base, overlay) == {
        "a": {"x": 1, "y": 20, "z": 30},
        "b": 3,
        "c": 4,
    }


def test_deep_merge_does_not_mutate_base():
    base = {"a": {"x": 1}}
    overlay = {"a": {"x": 2}}
    _deep_merge(base, overlay)
    assert base == {"a": {"x": 1}}


def test_deep_merge_result_shares_no_nested_dict_with_base():
    # The editor edits the merged result in place; it must not reach back into
    # the base. This includes base-only nested dicts absent from the overlay —
    # the bug that made the editor save an empty diff (it mutated the base, so
    # diff-against-base came out empty).
    base = {"section": {"leaf": "orig"}, "other": {"k": 1}}
    merged = _deep_merge(base, {})
    merged["section"]["leaf"] = "edited"
    merged["other"]["k"] = 99
    assert base == {"section": {"leaf": "orig"}, "other": {"k": 1}}


def test_deep_merge_scalar_and_list_replace_whole():
    # A non-dict overlay value replaces the base value outright (no list-merge).
    base = {"k": [1, 2, 3], "s": {"deep": 1}}
    overlay = {"k": [9], "s": "scalar-now"}
    assert _deep_merge(base, overlay) == {"k": [9], "s": "scalar-now"}


def test_deep_merge_overlay_none_replaces():
    base = {"k": "value"}
    assert _deep_merge(base, {"k": None}) == {"k": None}


# ── _overlay_diff ─────────────────────────────────────────────────────────────

def test_overlay_diff_keeps_only_deviations():
    base = {"a": {"x": 1, "y": 2}, "b": 3}
    merged = {"a": {"x": 1, "y": 99}, "b": 3, "c": 5}
    assert dict(_overlay_diff(merged, base)) == {"a": {"y": 99}, "c": 5}


def test_overlay_diff_empty_when_identical():
    base = {"a": {"x": 1}, "b": 2}
    assert dict(_overlay_diff(dict(base), base)) == {}


def test_overlay_diff_does_not_represent_deletions():
    # A key present in base but absent in merged is NOT a diff entry — the
    # overlay only adds/overrides, it never deletes a base key.
    base = {"a": 1, "b": 2}
    merged = {"a": 1}
    assert dict(_overlay_diff(merged, base)) == {}


def test_overlay_diff_new_nested_key():
    base = {"section": {"x": 1}}
    merged = {"section": {"x": 1, "new": 2}}
    assert dict(_overlay_diff(merged, base)) == {"section": {"new": 2}}


# ── end-to-end _load_config merge ─────────────────────────────────────────────

def _minimal_base() -> str:
    return textwrap.dedent(
        """
        folders:
          base:
            home: /home/u
          infrastructure:
            scripts: <home>/my_scripts
          repositories: {}
        machine:
          username: ""
          email: ""
          env_manager:
            name: mamba
            init: []
          scheduler:
            name: none
            init: []
        environments:
          Foo: bar_env
        """
    )


@pytest.fixture
def base_with_overlay(monkeypatch, tmp_path):
    """Write a temp base config + overlay and point ConfigManager at them.

    Yields a function that writes overlay YAML text (or removes the overlay)
    and returns a freshly loaded ConfigManager.
    """
    base_path = tmp_path / "config.cluster.yaml"
    overlay_path = tmp_path / ".config.cluster.yaml"
    base_path.write_text(_minimal_base(), encoding="utf-8")

    monkeypatch.setattr(
        ConfigManager, "_get_config_path",
        classmethod(lambda cls, variant=None: str(base_path)),
    )
    monkeypatch.setattr(
        ConfigManager, "_get_overlay_path",
        classmethod(lambda cls, variant=None: str(overlay_path)),
    )

    def _load(overlay_text=None):
        if overlay_text is None:
            if overlay_path.exists():
                overlay_path.unlink()
        else:
            overlay_path.write_text(textwrap.dedent(overlay_text), encoding="utf-8")
        ConfigManager._instance = None
        ConfigManager._config = None
        ConfigManager._variant = "cluster"
        return ConfigManager(variant="cluster")

    yield _load

    ConfigManager._instance = None
    ConfigManager._config = None
    ConfigManager._variant = None


def test_load_without_overlay_returns_base(base_with_overlay):
    cm = base_with_overlay(None)
    assert cm.get_environment("Foo") == "bar_env"
    assert cm.get_email() == ""


def test_overlay_overrides_a_leaf(base_with_overlay):
    cm = base_with_overlay(
        """
        machine:
          email: me@example.com
        """
    )
    assert cm.get_email() == "me@example.com"
    # Untouched base keys survive.
    assert cm.get_env_manager() == "mamba"
    assert cm.get_environment("Foo") == "bar_env"


def test_overlay_adds_a_new_key_without_redefining_siblings(base_with_overlay):
    # The whole point: an overlay touching one tool env leaves the rest intact.
    cm = base_with_overlay(
        """
        environments:
          Baz: baz_env
        """
    )
    assert cm.get_environment("Baz") == "baz_env"
    assert cm.get_environment("Foo") == "bar_env"


def test_overlay_overrides_scripts_folder(base_with_overlay):
    cm = base_with_overlay(
        """
        folders:
          infrastructure:
            scripts: /custom/scripts
        """
    )
    assert cm.get_scripts_folder() == "/custom/scripts"


def test_get_scripts_folder_resolves_placeholders(base_with_overlay):
    cm = base_with_overlay(None)
    # <home> from folders.base is substituted.
    assert cm.get_scripts_folder() == "/home/u/my_scripts"


def test_non_mapping_overlay_raises(base_with_overlay):
    with pytest.raises(ValueError):
        base_with_overlay("- just\n- a\n- list\n")
