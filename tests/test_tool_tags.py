"""The tag vocabulary, and the lint that keeps the docs honest about it.

Tags are the only part of the tool documentation that exists purely so a tool can be *found*
without knowing its name, which means nothing in normal use breaks when one is wrong. A tag
outside the vocabulary, a tool that carries none, or a section whose body argues for a
capability its Tags line omits are all invisible until someone searches and gets the wrong
answer — so they are checked here instead.

The load-bearing test is `test_a_section_declares_the_capabilities_it_describes`: CI can
trivially check that a declared tag is spelled correctly, but only keyword evidence can notice
a *missing* one, and a missing capability tag is the failure that started this whole feature.
"""

import io
import re

import pytest

from biopipelines import tool_docs, tool_tags

TOOLS = [(e["name"], tool_docs.section(e["name"]) or "") for e in tool_docs.collect()]
TAGGED = [(name, body, tool_tags.parse_tags_line(body) or []) for name, body in TOOLS
          if name not in tool_tags.UNTAGGED]


class TestEveryToolIsTagged:
    @pytest.mark.parametrize("name,body", TOOLS)
    def test_a_tool_declares_tags_or_is_a_named_exception(self, name, body):
        declared = tool_tags.parse_tags_line(body)
        if name in tool_tags.UNTAGGED:
            assert declared is None, f"{name} is listed as untagged but carries a Tags line"
        else:
            assert declared, f"{name} has no **Tags**: line in {tool_docs.DOCS}"

    @pytest.mark.parametrize("name,_body,tags", TAGGED)
    def test_every_declared_tag_is_in_the_vocabulary(self, name, _body, tags):
        unknown = tool_tags.unknown_tags(tags)
        assert not unknown, f"{name} declares tags outside the vocabulary: {unknown}"

    @pytest.mark.parametrize("name,_body,tags", TAGGED)
    def test_every_tool_says_what_it_does(self, name, _body, tags):
        """A tool with no ACTION tag is unreachable by the query anyone actually writes."""
        action = set(tool_tags.VOCABULARY["action"])
        assert set(tags) & action, f"{name} carries no ACTION tag: {tags}"

    @pytest.mark.parametrize("name,_body,tags", TAGGED)
    def test_a_tag_is_not_listed_twice(self, name, _body, tags):
        dupes = sorted({t for t in tags if tags.count(t) > 1})
        assert not dupes, f"{name} repeats {dupes}"


class TestTheLint:
    @pytest.mark.parametrize("name,body,tags", TAGGED)
    def test_a_section_declares_the_capabilities_it_describes(self, name, body, tags):
        """The one defect CI cannot see any other way: a capability in the body, not in the tags."""
        gaps = tool_tags.missing_capability_tags(body, tags, name)
        assert not gaps, (
            f"{name}'s documentation argues for {gaps} but its Tags line omits them. Add the tag, "
            f"or add an entry to tool_tags.LINT_EXEMPT saying why the keyword means something else."
        )

    @pytest.mark.parametrize("entry", sorted(tool_tags.LINT_EXEMPT))
    def test_an_exemption_still_refers_to_something_real(self, entry):
        tool, tag = entry
        assert tool_docs.section(tool), f"LINT_EXEMPT names {tool}, which is not a documented tool"
        assert tag in tool_tags.ALL_TAGS, f"LINT_EXEMPT names {tag}, which is not a tag"

    @pytest.mark.parametrize("entry", sorted(tool_tags.LINT_EXEMPT))
    def test_an_exemption_that_no_longer_fires_is_removed(self, entry):
        """A stale exemption silently disables the lint for that tool the next time the docs change."""
        tool, tag = entry
        body = tool_docs.section(tool) or ""
        strong, _weak = tool_tags.propose(body)
        declared = tool_tags.parse_tags_line(body) or []
        assert tag in strong and tag not in declared, (
            f"LINT_EXEMPT[{tool}, {tag}] no longer suppresses anything — delete it."
        )

    def test_every_pattern_names_a_real_tag(self):
        """A pattern for a retired tag proposes something `unknown_tags()` then rejects."""
        orphans = sorted(set(tool_tags.PATTERNS) - tool_tags.ALL_TAGS)
        assert not orphans, f"PATTERNS keys that are not tags: {orphans}"

    def test_the_tags_line_is_not_evidence_for_itself(self):
        """Without stripping it, every declared tag justifies itself and the lint can never fire."""
        section = "### Fake\n\n**Tags**: covalent\n\n**Parameters**:\n- `x`: int - nothing special\n"
        assert "covalent" not in tool_tags.evidence(section)


class TestTheVocabulary:
    @pytest.mark.parametrize("facet", tool_tags.FACETS)
    def test_no_tag_appears_in_two_facets(self, facet):
        others = {t for name, tags in tool_tags.VOCABULARY.items() if name != facet for t in tags}
        overlap = sorted(set(tool_tags.VOCABULARY[facet]) & others)
        assert not overlap, f"{facet} shares tags with another facet: {overlap}"

    @pytest.mark.parametrize("tag", sorted(tool_tags.ALL_TAGS))
    def test_every_tag_is_carried_by_at_least_one_tool(self, tag):
        """A tag nobody carries returns an empty list, which reads as 'we have no such tool'."""
        carriers = [name for name, _b, tags in TAGGED if tag in tags]
        assert carriers, f"no tool carries {tag!r}; drop it from VOCABULARY or apply it"

    @pytest.mark.parametrize("tag", sorted(tool_tags.ALL_TAGS))
    def test_every_tag_is_documented_for_humans(self, tag):
        text = io.open(tool_docs.ROOT / "docs" / "tool_tags.md", encoding="utf-8").read()
        assert f"`{tag}`" in text, f"{tag} is in VOCABULARY but not in docs/tool_tags.md"

    def test_the_doc_invents_no_tags_of_its_own(self):
        """Only the definitions are checked — the prose names rejected terms on purpose."""
        text = io.open(tool_docs.ROOT / "docs" / "tool_tags.md", encoding="utf-8").read()
        listed = set(re.findall(r"^\|\s*`([a-z][a-z-]+)`\s*\|", text, re.M))
        listed |= {t for line in text.splitlines() if re.fullmatch(r"(`[a-z-]+`)( · `[a-z-]+`)+", line)
                   for t in re.findall(r"`([a-z-]+)`", line)}
        invented = sorted(listed - tool_tags.ALL_TAGS)
        assert not invented, f"docs/tool_tags.md defines tags that do not exist: {invented}"


class TestVariantsAgree:
    @pytest.mark.parametrize("variant,parent", sorted(tool_tags.TAG_ALIAS.items()))
    def test_a_variant_carries_its_parents_tags(self, variant, parent):
        """SolubleMPNN runs what ProteinMPNN runs; a search must not find one and miss the other."""
        got = tool_tags.parse_tags_line(tool_docs.section(variant) or "") or []
        want = tool_tags.parse_tags_line(tool_docs.section(parent) or "") or []
        assert set(got) == set(want), f"{variant} and {parent} have drifted apart"


class TestFiltering:
    def test_terms_in_one_facet_are_or_ed(self):
        """Otherwise `["dock", "predict-structure"]` — two actions — can only return nothing."""
        assert tool_tags.matches(["dock", "protein"], ["dock", "predict-structure"])

    def test_terms_in_different_facets_are_and_ed(self):
        assert not tool_tags.matches(["dock", "protein"], ["dock", "covalent"])
        assert tool_tags.matches(["dock", "covalent"], ["dock", "covalent"])

    def test_exclusion_beats_a_match(self):
        assert not tool_tags.matches(["protein", "small-molecule"], ["protein"],
                                     exclude=["small-molecule"])

    def test_an_empty_query_matches_everything(self):
        """`exclude=["fetch"]` alone must be a valid question: what runs with no network."""
        assert tool_tags.matches(["protein"], [])
        assert not tool_tags.matches(["fetch"], [], exclude=["fetch"])

    def test_the_documented_negative_query_returns_prodigy(self):
        found = tool_tags.select(tool_docs.collect(), ["protein", "binding"], ["small-molecule"])
        assert [e["name"] for e in found] == ["Prodigy"]

    def test_excluding_fetch_answers_what_runs_offline(self):
        offline = {e["name"] for e in tool_tags.select(tool_docs.collect(), [], ["fetch"])}
        assert "AlphaFold" in offline, "AlphaFold accepts a supplied MSA and is not tagged fetch"
        assert "RCSB" not in offline and "UniProt" not in offline

    @pytest.mark.parametrize("typed,meant", sorted(tool_tags.QUERY_ALIAS.items()))
    def test_every_alias_points_at_a_real_tag(self, typed, meant):
        assert meant in tool_tags.ALL_TAGS, f"{typed} resolves to {meant}, which is not a tag"
        assert typed not in tool_tags.ALL_TAGS, f"{typed} is itself a tag; the alias shadows it"

    def test_an_alias_is_reported_not_silently_applied(self):
        resolved, applied = tool_tags.resolve_query(["affinity", "protein"])
        assert resolved == ["binding", "protein"]
        assert applied == {"affinity": "binding"}


def test_tags_reach_the_generated_index():
    """The index is what an agent greps; tags that stop at the docs help nobody."""
    rendered = tool_docs.render(tool_docs.collect())
    assert "tags: predict-structure" in rendered
    assert "covalent" in rendered
