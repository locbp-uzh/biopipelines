"""The tag vocabulary, and evidence for which tags a tool should carry.

An agent searching the index for "covalent" gets nothing today: the index carries one truncated
sentence per tool, and a capability like covalent linkage lives in a *parameter*, never in the
headline. Tags fix that at the source — in `docs/tool/*.md`, beside the tool they describe.

This module does two jobs with one set of rules:

* **Drafting.** Hand-writing ~76 tag lines invites a hasty pass, and a wrong tag is worse than a
  missing blurb: it turns an honest miss into a confident wrong answer. `evidence()` proposes
  tags from what each section actually says and shows the line it matched, so the first pass is
  a review rather than a blank page.
* **Linting, forever.** CI can check that every tag is in the vocabulary, but nothing can notice
  a *missing* one — except this: a section whose body talks about covalent bonds while its Tags
  line does not say `covalent` is a defect, and it stays detectable as the docs grow.

Evidence is weighted by where it appears. A keyword in a parameter name is strong: the tool
takes an argument for it. A keyword in prose or an example is weak: it may be describing
someone else's output, or a caveat about what the tool cannot do.
"""

import re

FACETS = ("action", "subject", "readout", "capability")

# The reconciled vocabulary. Kept closed on purpose: the discipline is that a tag earns its
# place only if it changes the answer to at least two distinct questions, and an open list
# drifts into `dock` / `docking` / `ligand-docking` within a month.
VOCABULARY = {
    "action": [
        "generate-backbone", "inverse-folding", "design-sequence", "predict-structure",
        "dock", "sample-ensemble", "refine-structure", "predict-property", "measure",
        "validate", "detect", "build-msa", "fetch", "data", "visualize",
    ],
    "subject": [
        "protein", "small-molecule", "nucleic-acid", "complex", "pocket", "residues",
        "msa", "ensemble",
    ],
    "readout": [
        "binding", "stability", "solubility", "fitness", "sasa", "interactions",
        "flexibility", "energy",
    ],
    "capability": [
        "covalent", "symmetry", "motif-scaffolding", "binder-design", "sequence-only",
        "all-atom", "sidechain",
    ],
}

ALL_TAGS = {tag for tags in VOCABULARY.values() for tag in tags}

# Variants that share a parent's implementation and must not drift apart, mirroring the
# BADGE_ALIAS that already exists for README badges.
TAG_ALIAS = {"SolubleMPNN": "ProteinMPNN", "LoadMultiple": "Load"}

# Scripting is the escape hatch: it can run anything, so every tag is true of it and none helps
# anyone find it. An empty tag set is the honest answer, and the lint exempts exactly this one.
UNTAGGED = {"Scripting"}

# Patterns are deliberately narrow for the capability facet, which is the facet that fixes the
# reported bug and the one where a false positive costs the most.
#
# Every key here must be a real tag: a pattern for a tag that no longer exists proposes
# something `unknown_tags()` then rejects, which wastes a reviewer's time. Not every tag needs a
# pattern — ACTION and SUBJECT come from the category seed, not from keywords.
PATTERNS = {
    "covalent": r"\bcovalent",
    "symmetry": r"\bsymmetr|\bcyclic\b",
    # `contigs?` bounded, not the bare prefix: `\bcontig` matched "contiguous segment".
    "motif-scaffolding": r"\bcontigs?\b|\bmotif[- ]scaffold|scaffold(ing)? a motif",
    "binder-design": r"\bbinder\b",
    "sequence-only": r"single[_ ]sequence|no MSA|without an MSA|sequence[- ]only|from sequence alone",
    "all-atom": r"all[- ]atom",
    "dock": r"\bdock(s|ing|ed)?\b",
    "build-msa": r"\bMSA\b.*\b(build|generate|search)|multiple sequence alignment",
    "binding": r"\baffinit|\bKd\b|\bpKd\b|binding free energy|\bddG\b.*bind",
    "stability": r"\bstabilit|\bddG\b|thermostab",
    "solubility": r"\bsolubilit|\baggregat",
    "pocket": r"\bpocket|binding site",
    "interactions": r"interaction fingerprint|\bhydrogen bond|contact map|\bcontacts\b",
    "sasa": r"\bSASA\b|solvent[- ]accessib",
    "flexibility": r"\bflexibilit|\bRMSF\b|conformational ensemble",
    "energy": r"\benergy\b|\benergies\b|minimi[sz]",
    "nucleic-acid": r"\bDNA\b|\bRNA\b|nucleic acid",
    "small-molecule": r"\bligand|small[- ]molecule|\bSMILES\b|\bcompound",
    "complex": r"\bcomplex\b|\binterface\b|protein[-–]protein|protein[-–]ligand",
    "ensemble": r"\bensemble\b|\bconformer",
}

# ACTION and SUBJECT are not keyword-detectable — "this is an inverse-folding tool for proteins"
# is never written in those words. They follow from the category the docs already assign, so the
# draft seeds them from there and the reviewer adjusts. CAPABILITY and READOUT are the facets
# keywords can genuinely propose, and the ones the blurb misses.
CATEGORY_SEED = {
    "Structure Generation": (["generate-backbone"], ["protein"]),
    "Sequence Design": (["inverse-folding"], ["protein"]),
    "Structure Prediction & Docking": (["predict-structure"], ["protein"]),
    "Analysis": (["measure"], ["protein"]),
    "Cheminformatics": (["measure"], ["small-molecule"]),
    "Sequence Statistics": (["predict-property"], ["protein"]),
    "MSAs": (["build-msa"], ["msa", "protein"]),
    "Data Management": (["data"], []),
    "Inputs & I/O": (["fetch"], []),
}

# `fetch` carries the runtime-network fact: these tools fail on an isolated compute node, so
# `exclude=["fetch"]` answers "what runs offline". A separate `network-required` tag was dropped
# because install-time weight downloads are near-universal and runtime network belongs to
# exactly this set. AlphaFold and Boltz2 are deliberately NOT fetch: querying an MSA server is
# incidental to what they do, and both accept a supplied MSA instead.
NOT_FETCH = {"AlphaFold", "Boltz2", "ESMFold", "ESMFold2"}


# Keyword evidence a reviewer has judged and rejected, with the reason. Without this list the
# lint is either noisy enough to be muted or quietly weakened; four named exceptions are cheaper
# than either. Every entry is a real sentence in the docs that uses the word in another sense.
LINT_EXEMPT = {
    ("ConformationalChange", "all-atom"):
        'atoms="all" picks which atoms the RMSD covers — a measurement mode, not all-atom modelling',
    ("LigandAtomSelector", "symmetry"):
        '"kept for API symmetry with DistanceSelector" is about the signature, not molecular symmetry',
    ("RCSB", "symmetry"):
        "SymmetryType is a search filter over curated PDB annotations, not a capability of the tool",
    ("LigandMPNN", "sequence-only"):
        '"emits sequence only" means no structure is written; sequence-only means runs without an MSA',
}


def seed_from_category(category):
    """The ACTION and SUBJECT tags a tool's category implies, before any review."""
    action, subject = CATEGORY_SEED.get(category, ([], []))
    return list(action), list(subject)


# Sections of a tool's doc entry, strongest evidence first.
ZONES = ("parameters", "streams", "tables", "prose", "example")
ZONE_WEIGHT = {"parameters": 3, "streams": 2, "tables": 2, "prose": 1, "example": 1}
STRONG = 3


def split_zones(section):
    """Break a tool's markdown section into the zones evidence is weighted by."""
    if not section:
        return {z: "" for z in ZONES}
    # The declared tags are not evidence for themselves: leaving the line in makes every tag
    # self-justifying and the missing-tag lint can never fire.
    section = re.sub(r"^\*\*Tags\*\*:.*$", "", section, flags=re.M)
    example_at = section.find("**Example**")
    example = section[example_at:] if example_at != -1 else ""
    body = section[:example_at] if example_at != -1 else section

    zones = {"example": example}
    for name, label in (("parameters", "**Parameters**"), ("streams", "**Streams**"),
                        ("tables", "**Tables**")):
        start = body.find(label)
        if start == -1:
            zones[name] = ""
            continue
        rest = body[start + len(label):]
        ends = [rest.find(other) for other in ("**Parameters**", "**Streams**", "**Tables**",
                                               "**References**", "**Environment**", "**WARNING")
                if rest.find(other) > 0]
        zones[name] = rest[:min(ends)] if ends else rest

    first_marker = min([p for p in (body.find("**Installation**"), body.find("**Environment**"),
                                    body.find("**Parameters**"), body.find("**References**"))
                        if p > 0] or [len(body)])
    zones["prose"] = body[:first_marker]
    return zones


def _excerpt(line, match, width=96):
    """The text around the match, not the start of the line.

    A parameter description runs long; quoting its first 110 characters usually shows the
    reviewer a sentence that does not contain the word being justified.
    """
    line = line.strip()
    if len(line) <= width:
        return line
    start = max(match.start() - width // 3, 0)
    end = min(start + width, len(line))
    return ("..." if start else "") + line[start:end] + ("..." if end < len(line) else "")


def evidence(section):
    """{tag: [(zone, weight, excerpt)]} for every tag the section gives a reason to carry."""
    zones = split_zones(section)
    found = {}
    for tag, pattern in PATTERNS.items():
        for zone in ZONES:
            hit = None
            for line in zones.get(zone, "").splitlines():
                m = re.search(pattern, line, re.I)
                if m:
                    hit = (zone, ZONE_WEIGHT[zone], _excerpt(line, m))
                    break
            if hit:
                found.setdefault(tag, []).append(hit)
    return found


def propose(section):
    """Tags with strong evidence, and tags worth a look — separated, never merged.

    The split is the point: a parameter named `covalent_linkage` is a fact about the tool, while
    the word "covalent" in a caveat may be saying the tool does *not* handle it.
    """
    found = evidence(section)
    strong, weak = [], []
    for tag, matches in found.items():
        (strong if max(w for _z, w, _l in matches) >= STRONG else weak).append(tag)
    return sorted(strong), sorted(weak)


def parse_tags_line(section):
    """The tags a section already declares, or None if it has no `**Tags**:` line."""
    match = re.search(r"^\*\*Tags\*\*:?\s*(.+)$", section or "", re.M | re.I)
    if not match:
        match = re.search(r"^\*\*Tags:\*\*\s*(.+)$", section or "", re.M | re.I)
    if not match:
        return None
    return [t.strip() for t in match.group(1).split(",") if t.strip()]


def unknown_tags(tags):
    return sorted(set(tags or []) - ALL_TAGS)


FACET_OF = {tag: facet for facet, tags in VOCABULARY.items() for tag in tags}


def matches(tool_tags_, wanted, exclude=()):
    """Does a tool's tag set satisfy the query?

    OR within a facet, AND across facets. `["dock", "predict-structure", "covalent"]` reads as
    *(docking or structure prediction) and covalent*, which is what someone typing three terms
    means; AND-ing all three would return nothing and OR-ing them would return the catalog.
    """
    have = set(tool_tags_ or [])
    if have & set(exclude or ()):
        return False
    by_facet = {}
    for tag in wanted or ():
        by_facet.setdefault(FACET_OF.get(tag), set()).add(tag)
    return all(have & group for group in by_facet.values())


def select(entries, tags=(), exclude=()):
    """The entries matching a tag query, in the order they were given."""
    return [e for e in entries if matches(e.get("tags"), tags, exclude)]


# Terms an agent reasonably types that are not the tag. Spelling suggestions cannot bridge these
# — nothing in "affinity" resembles "binding" — and each one was a real name during the design.
QUERY_ALIAS = {
    "affinity": "binding", "binding-affinity": "binding", "kd": "binding", "docking": "dock",
    "select-residues": "residues", "residue": "residues", "per-residue": "residues",
    "transform-data": "data", "table": "data", "network-required": "fetch", "download": "fetch",
    "structure-prediction": "predict-structure", "folding": "predict-structure",
    "sidechain-packing": "sidechain", "scaffolding": "motif-scaffolding",
    "aggregation": "solubility", "ddg": "stability", "rmsf": "flexibility",
}


def resolve_query(tags):
    """(resolved tags, {typed: meant}) — the aliases applied, so the answer can say what it read."""
    out, applied = [], {}
    for tag in tags or ():
        real = QUERY_ALIAS.get(tag.strip().lower(), tag)
        if real != tag:
            applied[tag] = real
        if real not in out:
            out.append(real)
    return out, applied


def close_tags(tag, limit=4):
    """Vocabulary entries nearest a term that is not one, for an error worth reading."""
    import difflib
    close = difflib.get_close_matches(tag, sorted(ALL_TAGS), n=limit, cutoff=0.4)
    lowered = tag.lower()
    return close + [t for t in sorted(ALL_TAGS) if lowered in t and t not in close][:limit]


def vocabulary_text():
    """The whole vocabulary, by facet — what an agent needs after getting a tag wrong."""
    return "\n".join(f"  {facet:<11} {', '.join(VOCABULARY[facet])}" for facet in FACETS)


def missing_capability_tags(section, declared, name=None):
    """Capability tags the section argues for but does not declare — what CI cannot otherwise see."""
    strong, _weak = propose(section)
    capability = set(VOCABULARY["capability"])
    gaps = (set(strong) & capability) - set(declared or [])
    return sorted(gaps - {tag for tool, tag in LINT_EXEMPT if tool == name})
