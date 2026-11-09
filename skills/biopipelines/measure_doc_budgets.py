#!/usr/bin/env python3
"""Measure the token budgets the agent prompts quote, so the figures can be refreshed from a measurement.

`skills/biopipelines/SKILL.md` and `llm/pipelines.md` tell an assistant to read the tool index instead of the
prose docs, and they justify that with token costs. Those costs go stale every time the catalog grows. Run this
after adding tools, and copy the numbers it prints into the two files.

Estimator: `tiktoken` (cl100k_base) when it is installed, otherwise characters / 4, the usual approximation for
English prose with markdown punctuation. The character estimate is within a few percent of cl100k_base for this
repo's documentation, and the printed output always names the estimator used.

Tool counts per category file are derived from the generated index, not from the prose: an index entry ends with
`docs/tool/<file>.md#<anchor>`, so counting those pointers per file is the same count the index itself presents.
"""

import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
INDEX = os.path.join("skills", "biopipelines", "references", "tool_index.md")
CHARS_PER_TOKEN = 4


def count_tokens(text):
    """(tokens, estimator name) for one string."""
    try:
        import tiktoken
    except ImportError:
        return len(text) // CHARS_PER_TOKEN, f"characters / {CHARS_PER_TOKEN} (tiktoken is not installed)"
    encoding = tiktoken.get_encoding("cl100k_base")
    return len(encoding.encode(text)), "tiktoken cl100k_base"


def read(relative_path):
    with open(os.path.join(REPO, relative_path), encoding="utf-8") as handle:
        return handle.read()


def tools_per_category(index_text):
    """{category file name: entries pointing at it}, read off the index's own pointers."""
    counts = {}
    for name in re.findall(r"`docs/tool/([a-z_]+\.md)#", index_text):
        counts[name] = counts.get(name, 0) + 1
    return counts


def thousands(tokens):
    """A token count as the prose writes it."""
    return f"~{tokens / 1000:.0f}k" if tokens >= 10000 else f"~{tokens / 1000:.1f}k"


def main():
    index_text = read(INDEX)
    per_category = tools_per_category(index_text)

    targets = [(INDEX, index_text), ("README.md", read("README.md"))]
    category_files = sorted(
        name for name in os.listdir(os.path.join(REPO, "docs", "tool")) if name.endswith(".md")
    )
    prose_total = 0
    for name in category_files:
        text = read(os.path.join("docs", "tool", name))
        targets.append((os.path.join("docs", "tool", name), text))
        prose_total += len(text)

    estimator = None
    rows = []
    for path, text in targets:
        tokens, estimator = count_tokens(text)
        rows.append((path.replace("\\", "/"), len(text), tokens, per_category.get(os.path.basename(path), "")))

    width = max(len(row[0]) for row in rows)
    print(f"Estimator: {estimator}\n")
    print(f"{'file'.ljust(width)}  {'bytes':>9}  {'tokens':>9}  tools")
    print(f"{'-' * width}  {'-' * 9}  {'-' * 9}  -----")
    for path, size, tokens, tools in rows:
        print(f"{path.ljust(width)}  {size:>9,}  {tokens:>9,}  {tools}")

    prose_tokens, _ = count_tokens("".join(read(os.path.join("docs", "tool", n)) for n in category_files))
    index_tokens = rows[0][2]
    readme_tokens = rows[1][2]
    print(f"\n{'docs/tool/*.md (all)'.ljust(width)}  {prose_total:>9,}  {prose_tokens:>9,}  "
          f"{sum(per_category.values())}")

    print("\nFigures for SKILL.md and llm/pipelines.md:")
    print(f"  tool index          {thousands(index_tokens)} tokens")
    print(f"  docs/tool/*.md      {thousands(prose_tokens)} tokens")
    print(f"  README.md           {thousands(readme_tokens)} tokens")
    analysis_tokens, _ = count_tokens(read(os.path.join("docs", "tool", "analysis.md")))
    print(f"  analysis.md         {thousands(analysis_tokens)} tokens for "
          f"{per_category.get('analysis.md', 0)} tools")
    return 0


if __name__ == "__main__":
    sys.exit(main())
