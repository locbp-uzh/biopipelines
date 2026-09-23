#!/usr/bin/env python3
"""Write the agent-facing tool index that ships with the skill.

All the parsing lives in `biopipelines/tool_docs.py`, which the MCP server reads through too,
so the index and `bp_tools` can never disagree. This script only decides where the file goes.

Run from the repo root: `python skills/biopipelines/build_tool_index.py`
"""

import io
import pathlib
import sys

# Running a script puts its own directory on sys.path, not the cwd, so an uninstalled
# checkout cannot import the package without this.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent.parent))

from biopipelines.tool_docs import ROOT, collect, render  # noqa: E402

OUT = ROOT / "skills" / "biopipelines" / "references" / "tool_index.md"


def main():
    entries = collect()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    io.open(OUT, "w", encoding="utf-8").write(render(entries))
    size = OUT.stat().st_size
    print(f"{len(entries)} tools -> {OUT.relative_to(ROOT)} ({size} bytes, ~{size // 4} tokens)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
