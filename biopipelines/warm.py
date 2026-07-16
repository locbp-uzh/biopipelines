# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""One-shot environment warm-up for BioPipelines.

`Tool.install()` only registers an install step *inside* an active Pipeline
context (see base_config.py). A bare `Boltz2.install()` at module scope is a
silent no-op. This wrapper opens a Pipeline, sets Resources, and calls
`.install()` for each named tool, so a container/cluster user can pre-build the
per-tool micromamba envs and download weights ONCE onto a persistent mount:

    bp-warm Boltz2 ProteinMPNN
    bp-warm --gpu A100 --time 2:00:00 Boltz2

All installs are idempotent; re-running skips tools already present. Selection
of config variant / OTF / output routing is via the usual env vars
(BIOPIPELINES_CONFIG_VARIANT, BIOPIPELINES_OTF, BIOPIPELINES_LOCAL_OUTPUT).
"""

import argparse
import sys

import biopipelines as bp
from biopipelines import Pipeline, Resources


def main():
    parser = argparse.ArgumentParser(
        prog="bp-warm",
        description="Pre-install BioPipelines tool environments and weights "
                    "(runs each Tool.install() inside a Pipeline context).",
    )
    parser.add_argument("tools", nargs="+",
                        help="Tool names to install, e.g. Boltz2 ProteinMPNN RFdiffusion")
    parser.add_argument("--gpu", default="A100",
                        help="GPU type for the install Resources (default: A100)")
    parser.add_argument("--memory", default="32GB", help="Memory (default: 32GB)")
    parser.add_argument("--time", default="2:00:00", help="Wallclock (default: 2:00:00)")
    parser.add_argument("--cpus", type=int, default=8, help="CPUs (default: 8)")
    parser.add_argument("--force-reinstall", action="store_true",
                        help="Reinstall even if the tool is already present")
    args = parser.parse_args()

    missing = [t for t in args.tools if not hasattr(bp, t)]
    if missing:
        print(f"ERROR: unknown tool(s): {', '.join(missing)}", file=sys.stderr)
        print("Run `bp-config list` or see docs/tool_index.md for valid names.", file=sys.stderr)
        sys.exit(1)

    with Pipeline("BioPipelines", "warm_up",
                  description="Environment warm-up (install-only)"):
        Resources(gpu=args.gpu, memory=args.memory, time=args.time, cpus=args.cpus)
        for name in args.tools:
            getattr(bp, name).install(force_reinstall=args.force_reinstall)
            print(f"queued install: {name}")
    print("warm-up complete.")


if __name__ == "__main__":
    main()
