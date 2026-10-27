# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Helpers shared by PPIDdG and LigandDdG: PyRosetta initialization, the bound-minus-unbound energy, replicate statistics and the per-id result loop (runtime, used by the pipe scripts), plus the install script and argument validators the two wrappers share (configuration time).

Imported by both sides, so it has no intra-package imports and never imports pyrosetta at module level. The protocol steps follow the Chaperone blocks `ppi_ddg`, `ligand_ddg` and `ligand_ddg_soft` (RosettaCommons/chaperone, MIT license).

PyRosetta is free for academic and non-profit use and needs a paid commercial license from the University of Washington. Installing it accepts that license.
"""

import math
import os
import re
import sys
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

# Copied from Chaperone's blocks.utils.pyrosetta.init_pyrosetta; used by PPIDdG and by LigandDdG(scoring="soft").
RELAX_INIT_OPTIONS = (
    "-run:preserve_header",
    "-use_input_sc",
    "-fast_restyping",
    "-ignore_waters true",
    "-beta",
    "-corrections::beta_nov16",
    "-skip_connect_info",
    "-mute all",
)

# The scoring-relevant lines of Chaperone's flags_rescore_genpot (the output-only ones are dropped); used by LigandDdG(scoring="hard").
HARD_INIT_OPTIONS = (
    "-use_input_sc",
    "-no_his_his_pairE",
    "-in:file:fullatom",
    "-no_optH false",
    "-flip_HNQ",
    "-beta_nov16",
    "-dunbrack_prob_buried 0.8",
    "-dunbrack_prob_nonburied 0.8",
    "-dunbrack_prob_buried_semi 0.8",
    "-dunbrack_prob_nonburied_semi 0.8",
)

UNBOUND_DISTANCE_A = 10000.0
PACK_RADIUS_A = 10.0

# Where the install script puts Rosetta's public python scripts, relative to the environment prefix (sys.prefix); the wheel does not ship molfile_to_params.py. Forward slashes: the install bash is generated on any OS.
ROSETTA_PUBLIC_DIR = "share/rosetta_public"
ROSETTA_PUBLIC_SPARSE_PATH = "source/scripts/python/public"
MOLFILE_TO_PARAMS = f"{ROSETTA_PUBLIC_DIR}/{ROSETTA_PUBLIC_SPARSE_PATH}/molfile_to_params.py"

_INITIALIZED = False


def molfile_to_params_path() -> str:
    """Absolute path of molfile_to_params.py inside the active environment."""
    return os.path.join(sys.prefix, *MOLFILE_TO_PARAMS.split("/"))


def init_pyrosetta(options: Sequence[str]) -> None:
    """Initialize PyRosetta; options bind at the first call and later calls are no-ops, since Rosetta initializes once per process."""
    global _INITIALIZED
    if _INITIALIZED:
        return
    import pyrosetta
    pyrosetta.init(" ".join(options), silent=True)
    _INITIALIZED = True


def move_chain_far_away(pose, chain: str, distance: float = UNBOUND_DISTANCE_A):
    """A clone of the pose with every residue of `chain` translated `distance` A along +x."""
    import pyrosetta as pyr
    ros = pyr.rosetta
    moved = pose.clone()
    subset = ros.core.select.residue_selector.ChainSelector(chain).apply(moved)
    x_unit = ros.numeric.xyzVector_double_t(1, 0, 0)
    shift = ros.numeric.xyzVector_double_t(distance, 0, 0)
    ros.protocols.toolbox.pose_manipulation.rigid_body_move(x_unit, 0, shift, moved, subset)
    return moved


def unbound_ddg(pose, scorefxn, move_chain: str) -> float:
    """score(bound) - score(`move_chain` translated far away), with no repacking."""
    pose = pose.clone()
    bound = scorefxn(pose)
    unbound = scorefxn(move_chain_far_away(pose, move_chain))
    return bound - unbound


def repack_around_chain(pose, lock_chain: str, scorefxn, radius: float = PACK_RADIUS_A) -> None:
    """Repack, in place, the residues within `radius` of `lock_chain`; the chain itself and everything farther stay fixed."""
    import pyrosetta as pyr
    ros = pyr.rosetta
    select = ros.core.select.residue_selector
    task = ros.core.pack.task
    lock = select.ChainSelector(lock_chain)
    nearby = select.NeighborhoodResidueSelector(lock, radius, False)
    frozen = select.OrResidueSelector(lock, select.NotResidueSelector(nearby))
    factory = task.TaskFactory()
    factory.push_back(task.operation.RestrictToRepacking())
    factory.push_back(task.operation.OperateOnResidueSubset(task.operation.PreventRepackingRLT(), frozen, False))
    factory.push_back(task.operation.IncludeCurrent())
    packer = ros.protocols.minimization_packing.PackRotamersMover()
    packer.score_function(scorefxn)
    packer.task_factory(factory)
    packer.apply(pose)


def mean_std(values: Sequence[float]) -> Tuple[float, float]:
    """Mean and population standard deviation (ddof = 0, as numpy's default)."""
    n = len(values)
    mean = sum(values) / n
    return mean, math.sqrt(sum((v - mean) ** 2 for v in values) / n)


def trimmed_mean_std(values: Sequence[float]) -> Tuple[float, float, int]:
    """Mean, population std and count after dropping the best and the worst when there are at least three values."""
    ordered = sorted(values)
    kept = ordered[1:-1] if len(ordered) >= 3 else ordered
    mean, std = mean_std(kept)
    return mean, std, len(kept)


def is_missing(value) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value))


def collect_rows(items: Iterable[Tuple[str, str]], score: Callable[[str, str], Dict],
                 columns: Sequence[str], step_id: str) -> Tuple[List[Dict], List[Dict], int, int]:
    """Score every (id, path); one result row per id that entered, NaN where scoring raised.

    Returns (rows, local_missing_rows, attempted, failures). A failure leaves the id in the result table with NaN values and adds a kind="failure" row naming the cause.
    """
    rows, missing = [], []
    attempted = failures = 0
    for sid, path in items:
        attempted += 1
        try:
            result = score(sid, path)
            rows.append({"id": sid, **{c: result.get(c, float("nan")) for c in columns if c != "id"}})
        except Exception as e:
            print(f"WARNING: {sid} failed: {e}", file=sys.stderr)
            failures += 1
            rows.append({"id": sid, **{c: float("nan") for c in columns if c != "id"}})
            missing.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})
    return rows, missing, attempted, failures


def write_outputs(table_csv: str, rows: List[Dict], columns: Sequence[str],
                  local_missing_csv: str, missing_rows: List[Dict]) -> None:
    import pandas as pd
    for path in (table_csv, local_missing_csv):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(rows, columns=list(columns)).to_csv(table_csv, index=False)
    pd.DataFrame(missing_rows, columns=["id", "removed_by", "kind", "cause"]).to_csv(local_missing_csv, index=False)


def step_exit_status(attempted: int, failures: int) -> Optional[str]:
    """The error message when the whole step failed, None when it did not.

    Nothing attempted means the upstream stream was empty or its files absent; every attempt raising means the tool itself is broken. NaN rows alone never fail the step.
    """
    if not attempted:
        return "no input entities to process (upstream stream is empty or its files are absent)"
    if failures == attempted:
        return f"all {attempted} attempted id(s) failed"
    return None


# ---------------------------------------------------------------------------
# Wrapper-side helpers (configuration time)
# ---------------------------------------------------------------------------

_CHAIN_RE = re.compile(r"^[A-Za-z0-9]$")


def validate_chain(field: str, value) -> None:
    if not isinstance(value, str) or not _CHAIN_RE.match(value):
        raise ValueError(f"{field} must be a single-character PDB chain id, got {value!r}")


def validate_positive_int(field: str, value) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{field} must be a positive integer, got {value!r}")


def pyrosetta_install_script(cls, folders, env_manager, force_reinstall, env_name=None) -> str:
    """Install the PyRosetta environment shared by PPIDdG and LigandDdG.

    Creates the environment from environments/PyRosetta.yaml + PyRosetta.pip.txt (the wheel comes from the RosettaCommons index, which also carries the license text), then fetches Rosetta's public python scripts with a sparse clone into <env prefix>/share/rosetta_public: the wheel does not ship molfile_to_params.py, which LigandDdG needs to parameterize a ligand.
    """
    env = cls._install_env(env_manager, env_name)
    biopipelines = folders.get("biopipelines", "")
    run = cls._env_run(env, env_manager)

    remove_block = cls._env_remove_block(env, env_manager) if force_reinstall else ""
    env_block = cls._env_install_block(env, env_manager, biopipelines)
    install_block = f"{remove_block}\n{env_block}"
    if not force_reinstall:
        # An existing env is still verified below, so a half-built one is not reported as installed.
        install_block = f"""if {cls._env_exists_check(env, env_manager)}; then
    echo "{env} env already present, verifying."
else
{install_block}
fi
"""
    verify = "import pyrosetta, rdkit, dimorphite_dl, pandas, yaml"
    return f"""echo "=== Installing PyRosetta ({env}) ==="
echo "PyRosetta is free for academic use and needs a paid commercial license; installing it accepts the Rosetta license."
{install_block}
PREFIX="$({run}python -c "import sys; print(sys.prefix)")"
ROSETTA_PUBLIC="$PREFIX/{ROSETTA_PUBLIC_DIR}"
if [ ! -d "$ROSETTA_PUBLIC/.git" ]; then
    mkdir -p "$(dirname "$ROSETTA_PUBLIC")"
    git clone --depth 1 --filter=blob:none --sparse https://github.com/RosettaCommons/rosetta.git "$ROSETTA_PUBLIC"
fi
git -C "$ROSETTA_PUBLIC" sparse-checkout set {ROSETTA_PUBLIC_SPARSE_PATH}

if {run}python -c "{verify}" >/dev/null 2>&1 && [ -f "$PREFIX/{MOLFILE_TO_PARAMS}" ]; then
    touch "$INSTALL_SUCCESS"
    echo "=== PyRosetta installation complete ==="
else
    echo "ERROR: PyRosetta verification failed (import of {verify.split(' ', 1)[1]} or molfile_to_params.py missing)"
    {run}python -c "{verify}" 2>&1 | tail -5
    exit 1
fi
"""
