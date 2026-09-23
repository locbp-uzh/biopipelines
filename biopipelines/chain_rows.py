# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""The `chains=` contract shared by the inverse-folding wrappers.

An MPNN writes every chain of its backbone into one record. A `sequences` row is one polymer chain, so the wrapper has to say which chains become rows and what their ids are. Both questions are answered here so ProteinMPNN and LigandMPNN cannot drift apart on them.
"""

from typing import List, Optional, Union

try:
    from .combinatorics import generate_multiplied_ids_pattern
except ImportError:  # standalone import, see the dual-import note in developer_manual.md
    from combinatorics import generate_multiplied_ids_pattern

ChainsSpec = Union[str, List[str], None]


def normalize_chains(chains: ChainsSpec) -> ChainsSpec:
    """None, "all", or a list — a bare chain id is promoted to a one-element list.

    ``"auto"`` is the spelling the retired ``chain`` parameter used for "work it out from
    the structure", which is what ``None`` means here, so it maps onto it.
    """
    if not isinstance(chains, str):
        return chains
    if chains.lower() == "all":
        return "all"
    if chains.lower() == "auto":
        return None
    return [chains]


def validate_chains(chains: ChainsSpec) -> None:
    """Raise if `chains` is not None, "all", or a list of distinct single-character ids."""
    if chains is None or chains == "all":
        return
    if not isinstance(chains, list) or not chains:
        raise ValueError('chains must be None, "all", a chain id, or a list of chain ids')
    seen = set()
    for c in chains:
        if not isinstance(c, str) or len(c) != 1 or not c.isalnum():
            raise ValueError(f"chains entries must be single-character chain ids, got {c!r}")
        if c in seen:
            raise ValueError(f"chains lists {c!r} twice")
        seen.add(c)


def positions_chain(chains: ChainsSpec) -> str:
    """The chain a chainless position selection (``fixed="10-20"``) attaches to.

    One named chain answers it outright. Otherwise the runtime works it out from the
    structure, which is also where the ambiguous case is caught. Several named chains are
    that ambiguous case, so they must not answer with the first of them.
    """
    if isinstance(chains, list) and len(chains) == 1:
        return chains[0]
    return "auto"


def accepts_multiple(chains: ChainsSpec) -> bool:
    """Whether the caller has said this step may span several chains.

    Answers two questions at once, which is why one parameter can serve both: it is what
    makes a chain suffix necessary on the row ids, and what decides whether a structure
    with several chains is an error or expected. When the caller has not said so, several
    chains is an error rather than a silent pick of the first — that guess is the failure
    mode this module exists to remove.
    """
    return chains == "all" or (isinstance(chains, list) and len(chains) > 1)


def chains_arg(chains: ChainsSpec) -> str:
    """The `--chains` value `pipe_fa_to_csv_fasta.py` parses."""
    if chains == "all":
        return "all"
    return ",".join(chains or [])


def chain_row_ids(design_ids: List[str], chains: ChainsSpec) -> List[str]:
    """Chain-row ids for `design_ids` under `chains`.

    Only a setting that can emit more than one chain per design adds a suffix, so a single-chain run keeps the ids it had before `chains` existed. "all" is the one spelling whose chain letters are unknown until the backbone is read, hence the lazy bracket.
    """
    if chains == "all":
        return [f"{d}[_<?>]" for d in design_ids]
    if isinstance(chains, list) and len(chains) > 1:
        return generate_multiplied_ids_pattern(
            design_ids, f"<{' '.join(chains)}>", input_stream_name="designs"
        )
    return list(design_ids)
