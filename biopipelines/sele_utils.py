# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Shared selection-string utilities.

Three formats coexist in the pipeline:

* **Chain-aware** – ``A1-50+B10``  (used by RFdiffusion TRB tables)
* **Chainless**   – ``1-50+10``    (used by ProteinMPNN jsonl output, pLDDT paths)
* **Legacy space-separated** – ``1-50 10``  (old tables, still accepted on input)

All public functions live here so that HelpScripts can import them without
duplicating logic.
"""


def chain_aware_sele(residues):
    """Convert list of (chain, resnum) tuples to compact chain-aware selection string.

    Examples:
        [('A',1),('A',2),('A',3),('B',10)] -> 'A1-3+B10'
        [('A',1),('A',3),('B',5)]          -> 'A1+A3+B5'
    """
    if not residues:
        return ""

    # Group by chain, preserving order of first appearance
    chain_order = []
    by_chain = {}
    for chain, resnum in residues:
        if chain not in by_chain:
            chain_order.append(chain)
            by_chain[chain] = []
        by_chain[chain].append(resnum)

    parts = []
    for chain in chain_order:
        nums = sorted(by_chain[chain])
        i = 0
        while i < len(nums):
            start = nums[i]
            j = i + 1
            while j < len(nums) and nums[j] == nums[j - 1] + 1:
                j += 1
            end = nums[j - 1]
            if end > start:
                parts.append(f"{chain}{start}-{end}")
            else:
                parts.append(f"{chain}{start}")
            i = j

    return "+".join(parts)


def sele_to_list(s):
    """Parse a selection string into a list of (chain, resnum) tuples.

    Accepts both chain-aware (``A1-50+B10``) and legacy chainless (``1-50+10``)
    formats, as well as the old space-separated format (``1-50 10``).

    Returns:
        Sorted list of (chain, resnum) tuples.  Chain is ``''`` for chainless input.
    """
    result = []
    if not s or str(s) == "nan":
        return result

    s = str(s)

    # Split on '+' or whitespace (legacy)
    if '+' in s:
        parts = s.split('+')
    else:
        parts = s.split()

    for part in parts:
        part = part.strip()
        if not part:
            continue

        # Detect leading chain letter: uppercase letter followed by a digit
        chain = ''
        body = part
        if len(part) >= 2 and part[0].isalpha() and part[0].isupper() and part[1].isdigit():
            chain = part[0]
            body = part[1:]

        if '-' in body and not body.startswith('-'):
            range_parts = body.split('-')
            if len(range_parts) == 2:
                start, end = int(range_parts[0]), int(range_parts[1])
                for r in range(start, end + 1):
                    result.append((chain, r))
            else:
                print(f"Warning: Malformed range '{part}', skipping")
        else:
            try:
                result.append((chain, int(body)))
            except ValueError:
                print(f"Warning: Could not parse '{part}', skipping")

    return sorted(result, key=lambda x: (x[0], x[1]))


def list_to_sele(a):
    """Convert a list of residues to a selection string.

    Accepts three input forms (may be mixed):
        * plain ints:             ``[1, 2, 3, 7]``
        * ``(chain, resnum)`` tuples: ``[('A',1), ('A',2), ('B',10)]``
        * chain-prefixed strings: ``['A1', 'A2', 'B10']``

    When chain information is present the output is chain-aware
    (``'A1-3+B10'``); otherwise chainless (``'1-3+7'``).

    Examples:
        [1,2,3,7]                          -> '1-3+7'
        [('A',1),('A',2),('A',3),('B',10)] -> 'A1-3+B10'
        ['A1','A2','A3','B10']             -> 'A1-3+B10'
    """
    if not a:
        return ""

    # Normalise every element to (chain, resnum)
    normalised = []
    for item in a:
        if isinstance(item, (tuple, list)):
            normalised.append((str(item[0]), int(item[1])))
        elif isinstance(item, str):
            if len(item) >= 2 and item[0].isalpha() and item[0].isupper() and item[1:].lstrip('-').isdigit():
                normalised.append((item[0], int(item[1:])))
            else:
                normalised.append(('', int(item)))
        else:
            normalised.append(('', int(item)))

    # Group by chain, preserving order of first appearance
    chain_order = []
    by_chain = {}
    for chain, resnum in normalised:
        if chain not in by_chain:
            chain_order.append(chain)
            by_chain[chain] = []
        by_chain[chain].append(resnum)

    parts = []
    for chain in chain_order:
        nums = sorted(set(by_chain[chain]))
        i = 0
        while i < len(nums):
            start = nums[i]
            j = i + 1
            while j < len(nums) and nums[j] == nums[j - 1] + 1:
                j += 1
            end = nums[j - 1]
            if end > start:
                parts.append(f"{chain}{start}-{end}")
            else:
                parts.append(f"{chain}{start}")
            i = j

    return "+".join(parts)
