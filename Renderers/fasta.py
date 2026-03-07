# Renderer for fasta/fa streams: inline colored sequence alignment

import os
import html as html_module


# ClustalX-style amino acid color groups
_AA_COLORS = {
    # Hydrophobic (blue)
    'A': '#80a0f0', 'V': '#80a0f0', 'I': '#80a0f0', 'L': '#80a0f0',
    'M': '#80a0f0', 'F': '#80a0f0', 'W': '#80a0f0', 'C': '#80a0f0',
    # Positive (red)
    'K': '#f01505', 'R': '#f01505', 'H': '#f01505',
    # Negative (magenta)
    'D': '#c048c0', 'E': '#c048c0',
    # Polar (green)
    'S': '#15c015', 'T': '#15c015', 'N': '#15c015', 'Q': '#15c015',
    # Glycine (orange)
    'G': '#f09048',
    # Proline (yellow)
    'P': '#c0c000',
    # Aromatic (cyan)
    'Y': '#15a8a8',
    # Nucleotides
    'U': '#f01505',
    # DNA/RNA
}
_NT_COLORS = {
    'A': '#64f073', 'T': '#f07268', 'U': '#f07268',
    'G': '#f0b264', 'C': '#73c2f5',
}
_GAP_COLOR = '#e8e8e8'
_DEFAULT_COLOR = '#dddddd'


def _parse_fasta(text):
    """Parse FASTA text into list of (id, sequence) tuples."""
    entries = []
    current_id = None
    parts = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_id is not None:
                entries.append((current_id, ''.join(parts)))
            current_id = line[1:].split()[0]
            parts = []
        else:
            parts.append(line)
    if current_id is not None:
        entries.append((current_id, ''.join(parts)))
    return entries


def _is_nucleotide(seq):
    return len(seq) > 0 and all(c in 'ACGTUacgtu-. ' for c in seq)


def render(stream, output):
    """Render fasta stream as inline colored sequence blocks."""
    # Find the first valid fasta file in the stream
    fasta_path = None
    if stream.files:
        for p in stream.files:
            if p and os.path.isfile(p):
                fasta_path = p
                break

    if fasta_path is None:
        return None

    try:
        with open(fasta_path, 'r') as f:
            text = f.read()
    except Exception:
        return None

    entries = _parse_fasta(text)
    if not entries:
        return None

    is_nt = _is_nucleotide(entries[0][1])
    color_map = _NT_COLORS if is_nt else _AA_COLORS

    MAX_SEQS = 50
    MAX_RES = 300
    truncated_seqs = len(entries) > MAX_SEQS
    entries = entries[:MAX_SEQS]
    max_id_len = max(len(e[0]) for e in entries)
    max_seq_len = max(len(e[1]) for e in entries)
    truncated_res = max_seq_len > MAX_RES

    parts = []
    parts.append(
        f'<div class="bp-section">'
        f'<div class="bp-section-title">{html_module.escape(stream.name)} '
        f'<span style="font-weight:normal;color:#666;">({stream.format}, '
        f'{len(stream)} item{"s" if len(stream) != 1 else ""})</span></div>'
    )

    # Monospace block, one row per sequence
    parts.append(
        '<div style="font-family:monospace;font-size:11px;line-height:1.4;'
        'overflow-x:auto;white-space:nowrap;border:1px solid #ddd;'
        'border-radius:4px;padding:6px;background:#fff;">'
    )

    for seq_id, seq in entries:
        display_seq = seq[:MAX_RES]
        # ID label
        padded_id = html_module.escape(seq_id.ljust(max_id_len))
        parts.append(
            f'<span style="color:#666;margin-right:8px;">{padded_id}</span>'
        )
        # Colored residues
        for ch in display_seq:
            upper = ch.upper()
            if upper in ('-', '.'):
                color = _GAP_COLOR
            else:
                color = color_map.get(upper, _DEFAULT_COLOR)
            parts.append(
                f'<span style="background:{color};padding:0 1px;">'
                f'{html_module.escape(ch)}</span>'
            )
        if truncated_res:
            parts.append('<span style="color:#aaa;"> …</span>')
        parts.append('<br>')

    if truncated_seqs:
        parts.append(
            f'<span style="color:#aaa;font-style:italic;">'
            f'… {len(entries)} of {stream._count if hasattr(stream, "_count") else "?"} sequences shown</span><br>'
        )

    parts.append('</div></div>')
    return '\n'.join(parts)
