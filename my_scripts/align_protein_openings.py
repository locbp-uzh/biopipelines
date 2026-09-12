# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""Face two protein openings across a plane and emit the posed complex, ready for linker design.

An *opening* is the pair of backbone atoms a linker has to bridge: for a scaffold with a domain cut out, the C of the residue before the gap and the N of the residue after it; for the domain going in, the N of its first residue and the C of its last. Each opening plus its own centroid defines a triangle. This step puts both triangles in the z=0 plane, rotates the second body so the two openings face each other, and separates them along the shared bisector by ``d``.

That fixes five of the six rigid-body degrees of freedom from geometry alone. The rest are scanned: ``d``, the tilt of each body about its own triangle normal (``a1``, ``a2``), and the roll of the mobile body about its N-C axis (``b2``). Among poses with no more than ``clash_max`` heavy-atom contacts below ``clash_hard``, the winner is the one whose *longer* linker is shortest, ties broken by how many atom pairs sit inside ``pack_cutoff``.

Minimizing the longer linker rather than the sum matters: the sum is happily paid down by making one linker very short and the other enormous, which is the opposite of what a sensor wants. ``objective="sum"`` restores the old behavior.

``d`` starts negative because the two bodies should end up packed against each other, not merely reachable. In jGCaMP8f the transducer touches the scaffold at 2.45 A with 60 heavy-atom pairs inside 4 A and only 2 inside 2.6 A -- which is what ``clash_max`` is calibrated to.

Ported from ``segments_align`` in the group's PyMOL scripts so it runs on a compute node — the PyMOL build there is missing libGL.

Parameters come from a table rather than arguments, because Scripting passes no free-form scalars and because per-design scan ranges are more useful than one global setting. Every column is optional and falls back to the default below; a single-row table is broadcast to every structure.

    Scripting("align_protein_openings.py", env="ProteinEnv",
              inputs={"structures": designs, "template": scaffold, "params": Table("scan.csv")})

``env="ProteinEnv"`` is required, not a preference: the neighbour search is ``scipy.spatial.cKDTree`` and the default ``biopipelines`` env carries no scipy, so the step raises ``ModuleNotFoundError`` after the upstream folds have already been paid for. ProteinEnv has scipy 1.13.0 and, since conda-forge ``libgl`` was added to it, imports on a CPU compute node — so the step still asks for ``gpu="none"``.

``scan.csv`` columns: ``id, d, a1, a2, b2, keep_residues, anchor_c, anchor_n, clash_hard, clash_max, keep_het, keep_het_mobile``.
``objective`` is ``balanced`` (default) or ``sum``. ``keep_het`` lists the HETATM residue names to retain (``CA`` for the calciums, ``*`` for all but water, empty for none); a residue-number range cannot do this job, since additives sit outside the protein range and a chromophore belonging to the excised domain sits inside it.
The scan columns accept a single number (``12``) or a comma-separated list (``8,12,16,20``).

Standalone, outside a pipeline:

    python align_protein_openings.py --template 7st4.pdb --designs 'ORDER_14/structures/*.pdb' \\
        --keep 8-29,278-422 --anchor-c 29 --anchor-n 278 --out poses
"""

DEFAULTS = {"d": "-8,-4,0,4,8,12", "a1": "-40,-30,-20,-10,0,10,20,30,40",
            "a2": "-40,-30,-20,-10,0,10,20,30,40", "b2": "0,45,90,135,180,225,270,315",
            "clash_hard": 2.6, "clash_max": 2, "pack_cutoff": 4.0, "objective": "balanced",
            "keep_het": "", "keep_het_mobile": "*"}


def _numbers(spec, fallback):
    text = str(spec if spec not in (None, "") else fallback)
    return [float(v) for v in text.replace(";", ",").split(",") if v.strip() != ""]


def _ranges(spec):
    """'8-29,278-422' -> [8..29, 278..422]"""
    out = []
    for part in str(spec).replace(";", ",").split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part.lstrip("-"):
            lo, hi = part.split("-")
            out.extend(range(int(lo), int(hi) + 1))
        else:
            out.append(int(part))
    return out


def _het_set(spec, fallback):
    """'' keeps no HETATM, '*' keeps all but water, 'CA,ZN' keeps those residue names."""
    text = fallback if spec in (None, "") and spec is not fallback else spec
    text = "" if text is None else str(text)
    if text.strip() == "*":
        return None
    return {t.strip() for t in text.replace(";", ",").split(",") if t.strip()}


def _geometry():
    """Everything numpy-dependent, imported inside so config time stays light."""
    import numpy as np

    def rot(axis, deg):
        a = np.asarray(axis, float)
        n = np.linalg.norm(a)
        if n < 1e-9:
            return np.eye(3)
        a, t = a / n, np.radians(deg)
        K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
        return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * (K @ K)

    def planify(n, c):
        """Matrix putting n on +x and c in the z=0 plane. Body must already be centered."""
        M, n2, c2 = np.eye(3), n.copy(), c.copy()
        for axis, vec, idx in ((( 0, 0, 1), "n", 1), ((0, 1, 0), "n", 2), ((1, 0, 0), "c", 2)):
            v = n2 if vec == "n" else c2
            if abs(v[idx]) < 1e-6:
                continue
            if axis == (0, 0, 1):
                R = rot(axis, -np.degrees(np.arctan2(n2[1], n2[0])))
            elif axis == (0, 1, 0):
                R = rot(axis, np.degrees(np.arctan2(n2[2], n2[0])))
            else:
                R = rot(axis, -np.degrees(np.arctan2(c2[2], c2[1])))
            M, n2, c2 = R @ M, R @ n2, R @ c2
            if axis == (0, 0, 1) and n2[0] < 0:
                R = rot(axis, 180)
                M, n2, c2 = R @ M, R @ n2, R @ c2
        return M

    return np, rot, planify


class _Body:
    """A rigid set of atoms plus the two backbone atoms its linkers attach to."""

    def __init__(self, records, coords, n_xyz, c_xyz):
        self.records, self.X, self.n, self.c = records, coords.copy(), n_xyz.copy(), c_xyz.copy()

    @property
    def F(self):
        return self.X.mean(0)

    def moved(self, M=None, t=None, origin=None):
        import numpy as np
        b = _Body(self.records, self.X, self.n, self.c)
        o = np.zeros(3) if origin is None else origin
        if M is not None:
            b.X = (b.X - o) @ M.T + o
            b.n = (b.n - o) @ M.T + o
            b.c = (b.c - o) @ M.T + o
        if t is not None:
            b.X, b.n, b.c = b.X + t, b.n + t, b.c + t
        return b


def _read_pdb(path, keep=None, keep_het=None, drop_first_residue=False):
    """Returns (records, coords, ordered protein resi) where a record is (kind, atom_name, resname, resi).

    keep_het is a set of HETATM residue names to retain, or None for all but water. A residue number
    filter cannot police HETATMs: cofactors and crystallization additives are numbered outside the
    protein range, and a chromophore belonging to a domain being cut out is numbered inside it.
    """
    import numpy as np
    rows, seen = [], []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[16] not in " A":
            continue
        resi, het = int(line[22:26]), line.startswith("HETATM")
        resn = line[17:20].strip()
        if het:
            if resn == "HOH":
                continue
            if keep_het is not None and resn not in keep_het:
                continue
        elif keep is not None and resi not in keep:
            continue
        elem = line[76:78].strip() or line[12:16].strip().lstrip("0123456789")[:1]
        rows.append((("HETATM" if het else "ATOM  "), line[12:16], line[17:20], resi, elem,
                     float(line[30:38]), float(line[38:46]), float(line[46:54])))
        if not het and resi not in seen:
            seen.append(resi)
    if drop_first_residue and seen:
        rows = [r for r in rows if not (r[0] == "ATOM  " and r[3] == seen[0])]
        seen = seen[1:]
    recs = [r[:5] for r in rows]
    return recs, np.array([r[5:] for r in rows], float), seen


def _atom(recs, X, resi, name):
    for r, x in zip(recs, X):
        if r[3] == resi and r[1].strip() == name:
            return x.copy()
    raise KeyError(f"atom {name} of residue {resi} not found")


def _face(scaf, mob, d):
    """Both openings into z=0, faced, separated by d."""
    np, rot, planify = _geometry()
    s = scaf.moved(t=-scaf.F)
    g = mob.moved(t=-mob.F)
    s = s.moved(M=planify(s.n, s.c))
    g = g.moved(M=planify(g.n, g.c))
    ang = lambda v: np.degrees(np.arctan2(v[1], v[0]))
    g = g.moved(M=rot([0, 0, 1], 180 + 0.5 * ang(s.c) - 0.5 * ang(g.c)))
    u = 0.5 * (s.n + s.c)
    u = u / np.linalg.norm(u)
    l1 = 0.5 * (abs(s.n @ u) + abs(s.c @ u))
    l2 = 0.5 * (abs(g.n @ u) + abs(g.c @ u))
    return s.moved(t=-(l1 + 0.5 * d) * u), g.moved(t=(l2 + 0.5 * d) * u)


def _scan(scaf, mob, p):
    """Best pose over the d/a1/a2/b2 grid. Returns (row, scaffold, mobile) or None."""
    np, rot, _ = _geometry()
    from scipy.spatial import cKDTree
    hard = float(p.get("clash_hard") or DEFAULTS["clash_hard"])
    cmax = int(float(p.get("clash_max") or DEFAULTS["clash_max"]))
    pack = float(p.get("pack_cutoff") or DEFAULTS["pack_cutoff"])
    balanced = str(p.get("objective") or DEFAULTS["objective"]).strip() != "sum"
    best = None
    for d in _numbers(p.get("d"), DEFAULTS["d"]):
        s0, g0 = _face(scaf, mob, d)
        n1 = np.cross(s0.n - s0.F, s0.c - s0.F)
        n2 = np.cross(g0.n - g0.F, g0.c - g0.F)
        for a1 in _numbers(p.get("a1"), DEFAULTS["a1"]):
            s = s0.moved(M=rot(n1, a1), origin=s0.F)
            tree = cKDTree(s.X)
            for a2 in _numbers(p.get("a2"), DEFAULTS["a2"]):
                gt = g0.moved(M=rot(n2, a2), origin=g0.F)
                for b2 in _numbers(p.get("b2"), DEFAULTS["b2"]):
                    g = gt.moved(M=rot(gt.n - gt.c, b2), origin=0.5 * (gt.n + gt.c))
                    dist = tree.query(g.X, workers=-1)[0]
                    nclash = int((dist < hard).sum())
                    if nclash > cmax:
                        continue
                    r1 = float(np.linalg.norm(g.n - s.c))
                    r2 = float(np.linalg.norm(s.n - g.c))
                    packed = int((dist < pack).sum())
                    # the longer linker is what a sensor pays for, so minimize the max rather than the sum: minimizing the sum buys a short L2 by making L1 enormous
                    key = (round(max(r1, r2), 1), -packed) if balanced else (round(r1 + r2, 1), -packed)
                    if best is None or key < best[0]:
                        best = (key,
                                dict(d=d, a1=a1, a2=a2, b2=b2, clashes=nclash, packed=packed,
                                     reach_L1=round(r1, 2), reach_L2=round(r2, 2)), s, g)
    return None if best is None else (best[1], best[2], best[3])


def _write(path, scaf, mob):
    with open(path, "w") as fh:
        i = 0
        for body, chain in ((scaf, "A"), (mob, "B")):
            for (kind, name, resn, resi, elem), xyz in zip(body.records, body.X):
                i += 1
                # columns 77-78 must carry the element: RFdiffusion's parser indexes l[77] unguarded
                fh.write(f"{kind}{i:5d} {name} {resn:>3s} {chain}{resi:4d}    "
                         f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00"
                         f"{'':10s}{elem:>2s}\n")
            fh.write("TER\n")
        fh.write("END\n")


def configuration(inputs):
    from biopipelines.scripting_api import Stream, Table
    ids = inputs["structures"].ids
    return {"structures": Stream("pdb", ids),
            "poses": Table(columns=["id", "d", "a1", "a2", "b2", "clashes", "packed",
                                    "reach_L1", "reach_L2", "L1_min", "L2_min"])}


def execution(inputs, outputs):
    import math

    params = inputs.get("params")
    template = next(iter(inputs["template"].iterate()))[1]

    def row_for(sid):
        if params is None:
            return {}
        r = params.row(sid)
        if r is None:
            rows = params.rows()
            r = rows[0] if len(rows) == 1 else {}
        return dict(r or {})

    for sid, path in inputs["structures"].iterate():
        p = row_for(sid)
        keep = _ranges(p["keep_residues"]) if p.get("keep_residues") else None
        anchor_c, anchor_n = int(p.get("anchor_c") or 0), int(p.get("anchor_n") or 0)
        het_t = _het_set(p.get("keep_het"), DEFAULTS["keep_het"])
        het_m = _het_set(p.get("keep_het_mobile"), DEFAULTS["keep_het_mobile"])

        srecs, sX, sres = _read_pdb(template, keep=keep, keep_het=het_t)
        scaf = _Body(srecs, sX, _atom(srecs, sX, anchor_n or sres[-1], "N"),
                     _atom(srecs, sX, anchor_c or sres[0], "C"))
        mrecs, mX, mres = _read_pdb(path, keep_het=het_m, drop_first_residue=True)
        mob = _Body(mrecs, mX, _atom(mrecs, mX, mres[0], "N"), _atom(mrecs, mX, mres[-1], "C"))

        hit = _scan(scaf, mob, p)
        if hit is None:
            outputs.drop(sid, cause="no pose below the clash ceiling")
            continue
        info, s, g = hit
        _write(outputs["structures"].file(sid, f"{sid}.pdb"), s, g)
        # rise measured on jGCaMP8f's own linkers: 3.01 A per residue for L1, 2.54 for L2
        info.update(id=sid,
                    L1_min=max(1, math.ceil(info["reach_L1"] / 3.01) - 1),
                    L2_min=max(1, math.ceil(info["reach_L2"] / 2.54) - 1))
        outputs["poses"].row(info)


if __name__ == "__main__":
    import argparse, csv, glob, math, os

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--template", required=True)
    ap.add_argument("--designs", required=True, help="glob for the mobile bodies")
    ap.add_argument("--keep", default=None, help="template residues to keep, e.g. 8-29,278-422")
    ap.add_argument("--anchor-c", type=int, default=0, help="template residue whose C the linker leaves")
    ap.add_argument("--anchor-n", type=int, default=0, help="template residue whose N the linker enters")
    for k, v in DEFAULTS.items():
        ap.add_argument(f"--{k.replace('_', '-')}", default=v)
    ap.add_argument("--out", default="poses")
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    p = {k: getattr(a, k) for k in DEFAULTS}
    keep = _ranges(a.keep) if a.keep else None
    srecs, sX, sres = _read_pdb(a.template, keep=keep, keep_het=_het_set(a.keep_het, DEFAULTS["keep_het"]))
    scaf = _Body(srecs, sX, _atom(srecs, sX, a.anchor_n or sres[-1], "N"),
                 _atom(srecs, sX, a.anchor_c or sres[0], "C"))

    rows = []
    for f in sorted(glob.glob(a.designs)):
        sid = os.path.splitext(os.path.basename(f))[0]
        mrecs, mX, mres = _read_pdb(f, keep_het=_het_set(a.keep_het_mobile, DEFAULTS["keep_het_mobile"]),
                                    drop_first_residue=True)
        mob = _Body(mrecs, mX, _atom(mrecs, mX, mres[0], "N"), _atom(mrecs, mX, mres[-1], "C"))
        hit = _scan(scaf, mob, p)
        if hit is None:
            print(f"{sid}: no pose below the clash ceiling")
            continue
        info, s, g = hit
        _write(os.path.join(a.out, f"{sid}.pdb"), s, g)
        info.update(id=sid,
                    L1_min=max(1, math.ceil(info["reach_L1"] / 3.01) - 1),
                    L2_min=max(1, math.ceil(info["reach_L2"] / 2.54) - 1))
        rows.append(info)
        print(f"{sid}: d={info['d']:.0f} a=({info['a1']:.0f},{info['a2']:.0f}) b={info['b2']:.0f} "
              f"clash={info['clashes']} packed={info['packed']:3d} "
              f"reach {info['reach_L1']:5.1f}/{info['reach_L2']:5.1f} A "
              f"-> L1>={info['L1_min']} L2>={info['L2_min']}")

    if rows:
        with open(os.path.join(a.out, "poses.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    print(f"\n{len(rows)} poses -> {a.out}/poses.csv")
