#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""PLIP runner. For each complex PDB, emits two CSVs (per-interaction
long table and per-complex summary) and optionally one .pse PyMOL session
per binding site.

Three analysis modes, selected by --mode:
  * ligand (default): protein-ligand interactions (auto-detected HETATM
    ligands, optionally restricted by --ligand-json codes).
  * peptide: protein-peptide / protein-protein. The --chains chain(s) are
    treated as the peptide/partner ligand (PLIP --peptides / config.PEPTIDES).
  * intra: intra-chain interactions within a single --chains chain
    (PLIP --intra / config.INTRA).

Two backends are supported, selected by --container-prefix:
  * empty -> use the PLIP Python API (`plip.structure.preparation.PDBComplex`)
    inside the active conda env; the env must provide the plip package.
  * set   -> shell `<prefix> -f <pdb> -o <out> -x [-y]` against the
    PharmAI PLIP apptainer image, then parse report.xml. The prefix must
    end in a trailing space (it's exactly what BaseConfig.container_prefix()
    returns)."""

import argparse
import glob
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from biopipelines.biopipelines_io import load_datastream, iterate_files, read_upstream_missing, step_id_from_table_path  # noqa: E402
from biopipelines.ligand_utils import resolve_ligand_codes  # noqa: E402


INT_COLS = ["id", "ligand", "interaction_type", "residue", "chain", "resnum", "distance", "details"]
SUM_COLS = ["id", "n_hbonds", "n_hydrophobic", "n_pi_stacking", "n_salt_bridges", "n_halogen", "n_water_bridges", "session_files"]


# ---------------------------------------------------------------------------
# Python-API backend (conda env path)
# ---------------------------------------------------------------------------

def collect_api(interactions, complex_id, ligand_label):
    """Flatten PLIP interaction objects from the Python API."""
    rows = []
    for hb in interactions.hbonds_pdon + interactions.hbonds_ldon:
        rows.append({
            "id": complex_id, "ligand": ligand_label,
            "interaction_type": "hbond",
            "residue": hb.restype, "chain": hb.reschain, "resnum": hb.resnr,
            "distance": round(hb.distance_ad, 3),
            "details": f"donor={hb.dtype} acceptor={hb.atype}",
        })
    for hc in interactions.hydrophobic_contacts:
        rows.append({
            "id": complex_id, "ligand": ligand_label,
            "interaction_type": "hydrophobic",
            "residue": hc.restype, "chain": hc.reschain, "resnum": hc.resnr,
            "distance": round(hc.distance, 3),
            "details": "",
        })
    for ps in interactions.pistacking:
        rows.append({
            "id": complex_id, "ligand": ligand_label,
            "interaction_type": "pi_stacking",
            "residue": ps.restype, "chain": ps.reschain, "resnum": ps.resnr,
            "distance": round(ps.distance, 3),
            "details": f"type={ps.type}",
        })
    for sb in interactions.saltbridge_lneg + interactions.saltbridge_pneg:
        rows.append({
            "id": complex_id, "ligand": ligand_label,
            "interaction_type": "salt_bridge",
            "residue": sb.restype, "chain": sb.reschain, "resnum": sb.resnr,
            "distance": round(sb.distance, 3),
            "details": "",
        })
    for hl in interactions.halogen_bonds:
        rows.append({
            "id": complex_id, "ligand": ligand_label,
            "interaction_type": "halogen",
            "residue": hl.restype, "chain": hl.reschain, "resnum": hl.resnr,
            "distance": round(hl.distance_ax, 3),
            "details": "",
        })
    for wb in interactions.water_bridges:
        rows.append({
            "id": complex_id, "ligand": ligand_label,
            "interaction_type": "water_bridge",
            "residue": wb.restype, "chain": wb.reschain, "resnum": wb.resnr,
            "distance": round(wb.distance_aw, 3),
            "details": "",
        })
    return rows


def _set_api_config(mode, chains):
    """Set PLIP's module-global analysis mode before PDBComplex.analyze().

    config.PEPTIDES (list of chains) drives --peptides/--inter; config.INTRA
    (single chain) drives --intra. Reset both each call so per-id state cannot
    leak across the loop."""
    from plip.basic import config as plip_config

    plip_config.PEPTIDES = list(chains) if mode == "peptide" else []
    plip_config.INTRA = chains[0] if mode == "intra" else None


def analyze_api(sid, pdb_path, filt, mode, chains):
    """Run PLIP via its Python API. Returns (int_rows, totals, has_sites)."""
    from plip.structure.preparation import PDBComplex

    _set_api_config(mode, chains)
    mol = PDBComplex()
    mol.load_pdb(pdb_path)
    mol.analyze()
    if not mol.interaction_sets:
        return [], _empty_totals(), False

    int_rows = []
    totals = _empty_totals()
    for site_id, site in mol.interaction_sets.items():
        resname = site_id.split(":")[0]
        if filt and resname not in filt:
            continue
        int_rows.extend(collect_api(site, sid, site_id))
        totals["n_hbonds"] += len(site.hbonds_pdon) + len(site.hbonds_ldon)
        totals["n_hydrophobic"] += len(site.hydrophobic_contacts)
        totals["n_pi_stacking"] += len(site.pistacking)
        totals["n_salt_bridges"] += len(site.saltbridge_lneg) + len(site.saltbridge_pneg)
        totals["n_halogen"] += len(site.halogen_bonds)
        totals["n_water_bridges"] += len(site.water_bridges)
    return int_rows, totals, True


def write_pses_api(pdb_path, sid, sessions_dir, mode, chains):
    """Shell `plip -y` on the host. Conda-env path."""
    os.makedirs(sessions_dir, exist_ok=True)
    with tempfile.TemporaryDirectory() as work:
        cmd = ["plip", "-f", pdb_path, "-o", work, "-y"] + _mode_cli_flags(mode, chains)
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"plip -y failed: {res.stderr.strip() or res.stdout.strip()}")
        return _stage_pses(work, sid, sessions_dir)


# ---------------------------------------------------------------------------
# Container backend (apptainer CLI + XML parser)
# ---------------------------------------------------------------------------

# Map PLIP XML <interaction-type> tags to (output_label, distance_field, details_fn)
XML_INTERACTIONS = {
    "hydrogen_bond":           ("hbond",        "dist_h-a",
                                lambda e: "donor=" + (e.findtext("donortype") or "?") +
                                          " acceptor=" + (e.findtext("acceptortype") or "?")),
    "hydrophobic_interaction": ("hydrophobic",  "dist",      lambda e: ""),
    "pi_stacking":             ("pi_stacking",  "centdist",
                                lambda e: "type=" + (e.findtext("type") or "?")),
    "pi_cation_interaction":   ("pi_cation",    "dist",      lambda e: ""),
    "salt_bridge":             ("salt_bridge",  "dist",      lambda e: ""),
    "halogen_bond":            ("halogen",      "dist",      lambda e: ""),
    "water_bridge":            ("water_bridge", "dist_a-w",  lambda e: ""),
    "metal_complex":           ("metal",        "dist",      lambda e: ""),
}

# Map XML section -> summary counter
XML_SUMMARY = {
    "hydrogen_bond":           "n_hbonds",
    "hydrophobic_interaction": "n_hydrophobic",
    "pi_stacking":             "n_pi_stacking",
    "salt_bridge":             "n_salt_bridges",
    "halogen_bond":            "n_halogen",
    "water_bridge":            "n_water_bridges",
}


def _empty_totals():
    return {"n_hbonds": 0, "n_hydrophobic": 0, "n_pi_stacking": 0,
            "n_salt_bridges": 0, "n_halogen": 0, "n_water_bridges": 0}


def _parse_distance(elem, field):
    val = elem.findtext(field)
    if val in (None, ""):
        return None
    try:
        return round(float(val), 3)
    except ValueError:
        return None


def parse_plip_xml(xml_path, sid, filt):
    """Parse PLIP report.xml into (int_rows, totals, has_sites)."""
    tree = ET.parse(xml_path)
    root = tree.getroot()
    sites = root.findall(".//bindingsite")
    if not sites:
        return [], _empty_totals(), False

    int_rows = []
    totals = _empty_totals()
    for site in sites:
        hetid = site.findtext(".//hetid", "") or "UNK"
        if filt and hetid.upper() not in filt:
            continue
        chain = site.findtext(".//chain", "") or ""
        position = site.findtext(".//position", "") or ""
        ligand_label = f"{hetid}:{chain}:{position}"

        for tag, (label, dist_field, details_fn) in XML_INTERACTIONS.items():
            for elem in site.findall(f".//{tag}"):
                int_rows.append({
                    "id": sid, "ligand": ligand_label,
                    "interaction_type": label,
                    "residue": elem.findtext("restype", ""),
                    "chain": elem.findtext("reschain", ""),
                    "resnum": elem.findtext("resnr", ""),
                    "distance": _parse_distance(elem, dist_field),
                    "details": details_fn(elem),
                })
                counter = XML_SUMMARY.get(tag)
                if counter:
                    totals[counter] += 1
    return int_rows, totals, True


def _mode_cli_flags(mode, chains):
    """PLIP CLI flags for the analysis mode (shared by container + PSE shell)."""
    if mode == "peptide":
        return ["--peptides"] + list(chains)
    if mode == "intra":
        return ["--intra", chains[0]]
    return []


def run_container(prefix_tokens, pdb_path, out_dir, generate_pse, mode, chains):
    """Invoke `<prefix_tokens> -f <pdb_path> -o <out_dir> -x [-y]` and wait."""
    cmd = prefix_tokens + ["-f", pdb_path, "-o", out_dir, "-x"] + _mode_cli_flags(mode, chains)
    if generate_pse:
        cmd.append("-y")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"plip container failed (exit {res.returncode}): "
                           f"{res.stderr.strip() or res.stdout.strip()}")
    return res


def _stage_pses(src_dir, sid, sessions_dir):
    """Move .pse files from src_dir into sessions_dir, prefixed by sid."""
    pses = sorted(glob.glob(os.path.join(src_dir, "*.pse")))
    out_paths = []
    for i, src in enumerate(pses):
        suffix = os.path.basename(src).replace(".pse", "")
        dst_name = f"{sid}.pse" if i == 0 else f"{sid}__{suffix}.pse"
        dst = os.path.join(sessions_dir, dst_name)
        shutil.copyfile(src, dst)
        out_paths.append(dst)
    return out_paths


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--structures-json", required=True)
    p.add_argument("--interactions-csv", required=True)
    p.add_argument("--summary-csv", required=True)
    p.add_argument("--missing-csv", required=True)
    p.add_argument("--ligand-json", default="",
                   help="compounds-stream JSON; restrict to its `code`(s)")
    p.add_argument("--mode", default="ligand", choices=["ligand", "peptide", "intra"])
    p.add_argument("--chains", nargs="*", default=[],
                   help="peptide/partner chain(s) (mode=peptide) or the single "
                        "chain (mode=intra)")
    p.add_argument("--generate-pse", action="store_true")
    p.add_argument("--sessions-dir", default="")
    p.add_argument("--sessions-map-csv", default="")
    p.add_argument("--container-prefix", default="",
                   help="When set, run PLIP via this shell prefix (e.g. "
                        "'apptainer exec -B ... <image> ') and parse report.xml. "
                        "Otherwise use the in-env Python API.")
    p.add_argument("--upstream-missing", nargs="*", default=None)
    args = p.parse_args()

    # Restrict to the ligand codes carried by the compounds stream, if given.
    filt = {c.upper() for c in resolve_ligand_codes(args.ligand_json)} if args.ligand_json else None
    prefix_tokens = shlex.split(args.container_prefix) if args.container_prefix.strip() else []
    use_container = bool(prefix_tokens)

    ds = load_datastream(args.structures_json)
    int_rows, sum_rows, sess_rows, missing_rows = [], [], [], []
    step_id = step_id_from_table_path(args.missing_csv)

    for sid, pdb_path in iterate_files(ds):
        try:
            if use_container:
                with tempfile.TemporaryDirectory() as work:
                    run_container(prefix_tokens, pdb_path, work, args.generate_pse,
                                  args.mode, args.chains)
                    xml_candidates = sorted(glob.glob(os.path.join(work, "*report.xml")))
                    if not xml_candidates:
                        missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": "container produced no report.xml"})
                        continue
                    rows, totals, has_sites = parse_plip_xml(xml_candidates[0], sid, filt)
                    if not has_sites:
                        print(f"  {sid}: no ligand sites detected", file=sys.stderr)
                        sum_rows.append({"id": sid, **_empty_totals(), "session_files": ""})
                        if args.generate_pse:
                            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": "no ligand sites detected"})
                        continue
                    int_rows.extend(rows)
                    session_paths = (_stage_pses(work, sid, args.sessions_dir)
                                     if (args.generate_pse and args.sessions_dir) else [])
            else:
                rows, totals, has_sites = analyze_api(sid, pdb_path, filt, args.mode, args.chains)
                if not has_sites:
                    print(f"  {sid}: no ligand sites detected", file=sys.stderr)
                    sum_rows.append({"id": sid, **_empty_totals(), "session_files": ""})
                    if args.generate_pse:
                        missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": "no ligand sites detected"})
                    continue
                int_rows.extend(rows)
                session_paths = []
                if args.generate_pse and args.sessions_dir:
                    try:
                        session_paths = write_pses_api(pdb_path, sid, args.sessions_dir,
                                                       args.mode, args.chains)
                    except Exception as e:
                        print(f"WARNING: {sid} PSE generation failed: {e}", file=sys.stderr)
                        missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": f"PSE generation failed: {str(e)[:160]}"})
                        session_paths = []
        except Exception as e:
            print(f"WARNING: {sid} PLIP analysis failed: {e}", file=sys.stderr)
            missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": str(e)[:200]})
            continue

        if args.generate_pse and args.sessions_dir:
            if session_paths:
                sess_rows.append({"id": sid, "file": session_paths[0]})
            else:
                missing_rows.append({"id": sid, "removed_by": step_id, "kind": "failure", "cause": "plip produced no .pse"})
        sum_rows.append({"id": sid, **totals, "session_files": ";".join(session_paths)})
        print(f"  {sid}: {sum(totals.values())} total interactions" +
              (f" + {len(session_paths)} PSE" if session_paths else ""))

    all_missing = read_upstream_missing(args.upstream_missing) + missing_rows

    for d in (args.interactions_csv, args.summary_csv, args.missing_csv):
        os.makedirs(os.path.dirname(d), exist_ok=True)
    pd.DataFrame(int_rows, columns=INT_COLS).to_csv(args.interactions_csv, index=False)
    pd.DataFrame(sum_rows, columns=SUM_COLS).to_csv(args.summary_csv, index=False)
    pd.DataFrame(all_missing, columns=["id", "removed_by", "kind", "cause"]).to_csv(args.missing_csv, index=False)
    print(f"Interactions: {args.interactions_csv} ({len(int_rows)} rows)")
    print(f"Summary: {args.summary_csv} ({len(sum_rows)} rows)")
    print(f"Missing: {args.missing_csv} ({len(all_missing)} rows)")
    if args.generate_pse and args.sessions_map_csv:
        os.makedirs(os.path.dirname(args.sessions_map_csv), exist_ok=True)
        pd.DataFrame(sess_rows, columns=["id", "file"]).to_csv(args.sessions_map_csv, index=False)
        print(f"Sessions: {args.sessions_map_csv} ({len(sess_rows)} rows)")

    if missing_rows:
        print(f"Failed: {len(missing_rows)}/{len(missing_rows)+len(sum_rows)}", file=sys.stderr)
    if not sum_rows:
        sys.exit(1)


if __name__ == "__main__":
    main()
