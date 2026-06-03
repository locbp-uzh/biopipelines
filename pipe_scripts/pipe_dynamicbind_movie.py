#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""DynamicBind movie renderer (runs under ProteinEnv / PyMOL).

Reads the movie-jobs manifest written by pipe_dynamicbind_postprocess.py
(columns: id, receptor_pdb, ligand_sdf). For each kept pose it loads the
relaxed receptor + docked ligand into headless PyMOL, renders a spinning
turntable as PNG frames, and encodes one <id>.mp4 per pose (falling back to a
.gif when ffmpeg is unavailable). Writes a movies map_table (id, file)."""

import argparse
import csv
import glob
import os
import sys

# Headless PyMOL must be configured before import (matches pipe_pymol.py).
os.environ["PYMOL_SYMMETRY_VIEWER"] = "0"
os.environ["DISPLAY"] = ""
os.environ["QT_QPA_PLATFORM"] = "offscreen"

import pymol  # noqa: E402
pymol.pymol_argv = ["pymol", "-cqp"]
pymol.finish_launching(["pymol", "-cqp"])
from pymol import cmd, movie  # noqa: E402

N_FRAMES = 60  # one full 360-degree turn


def render_pose(out_id, receptor_pdb, ligand_sdf, frames_dir, out_mp4):
    """Render a spinning movie of one pose. Returns the produced file path or ''."""
    cmd.reinitialize()
    cmd.load(receptor_pdb, "receptor")
    cmd.load(ligand_sdf, "ligand")
    cmd.hide("everything")
    cmd.show("cartoon", "receptor")
    cmd.color("gray80", "receptor")
    cmd.show("sticks", "ligand")
    cmd.util.cbag("ligand")
    cmd.orient("ligand")
    cmd.zoom("ligand", 8)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)

    # Build a 360-degree turntable across N_FRAMES.
    cmd.mset(f"1 x{N_FRAMES}")
    cmd.util.mroll(1, N_FRAMES, 1)

    os.makedirs(frames_dir, exist_ok=True)
    frame_prefix = os.path.join(frames_dir, f"{out_id}_")
    cmd.set("ray_trace_frames", 0)  # plain render is far faster than ray per frame
    cmd.mpng(frame_prefix, width=640, height=480)

    # Prefer an mp4 via PyMOL's movie.produce (uses ffmpeg when present).
    try:
        movie.produce(out_mp4, mode="mpeg", quality=80, quiet=1)
        if os.path.exists(out_mp4) and os.path.getsize(out_mp4) > 0:
            return out_mp4
    except Exception as exc:  # ffmpeg missing or encode failure
        print(f"  [warn] {out_id}: mp4 encode failed ({exc}); trying gif", file=sys.stderr)

    # Fallback: stitch the PNG frames into an animated gif with Pillow.
    frames = sorted(glob.glob(f"{frame_prefix}*.png"))
    if not frames:
        return ""
    try:
        from PIL import Image
        imgs = [Image.open(f).convert("RGB") for f in frames]
        out_gif = os.path.splitext(out_mp4)[0] + ".gif"
        imgs[0].save(out_gif, save_all=True, append_images=imgs[1:],
                     duration=60, loop=0)
        return out_gif if os.path.exists(out_gif) else ""
    except Exception as exc:
        print(f"  [warn] {out_id}: gif fallback failed ({exc})", file=sys.stderr)
        return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--movie-jobs-csv", required=True)
    ap.add_argument("--movies-folder", required=True)
    ap.add_argument("--movies-map", required=True)
    args = ap.parse_args()

    os.makedirs(args.movies_folder, exist_ok=True)
    os.makedirs(os.path.dirname(args.movies_map), exist_ok=True)
    frames_dir = os.path.join(args.movies_folder, "_frames")

    if not os.path.exists(args.movie_jobs_csv):
        print(f"ERROR: movie jobs manifest missing: {args.movie_jobs_csv}", file=sys.stderr)
        sys.exit(1)

    map_rows, failed = [], []
    with open(args.movie_jobs_csv, newline="") as f:
        jobs = list(csv.DictReader(f))

    for job in jobs:
        out_id = job["id"]
        receptor_pdb = job["receptor_pdb"]
        ligand_sdf = job["ligand_sdf"]
        out_mp4 = os.path.join(args.movies_folder, f"{out_id}.mp4")
        try:
            produced = render_pose(out_id, receptor_pdb, ligand_sdf, frames_dir, out_mp4)
            if produced:
                map_rows.append({"id": out_id, "file": produced})
                print(f"  {out_id}: {produced}")
            else:
                failed.append(out_id)
        except Exception as exc:
            print(f"WARNING: {out_id} movie failed: {exc}", file=sys.stderr)
            failed.append(out_id)

    with open(args.movies_map, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "file"])
        w.writeheader()
        w.writerows(map_rows)
    print(f"Movies: {args.movies_map} ({len(map_rows)} rendered, {len(failed)} failed)")


if __name__ == "__main__":
    main()
