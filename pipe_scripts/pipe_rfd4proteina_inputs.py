#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Build the RFD4-Proteina condition spec at execution time.

`rfd4-generate` reads one YAML whose top-level keys are design names and whose values hold `input` (a CIF), `contig` and `conditions`. The input paths are only known at run time, so this script writes the final spec from the config-time design entry and the input structures stream: one entry per input id, each pointing at a CIF. A CIF is used as it is; any other format (PDB) is converted with AtomWorks first, because upstream rejects PDB input.

Each prepared input is also built once, the way upstream builds it, so a contig that does not fit that input (a residue range the chain does not have) is recorded as a failure of that id with upstream's own message, instead of stopping the run at that example and losing every design after it.

It writes a manifest of which input ids were attempted, built and failed, which the post-processing step reads.

Usage:
    python pipe_rfd4proteina_inputs.py \\
        --config-yaml <config.yaml> --structures-json <ds.json> --inputs-dir <dir> \\
        --spec-yaml <condition_spec.yaml> --manifest-json <inputs_manifest.json>
"""

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Dict

import yaml

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files  # noqa: E402

# The structure suffixes upstream's `is_structure_file` accepts without conversion.
CIF_SUFFIXES = (".cif", ".cif.gz", ".cif.zst", ".bcif", ".bcif.gz", ".bcif.zst")


def is_cif_like(path: str) -> bool:
    return str(path).lower().endswith(CIF_SUFFIXES)


def valid_design_id(design_id: str) -> bool:
    """An id names one output directory upstream, so it cannot be blank, '.', '..' or carry a separator."""
    return design_id not in ("", ".", "..") and not any(c in design_id for c in ("/", "\\", "\0"))


def build_spec(template: Dict, inputs: Dict[str, str]) -> Dict[str, Dict]:
    """One design entry per input id: `input` first, then the shared contig and conditions."""
    return {sid: {"input": path, **template} for sid, path in inputs.items()}


def convert_to_cif(structure_path: str, out_path: str) -> None:
    """PDB to CIF with AtomWorks `parse` + `to_cif_file`, the conversion upstream's FAQ prescribes.

    Uses the DISTILLATION preset, which upstream built for PDB-format inputs (it does not use the altloc mode a PDB cannot carry); the condition spec re-parses the CIF with its own preset afterwards. The preset loads standard annotations, which AtomWorks refuses for a PDB (a PDB carries none), so that one flag is turned off.
    """
    from atomworks.io import parse
    from atomworks.io.utils.io_utils import to_cif_file
    from rfproteina.utils.parse_presets import DISTILLATION

    result = parse(structure_path, config=DISTILLATION.replace(build_assembly=None, load_standard_annotations=False))
    to_cif_file(result["asym_unit"][0], out_path)


def validate_example(sid: str, entry: Dict, spec_dir: str) -> None:
    """Build one design problem the way upstream does, so a contig that does not fit this input fails here, named, and not in the middle of the whole run.

    `rfd4-generate` builds its examples inside a data loader and stops at the first one that raises, which loses every design after it. The entry goes through YAML first so it carries the types upstream reads back.
    """
    try:
        from rfproteina.inference_utils.condition_spec import build_annotated_atom_array_from_condition_spec
    except ImportError:
        print("WARNING: rfproteina is not importable here, so the inputs are not checked before the run", file=sys.stderr)
        return
    try:
        build_annotated_atom_array_from_condition_spec(
            yaml.safe_load(yaml.safe_dump(entry)), spec_dir=Path(spec_dir), example_id=sid)
    except Exception as e:
        raise ValueError(" ".join(str(e).split())) from e


def main():
    ap = argparse.ArgumentParser(description="Build the RFD4-Proteina condition spec")
    ap.add_argument("--config-yaml", required=True, help="Config-time design entry (spec) and num_designs")
    ap.add_argument("--structures-json", required=True, help="Input structures DataStream JSON")
    ap.add_argument("--inputs-dir", required=True, help="Where converted CIFs are written")
    ap.add_argument("--spec-yaml", required=True, help="Condition spec to write")
    ap.add_argument("--manifest-json", required=True, help="Manifest of attempted, built and failed input ids")
    args = ap.parse_args()

    with open(args.config_yaml) as f:
        config = yaml.safe_load(f) or {}
    template = config["spec"]

    os.makedirs(args.inputs_dir, exist_ok=True)
    ds = load_datastream(args.structures_json)

    attempted, inputs, failed = [], {}, {}
    for sid, path in iterate_files(ds):
        attempted.append(sid)
        try:
            if not valid_design_id(sid):
                raise ValueError(f"{sid!r} cannot name an output directory")
            if is_cif_like(path):
                prepared = os.path.abspath(path)
            else:
                prepared = os.path.join(args.inputs_dir, f"{sid}.cif")
                convert_to_cif(path, prepared)
            validate_example(sid, {"input": prepared, **template}, os.path.dirname(os.path.abspath(args.spec_yaml)))
            inputs[sid] = prepared
        except Exception as e:
            print(f"WARNING: {sid} could not be prepared: {e}", file=sys.stderr)
            failed[sid] = str(e)[:200]

    for path in (args.spec_yaml, args.manifest_json):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(args.manifest_json, "w") as f:
        json.dump({"attempted": attempted, "built": list(inputs), "failed": failed}, f, indent=2)
    if inputs:
        with open(args.spec_yaml, "w") as f:
            yaml.safe_dump(build_spec(template, inputs), f, sort_keys=False, width=1000)

    print(f"Condition spec: {args.spec_yaml} ({len(inputs)} of {len(attempted)} input(s) prepared)")

    # Nothing attempted means the upstream stream was empty or its files are absent; a pipeline must not pass on that.
    if not attempted:
        print("ERROR: no input structures to process (stream is empty or its files are absent)", file=sys.stderr)
        sys.exit(1)
    if failed and len(failed) == len(attempted):
        print(f"ERROR: all {len(attempted)} input(s) failed", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
