#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Runtime helper script for Angle analysis.

This script analyzes protein structures to calculate angles between atoms.
Three modes are supported (determined by `mode` in the config):

- "bond"   : (a, o, b)           Bond angle at o (0-180 degrees)
- "torsion": (a, x1, x2, b)      Torsional/dihedral angle (-180 to 180 degrees)
- "vector" : [[a1, a2], [b1, b2]] Angle between vectors a1->a2 and b1->b2 (0-180 degrees)

Parses PDB files directly as text, using the same selection syntax as Distance.
"""

import os
import sys
import argparse
import json
import math
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple

# Import unified I/O utilities
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from biopipelines.biopipelines_io import load_datastream, iterate_files
from biopipelines.pdb_parser import parse_pdb_file, Atom, resolve_selection


def get_atom_centroid(atoms: List[Atom]) -> Tuple[float, float, float]:
    """
    Calculate the centroid of a list of atoms.

    Args:
        atoms: List of Atom objects

    Returns:
        Tuple of (x, y, z) coordinates

    Raises:
        ValueError: If atoms list is empty
    """
    if not atoms:
        raise ValueError("Cannot calculate centroid of empty atom list")

    x = sum(a.x for a in atoms) / len(atoms)
    y = sum(a.y for a in atoms) / len(atoms)
    z = sum(a.z for a in atoms) / len(atoms)

    return (x, y, z)


def calculate_bond_angle(p1: Tuple[float, float, float],
                         p2: Tuple[float, float, float],
                         p3: Tuple[float, float, float]) -> float:
    """
    Calculate bond angle at p2 (angle between p1-p2-p3).

    Args:
        p1: Coordinates of first atom
        p2: Coordinates of middle atom (vertex)
        p3: Coordinates of third atom

    Returns:
        Angle in degrees (0-180)
    """
    # Vector from p2 to p1
    v1 = (p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2])
    # Vector from p2 to p3
    v2 = (p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2])

    # Dot product
    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

    # Magnitudes
    mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)

    if mag1 == 0 or mag2 == 0:
        raise ValueError("Cannot calculate angle: zero-length vector")

    # Clamp to avoid numerical issues with acos
    cos_angle = max(-1.0, min(1.0, dot / (mag1 * mag2)))
    angle_rad = math.acos(cos_angle)

    return math.degrees(angle_rad)


def calculate_torsion_angle(p1: Tuple[float, float, float],
                            p2: Tuple[float, float, float],
                            p3: Tuple[float, float, float],
                            p4: Tuple[float, float, float]) -> float:
    """
    Calculate torsional/dihedral angle between planes p1-p2-p3 and p2-p3-p4.

    Uses the standard IUPAC convention for dihedral angles.

    Args:
        p1: Coordinates of first atom
        p2: Coordinates of second atom
        p3: Coordinates of third atom
        p4: Coordinates of fourth atom

    Returns:
        Torsion angle in degrees (-180 to 180)
    """
    # Vector from p1 to p2
    b1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    # Vector from p2 to p3 (central bond)
    b2 = (p3[0] - p2[0], p3[1] - p2[1], p3[2] - p2[2])
    # Vector from p3 to p4
    b3 = (p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2])

    # Normal to plane 1 (b1 x b2)
    n1 = (
        b1[1]*b2[2] - b1[2]*b2[1],
        b1[2]*b2[0] - b1[0]*b2[2],
        b1[0]*b2[1] - b1[1]*b2[0]
    )

    # Normal to plane 2 (b2 x b3)
    n2 = (
        b2[1]*b3[2] - b2[2]*b3[1],
        b2[2]*b3[0] - b2[0]*b3[2],
        b2[0]*b3[1] - b2[1]*b3[0]
    )

    # m1 = n1 x b2/|b2| (for determining sign)
    b2_mag = math.sqrt(b2[0]**2 + b2[1]**2 + b2[2]**2)
    if b2_mag == 0:
        raise ValueError("Cannot calculate torsion: central bond has zero length")

    b2_unit = (b2[0]/b2_mag, b2[1]/b2_mag, b2[2]/b2_mag)
    m1 = (
        n1[1]*b2_unit[2] - n1[2]*b2_unit[1],
        n1[2]*b2_unit[0] - n1[0]*b2_unit[2],
        n1[0]*b2_unit[1] - n1[1]*b2_unit[0]
    )

    # x = n1 . n2
    x = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]
    # y = m1 . n2
    y = m1[0]*n2[0] + m1[1]*n2[1] + m1[2]*n2[2]

    return math.degrees(math.atan2(y, x))


def calculate_vector_angle(p1: Tuple[float, float, float],
                           p2: Tuple[float, float, float],
                           p3: Tuple[float, float, float],
                           p4: Tuple[float, float, float]) -> float:
    """
    Calculate the angle between vectors p1->p2 and p3->p4.

    Args:
        p1, p2: Start and end of the first vector (direction p1 to p2)
        p3, p4: Start and end of the second vector (direction p3 to p4)

    Returns:
        Angle in degrees (0-180)
    """
    v1 = (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])
    v2 = (p4[0] - p3[0], p4[1] - p3[1], p4[2] - p3[2])

    dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    mag1 = math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    mag2 = math.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)

    if mag1 == 0 or mag2 == 0:
        raise ValueError("Cannot calculate vector angle: zero-length vector")

    cos_angle = max(-1.0, min(1.0, dot / (mag1 * mag2)))
    return math.degrees(math.acos(cos_angle))


def calculate_angle(structure_path: str, atom_selections,
                    mode: str) -> Optional[float]:
    """
    Calculate angle between selected atoms by parsing PDB file.

    Args:
        structure_path: Path to structure file
        atom_selections: For "bond"/"torsion": list of 3 or 4 selection strings.
                         For "vector": list of two 2-element lists [[a1, a2], [b1, b2]].
        mode: "bond", "torsion", or "vector"

    Returns:
        Calculated angle in degrees, or None if failed
    """
    try:
        # Parse PDB file
        atoms = parse_pdb_file(structure_path)

        if not atoms:
            print(f"  - Warning: No atoms found in {structure_path}")
            return None

        print(f"  - Total atoms in structure: {len(atoms)}")

        def resolve_point(selection):
            selected = resolve_selection(selection, atoms)
            print(f"  - Selection '{selection}': {len(selected)} atoms")
            if not selected:
                raise ValueError(f"No atoms found for selection '{selection}'")
            centroid = get_atom_centroid(selected)
            if len(selected) == 1:
                print(f"    Atom: {selected[0].atom_name} at ({centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f})")
            else:
                print(f"    Centroid of {len(selected)} atoms: ({centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f})")
            return centroid

        if mode == "vector":
            (sel_a1, sel_a2), (sel_b1, sel_b2) = atom_selections
            p_a1 = resolve_point(sel_a1)
            p_a2 = resolve_point(sel_a2)
            p_b1 = resolve_point(sel_b1)
            p_b2 = resolve_point(sel_b2)
            angle = calculate_vector_angle(p_a1, p_a2, p_b1, p_b2)
            print(f"  - Vector-vector angle: {angle:.2f}°")
        elif mode == "torsion":
            points = [resolve_point(s) for s in atom_selections]
            angle = calculate_torsion_angle(points[0], points[1], points[2], points[3])
            print(f"  - Torsional angle: {angle:.2f}°")
        else:  # bond
            points = [resolve_point(s) for s in atom_selections]
            angle = calculate_bond_angle(points[0], points[1], points[2])
            print(f"  - Bond angle: {angle:.2f}°")

        return angle

    except Exception as e:
        print(f"  - Error calculating angle: {e}")
        import traceback
        traceback.print_exc()
        return None


def analyze_angles(config_data: Dict[str, Any]) -> None:
    """
    Analyze angles in structures by parsing PDB files.

    Args:
        config_data: Configuration dictionary with analysis parameters
    """
    # Load structures DataStream
    structures_ds = load_datastream(config_data['structures_json'])

    atom_selections = config_data['atom_selections']
    mode = config_data['mode']
    metric_name = config_data['metric_name']
    unit = config_data.get('unit', 'degrees')
    output_csv = config_data['output_csv']

    angle_type = {"bond": "bond", "torsion": "torsional", "vector": "vector-vector"}[mode]
    print(f"Analyzing {angle_type} angles in {len(structures_ds.ids)} structures")
    if mode == "vector":
        (a1, a2), (b1, b2) = atom_selections
        print(f"Vector 1: {a1} -> {a2}  |  Vector 2: {b1} -> {b2}")
    else:
        print(f"Atom selections: {' -> '.join(atom_selections)}")
    print(f"Output column: {metric_name}")
    print(f"Output unit: {unit}")

    # Process each structure
    results = []
    structure_items = list(iterate_files(structures_ds))
    total = len(structure_items)

    for i, (structure_id, structure_path) in enumerate(structure_items):
        if not os.path.exists(structure_path):
            print(f"Warning: Structure file not found: {structure_path}")
            continue

        print(f"\nProcessing structure {i+1}/{total}: {structure_path}")
        print(f"  - ID: {structure_id}")

        # Calculate angle (always in degrees internally)
        angle = calculate_angle(structure_path, atom_selections, mode)

        # Convert to radians if requested
        if angle is not None and unit == 'radians':
            angle = math.radians(angle)

        # Store result
        result = {
            'id': structure_id,
            'source_structure': structure_path,
            metric_name: angle,
            'unit': unit
        }
        results.append(result)

        print(f"  - Result: {metric_name} = {angle}")

    # Create DataFrame and save
    if results:
        df = pd.DataFrame(results)

        # Create output directory
        output_dir = os.path.dirname(output_csv)
        print(f"\nCreating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)

        # Save results
        print(f"Writing CSV to: {output_csv}")
        df.to_csv(output_csv, index=False)

        # Verify file exists
        if os.path.exists(output_csv):
            file_size = os.path.getsize(output_csv)
            print(f"File verified, size: {file_size} bytes")
        else:
            print("ERROR: File does not exist after write!")
            raise RuntimeError("Failed to write output file")

        print(f"\nAngle analysis completed successfully!")
        print(f"Analyzed {len(results)} structures")
        print(f"Results saved to: {output_csv}")
        print(f"\nResults summary:")
        print(df)

        # Statistics
        angles = [r[metric_name] for r in results if r[metric_name] is not None]
        if angles:
            print(f"\nAngle statistics:")
            print(f"  Min: {min(angles):.2f}°")
            print(f"  Max: {max(angles):.2f}°")
            print(f"  Mean: {np.mean(angles):.2f}°")
            print(f"  Std: {np.std(angles):.2f}°")
        else:
            print(f"\nError: No valid angle measurements found!")
            print(f"All {len(results)} structures returned None for angle calculations.")
            print(f"Check that atom selections match atoms in structures:")
            for sel in atom_selections:
                print(f"  - '{sel}'")
            raise ValueError("No valid angle measurements - all results were None")
    else:
        raise ValueError("No valid results generated - check structure files and selections")


def main():
    parser = argparse.ArgumentParser(description='Analyze angles between atoms in structures')
    parser.add_argument('--config', required=True, help='JSON config file with analysis parameters')

    args = parser.parse_args()

    # Load configuration
    if not os.path.exists(args.config):
        print(f"Error: Config file not found: {args.config}")
        sys.exit(1)

    try:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)

    # Validate required parameters
    required_params = ['structures_json', 'atom_selections', 'mode',
                       'metric_name', 'output_csv']
    for param in required_params:
        if param not in config_data:
            print(f"Error: Missing required parameter: {param}")
            sys.exit(1)

    try:
        analyze_angles(config_data)

    except Exception as e:
        print(f"Error analyzing angles: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
