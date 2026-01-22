#!/usr/bin/env python3
"""
Custom BoltzGen merge script with ID renaming.

This script merges multiple BoltzGen output directories while renaming
design IDs to avoid collisions. It maintains the folder structure and
updates all file references in CSV files.

Usage:
    python pipe_boltzgen_merge.py \
        --sources /path/to/run1 /path/to/run2 ... \
        --output /path/to/merged \
        --id_template "batch{i:03d}_"
"""

import os
import sys
import argparse
import shutil
import re
from pathlib import Path


def get_source_prefix(source_path: str, index: int, id_template: str) -> str:
    """
    Generate a unique prefix for a source directory.

    Args:
        source_path: Path to source BoltzGen output
        index: Index of this source (0-based)
        id_template: Template string for prefix generation.
                    Supports {i} for index, {name} for folder name.

    Returns:
        Prefix string to prepend to design IDs
    """
    source_name = os.path.basename(source_path.rstrip('/\\'))
    return id_template.format(i=index, name=source_name)


def find_design_files(source_dir: str) -> dict:
    """
    Find all design-related files in a BoltzGen output directory.

    Returns:
        Dictionary mapping relative paths to file info
    """
    files = {
        'cif': [],      # Structure files
        'npz': [],      # Metadata files
        'csv': [],      # Metrics files
        'other': []     # Other files to copy as-is
    }

    # Directories to scan for design files
    design_dirs = [
        'intermediate_designs',
        'intermediate_designs_inverse_folded',
        'intermediate_designs_inverse_folded/refold_cif',
        'intermediate_designs_inverse_folded/refold_design_cif',
    ]

    for design_dir in design_dirs:
        full_dir = os.path.join(source_dir, design_dir)
        if not os.path.exists(full_dir):
            continue

        for entry in os.listdir(full_dir):
            full_path = os.path.join(full_dir, entry)
            if not os.path.isfile(full_path):
                continue

            rel_path = os.path.join(design_dir, entry)

            if entry.endswith('.cif'):
                files['cif'].append(rel_path)
            elif entry.endswith('.npz'):
                files['npz'].append(rel_path)
            elif entry.endswith('.csv'):
                files['csv'].append(rel_path)

    # CSV files in intermediate_designs_inverse_folded root
    inv_folded_dir = os.path.join(source_dir, 'intermediate_designs_inverse_folded')
    if os.path.exists(inv_folded_dir):
        for entry in os.listdir(inv_folded_dir):
            if entry.endswith('.csv'):
                rel_path = os.path.join('intermediate_designs_inverse_folded', entry)
                if rel_path not in files['csv']:
                    files['csv'].append(rel_path)

    return files


def extract_design_id(filename: str) -> str:
    """
    Extract design ID from a filename.

    Examples:
        design_spec_001.cif -> design_spec_001
        design_spec_001_model_0.cif -> design_spec_001_model_0
    """
    # Remove extension
    name = os.path.splitext(filename)[0]
    return name


def rename_design_id(old_id: str, prefix: str) -> str:
    """
    Rename a design ID by prepending a prefix.

    Args:
        old_id: Original design ID (e.g., "design_spec_001")
        prefix: Prefix to prepend (e.g., "batch001_")

    Returns:
        New design ID (e.g., "batch001_design_spec_001")
    """
    return f"{prefix}{old_id}"


def copy_and_rename_file(source_path: str, dest_dir: str, old_id: str, new_id: str) -> str:
    """
    Copy a file while renaming based on ID change.

    Returns:
        Path to the new file
    """
    filename = os.path.basename(source_path)
    new_filename = filename.replace(old_id, new_id)
    dest_path = os.path.join(dest_dir, new_filename)

    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
    shutil.copy2(source_path, dest_path)

    return dest_path


def update_csv_ids(csv_path: str, id_mappings: dict) -> None:
    """
    Update all ID references in a CSV file.

    Args:
        csv_path: Path to CSV file
        id_mappings: Dictionary mapping old_id -> new_id
    """
    with open(csv_path, 'r') as f:
        content = f.read()

    # Sort by length (longest first) to avoid partial replacements
    sorted_mappings = sorted(id_mappings.items(), key=lambda x: len(x[0]), reverse=True)

    for old_id, new_id in sorted_mappings:
        content = content.replace(old_id, new_id)

    with open(csv_path, 'w') as f:
        f.write(content)


def merge_csv_files(csv_files: list, output_path: str) -> None:
    """
    Merge multiple CSV files into one, keeping header from first file.

    Args:
        csv_files: List of CSV file paths
        output_path: Output merged CSV path
    """
    if not csv_files:
        return

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, 'w') as out_f:
        for i, csv_file in enumerate(csv_files):
            with open(csv_file, 'r') as in_f:
                lines = in_f.readlines()
                if i == 0:
                    # Write header and all data from first file
                    out_f.writelines(lines)
                else:
                    # Skip header for subsequent files
                    if len(lines) > 1:
                        out_f.writelines(lines[1:])


def merge_boltzgen_outputs(sources: list, output_dir: str, id_template: str,
                           verbose: bool = True) -> dict:
    """
    Merge multiple BoltzGen outputs with ID renaming.

    Args:
        sources: List of source directory paths
        output_dir: Output directory for merged results
        id_template: Template for generating ID prefixes
        verbose: Print progress information

    Returns:
        Dictionary with merge statistics
    """
    stats = {
        'sources': len(sources),
        'files_copied': 0,
        'ids_renamed': 0,
        'csvs_merged': {}
    }

    # Create output directory structure
    os.makedirs(output_dir, exist_ok=True)
    for subdir in ['intermediate_designs',
                   'intermediate_designs_inverse_folded',
                   'intermediate_designs_inverse_folded/refold_cif',
                   'intermediate_designs_inverse_folded/refold_design_cif',
                   'config']:
        os.makedirs(os.path.join(output_dir, subdir), exist_ok=True)

    # Track all ID mappings and CSV files to merge
    all_id_mappings = {}
    csv_files_to_merge = {}  # csv_name -> list of paths

    for i, source in enumerate(sources):
        prefix = get_source_prefix(source, i, id_template)
        if verbose:
            print(f"\nProcessing source {i+1}/{len(sources)}: {os.path.basename(source)}")
            print(f"  Prefix: {prefix}")

        files = find_design_files(source)
        source_id_mappings = {}

        # Process structure files (cif, npz)
        for file_type in ['cif', 'npz']:
            for rel_path in files[file_type]:
                source_path = os.path.join(source, rel_path)
                filename = os.path.basename(rel_path)
                old_id = extract_design_id(filename)
                new_id = rename_design_id(old_id, prefix)

                source_id_mappings[old_id] = new_id

                # Determine destination directory
                dest_subdir = os.path.dirname(rel_path)
                dest_dir = os.path.join(output_dir, dest_subdir)

                copy_and_rename_file(source_path, dest_dir, old_id, new_id)
                stats['files_copied'] += 1

        stats['ids_renamed'] += len(source_id_mappings)
        all_id_mappings.update(source_id_mappings)

        # Process CSV files - copy with renamed IDs, then merge later
        for rel_path in files['csv']:
            source_path = os.path.join(source, rel_path)
            csv_name = os.path.basename(rel_path)

            # Create temp copy with renamed IDs
            temp_dir = os.path.join(output_dir, '_temp_csv', f'source_{i}')
            os.makedirs(temp_dir, exist_ok=True)
            temp_path = os.path.join(temp_dir, csv_name)

            shutil.copy2(source_path, temp_path)
            update_csv_ids(temp_path, source_id_mappings)

            # Track for merging
            if csv_name not in csv_files_to_merge:
                csv_files_to_merge[csv_name] = []
            csv_files_to_merge[csv_name].append((rel_path, temp_path))

        if verbose:
            print(f"  Copied {len(files['cif'])} CIF files, {len(files['npz'])} NPZ files")
            print(f"  Renamed {len(source_id_mappings)} unique IDs")

    # Merge CSV files
    if verbose:
        print(f"\nMerging CSV files...")

    for csv_name, file_list in csv_files_to_merge.items():
        if file_list:
            # Use the relative path from first file
            rel_path = file_list[0][0]
            output_path = os.path.join(output_dir, rel_path)
            temp_paths = [fp[1] for fp in file_list]

            merge_csv_files(temp_paths, output_path)
            stats['csvs_merged'][csv_name] = len(file_list)

            if verbose:
                print(f"  {csv_name}: merged {len(file_list)} files")

    # Clean up temp directory
    temp_dir = os.path.join(output_dir, '_temp_csv')
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)

    # Copy config and other files from first source
    first_source = sources[0]
    config_dir = os.path.join(first_source, 'config')
    if os.path.exists(config_dir):
        for entry in os.listdir(config_dir):
            src = os.path.join(config_dir, entry)
            dst = os.path.join(output_dir, 'config', entry)
            if os.path.isfile(src):
                shutil.copy2(src, dst)

    # Copy design_spec.yaml if exists
    design_spec = os.path.join(first_source, 'design_spec.yaml')
    if os.path.exists(design_spec):
        shutil.copy2(design_spec, os.path.join(output_dir, 'design_spec.yaml'))

    # Copy steps.yaml if exists
    steps_yaml = os.path.join(first_source, 'steps.yaml')
    if os.path.exists(steps_yaml):
        shutil.copy2(steps_yaml, os.path.join(output_dir, 'steps.yaml'))

    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Merge BoltzGen outputs with ID renaming to avoid collisions'
    )
    parser.add_argument(
        '--sources',
        nargs='+',
        required=True,
        help='Source BoltzGen output directories to merge'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output directory for merged results'
    )
    parser.add_argument(
        '--id_template',
        default='batch{i:03d}_',
        help='Template for ID prefixes. Supports {i} for index (0-based), '
             '{name} for source folder name. Default: "batch{i:03d}_"'
    )
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Suppress progress output'
    )

    args = parser.parse_args()

    # Validate sources exist
    for source in args.sources:
        if not os.path.exists(source):
            print(f"Error: Source directory not found: {source}")
            sys.exit(1)

    print("=" * 60)
    print("BoltzGen Merge with ID Renaming")
    print("=" * 60)
    print(f"Sources: {len(args.sources)} directories")
    print(f"Output: {args.output}")
    print(f"ID template: {args.id_template}")

    stats = merge_boltzgen_outputs(
        sources=args.sources,
        output_dir=args.output,
        id_template=args.id_template,
        verbose=not args.quiet
    )

    print("\n" + "=" * 60)
    print("Merge Summary")
    print("=" * 60)
    print(f"Sources processed: {stats['sources']}")
    print(f"Files copied: {stats['files_copied']}")
    print(f"IDs renamed: {stats['ids_renamed']}")
    print(f"CSV files merged: {len(stats['csvs_merged'])}")
    for csv_name, count in stats['csvs_merged'].items():
        print(f"  - {csv_name}: {count} sources")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
