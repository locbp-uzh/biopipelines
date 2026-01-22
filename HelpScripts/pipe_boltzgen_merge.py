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


def update_csv_ids_regex(content: str, prefix: str) -> str:
    """
    Update all design IDs in CSV content using a single regex replacement.

    Since CSV files contain only filenames (not paths), we can use a regex
    pattern to match all design IDs and prepend the prefix in one pass.

    Args:
        content: CSV file content as string
        prefix: Prefix to prepend to all design IDs

    Returns:
        Updated content with renamed IDs
    """
    # Pattern matches design IDs like: design_spec_001, design_spec_001_model_0, etc.
    # The pattern captures the full ID to preserve suffixes like _model_0
    pattern = r'(design_spec_\d+(?:_model_\d+)?)'
    return re.sub(pattern, prefix + r'\1', content)


def collect_csv_headers(csv_files: list) -> list:
    """
    Collect all unique headers from multiple CSV files, preserving order.

    Args:
        csv_files: List of (source_path, prefix) tuples

    Returns:
        List of all unique column names in order of first appearance
    """
    all_columns = []
    seen = set()

    for csv_path, _ in csv_files:
        with open(csv_path, 'r') as f:
            header_line = f.readline().strip()
            if header_line:
                columns = header_line.split(',')
                for col in columns:
                    if col not in seen:
                        seen.add(col)
                        all_columns.append(col)

    return all_columns


def merge_csv_files_with_headers(csv_files: list, output_path: str) -> None:
    """
    Merge multiple CSV files into one, handling different headers.

    Collects all unique headers across files, then merges rows with
    proper column alignment. Missing columns get empty values.

    Args:
        csv_files: List of (source_csv_path, prefix) tuples
        output_path: Output merged CSV path
    """
    if not csv_files:
        return

    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # First pass: collect all unique headers
    all_columns = collect_csv_headers(csv_files)

    with open(output_path, 'w') as out_f:
        # Write unified header
        out_f.write(','.join(all_columns) + '\n')

        # Second pass: merge data with proper column alignment
        for csv_path, prefix in csv_files:
            with open(csv_path, 'r') as in_f:
                content = in_f.read()

            # Apply ID renaming with single regex
            content = update_csv_ids_regex(content, prefix)

            lines = content.strip().split('\n')
            if not lines:
                continue

            # Parse header to get column indices
            file_columns = lines[0].split(',')
            col_to_idx = {col: idx for idx, col in enumerate(file_columns)}

            # Process data rows
            for line in lines[1:]:
                if not line.strip():
                    continue

                values = line.split(',')
                # Build row with proper column alignment
                row = []
                for col in all_columns:
                    if col in col_to_idx and col_to_idx[col] < len(values):
                        row.append(values[col_to_idx[col]])
                    else:
                        row.append('')  # Empty for missing columns

                out_f.write(','.join(row) + '\n')


def rename_file_with_prefix(source_path: str, dest_dir: str, prefix: str) -> str:
    """
    Copy a file while renaming design IDs in the filename using regex.

    Args:
        source_path: Path to source file
        dest_dir: Destination directory
        prefix: Prefix to prepend to design IDs

    Returns:
        Path to the new file
    """
    filename = os.path.basename(source_path)
    # Apply same regex pattern to filename
    new_filename = re.sub(r'(design_spec_\d+(?:_model_\d+)?)', prefix + r'\1', filename)
    dest_path = os.path.join(dest_dir, new_filename)

    os.makedirs(dest_dir, exist_ok=True)
    shutil.copy2(source_path, dest_path)

    return dest_path


def merge_boltzgen_outputs(sources: list, output_dir: str, id_template: str,
                           verbose: bool = True) -> dict:
    """
    Merge multiple BoltzGen outputs with ID renaming.

    Optimized implementation:
    - Uses single regex for ID renaming instead of per-ID string replacement
    - Merges CSV files directly without temp files
    - Handles different CSV headers by collecting all columns first

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

    # Track CSV files to merge: csv_name -> list of (source_path, prefix)
    csv_files_to_merge = {}

    for i, source in enumerate(sources):
        prefix = get_source_prefix(source, i, id_template)
        if verbose:
            print(f"\nProcessing source {i+1}/{len(sources)}: {os.path.basename(source)}")
            print(f"  Prefix: {prefix}")

        files = find_design_files(source)
        unique_ids = set()

        # Process structure files (cif, npz) with regex-based renaming
        for file_type in ['cif', 'npz']:
            for rel_path in files[file_type]:
                source_path = os.path.join(source, rel_path)
                filename = os.path.basename(rel_path)

                # Extract ID for counting (optional, for stats)
                match = re.search(r'design_spec_\d+(?:_model_\d+)?', filename)
                if match:
                    unique_ids.add(match.group(0))

                # Determine destination directory
                dest_subdir = os.path.dirname(rel_path)
                dest_dir = os.path.join(output_dir, dest_subdir)

                rename_file_with_prefix(source_path, dest_dir, prefix)
                stats['files_copied'] += 1

        stats['ids_renamed'] += len(unique_ids)

        # Collect CSV files with their prefixes (no temp files needed)
        for rel_path in files['csv']:
            source_path = os.path.join(source, rel_path)
            csv_name = os.path.basename(rel_path)

            if csv_name not in csv_files_to_merge:
                csv_files_to_merge[csv_name] = {'rel_path': rel_path, 'files': []}
            csv_files_to_merge[csv_name]['files'].append((source_path, prefix))

        if verbose:
            print(f"  Copied {len(files['cif'])} CIF files, {len(files['npz'])} NPZ files")
            print(f"  Found {len(unique_ids)} unique design IDs")

    # Merge CSV files with header alignment
    if verbose:
        print(f"\nMerging CSV files...")

    for csv_name, csv_info in csv_files_to_merge.items():
        file_list = csv_info['files']
        if file_list:
            output_path = os.path.join(output_dir, csv_info['rel_path'])
            merge_csv_files_with_headers(file_list, output_path)
            stats['csvs_merged'][csv_name] = len(file_list)

            if verbose:
                print(f"  {csv_name}: merged {len(file_list)} files")

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
