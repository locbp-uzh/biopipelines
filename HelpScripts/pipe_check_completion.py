#!/usr/bin/env python3
# Authors (2026): G. Quargnali & P. Rivera-Fuentes @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# Licensed under the MIT License. See LICENSE file in the project root for details.

"""
Generic completion checker for pipeline tools.

Checks if expected output files exist and creates COMPLETED/FAILED status files accordingly.
This allows bash scripts to skip already completed steps and provides clear status indication.
"""

import os
import sys
import argparse
import json
from pathlib import Path
from typing import Dict, List, Any, Optional

def check_file_exists(file_path: str) -> bool:
    """
    Check if a file or directory exists.
    
    Args:
        file_path: Path to check
        
    Returns:
        True if exists, False otherwise
    """
    return os.path.exists(file_path) and (os.path.isfile(file_path) or os.path.isdir(file_path))

def load_expected_missing_ids(missing_csv_path: str) -> List[str]:
    """
    Load expected missing IDs from missing CSV file if it exists.

    Args:
        missing_csv_path: Path to the missing CSV file

    Returns:
        List of IDs that are expected to be missing
    """
    missing_ids = []

    if os.path.exists(missing_csv_path):
        try:
            import pandas as pd
            df = pd.read_csv(missing_csv_path)

            if 'id' in df.columns:
                missing_ids = df['id'].astype(str).tolist()

            print(f"Found {len(missing_ids)} IDs expected to be missing from {os.path.basename(missing_csv_path)}")
        except Exception as e:
            print(f"Warning: Could not read {os.path.basename(missing_csv_path)}: {e}")

    return missing_ids

def check_files_exist(file_list: List[str]) -> tuple[bool, List[str]]:
    """
    Check if all files in a list exist.
    
    Args:
        file_list: List of file paths to check
        
    Returns:
        Tuple of (all_exist: bool, missing_files: List[str])
    """
    missing_files = []
    for file_path in file_list:
        if not check_file_exists(file_path):
            missing_files.append(file_path)
    
    return len(missing_files) == 0, missing_files

def extract_file_list(category_data) -> List[str]:
    """
    Extract file paths from a category entry in expected outputs.

    The category data may be:
    - A list of file paths (legacy format)
    - A dict (serialized DataStream) with a 'files' key containing the list

    Args:
        category_data: List or dict from expected_outputs[category]

    Returns:
        List of file paths
    """
    if isinstance(category_data, dict):
        return category_data.get('files', [])
    elif isinstance(category_data, list):
        return category_data
    return []


def is_filter_output(expected_outputs: Dict[str, Any]) -> bool:
    """
    Check if the expected outputs are from a Filter tool (not RemoveDuplicates).
    
    Args:
        expected_outputs: Dictionary with standardized output format
        
    Returns:
        True if this is Filter tool output
    """
    # Check if tables contain missing AND there are structure files expected
    # This distinguishes Filter (works with structures) from RemoveDuplicates (sequences only)
    if 'tables' in expected_outputs and isinstance(expected_outputs['tables'], dict):
        has_missing_table = 'missing' in expected_outputs['tables']
        has_structures = bool(extract_file_list(expected_outputs.get('structures', [])))

        # Only treat as filter output if it has missing table AND expects structures
        return has_missing_table and has_structures
    
    return False




def check_expected_outputs_filter_aware(expected_outputs: Dict[str, Any], 
                                       output_folder: str,
                                       job_name: str) -> tuple[bool, Dict[str, List[str]], Dict[str, Any]]:
    """
    Check expected outputs with filter-aware validation.
    
    For filter tools:
    - Tables and manifests are required (critical)
    - Content files (structures, sequences) can be partially missing (warning)
    
    Args:
        expected_outputs: Dictionary with standardized output format
        output_folder: Tool's output folder
        job_name: Job name
        
    Returns:
        Tuple of (success: bool, missing_by_category: Dict, filter_info: Dict)
    """
    missing_by_category = {}
    warnings_by_category = {}
    filter_info = {}
    
    # Check if this is filter output
    is_filter = is_filter_output(expected_outputs)
    filter_info['is_filter'] = is_filter
    
    if is_filter:
        print("Filter-aware validation: This appears to be filter output")
        
        # For filters, categorize checks differently
        critical_files = []
        content_files = []
        
        # Tables are critical
        if 'tables' in expected_outputs:
            tables = expected_outputs['tables']
            if isinstance(tables, dict):
                for name, info in tables.items():
                    if isinstance(info, dict) and 'path' in info:
                        critical_files.append(info['path'])
                    else:
                        critical_files.append(str(info))
            elif isinstance(tables, list):
                critical_files.extend(tables)
        
        # Content files are non-critical for filters
        standard_categories = ['structures', 'compounds', 'sequences']
        for category in standard_categories:
            if category in expected_outputs and expected_outputs[category]:
                content_files.extend(extract_file_list(expected_outputs[category]))
        
        # Check critical files (must exist)
        if critical_files:
            exists, missing = check_files_exist(critical_files)
            if not exists:
                missing_by_category['critical'] = missing
        
        # Check content files (can be partially missing)
        if content_files:
            # Find missing CSV path from expected outputs
            missing_csv_path = None
            if 'tables' in expected_outputs and isinstance(expected_outputs['tables'], dict):
                if 'missing' in expected_outputs['tables']:
                    missing_info = expected_outputs['tables']['missing']
                    if isinstance(missing_info, dict) and 'path' in missing_info:
                        missing_csv_path = missing_info['path']
                    else:
                        missing_csv_path = str(missing_info)
            
            # Load expected missing IDs from missing CSV (if path found)
            missing_ids = []
            if missing_csv_path:
                missing_ids = load_expected_missing_ids(missing_csv_path)
            else:
                print("Warning: No missing CSV path found in expected outputs - cannot determine expected missing files")

            # Build set of expected-missing files by matching IDs to content file basenames
            expected_missing_files = set()
            for f in content_files:
                basename = os.path.splitext(os.path.basename(f))[0]
                for mid in missing_ids:
                    if mid == basename or basename.startswith(f"{mid}_") or basename.startswith(f"{mid}-"):
                        expected_missing_files.add(f)
                        break

            exists, missing = check_files_exist(content_files)
            if not exists:
                # Filter out expected missing files
                unexpected_missing = [f for f in missing if f not in expected_missing_files]
                expected_missing_found = [f for f in missing if f in expected_missing_files]
                
                if unexpected_missing:
                    warnings_by_category['content'] = unexpected_missing
                    print(f"Warning: {len(unexpected_missing)} content files unexpectedly missing")
                
                if expected_missing_found:
                    print(f"Info: {len(expected_missing_found)} content files missing as expected (filtered out)")
        
        # Success if all critical files exist
        success = 'critical' not in missing_by_category
        
        # Add warnings to result for reporting
        if warnings_by_category:
            filter_info['warnings'] = warnings_by_category
    
    else:
        # Standard validation for non-filter tools
        success, missing_by_category = check_expected_outputs(expected_outputs)
        filter_info['is_filter'] = False
    
    return success, missing_by_category, filter_info


def check_expected_outputs(expected_outputs: Dict[str, Any]) -> tuple[bool, Dict[str, List[str]]]:
    """
    Check if expected output files exist based on standardized output format.
    
    Args:
        expected_outputs: Dictionary with standardized output format
        
    Returns:
        Tuple of (all_exist: bool, missing_by_category: Dict[str, List[str]])
    """
    missing_by_category = {}
    all_exist = True
    
    # Check standard categories
    standard_categories = ['structures', 'compounds', 'sequences']
    
    for category in standard_categories:
        if category in expected_outputs:
            files = extract_file_list(expected_outputs[category])
            if files:
                exists, missing = check_files_exist(files)
                if not exists:
                    missing_by_category[category] = missing
                    all_exist = False
    
    # Check tables (handle both old and new format)
    if 'tables' in expected_outputs:
        tables = expected_outputs['tables']
        table_files = []
        
        if isinstance(tables, dict):
            # New format with named tables
            for name, info in tables.items():
                if isinstance(info, dict) and 'path' in info:
                    table_files.append(info['path'])
                else:
                    # Fallback - treat as path
                    table_files.append(str(info))
        elif isinstance(tables, list):
            # Old format - list of paths
            table_files = tables
        
        if table_files:
            exists, missing = check_files_exist(table_files)
            if not exists:
                missing_by_category['tables'] = missing
                all_exist = False
    
    return all_exist, missing_by_category

def create_status_file(output_folder: str, tool_name: str, status: str, details: Dict[str, Any] = None) -> str:
    """
    Create a status file indicating completion or failure.
    
    Args:
        output_folder: Tool's output folder
        tool_name: Name of the tool (e.g., "RFdiffusion", "ProteinMPNN")
        status: "COMPLETED" or "FAILED"
        details: Optional details about the status
        
    Returns:
        Path to created status file
    """
    # Get step number from output folder name (e.g., "1_RFdiffusion" -> "1")
    folder_name = os.path.basename(output_folder)
    if '_' in folder_name and folder_name.split('_')[0].isdigit():
        step_number = folder_name.split('_')[0]
        status_filename = f"{step_number}_{tool_name}_{status}"
    else:
        status_filename = f"{tool_name}_{status}"
    
    # Create status file in the parent directory (pipeline level)
    parent_dir = os.path.dirname(output_folder)
    status_file_path = os.path.join(parent_dir, status_filename)
    
    # Create the status file with optional details
    with open(status_file_path, 'w') as f:
        f.write(f"Status: {status}\n")
        f.write(f"Tool: {tool_name}\n")
        f.write(f"Output folder: {output_folder}\n")
        f.write(f"Timestamp: {__import__('datetime').datetime.now().isoformat()}\n")
        
        if details:
            f.write(f"\nDetails:\n")
            if 'missing_files' in details:
                f.write(f"Missing files:\n")
                for category, files in details['missing_files'].items():
                    f.write(f"  {category}:\n")
                    for file_path in files:
                        f.write(f"    - {file_path}\n")
            
            if 'error_message' in details:
                f.write(f"Error: {details['error_message']}\n")
    
    return status_file_path

def check_completion_status(output_folder: str, tool_name: str) -> Optional[str]:
    """
    Check if completion status file already exists.
    
    Args:
        output_folder: Tool's output folder
        tool_name: Name of the tool
        
    Returns:
        Status ("COMPLETED" or "FAILED") if exists, None otherwise
    """
    parent_dir = os.path.dirname(output_folder)
    folder_name = os.path.basename(output_folder)
    
    # Try step-numbered format first
    if '_' in folder_name and folder_name.split('_')[0].isdigit():
        step_number = folder_name.split('_')[0]
        completed_file = os.path.join(parent_dir, f"{step_number}_{tool_name}_COMPLETED")
        failed_file = os.path.join(parent_dir, f"{step_number}_{tool_name}_FAILED")
    else:
        completed_file = os.path.join(parent_dir, f"{tool_name}_COMPLETED")
        failed_file = os.path.join(parent_dir, f"{tool_name}_FAILED")
    
    if os.path.exists(completed_file):
        return "COMPLETED"
    elif os.path.exists(failed_file):
        return "FAILED"
    
    return None

def clean_old_status_files(output_folder: str, tool_name: str) -> None:
    """
    Remove old status files (both COMPLETED and FAILED) to allow fresh evaluation.
    
    Args:
        output_folder: Tool's output folder
        tool_name: Name of the tool
    """
    parent_dir = os.path.dirname(output_folder)
    folder_name = os.path.basename(output_folder)
    
    # Determine status file patterns
    if '_' in folder_name and folder_name.split('_')[0].isdigit():
        step_number = folder_name.split('_')[0]
        completed_file = os.path.join(parent_dir, f"{step_number}_{tool_name}_COMPLETED")
        failed_file = os.path.join(parent_dir, f"{step_number}_{tool_name}_FAILED")
    else:
        completed_file = os.path.join(parent_dir, f"{tool_name}_COMPLETED")
        failed_file = os.path.join(parent_dir, f"{tool_name}_FAILED")
    
    # Remove old status files if they exist
    for status_file in [completed_file, failed_file]:
        if os.path.exists(status_file):
            try:
                os.remove(status_file)
                print(f"Removed old status file: {os.path.basename(status_file)}")
            except OSError as e:
                print(f"Warning: Could not remove {os.path.basename(status_file)}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Check pipeline tool completion status")
    parser.add_argument("output_folder", help="Tool's output folder path")
    parser.add_argument("tool_name", help="Name of the tool (e.g., RFdiffusion, ProteinMPNN)")
    parser.add_argument("expected_outputs", help="JSON file or string with expected outputs")
    parser.add_argument("--check-only", action="store_true", 
                       help="Only check status, don't create status files")
    parser.add_argument("--force", action="store_true",
                       help="Force check even if status file exists")
    parser.add_argument("--job-name", 
                       help="Job name for filter manifest lookup (auto-detected if not provided)")
    
    args = parser.parse_args()
    
    # Auto-detect job name from output folder if not provided
    job_name = args.job_name
    if not job_name:
        # Try to extract from folder structure or use tool name as fallback
        folder_name = os.path.basename(args.output_folder)
        if '_' in folder_name and len(folder_name.split('_')) > 1:
            # Format like "1_ToolName" - use ToolName as job name
            job_name = '_'.join(folder_name.split('_')[1:])
        else:
            job_name = args.tool_name
    
    # Clean up old status files to allow fresh evaluation
    clean_old_status_files(args.output_folder, args.tool_name)
    
    # Check if we already have a status (unless forced)
    if not args.force:
        existing_status = check_completion_status(args.output_folder, args.tool_name)
        if existing_status == "COMPLETED":
            print(f"Tool {args.tool_name} already completed")
            sys.exit(0)
        elif existing_status == "FAILED":
            print(f"Tool {args.tool_name} previously failed")
            sys.exit(1)
    
    # Load expected outputs
    try:
        if os.path.isfile(args.expected_outputs):
            # Load from JSON file
            with open(args.expected_outputs, 'r') as f:
                expected_outputs = json.load(f)
        else:
            # Parse as JSON string
            expected_outputs = json.loads(args.expected_outputs)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"Error parsing expected outputs: {e}", file=sys.stderr)
        sys.exit(2)

    # Unwrap envelope format if present (tool_name, tool_class, output_structure wrapper)
    if 'output_structure' in expected_outputs and 'tool_name' in expected_outputs:
        expected_outputs = expected_outputs['output_structure']
    
    # Use filter-aware validation
    success, missing_by_category, filter_info = check_expected_outputs_filter_aware(
        expected_outputs, args.output_folder, job_name
    )
    
    if success:
        print(f"Required outputs found for {args.tool_name}")
        
        # Report any warnings for filter tools
        if filter_info.get('warnings'):
            print("Warnings (non-critical for filters):")
            for category, files in filter_info['warnings'].items():
                print(f"  {category}: {len(files)} files missing")
                # Optionally list first few missing files
                for file_path in files[:3]:
                    print(f"    - {os.path.basename(file_path)}")
                if len(files) > 3:
                    print(f"    ... and {len(files) - 3} more")
        
        if not args.check_only:
            # Include filter info in completion details
            details = {}
            if filter_info.get('is_filter'):
                details['filter_info'] = filter_info
            
            status_file = create_status_file(args.output_folder, args.tool_name, "COMPLETED", details)
            print(f"Created completion status file: {os.path.basename(status_file)}")
        sys.exit(0)
    else:
        print(f"Missing critical outputs for {args.tool_name}:")
        for category, files in missing_by_category.items():
            print(f"  {category}:")
            for file_path in files:
                print(f"    - {file_path}")
        
        if not args.check_only:
            details = {"missing_files": missing_by_category}
            if filter_info.get('is_filter'):
                details['filter_info'] = filter_info
            
            status_file = create_status_file(args.output_folder, args.tool_name, "FAILED", details)
            print(f"Created failure status file: {os.path.basename(status_file)}")
        sys.exit(1)

if __name__ == "__main__":
    main()