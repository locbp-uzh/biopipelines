#!/usr/bin/env python3
"""
Global MSA cache management for BioPipelines Boltz2.

This script manages a global cache of MSAs to avoid recalculating them across different jobs.
Uses a database.csv file to map sequences to MSA file names, with MSA files named msa_NNNNNN.csv.
"""

import os
import sys
import pandas as pd
import argparse
import hashlib
import shutil
from datetime import datetime
from typing import Optional, Tuple, List

def get_global_cache_folder(boltz_cache_folder: str) -> str:
    """Get path to global MSA cache folder."""
    return os.path.join(boltz_cache_folder, "GlobalMSAs")

def get_database_path(cache_folder: str) -> str:
    """Get path to MSA database CSV file."""
    return os.path.join(cache_folder, "msa_database.csv")

def ensure_cache_folder(cache_folder: str):
    """Ensure global MSA cache folder and database exist."""
    os.makedirs(cache_folder, exist_ok=True)
    
    database_path = get_database_path(cache_folder)
    if not os.path.exists(database_path):
        # Create empty database with headers
        df = pd.DataFrame(columns=['sequence', 'msa_file', 'date_created', 'sequence_hash'])
        df.to_csv(database_path, index=False)
        print(f"Created MSA database: {database_path}")

def get_sequence_hash(sequence: str) -> str:
    """Generate hash for sequence for faster lookup and verification."""
    return hashlib.md5(sequence.encode()).hexdigest()

def get_next_msa_filename(cache_folder: str) -> str:
    """Get next sequential MSA filename."""
    database_path = get_database_path(cache_folder)
    
    if os.path.exists(database_path):
        try:
            df = pd.read_csv(database_path)
            if not df.empty:
                # Extract numbers from existing msa files
                existing_numbers = []
                for msa_file in df['msa_file']:
                    if msa_file.startswith('msa_') and msa_file.endswith('.csv'):
                        try:
                            num = int(msa_file[4:10])  # Extract NNNNNN part
                            existing_numbers.append(num)
                        except ValueError:
                            continue
                
                if existing_numbers:
                    next_num = max(existing_numbers) + 1
                else:
                    next_num = 1
            else:
                next_num = 1
        except Exception:
            next_num = 1
    else:
        next_num = 1
    
    return f"msa_{next_num:06d}.csv"

def lookup_msa(cache_folder: str, sequence: str) -> Optional[str]:
    """
    Look up MSA file for given sequence.
    
    Args:
        cache_folder: Global MSA cache folder
        sequence: Protein sequence to look up
        
    Returns:
        Path to MSA file if found, None otherwise
    """
    database_path = get_database_path(cache_folder)
    
    if not os.path.exists(database_path):
        return None
    
    try:
        df = pd.read_csv(database_path)
        sequence_hash = get_sequence_hash(sequence)
        
        # Look up by hash first (faster), then by sequence
        matches = df[df['sequence_hash'] == sequence_hash]
        if matches.empty:
            matches = df[df['sequence'] == sequence]
        
        if not matches.empty:
            msa_file = matches.iloc[0]['msa_file']
            msa_path = os.path.join(cache_folder, msa_file)
            
            if os.path.exists(msa_path):
                print(f"Found cached MSA: {msa_file} for sequence hash {sequence_hash[:8]}...")
                return msa_path
            else:
                print(f"Warning: Database entry exists but MSA file not found: {msa_path}")
                # Could remove stale entries here
        
        return None
        
    except Exception as e:
        print(f"Error looking up MSA: {e}")
        return None

def add_msa_to_cache(cache_folder: str, sequence: str, msa_source_path: str) -> str:
    """
    Add new MSA to global cache.
    
    Args:
        cache_folder: Global MSA cache folder
        sequence: Protein sequence
        msa_source_path: Path to source MSA file to copy
        
    Returns:
        Path to cached MSA file
    """
    ensure_cache_folder(cache_folder)
    
    # Get next filename and copy MSA file
    msa_filename = get_next_msa_filename(cache_folder)
    msa_target_path = os.path.join(cache_folder, msa_filename)
    
    try:
        shutil.copy2(msa_source_path, msa_target_path)
        print(f"Copied MSA to cache: {msa_filename}")
        
        # Update database
        database_path = get_database_path(cache_folder)
        df = pd.read_csv(database_path)
        
        new_entry = {
            'sequence': sequence,
            'msa_file': msa_filename,
            'date_created': datetime.now().isoformat(),
            'sequence_hash': get_sequence_hash(sequence)
        }
        
        df = pd.concat([df, pd.DataFrame([new_entry])], ignore_index=True)
        df.to_csv(database_path, index=False)
        
        print(f"Added entry to MSA database: {msa_filename}")
        return msa_target_path
        
    except Exception as e:
        print(f"Error adding MSA to cache: {e}")
        raise

def get_cache_stats(cache_folder: str) -> dict:
    """Get statistics about the global MSA cache."""
    stats = {'total_msas': 0, 'database_entries': 0, 'orphaned_files': 0}
    
    if not os.path.exists(cache_folder):
        return stats
    
    # Count database entries
    database_path = get_database_path(cache_folder)
    if os.path.exists(database_path):
        try:
            df = pd.read_csv(database_path)
            stats['database_entries'] = len(df)
        except Exception:
            pass
    
    # Count MSA files
    try:
        msa_files = [f for f in os.listdir(cache_folder) if f.startswith('msa_') and f.endswith('.csv')]
        stats['total_msas'] = len(msa_files)
    except Exception:
        pass
    
    stats['orphaned_files'] = max(0, stats['total_msas'] - stats['database_entries'])
    
    return stats

def main():
    parser = argparse.ArgumentParser(description='Manage global MSA cache')
    parser.add_argument('action', choices=['lookup', 'add', 'stats'], help='Action to perform')
    parser.add_argument('--cache-folder', required=True, help='Global MSA cache folder')
    parser.add_argument('--sequence', help='Protein sequence (for lookup/add)')
    parser.add_argument('--msa-file', help='MSA file to add to cache (for add action)')
    parser.add_argument('--output', help='Output file path for lookup result')
    
    args = parser.parse_args()
    
    cache_folder = get_global_cache_folder(args.cache_folder)
    
    if args.action == 'lookup':
        if not args.sequence:
            print("Error: --sequence required for lookup action")
            sys.exit(1)
        
        msa_path = lookup_msa(cache_folder, args.sequence)
        if msa_path:
            print(f"Found MSA: {msa_path}")
            if args.output:
                shutil.copy2(msa_path, args.output)
                print(f"Copied to: {args.output}")
        else:
            print("MSA not found in cache")
            sys.exit(1)
    
    elif args.action == 'add':
        if not args.sequence or not args.msa_file:
            print("Error: --sequence and --msa-file required for add action")
            sys.exit(1)
        
        if not os.path.exists(args.msa_file):
            print(f"Error: MSA file not found: {args.msa_file}")
            sys.exit(1)
        
        cached_path = add_msa_to_cache(cache_folder, args.sequence, args.msa_file)
        print(f"MSA cached at: {cached_path}")
    
    elif args.action == 'stats':
        stats = get_cache_stats(cache_folder)
        print("Global MSA Cache Statistics:")
        print(f"  Database entries: {stats['database_entries']}")
        print(f"  MSA files: {stats['total_msas']}")
        print(f"  Orphaned files: {stats['orphaned_files']}")
        print(f"  Cache folder: {cache_folder}")

if __name__ == "__main__":
    main()