#!/usr/bin/env python3
"""
Script to extract sequences from MSA files and update global cache.

This script processes MSA files in a folder, extracts the query sequence (first sequence)
from each MSA file, and updates the global MSA cache if the sequence is not already cached.
"""

import os
import sys
import pandas as pd
import argparse
from pipe_manage_global_msa_cache import (
    get_global_cache_folder, lookup_msa, add_msa_to_cache, ensure_cache_folder
)

def extract_query_sequence_from_msa(msa_file_path: str) -> str:
    """
    Extract the query sequence (first sequence) from an MSA file.
    
    Args:
        msa_file_path: Path to MSA file (CSV format with key,sequence structure)
        
    Returns:
        Query sequence string, or empty string if extraction fails
    """
    try:
        df = pd.read_csv(msa_file_path)
        if df.empty:
            print(f"Warning: Empty MSA file: {msa_file_path}")
            return ""
        
        # First row should contain the query sequence
        if 'sequence' not in df.columns:
            print(f"Warning: No 'sequence' column in MSA file: {msa_file_path}")
            return ""
        
        query_sequence = df.iloc[0]['sequence'].strip()
        # Remove any gap characters (dashes) that might be in the query sequence
        query_sequence = query_sequence.replace('-', '')
        
        return query_sequence
        
    except Exception as e:
        print(f"Error extracting sequence from {msa_file_path}: {e}")
        return ""

def process_msa_folder(msas_folder: str, boltz_cache_folder: str) -> dict:
    """
    Process all MSA files in a folder and update global cache.
    
    Args:
        msas_folder: Path to folder containing MSA files
        boltz_cache_folder: Path to Boltz cache folder (for global cache)
        
    Returns:
        Dictionary with processing statistics
    """
    stats = {
        'total_files': 0,
        'sequences_extracted': 0,
        'already_cached': 0,
        'newly_cached': 0,
        'errors': 0
    }
    
    if not os.path.exists(msas_folder):
        print(f"MSA folder does not exist: {msas_folder}")
        return stats
    
    # Get global cache folder
    cache_folder = get_global_cache_folder(boltz_cache_folder)
    ensure_cache_folder(cache_folder)
    
    # Process each MSA file
    msa_files = [f for f in os.listdir(msas_folder) if f.endswith('.csv') or f.endswith('.a3m')]
    stats['total_files'] = len(msa_files)
    
    print(f"Processing {len(msa_files)} MSA files in: {msas_folder}")
    
    for msa_file in msa_files:
        msa_file_path = os.path.join(msas_folder, msa_file)
        sequence_id = os.path.splitext(msa_file)[0]  # Use filename without extension as sequence ID
        
        try:
            # Extract query sequence from MSA
            query_sequence = extract_query_sequence_from_msa(msa_file_path)
            if not query_sequence:
                print(f"Could not extract sequence from: {msa_file}")
                stats['errors'] += 1
                continue
            
            stats['sequences_extracted'] += 1
            
            # Check if sequence is already in global cache
            cached_msa_path = lookup_msa(cache_folder, query_sequence)
            if cached_msa_path:
                print(f"Sequence from {msa_file} already cached")
                stats['already_cached'] += 1
            else:
                # Add MSA to global cache
                try:
                    cached_path = add_msa_to_cache(cache_folder, query_sequence, msa_file_path)
                    print(f"Added {msa_file} to global cache: {os.path.basename(cached_path)}")
                    stats['newly_cached'] += 1
                except Exception as e:
                    print(f"Failed to cache MSA from {msa_file}: {e}")
                    stats['errors'] += 1
                    
        except Exception as e:
            print(f"Error processing {msa_file}: {e}")
            stats['errors'] += 1
    
    return stats

def main():
    parser = argparse.ArgumentParser(description='Update global MSA cache from MSA folder')
    parser.add_argument('MSA_FOLDER', help='Folder containing MSA files to process')
    parser.add_argument('BOLTZ_CACHE_FOLDER', help='Boltz cache folder (contains GlobalMSAs)')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Process the MSA folder
    stats = process_msa_folder(args.MSA_FOLDER, args.BOLTZ_CACHE_FOLDER)
    
    # Print summary
    print("\n=== Global MSA Cache Update Summary ===")
    print(f"Total MSA files processed: {stats['total_files']}")
    print(f"Sequences successfully extracted: {stats['sequences_extracted']}")
    print(f"Already in cache: {stats['already_cached']}")
    print(f"Newly added to cache: {stats['newly_cached']}")
    print(f"Errors: {stats['errors']}")
    
    if stats['newly_cached'] > 0:
        print(f"\nSuccessfully added {stats['newly_cached']} new MSAs to global cache")
    
    if stats['errors'] > 0:
        print(f"Warning: {stats['errors']} files had processing errors")
        sys.exit(1)
    else:
        print("Global MSA cache update completed successfully")

if __name__ == "__main__":
    main()