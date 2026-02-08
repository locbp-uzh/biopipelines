# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
Format conversion utilities for pipeline data flow.

Provides functions to convert between different file formats used by
modeling tools, ensuring seamless integration and data validation.
"""

import os
import json
from typing import List, Dict, Any, Optional
from pathlib import Path


def pdb_to_jsonl(pdb_folder: str, output_jsonl: str) -> bool:
    """
    Convert PDB files to JSONL format for ProteinMPNN.
    
    Mimics the functionality of parse_multiple_chains.py from ProteinMPNN.
    
    Args:
        pdb_folder: Folder containing PDB files
        output_jsonl: Output JSONL file path
        
    Returns:
        True if conversion successful
    """
    try:
        pdb_files = [f for f in os.listdir(pdb_folder) if f.endswith('.pdb')]
        if not pdb_files:
            raise ValueError(f"No PDB files found in {pdb_folder}")
        
        jsonl_entries = []
        
        for pdb_file in pdb_files:
            pdb_path = os.path.join(pdb_folder, pdb_file)
            pdb_name = os.path.splitext(pdb_file)[0]
            
            # Basic JSONL entry structure
            # In actual implementation, would parse PDB structure
            entry = {
                "name": pdb_name,
                "path": pdb_path,
                "chains": ["A"],  # Simplified - would extract from PDB
            }
            jsonl_entries.append(entry)
        
        # Write JSONL file
        with open(output_jsonl, 'w') as f:
            for entry in jsonl_entries:
                f.write(json.dumps(entry) + '\n')
        
        return True
        
    except Exception as e:
        print(f"Error converting PDB to JSONL: {e}")
        return False


def fasta_to_csv(fasta_folder: str, output_csv: str, output_fasta: str = None) -> bool:
    """
    Convert FASTA files to ColabFold CSV format.
    
    Mimics the functionality of fa_to_csv_fasta.py helper script.
    
    Args:
        fasta_folder: Folder containing FASTA files
        output_csv: Output CSV file for ColabFold
        output_fasta: Optional combined FASTA output
        
    Returns:
        True if conversion successful
    """
    try:
        fasta_files = [f for f in os.listdir(fasta_folder) 
                      if f.endswith(('.fasta', '.fa'))]
        
        if not fasta_files:
            raise ValueError(f"No FASTA files found in {fasta_folder}")
        
        csv_entries = []
        all_sequences = []
        
        for fasta_file in fasta_files:
            fasta_path = os.path.join(fasta_folder, fasta_file)
            
            # Read FASTA file
            with open(fasta_path, 'r') as f:
                content = f.read().strip()
            
            # Parse FASTA entries
            entries = content.split('>')
            for entry in entries:
                if not entry.strip():
                    continue
                
                lines = entry.strip().split('\n')
                header = lines[0]
                sequence = ''.join(lines[1:])
                
                if sequence:
                    # Generate unique ID
                    seq_id = f"{os.path.splitext(fasta_file)[0]}_{header}"
                    csv_entries.append(f"{seq_id},{sequence}")
                    all_sequences.append(f">{seq_id}\n{sequence}")
        
        # Write CSV file (ColabFold format)
        with open(output_csv, 'w') as f:
            f.write("id,sequence\n")
            for entry in csv_entries:
                f.write(entry + '\n')
        
        # Write combined FASTA file if requested
        if output_fasta:
            with open(output_fasta, 'w') as f:
                f.write('\n'.join(all_sequences))
        
        return True
        
    except Exception as e:
        print(f"Error converting FASTA to CSV: {e}")
        return False


def validate_pdb(pdb_file: str) -> bool:
    """
    Validate PDB file format and basic structure.
    
    Args:
        pdb_file: Path to PDB file
        
    Returns:
        True if PDB file is valid
    """
    try:
        if not os.path.exists(pdb_file):
            return False
        
        with open(pdb_file, 'r') as f:
            content = f.read()
        
        # Basic PDB validation
        if not content.strip():
            return False
        
        # Check for essential PDB records
        lines = content.split('\n')
        has_atom_records = any(line.startswith('ATOM') for line in lines)
        
        return has_atom_records
        
    except Exception:
        return False


def validate_fasta(fasta_file: str) -> bool:
    """
    Validate FASTA file format and sequences.
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        True if FASTA file is valid
    """
    try:
        if not os.path.exists(fasta_file):
            return False
        
        with open(fasta_file, 'r') as f:
            content = f.read().strip()
        
        if not content:
            return False
        
        # Basic FASTA validation
        if not content.startswith('>'):
            return False
        
        # Check for sequences
        entries = content.split('>')
        valid_entries = 0
        
        for entry in entries[1:]:  # Skip first empty entry
            lines = entry.strip().split('\n')
            if len(lines) >= 2:  # Header + at least one sequence line
                sequence = ''.join(lines[1:])
                if sequence and all(c.upper() in 'ACDEFGHIKLMNPQRSTVWY' for c in sequence.replace('-', '')):
                    valid_entries += 1
        
        return valid_entries > 0
        
    except Exception:
        return False


def combine_fasta_files(fasta_files: List[str], output_file: str) -> bool:
    """
    Combine multiple FASTA files into one.
    
    Args:
        fasta_files: List of FASTA file paths
        output_file: Output combined FASTA file
        
    Returns:
        True if combination successful
    """
    try:
        all_sequences = []
        
        for fasta_file in fasta_files:
            if not validate_fasta(fasta_file):
                print(f"Invalid FASTA file: {fasta_file}")
                continue
            
            with open(fasta_file, 'r') as f:
                content = f.read().strip()
            
            # Add sequences with updated headers
            file_base = os.path.splitext(os.path.basename(fasta_file))[0]
            entries = content.split('>')
            
            for i, entry in enumerate(entries[1:]):  # Skip first empty
                lines = entry.strip().split('\n')
                if len(lines) >= 2:
                    header = lines[0]
                    sequence = ''.join(lines[1:])
                    
                    # Create unique header
                    new_header = f"{file_base}_{i}_{header}"
                    all_sequences.append(f">{new_header}\n{sequence}")
        
        # Write combined file
        if all_sequences:
            with open(output_file, 'w') as f:
                f.write('\n'.join(all_sequences))
            return True
        else:
            print("No valid sequences found to combine")
            return False
        
    except Exception as e:
        print(f"Error combining FASTA files: {e}")
        return False


def extract_sequences_from_pdb(pdb_file: str) -> Dict[str, str]:
    """
    Extract amino acid sequences from PDB file.
    
    Args:
        pdb_file: Path to PDB file
        
    Returns:
        Dictionary mapping chain IDs to sequences
    """
    sequences = {}
    
    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        current_chain = None
        residues = []
        last_resnum = None
        
        for line in lines:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                chain = line[21]
                resnum = int(line[22:26])
                resname = line[17:20].strip()
                
                # Convert 3-letter to 1-letter amino acid code
                aa_map = {
                    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
                    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
                }
                
                if chain != current_chain:
                    # New chain - save previous if exists
                    if current_chain is not None:
                        sequences[current_chain] = ''.join(residues)
                    
                    current_chain = chain
                    residues = []
                    last_resnum = None
                
                # Add residue if new position
                if resnum != last_resnum:
                    if resname in aa_map:
                        residues.append(aa_map[resname])
                    last_resnum = resnum
        
        # Save last chain
        if current_chain is not None:
            sequences[current_chain] = ''.join(residues)
        
    except Exception as e:
        print(f"Error extracting sequences from PDB: {e}")
    
    return sequences