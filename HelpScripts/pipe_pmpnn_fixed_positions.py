#Copyright © 2024 LOCBP @ University of Zürich
#Distributed under MIT license
import argparse

parser = argparse.ArgumentParser(description='Establish residues for which to generate a sequence with ProteinMPNN. Can read from RFdiffusion datasheet, use direct selections, or apply pLDDT threshold')
parser.add_argument('JOB_FOLDER', type=str, help="Directory containing PDB files")
parser.add_argument('input_source', type=str, help = "'datasheet' for reading from datasheet, 'selection' for direct input, 'plddt' for threshold")
parser.add_argument('input_datasheet', type=str, help='Path to input datasheet (e.g., tool_results.csv) or "-" if not used')
parser.add_argument('pLDDT_thr', type=float, help="Will consider residues with pLDDT > threshold as fixed ")
parser.add_argument('FIXED', type=str, help="PyMOL selection for fixed positions or '-'")
parser.add_argument('DESIGNED', type=str, help="PyMOL selection for designed positions or '-'")
parser.add_argument('FIXED_CHAIN', type=str, help="A or B or whatever")
parser.add_argument('fixed_jsonl_file', type=str, help="output file")
parser.add_argument('sele_csv_file', type=str, help="output file with selections of fixed and mobile parts")
# Parse the arguments
args = parser.parse_args()
JOB_FOLDER=args.JOB_FOLDER
input_source=args.input_source
input_datasheet=args.input_datasheet
pLDDT_thr=args.pLDDT_thr
FIXED=args.FIXED
DESIGNED=args.DESIGNED
FIXED_CHAIN=args.FIXED_CHAIN
fixed_jsonl_file=args.fixed_jsonl_file
sele_csv_file=args.sele_csv_file

#Converts a pymol selection into an array
def sele_to_list(s):
    a = []
    if s == "": return a
    elif '+' in s:
        plus_parts = s.split('+')
        for pp in plus_parts:
            if '-' in pp:
                min,max = pp.split('-')
                for ri in range(int(min),int(max)+1):
                    a.append(ri)
            else:
                a.append(int(pp))
    else:
        if '-' in s:
            min,max = s.split('-')
            for ri in range(int(min),int(max)+1):
                a.append(ri)
        else:
            a.append(int(s))     
    return a

def list_to_sele(a):
    s = ""
    i = 0
    while i < len(a):   
        if i > 0: s += "+"
        s += f"{a[i]}"   
        #represent consecutive indeces with a dash
        if i < len(a) - 1:
            if int(a[i])+1 == int(a[i+1]):
                s += "-"
                j = i + 2
                while j < len(a):
                    if int(a[j]) != int(a[j-1])+1: break
                    j += 1
                i = j - 1
                s += f"{a[i]}" 
        i += 1        
    return s

import os
import pandas as pd
import pymol
from pymol import cmd

design_names = [d[:-4] for d in os.listdir(JOB_FOLDER) if d.endswith(".pdb")]
fixed_dict = dict()
mobile_dict = dict()

if input_source == "datasheet": #derive from datasheet
    try:
        # Read input datasheet
        df = pd.read_csv(input_datasheet)
        print(f"Reading from datasheet: {input_datasheet}")
        
        for _, row in df.iterrows():
            pdb_name = os.path.splitext(os.path.basename(row['pdb']))[0]
            if pdb_name in design_names:
                fixed_dict[pdb_name] = dict()
                mobile_dict[pdb_name] = dict()
                
                # Parse fixed and designed selections from datasheet
                if 'fixed' in df.columns and pd.notna(row['fixed']) and row['fixed'] != '':
                    fixed_dict[pdb_name][FIXED_CHAIN] = sele_to_list(str(row['fixed']))
                else:
                    fixed_dict[pdb_name][FIXED_CHAIN] = []
                    
                if 'designed' in df.columns and pd.notna(row['designed']) and row['designed'] != '':
                    mobile_dict[pdb_name][FIXED_CHAIN] = sele_to_list(str(row['designed']))
                else:
                    mobile_dict[pdb_name][FIXED_CHAIN] = []
                    
                print(f"Design: {pdb_name}, Fixed: {fixed_dict[pdb_name][FIXED_CHAIN]}, Designed: {mobile_dict[pdb_name][FIXED_CHAIN]}")
                
    except Exception as e:
        print(f"Error reading datasheet {input_datasheet}: {e}")
        print("Falling back to pLDDT threshold method...")
        input_source = "plddt"

if input_source == "selection":
    # Use direct selections provided
    FIXED = FIXED if FIXED != '-' else ''
    DESIGNED = DESIGNED if DESIGNED != '-' else ''
    
    for name in design_names:
        fixed_dict[name] = dict()
        mobile_dict[name] = dict()
        fixed_dict[name][FIXED_CHAIN] = sele_to_list(FIXED)
        mobile_dict[name][FIXED_CHAIN] = sele_to_list(DESIGNED)
        print(f"Design: {name}, Fixed: {FIXED}, Designed: {DESIGNED}")
        
elif input_source == "plddt" or input_source == "datasheet":  # datasheet fallback
    # Use pLDDT threshold method
    FIXED = FIXED if FIXED != '-' else ''
    if pLDDT_thr < 100:
        try:
            # Initialize PyMOL in headless mode (no GUI)
            pymol.pymol_argv = ['pymol', '-c']  # -q for quiet, -c for no GUI
            pymol.finish_launching()
        except Exception as e:
            print("Error while initializing pymol")
            print(str(e))
    
    for name in design_names:
        if name not in fixed_dict:  # Skip if already processed from datasheet
            fixed_dict[name] = dict()
            mobile_dict[name] = dict()
            
            if pLDDT_thr < 100:
                pdb_file = os.path.join(JOB_FOLDER,name+".pdb")
                try:
                    cmd.load(pdb_file,"prot")
                    fixed_residues = []
                    mobile_residues = []
                    atom_iterator = cmd.get_model("prot and name CA")
                    parfixed = sele_to_list(FIXED)
                    for atom in atom_iterator.atom:
                        resi = int(atom.resi)
                        if atom.b < pLDDT_thr and not resi in parfixed:
                            if not resi in mobile_residues:
                                mobile_residues.append(int(atom.resi))
                        else:
                            if not resi in fixed_residues:
                                fixed_residues.append(int(atom.resi))
                    cmd.delete("prot")
                    fixed_dict[name][FIXED_CHAIN] = fixed_residues[:]
                    mobile_dict[name][FIXED_CHAIN] = mobile_residues[:]
                    print(f"Design: {name}, pLDDT-based - Fixed: {len(fixed_residues)}, Mobile: {len(mobile_residues)}")
                except Exception as e:
                    print("Error while calculating fixed positions")
                    print(str(e))
                    fixed_dict[name][FIXED_CHAIN] = sele_to_list(FIXED)
                    mobile_dict[name][FIXED_CHAIN] = []
            else:
                fixed_dict[name][FIXED_CHAIN] = sele_to_list(FIXED)
                mobile_dict[name][FIXED_CHAIN] = sele_to_list(DESIGNED)
                print(f"Design: {name}, Selection-based - Fixed: {FIXED}, Redesigned: {DESIGNED}")

with open(fixed_jsonl_file,"w") as jsonl_file:
    #Python converts dictionaries to string having keys inside '', json only recognises ""
    jsonl_file.write(str(fixed_dict).replace("\'","\""))

with open(sele_csv_file,"w") as csv_file:
    csv_file.write("id,fixed,mobile")
    for id in fixed_dict.keys():
        fixed = list_to_sele(fixed_dict[id][FIXED_CHAIN]) if FIXED_CHAIN in fixed_dict[id].keys() else ""
        mobile = list_to_sele(mobile_dict[id][FIXED_CHAIN]) if FIXED_CHAIN in fixed_dict[id].keys() else ""
        csv_file.write("\n")
        csv_file.write(f"{id},{fixed},{mobile}")
        
print(f"Fixed positions written to: {fixed_jsonl_file}")
print(f"Selections summary written to: {sele_csv_file}")
        
#Dictionary of fixed positions looks like this
#{"design_0": {"A": [1, 2, 3, 7, 8, 9, 22, 25, 33], "B": []}, "design_1": {"A": [], "B": []}}
#Input modes:
# - datasheet: Read fixed/designed from input datasheet (e.g., from RFdiffusion or other tools)
# - selection: Use FIXED and DESIGNED parameters directly
# - plddt: Use pLDDT threshold to determine positions

if pLDDT_thr < 100:    
    cmd.quit() #MUST ALWAYS BE AT THE END OF THE SCRIPT