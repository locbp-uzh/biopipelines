# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.
"""
Script used to generate a library of compounds for Boltz 
Example from a library csv file of the kind
PRIMARY_KEY,a_b_c
a,CCC,CCCC
b,CN(=O),CN(C)
c,CCF
it will generate combinations
1,PRIM1,CCCCN(=O)CCF,...
2,PRIM2,CCCCN(=C)CCF,...
3,PRIM3,CCCCCN(=O)CCF,...
4,PRIM4,CCCCCN(=C)CCF,...
"""

import argparse
parser = argparse.ArgumentParser(description='Script used to generate a library of compounds for Boltz')
parser.add_argument('library_in_csv', type=str, help='csv library input')
parser.add_argument('primary_key', type=str, help='csv library input')
parser.add_argument('library_out_csv', type=str, help='csv library input')
args = parser.parse_args()

### Generate library
library = dict()
with open(args.library_in_csv,"r") as library_csv:
    library_csv.readline() #header
    for line in library_csv:
        line = line.strip()
        if ',' in line:
            vals = line.split(',')
            library[vals[0]] = [v for v in vals[1:] if v and not v.startswith('_')]
# remove any entries with empty keys
library = {k: v for k, v in library.items() if k}
#print('----Library-----')
#print(library)
# determine which keys are actually reachable from the primary key
no_more_branching = False
reachable_keys = [args.primary_key]
while not no_more_branching:
    no_more_branching = True
    for current_key in reachable_keys:
        for candidate_key in library.keys():
            if candidate_key in reachable_keys:
                continue
            # get the list of SMILES options for this key
            vals = library[current_key]
            if isinstance(vals, str):
                vals = [vals]
            # if the candidate_key appears in any of those SMILES, it’s reachable
            for val in vals:
                if candidate_key in val:
                    no_more_branching = False
                    reachable_keys.append(candidate_key)
                    break

print('-----Reachable keys-----')
print(' | '.join(reachable_keys))

# now prune the library—and recompute the master key list—so only reachable entries remain
library = {k: v for k, v in library.items() if k in reachable_keys}
library_keys = list(library.keys())

branch_keys = [key for key, options in library.items() if len(options) > 1] # Identify branching keys (i.e. those with more than one option)
print('-----Branching keys-----')
print(' | '.join(branch_keys))

# Expand each base compound from the primary key
final_compounds = []
for base in library[args.primary_key]:
    final_compounds.append({'SMILES':base,'BRANCHING':{}})
no_new_branching = False
while not no_new_branching:
    no_new_branching = True
    updated_compounds = []
    # Process every compound in current list
    for compound in final_compounds:
        key_found = False
        # Check for every library key in the current compound's SMILES
        for key in reachable_keys:
            if key in compound['SMILES']:
                key_found = True
                no_new_branching = False
                # For every possible substitution for the key, create a new compound
                for option in library[key]:
                    new_smiles = compound['SMILES'].replace(key, option, 1)
                    new_branching = compound['BRANCHING'].copy()
                    new_branching[key] = option
                    updated_compounds.append({'SMILES': new_smiles, 'BRANCHING': new_branching})
                # Process one key per compound per iteration
                break
        if not key_found:
            updated_compounds.append(compound)
    final_compounds = updated_compounds
num_compounds = len(final_compounds)
print(f"Number of compounds in the library: {num_compounds}")

for u_l_n in range(len(final_compounds)):
    ### Visualize molecule and atom indices for sanity check
    characters = 4
    if num_compounds>9: characters=3
    if num_compounds>99: characters=2
    if num_compounds>999: characters=1
    if num_compounds>99999: characters=0
    u_l_n_str = str(u_l_n)
    n0 = 5-characters-len(u_l_n_str)
    zeros_str = '0'*n0
    compound_name = args.primary_key if num_compounds == 1 else args.primary_key[:characters]+ zeros_str + u_l_n_str
    final_compounds[u_l_n]['ID']=str(u_l_n)
    final_compounds[u_l_n]['NAME']=compound_name

all_branch_keys = set()
for compound in final_compounds:
    record = compound['BRANCHING']
    all_branch_keys.update(record.keys())
all_branch_keys = sorted(list(all_branch_keys))

header = ["INDEX", "NAME", "SMILES"] + [x.replace('*','') for x in all_branch_keys]
with open(args.library_out_csv, 'w', newline='') as out:
    out.write(','.join(header))
    for compound in final_compounds:
        out.write('\n')
        row = [compound['ID'], compound['NAME'], compound['SMILES']]
        for key in all_branch_keys:
            row.append(compound['BRANCHING'].get(key, ""))
        out.write(','.join(row))