# Copyright (c) 2026 Gianluca Quargnali @ LOCBP (https://www.locbp.com/) University of Zurich Switzerland
#
# This software is freely available for use, modification, and redistribution.
# If you use this software or any derivative of it in published work,
# you must cite the original author and this repository.

"""
Script used by Boltz to rearrange all the outputs in one place and give a table for ranking

"""

import argparse

parser = argparse.ArgumentParser(description='Slices the queries for alphafold refolding of the best ones')
parser.add_argument('PREDICTION_FOLDER', type=str)
parser.add_argument('library_scores_csv_file', type=str)
parser.add_argument('RESULTS_FOLDER', type=str)
parser.add_argument('metric', type=str, help = "see JSON structure/scores explanation")
parser.add_argument('pymol_best_pse', type=int, help = "Create a pymol session contaning the N best models aligned with the original protein and colored by pLDDT")
parser.add_argument('pymol_pse_file', type=str, help = "Path to pymol session to be created")
parser.add_argument('highlight', type=str, help = "")

# Parse the arguments
args = parser.parse_args()

import os
import json
import shutil
import csv

# ----- Helper function to flatten nested confidence JSON -----
def flatten_confidence(conf, parent_key=""):
    """
    Flattens the confidence dictionary. Nested keys are concatenated.
    For a dict value, keys are appended with ':' (or '-' for the second-level nested dict)
    e.g. {'chains_ptm': {'0': 0.85, '1': 0.83}} becomes {'chains_ptm:0': 0.85, 'chains_ptm:1': 0.83}
         {'pair_chains_iptm': {'0': {'0': 0.85, '1': 0.81}, '1': {'0': 0.82, '1': 0.83}}}
             becomes {'pair_chains_iptm:0-0': 0.85, 'pair_chains_iptm:0-1': 0.81, 'pair_chains_iptm:1-0': 0.82, 'pair_chains_iptm:1-1': 0.83}
    """
    items = {}
    for key, value in conf.items():
        new_key = f"{parent_key}:{key}" if parent_key else key
        if isinstance(value, dict):
            # For second-level dictionaries, use a different delimiter (dash) if needed.
            for subkey, subvalue in value.items():
                if isinstance(subvalue, dict):
                    for subsubkey, subsubvalue in subvalue.items():
                        final_key = f"{new_key}-{subkey}-{subsubkey}"
                        items[final_key] = subsubvalue
                else:
                    final_key = f"{new_key}:{subkey}"
                    items[final_key] = subvalue
        else:
            items[new_key] = value
    return items

PREDICTION_FOLDER = args.PREDICTION_FOLDER  
library_scores_csv_file = args.library_scores_csv_file 
results_folder = args.RESULTS_FOLDER
predictions_folder = os.path.join(PREDICTION_FOLDER, "predictions")

confidence_data = {}  # dictionary mapping NAME -> flattened confidence dict
affinity_data = {}  # dictionary mapping NAME -> flattened affinity dict
structural_files = {} # dictionary mapping NAME -> pdb or cif file

# ----- Process each prediction subfolder -----
for folder in os.listdir(predictions_folder):
    folder_path = os.path.join(predictions_folder, folder)
    if os.path.isdir(folder_path):
        NAME = folder
        pdb_filename = f"{folder}_model_0.pdb"
        pdb_filepath = os.path.join(folder_path, pdb_filename)
        conf_filename = f"confidence_{folder}_model_0.json"
        conf_filepath = os.path.join(folder_path, conf_filename)
        aff_filename = f"affinity_{folder}.json"
        aff_filepath = os.path.join(folder_path, aff_filename)
        
        target_pdb_filepath = os.path.join(results_folder, f"{NAME}.pdb")
        try:
            if os.path.exists(pdb_filepath):
                shutil.move(pdb_filepath, target_pdb_filepath)
            else:
                pdb_filepath = pdb_filepath[:-3]+"cif"
                target_pdb_filepath = target_pdb_filepath[:-3]+"cif"
                shutil.move(pdb_filepath, target_pdb_filepath)
            print(f"Moved {pdb_filename} to {target_pdb_filepath}")
        except Exception as e:
            print(f"Error moving {pdb_filepath}: {e}")
        
        # Initialize data dictionaries for this compound
        confidence_data[NAME] = {}
        affinity_data[NAME] = {}
        
        # Load confidence data
        try:
            with open(conf_filepath, "r") as f:
                conf = json.load(f)
            flattened_conf = flatten_confidence(conf)
            confidence_data[NAME] = flattened_conf
            print(f"Loaded confidence for {NAME}: {len(flattened_conf)} keys")
        except Exception as e:
            print(f"Error reading confidence file {conf_filepath}: {e}")
            confidence_data[NAME] = {}
        
        # Load affinity data
        try:
            if os.path.exists(aff_filepath):
                with open(aff_filepath, "r") as f:
                    aff = json.load(f)
                flattened_aff = flatten_confidence(aff)
                affinity_data[NAME] = flattened_aff
                print(f"Loaded affinity for {NAME}: {len(flattened_aff)} keys")
            else:
                print(f"No affinity file found for {NAME}")
                affinity_data[NAME] = {}
        except Exception as e:
            print(f"Error reading affinity file {aff_filepath}: {e}")
            affinity_data[NAME] = {}
        
        # Store structural file path
        structural_files[NAME] = target_pdb_filepath

# Check if we have any data
print(f"> Total compounds processed: {len(confidence_data)}")
print(f"> Compounds with confidence data: {[name for name, data in confidence_data.items() if data]}")
print(f"> Compounds with affinity data: {[name for name, data in affinity_data.items() if data]}")

# Early exit check - moved after data processing
if library_scores_csv_file == '-':
    print("Library scores file set to '-', exiting without CSV creation")
    exit()

current_headers = []
existing_data = []

if not os.path.exists(library_scores_csv_file):
    current_headers = ["NAME"]
else:
    try:
        with open(library_scores_csv_file, 'r', newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile)
            rows = list(reader)
            if rows:
                current_headers = rows[0]
                existing_data = rows[1:]  # All rows except header
            else:
                current_headers = ["NAME"]
    except Exception as e:
        print(f"Error opening {library_scores_csv_file}: {e}")
        current_headers = ["NAME"]
        existing_data = []

if "NAME" not in current_headers:
    print("Library scores file must have a column 'NAME' for matching the compounds.")
    exit(1)

# Collect all unique keys from confidence data
new_conf_keys = set()
for conf_dict in confidence_data.values():
    new_conf_keys.update(conf_dict.keys())
new_conf_keys = sorted(list(new_conf_keys))

# Add confidence keys to headers
for key in new_conf_keys:
    if key not in current_headers:
        current_headers.append(key)

# Collect all unique keys from affinity data
new_aff_keys = set()
for aff_dict in affinity_data.values():
    new_aff_keys.update(aff_dict.keys())
new_aff_keys = sorted(list(new_aff_keys))

# Add affinity keys to headers
for key in new_aff_keys:
    if key not in current_headers:
        current_headers.append(key)

header_to_col = {header: idx for idx, header in enumerate(current_headers)}
name_col_idx = header_to_col["NAME"]

# Debug: show headers and existing rows before merging
#print(f"DEBUG: Headers: {current_headers}")
#print(f"DEBUG: Existing_data before merging: {len(existing_data)} rows")
#if existing_data:
#    print(f"DEBUG: Existing compound names: {[row[name_col_idx] for row in existing_data if len(row) > name_col_idx]}")

# If no existing rows, initialize one row per compound name
if not existing_data:
    all_compound_names = set(confidence_data.keys()) | set(affinity_data.keys())
    for name in sorted(all_compound_names):
        new_row = [''] * len(current_headers)
        new_row[name_col_idx] = name
        existing_data.append(new_row)

#print(f"DEBUG: Existing_data after initialization: {len(existing_data)} rows")
#print(f"DEBUG: Compound names after init: {[row[name_col_idx] for row in existing_data if len(row) > name_col_idx]}")

updated_data = []
for row in existing_data:
    # Extend row to match header length if needed
    while len(row) < len(current_headers):
        row.append('')
    
    # Get compound name from this row
    if len(row) > name_col_idx:
        compound_name = row[name_col_idx]
        #print(f"DEBUG: Processing compound: {compound_name}")
        
        # Update row with confidence data if available
        if compound_name in confidence_data and confidence_data[compound_name]:
            conf_flat = confidence_data[compound_name]
            #print(f"DEBUG: Found confidence data for {compound_name}: {list(conf_flat.keys())}")
            for key, value in conf_flat.items():
                if key in header_to_col:
                    col_idx = header_to_col[key]
                    row[col_idx] = str(value)
                    #print(f"DEBUG: Set {key} = {value} for {compound_name}")
        #else:
        #    print(f"DEBUG: No confidence data for {compound_name}")
            
        # Update row with affinity data if available
        if compound_name in affinity_data and affinity_data[compound_name]:
            aff_flat = affinity_data[compound_name]
            #print(f"DEBUG: Found affinity data for {compound_name}: {list(aff_flat.keys())}")
            for key, value in aff_flat.items():
                if key in header_to_col:
                    col_idx = header_to_col[key]
                    row[col_idx] = str(value)
                    #print(f"DEBUG: Set {key} = {value} for {compound_name}")
        #else:
        #    print(f"DEBUG: No affinity data for {compound_name}")
    
    updated_data.append(row)

# Write the updated CSV
with open(library_scores_csv_file, 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(current_headers)  # Write header
    writer.writerows(updated_data)    # Write all data rows

print(f"Updated library scores file saved as {library_scores_csv_file}")

scores_explanation="""
"confidence_score": Aggregated score used to sort the predictions, corresponds to 0.8 * complex_plddt + 0.2 * iptm (ptm for single chains)
"ptm": # Predicted TM score for the complex
"iptm": Predicted TM score when aggregating at the interfaces
"ligand_iptm": ipTM but only aggregating at protein-ligand interfaces
"protein_iptm": ipTM but only aggregating at protein-protein interfaces
"complex_plddt": Average pLDDT score for the complex
"complex_iplddt": Average pLDDT score when upweighting interface tokens
"complex_pde": Average PDE score for the complex
"complex_ipde": Average PDE score when aggregating at interfaces  
"chains_ptm": Predicted TM score within each chain
"pair_chains_iptm": Predicted (interface) TM score between each pair of chains
"""
with open(os.path.join(results_folder,"scores_info.txt"),'w') as scores_exp:
    scores_exp.write(scores_explanation)

"""PSE CREATION"""

metric = args.metric
pymol_best_pse = args.pymol_best_pse #int
pymol_pse_file = args.pymol_pse_file

if pymol_best_pse > 0:
    import pymol
    from pymol import cmd
    try:
        # Initialize PyMOL in headless mode (no GUI)
        pymol.pymol_argv = ['pymol', '-c']  # -q for quiet, -c for no GUI
        pymol.finish_launching()
        #Some settings for the session to have good pictures
        cmd.do("show cartoon")
        cmd.set("seq_view", 1)
        cmd.set("cartoon_gap_cutoff", 0)
        cmd.set("sphere_scale", 0.2)
        cmd.set("ray_trace_mode", 1)
        cmd.set("ray_shadows", 0)
        cmd.set("spec_reflect", 0)
        cmd.set("ray_trace_frames", 1)
        cmd.set("ray_trace_color", "gray20")
    except Exception as e:
        print("Error while initializing pymol")
        print(str(e))
 
    blue_rgb = [0,76,202]
    blue = []
    for c in blue_rgb:
        blue.append(c/255.0)
    lightblue_rgb = [73, 196, 238]
    lightblue = []
    for c in lightblue_rgb:
        lightblue.append(c/255.0)
    yellow_rgb = [255, 213, 57]
    yellow = []
    for c in yellow_rgb:
        yellow.append(c/255.0)
    orange_rgb = [255, 113, 67]
    orange = []
    for c in orange_rgb:
        orange.append(c/255.0)     
    cmd.set_color('blue_plddt', blue)
    cmd.set_color('lightblue_plddt', lightblue)
    cmd.set_color('yellow_plddt', yellow)
    cmd.set_color('orange_plddt', orange)

    def plddt(selection="all"):    
        #select and color blue
        blue_upper = 1.0
        blue_lower = 0.9
        blue_sel_str = selection + " & ! b < " + str(blue_lower) + " & ! b > " + str(blue_upper)
        cmd.color('blue_plddt', blue_sel_str)
        #select and color lightblue
        lightblue_upper = 0.9
        lightblue_lower = 0.7
        lightblue_sel_str = selection + " & ! b < " + str(lightblue_lower) + " & ! b > " + str(lightblue_upper)
        cmd.color('lightblue_plddt', lightblue_sel_str)
        #select and color yellow
        yellow_upper = 0.7
        yellow_lower = 0.5
        yellow_sel_str = selection + " & ! b < " + str(yellow_lower) + " & ! b > " + str(yellow_upper)
        cmd.color('yellow_plddt', yellow_sel_str)
        #select and color orange
        orange_upper = 0.5
        orange_lower = 0.0
        orange_sel_str = selection + " & ! b < " + str(orange_lower) + " & ! b > " + str(orange_upper)
        cmd.color('orange_plddt', orange_sel_str)
    cmd.extend('plddt', plddt)

    def near_selection(source_sel, target_sel, sele_name, cutoff="4.0"):
        cutoff = float(cutoff)
        """
        Selects residues from `source_sel` that have at least one atom within `cutoff` Å
        of any atom in `target_sel`, based on atom-by-atom distance.
        """
        # Helper to build a unique key for a residue (using segi, chain, resi).
        def res_key(atom):
            return (atom.segi if atom.segi.strip() else None, atom.chain, atom.resi)

        # Get target atoms and build their coordinate list and residue keys.
        target_atoms = cmd.get_model(target_sel).atom
        target_coords = [a.coord for a in target_atoms]
        target_keys = { res_key(a) for a in target_atoms }

        # Function to decide if an atom is within the cutoff of any target atom.
        def is_within_cutoff(atom):
            ax, ay, az = atom.coord
            for tx, ty, tz in target_coords:
                dx = ax - tx
                dy = ay - ty
                dz = az - tz
                if (dx*dx + dy*dy + dz*dz)**0.5 < cutoff:
                    return True
            return False

        # Get atoms from the source selection.
        source_atoms = cmd.get_model(source_sel).atom
        result_keys = set()

        for atom in source_atoms:
            key = res_key(atom)
            # Skip if this residue is part of the target.
            if key in target_keys:
                continue
            if is_within_cutoff(atom):
                result_keys.add(key)

        # Construct the selection string for whole residues.
        sele_parts = []
        for segi, chain, resi in result_keys:
            if segi:
                sele_parts.append(f"(segi {segi} and chain {chain} and resi {resi})")
            else:
                sele_parts.append(f"(chain {chain} and resi {resi})")

        if sele_parts:
            cmd.select(sele_name, " or ".join(sele_parts))
        else:
            cmd.select(sele_name, "none")
    cmd.extend("near_selection",near_selection)
    
    N_best = min(pymol_best_pse,len(list(confidence_data.keys())))
    metric_array = [(N,confidence_data[N][metric]) for N in confidence_data.keys()]
    sorted_metric_array = sorted(metric_array,key=lambda x: float(x[1]),reverse=True)
    best_names = [x[0] for x in sorted_metric_array[:N_best]]

    cmd.load(structural_files[best_names[0]],best_names[0])
    for BEST_NAME in best_names[1:]:
        cmd.load(structural_files[BEST_NAME],BEST_NAME)
        cmd.align(BEST_NAME,best_names[0])

    cmd.do(f"plddt")
    
    for BEST_NAME in best_names:
        cmd.do(f"show sticks, {BEST_NAME} and resn {BEST_NAME} and not name N+CA+C+O")
        cmd.do(f"color wheat, {BEST_NAME} and resn {BEST_NAME}")
    
    if args.highlight=='yes':
        for BEST_NAME in best_names:
            sele_name = f"{BEST_NAME}_contacts"
            cmd.do(f"near {BEST_NAME},{BEST_NAME} and resn {BEST_NAME}, {sele_name}")
            cmd.do(f"show sticks, {sele_name} and not name N+CA+C+O")
            cmd.do(f"delete {sele_name}")

    cmd.do(f"color atomic, (not elem C)")

    cmd.save(pymol_pse_file)

    print(">PSE file succesfully created")