#Copyright © 2024 LOCBP @ University of Zürich
#Distributed under MIT license
import argparse

parser = argparse.ArgumentParser(description='Makes csv files from ProteinMPNN.fa output ')
parser.add_argument('FA_FOLDER', type=str)
parser.add_argument('queries_csv_file', type=str)
parser.add_argument('queries_fasta_file', type=str)
parser.add_argument('-d', '--duplicates', action='store_true',
                    help='Allow duplicate sequences in output')

# Parse the arguments
args = parser.parse_args()

import os
fa_files = os.listdir(args.FA_FOLDER)
old_sequences = []
for fa in fa_files:
    if fa.endswith(".fa"):
        data = []
        with open(os.path.join(args.FA_FOLDER,fa),"r") as file:
            lines = [line.strip() for line in file.readlines()]
            #Starts from 2 to skip the original protein
            for i in range(2,len(lines),2):
                seq_data = dict()
                seq_data["id"] = "" #put here just to preserve the order
                seq_data["sequence"] = lines[i+1]
                params = lines[i][1:].split(", ")
                for p in params:
                    if not '=' in p: 
                        continue
                    key,value = p.split("=")
                    if key == 'id': key = 'sample' #ligandmpnn is structured differently from proteinmpnn
                    seq_data[key] = value
                seq_data["id"] = fa[:-3] + "_" + seq_data["sample"]
                if seq_data["sequence"] in old_sequences:
                    if args.duplicates:
                        data.append(seq_data)
                        old_sequences.append(seq_data["sequence"])
                    else:
                        print(f"Skipped duplicate: {seq_data['id']}")
                else:
                    data.append(seq_data)
                    old_sequences.append(seq_data["sequence"])
        if len(data) == 0: continue
        #Update CSV
        old_data = []
        if os.path.exists(args.queries_csv_file):
            with open(args.queries_csv_file,"r") as old_queries:
                old_data = old_queries.readlines()
        with open(args.queries_csv_file,"w") as queries:
            keys = list(data[0].keys())
            if len(old_data) != 0:
                queries.writelines(old_data)
            #Header is already part of old_data
            else:
                for i,k in enumerate(keys):
                    if i>0: queries.write(',')
                    queries.write(k)
            for seq_data in data:
                queries.write("\n")
                for i,k in enumerate(keys):
                    if i>0: queries.write(',')
                    queries.write(seq_data[k])
        #Update FASTA
        old_data = []
        if os.path.exists(args.queries_fasta_file):
            with open(args.queries_fasta_file,"r") as old_queries:
                old_data = old_queries.readlines()
        with open(args.queries_fasta_file,"w") as queries:
            if len(old_data) != 0:
                queries.writelines(old_data)
                queries.write("\n")
            for s_d,seq_data in enumerate(data):
                if s_d>0: queries.write('\n')
                queries.write(">")
                queries.write(seq_data["id"])
                queries.write('\n')
                queries.write(seq_data["sequence"])