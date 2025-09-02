"""
Script used by Boltz to generate a folder of yaml config files given a base yaml
"""
import os
### Arguments
import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument('base_yaml_file', type=str, help='full path to txt file containing basic configuration (txt or yaml file)')
parser.add_argument('library_csv_file', type=str, help='full path to csv file containing columns ID NAME SMILES ...')
parser.add_argument('primary_key', type=str, help='primary key to replace in the yaml')
parser.add_argument('id',type=str,choices=['CCD','SMILES'],default='SMILES')
parser.add_argument('CONFIG_FOLDER', type=str, help='folder where to store smiles and pictures')
args = parser.parse_args()

with open(args.base_yaml_file,"r") as base_yaml:
    base_config = base_yaml.read()

### Generate library
library = dict()
with open(args.library_csv_file,"r") as library_csv:
    library_csv.readline()
    for line in library_csv:
        line = line.strip()
        if ',' in line:
            vals = line.split(',') #ID,NAME,SMILES;BRANCHING
            library[vals[1]] = vals[2] #NAME=SMILES
NAMEs = list(library.keys())

for NAME in NAMEs:
    compound_name_config_file = os.path.join(args.CONFIG_FOLDER,f"{NAME}.yaml")
    #cmp_config=base_config.replace(f"/{primary_lib_key}/","***FOLDER***")
    #mp_config=cmp_config.replace(primary_lib_key,compound_name)
    #cmp_config=cmp_config.replace("***FOLDER***",f"/{primary_lib_key}/")
    with open(compound_name_config_file,'w') as config:
        if args.id == 'CCD': config.write(base_config.replace('\''+args.primary_key+'\'','\''+NAME+'\''))
        elif args.id == 'SMILES':
            config.write(base_config.replace('\''+args.primary_key+'\'','\''+library[NAME]+'\''))

