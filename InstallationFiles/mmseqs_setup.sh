#!/bin/bash

source /apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/mamba-23.11.0-0-334ztq7i4mzu762ew2x3kbbrrorhe6eg/etc/profile.d/conda.sh
module load mamba

env_exists=$(conda env list | grep 'MMseqs2')

if [ -z "$env_exists" ]; then
    echo "Creating conda environment..."
    conda create -n MMseqs2
else
    echo "MMseqs2 found"
fi
echo
conda activate MMseqs2 #Common to all models, contains pymol as well

if [ "$CONDA_DEFAULT_ENV" = "MMseqs2" ]; then
  echo "MMseqs2 activated"
else
  echo "Could not activate MMseqs2 environment. Check location of conda is /apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/mamba-23.11.0-0-334ztq7i4mzu762ew2x3kbbrrorhe6eg/etc/profile.d/conda.sh by typing: which conda (from cluster shell). Contact administrator otherwise"
  exit
fi

conda_list_output=$(conda list)
if echo "$conda_list_output" | grep -q "mmseqs2"; then
  echo "MMseqs2 is already installed"
else
  echo "Installing MMseqs2..."
  conda install -c conda-forge -c bioconda mmseqs2
fi

echo
echo "Setting up databases..."
cd /home/$USER/scratch
if [ ! -d "MMseqs2DBs" ]; then
    mkdir MMseqs2DBs
fi
cd MMseqs2DBs

# Set MMSEQS_NO_INDEX to skip the index creation step (not useful for colabfold_search in most cases)
ARIA_NUM_CONN=8

PDB_SERVER="${2:-"rsync.wwpdb.org::ftp"}"
PDB_PORT="${3:-"33444"}"
# do initial download of the PDB through aws?
# still syncs latest structures through rsync
PDB_AWS_DOWNLOAD="${4:-}"
PDB_AWS_SNAPSHOT="20240101"

UNIREF30DB="uniref30_2302"
MMSEQS_NO_INDEX=${MMSEQS_NO_INDEX:-}

# set which commits to use
MMSEQS_COMMIT=${1:-4589151554eb83a70ff0c4d04d21b83cabc203e4}
BACKEND_COMMIT=${2:-14e087560f309f989a5e1feb54fd1f9c988076d5}

# check that the correct mmseqs commit is there
downloadFile() {
    URL="$1"
    OUTPUT="$2"
    curl -L -o "$OUTPUT" "$URL"
}

# Make MMseqs2 merge the databases to avoid spamming the folder with files
export MMSEQS_FORCE_MERGE=1

#curl -s -o- https://mmseqs.com/archive/${1:-4589151554eb83a70ff0c4d04d21b83cabc203e4}/mmseqs-linux-avx2.tar.gz | tar -xzf - -C /scratch/gquarg/ mmseqs

if [ ! -f UNIREF30_READY ]; then
  downloadFile "https://wwwuser.gwdg.de/~compbiol/colabfold/${UNIREF30DB}.tar.gz" "${UNIREF30DB}.tar.gz"
  tar xzvf "${UNIREF30DB}.tar.gz"
  mmseqs tsv2exprofiledb uniref30_2302 uniref30_2302
  MMSEQS_NO_INDEX=1 MMSEQS_FORCE_MERGE=1 mmseqs tsv2exprofiledb uniref30_2302 uniref30_2302
  if [ -z "$MMSEQS_NO_INDEX" ]; then
    mmseqs createindex "${UNIREF30DB}_db" tmp1 --remove-tmp-files 1
  fi
  if [ -e ${UNIREF30DB}_db_mapping ]; then
    ln -sf ${UNIREF30DB}_db_mapping ${UNIREF30DB}_db.idx_mapping
  fi
  if [ -e ${UNIREF30DB}_db_taxonomy ]; then
    ln -sf ${UNIREF30DB}_db_taxonomy ${UNIREF30DB}_db.idx_taxonomy
  fi
  touch UNIREF30_READY
fi

if [ ! -f COLABDB_READY ]; then
  downloadFile "https://wwwuser.gwdg.de/~compbiol/colabfold/colabfold_envdb_202108.tar.gz" "colabfold_envdb_202108.tar.gz"
  tar xzvf "colabfold_envdb_202108.tar.gz"
  mmseqs tsv2exprofiledb "colabfold_envdb_202108" "colabfold_envdb_202108_db"
  # TODO: split memory value for createindex?
  if [ -z "$MMSEQS_NO_INDEX" ]; then
    mmseqs createindex "colabfold_envdb_202108_db" tmp2 --remove-tmp-files 1
  fi
  touch COLABDB_READY
fi

if [ ! -f PDB_READY ]; then
  downloadFile "https://wwwuser.gwdg.de/~compbiol/colabfold/pdb100_230517.fasta.gz" "pdb100_230517.fasta.gz"
  mmseqs createdb pdb100_230517.fasta.gz pdb100_230517
  if [ -z "$MMSEQS_NO_INDEX" ]; then
    mmseqs createindex pdb100_230517 tmp3 --remove-tmp-files 1
  fi
  touch PDB_READY
fi

if [ ! -f PDB100_READY ]; then
  downloadFile "https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pdb100_foldseek_230517.tar.gz" "pdb100_foldseek_230517.tar.gz"
  tar xzvf pdb100_foldseek_230517.tar.gz pdb100_a3m.ffdata pdb100_a3m.ffindex
  touch PDB100_READY
fi

if [ ! -f PDB_MMCIF_READY ]; then
  mkdir -p pdb/divided
  mkdir -p pdb/obsolete
  if [ -n "${PDB_AWS_DOWNLOAD}" ]; then
    aws s3 cp --no-sign-request --recursive s3://pdbsnapshots/${PDB_AWS_SNAPSHOT}/pub/pdb/data/structures/divided/mmCIF/ pdb/divided/
    aws s3 cp --no-sign-request --recursive s3://pdbsnapshots/${PDB_AWS_SNAPSHOT}/pub/pdb/data/structures/obsolete/mmCIF/ pdb/obsolete/
  fi
  rsync -rlpt -v -z --delete --port=${PDB_PORT} ${PDB_SERVER}/data/structures/divided/mmCIF/ pdb/divided
  rsync -rlpt -v -z --delete --port=${PDB_PORT} ${PDB_SERVER}/data/structures/obsolete/mmCIF/ pdb/obsolete
  touch PDB_MMCIF_READY
fi