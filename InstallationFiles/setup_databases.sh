#!/bin/bash -ex
# Setup everything for using mmseqs locally
# Set MMSEQS_NO_INDEX to skip the index creation step (not useful for colabfold_search in most cases)
ARIA_NUM_CONN=8
WORKDIR="${1:-$(pwd)}"

PDB_SERVER="${2:-"rsync.wwpdb.org::ftp"}"
PDB_PORT="${3:-"33444"}"
# do initial download of the PDB through aws?
# still syncs latest structures through rsync
PDB_AWS_DOWNLOAD="${4:-}"
PDB_AWS_SNAPSHOT="20240101"

UNIREF30DB="uniref30_2302"
MMSEQS_NO_INDEX=${MMSEQS_NO_INDEX:-}
cd "${WORKDIR}"


# set which commits to use
MMSEQS_COMMIT=${1:-4589151554eb83a70ff0c4d04d21b83cabc203e4}
BACKEND_COMMIT=${2:-14e087560f309f989a5e1feb54fd1f9c988076d5}

# check if all dependencies are there
for i in go curl git aria2c rsync; do
  if ! command -v "${i}" > /dev/null 2>&1; then
    echo "${i} is not installed, please install it first"
    exit 1
  fi 
done

# check that the correct mmseqs commit is there
if [ -x ./mmseqs/bin/mmseqs ]; then
  # download it again if its a different commit
  if [ $(./mmseqs/bin/mmseqs version) != ${MMSEQS_COMMIT} ]; then 
    rm -rf -- mmseqs
    curl -s -o- https://mmseqs.com/archive/${MMSEQS_COMMIT}/mmseqs-linux-avx2.tar.gz | tar -xzf - mmseqs/bin/mmseqs
  fi
else
  curl -s -o- https://mmseqs.com/archive/${MMSEQS_COMMIT}/mmseqs-linux-avx2.tar.gz | tar -xzf - mmseqs/bin/mmseqs
fi

# mmseqs needs to be in PATH for the setup_databases script to work
PATH="${WORKDIR}/mmseqs/bin:$PATH"


hasCommand () {
    command -v "$1" >/dev/null 2>&1
}

STRATEGY=""
if hasCommand aria2c; then STRATEGY="$STRATEGY ARIA"; fi
if hasCommand curl;   then STRATEGY="$STRATEGY CURL"; fi
if hasCommand wget;   then STRATEGY="$STRATEGY WGET"; fi
if [ "$STRATEGY" = "" ]; then
	    fail "No download tool found in PATH. Please install aria2c, curl or wget."
fi

if [ -n "${PDB_AWS_DOWNLOAD}" ]; then
    if ! hasCommand aws; then
        fail "aws command not found in PATH. Please install the aws command line tool."
    fi
fi

downloadFile() {
    URL="$1"
    OUTPUT="$2"
    set +e
    for i in $STRATEGY; do
        case "$i" in
        ARIA)
            FILENAME=$(basename "${OUTPUT}")
            DIR=$(dirname "${OUTPUT}")
            aria2c --max-connection-per-server="$ARIA_NUM_CONN" --allow-overwrite=true -o "$FILENAME" -d "$DIR" "$URL" && set -e && return 0
            ;;
        CURL)
            curl -L -o "$OUTPUT" "$URL" && set -e && return 0
            ;;
        WGET)
            wget -O "$OUTPUT" "$URL" && set -e && return 0
            ;;
        esac
    done
    set -e
    fail "Could not download $URL to $OUTPUT"
}

# Make MMseqs2 merge the databases to avoid spamming the folder with files
export MMSEQS_FORCE_MERGE=1

curl -s -o- https://mmseqs.com/archive/${1:-4589151554eb83a70ff0c4d04d21b83cabc203e4}/mmseqs-linux-avx2.tar.gz | tar -xzf - -C /scratch/gquarg/ mmseqs/bin/mmseqs

scratch/gquarg/mmseqs/bin/mmseqs

if [ ! -f UNIREF30_READY ]; then
  downloadFile "https://wwwuser.gwdg.de/~compbiol/colabfold/${UNIREF30DB}.tar.gz" "${UNIREF30DB}.tar.gz"
  tar xzvf "${UNIREF30DB}.tar.gz"
  /scratch/gquarg/mmseqs/bin/mmseqs /scratch/gquarg/MMSeqDB/tsv2exprofiledb /scratch/gquarg/MMSeqDB/uniref30_2302 /scratch/gquarg/MMSeqDB/uniref30_2302
  MMSEQS_NO_INDEX=1 MMSEQS_FORCE_MERGE=1 /scratch/gquarg/mmseqs/bin/mmseqs tsv2exprofiledb uniref30_2302 uniref30_2302
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