#!/bin/bash

#curl -o ProteinEnv.yml https://gitlab.uzh.ch/locbp/public/ProteinNotebooks/-/raw/main/InstallationFiles/ProteinEnv.yml

source /apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/mamba-23.11.0-0-334ztq7i4mzu762ew2x3kbbrrorhe6eg/etc/profile.d/conda.sh
module load mamba

notebooks=false
models=false

#Pre-process long options and convert them to short options
for arg in "$@"; do
  shift
  case "$arg" in
    "--notebooks") set -- "$@" "-n" ;;
    "--models")    set -- "$@" "-m" ;;
    *)            set -- "$@" "$arg"
  esac
done

# Process command-line options
while getopts ":nm" opt; do
  case $opt in
    n)
      notebooks=true
      ;;
    m)
      models=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

if $notebooks; then
  echo
  echo "Setting up notebooks..."
  if [ -f "protein_notebooks.sh" ]; then
  rm protein_notebooks.sh
  fi
  curl -o protein_notebooks.sh https://gitlab.uzh.ch/locbp/public/ProteinNotebooks/-/raw/main/InstallationFiles/protein_notebooks.sh
  chmod +x protein_notebooks.sh
  bash protein_notebooks.sh
  rm protein_notebooks.sh
fi


if $models; then
  # Install ProteinEnv if not present
  env_exists=$(conda env list | grep 'ProteinEnv')

  if [ -z "$env_exists" ]; then
  echo "Creating conda environment..."
  cd ProteinNotebooks/InstallationFiles
  mamba env create -f ProteinEnv.yml 
  cd ../..
  else
  echo "ProteinEnv found"
  fi
  echo

  echo
  conda activate ProteinEnv #Common to all models, contains pymol as well

  if [ "$CONDA_DEFAULT_ENV" = "ProteinEnv" ]; then
    echo
    echo "Setting up models..."
    cd /home/$USER/data
    if [ -f "protein_models.sh" ]; then
    rm protein_models.sh
    fi
    wget https://raw.githubusercontent.com/GQChem/ProteinNotebooks/main/InstallationFiles/protein_models.sh
    curl -o protein_models.sh https://gitlab.uzh.ch/locbp/public/ProteinNotebooks/-/raw/main/InstallationFiles/protein_models.sh
    chmod +x protein_models.sh 
    bash protein_models.sh
    rm protein_models.sh
  else
    echo "Could not activate ProteinEnv environment. Check location of conda is /apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/mamba-23.11.0-0-334ztq7i4mzu762ew2x3kbbrrorhe6eg/etc/profile.d/conda.sh by typing: which conda (from cluster shell). Contact administrator otherwise"
  fi
fi

echo
echo "Setup completed"
