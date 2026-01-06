#!/bin/bash
set -e

ENV_FILE="environment.yml"

# carregar conda
source $(conda info --base)/etc/profile.d/conda.sh

# criar ambiente
conda env create -f $ENV_FILE || echo "Ambiente já existe"

conda activate ngs_env

echo "✅ Ambiente pronto"
