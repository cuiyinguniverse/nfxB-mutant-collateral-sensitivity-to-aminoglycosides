#!/bin/bash
. /home/houj21/program/af2/localcolabfold/conda/etc/profile.d/conda.sh
conda activate /home/houj21/program/af2/localcolabfold/colabfold-conda
CURRENTPATH=`pwd`
colabfold_search "${CURRENTPATH}/input_sequences.fasta" "/data_1/houj21/af2" "${CURRENTPATH}/msas" --threads 40
