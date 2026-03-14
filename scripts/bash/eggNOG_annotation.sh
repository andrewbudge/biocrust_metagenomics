#!/bin/bash

cd ~/projects/biocrust_metagenomics/

source ~/miniconda3/etc/profile.d/conda.sh
conda activate eggnog

mkdir -p data/assembly_COG/eggNOG-mapper_output/

for f in data/assembly_COG/prodigal_output/*/*_proteins.faa; do
    base=$(basename "$f" _proteins.faa)
    mkdir -p "data/assembly_COG/eggNOG-mapper_output/$base"
    emapper.py -i "$f" \
        --output "$base" \
        --output_dir "data/assembly_COG/eggNOG-mapper_output/$base" \
        --data_dir ~/databases/eggnog_db \
        -m diamond \
        --cpu 8
done

conda deactivate
