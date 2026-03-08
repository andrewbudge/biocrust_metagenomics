#!/bin/bash

cd ~/projects/biocrust_metagenomics/

source ~/miniconda3/etc/profile.d/conda.sh
conda activate eggnog

mkdir -p data/gene_annotation/eggNOG-mapper_output/

for f in data/gene_annotation/prodigal_output/*/*_proteins.faa; do
    base=$(basename "$f" _proteins.faa)
    mkdir -p "data/gene_annotation/eggNOG-mapper_output/$base"
    emapper.py -i "$f" \
        --output "$base" \
        --output_dir "data/gene_annotation/eggNOG-mapper_output/$base" \
        --data_dir ~/databases/eggnog_db \
        -m diamond \
        --cpu 8
done

conda deactivate
