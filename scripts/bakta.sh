#!/bin/bash

cd ~/projects/biocrust_metagenomics/

source ~/miniconda3/etc/profile.d/conda.sh
conda activate bakta

mkdir -p data/genome_classification/genomes/Brevundimonas_sp_M20/bakta

# annotate Brevundimonas sp. M20 (4 proviruses, found in 5/6 samples)
bakta --db ~/databases/bakta_db/db-light/ --complete --threads 6 \
    --genus Brevundimonas --species "sp. M20" \
    --prefix Brevundimonas_sp_M20 \
    --output data/genome_classification/genomes/Brevundimonas_sp_M20/bakta \
    --keep-contig-headers \
    data/genome_classification/checkm2/input/SRR35809155_ctg370.fasta

conda deactivate
