#!/bin/bash

cd ~/projects/biocrust_metagenomics/

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pharokka

# extract proviruses for Brevundimonas sp. M20
bash scripts/extract_proviruses.sh SRR35809155_ctg370 \
    data/genome_classification/genomes/Brevundimonas_sp_M20/pharokka

# annotate proviruses with pharokka
pharokka.py \
    -i data/genome_classification/genomes/Brevundimonas_sp_M20/pharokka/proviruses.fasta \
    -o data/genome_classification/genomes/Brevundimonas_sp_M20/pharokka/output \
    -d ~/databases/pharokka_db/ \
    -t 6 -m -s \
    --prefix Brevundimonas_proviruses

conda deactivate
