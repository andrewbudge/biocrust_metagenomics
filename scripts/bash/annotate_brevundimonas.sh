#!/bin/bash

cd ~/projects/biocrust_metagenomics/

source ~/miniconda3/etc/profile.d/conda.sh

GENOME_DIR=data/genomes/bacterial/annotations/Brevundimonas_sp_M20
INPUT_FASTA=data/genomes/bacterial/checkm2/input/SRR35809155_ctg370.fasta

conda activate bakta

mkdir -p $GENOME_DIR/bakta

bakta --db ~/databases/bakta_db/db-light/ --complete --threads 6 \
    --genus Brevundimonas --species "sp. M20" \
    --prefix Brevundimonas_sp_M20 \
    --output $GENOME_DIR/bakta \
    --keep-contig-headers \
    $INPUT_FASTA

conda deactivate

conda activate pharokka

bash scripts/extract_proviruses.sh SRR35809155_ctg370 \
    $GENOME_DIR/pharokka

pharokka.py \
    -i $GENOME_DIR/pharokka/proviruses.fasta \
    -o $GENOME_DIR/pharokka/output \
    -d ~/databases/pharokka_db/ \
    -t 6 -m -s \
    --prefix Brevundimonas_proviruses

conda deactivate
