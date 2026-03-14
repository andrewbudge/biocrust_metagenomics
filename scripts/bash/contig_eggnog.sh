#!/bin/bash

# Predict genes on HQ bacterial chromosomes and annotate with eggNOG-mapper
# Input: labeled fasta with species names in headers
# Output: prodigal gene predictions + eggNOG-mapper COG annotations

cd ~/projects/biocrust_metagenomics/

BASE_DIR=data/genomes/bacterial
INPUT_FASTA=$BASE_DIR/taxonomy/hq_chromosomes_labeled.fasta
PRODIGAL_DIR=$BASE_DIR/contig_COG/prodigal_output
EGGNOG_DIR=$BASE_DIR/contig_COG/eggNOG_map_output

# Step 1: Predict genes with Prodigal (metagenome mode)
prodigal -i $INPUT_FASTA \
  -a $PRODIGAL_DIR/hq_chromosomes_proteins.faa \
  -o $PRODIGAL_DIR/hq_chromosomes_genes.gff \
  -f gff \
  -p meta

# Step 2: Run eggNOG-mapper on predicted proteins
source ~/miniconda3/etc/profile.d/conda.sh
conda activate eggnog

emapper.py -i $PRODIGAL_DIR/hq_chromosomes_proteins.faa \
  --output contig_eggnog \
  --output_dir $EGGNOG_DIR \
  --data_dir ~/databases/eggnog_db \
  -m diamond \
  --cpu 8

conda deactivate
