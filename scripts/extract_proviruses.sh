#!/bin/bash

cd ~/projects/biocrust_metagenomics/

# extract proviruses for a specific contig from genomad output
# usage: bash scripts/extract_proviruses.sh <contig_id> <output_dir>
# example: bash scripts/extract_proviruses.sh SRR35809155_ctg370 data/genome_classification/genomes/Brevundimonas_sp_M20/pharokka

CONTIG_ID=$1
OUTPUT_DIR=$2

mkdir -p "$OUTPUT_DIR"

# extract provirus sequences matching the contig ID
seqkit grep -r -p "$CONTIG_ID" \
    data/genome_classification/genomad_output/circular_contigs_find_proviruses/circular_contigs_provirus.fna \
    > "$OUTPUT_DIR/proviruses.fasta"

echo "Extracted $(grep -c '^>' "$OUTPUT_DIR/proviruses.fasta") proviruses for $CONTIG_ID"
seqkit stats "$OUTPUT_DIR/proviruses.fasta"
