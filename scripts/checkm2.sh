#!/bin/bash

cd ~/projects/biocrust_metagenomics/

source ~/miniconda3/etc/profile.d/conda.sh
conda activate checkm2

# extract chromosome contig IDs from genomad classification
awk -F'\t' 'NR>1 && $2 > 0.9 {print $1}' \
    data/genome_classification/genomad_output/circular_contigs_aggregated_classification/circular_contigs_aggregated_classification.tsv \
    > data/genome_classification/chromosome_ids.txt

# extract each chromosome contig into its own fasta
mkdir -p data/genome_classification/checkm2/input
while read id; do
    seqkit grep -p "$id" data/genome_classification/circular_contigs.fasta \
        > "data/genome_classification/checkm2/input/${id}.fasta"
done < data/genome_classification/chromosome_ids.txt

# run checkm2
mkdir -p data/genome_classification/checkm2/output
checkm2 predict --threads 6 -x .fasta \
    --database_path ~/databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
    -i data/genome_classification/checkm2/input/ \
    -o data/genome_classification/checkm2/output/

conda deactivate
