#!/bin/bash

cd ~/projects/biocrust_metagenomics/

source ~/miniconda3/etc/profile.d/conda.sh
conda activate kraken2

mkdir -p data/genome_classification/kraken2

# combine all high-quality genomes into one fasta
cat data/genome_classification/checkm2/input/*.fasta \
    > data/genome_classification/kraken2/hq_chromosomes.fasta

# run kraken2
kraken2 --db ~/databases/k2_pluspf_08gb --threads 6 \
    --output data/genome_classification/kraken2/hq_chromosomes.out \
    --report data/genome_classification/kraken2/hq_chromosomes.report \
    data/genome_classification/kraken2/hq_chromosomes.fasta

conda deactivate
