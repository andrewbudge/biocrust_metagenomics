#!/bin/bash

cd ~/projects/biocrust_metagenomics/

source ~/miniconda3/etc/profile.d/conda.sh
conda activate kraken2

mkdir -p data/genomes/bacterial/taxonomy

# combine all high-quality genomes into one fasta
cat data/genomes/bacterial/checkm2/input/*.fasta \
    > data/genomes/bacterial/taxonomy/hq_chromosomes.fasta

# run kraken2
kraken2 --db ~/databases/k2_pluspf_08gb --threads 6 \
    --output data/genomes/bacterial/taxonomy/hq_chromosomes.out \
    --report data/genomes/bacterial/taxonomy/hq_chromosomes.report \
    data/genomes/bacterial/taxonomy/hq_chromosomes.fasta

conda deactivate
