#!/bin/bash
cd ~/projects/biocrust_metagenomics

# conda, because we can't have nice things
source ~/miniconda3/etc/profile.d/conda.sh
conda activate kraken2

# Make output dirs
mkdir -p data/kraken2/reports
mkdir -p data/kraken2/outputs

# Run Kraken2 on the QC reads
echo "Running Kraken2..."
for file in data/reads/QC/*.fastq.gz; do
    base=$(basename "$file" _QC.fastq.gz)
    echo "Classifying: $base"
    kraken2 --db ~/databases/k2_pluspf_08gb \
        --threads 8 \
        --output data/kraken2/outputs/"${base}".out \
        --report data/kraken2/reports/"${base}".report \
        "$file"
    echo "Done: $base"
done

echo "All samples classified!"
conda deactivate
