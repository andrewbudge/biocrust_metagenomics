#!/bin/bash

cd ~/projects/biocrust_metagenomics

# use conda, the enemy of all
# If you want to use this, replace with whatever env has chopper on it
source ~/miniconda3/etc/profile.d/conda.sh
conda activate base

# Make QC dir
mkdir -p data/reads/QC

# run QC with chopper
echo "Starting QC with chopper..."
for file in data/reads/raw/*.fastq.gz; do
    base=$(basename "$file" .fastq.gz)
    echo "Processing: $base"
    pigz -p 6 -dc "$file" | chopper -t 8 -q 20 -l 1000 | pigz -p 6 > data/reads/QC/"${base}"_QC.fastq.gz
    echo "Done: $base"
done

echo "All samples complete!"
conda deactivate
