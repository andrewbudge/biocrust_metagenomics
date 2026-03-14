#!/bin/bash

cd ~/projects/biocrust_metagenomics/

# make dir for assemblies
mkdir -p data/assemblies

for file in data/reads/QC/*.fastq.gz; do
    base=$(basename "$file" _QC.fastq.gz)
    echo "Assembling: $base"
    mkdir -p "data/assemblies/$base"
    metaMDBG asm --out-dir "data/assemblies/$base" --in-hifi "$file" --threads 6
    echo "Done: $base"
done

echo "All assemblies complete"
