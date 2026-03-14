#!/bin/bash

cd ~/projects/biocrust_metagenomics/

# make dir for for prodigal annotations
mkdir -p data/assembly_COG

# clean up the tmp dirs in the assemble dir
rm -rf data/assemblies/*/tmp

# rename contigs and label headers with sample ID
for f in data/assemblies/*/contigs.fasta.gz; do
    base=$(basename $(dirname "$f"))
    pigz -dc "$f" | sed "s/^>/>${base}_/" | pigz > "data/assemblies/$base/${base}_contigs.fasta.gz.tmp"
    mv "data/assemblies/$base/${base}_contigs.fasta.gz.tmp" "data/assemblies/$base/${base}_contigs.fasta.gz"
done

# run prodigal in metagenomic mode
for f in data/assemblies/*/*_contigs.fasta.gz; do
    base=$(basename "$f" _contigs.fasta.gz)
    mkdir -p "data/assembly_COG/$base"
    pigz -dc "$f" | prodigal -p meta -i /dev/stdin \
        -a "data/assembly_COG/$base/${base}_proteins.faa" \
        -o "data/assembly_COG/$base/${base}_genes.gff" \
        -f gff
done

