#!/bin/bash

cd ~/projects/biocrust_metagenomics/

mkdir -p data/genome_classification

# extract circular contig IDs filtered by length >= 5kb and coverage >= 3x
pigz -dc data/assemblies/*/*_contigs.fasta.gz \
    | rg "circular=yes" \
    | awk '{
        for(i=1; i<=NF; i++) {
            if($i ~ /^length=/) len = substr($i, 8) + 0
            if($i ~ /^coverage=/) cov = substr($i, 10) + 0
        }
        if(len >= 5000 && cov >= 3) print $1
    }' \
    | sed 's/>//' > data/metadata/circular_ids.txt

# extract the sequences using the ID list
seqkit grep -f data/metadata/circular_ids.txt \
    data/assemblies/*/*_contigs.fasta.gz \
    > data/genome_classification/circular_contigs.fasta
