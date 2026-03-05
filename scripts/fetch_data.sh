#!/bin/bash
cd ~/projects/biocrust_metagenomics

# make raw data dir
mkdir -p data/reads/raw

# prefetch and download data using sratoolkit
for acc in $(cut -d',' -f10 data/metadata/biocrust_metadata.csv | tail -n +2); do
    prefetch -p -v "$acc" && \
    fastq-dump -v "$acc" -O data/reads/raw && \
    pigz -v data/reads/raw/"$acc".fastq && \
    rm -rf "$acc"
done
