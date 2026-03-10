#!/bin/bash

cd ~/projects/biocrust_metagenomics/

mkdir -p data/genome_classification/genomad_output

genomad end-to-end --cleanup --threads 6 --verbose --splits 6 \
    data/genome_classification/circular_contigs.fasta \
    data/genome_classification/genomad_output/ \
    ~/databases/genomad_db/
