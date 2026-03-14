#!/bin/bash

cd ~/projects/biocrust_metagenomics/

mkdir -p data/genomes/genomad_output

genomad end-to-end --cleanup --threads 6 --verbose --splits 6 \
    data/genomes/circular_contigs.fasta \
    data/genomes/genomad_output/ \
    ~/databases/genomad_db/
