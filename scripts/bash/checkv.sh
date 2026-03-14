#!/bin/bash

conda activate checkv

# cp viral hits from genomad over to viral dir
cp data/genomes/genomad_output/circular_contigs_summary/circular_contigs_virus.fna data/genomes/viral

# run whole checkv pipline
checkv end_to_end data/genomes/viral/circular_contigs_virus.fna data/genomes/viral/checkv/ -d ~/databases/checkv-db-v1.5/ --remove_tmp -t 6
