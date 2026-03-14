#!/bin/bash

cd ~/projects/biocrust_metagenomics/

source ~/miniconda3/etc/profile.d/conda.sh
conda activate pharokka

CHECKV_DIR=data/genomes/viral/checkv
VIRAL_DIR=data/genomes/viral
PHAROKKA_DIR=$VIRAL_DIR/pharokka

mkdir -p $PHAROKKA_DIR

# grab best hits, and the one non-provirus sample
awk -F'\t' 'NR==1 || ($10 == 100.0 && $12 == 0.0 && $14 == "") || /ctg3305/' \
    $CHECKV_DIR/quality_summary.tsv > $PHAROKKA_DIR/top_hits_for_pharroka.tsv

# get ids for seqkit grep
cut -f1 $PHAROKKA_DIR/top_hits_for_pharroka.tsv | grep -v "contig" > $PHAROKKA_DIR/pharroka_virus_ids.txt

# extract best hits
seqkit grep -f $PHAROKKA_DIR/pharroka_virus_ids.txt $VIRAL_DIR/circular_contigs_virus.fna \
    > $PHAROKKA_DIR/hq_viral_genomes.fna

# annotate with pharokka
pharokka.py -i $PHAROKKA_DIR/hq_viral_genomes.fna \
    -o $PHAROKKA_DIR/annotations \
    -d ~/databases/pharokka_db/ \
    -t 6 -m -s

conda deactivate
