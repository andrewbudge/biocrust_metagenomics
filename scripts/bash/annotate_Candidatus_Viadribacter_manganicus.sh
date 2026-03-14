#!/bin/bash

cd ~/projects/biocrust_metagenomics/

# Conda Pain
source ~/miniconda3/etc/profile.d/conda.sh

mkdir -p data/genomes/bacterial/annotations/Candidatus_Viadribacter_manganicus

# store path in a var since I am not writing all that
CVm_annotation_dir=data/genomes/bacterial/annotations/Candidatus_Viadribacter_manganicus

# Duct tape and glue never looked so good
cat data/genomes/bacterial/taxonomy/hq_chromosomes_labeled.fasta | seqkit seq -w0 | rg "SRR35809336_ctg8014" -A1 | sed 's/ /_/g ' > "$CVm_annotation_dir/Candidatus_Viadribacter_manganicus.fna"

# Run bakta on the genome
conda activate bakta

bakta --db ~/databases/bakta_db/db-light/ --complete --threads 8 \
    --genus Viadribacter --species manganicus \
    --prefix Candidatus_Viadribacter_manganicus \
    --output "$CVm_annotation_dir/bakta_results" \
    --keep-contig-headers \
    "$CVm_annotation_dir/Candidatus_Viadribacter_manganicus.fna"

conda deactivate
