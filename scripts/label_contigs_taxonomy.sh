#!/bin/bash

cd ~/projects/biocrust_metagenomics/

# build taxID -> species name lookup from kraken2 report
awk -F'\t' '{gsub(/^ +/, "", $6); print $5 "\t" $6}' \
    data/genome_classification/kraken2/hq_chromosomes.report \
    > data/genome_classification/kraken2/taxid_to_name.tsv

# build contig -> taxID lookup from kraken2 output
awk -F'\t' '{print $2 "\t" $3}' \
    data/genome_classification/kraken2/hq_chromosomes.out \
    > data/genome_classification/kraken2/contig_to_taxid.tsv

# join to create contig -> species name mapping
awk -F'\t' '
    NR==FNR {name[$1] = $2; next}
    {print $1 "\t" ($2 in name ? name[$2] : "unknown") "\t" $2}
' data/genome_classification/kraken2/taxid_to_name.tsv \
  data/genome_classification/kraken2/contig_to_taxid.tsv \
  > data/genome_classification/kraken2/contig_taxonomy.tsv

# rename fasta headers to include species
awk '
    NR==FNR {split($0, a, "\t"); gsub(/ /, "_", a[2]); species[a[1]] = a[2]; next}
    /^>/ {
        split($1, h, ">")
        id = h[2]
        if (id in species) {
            sub(id, id "_" species[id])
        }
        print
        next
    }
    {print}
' data/genome_classification/kraken2/contig_taxonomy.tsv \
  data/genome_classification/kraken2/hq_chromosomes.fasta \
  > data/genome_classification/kraken2/hq_chromosomes_labeled.fasta
