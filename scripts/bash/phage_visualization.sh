#!/bin/bash

# Visualize HQ phage genomes with LoVis4u
# Grouped by shared structure, plus unique singletons

cd ~/projects/biocrust_metagenomics/

source ~/miniconda3/etc/profile.d/conda.sh
conda activate lovis4u

GFF_DIR=data/genomes/viral/pharokka/annotations/single_gffs
OUT_DIR=data/genomes/viral/lovis4u

mkdir -p $OUT_DIR/group1 $OUT_DIR/group2 $OUT_DIR/group3 $OUT_DIR/unique

# --- Group 1: 39,631 bp phage (4 samples) ---
mkdir -p /tmp/lovis4u_group1
cp "$GFF_DIR/SRR35809155_ctg108|provirus_3612395_3652025.gff" /tmp/lovis4u_group1/
cp "$GFF_DIR/SRR35808611_ctg442|provirus_3721232_3760862.gff" /tmp/lovis4u_group1/
cp "$GFF_DIR/SRR35809332_ctg1415|provirus_2633736_2673366.gff" /tmp/lovis4u_group1/
cp "$GFF_DIR/SRR35808553_ctg921|provirus_1054300_1093930.gff" /tmp/lovis4u_group1/

lovis4u -gff /tmp/lovis4u_group1/ \
    -o $OUT_DIR/group1 \
    -hl --set-category-colour
# --- Group 2: ~45,334 bp phage (4 samples) ---
mkdir -p /tmp/lovis4u_group2
cp "$GFF_DIR/SRR35809155_ctg108|provirus_3243057_3288390.gff" /tmp/lovis4u_group2/
cp "$GFF_DIR/SRR35809332_ctg1415|provirus_2264398_2309731.gff" /tmp/lovis4u_group2/
cp "$GFF_DIR/SRR35808553_ctg921|provirus_684962_730295.gff" /tmp/lovis4u_group2/
cp "$GFF_DIR/SRR35808611_ctg442|provirus_4085068_4130200.gff" /tmp/lovis4u_group2/

lovis4u -gff /tmp/lovis4u_group2/ \
    -o $OUT_DIR/group2 \
    -hl --set-category-colour
# --- Group 3: 42,249 bp phage (2 samples) ---
mkdir -p /tmp/lovis4u_group3
cp "$GFF_DIR/SRR35809332_ctg866|provirus_3351243_3393491.gff" /tmp/lovis4u_group3/
cp "$GFF_DIR/SRR35808611_ctg1214|provirus_771119_813367.gff" /tmp/lovis4u_group3/

lovis4u -gff /tmp/lovis4u_group3/ \
    -o $OUT_DIR/group3 \
    -hl --set-category-colour
# --- Unique singletons (all in one image) ---
mkdir -p /tmp/lovis4u_unique
cp "$GFF_DIR/SRR35808900_ctg3305.gff" /tmp/lovis4u_unique/
cp "$GFF_DIR/SRR35808611_ctg78|provirus_2163821_2228021.gff" /tmp/lovis4u_unique/
cp "$GFF_DIR/SRR35808611_ctg514|provirus_397177_435424.gff" /tmp/lovis4u_unique/
cp "$GFF_DIR/SRR35808611_ctg707|provirus_5227802_5272024.gff" /tmp/lovis4u_unique/
cp "$GFF_DIR/SRR35809336_ctg10052|provirus_1805875_1849128.gff" /tmp/lovis4u_unique/
cp "$GFF_DIR/SRR35809336_ctg10687|provirus_4211115_4255250.gff" /tmp/lovis4u_unique/

lovis4u -gff /tmp/lovis4u_unique/ \
    -o $OUT_DIR/unique \
    -hl --set-category-colour
# cleanup
rm -rf /tmp/lovis4u_group1 /tmp/lovis4u_group2 /tmp/lovis4u_group3 /tmp/lovis4u_unique

conda deactivate
