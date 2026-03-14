# shared.R — common lookup tables, data loading, and helper functions
# Source this file at the top of any R analysis script

library(tidyverse)

# Read in metadata for sample info
meta_data <- read_csv("data/metadata/biocrust_metadata.csv")

# COG single letter code to description lookup
cog_names <- c(
  J = "Translation", A = "RNA processing", K = "Transcription",
  L = "Replication/repair", B = "Chromatin",
  D = "Cell cycle", V = "Defense", T = "Signal transduction",
  M = "Cell wall/membrane", N = "Cell motility", U = "Secretion/trafficking",
  O = "Protein turnover", W = "Extracellular", Z = "Cytoskeleton",
  C = "Energy production", G = "Carbohydrate metabolism",
  E = "Amino acid metabolism", F = "Nucleotide metabolism",
  H = "Coenzyme metabolism", I = "Lipid metabolism",
  P = "Inorganic ion transport", Q = "Secondary metabolites"
)

# genus to phylum lookup
phylum_map <- c(
  "Acidovorax" = "Pseudomonadota",
  "Allocoleopsis" = "Cyanobacteriota",
  "Altererythrobacter" = "Pseudomonadota",
  "Anatilimnocola" = "Planctomycetota",
  "Bosea" = "Pseudomonadota",
  "Bradyrhizobium" = "Pseudomonadota",
  "Brevundimonas" = "Pseudomonadota",
  "Candidatus" = "Pseudomonadota",
  "Chitinophaga" = "Bacteroidota",
  "Chryseolinea" = "Bacteroidota",
  "Devosia" = "Pseudomonadota",
  "Ferrovibrio" = "Pseudomonadota",
  "Ferruginibacter" = "Bacteroidota",
  "Flavisolibacter" = "Bacteroidota",
  "Flavobacterium" = "Bacteroidota",
  "Humisphaera" = "Planctomycetota",
  "Hydrogenophaga" = "Pseudomonadota",
  "Hyphomicrobiales" = "Pseudomonadota",
  "Lysobacter" = "Pseudomonadota",
  "Mesorhizobium" = "Pseudomonadota",
  "Methyloversatilis" = "Pseudomonadota",
  "Miltoncostaea" = "Actinomycetota",
  "Nostoc" = "Cyanobacteriota",
  "Phreatobacter" = "Pseudomonadota",
  "Pseudoxanthomonas" = "Pseudomonadota",
  "Rhizobium" = "Pseudomonadota",
  "Rhodocytophaga" = "Bacteroidota",
  "Roseateles" = "Pseudomonadota",
  "Roseomonas" = "Pseudomonadota",
  "Sediminibacterium" = "Bacteroidota",
  "Sphingomonas" = "Pseudomonadota",
  "Sphingopyxis" = "Pseudomonadota",
  "Tellurirhabdus" = "Bacteroidota",
  "Terricaulis" = "Pseudomonadota",
  "Terrirubrum" = "Pseudomonadota",
  "unclassified" = "Pseudomonadota"
)

# shared COG filtering: remove low confidence, unknowns, and ambiguous categories
filter_cog <- function(df) {
  df %>%
    filter(evalue < 1e-5, COG_category != "-") %>%
    separate_longer_position(COG_category, width = 1) %>%
    filter(!COG_category %in% c("S", "R", "Y")) %>%
    mutate(COG_description = cog_names[COG_category])
}

# clean geo_loc_name by removing "USA: " prefix
clean_geo_loc <- function(x) gsub("USA: ", "", x)
