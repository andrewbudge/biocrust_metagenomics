#setwd("projects/biocrust_metagenomics/")

library(tidyverse)
library(viridis)

# Read in metadata for sample info
meta_data <- read_csv("data/metadata/biocrust_metadata.csv")

# Read in all eggNOG annotation files
eggnog_data <- list.files("data/gene_annotation/eggNOG-mapper_output",
                          pattern = "*.emapper.annotations",
                          recursive = TRUE, full.names = TRUE) %>%
  map_dfr(function(f) {
    srr <- basename(f) %>% gsub(".emapper.annotations", "", .)
    read_tsv(f, comment = "##") %>%
      mutate(srr = srr)
  })

# join with metadata for site info
eggnog_data <- eggnog_data %>%
  left_join(meta_data, by = c("srr" = "srr"))

# filter data
COG_data <- eggnog_data %>%
  mutate(sample = str_extract(.data$title, "[^-]+$") %>% trimws(),
         geo_loc_name = gsub("USA: ", "", geo_loc_name)) %>%
  select(`#query`, COG_category, evalue, srr, sample, geo_loc_name) %>%
  filter(evalue < 1e-5) %>%
  filter(COG_category != "-") %>%
  separate_longer_position(COG_category, width = 1) %>% # split cog ID's into single letters
  filter(!COG_category %in% c("S", "R", "Y")) # remove the idk categories


# translation table from single letter codes to description 
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

# append description
COG_data <- COG_data %>%
  mutate(COG_description = cog_names[COG_category])

# normalize to proportions within each sample
COG_proportions <- COG_data %>%
  group_by(sample) %>%
  mutate(total = n()) %>%
  ungroup() %>%
  group_by(sample, COG_description, geo_loc_name) %>%
  summarize(proportion = n() / first(total), .groups = "drop")

# plot COG data by sample, faceted by state
ggplot(COG_proportions, aes(x = reorder(COG_description, proportion), y = proportion, fill = sample)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~geo_loc_name, scales = "free_x") +
  scale_fill_viridis_d(option = "A", begin = 0.1, end = 0.85) +
  labs(
    title = "Functional Profile of Desert Biocrust Metagenomes",
    x = "Predicted Gene Function",
    y = "Proportion",
    fill = "Sample") +
  theme(
    plot.title = element_text(size = 16, hjust = .5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title.position = "plot"
  )
