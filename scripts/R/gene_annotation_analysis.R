setwd("~/projects/biocrust_metagenomics/")
source("scripts/R/shared.R")

# assembly Annotation ----

eggnog_data <- list.files("data/assembly_COG/eggNOG-mapper_output",
                          pattern = "*.emapper.annotations",
                          recursive = TRUE, full.names = TRUE) %>%
  map_dfr(function(f) {
    srr <- basename(f) %>% gsub(".emapper.annotations", "", .)
    read_tsv(f, comment = "##") %>%
      mutate(srr = srr)
  }) %>%
  left_join(meta_data, by = "srr")

COG_data <- eggnog_data %>%
  mutate(sample = str_extract(title, "[^-]+$") %>% trimws(),
         geo_loc_name = gsub("USA: ", "", geo_loc_name)) %>%
  select(`#query`, COG_category, evalue, srr, sample, geo_loc_name) %>%
  filter_cog()

COG_proportions <- COG_data %>%
  count(sample, COG_description, geo_loc_name) %>%
  group_by(sample) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup() %>%
  mutate(sample = factor(sample, levels = c("CYAN3", "SON57", "HS003", "SEV30", "SON60", "HSN023")))

# assembly based COG plot ----

ggplot(COG_proportions, aes(x = reorder(COG_description, proportion), y = proportion, fill = sample)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~geo_loc_name, scales = "free_x") +
  scale_fill_viridis_d(option = "A", begin = 0.1, end = 0.85,
                       breaks = c("SON57", "SON60", "CYAN3", "SEV30", "HS003", "HSN023")) +
  labs(title = "Functional Profile of Desert Biocrust Metagenomes",
       x = "Predicted Gene Function", y = "Proportion", fill = "Sample") +
  theme(
    plot.title = element_text(size = 16, hjust = .5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title.position = "plot"
  )

# Chromosome based eggNOG ----

chromosome_eggnog_data <- read_tsv(
  "data/genomes/bacterial/contig_COG/eggNOG_map_output/contig_eggnog.emapper.annotations",
  comment = "##"
) %>%
  rename(query = `#query`) %>%
  mutate(
    srr = str_extract(query, "^SRR\\d+"),
    contig = str_extract(query, "^SRR\\d+_ctg\\d+"),
    species = str_replace(query, "^SRR\\d+_ctg\\d+_(.+)_\\d+$", "\\1"),
    gene_num = str_extract(query, "\\d+$")
  ) %>%
  left_join(meta_data, by = "srr") %>%
  mutate(geo_loc_name = gsub("USA: ", "", geo_loc_name)) %>%
  filter_cog() %>%
  mutate(
    genus = str_extract(species, "^[A-Za-z]+"),
    phylum = phylum_map[genus]
  )

# filter to species with enough annotated genes for meaningful proportions
tile_data <- chromosome_eggnog_data %>%
  group_by(species) %>%
  filter(n() >= 50) %>%
  ungroup() %>%
  count(species, phylum, COG_description) %>%
  group_by(species) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  complete(nesting(species, phylum), COG_description, fill = list(prop = 0, n = 0)) %>%
  mutate(prop = if_else(prop == 0, NA_real_, prop)) %>%
  mutate(species = str_replace_all(species, "_", " "))

# Chromosome based COG heatmap ----

ggplot(tile_data, aes(x = COG_description, y = species, fill = prop)) +
  facet_grid(phylum ~ ., scales = "free_y", space = "free_y") +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "A", begin = 0.1, end = 0.85, na.value = "grey60") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0),
    legend.position = "right",
  ) +
  labs(
    title = "Functional Profiles of Desert Biocrust Bacteria", 
    x = "Predicted Gene Function",
    y = "Species",
    fill = "Proportion"
  )
