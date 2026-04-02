library(maps)
library(vegan)
library(ggrepel)
setwd("~/projects/biocrust_metagenomics/")
source("scripts/R/shared.R")

# ---- Create master kraken2 dataset ----

master_kraken2_data <- list.files("data/reads/taxonomy/reports", full.names = TRUE) %>%
  map_dfr(function(f) {
    srr_num <- basename(gsub(".report", "", f))
    read_tsv(f, col_names = c("pct", "reads_clade", "reads_direct", "rank", "taxid", "name")) %>%
      mutate(srr = srr_num)
  }) %>%
  left_join(meta_data, by = "srr") %>%
  mutate(sample = str_extract(title, "(?<=- ).*"))

# ---- Biocrust community composition at phylum level ----

phylum_data <- master_kraken2_data %>%
  filter(rank == "P", pct > 1) %>%
  select(sample, name, pct, geo_loc_name)

ggplot(phylum_data, aes(x = sample, y = pct, fill = fct_reorder(name, pct))) +
  geom_col() +
  facet_wrap(~ geo_loc_name, scales = "free_x") +
  scale_fill_viridis_d(option = "A", begin = 0.1, end = 0.85) +
  labs(
    title = "Desert Biocrust Community Composition",
    x = "Sample ID",
    y = "Relative Abundance (%)",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(title = NULL, nrow = 2))

# ---- Alpha and Beta diversity ----

species_count <- master_kraken2_data %>%
  filter(rank == "S", reads_direct > 10) %>%
  select(geo_loc_name, sample, name, reads_direct)

# Species count per sample (horizontal bar plot)
species_per_sample <- species_count %>%
  group_by(sample, geo_loc_name) %>%
  summarise(n_species = n_distinct(name), .groups = "drop") %>%
  mutate(geo_loc_name = str_replace(geo_loc_name, "USA: ", ""))

ggplot(species_per_sample, aes(x = n_species, y = reorder(sample, n_species), fill = geo_loc_name)) +
  geom_col() +
  scale_fill_viridis_d(option = "A", begin = 0.1, end = 0.85) +
  labs(
    title = "Species Detected per Sample",
    x = "Number of Species",
    y = NULL,
    fill = "Location"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

vegan_species_data <- species_count %>%
  pivot_wider(names_from = name, values_from = reads_direct)

sample_info <- vegan_species_data %>%
  select(sample, geo_loc_name) %>%
  mutate(geo_loc_name = str_replace(geo_loc_name, "USA: ", ""))

vegan_matrix <- vegan_species_data %>%
  select(-sample, -geo_loc_name) %>%
  replace(is.na(.), 0)

sample_info$shannon <- diversity(vegan_matrix, index = "shannon")

# Alpha diversity plot
ggplot(sample_info, aes(x = sample, y = shannon, fill = geo_loc_name)) +
  geom_col() +
  facet_wrap(~ geo_loc_name, scales = "free_x") +
  scale_fill_viridis_d(option = "A", begin = 0.1, end = 0.85) +
  labs(
    title = "Desert Biocrust Community Alpha Diversity",
    x = "Sample ID",
    y = "Shannon Diversity",
    fill = "Location"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

# Beta diversity (PCoA)
pcoa_df <- vegdist(vegan_matrix, method = "bray") %>%
  cmdscale(k = 2) %>%
  data.frame(sample = sample_info$sample, geo_loc_name = sample_info$geo_loc_name)

ggplot(pcoa_df, aes(x = X1, y = X2, color = geo_loc_name)) +
  geom_point(size = 4) +
  geom_text(aes(label = sample), nudge_y = 0.03) +
  labs(
    title = "PCoA of Biocrust Communities",
    x = "PCoA Axis 1",
    y = "PCoA Axis 2",
    color = NULL
  ) +
  scale_color_viridis_d(option = "A", begin = 0.1, end = 0.85) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(color = guide_legend(title = NULL, nrow = 1))

# ---- Sample map ----

map_info <- meta_data %>%
  separate(lat_lon, into = c("lat", "lon"), sep = " N ") %>%
  mutate(
    lat = as.numeric(lat),
    lon = -parse_number(lon),
    sample = str_extract(title, "(?<=- ).*")
  ) %>%
  select(sample, lat, lon)

us_map <- map_data("state") %>%
  filter(region %in% c("utah", "arizona", "new mexico"))

ggplot() +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
               fill = "white", color = "black") +
  geom_point(data = map_info, aes(x = lon, y = lat), size = 4) +
  geom_text_repel(data = map_info, aes(x = lon, y = lat, label = sample),
                  size = 5, box.padding = .5) +
  theme_minimal()

# ---- Kraken2 of contigs: Circular genome taxonomy ----

contig_taxonomy <- read_tsv(
  "data/genomes/bacterial/taxonomy/contig_taxonomy.tsv",
  col_names = c("contig", "species", "taxid")
)

checkm2 <- read_tsv("data/genomes/bacterial/checkm2/output/quality_report.tsv")

genome_summary <- contig_taxonomy %>%
  left_join(checkm2, by = c("contig" = "Name")) %>%
  filter(Completeness > 90, Contamination < 5) %>%
  mutate(sample = str_extract(contig, "^SRR\\d+")) %>%
  left_join(meta_data, by = c("sample" = "srr")) %>%
  mutate(sample_name = str_extract(title, "(?<=- ).*"))

# Extract coverage from fasta headers
fasta_headers <- read_lines(pipe("grep '^>' data/genomes/bacterial/taxonomy/hq_chromosomes.fasta")) %>%
  keep(~ str_starts(.x, ">")) %>%
  tibble(header = .) %>%
  mutate(
    contig = str_extract(header, "(?<=>)\\S+"),
    coverage = as.numeric(str_extract(header, "(?<=coverage=)[\\d.]+"))
  ) %>%
  select(contig, coverage)

genome_summary <- genome_summary %>%
  left_join(fasta_headers, by = "contig")

genomes_within_sample <- genome_summary %>%
  select(sample, sample_name, species, geo_loc_name, coverage) %>%
  mutate(
    geo_loc_name = str_replace(geo_loc_name, "USA: ", ""),
    genus = str_extract(species, "^\\S+"),
    phylum = phylum_map[genus]
  )

species_order <- genomes_within_sample %>%
  distinct(species, phylum) %>%
  arrange(phylum, species) %>%
  pull(species)

genomes_within_sample <- genomes_within_sample %>%
  mutate(species = factor(species, levels = rev(species_order)))

# Chromosome count per species (horizontal bar plot)
chromosome_counts <- genomes_within_sample %>%
  mutate(species = as.character(species)) %>%
  count(species, name = "n_chromosomes")

ggplot(chromosome_counts, aes(x = n_chromosomes, y = reorder(species, n_chromosomes))) +
  geom_col(fill = viridis::viridis(1, option = "A", begin = 0.7)) +
  labs(
    title = "Recovered Bacterial Chromosomes",
    x = "Number of Chromosomes",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggplot(genomes_within_sample, aes(x = sample_name, y = species, color = phylum, size = coverage)) +
  geom_point() +
  scale_color_viridis_d(option = "A", begin = 0.1, end = 0.85) +
  facet_wrap(~ geo_loc_name, scales = "free_x") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, color = "grey20"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = NA),
    panel.spacing = unit(1, "lines"),
    panel.grid.major.y = element_line(linewidth = 0.3),
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12)
  ) +
  labs(
    title = "Complete Genomes Recovered by Sample",
    x = "Sample",
    y = NULL,
    color = "Phylum"
  )
