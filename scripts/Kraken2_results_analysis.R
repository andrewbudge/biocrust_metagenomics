# load packages
library(tidyverse)
library(viridis)
library(maps)
library(vegan)
library(ggrepel)
library(patchwork)
#setwd("~/projects/biocrust_metagenomics/")

# ---- create master kraken2 dataset ----

# create master data frame
master_kraken2_data <- tibble()

# read in metadata
meta_data <- read.csv("data/metadata/biocrust_metadata.csv")

# Fill master Kraken df
master_kraken2_data <- list.files("data/kraken2/reports", full.names = TRUE) %>%                                                     
  map_dfr(function(f) {
    #get srr num    
    srr_num <- basename(gsub(".report", "", f))
    # read in and data and add srr num to it
    srr_kraken_data <- read_tsv(f, col_names = c("pct", "reads_clade", "reads_direct", "rank", "taxid", "name"))
    srr_kraken_data <- mutate(srr_kraken_data, srr = srr_num)
  })

# attach metadata
master_kraken2_data <- master_kraken2_data %>%
  left_join(meta_data, by = "srr") %>%
  mutate(sample = str_extract(.data$title, "(?<=- ).*"))                                                                                   

# ---- Biocrust community composition at phylum level ----

# filter data, phylum level, at least one percent abundance
phylum_data <- master_kraken2_data %>%
  filter(rank == "P") %>%
  select(sample, name, pct, geo_loc_name) %>%
  filter(pct > 1)

# Plot (very pretty)
ggplot(phylum_data, aes(x = sample, y = pct, fill = fct_reorder(name, pct))) +
  geom_col() +                                                                                                 
  facet_wrap(~ geo_loc_name, scales = "free_x") +
  scale_fill_viridis_d(
    option = "A",
    begin = 0.1, end = 0.85) +
  labs(
    title = "Desert Biocrust Community Composition",
    x = "Sample ID",
    y = "Relative Abundance (%)",
    fill = "Phylum"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title = NULL, nrow = 2))

# Save this john
#ggsave("../plots/phylum_community_barplot.png", width = 10, height = 6, dpi = 300)
#ggsave("../plots/phylum_community_barplot.pdf", width = 10, height = 6)


# ---- Alpha and Beta diversity calculation ----

# get species counts
species_count <-  master_kraken2_data %>%
  filter(rank == "S") %>%
  filter(reads_direct > 10) %>%
  select(geo_loc_name, sample, name, reads_direct)

# Clean up the data and get it in the format vegan likes
vegan_species_data <- species_count %>%
  pivot_wider(names_from = name, values_from = reads_direct)

sample_info <- vegan_species_data %>%
  select(sample, geo_loc_name) %>%
  mutate(geo_loc_name = str_replace(geo_loc_name, "USA: ", ""))

vegan_matrix <- vegan_species_data %>%
  select(-sample, -geo_loc_name)

vegan_matrix[is.na(vegan_matrix)] <- 0

# calculate the shannon diversity and add it back to the sample info
sample_info$shannon <- diversity(vegan_matrix, index = "shannon")

# plot this (alpha div)
ggplot(sample_info, aes(x = sample, y = shannon, fill = geo_loc_name)) +
  geom_col() +
  facet_wrap(~ geo_loc_name, scales = "free_x") +
  scale_fill_viridis_d(
    option = "A",
    begin = 0.1, end = 0.85) +
  labs(
    title = "Desert Biocrust Community Alpha Diversity",
    x = "Sample ID",
    y = "Shannon Diversity",
    fill = "Location"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))

 # beta diversity time
bray_dist <- vegdist(vegan_matrix, method = "bray")
bray_dist

pcoa <- cmdscale(bray_dist, k = 2)  

pcoa_df <- data.frame(                      
  pcoa,                                                 
  sample = sample_info$sample,
  geo_loc_name = sample_info$geo_loc_name                                                                                                                
) 

# plot for beta diversity
ggplot(pcoa_df, aes(x = X1, y = X2, color = geo_loc_name)) + 
  geom_point(size = 4) +
  geom_text(aes(label = sample), nudge_y = 0.03) +                                                                                                     
  labs(
      title = "PCoA of Biocrust Communities",
      x = "PCoA Axis 1",
      y = "PCoA Axis 2",
      color = NULL) +
  scale_color_viridis_d(
    option = "A",
    begin = 0.1, end = 0.85) +
  theme(legend.position = "bottom") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(title = NULL, nrow = 1))


# TODO: Make this good
#---- make a quick map so people can see where these samples are from ----

# Parse lat/lon from metadata 
map_info <- meta_data %>%
  separate(lat_lon, into = c("lat", "lon"), sep = " N ") %>%                                                                                             
  mutate(                                                                                                                                              
    lat = as.numeric(lat),
    lon = -parse_number(lon),
    sample = str_extract(title, "(?<=- ).*")
  ) %>%
  select(sample, lat, lon)

# Subsample map to show only states we want
us_map <- map_data("state") %>%
  filter(region %in% c("utah", "arizona", "new mexico"))

ggplot() +                                                                                                              
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group),                                                  
                fill = "white", color = "black") +                                                                       
  geom_point(data = map_info, aes(x = lon, y = lat), size = 4) +                                                        
  geom_text_repel(data = map_info, aes(x = lon, y = lat, label = sample),
                  size = 5, box.padding = .5)+
  theme_minimal()

# ---- Kraken2 of contigs: Circular genome taxonomy ----

# read in contig taxonomy from kraken2
contig_taxonomy <- read_tsv(
  "data/genome_classification/kraken2/contig_taxonomy.tsv",
  col_names = c("contig", "species", "taxid")
)

# read in checkm2 quality results
checkm2 <- read_tsv("data/genome_classification/checkm2/output/quality_report.tsv")

# join taxonomy with quality, filter to high quality
genome_summary <- contig_taxonomy %>%
  left_join(checkm2, by = c("contig" = "Name")) %>%
  filter(Completeness > 90, Contamination < 5) %>%
  mutate(sample = str_extract(contig, "^SRR\\d+")) %>%
  left_join(meta_data, by = c("sample" = "srr")) %>%
  mutate(sample_name = str_extract(title, "(?<=- ).*")) %>%
  mutate(sample = str_extract(contig, "^SRR\\d+"))

# Extract the coverage of the genomes so that we can use it in the plot
fasta_headers <- read_lines(pipe("grep '^>' data/genome_classification/kraken2/hq_chromosomes.fasta")) %>%
  keep(~ str_starts(.x, ">")) %>%                                                                                                           
  tibble(header = .) %>%                                                                                                                    
  mutate(                                                                                                                                   
    contig = str_extract(header, "(?<=>)\\S+"),                                                                                             
    coverage = as.numeric(str_extract(header, "(?<=coverage=)[\\d.]+"))) %>%
    select(contig, coverage) 

# Add coverage to the plot
genome_summary <- genome_summary %>%
  left_join(fasta_headers, by = "contig")

# map genus to phylum
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

# build plot data with phylum
genomes_within_sample <- genome_summary %>%
  select(sample, sample_name, species, geo_loc_name, coverage) %>%
  mutate(
    geo_loc_name = str_replace(geo_loc_name, "USA: ", ""),
    genus = str_extract(species, "^\\S+"),
    phylum = phylum_map[genus]
  )

# sort species by phylum then alphabetically
species_order <- genomes_within_sample %>%
  distinct(species, phylum) %>%
  arrange(phylum, species) %>%
  pull(species)


axis_colors <- genomes_within_sample %>%
  distinct(species, phylum) %>%
  arrange(match(species, species_order)) %>%
  mutate(color = phylum_colors[phylum])

genomes_within_sample <- genomes_within_sample %>%
  mutate(species = factor(species, levels = rev(species_order)))

ggplot(genomes_within_sample, aes(x = sample_name, y = species, color = phylum, size = coverage)) +
  geom_point() +
  scale_color_viridis_d(option = "A", begin = 0.1, end = 0.85) +
  facet_wrap(~ geo_loc_name, scales = "free_x") +
  theme_minimal() +
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

