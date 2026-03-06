# load packages
library(tidyverse)
library(viridis)
library(maps)
library(vegan)
#setwd("projects/biocrust_metagenomics/scripts/")

# create master data frame
master_kraken2_data <- tibble()

# read in metadata
meta_data <- read.csv("../data/metadata/biocrust_metadata.csv")

# Fill master Kraken df
master_kraken2_data <- list.files("../data/kraken2/reports", full.names = TRUE) %>%                                                     
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

# TODO: make this look good
ggplot() +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group),
                fill = "white", color = "black") +
  geom_point(data = map_info, aes(x = lon, y = lat), size = 4) +
  theme_minimal()