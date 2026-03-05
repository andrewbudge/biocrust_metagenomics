# load packages
library(tidyverse)
library(viridis)
library(gganimate)
setwd("projects/biocrust_metagenomics/scripts/")

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
ggsave("../plots/phylum_community_barplot.png", width = 10, height = 6, dpi = 300)
ggsave("../plots/phylum_community_barplot.pdf", width = 10, height = 6)

