#setwd("projects/biocrust_metagenomics/")
library(tidyverse)
library(ggrepel)

# read in extra
meta_data <- read.csv("data/metadata/biocrust_metadata.csv")
alpha_div <- read_tsv("data/metadata/alpha_diversity.tsv")

# ---- Phage presence and absence ----
# join metadata to viral data and filter
viral_data <- read_tsv("data/genomes/viral/hq_phage_groups.tsv") %>%                                                              
  left_join(meta_data, by = c("sample" = "srr")) %>% 
  mutate(
    geo_loc_name = str_replace(as.character(geo_loc_name), "USA: ", ""),                                                                                
    site = str_extract(title, "(?<=- ).*")) %>%
  select(group, sample, contig, geo_loc_name, site)

# get count for heatmap tiles
viral_counts <- viral_data %>% count(group, site, geo_loc_name)

# make plot showing what phages appear where
ggplot(viral_data, aes(x = site, y = group, fill = group)) +                                                                                     
  geom_tile(color = "white") +
  facet_wrap(~ geo_loc_name, scales = "free_x")+
  scale_fill_viridis_d(option = "A", begin = 0.1, end = 0.85) +
  labs(
    title = "Phage Presence in Desert Biocrust",
    x = "Sample site",
    y = NULL
  ) +
  theme(
    legend.position = "none",
    title = element_text(size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)   
    )


# ---- Phage and overall diversity comparison ----
phage_richness <- viral_data %>%
  count(site, name = "phage_count") %>%
  left_join(alpha_div, by = c("site" = "sample"))

ggplot(phage_richness, aes(x = shannon, y = phage_count, color = geo_loc_name)) +                                                                                             
  geom_point(size = 3) +                                                                                                                                
  geom_text_repel(aes(label = site), size = 5, show.legend = F) +                                                                                                         
  labs(     
    title = "Phage Richness and Bacterial Diversity",                                                                                                                                            
    x = "Shannon Diversity",
    y = "Phage Richness",
    color = "State"                                                                                                                                
  ) +
  scale_color_viridis_d(option = "A", begin = 0.1, end = 0.75) + 
  theme(legend.position = "bottom",
        legend.text = element_text(size = 14),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)
      )
 



