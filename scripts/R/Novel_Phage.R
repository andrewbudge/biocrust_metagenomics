library(tidyverse)
library(ape)
library(phangorn)

# ---- Read in data ----

viptree_lst <- read_tsv("data/genomes/viral/phylogenomics/genome_table.lst",
                        show_col_types = FALSE)
viridic_sim <- read_tsv("data/genomes/viral/phylogenomics/VIRIDIC_sim-dist_table.tsv",
                        show_col_types = FALSE)

# clean up HTML column names
viptree_lst <- viptree_lst %>%
  rename(
    SG_Lepyrophagus = `<i>S</i><sub>G</sub> to Lepyrophagus`,
    SG_Eremotheras = `<i>S</i><sub>G</sub> to Eremotheras`,
    SG_Ammolykos = `<i>S</i><sub>G</sub> to Ammolykos`
  )

novel_phages <- c("Eremotheras", "Ammolykos", "Lepyrophagus")

# ---- VIRIDIC Heatmap ----

# Pivot to long format for ggplot
viridic_long <- viridic_sim %>%
  pivot_longer(-genome, names_to = "genome2", values_to = "similarity") %>%
  mutate(
    group1 = str_extract(genome, "^[^_]+"),
    group2 = str_extract(genome2, "^[^_]+"),
    # Shorten labels for readability
    label1 = str_replace(genome, "_(SRR\\d+)_ctg(\\d+)_provirus_(\\d+)_(\\d+)", "\n\\1 ctg\\2"), # this is horrible but it works so I will not change it.
    label2 = str_replace(genome2, "_(SRR\\d+)_ctg(\\d+)_provirus_(\\d+)_(\\d+)", "\n\\1 ctg\\2")
  )

# Order by group
genome_order <- viridic_sim$genome
viridic_long <- viridic_long %>%
  mutate(
    label1 = factor(label1, levels = str_replace(genome_order, "_(SRR\\d+)_ctg(\\d+)_provirus_(\\d+)_(\\d+)", "\n\\1 ctg\\2")), # this is also horrible
    label2 = factor(label2, levels = str_replace(genome_order, "_(SRR\\d+)_ctg(\\d+)_provirus_(\\d+)_(\\d+)", "\n\\1 ctg\\2")) # I hate regex
  )

ggplot(viridic_long, aes(x = label1, y = label2, fill = similarity)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(similarity, 1), color = similarity > 50),
            size = 5, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white")) +
  scale_fill_viridis_c(option = 'A', begin = 0.1, end = 0.85) +
  labs(
    title = "VIRIDIC Intergenomic Similarity",
    x = NULL, y = NULL,
    fill = "Similarity (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

# ---- ViPTree Similarity Lollipop Plot ----

sg_long <- bind_rows(
  viptree_lst %>% filter(!ID %in% novel_phages) %>%
    slice_max(SG_Eremotheras, n = 10) %>%
    select(ID, name, Hgroup, SG = SG_Eremotheras) %>%
    mutate(query = "Eremotheras"),
  viptree_lst %>% filter(!ID %in% novel_phages) %>%
    slice_max(SG_Ammolykos, n = 10) %>%
    select(ID, name, Hgroup, SG = SG_Ammolykos) %>%
    mutate(query = "Ammolykos"),
  viptree_lst %>% filter(!ID %in% novel_phages) %>%
    slice_max(SG_Lepyrophagus, n = 10) %>%
    select(ID, name, Hgroup, SG = SG_Lepyrophagus) %>%
    mutate(query = "Lepyrophagus")
) %>%
  mutate(Hgroup = ifelse(is.na(Hgroup) | Hgroup == "-", "Unknown", Hgroup))


 ggplot(sg_long, aes(x = SG, y = reorder(name, SG), color = Hgroup)) +
    geom_point(size = 3) +
    geom_segment(aes(xend = 0, yend = reorder(name, SG)), linewidth = 0.5) +
    facet_wrap(~ query, scales = "free_y") +                                                            
    geom_vline(xintercept = 0.70, linetype = "dashed", alpha = 0.4) +
    annotate("text", x = 0.72, y = 1, label = "Genus\nthreshold", size = 2.5, hjust = 0, alpha = 0.5) + 
    scale_color_viridis_d(option = "A", begin = 0.1, end = 0.85) +                                                              
    labs(                                                                                               
      title = "Closest Known Relatives (ViPTree SG Scores)",                                            
      x = "Genome Similarity (SG)",                         
      y = NULL,
      color = "Host Phylum"
    ) +                                                                                                 
    theme_minimal(base_size = 12) +
    theme(                                                                                              
      plot.title = element_text(hjust = 0.5),               
      legend.position = "bottom",
      strip.text = element_text(face = "bold.italic", size = 12)
    ) 
