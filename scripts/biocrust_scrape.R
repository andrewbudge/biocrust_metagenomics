library(rentrez)
library(xml2)
library(tidyverse)

setwd(dir = "~/projects/biocrust_metagenomics")

# Selected SRX accessions: 2 per site (Utah, New Mexico, Arizona)
accessions <- c(
  "SRX30858622",  # HS003  - Utah
  "SRX30858445",  # HSN023 - Utah
  "SRX30858626",  # CYAN3  - New Mexico
  "SRX30857901",  # SEV30  - New Mexico
  "SRX30858190",  # SON57  - Arizona
  "SRX30857843"   # SON60  - Arizona
)

# This is to clean up the search terms so rentrez likes it
query <- paste(paste0(accessions, "[accn]"), collapse = " OR ")

# Search the SRA database
search_result <- entrez_search(db = "sra", term = query, retmax = 20)

# Get all the ids
uids <- search_result$ids

# Fetch all data and parse it  
sra_xml <- entrez_fetch(db = "sra", id = uids, rettype = "xml")
sra_doc <- read_xml(sra_xml)
experiments <- xml_find_all(sra_doc, "//EXPERIMENT_PACKAGE")

srr_map <- lapply(experiments, function(exp) {
  srx <- xml_attr(xml_find_first(exp, ".//EXPERIMENT"), "accession")
  srr <- xml_attr(xml_find_first(exp, ".//RUN"), "accession")
  biosample_id <- xml_text(xml_find_first(exp, ".//EXTERNAL_ID[@namespace='BioSample']"))
  tibble(srx = srx, srr = srr, biosample = biosample_id)
})

srr_df <- bind_rows(srr_map)

# Get Biosample data and clean up data
biosample_xml <- entrez_fetch(db = "biosample", id = srr_df$biosample, rettype = "xml")
doc <- read_xml(biosample_xml)
samples <- xml_find_all(doc, "//BioSample")

sample_list <- lapply(samples, function(s) {
  attrs <- xml_find_all(s, ".//Attribute")
  tibble(
    biosample = xml_attr(s, "accession"),
    title = xml_text(xml_find_first(s, ".//Title")),
    attribute = xml_attr(attrs, "attribute_name"),
    value = xml_text(attrs)
  )
})

biosample_long <- bind_rows(sample_list)
biosample_wide <- biosample_long %>%
  pivot_wider(names_from = attribute, values_from = value)

final_table <- biosample_wide %>%
  left_join(srr_df, by = "biosample")

write_csv(final_table, "data/metadata/biocrust_metadata.csv")                                                                                                                         
