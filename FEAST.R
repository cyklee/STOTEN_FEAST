# FEAST source tracking analysis

library("phyloseq")
library("tidyverse")
library("speedyseq")
library("FEAST") # remotes::install_github("cozygene/FEAST")
# https://www.nature.com/articles/s41592-019-0431-x - FEAST citation
# qs and writexl are used as well.

dir.create("1_FEAST")

# Load phyloseq object (including ASV table) serialised via qs
ps_2 <- qs::qread("ps_pub.qs") 

sample_data(ps_2)$CountryType <- paste(sample_data(ps_2)$Country, sample_data(ps_2)$Type, sep = "-")
ps_2_merged <- merge_samples(ps_2, group = "CountryType")

# Repopulate those variables by splitting the merged "sample names"
sample_data(ps_2_merged)$Type <- sample_names(ps_2_merged) %>% str_split("-") %>% sapply("[", 2)
sample_data(ps_2_merged)$Country <- sample_names(ps_2_merged) %>% str_split("-") %>% sapply("[", 1)

meta <- as_tibble(sample_data(ps_2_merged))
otu <- otu_table(ps_2_merged)

meta$SampleID <- meta$.sample
meta$SourceSink <- ifelse(meta$Type == "NG", "Sink", "Source")

# All source ID = NA
# Sink are ID'ed

# Previous master branch FEAST (pre-commit cf88568, Apr-2020) used rownames
# While newer FEAST can use SampleID column for indexing
# Doing both here for compatibility 

rownames(meta) <- meta$SampleID 

# Rework the id numbering
meta[meta$SourceSink == "Source", "id"] <- NA # The Sources are ungrouped (hence blank)
meta[meta$SourceSink == "Sink", "id"] <- 1:nrow(meta[meta$SourceSink == "Sink", "id"])

meta_FEAST <- meta %>% mutate(Env = SampleID)  %>% select(Env, SourceSink, id) %>% data.frame()
rownames(meta_FEAST) <- meta$SampleID

otu <- otu[rownames(otu) %in% meta$SampleID,] # Vestigial

# Many sinks, many individual sources
FEAST_out <- FEAST(
  C = as.matrix(otu),
  metadata = meta_FEAST,
  EM_iterations = 1000,
  dir_path = file.path(getwd(), "1_FEAST"), # Was undocumented, I think this changes wd
  outfile = "GAM3",                         # Was undocumented
  different_sources_flag = 0)               # Each sink can have the same source

str(FEAST_out)

#unknown_source & unknown_source_rarefy are per OTU/ASV
subset_var <- names(FEAST_out)[1:13]
FEAST_df <- as_tibble(FEAST_out[subset_var])

colnames(FEAST_df) <- str_split(colnames(FEAST_df), "_") %>% sapply(`[`, 1)
rownames(FEAST_df) <- colnames(FEAST_df)[seq_len(length(colnames(FEAST_df)) - 1)] %>%  str_replace("Soil", "NG")

# Output FEAST numeric results as Excel spreadsheet
writexl::write_xlsx(rownames_to_column(FEAST_df),
                    path = "1_FEAST/FEAST_output.xlsx",
                    format_headers = TRUE)