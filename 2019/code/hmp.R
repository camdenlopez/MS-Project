# Creates datasets from the Human Microbiome Project
# sample map and OTU table.

# Data source and documentation:
# https://www.hmpdacc.org/HMQCP/
# https://www.hmpdacc.org/hmp/resources/metagenomics_sequencing_analysis.php

library(tidyverse)
options(stringsAsFactors = FALSE)

# RSID = random subject identifier
map <- read_tsv("../data/v13_map_uniquebyPSN.txt") %>%
  select(SampleID = `#SampleID`, RSID, visitno, sex, RUNCENTER, HMPbodysubsite)

samples <- map %>%
  filter(HMPbodysubsite == "Tongue_dorsum") %>%
  group_by(RSID) %>%
  filter(visitno == min(visitno))

tbl <- read_tsv("../data/otu_table_psn_v13.txt", skip = 1) %>%
  rename(OTU = `#OTU ID`, Lineage = `Consensus Lineage`) %>%
  select(c("OTU", "Lineage", as.character(samples$SampleID))) %>%
  gather(key = "SampleID", value = "Count",
         -OTU, -Lineage)

otus <- tbl %>%
  group_by(OTU) %>%
  summarize(phylum = sub("p__", "", str_extract(Lineage[1], "p__[^;]+")),
            class = sub("c__", "", str_extract(Lineage[1], "c__[^;]+")),
            order = sub("o__", "", str_extract(Lineage[1], "o__[^;]+")),
            family = sub("f__", "", str_extract(Lineage[1], "f__[^;]+")),
            genus = sub("g__", "", str_extract(Lineage[1], "g__[^;]+")),
            presence = mean(Count > 0)) %>%
  arrange(desc(presence))

hmp.sm <- tbl %>%
  merge(otus %>%
          filter(presence >= 0.90)) %>%
  select(OTU, SampleID, Count)

hmp.md <- tbl %>%
  merge(otus %>%
          filter(presence >= 0.80)) %>%
  select(OTU, SampleID, Count)

hmp.lg <- tbl %>%
  merge(otus %>%
          filter(presence >= 0.70)) %>%
  select(OTU, SampleID, Count)

save(samples, otus,
     hmp.sm, hmp.md, hmp.lg,
     file = "../data/hmp.RData")
