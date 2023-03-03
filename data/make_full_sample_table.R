library(tidyverse)


lkup <- read_csv("SraRunTable.txt") %>%
  dplyr::select(Run,GSM=`GEO_Accession (exp)`,sample_preparation_protocol)


df <- read_csv("tmp_subsample_table.csv")


df <- df %>% 
  left_join(lkup, by = c(sample_name = "Run"))

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162192
geo_lkp <- tribble(~geo_sample_name, ~GSM,
              "adult1", "GSM4946261",
              "adult2", "GSM4946262",
              "adult3", "GSM4946263",
              "adult4", "GSM4946264",
              "adult5", "GSM4946265")

df <- df %>% 
  left_join(geo_lkp) %>%
  mutate(library_type = case_when(sample_preparation_protocol == "protocol A" ~ "10xV2",
                                  sample_preparation_protocol == "protocol B" ~ "10xV3"))


df <- df %>% 
  mutate(ALEVIN_MAP_PARAM = case_when(sample_preparation_protocol == "protocol A" ~ "--chromium -l IU",
                                       sample_preparation_protocol == "protocol B" ~ "--chromiumV3 -l IU"))

write_csv(df,"../config/subsample_table.csv")