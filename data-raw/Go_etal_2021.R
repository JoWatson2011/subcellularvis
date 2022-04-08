library(subcellularvis)
library(tidyverse)

dat <- data.frame(
  rank = c(1:20),
  compartment = c(
    "Cell junction",
    "Chromatin",
    "ER membrane",
    "Mitochondrial matrix",
    "Actin, cytosol",
    "ER lumen",
    "Endosome, lysosome",
    "Nucleolus",
    "Nucleoplasm",
    "Nuclear body",
    "Plasma membrane",
    "Centrosome",
    "Mitochondrial Membrane, Peroxisome",
    "Golgi",
    "Nuclear outer membrane-ER membrane network",
    "Cytoplasmic RNP granule",
    "Microtubules",
    "Mitochondrial inner membrane, Mitochondrial inter membrane space",
    "Misc.",
    "Early endosome, Recycling endosome"
  )
)


prey <- readxl::read_excel("~/Supplementary table 8.xlsx", sheet = 2)  %>% 
  pivot_longer(-gene) %>% 
  group_by(gene) %>% 
  filter(value== max(value)) %>% 
  mutate(rank = as.integer(gsub("rank ", "", name))) %>% 
  select(-c(name, value)) %>% 
  ungroup() %>% 
  left_join(dat)

dat_n <- prey %>%
  group_by(compartment) %>%
  summarise(n_all=n())

subCell_predictions <- prey %>% 
  group_by(rank) %>% 
  group_split() %>% 
  lapply(function(i){
    compartmentData(i$gene)$enrichment %>% 
        slice(1) %>% 
      #filter(`FDR < 0.05` == T) %>% 
      mutate(predComp = Compartment,
             trueComp = unique(i$compartment))
  }) %>% 
  bind_rows() %>% 
  left_join(dat_n, by = c("trueComp" = "compartment")) %>% 
  mutate(n = paste0(n, "/", n_all)) %>% 
  select(trueComp, predComp, FDR, n)

subCell_predictions_supp <- prey %>% 
  group_by(rank) %>% 
  group_split() %>% 
  lapply(function(i){
    compartmentData(i$gene)$enrichment %>% 
    #  slice(1) %>% 
      #filter(`FDR < 0.05` == T) %>% 
      mutate(predComp = Compartment,
        trueComp = unique(i$compartment))
    }) %>% 
  bind_rows() %>% 
  filter(FDR < 0.05) %>% 
  left_join(dat_n, by = c("trueComp" = "compartment")) %>% 
  mutate(n = paste0(n, "/", n_all)) %>% 
  select(trueComp, predComp, FDR, n)

readr::write_csv(subCell_predictions, "data-raw/subcellResults_Go_etal_2021.csv")
readr::write_csv(subCell_predictions_supp, "data-raw/subcellResults_Go_etal_2021_supp.csv")


