library(subcellularvis)
library(tidyverse)
#library(org.Sc.sgd.db)

dat <- readxl::read_excel("~/Nightingale2019_supp3.xlsx",
                          col_names = c("UNIPROT", "Protein name", 
                                        "Compartment", "SVM"))
# mapGenes <- AnnotationDbi::select(org.Sc.sgd.db,
#                                   dat$UNIPROT,
#                                   "GENENAME",
#                                   "UNIPROT")
# prey <- left_join(dat, mapGenes, by = "UNIPROT") %>% 
#   dplyr::select(gene = GENENAME, Compartment)

dat_n <- dat %>%
  select(UNIPROT, Compartment) %>% 
  group_by(Compartment) %>%
  summarise(n_all=n())

subCell_predictions <- dat %>% 
  group_by(Compartment) %>% 
  group_split() %>% 
  lapply(function(i){
    compartmentData(i$UNIPROT,
                    organism = "Yeast",
                    id_type = "UNIPROT")$enrichment %>% 
     # dplyr::slice(1) %>% 
     # filter(`FDR < 0.05` == T) %>% 
      mutate(trueComp = unique(i$Compartment))
  }) %>% 
  bind_rows() %>% 
  left_join(dat_n, by = c("trueComp" = "Compartment")) %>% 
  mutate(n = paste0(n, "/", n_all)) %>% 
  dplyr::select(trueComp, predComp = Compartment, FDR, n)

loc_supp <- dat %>% 
  group_by(Compartment) %>% 
  group_split() %>% 
  lapply(function(i)
  # compartmentData(i)[1:2,] ) %>% 
    compartmentData(i$UNIPROT,
                    organism = "Yeast",
                    id_type = "UNIPROT")$enrichment %>% 
    mutate(predComp = Compartment,
           trueComp = unique(i$Compartment))
) %>% 
  bind_rows() %>% 
  # mutate(trueComp = as.vector(sapply(names(loc_genes), rep, 2))) %>% 
  dplyr::select(trueComp, predComp, FDR, n) %>% 
  left_join(dat_n, c("trueComp" = "Compartment")) %>% 
  mutate(n = paste0(n, "/", n_all)) %>% 
  dplyr::select(-n_all) %>% 
  filter(FDR < 0.05)

readr::write_csv(subCell_predictions, 
                 "data-raw/subcellResults_Nightingale_etal_2019.csv")
readr::write_csv(loc_supp, 
                 "data-raw/subcellResults_Nightingale_etal_2019_supp.csv")
