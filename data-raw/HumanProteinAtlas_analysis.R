library(tidyverse)
library(subcellularvis)
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!group_keys(grouped), sep = " / "))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}
#PROvision, PROlocalise

loc <- readr::read_tsv("data-raw/proteinatlas_subcellular_location.tsv")

loc_tidy <- loc %>% 
  filter(Reliability != "Uncertain") %>% 
  select(`Gene name`,`Main location`, `Additional location`) %>%  
  separate(col = `Main location`, 
           into = c("Main location","addLoc"),
           extra = "merge", sep = ";") %>% 
  rowwise() %>% 
  mutate(`Other locations` = 
           paste(
             na.omit(c(`addLoc`, `Additional location`)), 
             collapse = ";"
           )
         ) %>% 
  select(-c(`addLoc`, `Additional location`))# %>% 
  # group_by(`Main location`, `Other locations`) %>% 
  # summarise(n = n(), .groups = "keep") %>% 
  # arrange(desc(n))

  
loc_genes <- loc_tidy %>% 
  filter(`Other locations` == "") %>% 
  distinct() %>%
  named_group_split(`Main location`) %>% 
  map(function(i){
    unique(i$`Gene name`)
  })
loc_genes <- loc_genes[sapply(loc_genes, length) > 5]

loc_n <- loc_tidy %>% 
  filter(`Other locations` == "") %>% 
  select(-`Other locations`) %>% 
  distinct() %>% 
  group_by(`Main location`) %>%
  summarise(nGenes = n(), .groups = "keep")

loc_res <- lapply(loc_genes, function(i)
 # compartmentData(i)[1:2,] ) %>% 
  compartmentData(i)$enrichment[1,] 
 ) %>% 
  bind_rows() %>% 
 mutate(trueComp = names(loc_genes)) %>% 
 # mutate(trueComp = as.vector(sapply(names(loc_genes), rep, 2))) %>% 
  select(trueComp, predComp = Compartment, FDR, n) %>% 
  left_join(loc_n, c("trueComp" = "Main location")) %>% 
  mutate(n = paste0(n, "/", nGenes)) %>% 
  select(-nGenes)

loc_supp <- lapply(1:length(loc_genes), function(i)
  # compartmentData(i)[1:2,] ) %>% 
  compartmentData(loc_genes[[i]])$enrichment %>% 
    mutate(predComp = Compartment,
           trueComp = names(loc_genes)[i])
) %>% 
  bind_rows() %>% 
  # mutate(trueComp = as.vector(sapply(names(loc_genes), rep, 2))) %>% 
  select(trueComp, predComp, FDR, n) %>% 
  left_join(loc_n, c("trueComp" = "Main location")) %>% 
  mutate(n = paste0(n, "/", nGenes)) %>% 
  select(-nGenes) %>% 
  filter(FDR < 0.05)

write_csv(loc_res, "data-raw/subcellResults_HPA.csv")
write_csv(loc_supp, "data-raw/subcellResults_HPA_supp.csv")