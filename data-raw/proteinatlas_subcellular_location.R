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

nuc_loc <- loc %>% 
  filter(Reliability != "Uncertain") %>% 
  select(Gene, `Gene name`, `Main location`,
         `Additional location`, `GO id`) %>% 
  pivot_longer(c(`Main location`, `Additional location`), 
              names_to = "type", values_to = "location") %>% 
  separate_rows(location, sep = ";") %>% 
  filter(grepl("[Nn]uc", location)) 

nuc_go <- nuc_loc %>% 
  separate_rows(`GO id`, sep = ";") %>%
  separate(`GO id`, c("GO term", "GO id"), sep = " \\(") %>% 
  mutate(`GO id` = gsub(")$", "", `GO id`)
  ) %>% 
  filter(`GO term` == location) %>% 
  select(`GO term`, `GO id`) %>% 
  unique()

nuc_loc %>% 
  group_by(type, location) %>% 
  summarise(n = n(), .groups = "keep") %>% 
  ggplot(aes(x = location, y = n, fill = type))+
  geom_col() +
  theme_bw() +
  xlab("Nuclear location") +
  scale_fill_viridis_d() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tmp <- nuc_loc %>% 
  group_by(Gene) %>% 
  mutate(type = paste0("node", 1:n())) %>% 
  select(`Gene`, location) %>% 
  group_split() %>% 
  lapply(function(i){
    data.frame(
      Gene = i$Gene[1],
      expand.grid(i$location, i$location)
    ) %>% 
      filter(Var1 != Var2)
  }) %>% 
  bind_rows() %>% 
  ungroup() %>% 
  group_by(Var1,Var2) %>% 
  summarise(n = n(), .groups = "keep")
  
library(igraph)  
library(ggraph)
mygraph <- graph_from_data_frame(tmp, directed = F)
E(mygraph)$weight <-  tmp$n
ggraph(mygraph, layout = 'linear') + 
  geom_edge_arc(aes(edge_width = n, alpha = n),
                end_cap = circle(2, 'mm'),
                start_cap = circle(2, 'mm')) +
  geom_node_label(aes(label = name))

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

tmp <- lapply(loc_genes, function(i)
 # compartmentData(i)[1:2,] ) %>% 
  compartmentData(i)[1,] ) %>% 
  bind_rows() %>% 
 mutate(trueComp = names(loc_genes)) %>% 
 # mutate(trueComp = as.vector(sapply(names(loc_genes), rep, 2))) %>% 
  select(trueComp, predComp = Compartment, FDR, n) %>% 
  left_join(loc_n, c("trueComp" = "Main location")) %>% 
  mutate(n = paste0(n, "/", nGenes)) %>% 
  select(-nGenes)




write_csv(tmp, "~/subcell.csv")

