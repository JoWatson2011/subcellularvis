library(org.Hs.eg.db)
library(GO.db)
library(dplyr)
library(tidyr)

# Child term list from GO.db package
child <- as.list(GOCCOFFSPRING)

####
# Endolysosome system
####

# Reference dataframe containg "top level" GO term and ID.
# Found on https://supfam.org/SUPERFAMILY/cgi-bin/dcgo.cgi
COMPARTMENTS_traffic <- data.frame(compartment = c("Cytoplasm", "Plasma membrane", 
                                                  "Intracellular vesicle", "Lysosome",
                                                  "Recycling Endosome", "Endolysosome",
                                                  "Late Endosome", "Early Endosome"
),
ID = c("GO:0005737", "GO:0005886", "GO:0097708",
       "GO:0005764", "GO:0055037", "GO:0036019",
       "GO:0005770", "GO:0005769"),
stringsAsFactors = F)

# Subset child terms for those in reference dataframe
COMPARTMENTS_traffic_offspring <- child[names(child) %in% COMPARTMENTS_traffic$ID] 

# Convert lists -> data.frame and merge with original reference data.frame.
# Then pivot_longer to keep the "top level" term name and all child terms
COMPARTMENTS_traffic_offspring <- lapply(names(COMPARTMENTS_traffic_offspring), 
                                         function(i){
  data.frame(ID = i, subid = COMPARTMENTS_traffic_offspring[[i]])
}) %>% 
  do.call(rbind, .) %>% 
  merge(COMPARTMENTS_traffic, by = "ID") %>% 
  pivot_longer(cols = -compartment, names_to = "level", values_to = "ID") %>% 
  dplyr::select(compartment, ID) %>% 
  unique()

# Find all the HGNC gene symbols associated to the GOID  
traffic_annots <- AnnotationDbi::select(org.Hs.eg.db, 
                                         COMPARTMENTS_traffic_offspring$ID,
                                         c("GO", "SYMBOL"),
                                         "GO") %>% 
  dplyr::select(GO, SYMBOL) %>% 
  unique() %>% 
  merge(COMPARTMENTS_traffic_offspring, by.x = "GO", by.y = "ID") %>% 
  dplyr::select(SYMBOL, compartment) %>% 
  na.omit() %>% 
  unique()

# Export
saveRDS(traffic_annots, "data/traffic_annots.rds")

####
# Whole cell
####

# Reference dataframe containg "top level" GO term and ID.
# Found on https://supfam.org/SUPERFAMILY/cgi-bin/dcgo.cgi
COMPARTMENTS_parent <- data.frame(compartment = c("Nucleus", "Cytoplasm", "Cytoskeleton",
                                                  "Peroxisome", "Vacuole", "Endoplasmic reticulum",
                                                  "Golgi apparatus", "Plasma membrane", "Endosome",
                                                  "Extracellular region", "Mitochondrion",
                                                  "Intracellular vesicle", "Ribosome", "Lysosome"
),
ID = c("GO:0005634", "GO:0005737", "GO:0005856",
       "GO:0005777", "GO:0005773", "GO:0005783",
       "GO:0005794", "GO:0005886", "GO:0005768",
       "GO:0005576", "GO:0005739", "GO:0097708",
       "GO:0005840", "GO:0005764"
),
stringsAsFactors = F)

# Subset child terms for those in reference dataframe
COMPARTMENTS_offspring <- child[names(child) %in% COMPARTMENTS_parent$ID] 

# Convert lists -> data.frame and merge with original reference data.frame.
# Then pivot_longer to keep the "top level" term name and all child terms
COMPARTMENTS_offspring <- lapply(names(COMPARTMENTS_offspring), function(i){
  data.frame(ID = i, subid = COMPARTMENTS_offspring[[i]])
}) %>% 
  do.call(rbind, .) %>% 
  merge(COMPARTMENTS_parent, by = "ID") %>% 
  pivot_longer(cols = -compartment, names_to = "level", values_to = "ID") %>% 
  dplyr::select(compartment, ID) %>% 
  unique()

# Find all the HGNC gene symbols associated to the GOID  
annots <- AnnotationDbi::select(org.Hs.eg.db, 
                                        COMPARTMENTS_offspring$ID,
                                        c("GO", "SYMBOL"),
                                        "GO") %>% 
  dplyr::select(GO, SYMBOL) %>% 
  unique() %>% 
  merge(COMPARTMENTS_offspring, by.x = "GO", by.y = "ID") %>% 
  dplyr::select(SYMBOL, compartment) %>% 
  na.omit() %>% 
  unique()

# Export
saveRDS(annots, "data/annots.rds")
