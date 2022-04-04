library(org.Hs.eg.db)
library(GO.db)
library(dplyr)
library(tidyr)
child <- as.list(GOCCOFFSPRING)

### #####################################
###"Top level"/summarised annotations,
### used for compartments visualisation
#########################################
#Trafficking
COMPARTMENTS_traffic <- data.frame(compartment = c("Cytoplasm", "Plasma membrane", 
                                                  "Intracellular vesicle", "Lysosome",
                                                  "Recycling Endosome", "Endolysosome",
                                                  "Late Endosome", "Early Endosome"
),
ID = c("GO:0005737", "GO:0005886", "GO:0097708",
       "GO:0005764", "GO:0055037", "GO:0036019",
       "GO:0005770", "GO:0005769"),
stringsAsFactors = F)

COMPARTMENTS_traffic_offspring <- child[names(child) %in% COMPARTMENTS_traffic$ID] 
COMPARTMENTS_traffic_offspring <- lapply(names(COMPARTMENTS_traffic_offspring), function(i){
  data.frame(ID = i, subid = COMPARTMENTS_traffic_offspring[[i]])
}) %>% 
  do.call(rbind, .) %>% 
  merge(COMPARTMENTS_traffic, by = "ID") %>% 
  pivot_longer(cols = -compartment, names_to = "level", values_to = "ID") %>% 
  dplyr::select(compartment, ID) %>% 
  unique()
  
Human_traffic_annots <- AnnotationDbi::select(org.Hs.eg.db, 
                                         COMPARTMENTS_traffic_offspring$ID,
                                         c("GO", "SYMBOL", "UNIPROT"),
                                         "GO") %>% 
  dplyr::select(GO, SYMBOL, UNIPROT) %>% 
  unique() %>% 
  merge(COMPARTMENTS_traffic_offspring, by.x = "GO", by.y = "ID") %>% 
  dplyr::select(compartment, SYMBOL, UNIPROT) %>% 
  distinct()

### Whole cell 
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

COMPARTMENTS_offspring <- child[names(child) %in% COMPARTMENTS_parent$ID] 
COMPARTMENTS_offspring <- lapply(names(COMPARTMENTS_offspring), function(i){
  data.frame(ID = i, subid = COMPARTMENTS_offspring[[i]])
}) %>% 
  do.call(rbind, .) %>% 
  merge(COMPARTMENTS_parent, by = "ID") %>% 
  pivot_longer(cols = -compartment, names_to = "level", values_to = "ID") %>% 
  dplyr::select(compartment, ID) %>% 
  unique()

Human_annots <- AnnotationDbi::select(org.Hs.eg.db, 
                                        COMPARTMENTS_offspring$ID,
                                        c("GO", "SYMBOL", "UNIPROT"),
                                        "GO") %>% 
  dplyr::select(GO, SYMBOL, UNIPROT) %>% 
  distinct() %>% 
  merge(COMPARTMENTS_offspring, by.x = "GO", by.y = "ID") %>% 
  dplyr::select(compartment, SYMBOL, UNIPROT) %>% 
  distinct()

#########################
### Detailed annotations
#########################

#Trafficking
traffic_subannots_terms <- AnnotationDbi::select(
  GO.db,
  COMPARTMENTS_traffic_offspring$ID,
  "TERM",
  "GOID"
) %>% 
  distinct() %>% 
  left_join(COMPARTMENTS_traffic_offspring, by = c("GOID"="ID")) %>% 
  distinct()

Human_traffic_subannots <- AnnotationDbi::select(
  org.Hs.eg.db,
  traffic_subannots_terms$GOID,
  c("SYMBOL","UNIPROT"),
  "GO"
) %>% 
  left_join(traffic_subannots_terms, by = c("GO" = "GOID")) %>% 
  dplyr::select(GOID = GO, SYMBOL, UNIPROT,
                compartment = TERM, group = compartment)  %>% 
  distinct()

# Whole cell
subannots_terms <- AnnotationDbi::select(
  GO.db,
  COMPARTMENTS_offspring$ID,
  "TERM",
  "GOID"
) %>% 
  distinct() %>% 
  left_join(COMPARTMENTS_offspring, by = c("GOID"="ID")) %>% 
  distinct()

Human_subannots <- AnnotationDbi::select(
  org.Hs.eg.db,
  subannots_terms$GOID,
  c("SYMBOL","UNIPROT"),
  "GO"
) %>% 
  left_join(subannots_terms, by = c("GO" = "GOID")) %>% 
  dplyr::select(GOID = GO, SYMBOL, UNIPROT,
                compartment = TERM, group = compartment)  %>% 
  distinct()

#################################
# Human Protein Atlas annotations
#################################

hpa <- readr::read_tsv("data-raw/proteinatlas_subcellular_location.tsv")

## Some of the GO terms in hpa aren't
# child terms to those specified above
# so I have manually assigned them / 
# used HPA website for categories below 

orphanHPAAnnots <- data.frame(
  stringsAsFactors = FALSE,
  GOTerm = c("Aggresome (GO:0016235)",
             "Cell Junctions (GO:0030054)",
             "Focal adhesion sites (GO:0005925)",
             "Kinetochore (GO:0000776)",
             "Lipid droplets (GO:0005811)","Midbody (GO:0030496)",
             "Midbody ring (GO:0090543)",
             "Mitotic chromosome (GO:0005694)",
             "Rods & Rings ()",
             "Vesicles (GO:0043231)"),
  # GOID = c("GO:0016235","GO:0030054",
  #          "GO:0005925","GO:0000776","GO:0005811","GO:0030496",
  #          "GO:0090543","GO:0005694",NA,"GO:0043231"),
  NewTerm = c("Cytoplasm", "Plasma membrane",
              "Plasma membrane", "Cytoplasm",
              NA, "Cytoskeleton", "Cytoplasm",
              "Cytoplasm", "Cytoplasm",
              "Intracellular vesicles"),
  n = c(19L, 313L, 129L, 5L, 39L, 53L, 30L, 69L, 19L, 1983L)
)

hpa_tidy <- hpa %>% 
  filter(Reliability != "Uncertain") %>% 
  select(Gene, `Gene name`, GOTerm = `GO id`) %>% 
  separate_rows(GOTerm, sep = ";") %>% 
  mutate(GOID = gsub("^.*\\(", "", gsub(")$", "", GOTerm))
         ) %>% 
  left_join(COMPARTMENTS_offspring, by = c("GOID" = "ID")) %>% 
  left_join(orphanHPAAnnots, by = c("GOTerm")) %>% 
  mutate(compartment = ifelse(is.na(compartment), NewTerm, compartment)) 

Human_HPA_annots <- hpa_tidy %>% 
  select(SYMBOL = `Gene name`, compartment)
Human_HPA_annots <- AnnotationDbi::select(org.Hs.eg.db, 
                                          Human_HPA_annots$SYMBOL,
                                            c("UNIPROT"),
                                            "SYMBOL") %>% 
  left_join(Human_HPA_annots, by = "SYMBOL") %>% 
  distinct()

