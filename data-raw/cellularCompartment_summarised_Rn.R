library(org.Rn.eg.db)
library(GO.db)
library(dplyr)
library(tidyr)
child <- as.list(GOCCOFFSPRING)

### #####################################
###"Top level"/summarised annotations,
### used for compartments visualisation
#########################################
#Trafficking
COMPARTMENTS_traffic <-
  data.frame(
    compartment = c(
      "Cytoplasm",
      "Plasma membrane",
      "Intracellular vesicle",
      "Lysosome",
      "Recycling Endosome",
      "Endolysosome",
      "Late Endosome",
      "Early Endosome"
    ),
    ID = c(
      "GO:0005737",
      "GO:0005886",
      "GO:0097708",
      "GO:0005764",
      "GO:0055037",
      "GO:0036019",
      "GO:0005770",
      "GO:0005769"
    ),
    stringsAsFactors = F
  )

COMPARTMENTS_traffic_offspring <-
  child[names(child) %in% COMPARTMENTS_traffic$ID]
COMPARTMENTS_traffic_offspring <-
  lapply(names(COMPARTMENTS_traffic_offspring), function(i) {
    data.frame(ID = i, subid = COMPARTMENTS_traffic_offspring[[i]])
  }) %>%
  do.call(rbind, .) %>%
  merge(COMPARTMENTS_traffic, by = "ID") %>%
  pivot_longer(cols = -compartment,
               names_to = "level",
               values_to = "ID") %>%
  dplyr::select(compartment, ID) %>%
  unique()

Rat_traffic_annots <- AnnotationDbi::select(org.Rn.eg.db,
                                            COMPARTMENTS_traffic_offspring$ID,
                                            c("GO", "SYMBOL"),
                                            "GO") %>%
  dplyr::select(GO, SYMBOL) %>%
  unique() %>%
  merge(COMPARTMENTS_traffic_offspring,
        by.x = "GO",
        by.y = "ID") %>%
  dplyr::select(SYMBOL, compartment) %>%
  na.omit() %>%
  distinct()


### Whole cell
COMPARTMENTS_parent <-
  data.frame(
    compartment = c(
      "Nucleus",
      "Cytoplasm",
      "Cytoskeleton",
      "Peroxisome",
      "Vacuole",
      "Endoplasmic reticulum",
      "Golgi apparatus",
      "Plasma membrane",
      "Endosome",
      "Extracellular region",
      "Mitochondrion",
      "Intracellular vesicle",
      "Ribosome",
      "Lysosome"
    ),
    ID = c(
      "GO:0005634",
      "GO:0005737",
      "GO:0005856",
      "GO:0005777",
      "GO:0005773",
      "GO:0005783",
      "GO:0005794",
      "GO:0005886",
      "GO:0005768",
      "GO:0005576",
      "GO:0005739",
      "GO:0097708",
      "GO:0005840",
      "GO:0005764"
    ),
    stringsAsFactors = F
  )

COMPARTMENTS_offspring <-
  child[names(child) %in% COMPARTMENTS_parent$ID]
COMPARTMENTS_offspring <-
  lapply(names(COMPARTMENTS_offspring), function(i) {
    data.frame(ID = i, subid = COMPARTMENTS_offspring[[i]])
  }) %>%
  do.call(rbind, .) %>%
  merge(COMPARTMENTS_parent, by = "ID") %>%
  pivot_longer(cols = -compartment,
               names_to = "level",
               values_to = "ID") %>%
  dplyr::select(compartment, ID) %>%
  unique()

Rat_annots <- AnnotationDbi::select(org.Rn.eg.db,
                                    COMPARTMENTS_offspring$ID,
                                    c("GO", "SYMBOL"),
                                    "GO") %>%
  dplyr::select(GO, SYMBOL) %>%
  unique() %>%
  merge(COMPARTMENTS_offspring, by.x = "GO", by.y = "ID") %>%
  dplyr::select(SYMBOL, compartment) %>%
  na.omit() %>%
  distinct()


########################################
### Detailed annotations
###
#########################################

#Trafficking
traffic_subannots_terms <- AnnotationDbi::select(GO.db,
                                                 COMPARTMENTS_traffic_offspring$ID,
                                                 "TERM",
                                                 "GOID") %>%
  distinct() %>%
  left_join(COMPARTMENTS_traffic_offspring, by = c("GOID" = "ID")) %>%
  distinct()

Rat_traffic_subannots <- AnnotationDbi::select(org.Rn.eg.db,
                                               traffic_subannots_terms$GOID,
                                               "SYMBOL",
                                               "GO") %>%
  left_join(traffic_subannots_terms, by = c("GO" = "GOID")) %>%
  dplyr::select(GOID = GO,
                SYMBOL,
                compartment = TERM,
                group = compartment) %>%
  distinct()


# Whole cell
subannots_terms <- AnnotationDbi::select(GO.db,
                                         COMPARTMENTS_offspring$ID,
                                         "TERM",
                                         "GOID") %>%
  distinct() %>%
  left_join(COMPARTMENTS_offspring, by = c("GOID" = "ID")) %>%
  distinct()

Rat_subannots <- AnnotationDbi::select(org.Rn.eg.db,
                                       subannots_terms$GOID,
                                       "SYMBOL",
                                       "GO") %>%
  left_join(subannots_terms, by = c("GO" = "GOID")) %>%
  dplyr::select(GOID = GO,
                SYMBOL,
                compartment = TERM,
                group = compartment) %>%
  distinct()