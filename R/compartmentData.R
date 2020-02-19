#' Annotate genes with filtered Gene Ontology (GO)
#' Cellular Compartment (CC) terms.
#'
#' @param genes A vector of HGNC gene symbols
#' @return A dataframe with `genes` annotated to cellular compartment
#' @export
#' @examples
#' genes <- c("MAPK1", "MAPK2", "MAPK3")
#' compartmentData(genes)

compartmentData <- function(genes){

  # Identify 'parent' GO terms + IDs.
  # http://www.supfam.org/SUPERFAMILY/cgi-bin/dcgo.cgi
  COMPARTMENTS_parent <- data.frame(compartment = c("Nucleus", "Cytoplasm", "Cytoskeleton",
                                                    "Peroxisome", "Vacuole", "Endoplasmic reticulum",
                                                    "Golgi apparatus", "Plasma membrane", "Endosome",
                                                    "Extracellular region", "Mitochondrion",
                                                    "Intracellular vesicle", "Ribosome", "Lysosome"
  ),
  GOID = c("GO:0005634", "GO:0005737", "GO:0005856",
           "GO:0005777", "GO:0005773", "GO:0005783",
           "GO:0005794", "GO:0005886", "GO:0005768",
           "GO:0005576", "GO:0005739", "GO:0097708",
           "GO:0005840", "GO:0005764"
  ),
  stringsAsFactors = F)

  # Annotate genes with GO CC terms and then filter
  # for those in lookup table
  annots <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = genes, columns = "GO",
                                keytype = "SYMBOL") %>%
    filter(ONTOLOGY == "CC") %>%
    dplyr::select(SYMBOL, GOID = GO) %>%
    filter(GOID %in% COMPARTMENTS_parent$GOID) %>%
    merge(COMPARTMENTS_parent, by = "GOID")

  return(annots)
}

