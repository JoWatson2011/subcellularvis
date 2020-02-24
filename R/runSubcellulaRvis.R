compartmentData <- function(genes, bkgdSize = 20367){
  library(GO.db)
  library(org.Hs.eg.db)
  library(dplyr)
  
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
  
  enrichment <- lapply(COMPARTMENTS_parent$GOID, function(x){
    annots <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = x, columns = "SYMBOL",
                                    keytype = "GO")

    successes.sample <- annots %>% 
      filter(SYMBOL %in% genes) %>% 
      nrow() # How many times term 'i' occurred in sample
    
    successes.bkgd <- nrow(annots) #How many times term 'i' occured in background
    
    failure <- bkgdSize - successes.bkgd  #Number of peptides in background minus the peptides with that term
    
    sampleSize <- length(genes)
    
    samplep <- phyper(successes.sample-1,
                            successes.bkgd,
                            failure,
                            sampleSize,
                            lower.tail = F)
      
      
  return(data.frame(GOID = x, p = samplep, 
                    stringsAsFactors = F))    
  }) %>% 
    do.call(rbind, .) %>% 
    merge(COMPARTMENTS_parent, by = "GOID")

  return(enrichment)
}
  
 
runSubcellulaRvis <- function(compartmentData, colScheme){
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(png)
  
  vis_map <- function(W, H){
    points <- data.frame(
      c(0, 0, W, W),
      c(0, H, H, 0)
    )
    colnames(points) <- c("x", "y")
    
    return(points)
  }
  
  vis_organelles <- function(compartment, W, H, X, Y, compartments_df = compartmentData){
    
    if(compartment %in% compartments_df$compartment){
      p <- compartments_df[compartments_df$compartment == compartment,]$p
    }else{
      p <- 1
    }
    points <- data.frame(
      c(X, X, X+W, X+W),
      c(Y, Y+H, Y+H, Y),
      rep(p, times = 4)
    )
    colnames(points) <- c(paste0(compartment, "_x"),
                          paste0(compartment, "_y"),
                          paste0(compartment, "_fill"))
    
    return(points)
  }
  
  df<- cbind(
    vis_map(1779.819, 1622.419),
    vis_organelles(compartment = "Cytoplasm", W = 1467, H = 1325, X = 153.5, Y = 1622.4 - 1467.5),
    vis_organelles("Nucleus", 431.7, 397, 867, 1622.4 - 1365.7),
    vis_organelles("Cytoskeleton", 52, 1255.5, 1510.5, 1622.4 - 1429.5),
    vis_organelles("Mitochondrion", 262.3, 132.445, 242.1, 1622.4 -  1164.8),
    vis_organelles("Endosome", 324.469, 493.7, 315.3, 1622.4 -  699.874),
    vis_organelles("Intracellular vesicle", 59.3, 59.3, 744.2, 1622.4 -  258.5),
    vis_organelles("Ribosome", 89.7, 109, 583.8, 1622.4 - 1321.5),
    vis_organelles("Plasma membrane", 1626.8, 1499, 72.4, 60),
    vis_organelles("Lysosome", 144.4, 140.9, 285.3, 1622.4 - 893.1),
    vis_organelles("Peroxisome", 143.5, 146.4, 572.9, 1622.4 - 921.8),
    vis_organelles("Extracellular region", 1779.8, 1622.4, 0, 0),
    vis_organelles("Endoplasmic reticulum", 532.2, 284.2, 814.6, 1622.4 - 637.6),
    vis_organelles("Golgi apparatus", 479.7, 208.4, 836.5, 1622.4 - 329.6)
  )
  
  img <- readPNG("data/CELL.png")
  
  g <- ggplot(df) +
    geom_polygon(aes(x = `Extracellular region_x`,
                     y = `Extracellular region_y`,
                     fill = `Extracellular region_fill`)) +
    geom_polygon(aes(x = `Plasma membrane_x`,
                     y = `Plasma membrane_y`,
                     fill = `Plasma membrane_fill`)) +
    geom_polygon(aes(x = Cytoplasm_x,
                     y = Cytoplasm_y,
                     fill = Cytoplasm_fill)) +
    geom_polygon(aes(x = Nucleus_x,
                     y = Nucleus_y,
                     fill = Nucleus_fill)) +
    geom_polygon(aes(x = Cytoskeleton_x,
                     y = Cytoskeleton_y,
                     fill = Cytoskeleton_fill)) +
    geom_polygon(aes(x = Mitochondrion_x,
                     y = Mitochondrion_y,
                     fill = Mitochondrion_fill)) +
    geom_polygon(aes(x = Endosome_x,
                     y = Endosome_y,
                     fill = Endosome_fill)) +
    geom_polygon(aes(x = `Intracellular vesicle_x`,
                     y = `Intracellular vesicle_y`,
                     fill = `Intracellular vesicle_fill`)) +
    geom_polygon(aes(x = Ribosome_x,
                     y = Ribosome_y,
                     fill = Ribosome_fill)) +
    geom_polygon(aes(x = Lysosome_x,
                     y = Lysosome_y,
                     fill = Lysosome_fill)) +
    geom_polygon(aes(x = Peroxisome_x,
                     y = Peroxisome_y,
                     fill = Peroxisome_fill)) +
    annotation_raster(img,
                      xmin=0, xmax=1777.411,
                      ymin=0, ymax= 1621.419) + #size of image
    scale_fill_gradientn(name = "Enrichment (p-value)",
                         colors = c("white", "grey", colScheme),
                         values = scales::rescale(c(1, 0.05, 0)),
                         guide = "colorbar") +
    theme(line = element_blank(),
          axis.text = element_blank(),
          title = element_blank(),
          panel.background = element_blank(), 
          legend.title = element_text(size = 6)
    ) +
    guides(fill = guide_colourbar(reverse = TRUE))
  
  return(g)
  
}
  