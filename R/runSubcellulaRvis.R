#' Calculate enrichment of subcellular compartment
#'
#' @param genes input genes (HGNC symbol) as character vector
#' @param bkgdSize size of background proteome. Default = 20367, HeLa proteome
#' @param trafficking Logical; Calculate enrichment within endolysosome system?
#'
#' @importFrom magrittr `%>%`
#' @importFrom dplyr filter
#' @importFrom phyper stats
#'
#' @return
#' @export
#'
#' @examples
#' 
compartmentData <- function(genes, bkgdSize = 20367, trafficking = F){
  # Identify 'parent' GO terms + IDs.
  # http://www.supfam.org/SUPERFAMILY/cgi-bin/dcgo.cgi
  if(trafficking){
    COMPARTMENTS_parent <- readRDS("data/traffic_annots.rds")
  }else{
    COMPARTMENTS_parent <- readRDS("data/annots.rds")
  }
  
  # Annotate genes with GO CC terms and then filter
  # for those in lookup table
  
  enrichment <- lapply(unique(COMPARTMENTS_parent$compartment), function(x){
    successes.sample <- COMPARTMENTS_parent %>% 
      filter(compartment == x) %>% 
      filter(SYMBOL %in% genes) %>% 
      unique() %>% 
      nrow() # How many times term 'i' occurred in sample
    
    successes.bkgd <- COMPARTMENTS_parent %>% 
      filter(compartment == x) %>% 
      unique() %>% 
      nrow()    #How many times term 'i' occured in background
    
    failure <- bkgdSize - successes.bkgd  #Number of peptides in background minus the peptides with that term
    
    sampleSize <- length(genes)
    
    samplep <- phyper(successes.sample-1,
                      successes.bkgd,
                      failure,
                      sampleSize,
                      lower.tail = F)
    
    return(data.frame(compartment = x, p = samplep))    
  }) %>% 
    do.call(rbind, .) 
  
  return(enrichment)
}


#' Visualise subcellular enrichment
#'
#' @param compartmentData Output of compartmentData() 
#' @param colScheme_low Low value (i.e. most statistically ignificant) of colour scheme 
#' @param colScheme_high High value (i.e. most statistically significant) of colour scheme
#' @param trafficking Logical; Visualise endolysosome system? 
#' @param text_size Size of text if using 
#' @param legend Logical; include legend?
#' @param legend.pos passes to ggplot::theme(legend.postion = legend.pos). One of "right", "left", "bottom", "top"

#' @importFrom ggplot2 ggplot aes geom_polygon annotation_raster scale_fill_gradientn theme element_blank element_text unit margin guides guide_colourbar
#' @importFrom png readPNG

#' @return
#' @export
#'
#' @examples
runSubcellulaRvis <- function(compartmentData, colScheme_low, 
                              colScheme_high, trafficking = F, 
                              text_size = 6,
                              legend = T,
                              legend.pos = "right"){

  # library(ggplot2) 
  # library(png)
  
  vis_map <- function(W, H){
    points <- data.frame(
      c(0, 0, W, W),
      c(0, H, H, 0)
    )
    colnames(points) <- c("x", "y")
    
    return(points)
  }
  
  vis_organelles <- function(compartment, W, H, X, Y, compartments_df = compartmentData){
    
  #  if(compartment %in% compartments_df$compartment){
      p <- compartments_df[compartments_df$compartment == compartment,]$p
 #   }else{
 #     p <- 1
 #   }
    points <- data.frame(
      c(X, X, X+W, X+W),
      c(Y, Y+H, Y+H, Y),
      rep(p, times = 4)
    )
    colnames(points) <- c(paste0(compartment, "_x"),
                          paste0(compartment, "_y"),
                          compartment)
    
    return(points)
  }
  
  if(trafficking){
    df<- cbind(
      vis_map(210.4, 189.85),
      vis_organelles(compartment = "Plasma membrane", W = 210.4, H = 15.663, X = 0, Y = 164.278),
      vis_organelles("Intracellular vesicle", 19.082, 18.552, 129.9, 135.2),
      vis_organelles("Lysosome", 39.5, 40.46, 149.186, 15.059),
      vis_organelles("Recycling Endosome", 67.8, 37.2, 122,  87.6),
      vis_organelles("Late Endosome", 92.709, 42.438, 34.619,  5.515),
      vis_organelles("Early Endosome", 75.373,  87.449, 17.575, 61.062)
    )
    
    # compartmentLabels <- data.frame(
    #   label = c("Plasma membrane", "Intracellular vesicle", "Lysosome",
    #             "Recycling Endosome", "Late Endosome", "Early Endosome")
    #   x = 
    # )
    img <- png::readPNG("data/endosomalCell.png")
    
    g <- ggplot2::ggplot(df) +
      ggplot2::geom_polygon(ggplot2::aes(x = `Plasma membrane_x`,
                       y = `Plasma membrane_y`,
                       fill = `Plasma membrane`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = `Intracellular vesicle_x`,
                       y = `Intracellular vesicle_y`,
                       fill = `Intracellular vesicle`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = `Lysosome_x`,
                       y = `Lysosome_y`,
                       fill = `Lysosome`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = `Recycling Endosome_x`,
                       y = `Recycling Endosome_y`,
                       fill = `Recycling Endosome`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = `Late Endosome_x`,
                       y = `Late Endosome_y`,
                       fill = `Late Endosome`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = `Early Endosome_x`,
                       y = `Early Endosome_y`,
                       fill = `Early Endosome`)) +
      ggplot2::annotation_raster(img,
                        xmin=0, xmax=210.4,
                        ymin=0, ymax= 189.85)   #size of image
    
  }else{
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
      vis_organelles("Endoplasmic reticulum", 532.2, 284.2, 814.6, 696.72),
      vis_organelles("Golgi apparatus", 479.7, 208.4, 836.5, 1083.75)
    ) 
    
    
    img <- pngreadPNG("data/CELL.png")
    
    g <- ggplot2::ggplot(df) +
      ggplot2::geom_polygon(ggplot2::aes(x = `Extracellular region_x`,
                       y = `Extracellular region_y`,
                       fill = `Extracellular region`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = `Plasma membrane_x`,
                       y = `Plasma membrane_y`,
                       fill = `Plasma membrane`)) +
      ggplot2::eom_polygon(ggplot2::aes(x = Cytoplasm_x,
                       y = Cytoplasm_y,
                       fill = Cytoplasm)) +
      ggplot2::geom_polygon(ggplot2::aes(x = Nucleus_x,
                       y = Nucleus_y,
                       fill = Nucleus)) +
      ggplot2::geom_polygon(ggplot2::aes(x = Cytoskeleton_x,
                       y = Cytoskeleton_y,
                       fill = Cytoskeleton)) +
      ggplot2::geom_polygon(ggplot2::aes(x = Mitochondrion_x,
                       y = Mitochondrion_y,
                       fill = Mitochondrion)) +
      ggplot2::geom_polygon(ggplot2::aes(x = Endosome_x,
                       y = Endosome_y,
                       fill = Endosome)) +
      ggplot2::geom_polygon(ggplot2::aes(x = `Intracellular vesicle_x`,
                       y = `Intracellular vesicle_y`,
                       fill = `Intracellular vesicle`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = Ribosome_x,
                       y = Ribosome_y,
                       fill = Ribosome)) +
      ggplot2::geom_polygon(ggplot2::aes(x = Lysosome_x,
                       y = Lysosome_y,
                       fill = Lysosome)) +
      ggplot2::geom_polygon(ggplot2::aes(x = Peroxisome_x,
                       y = Peroxisome_y,
                       fill = Peroxisome)) +
      ggplot2::geom_polygon(ggplot2::aes(x = `Endoplasmic reticulum_x`,
                       y = `Endoplasmic reticulum_y`,
                       fill = `Endoplasmic reticulum`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = `Golgi apparatus_x`,
                       y = `Golgi apparatus_y`,
                       fill = `Golgi apparatus`)) +
      ggplot2::annotation_raster(img,
                        xmin=0, xmax=1777.411,
                        ymin=0, ymax= 1621.419)   #size of image
  }
  
    g <- g +
      ggplot2::scale_fill_gradientn(name = "FDR",
                           breaks = seq(0, 0.05, 0.01), 
                           limits = c(0, 0.05),
                           colors = c(colScheme_low, colScheme_high),
                           guide = "colorbar",
                           na.value = "white")+
      ggplot2::theme(line = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.ticks.length = unit(0, "pt"),
            axis.line = element_blank(),
            title = element_blank(),
            panel.background = element_blank(), 
            panel.border=element_blank(), 
            panel.spacing = unit(0, "cm"),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            legend.title = element_text(size = text_size),
            legend.text = element_text(size = text_size,
                                       angle = ifelse(legend.pos %in% c("right", "left"),
                                                      0,
                                                      45),
                                       hjust = ifelse(legend.pos %in% c("right", "left"),
                                                      0,
                                                      1)),
            legend.key.height = unit(0.1, "cm"),
            legend.key.width = unit(ifelse(legend.pos %in% c("right", "left"),
                                           0.1,
                                           0.5), 
                                    "cm"),
            legend.margin = margin(0, 0, 0, 0, "mm"),
            legend.box.background = element_blank(),
            legend.box.spacing = margin(0, 0, 0, 0, "mm"),
            legend.position = legend.pos,
            legend.direction = ifelse(legend.pos %in% c("right", "left"), "vertical", "horizontal"),
            plot.margin = margin(0, 0, 0, 0, "mm"),
            plot.background = element_blank()
      )+
      guides(fill = guide_colourbar(reverse = TRUE))
    
    if(legend == F){
      g <- g +
        guides(fill = "none")
    }
      
  return(g)
  
}

