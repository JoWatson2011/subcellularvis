#' Calculate enrichment of subcellular compartment
#'
#' @param genes input genes (HGNC symbol) as character vector
#' @param bkgd background population. Defaults to size of species library
#' @param trafficking Logical; Calculate enrichment within endolysosome system?
#' @param subAnnots Logical; Calculate enrichment on all terms, not summarised
#' @param organism One of "Human", "Mouse",  "Drosophila", "Yeast", "Rat", "Xenopus"
#'
#'
#' @importFrom dplyr mutate filter select arrange bind_rows distinct group_by summarise
#' @importFrom stats phyper p.adjust
#' @importFrom rlang .data
#' @importFrom magrittr `%>%`
#'
#' @return data frame
#' @export


compartmentData <- function(genes, bkgd = NULL,
                            trafficking = F, subAnnots = F,
                            organism = c("Human", "Mouse", 
                                         "Drosophila", "Yeast", 
                                         "Rat", "Xenopus")){
  # Identify 'parent' GO terms + IDs.
  # http://www.supfam.org/SUPERFAMILY/cgi-bin/dcgo.cgi
  
  if(any(genes == "" | length(genes) < 1)){
    stop("No genes entered")
  }
  
  dat <- paste0(organism[1], 
                ifelse(trafficking, "_traffic_", "_"),
                ifelse(subAnnots, "subAnnots", "annots"))
  COMPARTMENTS_parent <- eval(as.name(dat))
  
  if(is.null(bkgd)){
    bkgd_size  <-  nrow(COMPARTMENTS_parent)
  } else {
    COMPARTMENTS_parent <- dplyr::filter(COMPARTMENTS_parent, .data$SYMBOL %in% bkgd)
    bkgd_size  <-  length(bkgd)
  }
  
  if(subAnnots == T) {
    COMPARTMENTS_parent$calc <- ifelse(COMPARTMENTS_parent$SYMBOL %in% genes, T, F)
  } else {
    COMPARTMENTS_parent$calc <- T
  }
  # Annotate genes with GO CC terms and then filter
  # for those in lookup table
  
  if(sum(COMPARTMENTS_parent$SYMBOL %in% genes) == 0){
    return(NULL)
  } else {
    enrichment <- lapply(
      unique(
        COMPARTMENTS_parent$compartment[COMPARTMENTS_parent$calc]),
      function(x){
        successes.sample <- COMPARTMENTS_parent %>% 
          dplyr::filter(.data$compartment == x) %>% 
          dplyr::filter(.data$SYMBOL %in% genes) %>% 
          unique() %>% 
          nrow() # How many times term 'i' occurred in sample
        
        successes.bkgd <- COMPARTMENTS_parent %>% 
          dplyr::filter(.data$compartment == x) %>% 
          unique() %>% 
          nrow()    #How many times term 'i' occurred in background
        
        failure <- bkgd_size - successes.bkgd  #Number of peptides in background minus the peptides with that term
        
        sampleSize <- length(genes)
        
        samplep <- stats::phyper(successes.sample,
                                 successes.bkgd,
                                 failure,
                                 sampleSize,
                                 lower.tail = F)
        
        symbols <- COMPARTMENTS_parent %>% 
          dplyr::filter(.data$compartment == x) %>% 
          dplyr::filter(.data$SYMBOL %in% genes) %>% 
          dplyr::select(.data$SYMBOL) %>% 
          unique()
        
        return(
          data.frame(
            Compartment = x,
            Symbol = paste(symbols$SYMBOL, collapse = ","),
            n = length(unique(symbols$SYMBOL)),
            p = samplep)
        )  
      }) %>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate(FDR = p.adjust(.data$p),
                    `FDR < 0.05` = ifelse(.data$FDR < 0.05, T, F)
      ) %>% 
      dplyr::arrange(.data$FDR) %>% 
      dplyr::select(.data$Compartment, .data$p,
                    .data$FDR, .data$`FDR < 0.05`,
                    .data$n, .data$Symbol)
    
    
    if(subAnnots){
      enrichment <- 
        merge(enrichment, 
              COMPARTMENTS_parent[,c("compartment", "group")], 
              by.x = "Compartment", by.y = "compartment",
              all.y = F) %>% 
        dplyr::distinct() %>% 
        dplyr::group_by(.data$Compartment, .data$FDR, .data$`FDR < 0.05`, .data$n, .data$Symbol) %>% 
        dplyr::summarise(group = paste(.data$group, collapse = ", "), 
                         .groups = "keep") %>% 
        dplyr::arrange(.data$FDR) %>% 
        dplyr::select(.data$Compartment,.data$group,
                      .data$FDR, .data$`FDR < 0.05`,
                      .data$n, .data$Symbol) %>% 
        as.data.frame()
    }
    
    return(enrichment)
  }
}


#' Visualise subcellular enrichment
#'
#' @param compsDat Output of compartmentData() 
#' @param colScheme_low Low value (i.e. most statistically ignificant) of colour scheme 
#' @param colScheme_high High value (i.e. most statistically significant) of colour scheme
#' @param trafficking Logical; Visualise endolysosome system? 
#' @param text_size Size of text if using 
#' @param legend Logical; include legend?
#' @param legend.pos passes to ggplot2::theme(legend.postion = legend.pos). One of "right", "left", "bottom", "top"
#' @param labels show labels on plot? Default is TRUE
#'
#' @importFrom ggplot2 ggplot aes geom_polygon annotation_raster scale_fill_gradientn theme element_blank element_text margin guides geom_label geom_text
#' @importFrom grid unit
#' @importFrom rlang .data
#' @importFrom magick image_read 
#' 
#' @return ggplot2 object
#' @export
runSubcellulaRvis <- function(compsDat, colScheme_low, 
                              colScheme_high, trafficking = F, 
                              text_size = 2,
                              legend = T,
                              legend.pos = "right",
                              labels = T){
  
  if(any(colnames(compsDat) != c("Compartment","p","FDR","FDR < 0.05","n","Symbol"))){
    stop("compsDat must be the output of compartmentData() function.")
  }
  
  vis_organelles <- function(compartment, W, H, X, Y, compartments_df = compsDat){
    
    if(!compartment %in% compsDat$Compartment){
      stop("compsDat must be the output of compartmentData() function.")
    }
    
    FDR <- compartments_df[compartments_df$Compartment == compartment,]$FDR
    
    if(length(FDR) == 0){
      FDR <-  1
    }
    
    points <- data.frame(
      c(X, X, X+W, X+W),
      c(Y, Y+H, Y+H, Y),
      rep(FDR, times = 4)
    )
    colnames(points) <- c(paste0(compartment, "_x"),
                          paste0(compartment, "_y"),
                          compartment)
    
    return(points)
  }
  
  
  if(trafficking){
    
    if("Golgi apparatus" %in% compsDat$Compartment){
      # Golgi selected randomly as an example not in trafficking subset
      stop("Trafficking variable = T, but compsDat not trafficking subset.\n
           Rerun compartmentsData with trafficking = T")
    }
    
    df<- cbind(
      vis_organelles(compartment = "Plasma membrane", W = 210.4, H = 15.663, X = 0, Y = 164.278),
      vis_organelles("Intracellular vesicle", 19.082, 18.552, 129.9, 135.2),
      vis_organelles("Lysosome", 39.5, 40.46, 149.186, 15.059),
      vis_organelles("Recycling Endosome", 67.8, 37.2, 122,  87.6),
      vis_organelles("Late Endosome", 92.709, 42.438, 34.619,  5.515),
      vis_organelles("Early Endosome", 75.373,  87.449, 17.575, 61.062)
    )
    
    img <- magick::image_read(system.file("extdata",
                                          "endosomalCell_scaled.png",
                                          package = "subcellularvis"))
    
    g <- ggplot2::ggplot(df
    ) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$`Plasma membrane_x`,
                                         y = .data$`Plasma membrane_y`,
                                         fill = .data$`Plasma membrane`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$`Intracellular vesicle_x`,
                                         y = .data$`Intracellular vesicle_y`,
                                         fill = .data$`Intracellular vesicle`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$`Lysosome_x`,
                                         y = .data$`Lysosome_y`,
                                         fill = .data$`Lysosome`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$`Recycling Endosome_x`,
                                         y = .data$`Recycling Endosome_y`,
                                         fill = .data$`Recycling Endosome`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$`Late Endosome_x`,
                                         y = .data$`Late Endosome_y`,
                                         fill = .data$`Late Endosome`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$`Early Endosome_x`,
                                         y = .data$`Early Endosome_y`,
                                         fill = .data$`Early Endosome`)) +
      ggplot2::annotation_raster(img,
                                 xmin=0, xmax=210.4,
                                 ymin=0, ymax= 189.85)   #size of image
    
  }else{
    if("Recycling Endosome" %in% compsDat$Compartment){
      # Recycling Endosome selected randomly as an example not in trafficking subset
      stop("Trafficking variable = F, but compsDat is the trafficking subset.\n
           Rerun compartmentsData with trafficking = F")
    }
    
    df<- cbind(
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
    
    img <- magick::image_read(system.file("extdata",
                                          "CELL_scaled.png",
                                          package = "subcellularvis")
                              )
    
    g <-  ggplot2::ggplot(df) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$`Extracellular region_x`,
                                         y = .data$`Extracellular region_y`,
                                         fill = .data$`Extracellular region`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$`Plasma membrane_x`,
                                         y = .data$`Plasma membrane_y`,
                                         fill = .data$`Plasma membrane`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$Cytoplasm_x,
                                         y = .data$Cytoplasm_y,
                                         fill = .data$Cytoplasm)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$Nucleus_x,
                                         y = .data$Nucleus_y,
                                         fill = .data$Nucleus)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$Cytoskeleton_x,
                                         y = .data$Cytoskeleton_y,
                                         fill = .data$Cytoskeleton)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$Mitochondrion_x,
                                         y = .data$Mitochondrion_y,
                                         fill = .data$Mitochondrion)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$Endosome_x,
                                         y = .data$Endosome_y,
                                         fill = .data$Endosome)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$`Intracellular vesicle_x`,
                                         y = .data$`Intracellular vesicle_y`,
                                         fill = .data$`Intracellular vesicle`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$Ribosome_x,
                                         y = .data$Ribosome_y,
                                         fill = .data$Ribosome)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$Lysosome_x,
                                         y = .data$Lysosome_y,
                                         fill = .data$Lysosome)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$Peroxisome_x,
                                         y = .data$Peroxisome_y,
                                         fill = .data$Peroxisome)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$`Endoplasmic reticulum_x`,
                                         y = .data$`Endoplasmic reticulum_y`,
                                         fill = .data$`Endoplasmic reticulum`)) +
      ggplot2::geom_polygon(ggplot2::aes(x = .data$`Golgi apparatus_x`,
                                         y = .data$`Golgi apparatus_y`,
                                         fill = .data$`Golgi apparatus`))  +
      ggplot2::annotation_raster(img,
                                 xmin=0, xmax=1777.411,
                                 ymin=0, ymax= 1621.419)   #size of image
  }
  
  g <- g +
    ggplot2::scale_fill_gradientn(name = "FDR",
                                  breaks = seq(0, max(compsDat$FDR[compsDat$FDR <= 0.05]), max(compsDat$FDR[compsDat$FDR <= 0.05])/5), 
                                  limits = c(0, max(compsDat$FDR[compsDat$FDR <= 0.05])),
                                  colors = c(colScheme_low, colScheme_high),
                                  guide = "colorbar",
                                  na.value = "white")+
    ggplot2::theme(line = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.ticks.length = grid::unit(0, "pt"),
                   axis.line = ggplot2::element_blank(),
                   title = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(), 
                   panel.border= ggplot2::element_blank(), 
                   panel.spacing = grid::unit(0, "cm"),
                   panel.grid.major= ggplot2::element_blank(),
                   panel.grid.minor= ggplot2::element_blank(),
                   legend.title = ggplot2::element_text(size = text_size+5),
                   legend.text = ggplot2::element_text(size = text_size+5,
                                                       angle = ifelse(
                                                         legend.pos %in% 
                                                           c("right", "left"),
                                                         0,
                                                         45),
                                              hjust = ifelse(legend.pos %in% c("right", "left"),
                                                             0,
                                                             1)),
                   # legend.key.height = grid::unit(0.1, "cm"),
                   # legend.key.width = grid::unit(ifelse(legend.pos %in% c("right", "left"),
                   #                                0.1,
                   #                                0.5), 
                   #                         "cm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0, "mm"),
                   legend.box.background = ggplot2::element_blank(),
                   legend.box.spacing = ggplot2::margin(0, 0, 0, 0, "mm"),
                   legend.position = legend.pos,
                   legend.direction = ifelse(legend.pos %in% c("right", "left"), "vertical", "horizontal"),
                   plot.margin = ggplot2::margin(0, 0, 0, 0, "mm"),
                   plot.background = ggplot2::element_blank()
    ) 
  
  
  
  if(legend == F){
    g <- g +
      ggplot2::guides(fill = "none")
  }
  
  if(labels == T){
    if(trafficking){
      labels_df <-  data.frame(
        label = c("Early Endosome", "Late Endosome",
                  "Vesicle", "Recycling Endosome",
                  "Lysosome", "Plasma Membrane"),
        x = c(55, 80, 140, 155, 170, 105),
        y = c(110, 25, 145, 105, 35, 170),
        angle = c(rep(0, 6))
      )
    }else{
      labels_df <-  data.frame(
        label = c("Vesicle", "Endosomal system", "Lysosome",
                  "Peroxisome", "Mitochondria", "Ribosome",
                  "Golgi apparatus", "Endoplasmic reticulum", 
                  "Nucleus", "Extracellular region", "Plasma membrane", 
                  "Cytosol", "Cytoskeleton"),
        x = c(760, 500, 350, 650, 375, 645, 1075, 1075, 1075, 200, 320, 450, 1470),
        y = c(1400, 1200, 800, 760, 520, 375, 1200, 865, 450, 40, 100, 168, 1200),
        rect = c(F, T, T, T, T, F, T, T, T, F, F, F, F),
        angle = c(rep(0, 12), 90)
      )
    }
    
    g <- g + 
      ggplot2::geom_label(data = labels_df[labels_df$rect == T,],
                          ggplot2::aes(x = .data$x, 
                                       y = .data$y, 
                                       label = .data$label),
                          alpha = 0.5, color = NA,
                          size = text_size) +
      ggplot2::geom_text(data = labels_df, 
                         ggplot2::aes(x = .data$x,
                                      y = .data$y,
                                      label = .data$label,
                                      angle = .data$angle),
                         fontface = "bold",
                         size = text_size
      )
  }
  
  return(g)
  
}