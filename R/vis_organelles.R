#' Generate dimensions of organelles in cell image to map to.
#' @param compartment The name of the organelle (i.e. the GO CC term) to be mapped to.
#' @param W The width of the organelle in the cell image, in pixels.
#' @param H The height of the organelle in the cell image, in pixels.
#' @param X The bottom left X coordinate of the organelle in the cell image, in pixels.
#' @param Y The bottom left Y coordinate of the organelle in the cell image, in pixels.
#' @param compartments_df The output of compartmentsData.
#' @examples
#' organelle_map("Cytoplasm", 1467, 1325, 153.5, 154.9)

vis_organelles <- function(compartment, W, H, X, Y, compartments_df = compartmentAnnots){

  freq <- compartments_df[compartments_df$compartment == compartment,]$freq

  points <- data.frame(
    c(X, X, X+W, X+W),
    c(Y, Y+H, Y+H, Y),
    rep(freq, times = 4)
  )
  colnames(points) <- c(paste0(compartment, "_x"),
                        paste0(compartment, "_y"),
                        paste0(compartment, "_fill"))

  return(points)
}
