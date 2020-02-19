#' Generate dimensions of cell image to map to
#'
#' @param W Width of image in pixels
#' @param H Height of image in pixels
#' @return Dimensions of image to pass to ggplot2
#' @examples
#' vis_map(500, 500)

vis_map <- function(W, H){
  points <- data.frame(
    c(0, 0, W, W),
    c(0, H, H, 0)
  )
  colnames(points) <- c("x", "y")

  return(points)
}
