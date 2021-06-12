cran <- c(
  "shiny",
  "shinythemes",
  "rmarkdown",
  "dplyr",
  "magick",
  "colourpicker"
)

lapply(cran, function(pkg){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
})