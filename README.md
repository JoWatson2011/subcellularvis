
# subcellularvis

<!-- badges: start -->
<!-- badges: end -->

The goal of subcellularvis is to simplify the interpretation of Gene Ontology Cellular Compartment enrichment analyses.
The app is hosted [here](phenome.manchester.ac.uk/subcellular/)

## Installation

You can install the released version of subcellularvis using [devtools](https://cran.r-project.org/web/packages/devtools/index.html) with:

``` r
devtools::install_github("jowatson2011/subcellularvis")
```

## Running the Shiny App from Rstudio

This example  shows you how to run a basic analysis.

Either through the Shiny app interface:

``` r
library(subcellularvis)

subcellularapp()
```

Or you can run analysis outside of the app:
``` r
genes <- c("MAPK1", "MAPK3")
comps <- compartmentData(genes)
runSubcellulaRvis(comps$enrichment, 
                  colScheme_low = "lightblue", 
                  colScheme_high = "darkblue")

# Use plotly for interactive visualisation
plotly::ggplotly(runSubcellulaRvis(comps$enrichment, 
                                   colScheme_low = "lightblue", 
                                   colScheme_high = "darkblue")
                                   )
```

