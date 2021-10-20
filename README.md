
# subcellularvis

<!-- badges: start -->
<!-- badges: end -->

The goal of subcellularvis is to ...
The app is hosted [here]()

## Installation

You can install the released version of subcellularvis using [devtools](https://cran.r-project.org/web/packages/devtools/index.html) with:

``` r
devtools::install_github("jowatson2011/subcellularvis")
```

## Running the Shiny App from Rstudio

This is a basic example which shows you how to solve a common problem:

``` r
library(subcellularvis)

subcellularapp()
```

To run a basic analysis outside of the app:
``` r
genes <- c("MAPK1", "MAPK3")
comps <- compartmentData(genes)
runSubcellulaRvis(comps, 
                  colScheme_low = "lightblue", 
                  colScheme_high = "darkblue")

# Use plotly for interactive visualisation
plotly::ggplotly(runSubcellulaRvis(comps, 
                                   colScheme_low = "lightblue", 
                                   colScheme_high = "darkblue")
                                   )
```

