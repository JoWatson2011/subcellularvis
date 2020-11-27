library(shiny)
library(visNetwork)
library(igraph)
library(dplyr)
library(tidyr)
simplifyData1 <<- readRDS("data/simplifyData1.rds")
STRING.expmt.gene <<- readRDS("data/STRINGexpmtgene.rds")
GOBP <<- list(
  bg = readRDS("data/GOBP_background.rds"),
  distr = readRDS("data/GOBPdistr.rds")
)
GOCC <<- list(
  bg = readRDS("data/GOCC_background.rds"),
  distr = readRDS("data/GOCCdistr.rds")
)
GOMF <<- list(
  bg = readRDS("data/GOMF_background.rds"),
  distr = readRDS("data/GOMFdistr.rds")
)
REACT <<-list(
  bg = readRDS("data/REACTOME_background.rds"),
  distr = readRDS("data/REACTdistr.rds")
)
databases <- list(
  GOBP = GOBP,
  GOCC = GOCC,
  GOMF = GOMF,
  REACT = REACT
)
source("R/runSubcellulaRvis.R")

### DEFINE UI FOR APP ###

ui <- fluidPage(

  #Title
  titlePanel("Visualise Subcellular Localisation of 'Omics Data"),

  sidebarLayout(
    sidebarPanel(textAreaInput(inputId = "input_genes_text",
                               label = "List of genes, \none per line"),
                 
                 # actionButton(inputId = "action_text",
                 #              label = "From textbox"),
                 # br(),br(),
                 fileInput(inputId = "input_genes_file",
                           label = "Select file",
                           multiple = F,
                           placeholder = ".csv file of genes"),
                 # actionButton(inputId = "action_file",
                 #              label = "From .csv file")
                 checkboxInput(inputId = "traffic",
                               label = "Intersted in trafficking?"),
                 br(),
                 actionButton(inputId = "action_button",
                              label = "Calculate Enrichment & Visualise")
    ),

    mainPanel(
      
      tabsetPanel(
        tabPanel("Table", 
                 tableOutput(outputId = "compartment_df")),
        
        tabPanel("Plot", 
                 plotOutput(outputId = "plot_cell"),
                 selectInput(inputId = "colScheme", 
                             label = "Plot colour scale",
                             choices = list("Red" = "#800000", 
                                            "Green" = "#065535",
                                            "Blue" = "#003366",
                                            "Purple" = "#8968cd",
                                            "Pink" = "#c71585",
                                            "Orange" = "#cd3700"),
                             selected = 1),
                 actionButton(inputId = "remakePlot",
                              label = "Remake Plot"),
                 downloadButton(outputId = "downloadPlot", 
                                label = "Download"))
      ),

      p("Text placeholder")
    )
  )
)

### DEFINE SERVER LOGIC ###

server <- function(input, output){
  
  plot_cell <- reactiveVal()
  
  observeEvent(input$action_button, {

    file <- input$input_genes_file
    
    if(!is.null(file)){
      genes_tidy <<- read.csv(file$datapath,
                             header = F,
                             stringsAsFactors = F)[,1]
    }else{
      genes_tidy <<- as.vector(stringi::stri_split(input$input_genes_text,
                                                  regex = "\n",
                                                  simplify=T))
    }
    
    withProgress(message = "calculating", value = 15, {
      comps <<- compartmentData(genes = genes_tidy, trafficking = input$traffic)
      
      plot_cell(runSubcellulaRvis(comps, colScheme = input$colScheme, trafficking = input$traffic))
      
    })
    
    output$compartment_df <- renderTable({
      comps
    })
    
    output$plot_cell <- renderPlot({
     plot_cell()
    })
  })
  
  observeEvent(input$remakePlot, {
    withProgress(message = "rendering", value = 15, {
      plot_cell(runSubcellulaRvis(comps, colScheme = input$colScheme))
    })
  })
  
  
  output$downloadPlot <- downloadHandler(
    filename = 'test.png',
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = width, height = height,
                       res = 300, units = "in")
      }
      ggsave(file, device = device, width = 15, height = 10, units = "cm")
    }
  )
  
}


### RUN THE APP ###
shinyApp(ui = ui, server = server)
