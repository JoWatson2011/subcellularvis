library(shiny)
library(ggplot2) # ggsave
library(stringi) # stri_split
library(plotly)

source("R/runSubcellulaRvis.R")

### DEFINE UI FOR APP ###

ui <- fluidPage(
  theme = "cosmo",
  #Title
  titlePanel("Visualise Subcellular Localisation of 'Omics Data"),
  
  sidebarLayout(
    sidebarPanel(textAreaInput(inputId = "input_genes_text",
                               label = "List of genes,\none per line"),
                 
                 actionLink(inputId = "exampleGenes",
                            label = "Example gene set"),
                 
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
                               label = "Interested in trafficking?"),
                 br(),
                 actionButton(inputId = "action_button",
                              label = "Calculate Enrichment & Visualise")
    ),
    
    mainPanel(
      
      tabsetPanel(
        tabPanel("Table", 
                 tableOutput(outputId = "compartment_df")),
        
        tabPanel("Plot", 
                 fluidRow(
                   column(width = 12,
                          plotlyOutput(outputId = "plot_cell")
                   )
                 ),
                 
                 fluidRow(
                   column(width = 6,
                          p(strong("Colour scale:")),
                          textInput(inputId = "colScheme_low", 
                                    label = "Colour scale low",
                                    value = "#800000"),
                          textInput(inputId = "colScheme_high", 
                                    label = "Colour scale high",
                                    value = "#d8b2b2"),
                          actionButton(inputId = "remakePlot",
                              label = "Remake Plot")
                   ),
                 column(width = 6,  
                        p(strong("Options for export")),
                        checkboxInput("includeLegend",
                                      "Include legend?",
                                      T),
                        selectInput("legendPos",
                                    "Legend positon",
                                    list("bottom", "right", "left", "top")),
                        numericInput("plotWidth", "Plot width, cm", 15),
                        numericInput("plotHeight", "Plot height, cm", 10),
                        numericInput("textSize", "Text size", 6),
                        selectInput("export_fileType", "Image file type",
                                    choices = c("pdf", "tiff", "png", "jpeg")),
                        downloadButton(outputId = "downloadPlot", 
                                       label = "Download")
                        )
                 ),
                 
                 fluidRow(
                   column(width = 12,
                          p(strong("Export log")), 
                          p("This can be used in manuscript methods session and aid reproducibility of your analysis"),
                          downloadButton(outputId = "downloadLog", 
                                         label = "Export log.txt")
                          )
                   
                 )
        ),
        
        tabPanel(
          
        )
      )
    )
  )
)

### DEFINE SERVER LOGIC ###

server <- function(input, output){
  
  observeEvent(input$exampleGenes, {
  updateTextAreaInput(inputId = "input_genes_text", 
                      value = "MAPK1\nMAPK3")
  })
  
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
      
      plot_cell(runSubcellulaRvis(comps,
                                  colScheme_low = input$colScheme_low,
                                  colScheme_high = input$colScheme_high, 
                                  trafficking = input$traffic,
                                  legend = input$includeLegend))
      
    })
    
    output$compartment_df <- renderTable({
      comps
    })
    
    output$plot_cell <- renderPlotly({
      print(ggplotly(plot_cell()))
    })
  })
  
  observeEvent(input$remakePlot, {
    withProgress(message = "rendering", value = 15, {
      plot_cell(runSubcellulaRvis(comps, 
                                  colScheme_low = input$colScheme_low,
                                  colScheme_high = input$colScheme_high, 
                                  text_size = input$textSize,
                                  trafficking = input$traffic,
                                  legend = input$includeLegend,
                                  legend.pos = input$legendPos))
    })
  })
  
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("CellVis_", Sys.Date(), input$export_fileType, sep="")
    },
    content = function(file) {
      ggplot2::ggsave(file, 
             plot = plot_cell(),
             device = input$export_fileType,
             width = input$plotWidth,
             height = input$plotHeight,
             dpi = 300,
             units = "cm")
      },
    contentType = ifelse(input$export_fileType == "pdf", 
                         "application/pdf",
                         paste0("image/", input$export_fileType))
  )
  
  output$downloadLog <- downloadHandler(
    filename = function() {
      paste("CellVis_log_", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
        cat(
          sprintf(
            "Analysis performed on %s using SubCellularVis tool\n(https://github.com/JoWatson2011/shiny_subcellularvis)\n",
            Sys.Date()
          ),
          file = file
        )
      cat("\nEnrichment FDR:\n", file = file, append = T)
      apply(comps, 1, 
            function(i){
              cat(c(paste(i, collapse = "\t"), "\n"),
                  file = file, append = T)
            })
      cat("\nGene list:\n", file = file, append = T)
      cat(genes_tidy, file = file, sep = "\n", append = T)
    },
    contentType = "text/plain"
  )

}


### RUN THE APP ###
shinyApp(ui = ui, server = server)
