library(shiny)
library(shinythemes)
library(ggplot2) # ggsave
library(stringi) # stri_split
library(plotly)
library(shinyFeedback)

source("R/runSubcellulaRvis.R")

### DEFINE UI FOR APP ###

ui <- fluidPage(
  theme = shinytheme("cosmo"),
  shinyFeedback::useShinyFeedback(),
  #Title
  titlePanel("SubcellulaRVis: Visualising enrichment of protein localisation"),
  
  sidebarLayout(
    sidebarPanel(textAreaInput(inputId = "input_genes_text",
                               label = "List of genes,\none per line"),
                 
                 actionLink(inputId = "exampleGenes",
                            label = "Example gene set"),
                 br(),
                 # actionButton(inputId = "action_text",
                 #              label = "From textbox"),
                 # br(),br(),
                 fileInput(inputId = "input_genes_file",
                           label = "Select file",
                           multiple = F,
                           placeholder = ".csv file of genes"),
                 # actionButton(inputId = "action_file",
                 #              label = "From .csv file")
                 
                 numericInput("bkgdSize",
                              "Background size",
                              20367),
                 
                 actionLink(inputId = "explainBkgd",
                            label = "What is the background size?",
                            icon = icon("question")),
                 
                 
                 checkboxInput(inputId = "traffic",
                               label = "Interested in trafficking?"),
                 br(),
                 actionButton(inputId = "action_button",
                              label = "Calculate Enrichment & Visualise"),
                 textOutput("checkInput")
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabsPanel",
        tabPanel(
          "About",
          
          p(strong("SubcellulaRVis a tool for visualising enrichment of\n
                   Gene Ontology Cellular Compartments within gene lists")),
          
         # p("Source code can be found at")
          tagList("Source code can be found at:", 
                  a("www.github.com/jowatson2011/shiny_subcellularvis",
                    href = "www.github.com/jowatson2011/shiny_subcellularvis")
                  ),
         br(),br(),
         p("Supported by:"),
         img(src="https://bbsrc.ukri.org/bbsrc/includes/themes/BBSRC/images/logo-1.png",
             width = "200px"),
         img(src="https://www.staffnet.manchester.ac.uk/brand/visual-identity/logo/logo_big.gif",
             width = "200px",)
         
          
          #       icon = "dna"
        ),
        
        tabPanel("Plot", 
                 fluidRow(
                   column(width = 12,
                          plotlyOutput(outputId = "plot_cell",
                                       width = "600px")
                   )
                 ),
                 
                 fluidRow(
                   column(width = 6,
                          p(strong("Colour scale:")),
                          textInput(inputId = "colScheme_low", 
                                    label = "Colour scale low",
                                    value = "#567ace"),
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
                                      list("bottom", "right", "left", "top"),
                                      "right"),
                          checkboxInput("includeLabels",
                                        "Label compartments?",
                                        T),
                          numericInput("plotWidth", "Plot width, cm", 15),
                          numericInput("plotHeight", "Plot height, cm", 10),
                          numericInput("textSize", "Text size", 6),
                          selectInput("export_fileType", "Image file type",
                                      choices = list("pdf", "tiff", "png", "jpeg")),
                          downloadButton(outputId = "downloadPlot", 
                                         label = "Download")
                   )
                 ),
                 
                 #              icon = "chart-bar"
        ),
        tabPanel("Table", 
                 fluidRow(width = 6,
                          tableOutput(outputId = "compartment_df")
                 ),
                 fluidRow(
                   column(width = 6,
                          p(strong("Export log")), 
                          p("This can be used in manuscript methods session and aid reproducibility of your analysis"),
                          downloadButton(outputId = "downloadLog", 
                                         label = "Export log.doc")
                   )
                   
                 ),
                 #             icon = "table"
        ),
        tabPanel(
          "Full enrichment",
          br(),
          
          fluidRow(
            column(width = 6,
                   tableOutput(outputId = "fullEnrichment_df")
            ),
            column(width = 6,
                   br(),
                   p("Top 50 terms shown"),
                   downloadButton(outputId = "downloadFullEnrichment", 
                                  label = "Export full results")
            ),
          ),
          #     icon = "table"
        ),
        tabPanel(
          "Help",
          br(),
          p(strong("How does SubCellulaRVis work?")),
          p("SubCellulaRVis calculates the enrichment of . The broad categories of organelles and compartments were drawn from the Gene Ontology Cellular Compartments (GOCC). Using the GO hierachy, all \'child\' terms of these categories are grouped together to calculate a single enrichment score for each category. We can then visualise this in a single graphic."),
          p(strong("What is the background size?")),
          p("The \"background size\" is the size of population being sampled (e.g. proteome size for a cell type or species). The default background size is the predicted number of proteins in the HeLa proteome."),
          p(strong("What is the difference between the simple and full enrichment?")),
          p("The simple enrichment is based on the categorisation as described above. The full enrichment is a standard enrichment analysis based on all the GOCC terms."),
          br(),
          p(strong("How can I export my results?")),
          p("The visualisation can be exported in multiple image formats from the \"Plot\" tab.\n
            Enrichment results can be exported as a .csv file from the \"Full Enrichment\" tab.\n
            A log of the analysis can be exported as a .doc file from the \"Table\" tab"),
          br(),
          br(),
          p("If you could not find your answer here, please contact: joanne.watson@manchester.ac.uk\n
            or raise an issue on the GitHub repository:" ),
          tagList(a("www.github.com/jowatson2011/shiny_subcellularvis",
                    href = "www.github.com/jowatson2011/shiny_subcellularvis")
          )
         
          #    icon = "question-circle"
        )
      )
    )
  )
)

### DEFINE SERVER LOGIC ###

server <- function(input, output){
  
  observeEvent(input$exampleGenes, {
    updateTextAreaInput(inputId = "input_genes_text", 
                        value = paste0(read.csv("data/exampleGenes.csv",
                                                header = F,
                                                stringsAsFactors = F)[,1], collapse = "\n"))
  })
  
  observeEvent(input$explainBkgd, {
    updateTabsetPanel(inputId = "tabsPanel", selected = "Help")
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
    output$checkInput <- renderText({  #genesNull())
      if(!isTruthy(genes_tidy)){
        validate("HGNC gene symbols need to be entered,\n one per line\nas .csv or in the text box.")
      }

    })

    output$downloadLog <- downloadHandler(
      filename = function() {
        paste("CellVis_log_", Sys.Date(), ".doc", sep="")
      },
      content = function(file) {
        tempReport <- file.path(tempdir(), "log.Rmd")
        file.copy("log.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(
          date = as.character(Sys.Date()),
          comps = comps,
          genes = genes_tidy
        )
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )}
    )
    if(isTruthy(genes_tidy)){
      updateTabsetPanel(inputId = "tabsPanel", selected = "Table")
    }

    output$compartment_df <- renderTable({
      withProgress(message = "Calculating", value = 15, {
        req(genes_tidy)
        comps <- compartmentData(genes = genes_tidy,
                                 bkgdSize = input$bkgdSize,
                                 trafficking = input$traffic)
        
        comps
      })
    }, 
    digits = 5)
    
    
    output$plot_cell <- renderPlotly({
      withProgress(message = "Visualising", value = 15, {
        req(genes_tidy,cancelOutput = T)
        comps <- compartmentData(genes = genes_tidy,
                                 bkgdSize = input$bkgdSize,
                                 trafficking = input$traffic)
        plot_cell(
          runSubcellulaRvis(comps,
                            colScheme_low = input$colScheme_low,
                            colScheme_high = input$colScheme_high,
                            trafficking = input$traffic,
                            text_size = input$textSize,
                            legend = input$includeLegend,
                            legend.pos = input$legendPos,
                            labels = F)
        )
        ggplotly(plot_cell())
      })
    })
  })
  
  
  observeEvent(input$remakePlot, {
    withProgress(message = "Revisualisng", value = 15, {
      plot_cell(
        runSubcellulaRvis(comps,
                          colScheme_low = input$colScheme_low,
                          colScheme_high = input$colScheme_high,
                          text_size = input$textSize,
                          trafficking = input$traffic,
                          legend = input$includeLegend,
                          legend.pos = input$legendPos,
                          labels = F)
      )
    })
  })
  
  
  # Export  visualisation
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
  
  observeEvent(input$export_fileType,{
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste("CellVis_", Sys.Date(),".", input$export_fileType, sep="")
      },
      content = function(file) {
        ggplot2::ggsave(file, 
                        plot = runSubcellulaRvis(comps,
                                                 colScheme_low = input$colScheme_low,
                                                 colScheme_high = input$colScheme_high,
                                                 text_size = input$textSize,
                                                 trafficking = input$traffic,
                                                 legend = input$includeLegend,
                                                 legend.pos = input$legendPos,
                                                 labels = input$includeLabels),
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
  })
  
  #####
  # Calculate "full" enrichment
  #####
  
  
  output$fullEnrichment_df <- renderTable({
    req(genes_tidy)
    withProgress(message = "Calculating", value = 15, {
      comps_full <<- compartmentData(genes = genes_tidy,
                                     bkgdSize = input$bkgdSize,
                                     trafficking = input$traffic,
                                     subAnnots = T)[1:50,]
      comps_full
    }) 
  },
  digits = 5)
  
  
  output$downloadFullEnrichment <- downloadHandler(
    filename = function() {
      paste("CellVis_fullEnrichment_", Sys.Date(), ".csv", sep="")
    },
    content = function(file){
      write.csv(comps_full, file)
    })
  
  
}


### RUN THE APP ###
shinyApp(ui = ui, server = server)
