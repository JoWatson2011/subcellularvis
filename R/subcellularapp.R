#' Run Shiny App for SubcellularVis
#'
#' @param ... empty
#'
#' @import shiny 
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importFrom shinythemes shinytheme
#' @importFrom htmltools br p tagList img a strong
#' @importFrom colourpicker colourInput
#' @importFrom utils read.csv write.csv
#' @importFrom stats na.omit
#' @importFrom plotly plotlyOutput renderPlotly ggplotly
#' @importFrom stringi stri_split
#'
#' @return
#' @export
 
subcellularapp <- function(...){
  ui <- shiny::fluidPage(
    theme = shinythemes::shinytheme("cosmo"),
    #Title
    shiny::titlePanel("SubcellulaRVis: Visualising enrichment of protein localisation"),
    
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::selectInput(inputId = "input_organism",
                           label = "Organism",
                           choices = c("Human", "Mouse", "Drosophila", "Yeast", "Rat", "Xenopus"),
                           selected = "Human",
                           multiple = F),
        shiny::textAreaInput(inputId = "input_genes_text",
                             label = "Input genes,\none per line"),
        
        shiny::actionLink(inputId = "exampleGenes",
                          label = "Example gene set"),
        htmltools::br(),
        shiny::fileInput(inputId = "input_genes_file",
                         accept = ".csv",
                         label = "Input genes (file)",
                         multiple = F,
                         placeholder = ".csv file of genes; 1 column with HGNC symbols"),
        htmltools::br(),
        shiny::fileInput(inputId = "input_bkgd_file",
                         #accept = ".csv",
                         label = "Background population (optional)",
                         multiple = F,
                         placeholder = ".csv file of genes; 1 column with HGNC symbols", ),
        shiny::actionLink(inputId = "explainBkgd",
                          label = "What is the background?",
                          icon = shiny::icon("question")),
        
        
        shiny::checkboxInput(inputId = "traffic",
                             label = "Interested in trafficking?"),
        shiny::numericInput(
          inputId = "significanceThresh",
          label = "Stastical Signficance Threshold",
          min = 0, max = 0.05,
          value = 0.05),
        htmltools::br(),
        shiny::actionButton(inputId = "action_button",
                            label = "Calculate Enrichment & Visualise"),
        shiny::textOutput("checkInput")
      ),
      
      shiny::mainPanel(
        shiny::tabsetPanel(
          id = "tabsPanel",
          shiny::tabPanel(
            "About",
            htmltools::p(htmltools::strong("SubcellulaRVis a tool for visualising enrichment of\n
                   Gene Ontology Cellular Compartments within gene lists")),
            htmltools::tagList("Source code can be found at:", 
                               htmltools::a("www.github.com/jowatson2011/shiny_subcellularvis",
                                            href = "www.github.com/jowatson2011/shiny_subcellularvis")
            ),
            htmltools::br(),htmltools::br(),
            htmltools::p("Supported by:"),
            htmltools::img(src="https://bbsrc.ukri.org/bbsrc/includes/themes/BBSRC/images/logo-1.png",
                           width = "200px"),
            htmltools::img(src="https://www.staffnet.manchester.ac.uk/brand/visual-identity/logo/logo_big.gif",
                           width = "200px",)
          ),
          
          shiny::tabPanel("Plot", 
                          shiny::fluidRow(
                            shiny::column(width = 12,
                                          plotly::plotlyOutput(outputId = "plot_cell",
                                                               width = "600px")
                            )
                          ),
                          
                          shiny::fluidRow(
                            shiny::column(width = 6,
                                          htmltools::p(htmltools::strong("Colour scale:")),
                                          colourpicker::colourInput(inputId = "colScheme_low", 
                                                                    label = "Colour scale low",
                                                                    value = "#609377"),
                                          colourpicker::colourInput(inputId = "colScheme_high", 
                                                                    label = "Colour scale high",
                                                                    value = "#d8b2b2"),
                                          # shiny::textInput(inputId = "colScheme_low", 
                                          #           label = "Colour scale low",
                                          #           value = "#609377"),
                                          # shiny::textInput(inputId = "colScheme_high", 
                                          #           label = "Colour scale high",
                                          #           value = "#d8b2b2"),
                                          shiny::actionButton(inputId = "remakePlot",
                                                              label = "Remake Plot")
                            ),
                            shiny::column(width = 6,  
                                          htmltools::p(htmltools::strong("Options for export")),
                                          shiny::checkboxInput("includeLegend",
                                                               "Include legend?",
                                                               T),
                                          shiny::selectInput("legendPos",
                                                             "Legend positon",
                                                             list("bottom", "right", "left", "top"),
                                                             "right"),
                                          shiny::checkboxInput("includeLabels",
                                                               "Label compartments?",
                                                               T),
                                          shiny::numericInput("plotWidth", "Plot width, cm", 15),
                                          shiny::numericInput("plotHeight", "Plot height, cm", 10),
                                          shiny::numericInput("textSize", "Text size", 3),
                                          shiny::selectInput("export_fileType", "Image file type",
                                                             choices = list("pdf", "tiff", "png", "jpeg")),
                                          shiny::downloadButton(outputId = "downloadPlot", 
                                                                label = "Download")
                            )
                          ),
          ),
          shiny::tabPanel("Table", 
                          shiny::fluidRow(width = 6,
                                          shiny::tableOutput(outputId = "compartment_df")
                          ),
                          shiny::fluidRow(
                            shiny::column(width = 6,
                                          # htmltools::p(htmltools::strong("Export log")), 
                                          # htmltools::p("This can be used in manuscript methods session and aid reproducibility of your analysis"),
                                          # shiny::downloadButton(outputId = "downloadLog", 
                                          # label = "Export log.doc")
                                          htmltools::br(),
                                          shiny::downloadButton(outputId = "downloadEnrichment", 
                                                                label = "Export results")      
                            )
                            
                          ),
          ),
          shiny::tabPanel(
            "Full enrichment",
            htmltools::br(),
            
            shiny::fluidRow(
              shiny::column(width = 6,
                            shiny::tableOutput(outputId = "fullEnrichment_df")
              ),
              shiny::column(width = 6,
                            htmltools::br(),
                            htmltools::p("Top 50 terms shown"),
                            shiny::downloadButton(outputId = "downloadFullEnrichment", 
                                                  label = "Export full results")
              ),
            ),
          ),
          shiny::tabPanel(
            "Help",
            htmltools::br(),
            htmltools::p(htmltools::strong("How does SubCellulaRVis work?")),
            htmltools::p("SubCellulaRVis calculates the enrichment of . The broad categories of organelles and compartments were drawn from the Gene Ontology Cellular Compartments (GOCC). Using the GO hierachy, all \'child\' terms of these categories are grouped together to calculate a single enrichment score for each category. We can then visualise this in a single graphic."),
            htmltools::p(htmltools::strong("What is the background size?")),
            htmltools::p("The \"background size\" is the size of population being sampled (e.g. proteome size for a cell type or species). The default background size is the predicted number of proteins in the HeLa proteome."),
            htmltools::p(htmltools::strong("What is the difference between the simple and full enrichment?")),
            htmltools::p("The simple enrichment is based on the categorisation as described above. The full enrichment is a standard enrichment analysis based on all the GOCC terms."),
            htmltools::br(),
            htmltools::p(htmltools::strong("How can I export my results?")),
            htmltools::p("The visualisation can be exported in multiple image formats from the \"Plot\" tab.\n
            Enrichment results can be exported as a .csv file from the \"Full Enrichment\" tab.\n
            A log of the analysis can be exported as a .doc file from the \"Table\" tab"),
            htmltools::br(),
            htmltools::br(),
            htmltools::p("If you could not find your answer here, please contact: joanne.watson@manchester.ac.uk\n
            or raise an issue on the GitHub repository:" ),
            htmltools::tagList(
              htmltools::a("www.github.com/jowatson2011/shiny_subcellularvis",
                           href = "www.github.com/jowatson2011/shiny_subcellularvis")
            )
          )
        )
      )
    )
  )
  
  ### DEFINE SERVER LOGIC ###
  
  server <- function(input, output){
    
    shiny::observeEvent(input$exampleGenes, {
      
      # exampleGenes <- read.csv(paste0("data/",
      #                                 input$input_organism,
      #                                 "/exampleGenes.csv"),
      #                          header = F,
      #                          stringsAsFactors = F
      # )[,1]
      if(input$input_organism == "Mouse"){
        exampleGenes <- subcellularvis::Mouse_exampleGenes[,1, drop = T]
      } else if(input$input_organism == "Human"){
        exampleGenes <- subcellularvis::Human_exampleGenes[,1, drop = T]
      } else {
        updateSelectInput(inputId = "input_organism",
                          selected = "Human")
        exampleGenes <- subcellularvis::Human_exampleGenes[,1, drop = T]
      }
      
      shiny::updateTextAreaInput(inputId = "input_genes_text", 
                                 value = paste(exampleGenes, collapse = "\n")
      )
    })
    
    shiny::observeEvent(input$explainBkgd, {
      shiny::updateTabsetPanel(inputId = "tabsPanel", selected = "Help")
    })
    
    
    plot_cell <- shiny::reactiveVal()
    
    
    
    shiny::observeEvent(input$action_button, {
      
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
      output$checkInput <- shiny::renderText({ 
        if(!shiny::isTruthy(genes_tidy)){
          shiny::validate("HGNC gene symbols need to be entered,\n one per line\nas .csv or in the text box.")
        }
      })
      
      if(shiny::isTruthy(genes_tidy)){
        
        # getBkgd <- reactive({
        #   if(is.null(input$input_bkgd_file)) return(NULL)
        #   ......
        # })
        # 
        # bkgd <- if(is.null(getBkgd)){
        #   return(NULL)
        # } else {
        #   return(vroom::vroom(input$input_bkgd_file$datapath,
        #                       delim = ",", 
        #                       col_select = 1))
        # }
        bkgd_file <- input$input_bkgd_file
        
        if(!is.null(bkgd_file)){
          bkgd <<- read.csv(bkgd_file$datapath,
                            header = F,
                            stringsAsFactors = F)[,1]
        }else{
          bkgd <<- NULL
        }
        
        
        comps <<- compartmentData(genes = genes_tidy,
                                  bkgd = bkgd,
                                  trafficking = input$traffic,
                                  organism = input$input_organism)
      }
      
      # output$downloadLog <- shiny::downloadHandler(
      #   filename = function() {
      #     paste("CellVis_log_", Sys.Date(), ".doc", sep="")
      #   },
      #   content = function(file) {
      #   #  tempReport <- file.path(tempdir(), "log.Rmd")
      #  #   file.copy("data/log.Rmd", tempReport, overwrite = TRUE)
      #     
      #     # Set up parameters to pass to Rmd document
      #     params <- list(
      #       date = as.character(Sys.Date()),
      #       comps = comps,
      #       genes = genes_tidy
      #     )
      #     
      #     # Knit the document, passing in the `params` list, and eval it in a
      #     # child of the global environment (this isolates the code in the document
      #     # from the code in this app).
      #     rmarkdown::render(#tempReport,
      #       "data/log.Rmd", output_file = file,
      #       params = params,
      #       envir = new.env(parent = globalenv())
      #     )}
      # )
      if(shiny::isTruthy(genes_tidy)){
        shiny::updateTabsetPanel(inputId = "tabsPanel", selected = "Table")
      }
      
      output$compartment_df <- shiny::renderTable({
        shiny::withProgress(message = "Calculating", value = 15, {
          dplyr::mutate(comps, Symbol = sapply(.data$Symbol, function(i) {
                     
                     vec <- na.omit(strsplit(i, ",")[[1]][1:7])
                     
                     if(length(vec) == 0){
                       vec <- ""
                     }else if(length(na.omit(strsplit(i, ",")[[1]])) < 8){
                       vec <- paste(vec, collapse = ", ")
                     }else{
                       vec <- paste(c(vec, "..."), collapse = ", ")
                     }
                     
                     return(vec)
                     
                   })
          )
          # comps[,c("compartment",
          # "p",
          # "FDR",
          # "FDR < 0.05")]
        })
      }, 
      digits = 5)
      
      output$downloadEnrichment <- shiny::downloadHandler(
        filename = function() {
          paste("CellVis_Enrichment_", Sys.Date(), ".csv", sep="")
        },
        content = function(file){
          write.csv(comps, file)
        })
      
      
      
      output$plot_cell <- plotly::renderPlotly({
        shiny::withProgress(message = "Visualising", value = 15, {
          plot_cell(
            runSubcellulaRvis(comps,
                              colScheme_low = input$colScheme_low,
                              colScheme_high = input$colScheme_high,
                              trafficking = input$traffic,
                              text_size = 6,
                              legend = input$includeLegend,
                              legend.pos = input$legendPos,
                              labels = F)
          )
          plotly::ggplotly(plot_cell())
        })
      })
    })
    
    
    shiny::observeEvent(input$remakePlot, {
      shiny::withProgress(message = "Revisualisng", value = 15, {
        
        plot_cell(
          runSubcellulaRvis(comps,
                            colScheme_low = input$colScheme_low,
                            colScheme_high = input$colScheme_high,
                            text_size = 6,
                            trafficking = input$traffic,
                            legend = input$includeLegend,
                            legend.pos = input$legendPos,
                            labels = F)
        )
      })
    })
    
    
    # Export  visualisation
    output$downloadPlot <- shiny::downloadHandler(
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
    
    shiny::observeEvent(input$export_fileType,{
      output$downloadPlot <- shiny::downloadHandler(
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
    
    
    output$fullEnrichment_df <- shiny::renderTable({
      shiny::req(genes_tidy)
      shiny::withProgress(message = "Calculating", value = 15, {
        comps_full <<- compartmentData(genes = genes_tidy,
                                       bkgd = bkgd,
                                       trafficking = input$traffic,
                                       subAnnots = T,
                                       organism = input$input_organism)[1:50,]
        dplyr::mutate(comps_full, Symbol = 
                 sapply(.data$Symbol, function(i) {
                   
                   vec <- na.omit(strsplit(i, ",")[[1]][1:7])
                   
                   if(length(vec) == 0){
                     vec <- ""
                   }else if(length(na.omit(strsplit(i, ",")[[1]])) < 8){
                     vec <- paste(vec, collapse = ", ")
                   }else{
                     vec <- paste(c(vec, "..."), collapse = ", ")
                   }
                   
                   return(vec)
                   
                 })
        )
        #comps_full[,c("compartment",
        # "FDR",
        # "FDR < 0.05",
        # "group")]
      }) 
    },
    digits = 5)
    
    
    output$downloadFullEnrichment <- shiny::downloadHandler(
      filename = function() {
        paste("CellVis_fullEnrichment_", Sys.Date(), ".csv", sep="")
      },
      content = function(file){
        write.csv(comps_full, file)
      })
    
    
  }
  
  
  ### RUN THE APP ###
  shiny::shinyApp(ui = ui, server = server)
}