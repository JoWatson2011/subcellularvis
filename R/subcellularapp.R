#' Run Shiny App for SubcellularVis
#'
#' @param ... empty
#'
#' @import shiny 
#' @importFrom dplyr mutate left_join select
#' @importFrom rlang .data
#' @importFrom shinythemes shinytheme
#' @importFrom htmltools br p tagList img a strong
#' @importFrom colourpicker colourInput
#' @importFrom utils read.csv write.csv
#' @importFrom stats na.omit
#' @importFrom plotly plotlyOutput renderPlotly ggplotly
#' @importFrom stringi stri_split
#' @importFrom purrr reduce splice
#' @importFrom formattable formattable color_tile formattableOutput renderFormattable
#' @importFrom grDevices dev.off jpeg png tiff
#' @return Shiny App
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
                           choices = c("Human", "Mouse", "Drosophila", 
                                       "Yeast", "Rat", "Xenopus"),
                           selected = "Human",
                           multiple = F),
        conditionalPanel(
          condition = "input.input_organism == 'Human'",
          shiny::selectInput(inputId = "annotationSource",
                             label = "Annotation Source",
                             choices = c("Gene Ontology", 
                                         "Human Protein Atlas"),
                             selected = "Gene Ontology",
                             multiple = F
          )
        ),
        shiny::selectInput(inputId = "identifierType",
                           label = "Identifier type",
                           choices = c("HGNC Symbol",
                                       "UniProt ID"),
                           selected = "HGNC Symbol",
                           multiple = F),
        shiny::actionLink(inputId = "explainIdentifier",
                          label = "What if my list uses different identifiers?",
                          icon = shiny::icon("question")),
        htmltools::br(),
        shiny::actionLink(inputId = "exampleGenes",
                          label = "Example gene set"),
        shiny::textAreaInput(inputId = "input_genes_text",
                             label = "Input genes,\none per line"),
        shiny::fileInput(inputId = "input_genes_file",
                         accept = ".csv",
                         label = "Input genes (file)",
                         multiple = F,
                         placeholder = ".csv file of genes; 1 column, one gene per row"),
        shiny::fileInput(inputId = "input_bkgd_file",
                         #accept = ".csv",
                         label = "Background population (optional)",
                         multiple = F,
                         placeholder = ".csv file of genes; 1 column, one gene per row", ),
        shiny::actionLink(inputId = "explainBkgd",
                          label = "What is the background?",
                          icon = shiny::icon("question")),
        
        shiny::selectInput(inputId = "traffic",
                           "Visualisation of...",
                           c("Whole cell","Endosomal system")),
        # shiny::checkboxInput(inputId = "traffic",
        #                      label = "Interested in trafficking?"),
        shiny::numericInput(
          inputId = "significanceThresh",
          label = "FDR Significance Threshold",
          min = 0, max = 0.05,
          value = 0.05),
        htmltools::br(),
        withBusyIndicatorUI(
          shiny::actionButton(inputId = "action_button",
                              label = "Calculate Enrichment & Visualise",
                              class = "btn-primary" )
        ),
        shiny::textOutput("checkInputEmpty"),
        tags$head(tags$style("#checkInputEmpty{color: red;
                                  }"
        )),
        shiny::textOutput("checkInputLength"),
        tags$head(tags$style("#checkInputLength{color: red;
                                  }"
        )),
        shiny::textOutput("checkInputFormat"),
        tags$head(tags$style("#checkInputLength{color: red;
                                  }"
        )
        )
      ),
      
      shiny::mainPanel(
        shiny::tabsetPanel(
          id = "tabs",
          shiny::tabPanel(
            align="center",
            "About",
            htmltools::p(htmltools::strong("SubcellulaRVis is a tool for visualising enrichment of\n
                   Gene Ontology Cellular Compartments within gene lists")),
            htmltools::img(src= "https://raw.githubusercontent.com/JoWatson2011/subcellularvis/master/data-raw/CELL.png",
                           width = "300px"
            ),
            htmltools::img(src= "https://raw.githubusercontent.com/JoWatson2011/subcellularvis/master/data-raw/endosomalCell.png",
                           width = "275px"
            ),
            htmltools::p("Input your gene or protein list as a text or .csv input in the box to the left.\n
                         Make sure to select the correct organism and gene/protein identifier type.\n
                         You can calculate enrichment based on the whole cell (left above) or the endosomal system (right above)"),
            htmltools::tagList("Source code for this app can be found at:", 
                               htmltools::a("www.github.com/jowatson2011/subcellularvis",
                                            href = "www.github.com/jowatson2011/subcellularvis")
            ),
            htmltools::br(),htmltools::br(),
            htmltools::p("Supported by:"),
            htmltools::img(src="https://www.ukri.org/wp-content/uploads/2020/06/our-council-logo-bbsrc.png",
                           width = "75px"),
            htmltools::img(src="https://www.staffnet.manchester.ac.uk/brand/visual-identity/logo/logo_big.gif",
                           width = "200px"),
            htmltools::img(src="https://upload.wikimedia.org/wikipedia/commons/thumb/5/58/Wellcome_Trust_logo.svg/800px-Wellcome_Trust_logo.svg.png",
                           width = "75px")
          ),
          shiny::tabPanel("Table",
                          shiny::fluidRow(
                            shiny::column(width = 6,
                                          textOutput("checkGenesMap1"),
                                          tags$head(tags$style("#checkGenesMap1{color: red;
                                  }"
                                          )
                                          )
                            )
                          ),
                          shiny::fluidRow(
                            shiny::column(width = 6,
                                          textOutput("noGenesMapped1"),
                                          tags$style(type="text/css", "#unmappedGenes1 {word-wrap: break-word;word-break: break-all;white-space: pre-wrap;}"),
                                          textOutput("unmappedGenes1"),
                                          # htmlOutput("unmappedGenes1")
                            )
                          ),
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
            "Plot", 
            shiny::fluidRow(
              shiny::column(width = 6,
                            textOutput("checkGenesMap3"),
                            tags$head(tags$style("#checkGenesMap3{color: red;
                                  }"
                            )
                            )
              )
            ),
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
                                                label = "Remake Plot",
                                                class = "btn-primary")
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
                            shiny::selectInput("export_fileType", 
                                               "Export image type",
                                               choices = list("pdf",
                                                              "svg",
                                                              "tiff", 
                                                              "png",
                                                              "jpeg")),
                            shiny::downloadButton(outputId = "downloadPlot", 
                                                  label = "Download")
              )
            ),
          ),
          
          shiny::tabPanel(
            "Full enrichment",
            htmltools::p("Note: This enrichment is performed using GO annotations,
                         even if Human Protein Atlas has been selected."),
            htmltools::br(),
            
            shiny::fluidRow(
              shiny::column(width = 6,
                            textOutput("checkGenesMap2"),
                            tags$head(tags$style("#checkGenesMap2{color: red;
                                  }"
                            )
                            ))
            ),
            shiny::fluidRow(
              shiny::column(width = 6,
                            textOutput("noGenesMapped2"),
                            # tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
                            tags$style(type="text/css", "#unmappedGenes2 {word-wrap: break-word;word-break: break-all;white-space: pre-wrap;}"),
                            textOutput("unmappedGenes2")
              )
            ),
            shiny::fluidRow(
              shiny::column(width = 6,
                            shiny::tableOutput(outputId = "fullEnrichment_df")
              ),
              shiny::column(width = 6,
                            htmltools::br(),
                            htmltools::p("Top 50 terms shown"),
                            shiny::downloadButton(outputId = "downloadFullEnrichment", 
                                                  label = "Export full results")
              )
            )
          ),
          shiny::tabPanel(
            "Overlap",
            shiny::fluidRow(
              column(
                width = 12,
                  htmltools::HTML("<p>More details on how to interpret Upset plots 
                                  can be found here:
                                  <a href = https://upset.app/> https://upset.app/</a> </p>")
              )
            ),
            shiny::fluidRow(
              shiny::column(width = 12, 
                            shiny::plotOutput("upsetPlot"),
                            shiny::textOutput("noIntersects"),
                            height = "600px"
              )
            ),
            shiny::fluidRow(
              shiny::column(width = 12,
                            htmltools::br(),
                            shiny::tableOutput("upsetData")
              )
            ),
            shiny::fluidRow(
              shiny::column(width = 3,
                            shiny::selectInput("exportUpset_fileType",
                                               "Image file type",
                                               choices = list("tiff", 
                                                              "png",
                                                              "jpeg")),
              ),
              shiny::column(width = 3,
                            shiny::downloadButton(outputId = "downloadUpsetPlot", 
                                                  label = "Download plot")
              ), 
              shiny::column(width = 3,
                            shiny::downloadButton(outputId = "downloadUpsetData", 
                                                  label = "Download data"))
            )
          ),
          shiny::tabPanel(
            "Multiple comparisons",
            fluidRow(
              htmltools::p("")
            ),
            fluidRow(
              column(
                width = 6,
                selectInput("inputType",
                            "Input type",
                            c("Text","File"),
                            "Text"),
                sliderInput("numInputs", "How many inputs do you want", 1, 10,1),
                # place to hold dynamic inputs
              ),
              column(
                width = 6,
                uiOutput("inputGroup"),
                #   withBusyIndicatorUI(
                actionButton("compareAction", 
                             "compare",
                             #                 class = "btn-primary"
                # )
              )
            )
          ),
          fluidRow(
            column(
              width = 12,
              tableOutput("compDatComparison"),
              shiny::downloadButton(outputId = "downloadMultipleComparison", 
                                    label = "Download table")
            )
          )
        ),
        shiny::tabPanel(
          "Help",
          htmltools::br(),
          htmltools::p("View our publication describing SubcellulaRVis:"),
          htmltools::tagList(
            htmltools::a("https://doi.org/10.1093/nar/gkac336",
                         href = "https://doi.org/10.1093/nar/gkac336")
          ),
          htmltools::br(),
          htmltools::br(),
          htmltools::p(htmltools::strong("How does SubCellulaRVis work?")),
          htmltools::p("SubCellulaRVis calculates the enrichment of proteins
                         for cellular compartments, inferred from the Gene
                         Ontology Cellular Compartment aspect. The broad
                         categories of organelles and compartments were drawn
                         from the Gene Ontology Cellular Compartments (GOCC).
                         Using the GO hierachy, all \'child\' terms of these
                         categories are grouped together to calculate a single
                         enrichment score for each category. We can then visualise
                         this in a single graphic."),
          htmltools::p(htmltools::strong("How is the gene list processed?")),
          htmltools::p("Duplicated genes are removed, and identifiers are then mapped to their
                       annotations in the Gene Ontology or Human Protein Atlas, as specified.
                       The number of successfully mapped genes and a list of unmapped genes can
                       be found at the top of the \'Table\' and \'Full Enrichment\' tabs"),
          htmltools::p(htmltools::strong("How do I map my gene list to HGNC symbols or Uniprot ID?")),
          htmltools::p("The BioMart service (https://m.ensembl.org/biomart/martview/) can
                         be used to map identifiers. Full instructions can be found on the
                         Ensembl website (https://m.ensembl.org/info/data/biomart/index.html)."),
          htmltools::p(htmltools::strong("What is the background population?")),
          htmltools::p("The \"background population\" is the total population
                         of expressed proteins in the sampled (e.g. proteome 
                         size for a cell type or species). The default 
                         background is the genes in the
                         reference genomes."),
          htmltools::p(htmltools::strong("What is the difference between the
                                           simple and full enrichment?")),
          htmltools::p("The simple enrichment is based on the categorisation
                         as described above. The full enrichment is a standard
                         enrichment analysis based on all the GOCC terms."),
          htmltools::p(htmltools::strong("How do I interpret the results?")),
          htmltools::p("SubcellulaRVis provides a simplified version of Gene 
                         Ontology Cellular Compartment (GOCC) enrichment results.
                         Enrichment analyses are based on over-representation
                         of proteins or genes annotated to annotations such as
                         GOCC terms. SubcellulaRvis projects these results onto 
                         a schematic of the cell (or subcellular system). The 
                         compartments are representative of a set of GOCC
                         terms that describe that compartment; enrichment for
                         that compartments suggests enrichment in your gene list
                         of some in that set of terms. The full set of GOCC terms
                         can be found in the \'Full Enrichment\' tab."),
          htmltools::p(htmltools::strong("How can I export my results?")),
          htmltools::p("The visualisation can be exported in multiple image
            formats from the \"Plot\" tab.\n
            Enrichment results can be exported as a .csv file from the
            \"Full Enrichment\" tab.\n
            A log of the analysis can be exported as a .doc file from the
                         \"Table\" tab"),
          htmltools::p(htmltools::strong("Can I run SubcellulaRVis locally or programmatically?")),
          htmltools::HTML("<p>An R package for SubcellulaRVis can be used for both of these purposes. The package
                       is simple to use and can be downloaded from GitHub, where details on how to use the 
                       package can be found:<a href =https://github.com/JoWatson2011/subcellularvis>
                          https://github.com/JoWatson2011/subcellularvis</a>.</p>"),
          htmltools::br(),
          htmltools::p("If you could not find your answer here, please 
            contact: jean-marc.schwartz@manchester.ac.uk\n
            or raise an issue on the GitHub repository:" ),
          htmltools::tagList(
            htmltools::a("www.github.com/jowatson2011/subcellularvis",
                         href = "www.github.com/jowatson2011/subcellularvis")
          )
        )
      )
    )
  )
  )

### DEFINE SERVER LOGIC ###

server <- function(input, output){
  
  shiny::observe({
    if(input$annotationSource == "Human Protein Atlas"){
      updateSelectInput(inputId = "traffic",
                        label = "Visualisation of...",
                        choices = "Whole cell",
                        selected ="Whole cell"
      )
    }
  })
  
  shiny::observeEvent(input$exampleGenes, {
    
    if(input$input_organism == "Mouse"){
      exampleGenes <- subcellularvis::Mouse_exampleGenes[,1, drop = T]
      updateSelectInput(inputId = "identifierType",
                        selected = "HGNC Symbol")
    } else if(input$input_organism == "Human"){
      exampleGenes <- subcellularvis::Human_exampleGenes[,1, drop = T]
      shiny::updateSelectInput(inputId = "identifierType",
                               selected = "HGNC Symbol")
    } else {
      shiny::updateSelectInput(inputId = "input_organism",
                               selected = "Human")
      shiny::updateSelectInput(inputId = "identifierType",
                               selected = "HGNC Symbol")
      exampleGenes <- subcellularvis::Human_exampleGenes[,1, drop = T]
    }
    
    shiny::updateTextAreaInput(inputId = "input_genes_text", 
                               value = paste(exampleGenes, collapse = "\n")
    )
  })
  
  shiny::observeEvent(input$explainBkgd, {
    shiny::updateTabsetPanel(inputId = "tabs", selected = "Help")
  })
  
  shiny::observeEvent(input$explainIdentifier,{
    shiny::updateTabsetPanel(inputId = "tabs", selected = "Help")
  })
  
  shiny::observeEvent(input$action_button, {
    withBusyIndicatorServer("action_button", {
      plot_cell <- shiny::reactiveVal()
      
      file <- input$input_genes_file
      
      genes_tidy <<- ""
      
      if(!is.null(file)){
        genes_tidy <<- read.csv(file$datapath,
                                header = F,
                                stringsAsFactors = F)[,1]
        genes_tidy <<- genes_tidy[genes_tidy != ""]
      }else{
        genes_tidy <<- as.vector(stringi::stri_split(input$input_genes_text,
                                                     regex = "\n",
                                                     simplify=T))
        genes_tidy <<- genes_tidy[genes_tidy != ""]
      }
      output$checkInputEmpty <- shiny::renderText({
        if(!shiny::isTruthy(genes_tidy)){
          shiny::validate("Genes need to be entered,\n one per line\nas .csv or in the text box.")
        }
      })
      output$checkInputLength <- shiny::renderText({
        if(length(genes_tidy) == 1){
          # as.character(icon("exclamation-circle"))
          shiny::validate("More than one gene needs to be entered.")
        }
      })
      output$checkInputFormat <- shiny::renderText({
        if(any(grepl("[,;.]", genes_tidy))){
          shiny::validate(paste(
            "Genes need to be entered,\n one per line\nas .csv or in the text box.
            Error at line:",
            which(sapply(c("MAPK1,MAPK3", "PI3K"), function(i) grepl(",", i)))
          ))
        }
      })
      
      if(shiny::isTruthy(genes_tidy) & length(genes_tidy) > 1){
     
        bkgd_file <- input$input_bkgd_file
        
        if(!is.null(bkgd_file)){
          bkgd <<- read.csv(bkgd_file$datapath,
                            header = F,
                            stringsAsFactors = F)[,1]
        }else{
          bkgd <<- NULL
        }
        
        shiny::withProgress(message = "Calculating", value = 10, {
          tr <- isolate(input$traffic)
          id <- isolate(input$identifierType)
          comps <<- subcellularvis::compartmentData(
            genes = genes_tidy,
            bkgd = bkgd,
            id_type = ifelse(id == "HGNC Symbol", "SYMBOL", "UNIPROT"),
            aspect = tr,
            organism = input$input_organism,
            annotationSource = input$annotationSource,
            significanceThresh = input$significanceThresh)
          
          comps_full <<- subcellularvis::compartmentData(
            genes = genes_tidy,
            bkgd = bkgd,
            id_type = ifelse(id == "HGNC Symbol", "SYMBOL", "UNIPROT"),
            #trafficking = input$traffic,
            aspect = tr,
            subAnnots = T,
            organism = input$input_organism,
            annotationSource = input$annotationSource,
            significanceThresh = input$significanceThresh)
        })
      }
      
      output$noGenesMapped1 <-  
        shiny::renderText({
          if(!is.null(comps)){
            paste("Number of genes mapped:", comps$nMapped)
          }
        })
      output$noGenesMapped2 <- shiny::renderText({
        if(!is.null(comps)){
          paste("Number of genes mapped:", comps_full$nMapped)
        }
      })
      output$unmappedGenes1 <- 
        shiny::renderText({
          if(!is.null(comps)){
            paste("Genes not mapped:", comps$unmapped)
          }
        })
      output$unmappedGenes2 <- 
        shiny::renderText({
          if(!is.null(comps)){
            paste("Genes not mapped:", comps_full$unmapped)
          }
        })
      output$checkGenesMap1 <- shiny::renderText({
        if(is.null(comps)){
          shiny::validate("Gene names did not map. Did you use
                          select correct organism and identifier type?")
        }
      })
      output$checkGenesMap2 <- shiny::renderText({
        if(is.null(comps)){
          shiny::validate("Gene names did not. Did you use
                          select correct organism and identifier type?")
        }
      })
      output$checkGenesMap3 <- shiny::renderText({
        if(is.null(comps)){
          shiny::validate("Gene names did not. Did you use
                          select correct organism and identifier type?")
        }
      })

      if(shiny::isTruthy(genes_tidy)){
        shiny::updateTabsetPanel(inputId = "tabs", selected = "Table")
      }
      
      output$compartment_df <- shiny::renderTable({
        req(comps)
        
        dplyr::mutate(comps$enrichment, 
                      Genes = sapply(comps$enrichment[,6], function(i) {
                        
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
      }, 
      digits = 5)
      
      output$downloadEnrichment <- shiny::downloadHandler(
        filename = function() {
          paste("CellVis_Enrichment_", Sys.Date(), ".csv", sep="")
        },
        content = function(file){
          write.csv(comps$enrichment, file)
        },
        contentType = "text/csv")
      
      output$plot_cell <- plotly::renderPlotly({
        req(comps)
        tr <- isolate(input$traffic)
        sig <- isolate(input$significanceThresh)
        shiny::withProgress(message = "Visualising", value = 15, {
          plot_cell(
            runSubcellulaRvis(comps$enrichment,
                              colScheme_low = input$colScheme_low,
                              colScheme_high = input$colScheme_high,
                              aspect = tr,
                              text_size = 6,
                              legend = input$includeLegend,
                              legend.pos = input$legendPos,
                              labels = F,
                              significanceThresh = sig)
          )
          plotly::ggplotly(plot_cell())
        })
      })
      
      output$fullEnrichment_df <- shiny::renderTable({
        req(comps_full)
        dplyr::mutate(
          comps_full$enrichment[1:50,],
          Genes =
            sapply(.data$Genes, function(i) {
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
      },
      digits = 5)
      
      if(!is.null(comps)){
        upset <<- subcellularvis::plotOverlap(comps$enrichment)
      }
      output$noIntersects <- renderText({
        paste(
          ifelse(nrow(upset$upsetDat) < 15, nrow(upset$upsetDat), 15),
          "/",
          nrow(upset$upsetDat), "intersects plotted")
      })
      output$upsetPlot <- renderPlot({
        req(upset)
        upset$upsetPlot
      })
      output$upsetData <- renderTable({
        req(upset)
        dplyr::mutate(
          upset$upsetDat,
          Genes =
            sapply(.data$Genes, function(i) {
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
      })
      
      # observe changes in "numInputs", and create corresponding number of inputs
      observeEvent(input$numInputs, {
        output$inputGroup = renderUI({
          input_list <- lapply(1:input$numInputs, function(i) {
            # for each dynamically generated input, give a different name
            inputName <- paste("Input genes", i, sep = "")
            if(input$inputType == "Text"){
              textAreaInput(inputName, 
                            label = inputName)
            } else {
              fileInput(inputId = inputName,
                        label =  inputName,
                        accept = ".csv")
            }
          })
          do.call(tagList, input_list)
        })
      })
      
      observeEvent(input$compareAction, {
        num <- isolate(input$numInputs)
        tr <- isolate(input$traffic)
        id <- isolate(input$identifierType)
        
        dat_list <- lapply(1:num, function(i){
          inputName <- paste("Input genes", i, sep = "")
          
          if(input$inputType == "Text"){
            compGenes_tidy <<- as.vector(stringi::stri_split(input[[inputName]],
                                                             regex = "\n",
                                                             simplify=T))
          }else{
            compGenes_tidy <<- read.csv(input[[inputName]]$datapath,
                                        header = F,
                                        stringsAsFactors = F)[,1]
          }
          compGenes_tidy <<- compGenes_tidy[compGenes_tidy != ""]
          
          compCompartmentDat <- compartmentData(compGenes_tidy,
                                                bkgd = bkgd,
                                                id_type = ifelse(id == "HGNC Symbol", 
                                                                 "SYMBOL", 
                                                                 "UNIPROT"),
                                                aspect = tr,
                                                organism = input$input_organism,
                                                annotationSource = input$annotationSource,
                                                significanceThresh = input$significanceThresh
          )$enrichment %>% 
            select(.data$Compartment, .data$FDR)
          colnames(compCompartmentDat)[2] <- paste(colnames(compCompartmentDat)[2], i)
          return(compCompartmentDat)
        })
        
        dat_joined <<- purrr::reduce(
          purrr::splice(
            comps$enrichment[,c("Compartment","FDR")],
            dat_list
          ),
          left_join, 
          by = "Compartment")
        
        output$compDatComparison <- shiny::renderTable({
          dat_joined
          #      })
        })
      })
    })
  })
  
  shiny::observeEvent(input$remakePlot, {
    shiny::withProgress(message = "Revisualising", value = 15, {
      plot_cell <- shiny::reactiveVal()
      tr <- isolate(input$traffic)
      plot_cell(
        runSubcellulaRvis(comps$enrichment,
                          colScheme_low = input$colScheme_low,
                          colScheme_high = input$colScheme_high,
                          text_size = 6,
                          aspect = input$traffic,
                          legend = input$includeLegend,
                          legend.pos = input$legendPos,
                          labels = F,
                          significanceThresh = input$significanceThresh)
      )
    })
  })
  
  
  # Export  visualisations & data
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
        paste("CellVis_", Sys.Date(),".", 
              input$export_fileType, sep="")
      },
      content = function(file) {
        ggplot2::ggsave(file, 
                        plot = runSubcellulaRvis(comps$enrichment,
                                                 colScheme_low = input$colScheme_low,
                                                 colScheme_high = input$colScheme_high,
                                                 text_size = input$textSize,
                                                 aspect =  input$traffic,
                                                 legend = input$includeLegend,
                                                 legend.pos = input$legendPos,
                                                 labels = input$includeLabels,
                                                 significanceThresh = input$significanceThresh),
                        device = input$export_fileType,
                        width = input$plotWidth,
                        height = input$plotHeight,
                        dpi = 300,
                        units = "cm")
      },
      contentType = ifelse(input$export_fileType == "pdf", 
                           "application/pdf",
                           ifelse(input$export_fileType == "svg", 
                                  "svg+xml", 
                                  paste0("image/", 
                                         input$export_fileType)
                           )
      )
    )
  })
  
  output$downloadUpsetPlot <- shiny::downloadHandler(
    filename = function() {
      paste("CellVis_Upset_", Sys.Date(), ".",
            input$exportUpset_fileType, sep="")
    },
    content = function(file) {
      
      # if(names(dev.cur()) != "null device"){
      #   dev.off()
      # }
      
      if(input$exportUpset_fileType == "png"){
        png(file, width = 600)
        print(upset$upsetPlot)
        dev.off()
      }else if(input$exportUpset_fileType == "tiff"){
        tiff(file, width = 600)
        print(upset$upsetPlot)
        dev.off()
      }else if(input$exportUpset_fileType == "jpeg"){
        jpeg(file, width = 600)
        print(upset$upsetPlot)
        dev.off()
      }
    },
    contentType = paste0("image/", input$exportUpset_fileType)
  )
  
  output$downloadUpsetData <- shiny::downloadHandler(
    filename = function() {
      paste("CellVis_UpsetData_", Sys.Date(), ".csv", sep="")
    },
    content = function(file){
      write.csv(upset$upsetDat, file)
    },
    contentType = "text/csv")
  
  output$downloadFullEnrichment <- shiny::downloadHandler(
    filename = function() {
      paste("CellVis_fullEnrichment_", Sys.Date(), ".csv", sep="")
    },
    content = function(file){
      write.csv(comps_full$enrichment, file)
    },
    contentType = "text/csv")
  
  output$downloadMultipleComparisons <- shiny::downloadHandler(
    filename = function() {
      paste("CellVis_multipleComparisons_", Sys.Date(), ".csv", sep="")
    },
    content = function(file){
      write.csv(dat_joined, file)
    },
    contentType = "text/csv")
  
}

### RUN THE APP ###
shiny::shinyApp(ui = ui, server = server)
}