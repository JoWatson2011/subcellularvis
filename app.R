library(shiny)
source("R/runSubcellulaRvis.R")

### DEFINE UI FOR APP ###

ui <- fluidPage(

  #Title
  titlePanel("Visualise Subcellular Localisation of Omics Data"),

  sidebarLayout(
    sidebarPanel(textInput(inputId = "input_genes",
                           label = NULL),
                 
                  fileInput(inputId = "input_genes",
                           label = "Select file",
                           multiple = F),
                 
                 actionButton(inputId = "action",
                              label = "Visualise"),
                 
                 downloadButton(outputId = "downloadPlot", 
                                label = "Download")
                 ),

    mainPanel(
      plotOutput(outputId = "plot_cell"),
      br(),br(),
      p("Text placeholder for link to metadata"),
    )
  )
)

### DEFINE SERVE LOGIC ###

server <- function(input, output){
  
  plot_cell <- reactiveVal()
  
  observeEvent(input$action, {

    # genes_tidy <- stringi::stri_split(genes, 
    #                                   ";",
    #                                   fixed = T,
    #                                   simplify=T)
    
    comps <- compartmentData(genes = input$input_genes)

    plot_cell(runSubcellulaRvis(comps))
  }
  )
 
  output$plot_cell <- renderPlot({
    plot_cell()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = 'test.png',
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(..., width = width, height = height,
                       res = 300, units = "in")
      }
      ggsave(file, plot = plotInput(), device = device)
    }
  )
  
}


### RUN THE APP ###
shinyApp(ui = ui, server = server)
