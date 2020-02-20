library(shiny)
source("R/runSubcellulaRvis.R")

### DEFINE UI FOR APP ###

ui <- fluidPage(

  #Title
  titlePanel("Visualise Subcellular Localisation of 'Omics Data"),

  sidebarLayout(
    sidebarPanel(textAreaInput(inputId = "input_genes_text",
                               label = "List of genes, \none per line"),
                 
                 actionButton(inputId = "action_text",
                              label = "From textbox"),
                 
                  fileInput(inputId = "input_genes_file",
                           label = "Select file",
                           multiple = F,
                           placeholder = ".txt file of genes"),
                 
                 actionButton(inputId = "action_file",
                              label = "From .txt file")
                 ),

    mainPanel(
      plotOutput(outputId = "plot_cell"),
      br(),br(),
      downloadButton(outputId = "downloadPlot", 
                     label = "Download"),
      p("Text placeholder"),
    )
  )
)

### DEFINE SERVER LOGIC ###

server <- function(input, output){
  
  plot_cell <- reactiveVal()
  
  # visualise from textbox
  observeEvent(input$action_text, {

    genes_tidy <- as.vector(stringi::stri_split(input$input_genes_text,
                                      regex = "\n",
                                      simplify=T))
    
    comps <- compartmentData(genes = genes_tidy)

    plot_cell(runSubcellulaRvis(comps))
  }
  )
  
  # Visualise from file
  observeEvent(input$action_file, {
    
    genes_df <- read.csv(input$input_genes_file$datapath,
                           header = F,
                           stringsAsFactors = F)

    comps <- compartmentData(genes = genes_df[,1])
    
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
      ggsave(file, device = device)
    }
  )
  
}


### RUN THE APP ###
shinyApp(ui = ui, server = server)
