library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  titlePanel("Real-time PCR data analysis tool"),
  sidebarLayout(
    sidebarPanel(
      helpText("Step 1: Please select the file with reference gene names"),
      fileInput("file2", label = "Choose the file", multiple =
                  T),
      helpText("Step 2: Please select any number of reference genes from the list"),
      uiOutput("refgen"),
      
      br(),
      helpText("Step 3: Please choose the Real-time PCR data"),
      fileInput("file1", label = "Choose the files", multiple =
                  T),
      p("Step 4: Press this button to perform the analysis"),
      actionButton("act","ANALYSE"),
      br(),
      br(),
      p("Press this button in order to download the results"),
      downloadButton("down","Download the results")
      
      ),
    
    
    mainPanel(
      tabsetPanel(
        tabPanel("Tool description", helpText("This is a tool which allows to performa a simple and quick Real-time PCR data analysis."),
                 helpText("The following files are necessary:"),
                 helpText("- the file with reference gene names,"),
                 img(src = "ref_names.png", height = 250, width = 150),
                 helpText("- the Real-time PCR data."),
                 img(src = "qPCR.png", height = 250, width = 400),
                 helpText("The file name must have the following format: gene_pow2_other"),
                 helpText("Where:"),
                 helpText("gene - gene name,"),
                 helpText("pow2 - for second repeat only, otherwise this part should be ommited,"),
                 helpText("other - any additional information."),
                 helpText("The output table will be visible in the Results tab, while the boxplots will be available in their corresponding parts."),
                 helpText("If you wish to save the results, press the 'Download the results' button."),
                 helpText("Note: The program assumes that each plate has the same samples."),
                 br(),
                 br()
                 
                 ),
        tabPanel("Results", tableOutput("tabelka")),
        tabPanel("Boxplots - genes", plotOutput("genplot")),
        tabPanel("Boxplots - patients", plotOutput("patplot")))
      

    )
  )
))