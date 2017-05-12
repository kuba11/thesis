library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  titlePanel("Simple qPCR"),
  sidebarLayout(
    sidebarPanel(
      helpText("Step 1: Please select the file with efficiencies"),
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
        tabPanel("Tool description", helpText("This is a tool which allows to perform a simple and quick Real-time PCR data analysis.
                 It uses the Pfaffl method for quantification and the efficiency values from the experiment."),
                 helpText("The following files are necessary:"),
                 helpText("- the file with gene list and efficiencies,"),
                 img(src = "eff.png", height = 210, width = 120),
                 helpText("- the Real-time PCR data in SDS 2.3 format."),
                 img(src = "qPCR.png", height = 300, width = 450),
                 helpText("The file names of the data must have the following format: gene_other"),
                 helpText("Where:"),
                 helpText("gene - gene name,"),
                 helpText("other - any additional information."),
                 helpText("The output table will be visible in the Results tab, while the boxplots will be available in the boxplots tab."),
                 helpText("If you wish to save the results, press the 'Download the results' button."),
                 helpText("Note: The program assumes that each plate has the same samples."),
                 br(),
                 br()
                 
                 ),
        tabPanel("Advanced options", radioButtons(
          "configtype",
          label = "File type:",
          choices = list("Efficiency file", "Full configuration file"), selected = "Full configuration file"
        )),
        tabPanel("Results", tableOutput("tabelka")),
        tabPanel("Boxplots - samples", uiOutput("patplot"))
        )
      

    )
  )
))