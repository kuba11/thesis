library(shiny)


shinyUI(fluidPage(
  titlePanel("Simple qPCR"),
  sidebarLayout(
    sidebarPanel(
      helpText("Step 1: Please select the configuration file with reaction efficiencies"),
      fileInput("file2", label = "Choose the file", multiple =
                  T),
      helpText("Step 2: Please select any number of reference genes from the list"),
      uiOutput("refgen"),
      
      br(),
      helpText("Step 3: Please choose the Real-time PCR data"),
      fileInput("file1", label = "Choose the files", multiple =
                  T),
      
      helpText("Step 4: Please select the calibration samples"),
      uiOutput("control"),
      
      p("Step 5: Press this button to perform the analysis"),
      actionButton("act","ANALYSE"),
      br(),
      br(),
      p("Press this button in order to download the results"),
      downloadButton("down","Download the results")
      
      ),
    
    
    mainPanel(
      tabsetPanel(
        tabPanel("Tool description",
helpText("Simple qPCR is a tool which allows to perform a simple and quick Real-time PCR data analysis.
                 It uses relative quantification method utilising the Pfaffl model and efficiency values from the experiment."),
helpText("The scheme of the used algorithm can be seen below:"),                 
img(src = "algorithm.png", height = 400, width = 450),
helpText("The following files are necessary:"),
                 helpText("- the Real-time PCR data in SDS 2.3 format."),
                 img(src = "qPCR.png", height = 300, width = 450),
                 br(),
                 downloadButton("downdata","Download sample data"),
                 helpText("- the configuraion file,"),
                 helpText('The configuration file should be a XLSX file and contain:'),
                 helpText('- efficiency values (sheet 1)'), 
                 helpText('- gene names (sheet 2)'),
                 helpText('- set number for each gene (sheet 3)'),
                 helpText('- information about sample placement in each set (sheet 4, picture below).'),
                 img(src = "conf.png", height = 300, width = 200),
                 br(),
                 downloadButton("downconf","Download a sample configuration file"),
                 downloadButton("downconfe","Download an empty configuration file"),
                 helpText("If all of your files have the same samples in the same order, you can provide a simple file with efficiencies instead. Check the 'Efficiency' tab for more info."),
                 br(),
                 helpText("You can also model the fluorescence curve and calculate the Ct value based on the raw fluorescence data. Check the 'Raw fluorescence data' tab for more info."),
                 br(),
                 helpText("The file names of the data must have the following format: gene_other"),
                 helpText("Where:"),
                 helpText("gene - gene name,"),
                 helpText("other - any additional information."),
                 helpText("The output table will be visible in the Results tab, while the boxplots will be available in the boxplots tab."),
                 helpText("If you wish to save the results, press the 'Download the results' button."),
                 br(),
                 br()
                 
                 ),
        tabPanel("Efficiency",
                 helpText('Here you can specify if you want to use a simple efficiency file as a configuration file.'),
                 helpText('The file should contain only one sheet with a list of gene names and efficiencies as seen below,
                 as in this case it is assumed that the same samples are present in each of the files and they are in the same order.'),
                 checkboxInput("eff.file", "Use a simple efficiency file", value = FALSE, width = NULL),
                 img(src = "eff.png", height = 210, width = 120),
                 br(),
                 downloadButton("downeff","Download a sample efficiency file"),
                 downloadButton("downeffe","Download an empty efficiency file"),
                 br()
                 ),

        tabPanel("Raw fluorescence data",
         helpText("Raw fluorescence values from the real-time PCR experiment can be used to create a regression function model that allow to determine precise Ct values.
                  If you don't have the files with Ct values or you wish to calculate it from the models, you can chack the box below upload the raw fluorescence files in the 'Step 3'."),
              checkboxInput("fluo.file", "Use raw fluorescence data", value = FALSE, width = NULL),
         helpText("Please select the threshold value of fluorescence that will be used used to determine the threshold cycle (Ct)."),
         numericInput('fluo.thresh', 'Fluorescence threshold', 0.5, min = 0, step = 0.1, max = 3),
         helpText("Please indicate the number of cycles/columns in the file."),
         numericInput('fluo.cols', 'Number of columns', 40, min = 0, max = 100, step = 1),
         br(),
         downloadButton("downfluo","Download sample raw fluorescence data"),
         br(),
         helpText("Please note that if you do provide a simple efficiency file, the well numbers will be used as the sample names - if they do not match between the files, this will cause errors and the analysis will fail.")
         ),
        tabPanel("Results", tableOutput("tabelka")),
        tabPanel("Plots", uiOutput("patplot"))
        )
      

    )
  )
))