library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  titlePanel("Uploading Files"),
  sidebarLayout(
    sidebarPanel(
      helpText("Krok 3: Prosze wybrac dane z eksperymentu RT-PCR"),
      fileInput("file1", label = "Wybierz pliki", multiple =
                  T),
      
      br(),
      p("Kliknij w ponizszy przycisk jesli chcesz pobrac uzyskane wyniki analizy"),
      downloadButton("down","Pobierz wyniki")
      
      ),
    
    
    mainPanel(
      helpText("Krok 1: Prosze wybrac plik konfiguracyjny"),
      fileInput("file2", label = "choose file", multiple =
                  T),
      br(),
      
      helpText("Krok 2: Prosze wybrac dowolna ilosc genow referencyjnych z listy"),
      uiOutput("refgen"),
      
      tableOutput('tabelka')
    )
  )
))