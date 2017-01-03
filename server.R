library(shiny)
library('xlsx')

shinyServer(function(input, output) {
source(paste(getwd(), 'funkcja1.R', sep = "/"))


  #Wybor genow referencyjnych
  output$refgen<-renderUI({
    
    if (is.null(input$file2))
      return(NULL)
    
    inFile2 <- input$file2
    selectizeInput("refergene",
  label="Wybierz geny referencyjne", choices=read.xlsx(inFile2$datapath, header = T, 2), multiple=T)
    
    })
  
  
  #Obliczanie Fold difference
  Fd <- reactive({

    Fd <- funkcja1(input$file1, input$file2, input$refergene)
  
    Fd})

  #Macierz wynikowa
  output$tabelka <- renderTable({
    
    if (is.null(input$file1) | is.null(input$file2))
      return(NULL)
    
    Fd()
    
  })

  
  output$down <- downloadHandler(
    filename=function() {
      paste('wynik','.xlsx', sep='')
    },
    content = function(file){
      write.xlsx(Fd(), file, row.names = F)
      
    })

})
