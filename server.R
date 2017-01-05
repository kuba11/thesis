library(shiny)
library('xlsx')

shinyServer(function(input, output) {
  source(paste(getwd(), 'R/funkcja1.R', sep = "/"))


  #Wybor genow referencyjnych
  output$refgen<-renderUI({
    
    if (is.null(input$file2))
      return(NULL)
    
    inFile2 <- input$file2
    selectizeInput("refergene",
  label="Reference genes", choices=read.xlsx(inFile2$datapath, header = T, 1), multiple=T)
    
    })
  
  
  #Obliczanie Fold difference
  Fd <- reactive({

    Fd <- funkcja1(input$file1, input$refergene)
  
    Fd})

  #Macierz wynikowa
  output$tabelka <- renderTable({
    
    if (is.null(input$file1) | is.null(input$file2) | is.null(input$refergene) | input$act == 0)
      return(NULL)
    
    Fd()
    
  })

  #Wykresy pudelkowe - geny
    
    output$genplot <- renderPlot({
      
      if (is.null(input$file1) | is.null(input$file2) | is.null(input$refergene) | input$act == 0)
        {return(NULL)}
      
      y <- matrix(ncol = (dim(Fd())[2]-1), nrow = dim(Fd())[1]  )
      for (i in 2:(dim(Fd())[2])){
      y[, i-1] <- as.numeric(data.matrix(Fd()[, i]))
      }
      
      x <- colnames(Fd())
      x <- x[-c(1)]
    colnames(y) <- x
    x <- rep(x, each = dim(Fd())[1])
    x <- matrix(x, ncol = dim(Fd())[2]-1)
      
      boxplot (y ~ x, outline = F)
    })

    #Wykresy pudelkowe - pacjenci
    
    output$patplot <- renderPlot({
      
      if (is.null(input$file1) | is.null(input$file2) | is.null(input$refergene) | input$act == 0)
      {return(NULL)}
      
      y <- matrix(ncol = (dim(Fd())[2]-1), nrow = dim(Fd())[1]  )
      for (i in 2:(dim(Fd())[2])){
        y[, i-1] <- as.numeric(data.matrix(Fd()[, i]))
      }
      x <- data.matrix(Fd()[, 1])
      x <- rep(x, each = dim(Fd())[2]-1)
      x <- matrix(x, ncol = dim(Fd())[2]-1, byrow = T)
      

      boxplot(y ~ x, outline = F)

      
    })
  
  output$down <- downloadHandler(
    filename=function() {
      paste('wynik','.xlsx', sep='')
    },
    content = function(file){
      write.xlsx(Fd(), file, row.names = F)
      
    })

})
