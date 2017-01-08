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
    
  output$genplot <- renderUI({
      
      if (is.null(input$file1) | is.null(input$file2) | is.null(input$refergene) | input$act == 0)
        {return(NULL)}
      
      #Boxplot zbiorczy
    plota <- renderPlot({
      y <- matrix(ncol = (dim(Fd())[2]-1), nrow = dim(Fd())[1])
      for (i in 2:(dim(Fd())[2])){
      y[, i-1] <- as.numeric(as.character(data.frame(lapply(Fd()[, i], as.character), stringsAsFactors=FALSE)))
      }
    
      x <- colnames(Fd())
      x <- x[-c(1)]
    colnames(y) <- x
    x <- rep(x, each = dim(Fd())[1])
    x <- matrix(x, ncol = dim(Fd())[2]-1)
      boxplot (y ~ x, outline = F, main = 'Boxplots for all genes')
      
    })
      
      #Boxploty dla genow

      plotb <- lapply(1:(dim(Fd())[2] -1), function(i) {
        output[[paste0('b', i)]] <- renderPlot({
          
          
          y <- as.numeric(as.character(data.frame(lapply(Fd()[, i+1], as.character), stringsAsFactors=FALSE)))
          
          x <- colnames(Fd())
          x <- x[-c(1)]
          x<- x[i]
          x <- rep(x, each = dim(Fd())[1])
          x <- matrix(x, ncol = 1)
          
          boxplot(y ~ x, outline = F, main = x[1])
          
        })
      })
      
      plot_output_list <- list(plota, plotb)
      do.call(tagList, plot_output_list)
    })
    
 
    


    output$patplot <- renderUI({
    plot_output_list <- lapply(1:dim(Fd())[1], function(i) {
  output[[paste0('b', i)]] <- renderPlot({
    
    
            y <- as.numeric(as.character(data.frame(lapply(Fd()[i, -c(1)], as.character), stringsAsFactors=FALSE)))
            
    x <- data.matrix(Fd()[i, 1])
    x <- rep(x, each = dim(Fd())[2]-1)
    
    
    boxplot(y ~ x, outline = F, main = paste(c('Sample'), x[1], sep = ' '))
    
  })
})
    do.call(tagList, plot_output_list)
    })

    

    
    output$down <- downloadHandler(
      filename=function() {
      paste('wynik','.xlsx', sep='')
    },
    content = function(file){
      write.xlsx(Fd(), file, row.names = F)
      
    })

})
