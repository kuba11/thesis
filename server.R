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


#Boxplots for samples
    output$patplot <- renderUI({
      
      if (is.null(input$file1) | is.null(input$file2) | is.null(input$refergene) | input$act == 0)
        return(NULL)
      
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
