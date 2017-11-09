library(shiny)
library('xlsx')

shinyServer(function(input, output) {
  source(paste(getwd(), 'R/funkcja1.R', sep = "/"))
  source(paste(getwd(), 'R/funkcja2.R', sep = "/"))


  ## Choosing the reference genes
  output$refgen<-renderUI({
    
    if (is.null(input$file2))
      return(NULL)
    
    inFile2 <- input$file2
    
    tryCatch(
    selectizeInput("refergene",
  label="Reference genes", choices=read.xlsx(inFile2$datapath, header = T, 1)[, 1], multiple=T),
  error = function(e) {print('Not a proper file with efficiencies:no gene names!')}
    )

    })
  
  ### Choosing the control samples
  output$control<-renderUI({
    
    if (is.null(input$file1) | is.null(input$file2))
      return(NULL)
    
    ## Reading the list of samples from the first loaded file (in case of efficiency file)
    if(input$eff.file == T){
      
      ## List of samples taken from the first file
      sample_list <- read.table(input$file1$datapath[1], sep = '\t', col.names = rep('a', 29), fill = T)
      sample_list <- sample_list[which(sample_list[, 1] == "Well"): which(sample_list[, 1] == "Slope"), 2]
      sample_list <- sample_list[-1]
      
      # The dropdown selection menu
      tryCatch(
        selectizeInput("control",
                       label="Control samples", choices = sample_list, multiple=T),
        error = function(e) {print('Not a proper file with efficiencies!')}
      )
    }else{
      
      ## List of samples
      # Information about sample placement for each set
      plates <- read.xlsx(input$file2$datapath, header = T, 4)
      # information which gene uses which set 
      set <- data.frame(read.xlsx(input$file2$datapath, header = T, 3))
      # Getting the list of all gene names
      names <- as.character(unlist(read.xlsx(input$file2$datapath, header = T, 2)))
      # Checking which genes are loaded and which sets to use
      sets <- unique(set[names %in% as.matrix(data.frame(strsplit(input$file1$name, '_'))[1, ]), 2])
      # List of all samples that appear in the sets (that are corresponding to the genes loaded)
      uniques <- unique(unlist(subset(plates, select=as.character(sets))))
      # Sortling sample names
      sample_list <- sort(as.character(uniques[is.finite(uniques)]))
      
      x <- match(names, as.matrix(data.frame(strsplit(input$file1$name, '_'))[1, ]))
      
      tryCatch(
        selectizeInput("control",
                       label="Control samples", choices = sample_list, multiple=T),
        error = function(e) {print('Not a proper file with efficiencies:2!')}
      )
    }
    

    
  })
  
  
  #Obliczanie Fold difference
  Fd <- reactive({
    
    if(!is.null(input$control)){
    if (input$eff.file == T){
      
    Fd <- funkcja1(input$file1, input$file2, input$refergene, input$control)
    
    }else{
      
    Fd <- funkcja2(input$file1, input$file2, input$refergene, input$control)
    
    }
    }else{Fd <- c("No control samples were chosen!")}
    
    Fd
    })

  #Macierz wynikowa
  output$tabelka <- renderTable({
    
    if (is.null(input$file1) | is.null(input$file2) | is.null(input$refergene) | input$act == 0)
      return(NULL)
    
    Fd()
    
  })


#Boxplots for samples
    output$patplot <- renderUI({
      
      if (is.null(input$file1) | is.null(input$file2) | is.null(input$refergene) |  is.null(input$control) | input$act == 0 | length(dim(Fd())) == 0)
        return(NULL)
    
      
      
      
### Single boxplots
    plot_output_1 <- lapply(1:dim(Fd())[1], function(i) {
  output[[paste0('b', i)]] <- renderPlot({

    ## Removing first (sample name) and last (geNorm) columns
    y <- as.numeric(as.character(data.frame(lapply(Fd()[i, -c(1, dim(Fd())[1])], as.character), stringsAsFactors=FALSE)))
    x <- data.matrix(Fd()[i, 1])
    x <- rep(x, each = dim(Fd())[2]-1)
    
    
      try(boxplot(y ~ x, outline = F, main = paste(c('Sample'), x[1], sep = ' ')), silent = T)

  })
})
    ### First plot for all samples
    plot_output_2 <- renderPlot({
      
      y <- Fd()[, -c(1, dim(Fd())[1])]
      y <- as.numeric(c(unlist(t(y))))
      x <- data.matrix(Fd()[, 1])
      x <- rep(x, each = dim(Fd())[2]-1)
      
      
      
      
      boxplot(y ~ x, outline = F, main = c("Boxplot for all samples"))
      # plot(length(dim(y)), 1)
      #!!!wydajnosci
    })
    
    
    plot_output_list <- c(plot_output_2, plot_output_1)
    do.call(tagList, plot_output_list)

   })

    

    ### Results table
    output$down <- downloadHandler(
      filename=function() {
      paste('wynik','.xlsx', sep='')
    },
    content = function(file){
      write.xlsx(Fd(), file, row.names = F)
    })
    
    ### Sample files
    ## Sample data
    output$downdata <- downloadHandler(
      filename = function() {
        paste('ACTB_sample','txt', sep='.')
      },
      content = function(file){
        file.copy("www/ACTB_sample.txt", file)
      })
    
    ## Sample efficiency
    output$downeff <- downloadHandler(
      filename=function() {
        paste('sample_efficiency','.xlsx', sep='')
      },
      content = function(file){
        file.copy("www/efficiency.xlsx", file)
      })
    
    ## Empty efficiency
    output$downeffe <- downloadHandler(
      filename=function() {
        paste('empty_efficiency','.xlsx', sep='')
      },
      content = function(file){
        file.copy("www/empty_efficiency.xlsx", file)
      })
    
    ## Sample cofniguration
    output$downconf <- downloadHandler(
      filename=function() {
        paste('sample_config','.xlsx', sep='')
      },
      content = function(file){
        file.copy("www/config file.xlsx", file)
      })
    
    ## Empty cofniguration
    output$downconfe <- downloadHandler(
      filename=function() {
        paste('empty_config','.xlsx', sep='')
      },
      content = function(file){
        file.copy("www/empty_config file.xlsx", file)
      })
    
    
    
})
