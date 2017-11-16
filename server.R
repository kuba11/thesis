library(shiny)
library('xlsx')
library('ggplot2')

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
        error = function(e) {print('Error while listing samples from the experimental data. Please check the input files.')}
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
        error = function(e) {print('Error while listing samples from the configuration file, please make sure it has the proper structure.')}
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




      boxplot(y ~ x, outline = F, main = c("Normalized and callibrated realive expression value"), col = heat.colors(length(unique(x))), xlab = "Samples")
      # plot(length(dim(y)), 1)
      #!!!wydajnosci

      
      })
    
    
    ### Barplot for all samples
    plot_output_3 <- renderPlot({
      

      y <- Fd()[, -c(1, dim(Fd())[1])]
      y2 <- as.numeric(c(unlist(t(y))))
      y2[y2>5] <- NA
      y3 <- matrix(y2, nrow = dim(y)[1], ncol = dim(y)[2], byrow = T)
      m <- apply(y3, 1, sd, na.rm = T)
      y4 <- apply(y3, 1, mean, na.rm = T)
      x <- data.matrix(Fd()[, 1])
      
      d <- data.frame(x, y = y4)
      f=ggplot(d, aes(x=x,y = m)) + geom_bar(stat = "identity",  fill = heat.colors(length(x)-1))
      f+geom_errorbar(aes(ymax=m+(y4)/2, ymin=m-(y4)/2), position="dodge")+
        ggtitle("Normalized and callibrated realive expression value") + theme(plot.title = element_text(lineheight=.8, face="bold")) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Samples")
      
      
    })
    
    
    
    plot_output_list <- c(plot_output_2, plot_output_3, plot_output_1)
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
        paste('sample_data','zip', sep='.')
      },
      content = function(file){
        file.copy("www/sample_data.zip", file)
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
