library(shiny)
library('xlsx')
library('ggplot2')
library(qpcR)

shinyServer(function(input, output) {
  source(paste(getwd(), 'R/funkcja1.R', sep = "/"))
  source(paste(getwd(), 'R/funkcja2.R', sep = "/"))
  source(paste(getwd(), 'R/fluorescence1.R', sep = "/"))
  source(paste(getwd(), 'R/fluorescence2.R', sep = "/"))

####### 1. DATA INPUT ---------------------------------------------------------------------------- 
  ## Choosing the reference genes (step 2)
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
  
  ### Choosing the control samples (step 4)
  output$control<-renderUI({
    
    if (is.null(input$file1) | is.null(input$file2))
      return(NULL)
    

    ## Reading the list of samples from the first loaded file (in case of efficiency file)
    if(input$eff.file == T){
      
      ## List of samples taken from the first file
      if(input$fluo.file == F){
        sample_list <- read.table(input$file1$datapath[1], sep = '\t', col.names = rep('a', 29), fill = T)
        sample_list <- sample_list[which(sample_list[, 1] == "Well"): which(sample_list[, 1] == "Slope"), 2]
        sample_list <- sample_list[-1]
      }else{
        sample_list <- read.table(input$file1$datapath[1], sep = '\t', header = F, skip = 2)
        sample_list <- sample_list[, 1]
        
      }

      
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
  
  
####### 2. CALCULATIONS------------------------------------------------------------------------
  
  Fd <- reactive({
    
    if(!is.null(input$control)){
    if (input$eff.file == T){
      
      # Inputs: Input files(raw or ct), configuration file, reference gene names, control sample names, is it raw fluorescence data? [T/F]
    Fd <- funkcja1(input$file1, input$file2, input$refergene, input$control, input$fluo.file, input$fluo.thresh, input$fluo.cols)
    
    }else{
        
        Fd <- funkcja2(input$file1, input$file2, input$refergene, input$control, input$fluo.file, input$fluo.thresh, input$fluo.cols)
      }
    }else{Fd <- c("No control samples were chosen!")}
    
    Fd
    })

  #Output matrix
  output$tabelka <- renderTable({
    
    if(input$act == 0){
      mes <- data.frame("Please upload the files, make relevant selections and press the 'ANALYSE' button")
      colnames(mes) <- 'Information'
      return(mes)
    }else {

      if (is.null(input$file1) | is.null(input$file2) | is.null(input$refergene) )
        return('One of the input files or information is missing!')

      Fd()
    }
  })


####### 3. GRAPHICAL RESULTS -------------------------------------------------------------
  
### Plots for samples
    output$patplot <- renderUI({
      
      if (is.null(input$file1) | is.null(input$file2) | is.null(input$refergene) |  is.null(input$control) | input$act == 0 | length(dim(Fd())) == 0)
        return(NULL)

  ##### 3.1 Single boxplots
    plot_output_1 <- lapply(1:dim(Fd())[1], function(i) {
  output[[paste0('b', i)]] <- renderPlot({

    ## Removing first (sample name) and last (geNorm) columns
    y <- as.numeric(as.character(data.frame(lapply(Fd()[i, -c(1, dim(Fd())[1])], as.character), stringsAsFactors=FALSE)))
    x <- data.matrix(Fd()[i, 1])
    x <- rep(x, each = dim(Fd())[2]-1)
    
    
      try(boxplot(y ~ x, outline = F, main = paste(c('Sample'), x[1], sep = ' '), col = heat.colors(dim(Fd())[1])[i], xlab = 'Samples', ylab = 'Fold difference'), silent = T)

  })
})
  ##### 3.2 First plot for all samples
    plot_output_2 <- renderPlot({
      
      y <- Fd()[, -c(1, dim(Fd())[2])]
      y <- as.numeric(c(unlist(t(y))))
      x <- data.matrix(Fd()[, 1])
      x <- rep(x, each = dim(Fd())[2]-2)



      # showNotification(paste('1', length(y)), duration = NULL)
      boxplot(y ~ x, outline = F, main = c("Boxplot - Normalized and callibrated relative expression value"), col = heat.colors(length(unique(x))), xlab = "Samples", ylab = 'Fold difference',  las = 2)
      # plot(length(dim(y)), 1)
      #!!!wydajnosci

      
      })
    
    
  ##### 3.3 Barplot for all samples
    plot_output_3 <- renderPlot({
      

      y <- Fd()[, -c(1, dim(Fd())[1])]
      y2 <- as.numeric(c(unlist(t(y))))
      y3 <- matrix(y2, nrow = dim(y)[1], ncol = dim(y)[2], byrow = T)
      sd <- apply(y3, 1, sd, na.rm = T)
      mean <- apply(y3, 1, mean, na.rm = T)
      sample_names <- data.matrix(Fd()[, 1])
      #showNotification(paste('', mean[18], mean[19], mean[20], sep = ' '), duration = NULL)

      names(mean) <- sample_names#[is.finite(mean)]
      mean <- mean[is.finite(mean)]
      
      barCenters <- barplot(mean, col = heat.colors(length(sample_names)), las = 2, ylim = c(0, max(mean + sd/2, na.rm = T)),
                            main = "Barplot - Normalized and callibrated relative expression value", ylab = 'Fold difference')

      segments(barCenters, mean - sd/2, barCenters,
               mean + sd/2, lwd = 1.5)

      arrows(barCenters, mean - sd/2, barCenters,
             mean + sd/2, lwd = 1.5, angle = 90,
             code = 3, length = 0.05)
    })
    
    
    ### Merging all plots into one object to display on the page
    
    plot_output_list <- c(plot_output_2, plot_output_3, plot_output_1)
    do.call(tagList, plot_output_list)

   })
    
    

####### 4. RESULTS TABLE -------------------------------------------------------------
  
    output$down <- downloadHandler(
      filename=function() {
      paste('Results','.xlsx', sep='')
    },
    content = function(file){
      write.xlsx(Fd(), file, row.names = F)
    })
    
####### 5. SAMPLE FILES ---------------------------------------------------------------------
   
    ## Sample data (Ct)
    output$downdata <- downloadHandler(
      filename = function() {
        paste('sample_data','zip', sep='.')
      },
      content = function(file){
        file.copy("www/sample_data.zip", file)
      })
  
  ## Sample fluorescence data
  output$downfluo <- downloadHandler(
    filename = function() {
      paste('fluo_data','zip', sep='.')
    },
    content = function(file){
      file.copy("www/fluo_data.zip", file)
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
