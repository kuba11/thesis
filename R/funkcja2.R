funkcja2 <- function (x, y, z, control){
  
  ## Removing '.txt' from the file names
  inFile <- x
  file1.name <- as.matrix(data.frame(strsplit(inFile$name, '[.]'))[1, ])
  
  ### Creating a list of matrices for each file

  a = list()
  for (i in 1:length(x[, 1])){ 
    a.temp <- read.table(inFile$datapath[i], sep = '\t', col.names = rep('a', 29), fill = T)
    a.temp <- a.temp[which(a.temp[, 1] == "Well"): which(a.temp[, 1] == "Slope"), 1:6]
    colnames(a.temp) <- as.character(unlist(a.temp[1,]))
    a[[i]] <- a.temp[-c(1, length(a.temp[, 1])), ]
  }
  
  ### Checking if the file is proper (has 4 sheets)
  inFile2 <- y
  plates <- try(read.xlsx(inFile2$datapath, header = T, 4))
  if("try-error" %in% class(plates)) {
      return ("This isn't a proper configuration file.")
    }else{
  
  ### Getting the efficiency from the file
  efficiency <- read.xlsx(inFile2$datapath, header = T, 1)
  eff <- efficiency[efficiency[, 1] %in% as.matrix(data.frame(strsplit(file1.name, '_'))[1, ]), ]
  
  ### Creating a matrix with all sample names for all sets based on the config file
  # information which gene uses which set 
  set <- data.frame(read.xlsx(inFile2$datapath, header = T, 3))
  # Checking which genes are loaded and which sets to use
  sets <- unique(set[efficiency[, 1] %in% as.matrix(data.frame(strsplit(file1.name, '_'))[1, ]), 2])
  # List of all samples that appear in the sets (that are corresponding to the genes loaded)

  uniques <- unique(unlist(subset(plates, select=as.character(sets))))
  # Sorting sample names
  uniques <- sort(as.character(uniques[is.finite(uniques)]))
  
  
  ### Reference gene selection
  k = 1  
  ref <- list()
  j <- 1
  liczba_usuniec <- 0
  ref_index <- c()
  
  # Creating an index for normal genes
  gene_index <- c(1:length(a))
  a.name <- lapply(file1.name, strsplit, '_')
  # Choosing first element (gene name) from the lists
  a.name <- unlist(lapply(a.name, function(l) l[[1]][1])) 
  ref_index <- which(a.name %in% z)
  ref_index <- ref_index[is.finite(ref_index)]
  gene_index <- gene_index[-ref_index]
  ref <- a[ref_index]
  a[ref_index] <- NULL
  
  
  ### Checking in gene names from the configuration file match the file names
  name_match <- match(a.name, eff[, 1])
  eff <- eff[name_match, ]
  eff <- eff[, 2]

  ### Checking, if normal and reference genes are present
  if (length(ref) > 0 & length(a) > 0){

  ### Reference genes calculations

    ## dCt for reference genes

    sample_no <- length(uniques)

    # Ct values for reference genes are matched to corresponding sample names

    refCt <- data.frame(matrix(nrow = length(ref), ncol = (sample_no - length(control))))
    refdCt <- refCt
    preQ <- refCt
    refnames <- refCt

    test_index <- which(uniques %in% unique(ref[[1]][, 2]) & !uniques %in% control) # Getting column numbers to consider based on the first gene
    
    for (i in 1:length(ref)){
      refCt[i, ] <- as.numeric(as.character(ref[[i]][test_index, 6]))
      refCal <- mean(as.numeric(as.character(ref[[i]][as.character(ref[[i]][, 2]) %in% control, 6])), na.rm = T)
      refdCt[i, ] <- refCal - refCt[i, ]
      refnames[i, ] <- ref[[i]][test_index, 2]
    }


    ## Calculating Q, efficiencies needed
    for (i in 1:length(ref)){
      preQ[i, ] <- eff[ref_index[i]]^refdCt[i, ]
    }


    refQ <- data.frame(matrix(nrow = 1, ncol = (sample_no - length(control))))
    colnames(refQ) <- uniques[!uniques %in% control]

    ### Geometric average. We create a list because we expect different number of data per sample.
    for (i in 1:length(uniques[!uniques %in% control])){
      x <- which(colnames(refQ)[i] ==  refnames, arr.ind = T)
      if(length(x) > 0)
      {
      ## If there is more than one nample of the same name
      refQ[i] <- (prod(preQ[x[,1],x[,2]], na.rm = T))^(1/sum(!is.na(preQ[x[,1],x[,2]])))
      }else{
        refQ[i] <- NA}
      }

    ### Normal genes calculations

    ## Calculating dCt
    Ct <- data.frame(matrix(nrow = length(a), ncol = (sample_no - length(control))))
    Q <- Ct
    samnames <- Ct
    Cal <- c()

    for (i in 1:length(a)){
      Ct[i, ] <- as.numeric(as.character(a[[i]][test_index, 6]))
      Cal[i] <- mean(as.numeric(as.character(a[[i]][as.character(a[[i]][, 2]) %in% control, 6])), na.rm = T)
      samnames[i, ] <- a[[i]][test_index, 2]
      }
    dCt <- Cal - Ct

    ## Calculating Q, efficiencies needed

    gene_index <- c(1:(length(a)+length(ref_index)))[-ref_index]
    for (i in 1:length(a)){
      Q[i, ] <- eff[gene_index[i]]^dCt[i, ]
    }



### Fold difference

Fd <- data.frame(matrix(ncol = length(gene_index), nrow = length(uniques[!uniques %in% control]))) #Tu bylo sample_no
rownames(Fd) <- uniques[!uniques %in% control]
x <- matrix()

for (i in 1:length(gene_index)){
  for (j in 1:length(uniques[!uniques %in% control])){

    ## We match data to corresponding names
    ## Each data matrix has corresponding matrix with their corresponding sample names
    x <- which(rownames(Fd)[j] == samnames[i, ], arr.ind = T)
    if (length(x) > 0){
    Fd[j, i] <- Q[i, x[, 2]]
    }else{
    Fd[j, i] <- NA
    }
    ###???!!!
    Fd[j, i] <- Fd[j, i]/refQ[which(rownames(Fd)[j] == colnames(refQ))] #?le, trzeba dopasowa? pr?bki z gen?w do ref.
  }
}

Fd <- cbind(uniques[!uniques %in% control], Fd)


## Row and column names
  k <- 1
  columns <- c()

  for (i in gene_index){
    columns[k] <- strsplit(file1.name[i], '_')[[1]][1]
    k <- k + 1

    Fd[is.na(Fd)] <- c('No data')
  }
  # Adding geometric average vals (ref genes) to the results
  Fd <- cbind(Fd, t(refQ))
  colnames(Fd) <- c('Sample Name', columns, "geNorm")


  } else if ( length(ref) == 0) {
    Fd <- c('No reference genes')
  } else{
    Fd <- c('No genes for the analysys (all genes are reference genes)')
  }

  if ((NA %in% match(z, a.name)) == T) {
    Fd <- c("Not all of the reference gene files are present!")
  }
  

  
  return(Fd)
    }    
}


# Wczytanie danych
# Zbada? to NA w match()

# Genenames[i, ] - sprawdzic ilocs probek dla kazdego !!! +
# Trzeba dopasowac probki z genow do ref. Trzeba to zrobic w petli, gdzie sie liczy srednia geometryczna

#samnames is the same size as Q, and hac names in the cells corresponding to each Q value
#Czy zmienily� sie wybór wydajności na zły?
