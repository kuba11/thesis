funkcja1 <- function (x, y, z, control, fluo, thresh = 1, cols = 1){
  
####### 1. File loading ---------------------------------------------------------------
  
  inFile <- x
  ## Removing '.txt' from the file names
  file1.name <- as.matrix(data.frame(strsplit(inFile$name, '[.]'))[1, ])
  
  ##### 1.1 Raw fluorescence files
  if(fluo == T){
    a <- fluorescence1(x, thresh, cols + 1)
    
  }else{
  ##### 1.2 Processed data (with Ct values)
  ### Creating a list of matrices for each file
  
  a = list()
  for (i in 1:length(inFile[, 1])){ 
    a.temp <- read.table(inFile$datapath[i], sep = '\t', col.names = rep('a', 29), fill = T)
    a.temp <- a.temp[which(a.temp[, 1] == "Well"): which(a.temp[, 1] == "Slope"), c(2, 6)]
    colnames(a.temp) <- as.character(unlist(a.temp[1,]))
    a[[i]] <- a.temp[-c(1, length(a.temp[, 1])), ]
  }
  }
  
  ### Getting the efficiency from the file
  inFile2 <- y
  efficiency <- read.xlsx(inFile2$datapath, header = T, 1)
  eff <- efficiency[efficiency[, 1] %in% as.matrix(data.frame(strsplit(file1.name, '_'))[1, ]), ]###!!!!!!!!!!!!!!!!!!!



####### 2. Separating samples into reference and test ---------------------------------
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


### Checking in gene names from the efficiency file match the file names
name_match <- match(a.name, eff[, 1])
eff <- eff[name_match, ]
eff <- eff[, 2]




####### 3. Calculation of Fold difference -------------------------------------------------------

### Checking, if normal and reference genes are present
  if (length(ref) > 0 & length(a) > 0){

  ##### 3.1 Reference genes calculations
  ## dCt for reference genes
  sample_no <- length(unique(ref[[1]][, 2]))
  refCt <- data.frame(matrix(nrow = length(ref), ncol = (sample_no - length(control))))
  refdCt <- refCt
  preQ <- refCt
  refQ <- data.frame(matrix(nrow = 1, ncol = (sample_no - length(control))))

  for (i in 1:length(ref)){
    refCt[i, ] <- as.numeric(as.character(ref[[i]][!as.character(ref[[i]][, 1]) %in% control, 2]))
    refCal <- mean(as.numeric(as.character(ref[[i]][as.character(ref[[i]][, 1]) %in% control, 2])), na.rm = T)
    refdCt[i, ] <- refCal - refCt[i, ]
  }

  ## Calculating Q, efficiencies needed
  for (i in 1:length(ref)){
    preQ[i, ] <- eff[ref_index[i]]^refdCt[i, ]
  }


  ## Geometric average
  for (i in 1:dim(preQ)[2]){
    refQ[i] <- (prod(preQ[, i], na.rm = T))^(1/sum(!is.na(preQ[, i])))
  }


  ##### 3.2 Normal genes calculations

  ## Calculating dCt
  Ct <- data.frame(matrix(nrow = length(a), ncol = (sample_no - length(control))))
  Q <- Ct
  Cal <- c()
  xx <- list()

  for (i in 1:length(a)){
    Ct[i, ] <- as.numeric(as.character(a[[i]][!as.character(a[[i]][, 1]) %in% control, 2]))
    Cal[i] <- mean(as.numeric(as.character(a[[i]][as.character(a[[i]][, 1]) %in% control, 2])), na.rm = T)
  #  xx[[i]] <- as.character(a[[i]][, 1])
  }
    dCt <- Cal - Ct

  ## Calculating Q, efficiencies needed
  for (i in 1:length(a)){
    Q[i, ] <- eff[gene_index[i]]^dCt[i, ]
  }


  ##### 3.3 Fold difference

Fd <- data.frame(matrix(nrow = length(gene_index), ncol = (sample_no - length(control))))

for (i in 1:length(gene_index)){
  Fd[i, ] <- Q[i, ]/refQ
}
Fd <- t(Fd)

### Checking if the number of samples is the same for normal and reference genes
if (dim(refQ)[2] != dim(Q)[2]){
  return(NULL)
}else{


    ## Row and column names
    Fd <- cbind(as.character(a[[1]][!as.character(a[[1]][, 1]) %in% control, 1]), Fd)
    k <- 1
    columns <- c()

    for (i in gene_index){
     columns[k] <- strsplit(file1.name[i], '_')[[1]][1]
    k <- k + 1

    Fd[is.na(Fd)] <- c('No data')
    }

    # Adding reference genes' geometric mean values
    Fd <- cbind(Fd, t(refQ))

    colnames(Fd) <- c('Sample Name', columns, "geNorm")
}

  } else if ( length(ref) == 0) {
    Fd <- c("No reference genes. The name of the reference gene chosen doesn't match the first part of any filename.")
  } else {
    Fd <- c('No genes for the analysys (all genes are reference genes)')
  }



### Checking if all reference gene files are present
if ((NA %in% match(z, a.name)) == T) {
  Fd <- c("Not all of the reference gene files are present!")
}



  return(Fd)
}
#eff!!!