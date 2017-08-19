funkcja2 <- function (x, y, z){
  
  # Tworzymy macierz, gdzie ka?da cz. listy to 1 plik
  inFile <- x
  a = list()
  for (i in 1:length(x[, 1])){ 
    a.temp <- read.table(inFile$datapath[i], sep = '\t', col.names = rep('a', 29), fill = T)
    a.temp <- a.temp[which(a.temp[, 1] == "Well"): which(a.temp[, 1] == "Slope"), 1:6]
    colnames(a.temp) <- as.character(unlist(a.temp[1,]))
    a[[i]] <- a.temp[-c(1, (length(a.temp[, 1])-2):length(a.temp[, 1])), ]
  }
  
  # Wydajno??
  inFile2 <- y
  plates <- read.xlsx(inFile2$datapath, header = T, 4)
  
  efficiency <- read.xlsx(inFile2$datapath, header = T, 1)
  eff <- efficiency[efficiency[, 1] %in% as.matrix(data.frame(strsplit(inFile$name, '_'))[1, ]), ]
  
  # Creating a  matrix with all sample names for all sets based on the config file
  
  set <- data.frame(read.xlsx(inFile2$datapath, header = T, 3))
  sets <- unique(set[efficiency[, 1] %in% as.matrix(data.frame(strsplit(inFile$name, '_'))[1, ]), 2])
  uniques <- unique(unlist(subset(plates, select=as.character(sets))))
  uniques <- sort(uniques[is.finite(uniques)])
  
  # Slecting reference genes

  #Selekcja genów referencyjnych
  
  k = 1  
  ref <- list()
  j <- 1
  liczba_usuniec <- 0
  ref_index <- c()
  
  
  gene_index <- c(1:length(a))
  a.name <- lapply(inFile$name, strsplit, '_')
  a.name <- unlist(lapply(a.name, function(l) l[[1]][1])) #Choosing firs element (gene name) from the lists
  ref_index <- which(z == a.name)
  ref_index <- ref_index[is.finite(ref_index)]
  gene_index <- gene_index[-ref_index]
  ref <- a[ref_index]
  a[ref_index] <- NULL
  
  
  ### Sprawdzanie, czy geny z pliku z wydajnoœciami pasuj¹ do nazw plików
  name_match <- match(a.name, eff[, 1])
  xxx = eff[,1]
  eff <- eff[name_match, ]
  eff <- eff[, 2]

  ### Sprawdzanie, czy s? geny referencyjne i zwyk?e
  if (length(ref) > 0 & length(a) > 0){

    ### Obliczenia dla gen?w referencyjnych
    #Obliczanie dCt dla ref

    sample_no <- length(ref[[1]][, 1])
    
    # Ct values for reference genes are matched to corresponding names
    
    refCt <- data.frame(matrix(nrow = length(ref), ncol = sample_no))
    refdCt <- refCt
    preQ <- refCt
    refnames <- refCt

    for (i in 1:length(ref)){
      refCt[i, ] <- as.numeric(as.character(ref[[i]][, 6]))
      refCal <- mean(as.numeric(as.character(ref[[i]][(sample_no-5):sample_no, 6])), na.rm = T)
      refdCt[i, ] <- refCal - refCt[i, ]
      refnames[i, ] <- ref[[i]][, 2]
    }

    
    # Obliczamy Q, czyli potrzebne s? wydajno?ci
    for (i in 1:length(ref)){
      preQ[i, ] <- eff[ref_index[i]]^refdCt[i, ]
    }

    
    refQ <- data.frame(matrix(nrow = 1, ncol = length(uniques)))
    colnames(refQ) <- uniques
    
    # Geometric average. We create a list because we expect different number of data per sample.
    for (i in 1:length(uniques)){

      x <- which(colnames(refQ)[i] ==  refnames, arr.ind = T)
      
      if(length(x) > 0)
      {

      refQ[i] <- (prod(preQ[x[,1],x[,2]], na.rm = T))^(1/sum(!is.na(preQ[x[,1],x[,2]])))

      }else{
        refQ[i] <- NA}

      }


    ##### Przygotowanie danych
    # Obliczanie dCt
    Ct <- data.frame(matrix(nrow = length(a), ncol = sample_no))
    Q <- Ct
    samnames <- Ct
    Cal <- c()

    for (i in 1:length(a)){
      Ct[i, ] <- as.numeric(as.character(a[[i]][, 6]))
      Cal[i] <- mean(as.numeric(as.character(a[[i]][(sample_no-5):sample_no, 6])), na.rm = T)
      samnames[i, ] <- a[[i]][, 2]
      }
    dCt <- Cal - Ct

    # Obliczamy Q, czyli potrzebne s? wydajno?ci

    gene_index <- c(1:(length(a)+length(ref_index)))[-ref_index]

    for (i in 1:length(a)){
      Q[i, ] <- eff[gene_index[i]]^dCt[i, ]
    }



    ### Fold difference

    Fd <- data.frame(matrix(ncol = length(gene_index), nrow = length(uniques))) #Tu bylo sample_no
    rownames(Fd) <- uniques
    x <- matrix()

    for (i in 1:length(gene_index)){
      for (j in 1:length(uniques)){

        # We match data to corresponding names
        # Each data matrix has corresponding matrix with their corresponding sample names
        x <- which(rownames(Fd)[j] == samnames[i, ], arr.ind = T)
        if (length(x) > 0){
        Fd[j, i] <- Q[i, x[, 2]]
        }else{
        Fd[j, i] <- NA
        }
        Fd[j, i] <- Fd[j, i]/refQ[which(rownames(Fd)[j] == colnames(refQ))] #?le, trzeba dopasowa? pr?bki z gen?w do ref.
      }
    }

    
    Fd <- cbind(uniques, Fd)
    
      #Nazwy wierszy i kolumn


      k <- 1
      columns <- c()

      for (i in gene_index){
        columns[k] <- strsplit(inFile$name[i], '_')[[1]][1]
        k <- k + 1

        Fd[is.na(Fd)] <- c('No data')
      }
      colnames(Fd) <- c('Sample Name', columns)


  } else if ( length(ref) == 0) {
    Fd <- c('No reference genes')
  } else{
    Fd <- c('No genes for the analysys (all genes are reference genes)')


  }

  return(Fd)
         
}

# Wczytanie danych
# Zapytac, czy o to chodzi i o dane z roznymi nazwami probek (raczej nie bedzie)
# Zbada? to NA w match()

# Genenames[i, ] - sprawdzic ilocs probek dla kazdego !!! +
# Trzeba dopasowac probki z genow do ref. Trzeba to zrobic w petli, gdzie sie liczy srednia geometryczna

#Nie dzila dla kilku genow ref!!!
#Wyniki sa w roznej kolejnosci

#samnames is the same size as Q, and hac names in the cells corresponding to each Q value
#Czy zmienily‚ sie wybÃ³r wydajnoÅ›ci na zÅ‚y?
