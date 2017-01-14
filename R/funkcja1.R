funkcja1 <- function (x, y, z){
  
  #Tworzymy macierz, gdzie ka¿da cz. listy to 1 plik
  inFile <- x
  a = list()
  for (i in 1:length(x[, 1])){ 
  a[[i]] <- read.table(textConnection(rev(rev(readLines(inFile$datapath[i]))[-(1:10)])),
                              skip = 10, header = T, sep = "\t", nrows = 68)
  }
  
  #Wydajnoœæ
  inFile2 <- y
  eff <- read.xlsx(inFile2$datapath, header = T, 1)[, 2]
  
  #Selekcja genów referencyjnych
  
  k = 1  
  ref <- list()
  j <- 1
  liczba_usuniec <- 0
  ref_index <- c()
  
  
  for (i in 1:length(z)){
    for (j in 1:length(a)){
  if (strsplit(inFile$name[j], '_')[[1]][1] == z[i]){
    #Przypisywanie plików z genami ref do nowej zmiennej
    ref[[k]] = a[[j]]
    


  k=k+1
  #Usuwanie z listy normalnych genów
  
  
  ref_index[k - 1] <- j 
  }}}
  
  k <- 1
for (i in ref_index){
  j <- i + 1 - k
  a[[j]] <- NULL
  k <- k + 1
  }


### Sprawdzanie, czy s¹ geny referencyjne i zwyk³e 
  if (length(ref) > 0 & length(a) > 0){
  
### Obliczenia dla genów referencyjnych
  #Obliczanie dCt dla ref
    
  sample_no <- length(ref[[1]][, 1])
  refCt <- data.frame(matrix(nrow = length(ref), ncol = sample_no))
  refdCt <- refCt
  preQ <- refCt
  refQ <- data.frame(matrix(nrow = 1, ncol = sample_no))
  
  for (i in 1:length(ref)){
    refCt[i, ] <- as.numeric(as.character(ref[[i]][, 6]))
    refCal <- mean(as.numeric(as.character(ref[[i]][(sample_no-5):sample_no, 6])), na.rm = T)
    refdCt[i, ] <- refCal - refCt[i, ]
  }
  
  # Obliczamy Q, czyli potrzebne s¹ wydajnoœci
  for (i in 1:length(ref)){
    preQ[i, ] <- eff[ref_index[i]]^refdCt[i, ]
  }
  
  
  #Œrednia geometryczna
  for (i in 1:dim(preQ)[2]){
    refQ[i] <- (prod(preQ[, i], na.rm = T))^(1/sum(!is.na(preQ[, i]))) 
  }

  
### Przygotowanie danych
  # Obliczanie dCt
  Ct <- data.frame(matrix(nrow = length(a), ncol = sample_no))
  Q <- Ct
  Cal <- c()

  for (i in 1:length(a)){
    Ct[i, ] <- as.numeric(as.character(a[[i]][, 6]))
    Cal[i] <- mean(as.numeric(as.character(a[[i]][(sample_no-5):sample_no, 6])), na.rm = T)
  }
    dCt <- Cal - Ct
  
  # Obliczamy Q, czyli potrzebne s¹ wydajnoœci
  
  gene_index <- c(1:(length(a)+length(ref_index)))[-ref_index]
  
  for (i in 1:length(a)){
    Q[i, ] <- eff[gene_index[i]]^dCt[i, ]
  }
  

### Fold difference

Fd <- data.frame(matrix(nrow = length(gene_index), ncol = sample_no))

for (i in 1:length(gene_index)){
  Fd[i, ] <- Q[i, ]/refQ
}
Fd <- t(Fd)

### Sprawdzamy, czy liczba próbek jest taka sama
if (dim(refQ)[2] != dim(Q)[2]){
  return(NULL)
}else{


    #Nazwy wierszy i kolumn
    
    Fd <- cbind(as.character(a[[1]][, 2]), Fd)
    k <- 1 
    columns <- c()
    
    for (i in gene_index){
     columns[k] <- strsplit(inFile$name[i], '_')[[1]][1]
    k <- k + 1
    
    Fd[is.na(Fd)] <- c('No data')
}
    colnames(Fd) <- c('Sample Name', columns)
}

  } else if ( length(ref) == 0) {
    Fd <- c('No reference genes')
  } else{
    Fd <- c('No genes for the analysys (all genes are reference genes)')
  
    
    }
  return(Fd)
}