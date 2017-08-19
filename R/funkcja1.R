funkcja1 <- function (x, y, z){
    
  #Tworzymy macierz, gdzie ka¿da cz. listy to 1 plik
  inFile <- x
  a = list()
  for (i in 1:length(x[, 1])){ 
    a.temp <- read.table(inFile$datapath[i], sep = '\t', col.names = rep('a', 29), fill = T)
    a.temp <- a.temp[which(a.temp[, 1] == "Well"): which(a.temp[, 1] == "Slope"), 1:6]
    colnames(a.temp) <- as.character(unlist(a.temp[1,]))
    a[[i]] <- a.temp[-c(1, (length(a.temp[, 1])-2):length(a.temp[, 1])), ]
  }
  
  #Wydajnoœæ
  inFile2 <- y
  efficiency <- read.xlsx(inFile2$datapath, header = T, 1)
  eff <- efficiency[efficiency[, 1] %in% as.matrix(data.frame(strsplit(inFile$name, '_'))[1, ]), ]###!!!!!!!!!!!!!!!!!!!
  
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
    Fd <- c("No reference genes. The name of the reference gene chosen doesn't match the first part of any filename.")
  } else{
    Fd <- c('No genes for the analysys (all genes are reference genes)')


}

Fd <- as.data.frame(Fd)

  return(Fd)}
#eff!!!