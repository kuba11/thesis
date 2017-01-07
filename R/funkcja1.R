funkcja1 <- function (x, z){
  
  #Tworzymy macierz, gdzie ka¿da cz. listy to 1 plik
  inFile <- x
  a = list()
  for (i in 1:length(x[, 1])){ 
  a[[i]] <- matrix(read.xlsx(inFile$datapath[i], header = F, 1))
  }
  
  #Œrednia z podwójnych prób
  
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
  refCt <- matrix(nrow = length(ref), ncol = (length(ref[[1]][[1]]))-1)
  refdCt <- refCt
  preQ <- refCt
  refQ <- matrix(nrow = 1, ncol = (length(ref[[1]][[1]]))-1)
  
  for (i in 1:length(ref)){
    refCt[i, ] <- c(as.numeric(sub(",", ".", (ref[[i]][[2]][-1]), fixed=TRUE)))
    refCal <- max(ref[[i]][[4]][!is.na(ref[[i]][[ 4]])])
    refdCt[i, ] <- refCal - refCt[i, ]
  }
  
  # Obliczamy Q, czyli potrzebne s¹ wydajnoœci
  for (i in 1:length(ref)){
    preQ[i, ] <- min(ref[[i]][[4]][!is.na(ref[[i]][[ 4]])])^refdCt[i, ]
  }
  
  # Sprawdzamy, czy mamy pliki dla wielu powtórzeñ - Uwaga na to!
  
  gene_remove <- c()
  remove_row <- c()
  
  for (i in ref_index){
    if (strsplit(inFile$name[i], '_')[[1]][2] == 'pow2'){
      for (j in ref_index){
        # Szukamy 2giego pliku bez pow2
        if (strsplit(inFile$name[i], '_')[[1]][1] == strsplit(inFile$name[j], '_')[[1]][1] & strsplit(inFile$name[j], '_')[[1]][2] != 'pow2'){
          # Dopasowanie indeksu genu z danych do indeksu wartoœci Q
          preQ[match(j, ref_index), ] <- (preQ[match(i, ref_index), ] + preQ[match(j, ref_index), ])/2
          gene_remove[k] <- i
          remove_row[k] <- match(i, ref_index)
          k <- k + 1
          
          
        }}}}
  remove_row <- remove_row[is.finite(remove_row)]
  if (length(remove_row) > 0){
    preQ <- preQ[-remove_row, ]
  }
  ref_index_old <- ref_index
  ref_index <- setdiff(ref_index, gene_remove)
  preQ <- matrix(preQ, nrow = length(ref_index), byrow = F)
  
  
  #Œrednia
  for (i in 1:dim(preQ)[2]){
    refQ[i] <- (prod(preQ[, i]))^(1/length(preQ[, i])) 
  }

  
### Przygotowanie danych
  # Obliczanie dCt
  Ct <- matrix( nrow = length(a), ncol = (length(a[[1]][[1]])-1))
  Q <- Ct
  Cal <- c()

  for (i in 1:length(a)){
    Ct[i, ] <- c(as.numeric(sub(",", ".", (a[[i]][[2]][-1]), fixed=TRUE)))
    Cal[i] <- max(a[[i]][[ 4]][!is.na(a[[i]][[4]])])
  }
  dCt <- Cal - Ct
  
  # Obliczamy Q, czyli potrzebne s¹ wydajnoœci
  for (i in 1:length(a)){
    Q[i, ] <- min(a[[i]][[4]][!is.na(a[[i]][[ 4]])])^dCt[i, ]
  }
  
  # Sprawdzamy, czy mamy pliki dla wielu powtórzeñ 

gene_index <- c(1:(length(a)+length(ref_index)))[-ref_index_old]
k <- 1
gene_remove <- c()
remove_row <- c()

for (i in gene_index){
  if (strsplit(inFile$name[i], '_')[[1]][2] == 'pow2'){
    for (j in gene_index){
      # Szukamy 2giego pliku bez pow2
      if (strsplit(inFile$name[i], '_')[[1]][1] == strsplit(inFile$name[j], '_')[[1]][1] & strsplit(inFile$name[j], '_')[[1]][2] != 'pow2'){
        # Dopasowanie indeksu genu z danych do indeksu wartoœci Q
        Q[match(j, gene_index), ] <- (Q[match(i, gene_index), ] + Q[match(j, gene_index), ])/2
        gene_remove[k] <- i
        remove_row[k] <- match(i, gene_index)
        k <- k + 1
        #Uaktualniæ gene index?
       
      }}}}
if (length(remove_row) > 0){
Q <- Q[-remove_row, ]
}
gene_index <- setdiff(gene_index, gene_remove)
Q <- matrix(Q, nrow = length(gene_index), byrow = F)
  
### Fold difference

Fd <- matrix( ncol = length(gene_index), nrow = (length(a[[1]][[1]])-1))

for (i in 1:length(gene_index)){
  Fd[, i] <- Q[i, ]/refQ
}
### Sprawdzamy, czy liczba próbek jest taka sama
if (dim(refQ)[2] != dim(Q)[2]){
  return(NULL)
}else{



    #Nazwy wierszy i kolumn
    
    Fd <- cbind(as.character(a[[1]][[1]][-1]), Fd)
    k <- 1
    columns <- c()
    
    for (i in gene_index){
     columns[k] <- strsplit(inFile$name[i], '_')[[1]][1]
    k <- k + 1
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