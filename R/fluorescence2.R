fluorescence2 <- function (inFile, inFile2, genes, threshold = 0.5){
  
  a <- list()
  plates <- data.frame(read.xlsx(inFile2$datapath, header = T, 3))
  sample <- data.frame(read.xlsx(inFile2$datapath, header = T, 4))
  sample[, 1] <- as.character(as.integer(sample[, 1]))
  
  for (i in 1:length(inFile[, 1])){

  data <- read.table(inFile$datapath[i], header = F, skip = 2, sep = "\t")
  names <- as.character(as.integer(data[, 1]))
  data <- data[, -c(1:3)]
  data <- t(data[, 1:40])
  data <- cbind(1:40, data)
  data <- data.frame(data)
  colnames(data) <- c('Cycle', names)
  
  # Changing well numbers to sample names

  sam_nam <- sapply(1:length(data[, 1]), function(j){ sample[data[, j], as.character(plates[genes[i], 2])]})
  #mode=modlist(data, cyc = 1, fluo = 2:40, model = l4)
  # model=pcrbatch(mode, plot=F) # Model is only here to check, if fitting failed
  # 
  # # Fitting model for not fitted stuff
  # not_fitted=grep("\\*",names(model[-1]),perl=T)
  # #not_fitted=c(1,2,3)
  # if(length(not_fitted)!=0){
  #   smoothPAR = list(span = 0.1)
  #   for(i in not_fitted){
  #     if(length(na.omit(data[, c(1,i+1)]))==0)
  #       mode[i]=modlist(na.omit(data[, c(1,i+1)]),1,2, model = l4,smooth = c("supsmu"), smoothPAR = smoothPAR)
  #   }
  #   model=pcrbatch(mode,plot=F)
  # }

  # Selecting fluorescence threshold to determine Ct value
  # Ct = unlist(lapply(1:length(mode), function(x){
  # unlist(predict(mode[[x]], newdata = data.frame(Fluo = threshold), which = "x"))}))
  # 
  # # Adding data with Ct values to the list
  # a[[i]] <- data.frame(cbind(sam_nam, Ct))
  }
  x=sapply(2:length(data[1, ]), function(j){ sample[colnames(data)[j] == sample[, 1], 1]})
  y=match(colnames(data), sample[, 1])
  #sample[y, 1]
  
  return(genes %in% plates[, 1])
}
# 
# #data
# setwd('C:/Users/jakpo_000/Desktop/Nauka/Inzynierka/Program/thesis/Dane')
# data=read.table('ATP6V1_jajniki_22_09_2010.txt', header = F, skip = 2, sep = "\t")
# #data=read.table('test.txt', header = T, sep = "\t")
# data <- data[, -c(1:3)]
# data <- t(data[, 1:40])
# data <- cbind(1:40, data)
# 
# # backsub - background substraction. More: https://www.rdocumentation.org/packages/qpcR/versions/1.3-1/topics/pcrbatch
# 
# mode=modlist(data, cyc = 1, fluo = 2:40, model = l4)
# model=pcrbatch(mode, plot=F) # Model is only here to check, if fitting failed
# 
# # Fitting model for not fitted stuff
# not_fitted=grep("\\*",names(model[-1]),perl=T)
# #not_fitted=c(1,2,3)
# if(length(not_fitted)!=0){
#   smoothPAR = list(span = 0.1)
#   for(i in not_fitted){
#     if(length(na.omit(data[, c(1,i+1)]))==0)
#       mode[i]=modlist(na.omit(data[, c(1,i+1)]),1,2, model = l4,smooth = c("supsmu"), smoothPAR = smoothPAR)
#   }
#   model=pcrbatch(mode,plot=F)
# }
# 
# ## Selecting fluorescence threshold to determine Ct value
# Ct = unlist(lapply(1:length(mode), function(x){
#   unlist(predict(mode[[x]], newdata = data.frame(Fluo = threshold), which = "x"))
# })) # works from fluo 0.29 for some