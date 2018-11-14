fluorescence1 <- function (inFile, threshold = 0.5){
  
  a <- list()
  for (i in 1:length(inFile[, 1])){

  data <- read.table(inFile$datapath[i], header = F, skip = 2, sep = "\t")
  data <- data[, -c(1:3)]
  data <- t(data[, 1:40])
  data <- cbind(1:40, data)
  data <- data.frame(data)

  # backsub - background substraction. More: https://www.rdocumentation.org/packages/qpcR/versions/1.3-1/topics/pcrbatch

  mode=modlist(data, cyc = 1, fluo = 2:40, model = l4)
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
  Ct = unlist(lapply(1:length(mode), function(x){
  unlist(predict(mode[[x]], newdata = data.frame(Fluo = threshold), which = "x"))}))
  
  # Adding data with Ct values to the list
  a[[i]] <- data.frame(cbind(data[, 1], Ct))
  }
  return(a)
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