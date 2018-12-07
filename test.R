library(qpcR)
#data
setwd('C:/Users/jakpo_000/Desktop/Nauka/Inzynierka/Program/thesis/Dane/test/fluorescence')
data=read.table('ACTB_pl2_07_10_2014_clipped.txt', header = F, skip = 2, sep = "\t")
#data=read.table('test.txt', header = T, sep = "\t")
data <- data[, -c(1:3)]
data <- t(data[, 1:40])
data <- cbind(1:40, data)
data <- data.frame(data)
# backsub - background substraction. More: https://www.rdocumentation.org/packages/qpcR/versions/1.3-1/topics/pcrbatch

mode=modlist(data, cyc = 1, fluo = 2:40, model = l4, opt = F) 
model=pcrbatch(mode, plot=F) # Model is only here to check, if fitting failed

# Fitting model for not fitted stuff
# While loop here?
not_fitted=grep("\\*", names(model[-1]), perl=T)
#not_fitted=c(1,2,3)
if(length(not_fitted)!=0){
  smoothPAR = list(span = 0.1)
  for(i in not_fitted){
    if(length(na.omit(data[, c(1,i+1)]))==0)
    mode[i]=modlist(na.omit(data[, c(1,i+1)]),1,2, model = l4,smooth = c("supsmu"), smoothPAR = smoothPAR)
  }
  model=pcrbatch(mode, plot=F)
}

## Selecting fluorescence threshold to determine Ct value
Ct = unlist(lapply(1:length(mode), function(x){
  unlist(predict(mode[[x]], newdata = data.frame(Fluo = 0.4), which = "x"))
  })) # works from fluo 0.29 for some

# #equation formula
# form = mode[[1]]$MODEL$expr
# # inverse function - PERFECT to calculate Ct!
# inv = mode[[1]]$MODEL$inv
# parm <- unlist(mode[[1]]$model[4:7])
# Ct=inv(0.5, parm)

# Sent code
mode=modlist(im[1:plateau.end,],1,2:ncol(im),modtype,backsub =1:20,opt=FALSE)
model=pcrbatch(mode,plot=T)
not_fitted=grep("\\*",names(model[-1]),perl=T)
if(length(not_fitted)!=0){
  smoothPAR = list(span = 0.1)
  for(i in not_fitted){
    if(length(na.omit(im[1:plateau.end,c(1,i+1)]))==0)
      im[1:plateau.end,c(i+1)]=0
    mode[i]=modlist(na.omit(im[1:plateau.end,c(1,i+1)]),1,2, model = l4,backsub =backsub,opt=opt,smooth = c("supsmu"),
                    smoothPAR = smoothPAR)
  }
  model=pcrbatch(mode,plot=F)
}
