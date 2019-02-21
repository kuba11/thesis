library(qpcR)
#data
setwd('C:/Users/jakpo_000/Desktop/Nauka/Inzynierka/Program/thesis/Dane/test/fluorescence')
data=read.table('CCL2_pl1_07_10_2014_clipped.txt', header = F, skip = 2, sep = "\t")
#data=read.table('test.txt', header = T, sep = "\t")
data <- data[, -c(1:3)]
data <- t(data[, 1:40])
data <- cbind(1:40, data)
data <- data.frame(data)
# backsub - background substraction. More: https://www.rdocumentation.org/packages/qpcR/versions/1.3-1/topics/pcrbatch

mode=modlist(data, cyc = 1, fluo = 2:69, model = l4, opt = F) 
model=pcrbatch(mode, plot=F) # Model is only here to check, if fitting failed

# Fitting model for not fitted stuff
star12 <- grep("\\*", names(model[-1]), perl=T) #fitting failed + no sigmoidal structure
star2 <- grep("\\*\\*",names(model[-1]), perl=T) #no sigmoidal structure
#not_fitted=c(1,2,3)
if(length(star2)!=0){
  
  for(i in not_fitted){
    if(length(na.omit(data[, c(1,i+1)]))==0)
    mode[i]
    }
  model=pcrbatch(mode, plot=F)
}

# ## Selecting fluorescence threshold to determine Ct value
# Ct = unlist(lapply(which(lapply(mode, function(y) length(y))==15), function(x){
#   unlist(predict(mode[[x]], newdata = data.frame(Fluo = 0.5), which = "x"))
#   })) # works from fluo 0.29 for some
# 

# Ct calculation. Needs to be in a loop to see if the modelling succeeded
Ct <- c()
for(x in 1:length(mode)){
  if(length(mode[[x]]) == 15) #if = 15 then ok, if 4 then modeling failed
  Ct <- c(Ct, predict(mode[[x]], newdata = data.frame(Fluo = 0.5), which = "x"))
  else
    Ct <- c(Ct, NA)
  }

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







for(x in 1:28){#length(mode)){
z=predict(mode[[28]], newdata = data.frame(Fluo = 0.5), which = "x")
}
x=which(lapply(mode, function(y) length(y))==15)

