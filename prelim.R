library("R.matlab")
library("glmnet")

# Convert matlab matrices to R data.frame
# Create a file with converted data so that we don't need to repeat this process
if (!file.exists("experiments.Rdata")){
  dr = "BiomechParam/"
  file.names <- dir(dr, pattern ="mat")
  
  ## BIG LIST
  biglist = list()
  for(i in 1:length(file.names))
  {
    print(paste("Processing",file.names[i]))
    ss = strsplit(file.names[i],split = ".",fixed=TRUE)[[1]][1]
    ss = strsplit(ss,split = "_",fixed=TRUE)[[1]][2]
    data = readMat(paste0(dr,"BiomechParam_R105.mat"))
    
    nobs = dim(data$BiomechParamMat)[3]
    for (j in 1:nobs)
    {
      newrow.lst = data$BiomechParamMat[,,j]
      newrow.lst$Subject = ss
      for (nm in names(newrow.lst)){
        biglist[[nm]] = rbind(biglist[[nm]], newrow.lst[[nm]])
      }
    }
  }
  experiments = biglist
  
  save(experiments, file="experiments.Rdata")
}

# Load the data
load("experiments.Rdata")
data = experiments

# remove KAM from the data 
KAM = data$KAM
data$KAM = NULL

# create a data.frame from the list of features
df = data.frame(data)

# split to training and testing sets
subjects = unique(df$Subject)
n = length(subjects)
test.subj = subjects[sample(n)[1:floor(n*0.10)]]
train.mask = !(df$Subject %in% test.subj)

# make factors from strings
df$Subject = factor(df$Subject)
df$TrialName = factor(df$TrialName)

# build baseline models
model = lm(KAM[train.mask,] ~ ., data = df[train.mask,-1])
preds = predict(model, newdata = df[!train.mask,-1])

# check how they perform
mean((preds - KAM[!train.mask,])**2) / mean((KAM[!train.mask,])**2)

# realize that this is too good.