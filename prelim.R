library("R.matlab")
library("glmnet")

# Convert matlab matrices to R data.frame
# Create a file with converted data so that we don't need to repeat this process
if (!file.exists("experiments.Rdata")){
  dr = "../BiomechParam/"
  file.names <- dir(dr, pattern ="mat")
  
  ## BIG LIST
  biglist = list()
  for(i in 1:length(file.names))
  {
    print(paste("Processing",file.names[i]))
    ss = strsplit(file.names[i],split = ".",fixed=TRUE)[[1]][1]
    ss = strsplit(ss,split = "_",fixed=TRUE)[[1]][2]
    data = readMat(paste0(dr,"BiomechParam_",ss,".mat"))
    
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

res = c()
res

dim(biglist[[1]])

# Load the data
load("experiments.Rdata")
data = experiments

# remove KAM from the data 
KAM = data$KAM
data$KAM = NULL
data$Fz = NULL
data$COPx = NULL
data$COPy = NULL
data$Alignment = NULL
data$KneeAlignment = NULL

nc = 1
data$TibiaAngle = prcomp(data$TibiaAngle)$x[,1:nc]
data$PelvicAxialRot = prcomp(data$PelvicAxialRot)$x[,1:nc]
data$PelvicList = prcomp(data$PelvicList)$x[,1:nc]
data$KneeFlexion = prcomp(data$KneeFlexion)$x[,1:nc]
data$VarusThrust = prcomp(data$VarusThrust)$x[,1:nc]
data$TrunkSway = prcomp(data$TrunkSway)$x[,1:nc]

# create a data.frame from the list of features
df = data.frame(data)
df$target = df$StepMaxRedux
df$StepMaxRedux = NULL
df$StepP1Redux = NULL
df$StepP2Redux = NULL
df$StepwiseBenefitP1 = NULL
df$StepwiseBenefitP2 = NULL
df$BinwiseBenefitP1 = NULL
df$BinwiseBenefitP2 = NULL
df$BinwiseP1Redux = NULL
df$BinwiseP2Redux = NULL
df$FPABin = NULL
df$TrialName = NULL

# split to training and testing sets
subjects = unique(df$Subject)
n = length(subjects)
test.subj = subjects[sample(n)[1]]
train.mask = !(df$Subject %in% test.subj)

# make factors from strings
df$Subject = factor(df$Subject)
#df$TrialName = factor(df$TrialName)
df$Subject = NULL

# Build models
model.lm = lm(target ~ ., data = df[train.mask,])
preds = predict(model.lm, newdata = df[!train.mask,])

library(ggplot2)
paper.theme = theme_set(theme_grey(base_size = 22)) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1))

FPA = c(df$FPA[!train.mask],df$FPA[!train.mask])
reduction = c(df$target[!train.mask], preds)
group = c(rep("true", sum(!train.mask)),rep("predicted", sum(!train.mask)))
df.plot = data.frame(FPA=FPA, reduction=reduction, group=group)
ggplot(df.plot, aes(x=FPA, y=reduction, group=group, color=group)) + paper.theme +
  geom_point()

### OLD STUFF
# build baseline models
model = lm(KAM[train.mask,] ~ ., data = df[train.mask,c(2:ncol(df))])
preds = predict(model, newdata = df[!train.mask,c(2:ncol(df))])

# check how they perform
var.exp = 1 - colMeans((preds - KAM[!train.mask,])**2) / colMeans((KAM[!train.mask,])**2)
var.exp[var.exp<0] = 0
res = rbind(res, var.exp)
plot(sqrt(var.exp))

# matplot(t(preds),t='l')
# matplot(t(KAM[!train.mask,]),t='l')

## Some plots
mm = colMeans(KAM[!train.mask,])
#mm[] = 0

sd.true = sqrt(diag(cov(KAM[!train.mask,])))
sd.pred = sqrt(diag(cov(preds - KAM[!train.mask,])))

plot(sd.true,ylim=c(0,max(sd.true)))
lines(sd.pred)




# matplot(t(preds - KAM[!train.mask,]) + mm,t='l',add=TRUE)
# lines(mm + sd.true, lwd=4, t='l', ylim=c(-max(sd.true) + min(mm),max(sd.true) + max(mm)))
# lines(mm - sd.true, lwd=4, t="l")
# 
# 
# lines(mm + sd.pred, lwd=4, t="l", col=1)
# lines(mm - sd.pred, lwd=4, t="l", col=1)
# 
# matplot(t(KAM) - mm,t='l')
#plot(sd.pred**2 / sd.true**2)
# realize that this is too good.

plot(colMeans(res))