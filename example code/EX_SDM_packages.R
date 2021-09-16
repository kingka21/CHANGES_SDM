#### SDMS ##### 

# RANDOM CODE AND PACKAGES # 

################################
#### GJAM CODE ####
#Gjam - Analyzes joint attribute data (e.g., species abundance) that are combinations of continuous and discrete data with Gibbs sampling. 
#I think for joint sp-dist modeling 
#Full model and computation details are described in Clark et al. (2018) <doi:10.1002/ecm.1241>.
install.packages("gjam")
library(gjam)

## discrete abundance with heterogeneous effort 
S   <- 5                             
n   <- 1000
eff <- list( columns = 1:S, values = round(runif(n,.5,5),1) )
f   <- gjamSimData(n, S, typeNames='DA', effort=eff)
ml  <- list(ng = 500, burnin = 50, typeNames = f$typeNames, effort = eff)
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
summary(out)

# repeat with ng = 2000, burnin = 500, then plot data:
pl  <- list(trueValues = f$trueValues)
gjamPlot(out, plotPars = pl)

########################################
#### code From class ####
### Run a full model 
sdmdata<-modeldata

#check for collinearity 
pairs(sdmdata[,2:5], cex=0.1, fig=TRUE)
ext <- extent(-90, -32, -33, 23)

#GAM model  (Generalized Additive Models ) 
install.packages('gam')
library(gam)
gam.model<-gam(abundance~forest_buf + wetland_buf + damdens_ws + agriculture_ws + area_ws + elev_ws, 
               data=MI_data_with_pred)
summary(gam.model)

### prediction 
full.pred<-predict(gam.model, all_vars) # dataset with all 51,000 lakes 

###reduced model 
red.gam.model<-gam(abundance~forest_buf + elev_ws, 
                   data=MI_data_with_pred)

summary(red.gam.model)

#reduced model prediction 

red.pred<-predict(red.gam.model, all_vars)

## map 
#get observations 
library(dismo)
bradypus <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
bradypus.obs <- read.csv(bradypus)
prj <- CRS("+proj=longlat +datum=NAD83")
brady.sp <- SpatialPoints(coords=bradypus.obs[,c(2,3)],
                          proj4string=prj)
# Reproject to matches 
crs(brady.sp)<-crs(predictors)

plot(full.pred, main= "Full model prediction")
points(brady.sp)
plot(red.pred, main= "Reduced model prediction")
points(brady.sp)

#### testing AUC  
samp <- sample(nrow(MI_data_with_pred), round(0.75 * nrow(MI_data_with_pred)))
traindata <- MI_data_with_pred[samp,]
traindata <- traindata[traindata[,1] == 1, 2:9]
testdata <- MI_data_with_pred[-samp,]

full <- evaluate(testdata[testdata==1,], testdata[testdata==0,], gam.model)
red <-evaluate(testdata[testdata==1,], testdata[testdata==0,], red.gam.model)

par(mfrow=c(1,2))
plot(full, 'ROC')
plot(red, "ROC")

#### GLMM ADMB package #### 
library(glmmADMB)
ln_effort=log(effort) 
nbinom=glmmadmb(catch~year+season+area+target+offset(ln_effort),family='nbinom',link='log')
summary(nbinom) 
ZI_nbinom=glmmadmb(catch~year+season+area+target+offset(ln_effort),family='nbinom',link='log',zeroInflation=TRUE)
summary(ZI_nbinom)
AICtab(nbinom,ZI_nbinom)