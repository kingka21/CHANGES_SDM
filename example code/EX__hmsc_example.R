#HMSC: Hierarchical Modelling of Species Communities (HMSC) ####
#is a model-based approach for analyzing community ecological data. 
#This package implements it in the Bayesian framework with Gibbs Markov chain Monte Carlo (MCMC) sampling (Tikhonov et al. (2020) <doi:10.1111/2041-210X.13345>).
#can be used for joint sp dist modeling but also single species 
#https://cran.r-project.org/web/packages/Hmsc/vignettes/vignette_1_univariate.pdf
library(Hmsc)

##### Introduction to linear models #### 
# example of HMSC to fit a univariate linear model.
#we also apply the basic lm function to the same data and compare the results.
set.seed(1)
n = 50 #number of data points 
x = rnorm(n) # 50 parameter estimates from a continuous covariate  
alpha = 0 # intercept 
beta = 1 #slope
sigma = 1 #sd
L = alpha + beta*x #linear predictor 
y = L + rnorm(n, sd = sigma) 
plot(x, y, las=1)

#using the lm function - constructs and fits the model at the same time 
df = data.frame(x,y)
m.lm = lm(y ~ x, data=df) 
summary(m.lm)

#check model assumptions of residuals being normally distributed and homoscedastic
nres.lm = rstandard(m.lm)  #residuals 
preds.lm = fitted.values(m.lm) 
par(mfrow=c(1,2))
hist(nres.lm, las = 1) 
plot(preds.lm,nres.lm, las = 1) 
abline(a=0,b=0)

#same thing as lm  but using the hmsc 
Y = as.matrix(y)
XData = data.frame(x = x)
m = Hmsc(Y = Y, XData = XData, XFormula = ~x, distr = "normal") # this constructs the model object # normal dist for continuous data

#To fit the HMSC model with Bayesian inference, we use the sampleMcmc function.
nChains = 2  #how many chains to sample (nChains)
thin=5 #how many samples to keep 
samples = 1000 #how many samples to obtain per chain (samples)
transient=500*thin #how long transient (also called burn-in) 
verbose = 500*thin #how frequently we wish to see the progress of the MCMC sampling (verbose). 
#fit the model 
m1 = sampleMcmc(m, thin = thin, samples = samples, transient = transient, nChains = nChains, verbose = verbose)

#To look at the parameter estimates, we extract the posterior distribution
mpost = convertToCodaObject(m1) 
summary(mpost$Beta)

#investigate MCMC convergence 
plot(mpost$Beta) #trace plot and density plot - two chains(red and black) look identical, the chains mix very well, i.e. they go fast up and down without any apparent autocorrelation, they seem to have reached a stationary distribution  
effectiveSize(mpost$Beta) #did we set enough samples? effective sample sizes are very close to the theoretical value of the actual number of samples, which is 2000 (1000 per chain)
gelman.diag(mpost$Beta,multivariate=FALSE)$psrf #scale reduction factors are very close to one, which indicates that two chains gave consistent results
par(mfrow=c(1,2)) #plot if you have a lot of parameters 
hist(effectiveSize(mpost$Beta), main="ess(beta)") 
hist(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")

# get RMSE and R2 by looking at the predicted dist
preds = computePredictedValues(m1) #posterior predictive distribution 
evaluateModelFit(hM=m1, predY=preds)

#diagnostic plots 
preds.mean = apply(preds, FUN=mean, MARGIN=1) #get the posterior means 
nres = scale(y-preds.mean) #standardize the residuals
par(mfrow=c(1,2))
hist(nres, las = 1)
plot(preds.mean,nres, las = 1) 
abline(a=0,b=0)


#### Presence-abs data probit model #### 
#probit is similar to logistic regression - uses a probit instead of logit link function
y = 1*(L + rnorm(n, sd = 1) > 0) #sets to 0s and 1s
plot(x,y, las = 1)

Y=as.matrix(y)
m = Hmsc(Y=Y, XData=XData, XFormula=~x, distr="probit")
verbose = 0
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)
mpost = convertToCodaObject(m)
summary(mpost$Beta)
effectiveSize(mpost$Beta)
gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
preds = computePredictedValues(m) 
evaluateModelFit(hM=m, predY=preds) #this gives you RMSE, AUC, and Tjur R2

#### Poisson model for count data #### 
#a non-negative integer 0,1,2,3,..., representing the count of individuals. 
y = rpois(n, lambda = exp(L)) 
plot(x,y,las=1)
Y=as.matrix(y)
m = Hmsc(Y=Y, XData=XData, XFormula=~x, distr="poisson")
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)
mpost = convertToCodaObject(m) 
summary(mpost$Beta)
effectiveSize(mpost$Beta)
gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
preds = computePredictedValues(m, expected = FALSE) #expected=FALSE, the predictions are not expected values (e.g. on average we expect to see 2.3 individuals), but a posterior predictive distribution of data, i.e., actual counts 0,1,2,3,that include the Poisson sampling error. We have selected expected = FALSE so that presence-absences can be inferred from the predictions
evaluateModelFit(hM=m, predY=preds) #get RMSE, SR2 (pseudo-R2), O.AUC, O.TjurR2, O.RMSE (O stand for how well sp occurrences are predicted), C.SR2, C.RMSE (conditional on presence, how well the abundance is predicted)

#### lognormal Poisson model for count data #### 
#To allow for more variation, fit a lognormal Poisson model (similar to neg binomial). zero inflated not an option 
y = rpois(n, lambda = exp(L+rnorm(n, sd=2)))
plot(x,y, las = 1)
Y=as.matrix(y)
m = Hmsc(Y=Y, XData=XData, XFormula=~x, distr="lognormal poisson")
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)
mpost = convertToCodaObject(m) 
summary(mpost$Beta)
effectiveSize(mpost$Beta)
gelman.diag(mpost$Beta, multivariate=FALSE)$psrf
preds = computePredictedValues(m, expected = FALSE) 
evaluateModelFit(hM=m, predY=preds)

#### hierarchical structure ####
# note this is how you would do it with max likelihood aproach: lmer(y~x+(1|plot.id),data=da)
n = 100
x = rnorm(n)
alpha = 0
beta = 1
sigma = 1
L = alpha + beta*x
np = 10 ### this is the grouping variable 
sigma.plot = 1
sample.id = 1:n
plot.id = sample(1:np, n, replace = TRUE) 
ap = rnorm(np, sd = sigma.plot)
a = ap[plot.id]
y = L + a + rnorm(n, sd = sigma)
plot.id = as.factor(plot.id)
plot(x,y,col = plot.id, las = 1)


XData = data.frame(x = x) 
Y = as.matrix(y)
studyDesign = data.frame(sample = as.factor(sample.id), plot = as.factor(plot.id)) #individual samples and their grouping factor
rL = HmscRandomLevel(units = studyDesign$plot) #default options for latent variables. this function is to reduce dimensionality of residual association. latent features should not exceed the number of species  
m = Hmsc(Y=Y, XData=XData, XFormula= ~x,
         studyDesign=studyDesign, ranLevels=list("plot"=rL)) #list of study design that will be random
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)
preds=computePredictedValues(m) 
MF=evaluateModelFit(hM=m, predY=preds)
MF #this R2 is explanatory and is missleading because of overfitting. better to evaluate with cross-validation predictive R2

#cross-validation 
#2-fold cross-validation - randomly sample into 2 groups, model fitting and predictions are made separately for each fold 
partition=createPartition(m, nfolds=2, column="sample") #this partitions based on the samples 
preds=computePredictedValues(m, partition = partition)
MF=evaluateModelFit(hM=m, predY=preds)
MF #predictive R2 is lower than the explanatory R2 above 

#### spatial structure: including coordinates instead of grouping by plot #### 
sigma.spatial = 2 
alpha.spatial = 0.5 
sample.id = rep(NA,n) 
for (i in 1:n){ 
  sample.id[i] = paste0("location_",as.character(i))
}
sample.id = as.factor(sample.id) 
xycoords = matrix(runif(2*n), ncol=2) 
rownames(xycoords) = sample.id 
colnames(xycoords) = c("x-coordinate","y-coordinate") 
a = MASS::mvrnorm(mu=rep(0,n), Sigma = sigma.spatial^2*exp(-as.matrix(dist(xycoords))/alpha.spatial))
y = L + a + rnorm(n, sd = sigma) 
Y=as.matrix(y) 
colfunc = colorRampPalette(c("cyan", "red")) 
ncols = 100 
cols = colfunc(100)
par(mfrow=c(1,2)) 
#example ploting the predictor (x) and the response (y) to look for spatial pattern
for (i in 1:2){ 
  if (i==1) value = x 
  if (i==2) value = y 
  value = value-min(value) 
  value = 1+(ncols-1)*value/max(value) 
  plot(xycoords[,1],xycoords[,2],col=cols[value],pch=16,main=c("x","y")[i], asp=1)
}

#fit model 
##As the model includes a spatially structured random effect, its predictive power is based on both the fixed and the random effects. 
#because of random effect, the model can utilize observed data from nearby sampling units included in model fitting when predicting the response for a focal sampling unit that is not included in model fitting
studyDesign = data.frame(sample = sample.id) 
rL = HmscRandomLevel(sData = xycoords) #use the sData argument to construct the random effect for spatially explicit
m = Hmsc(Y=Y, XData=XData, XFormula=~x, 
         studyDesign=studyDesign, ranLevels=list("sample"=rL))
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient, 
               nChains = nChains, verbose = verbose)
preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds) 
MF$R2
##cross-val
partition = createPartition(m, nfolds = 2, column = "sample") 
preds = computePredictedValues(m, partition=partition)
MF = evaluateModelFit(hM=m, predY=preds) 
MF$R2

#look at estimates of Beta parameters 
postBeta = getPostEstimate(m, parName = "Beta")
plotBeta(m, post = postBeta, param = "Support", supportLevel = 0.95, colors = colorRampPalette(c("darkgreen", "white", "purple"))) #95% credible interval

