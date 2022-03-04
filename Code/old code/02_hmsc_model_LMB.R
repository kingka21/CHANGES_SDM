#HMSC: Hierarchical Modelling of Species Communities (HMSC) ####
#is a model-based approach for analyzing community ecological data. 
#This package implements it in the Bayesian framework with Gibbs Markov chain Monte Carlo (MCMC) sampling (Tikhonov et al. (2020) <doi:10.1111/2041-210X.13345>).
#can be used for joint sp dist modeling but also single species 
#https://cran.r-project.org/web/packages/Hmsc/vignettes/vignette_1_univariate.pdf
library(Hmsc)
library(MCMCvis)
library(dotwhisker)
#note that I cant do zero-inflated
#note that as default, scaling is applied for X, but not for Y matrix

#### standardize and transform predictor values ####
lmb_count_dat<-lmb_count_dat %>%
  ungroup()

##subset into four different gear response variables and remove 0s, select subset of metrics not correlated for analysis and remove NAs #
dat_shock<-filter(lmb_count_dat, fish_count_SHOCK >0) %>%
  na.omit()
dat_fyke<-filter(lmb_count_dat, fish_count_FT_NET >0) %>%
  na.omit() 
dat_gill<-filter(lmb_count_dat, fish_count_GILL >0) %>% 
  na.omit() 
dat_seine<-filter(lmb_count_dat, fish_count_SEINE >0) %>% 
  na.omit() 

#your model works on the representation given by its input predictors, differently scaled test/prediction data would potentially lead to over/under-exaggeration of a feature.
#standardize and transform 
#shock 
dat_shock$z_order<-as.numeric(scale(log(dat_shock$lake_order+0.001))) # min is 0 because these are the isolated
dat_shock$z_lake_area<-as.numeric(scale(log(dat_shock$lake_area_m2))) #good 
dat_shock$z_perim_km<-as.numeric(scale(log(dat_shock$lake_perim_km))) 
dat_shock$z_ws_area_km2<-as.numeric(scale(log(dat_shock$ws_area_km2)))  
dat_shock$z_ws_urban<-as.numeric(scale(asin(sqrt(dat_shock$ws_urban_prop)))) #good
dat_shock$z_ws_agriculture<-as.numeric(scale(asin(sqrt(dat_shock$ws_agriculture_prop)))) #ok
dat_shock$z_ws_forest<-as.numeric(scale(asin(sqrt(dat_shock$ws_forest_prop)))) #good
dat_shock$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat_shock$ws_wetland_prop)))) #good
dat_shock$z_ws_shrub<-as.numeric(scale(asin(sqrt(dat_shock$ws_shrub_prop)))) #ok
dat_shock$z_ws_slope<-as.numeric(scale(log(dat_shock$ws_slope_deg + 0.001))) #has 0s so added 0.001 
dat_shock$z_ws_elevation<-as.numeric(scale(log(dat_shock$ws_mean_elevation_m))) #good 
dat_shock$z_dd_mean<-as.numeric(scale(log(dat_shock$dd_mean))) #good
dat_shock$z_houses<-as.numeric(scale(log(dat_shock$houses_km+ 0.001))) #ok 
dat_shock$z_secchi<-as.numeric(scale(log(dat_shock$secchi_m)))
dat_shock$z_depth<-as.numeric(scale(log(dat_shock$maxdepth_m)))
dat_shock$z_geom<-as.numeric(scale(log(dat_shock$geom_ratio)))

#fyke
dat_fyke$z_order<-as.numeric(scale(log(dat_fyke$lake_order+0.001))) # min is 0 because these are the isolated
dat_fyke$z_lake_area<-as.numeric(scale(log(dat_fyke$lake_area_m2))) #good 
dat_fyke$z_perim_km<-as.numeric(scale(log(dat_fyke$lake_perim_km))) 
dat_fyke$z_ws_area_km2<-as.numeric(scale(log(dat_fyke$ws_area_km2)))  
dat_fyke$z_ws_urban<-as.numeric(scale(asin(sqrt(dat_fyke$ws_urban_prop)))) #good
dat_fyke$z_ws_agriculture<-as.numeric(scale(asin(sqrt(dat_fyke$ws_agriculture_prop)))) #ok
dat_fyke$z_ws_forest<-as.numeric(scale(asin(sqrt(dat_fyke$ws_forest_prop)))) #good
dat_fyke$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat_fyke$ws_wetland_prop)))) #good
dat_fyke$z_ws_shrub<-as.numeric(scale(asin(sqrt(dat_fyke$ws_shrub_prop)))) #ok
dat_fyke$z_ws_slope<-as.numeric(scale(log(dat_fyke$ws_slope_deg + 0.001))) #has 0s so added 0.001 
dat_fyke$z_ws_elevation<-as.numeric(scale(log(dat_fyke$ws_mean_elevation_m))) #good 
dat_fyke$z_dd_mean<-as.numeric(scale(log(dat_fyke$dd_mean))) #good
dat_fyke$z_houses<-as.numeric(scale(log(dat_fyke$houses_km+ 0.001))) #ok 
dat_fyke$z_secchi<-as.numeric(scale(log(dat_fyke$secchi_m)))
dat_fyke$z_depth<-as.numeric(scale(log(dat_fyke$maxdepth_m)))
dat_fyke$z_geom<-as.numeric(scale(log(dat_fyke$geom_ratio)))

#gill 
dat_gill$z_order<-as.numeric(scale(log(dat_gill$lake_order+0.001))) # min is 0 because these are the isolated
dat_gill$z_lake_area<-as.numeric(scale(log(dat_gill$lake_area_m2))) #good 
dat_gill$z_perim_km<-as.numeric(scale(log(dat_gill$lake_perim_km))) 
dat_gill$z_ws_area_km2<-as.numeric(scale(log(dat_gill$ws_area_km2)))  
dat_gill$z_ws_urban<-as.numeric(scale(asin(sqrt(dat_gill$ws_urban_prop)))) #good
dat_gill$z_ws_agriculture<-as.numeric(scale(asin(sqrt(dat_gill$ws_agriculture_prop)))) #ok
dat_gill$z_ws_forest<-as.numeric(scale(asin(sqrt(dat_gill$ws_forest_prop)))) #good
dat_gill$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat_gill$ws_wetland_prop)))) #good
dat_gill$z_ws_shrub<-as.numeric(scale(asin(sqrt(dat_gill$ws_shrub_prop)))) #ok
dat_gill$z_ws_slope<-as.numeric(scale(log(dat_gill$ws_slope_deg + 0.001))) #has 0s so added 0.001 
dat_gill$z_ws_elevation<-as.numeric(scale(log(dat_gill$ws_mean_elevation_m))) #good 
dat_gill$z_dd_mean<-as.numeric(scale(log(dat_gill$dd_mean))) #good
dat_gill$z_houses<-as.numeric(scale(log(dat_gill$houses_km+ 0.001))) #ok 
dat_gill$z_secchi<-as.numeric(scale(log(dat_gill$secchi_m)))
dat_gill$z_depth<-as.numeric(scale(log(dat_gill$maxdepth_m)))
dat_gill$z_geom<-as.numeric(scale(log(dat_gill$geom_ratio)))

#seine 
dat_seine$z_order<-as.numeric(scale(log(dat_seine$lake_order+0.001))) # min is 0 because these are the isolated
dat_seine$z_lake_area<-as.numeric(scale(log(dat_seine$lake_area_m2))) #good 
dat_seine$z_perim_km<-as.numeric(scale(log(dat_seine$lake_perim_km))) 
dat_seine$z_ws_area_km2<-as.numeric(scale(log(dat_seine$ws_area_km2)))  
dat_seine$z_ws_urban<-as.numeric(scale(asin(sqrt(dat_seine$ws_urban_prop)))) #good
dat_seine$z_ws_agriculture<-as.numeric(scale(asin(sqrt(dat_seine$ws_agriculture_prop)))) #ok
dat_seine$z_ws_forest<-as.numeric(scale(asin(sqrt(dat_seine$ws_forest_prop)))) #good
dat_seine$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat_seine$ws_wetland_prop)))) #good
dat_seine$z_ws_shrub<-as.numeric(scale(asin(sqrt(dat_seine$ws_shrub_prop)))) #ok
dat_seine$z_ws_slope<-as.numeric(scale(log(dat_seine$ws_slope_deg + 0.001))) #has 0s so added 0.001 
dat_seine$z_ws_elevation<-as.numeric(scale(log(dat_seine$ws_mean_elevation_m))) #good 
dat_seine$z_dd_mean<-as.numeric(scale(log(dat_seine$dd_mean))) #good
dat_seine$z_houses<-as.numeric(scale(log(dat_seine$houses_km+ 0.001))) #ok 
dat_seine$z_secchi<-as.numeric(scale(log(dat_seine$secchi_m)))
dat_seine$z_depth<-as.numeric(scale(log(dat_seine$maxdepth_m)))
dat_seine$z_geom<-as.numeric(scale(log(dat_seine$geom_ratio)))


##### linear models #### 
# example of HMSC to fit a linear model.
#using the lm function - constructs and fits the model at the same time 
cpue_shock<-filter(lmb_count_dat, SHOCK >0) %>% 
  dplyr::select(new_key, FMU_Code, SHOCK, z_lake_area, z_perim_km, z_order, z_ws_area_km2, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation, z_dd_mean, z_houses, z_secchi) %>%
  na.omit()
m.lm_shock = lm(log(SHOCK) ~ z_lake_area+ z_perim_km+ z_order+ z_ws_area_km2+ z_ws_urban+ z_ws_forest+ z_ws_agriculture+ z_ws_shrub+ z_ws_wetland+ z_ws_slope+ z_ws_elevation+ z_dd_mean+ z_houses+ z_secchi,
              data=dat_shock) 
summary(m.lm_shock)
#for reduced model chose sig lake area, elevation, dd, secchi, ws area 
m.red_shock<-lm(log(SHOCK) ~ z_lake_area+ z_ws_area_km2+ z_ws_elevation+ z_dd_mean+ z_secchi,
                data=cpue_shock) 
summary(m.red_shock)

#check model assumptions of residuals being normally distributed and homoscedastic
nres.lm = rstandard(m.lm_shock)  #residuals 
preds.lm = fitted.values(m.lm_shock) 
par(mfrow=c(1,2))
hist(nres.lm, las = 1) 
plot(preds.lm,nres.lm, las = 1) 
abline(a=0,b=0)

#To fit the HMSC model with Bayesian inference, we use the sampleMcmc function.
nChains = 3  #how many chains to sample (nChains)
thin=1 #how many samples to keep 
samples = 10000 #how many samples to obtain per chain (samples)
transient=5000 #how long transient (also called burn-in) 
verbose = 5000 #how frequently we wish to see the progress of the MCMC sampling (verbose). 

#same thing as lm  but using the hmsc 
#shock model 
Y = as.matrix(log(dat_shock$SHOCK))
XData = dplyr::select(dat_shock, z_lake_area , z_geom, z_order, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation, z_dd_mean, z_houses, z_secchi)
#XData = data.frame(cpue_shock$z_lake_area, cpue_shock$z_perim_km, cpue_shock$z_order, cpue_shock$z_ws_area_km2, cpue_shock$z_ws_urban, cpue_shock$z_ws_forest, cpue_shock$z_ws_agriculture, cpue_shock$z_ws_shrub, cpue_shock$z_ws_wetland, cpue_shock$z_ws_slope, cpue_shock$z_ws_elevation, cpue_shock$z_dd_mean, cpue_shock$z_houses, cpue_shock$z_secchi)
m_shock = Hmsc(Y = Y, XData = XData, XFormula = ~z_lake_area + z_geom+ z_order+ z_ws_urban+ z_ws_forest+ z_ws_agriculture+ z_ws_shrub+ z_ws_wetland+ z_ws_slope+ z_ws_elevation+ z_dd_mean+ z_houses+ z_secchi, 
         XScale = FALSE, distr = "normal") # this constructs the model object # normal dist for continuous data #alrady scaled my parameters 
#fit the model
m_shock_l = sampleMcmc(m_shock, thin = thin, samples = samples, transient = transient, nChains = nChains, verbose = verbose)

#fyke model 
Y = as.matrix(log(dat_fyke$FT_NET))
XData = dplyr::select(dat_fyke, z_lake_area , z_geom, z_order, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation, z_dd_mean, z_houses, z_secchi)
#XData = data.frame(cpue_shock$z_lake_area, cpue_shock$z_perim_km, cpue_shock$z_order, cpue_shock$z_ws_area_km2, cpue_shock$z_ws_urban, cpue_shock$z_ws_forest, cpue_shock$z_ws_agriculture, cpue_shock$z_ws_shrub, cpue_shock$z_ws_wetland, cpue_shock$z_ws_slope, cpue_shock$z_ws_elevation, cpue_shock$z_dd_mean, cpue_shock$z_houses, cpue_shock$z_secchi)
m_fyke = Hmsc(Y = Y, XData = XData, XFormula = ~z_lake_area + z_geom+ z_order+ z_ws_urban+ z_ws_forest+ z_ws_agriculture+ z_ws_shrub+ z_ws_wetland+ z_ws_slope+ z_ws_elevation+ z_dd_mean+ z_houses+ z_secchi, 
               XScale = FALSE, distr = "normal") # this constructs the model object # normal dist for continuous data #alrady scaled my parameters 
#fit the model
m_fyke_l = sampleMcmc(m_fyke, thin = thin, samples = samples, transient = transient, nChains = nChains, verbose = verbose)

#gill model 
Y = as.matrix(log(dat_gill$GILL))
XData = dplyr::select(dat_gill, z_lake_area , z_geom, z_order, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation, z_dd_mean, z_houses, z_secchi)
#XData = data.frame(cpue_shock$z_lake_area, cpue_shock$z_perim_km, cpue_shock$z_order, cpue_shock$z_ws_area_km2, cpue_shock$z_ws_urban, cpue_shock$z_ws_forest, cpue_shock$z_ws_agriculture, cpue_shock$z_ws_shrub, cpue_shock$z_ws_wetland, cpue_shock$z_ws_slope, cpue_shock$z_ws_elevation, cpue_shock$z_dd_mean, cpue_shock$z_houses, cpue_shock$z_secchi)
m_gill = Hmsc(Y = Y, XData = XData, XFormula = ~z_lake_area + z_geom+ z_order+ z_ws_urban+ z_ws_forest+ z_ws_agriculture+ z_ws_shrub+ z_ws_wetland+ z_ws_slope+ z_ws_elevation+ z_dd_mean+ z_houses+ z_secchi, 
              XScale = FALSE, distr = "normal") # this constructs the model object # normal dist for continuous data #alrady scaled my parameters 
#fit the model
m_gill_l = sampleMcmc(m_gill, thin = thin, samples = samples, transient = transient, nChains = nChains, verbose = verbose)

#seine model 
Y = as.matrix(log(dat_seine$SEINE))
XData = dplyr::select(dat_seine, z_lake_area , z_geom, z_order, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation, z_dd_mean, z_houses, z_secchi)
#XData = data.frame(cpue_shock$z_lake_area, cpue_shock$z_perim_km, cpue_shock$z_order, cpue_shock$z_ws_area_km2, cpue_shock$z_ws_urban, cpue_shock$z_ws_forest, cpue_shock$z_ws_agriculture, cpue_shock$z_ws_shrub, cpue_shock$z_ws_wetland, cpue_shock$z_ws_slope, cpue_shock$z_ws_elevation, cpue_shock$z_dd_mean, cpue_shock$z_houses, cpue_shock$z_secchi)
m_seine = Hmsc(Y = Y, XData = XData, XFormula = ~z_lake_area + z_geom+ z_order+ z_ws_urban+ z_ws_forest+ z_ws_agriculture+ z_ws_shrub+ z_ws_wetland+ z_ws_slope+ z_ws_elevation+ z_dd_mean+ z_houses+ z_secchi, 
              XScale = FALSE, distr = "normal") # this constructs the model object # normal dist for continuous data #alrady scaled my parameters 
#fit the model
m_seine_l = sampleMcmc(m_seine, thin = thin, samples = samples, transient = transient, nChains = nChains, verbose = verbose)

#### investigating output ####
# extract the posterior distribution
shock_post = convertToCodaObject(m_shock_l) 
fyke_post = convertToCodaObject(m_fyke_l) 
gill_post = convertToCodaObject(m_gill_l) 
seine_post = convertToCodaObject(m_seine_l) 

#pull out means and 90% credible intervals 
shock_betas<-MCMCsummary(shock_post$Beta, 
            Rhat = TRUE, 
            n.eff = TRUE, 
            probs = c(0.05, 0.5, 0.95), 
            round = 2)
shock_betas= shock_betas[-1,]

fyke_betas<-MCMCsummary(fyke_post$Beta, 
                         Rhat = TRUE, 
                         n.eff = TRUE, 
                         probs = c(0.05, 0.5, 0.95), 
                         round = 2)
fyke_betas= fyke_betas[-1,]

gill_betas<-MCMCsummary(gill_post$Beta, 
                         Rhat = TRUE, 
                         n.eff = TRUE, 
                         probs = c(0.05, 0.5, 0.95), 
                         round = 2)
gill_betas= gill_betas[-1,]

seine_betas<-MCMCsummary(seine_post$Beta, 
                         Rhat = TRUE, 
                         n.eff = TRUE, 
                         probs = c(0.05, 0.5, 0.95), 
                         round = 2)
seine_betas= seine_betas[-1,]

#* Plot Posterior means and CIs for all parameters to plot ####
betaEsts <- matrix(NA, nrow=13,ncol=3)
colnames(betaEsts) <- c("estimate", "conf.low", "conf.high")

betaEsts[,1] <- seine_betas$mean
betaEsts[,2] <- seine_betas$`5%`
betaEsts[,3] <- seine_betas$`95%` 
betaEsts<-as.data.frame(betaEsts)

betaEsts$term <- c("lake_area", "geom_ratio", "order", "ws_urban", "ws_forest", "ws_agriculture",
                   "ws_shrub", "ws_wetland", "ws_slope", "ws_elevation", "dd_mean", "shoreline_houses", "Secchi")


#add colors for sig different than 0
betaEsts$color <- as.numeric(betaEsts[,2] * betaEsts[,3] > 0 )

# plot using the dotwhisker package 
lake_plot<-dwplot(betaEsts, style = "dotwhisker", 
                  dot_args = list(aes(colour = factor(color))),
                  whisker_args = list(aes(colour = factor(color))), 
                  vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("seine") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
lake_plot

# Plot response of species over the gradient of environmental variable x
Gradient = constructGradient(m_seine_l, focalVariable="z_dd_mean")
predY = predict(m_seine_l, Gradient=Gradient)
plotGradient(m_seine_l, Gradient, pred=predY, measure="Y", q = c(0.05, 0.5, 0.95), showData = TRUE) 

#*investigate MCMC convergence ####
plot(shock_post$Beta) #trace plot and density plot - two chains(red and black) look identical, the chains mix very well, i.e. they go fast up and down without any apparent autocorrelation, they seem to have reached a stationary distribution  
effectiveSize(shock_post$Beta) #did we set enough samples? effective sample sizes are very close to the theoretical value of the actual number of samples, which is 2000 (1000 per chain)
gelman.diag(shock_post$Beta,multivariate=FALSE)$psrf #scale reduction factors are very close to one, which indicates that two chains gave consistent results
par(mfrow=c(1,2)) #plot if you have a lot of parameters 
hist(effectiveSize(shock_post$Beta), main="ess(beta)") 
hist(gelman.diag(shock_post$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")

#* get RMSE and R2 by looking at the predicted dist distribution ####
#Lower values of RMSE indicate better fit.
#explanatory R2 
preds = computePredictedValues(m_shock_l) #posterior predictive distribution 
evaluateModelFit(hM=m_shock_l, predY=preds)

#diagnostic plots 
preds.mean = apply(preds, FUN=mean, MARGIN=1) #get the posterior means 
nres = scale(Y-preds.mean) #standardize the residuals
par(mfrow=c(1,2))
hist(nres, las = 1)
plot(preds.mean,nres, las = 1) 
abline(a=0,b=0)

#### hierarchical structure ####
#check ICC using lmer
m_region<-lmer(SHOCK ~ 1 + (1|FMU_Code), cpue_shock, REML=FALSE) #use maximum likelihood instead of resricted maximum likelihod
summary(m_region) # this gives you the between region variance as the intercept and the residual variance is the within region variance  
icc_n <- as.data.frame(VarCorr(m_region))[1,4]
icc_d <- as.data.frame(VarCorr(m_region))[1,4] + 
  as.data.frame(VarCorr(m_region))[2,4]
icc_n / icc_d 

# note this is how you would do it with max likelihood approach: lmer(y~x+(1|plot.id),data=da)
#still have Y and XData from above models 
#Y = as.matrix(log(cpue_shock$SHOCK))
#XData = data.frame(dat_net$lake_order,  dat_net$dd_mean, dat_net$secchi_m)
studyDesign = data.frame(sample = as.factor(cpue_shock$new_key), plot = as.factor(cpue_shock$FMU_Code)) #individual samples and grouping by FMU 
rL = HmscRandomLevel(units = studyDesign$plot) #default options for latent variables. this function is to reduce dimensionality of residual association. latent features should not exceed the number of species  
#construct model object
m = Hmsc(Y=Y, XData=XData, XFormula= ~z_lake_area + z_perim_km+ z_order+ z_ws_area_km2+ z_ws_urban+ z_ws_forest+ z_ws_agriculture+ z_ws_shrub+ z_ws_wetland+ z_ws_slope+ z_ws_elevation+ z_dd_mean+ z_houses+ z_secchi, distr = "normal",
         studyDesign=studyDesign, ranLevels=list("plot"=rL)) #list of study design that will be random

#fit the model 7mins
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient,
               nChains = nChains, verbose = verbose)
preds=computePredictedValues(m) 
MF=evaluateModelFit(hM=m, predY=preds)
MF #this R2 is explanatory and is misleading because of over fitting. better to evaluate with cross-validation predictive R2

#To look at the parameter estimates, we extract the posterior distribution
mpost = convertToCodaObject(m1) 
sum<-summary(mpost$Beta)
#how to pull these means to plot!? 
mean_values<-as.data.frame(sum[["statistics"]])
mean_values = mean_values[-1,]
quant_values<-as.data.frame(sum[["quantiles"]])
quant_values = quant_values[-1,]

#* Posterior means and CIs for all parameters to plot ####
betaEsts <- matrix(NA, nrow=14,ncol=3)
colnames(betaEsts) <- c("estimate", "conf.low", "conf.high")

betaEsts[,1] <- mean_values$Mean
betaEsts[,2] <- quant_values$`2.5%`
betaEsts[,3] <- quant_values$`97.5%` #this is if 95% credible intervals don't cross - so more than the 90%
betaEsts<-as.data.frame(betaEsts)

betaEsts$term <- c("lake_area", "perim", "order", "ws_area", "ws_urban", "ws_forest", "ws_agriculture",
                   "ws_shrub", "ws_wetland", "ws_slope", "ws_elevation", "dd_mean", "shoreline_houses", "Secchi")


#add colors for sig different than 0
betaEsts$color <- as.numeric(betaEsts[,2] * betaEsts[,3] > 0 )

# plot using the dotwhisker package 
library(dotwhisker)
lake_plot<-dwplot(betaEsts, style = "dotwhisker", 
                  dot_args = list(aes(colour = factor(color))),
                  whisker_args = list(aes(colour = factor(color))), 
                  vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("lakes") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
lake_plot

#investigate MCMC convergence 
plot(mpost$Beta) #trace plot and density plot - two chains(red and black) look identical, the chains mix very well, i.e. they go fast up and down without any apparent autocorrelation, they seem to have reached a stationary distribution  
effectiveSize(mpost$Beta) #did we set enough samples? effective sample sizes are very close to the theoretical value of the actual number of samples, which is 2000 (1000 per chain)
gelman.diag(mpost$Beta,multivariate=FALSE)$psrf #scale reduction factors are very close to one, which indicates that two chains gave consistent results

#cross-validation 
#2-fold cross-validation - randomly sample into 2 groups, model fitting and predictions are made separately for each fold 
partition=createPartition(m, nfolds=2, column="sample") #this assigns a partition based on the samples 
preds=computePredictedValues(m, partition = partition)
MF=evaluateModelFit(hM=m, predY=preds)
MF #predictive R2 is lower than the explanatory R2 above 



#### spatial structure: including coordinates instead of grouping by Fisheries Management Unit #### 

#fit model 
##As the model includes a spatially structured random effect, its predictive power is based on both the fixed and the random effects. 
#because of random effect, the model can utilize observed data from nearby sampling units included in model fitting when predicting the response for a focal sampling unit that is not included in model fitting
studyDesign = data.frame(sample = as.factor(dat_pres$new_key) )
xycoords<-dplyr::select(dat_pres, LONG_DD, LAT_DD)
qu25 <- xycoords %>% 
  ungroup() #dplyr add a grouping value if it has been used previously
xycoords<-select(qu25, -c(new_key))
rownames(xycoords) = dat_pres$new_key
rL = HmscRandomLevel(sData = xycoords, longlat= TRUE) #use the sData argument to construct the random effect for spatially explicit
m = Hmsc(Y=Y, XData=XData, XFormula= ~dat_pres.lake_area_ac + dat_pres.dd_mean, 
         studyDesign=studyDesign, ranLevels=list("sample"=rL))
# Set timer 
start.time = Sys.time()  
#run model 
m = sampleMcmc(m, thin = thin, samples = samples, transient = transient, 
               nChains = nChains, verbose = verbose)

#calculate run time 
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

#To look at the parameter estimates, we extract the posterior distribution
mpost = convertToCodaObject(m) 
summary(mpost$Beta)

#investigate MCMC convergence 
plot(mpost$Beta) #trace plot and density plot - two chains(red and black) look identical, the chains mix very well, i.e. they go fast up and down without any apparent autocorrelation, they seem to have reached a stationary distribution  
effectiveSize(mpost$Beta) #did we set enough samples? effective sample sizes are very close to the theoretical value of the actual number of samples, which is 2000 (1000 per chain)
gelman.diag(mpost$Beta,multivariate=FALSE)$psrf #scale reduction factors are very close to one, which indicates that two chains gave consistent results

#evaluate model fit 
preds = computePredictedValues(m)
MF = evaluateModelFit(hM=m, predY=preds) 
MF

#diagnostic plots 
preds.mean = apply(preds, FUN=mean, MARGIN=1) #get the posterior means 
nres = scale(Y-preds.mean) #standardize the residuals
par(mfrow=c(1,2))
hist(nres, las = 1)
plot(preds.mean,nres, las = 1) 
abline(a=0,b=0)

##cross-val - this takes a lot of time 
partition = createPartition(m, nfolds = 2, column = "sample") 
preds = computePredictedValues(m, partition=partition)
MF = evaluateModelFit(hM=m, predY=preds) 
MF$R2

#look at estimates of Beta parameters 
postBeta = getPostEstimate(m, parName = "Beta")
plotBeta(m, post = postBeta, param = "Support", supportLevel = 0.95, colors = colorRampPalette(c("darkgreen", "white", "purple"))) #95% credible interval

#### PREDICTION TO UNSAMPLED LAKES #### 
predict.Hmsc 