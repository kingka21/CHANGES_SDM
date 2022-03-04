#HMSC: Hierarchical Modelling of Species Communities (HMSC) ####

#final drivers of Hansen et al 2017 
#walleye- degree days, lake area, conductivity, and shoreline complexity 
#note no seine catches for walleye 

# Walleye model 
library(Hmsc)
library(MCMCvis)
library(dotwhisker)

#### standardize and transform predictor values ####
wae_count_dat<-wae_count_dat %>%
  ungroup()

#transform 
wae_count_dat$lake_area_km2<-wae_count_dat$lake_area_ac/247 #change from acre to km2

##subset into three different gear response variables and remove 0s, select subset of metrics not correlated for analysis and remove NAs #
dat_shock<-filter(wae_count_dat, fish_count_SHOCK >0) %>%
  na.omit()
dat_fyke<-filter(wae_count_dat, fish_count_FT_NET >0) %>%
  na.omit() 
dat_gill<-filter(wae_count_dat, fish_count_GILL >0) %>% 
  na.omit() 

#your model works on the representation given by its input predictors, differently scaled test/prediction data would potentially lead to over/under-exaggeration of a feature.
#standardize and transform 
#shock 
dat_shock$z_order<-as.numeric(scale(log(dat_shock$lake_order+0.001))) # min is 0 because these are the isolated
dat_shock$z_lake_area<-as.numeric(scale(log(dat_shock$lake_area_km2))) #good 
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

#fyke
dat_fyke$z_order<-as.numeric(scale(log(dat_fyke$lake_order+0.001))) # min is 0 because these are the isolated
dat_fyke$z_lake_area<-as.numeric(scale(log(dat_fyke$lake_area_km2))) #good 
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

#gill 
dat_gill$z_order<-as.numeric(scale(log(dat_gill$lake_order+0.001))) # min is 0 because these are the isolated
dat_gill$z_lake_area<-as.numeric(scale(log(dat_gill$lake_area_km2))) #good 
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

#### Bays linear model #### 
#To fit the HMSC model with Bayesian inference, we use the sampleMcmc function.
nChains = 3  #how many chains to sample (nChains)
thin=1 #how many samples to keep 
samples = 10000 #how many samples to obtain per chain (samples)
transient=5000 #how long transient (also called burn-in) 
verbose = 5000 #how frequently we wish to see the progress of the MCMC sampling (verbose). 

#shock model 
Y = as.matrix(log(dat_shock$SHOCK))
XData = dplyr::select(dat_shock, z_lake_area , z_perim_km, z_order, z_ws_area_km2, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation, z_dd_mean, z_houses, z_secchi)
#XData = data.frame(cpue_shock$z_lake_area, cpue_shock$z_perim_km, cpue_shock$z_order, cpue_shock$z_ws_area_km2, cpue_shock$z_ws_urban, cpue_shock$z_ws_forest, cpue_shock$z_ws_agriculture, cpue_shock$z_ws_shrub, cpue_shock$z_ws_wetland, cpue_shock$z_ws_slope, cpue_shock$z_ws_elevation, cpue_shock$z_dd_mean, cpue_shock$z_houses, cpue_shock$z_secchi)
m_shock = Hmsc(Y = Y, XData = XData, XFormula = ~z_lake_area + z_perim_km+ z_order+ z_ws_area_km2+ z_ws_urban+ z_ws_forest+ z_ws_agriculture+ z_ws_shrub+ z_ws_wetland+ z_ws_slope+ z_ws_elevation+ z_dd_mean+ z_houses+ z_secchi, 
               XScale = FALSE, distr = "normal") # this constructs the model object # normal dist for continuous data #alrady scaled my parameters 
#fit the model
m_shock_l = sampleMcmc(m_shock, thin = thin, samples = samples, transient = transient, nChains = nChains, verbose = verbose)

#fyke model 
Y = as.matrix(log(dat_fyke$FT_NET))
XData = dplyr::select(dat_fyke, z_lake_area , z_perim_km, z_order, z_ws_area_km2, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation, z_dd_mean, z_houses, z_secchi)
#XData = data.frame(cpue_shock$z_lake_area, cpue_shock$z_perim_km, cpue_shock$z_order, cpue_shock$z_ws_area_km2, cpue_shock$z_ws_urban, cpue_shock$z_ws_forest, cpue_shock$z_ws_agriculture, cpue_shock$z_ws_shrub, cpue_shock$z_ws_wetland, cpue_shock$z_ws_slope, cpue_shock$z_ws_elevation, cpue_shock$z_dd_mean, cpue_shock$z_houses, cpue_shock$z_secchi)
m_fyke = Hmsc(Y = Y, XData = XData, XFormula = ~z_lake_area + z_perim_km+ z_order+ z_ws_area_km2+ z_ws_urban+ z_ws_forest+ z_ws_agriculture+ z_ws_shrub+ z_ws_wetland+ z_ws_slope+ z_ws_elevation+ z_dd_mean+ z_houses+ z_secchi, 
              XScale = FALSE, distr = "normal") # this constructs the model object # normal dist for continuous data #alrady scaled my parameters 
#fit the model
m_fyke_l = sampleMcmc(m_fyke, thin = thin, samples = samples, transient = transient, nChains = nChains, verbose = verbose)

#gill model 
Y = as.matrix(log(dat_gill$GILL))
XData = dplyr::select(dat_gill, z_lake_area , z_perim_km, z_order, z_ws_area_km2, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation, z_dd_mean, z_houses, z_secchi)
#XData = data.frame(cpue_shock$z_lake_area, cpue_shock$z_perim_km, cpue_shock$z_order, cpue_shock$z_ws_area_km2, cpue_shock$z_ws_urban, cpue_shock$z_ws_forest, cpue_shock$z_ws_agriculture, cpue_shock$z_ws_shrub, cpue_shock$z_ws_wetland, cpue_shock$z_ws_slope, cpue_shock$z_ws_elevation, cpue_shock$z_dd_mean, cpue_shock$z_houses, cpue_shock$z_secchi)
m_gill = Hmsc(Y = Y, XData = XData, XFormula = ~z_lake_area + z_perim_km+ z_order+ z_ws_area_km2+ z_ws_urban+ z_ws_forest+ z_ws_agriculture+ z_ws_shrub+ z_ws_wetland+ z_ws_slope+ z_ws_elevation+ z_dd_mean+ z_houses+ z_secchi, 
              XScale = FALSE, distr = "normal") # this constructs the model object # normal dist for continuous data #alrady scaled my parameters 
#fit the model
m_gill_l = sampleMcmc(m_gill, thin = thin, samples = samples, transient = transient, nChains = nChains, verbose = verbose)


#### investigating output ####
# extract the posterior distribution
shock_post = convertToCodaObject(m_shock_l) 
fyke_post = convertToCodaObject(m_fyke_l) 
gill_post = convertToCodaObject(m_gill_l) 

#pull out means and 95% credible intervals 
shock_betas<-MCMCsummary(shock_post$Beta) #, 
                       #  Rhat = TRUE, 
                        # n.eff = TRUE, 
                        # probs = c(0.1, 0.5, 0.9), # if you want to change to 90% CIs
                        # round = 2)
shock_betas= shock_betas[-1,]

fyke_betas<-MCMCsummary(fyke_post$Beta) #, 
                      #  Rhat = TRUE, 
                      #  n.eff = TRUE, 
                      #  probs = c(0.1, 0.5, 0.9), 
                      #  round = 2)
fyke_betas= fyke_betas[-1,]

gill_betas<-MCMCsummary(gill_post$Beta , 
                       Rhat = TRUE, 
                        n.eff = TRUE, 
                        probs = c(0.05, 0.5, 0.95), 
                        round = 2)
gill_betas= gill_betas[-1,]

#* Plot Posterior means and CIs for all parameters to plot ####
betaEsts <- matrix(NA, nrow=14,ncol=3)
colnames(betaEsts) <- c("estimate", "conf.low", "conf.high")

betaEsts[,1] <- gill_betas$mean
betaEsts[,2] <- gill_betas$`5%`
betaEsts[,3] <- gill_betas$`95%` 
betaEsts<-as.data.frame(betaEsts)

betaEsts$term <- c("lake_area", "perim", "order", "ws_area", "ws_urban", "ws_forest", "ws_agriculture",
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
  ggtitle("gill") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
lake_plot
