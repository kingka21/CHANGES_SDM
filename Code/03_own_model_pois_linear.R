# modeling linear model (with varying intercepts only (no grouping factor)
#written by Katelyn King 
#code adapted from Tyler Wagner 
#using catch as response with an offset value because you can use Poisson distribution 
#also should look into zero-inflated 

# rm(list=ls())
library(R2jags)
library(lme4)
library(MCMCpack)
library(dotwhisker) #https://cran.r-project.org/web/packages/dotwhisker/vignettes/dotwhisker-vignette.html
library(dplyr)

#### standardize and transform predictor values ####
lmb_count_dat<-lmb_count_dat %>%
  ungroup()

#transform 
lmb_count_dat$lake_area_km2<-lmb_count_dat$lake_area_ac/247 #change from acre to km2

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
#seine 
dat_seine$z_order<-as.numeric(scale(log(dat_seine$lake_order+0.001))) # min is 0 because these are the isolated
dat_seine$z_lake_area<-as.numeric(scale(log(dat_seine$lake_area_km2))) #good 
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

# select variables for analysis 
dat_shock<-dplyr::select(dat_shock, new_key, FMU_Code, fish_count_SHOCK, effort_SHOCK, z_lake_area, z_perim_km, z_order, z_ws_area_km2, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation, z_dd_mean, z_houses, z_secchi) 
dat_fyke<-dplyr::select(dat_fyke, new_key, FMU_Code, fish_count_FT_NET, effort_FT_NET, z_lake_area, z_perim_km, z_order, z_ws_area_km2, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation, z_dd_mean, z_houses, z_secchi) 
dat_gill<-dplyr::select(dat_gill, new_key, FMU_Code, fish_count_GILL, effort_GILL, z_lake_area, z_perim_km, z_order, z_ws_area_km2, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation, z_dd_mean, z_houses, z_secchi) 
dat_seine<-dplyr::select(dat_seine, new_key, FMU_Code, fish_count_SEINE, effort_SEINE, z_lake_area, z_perim_km, z_order, z_ws_area_km2, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation, z_dd_mean, z_houses, z_secchi) 

#view response - sort of normal. could try zero inflation model and compare 
hist(dat_shock$fish_count_SHOCK)
hist(log(dat_shock$fish_count_SHOCK))
hist(dat_fyke$fish_count_FT_NET)
hist(log(dat_fyke$fish_count_FT_NET))
hist(dat_gill$fish_count_GILL)
hist(log(dat_gill$fish_count_GILL))
hist(dat_seine$fish_count_SEINE)
hist(log(dat_seine$fish_count_SEINE))

# Remove observations with no effort
#dat_shock <- dat_shock[dat_shock$effort_SHOCK != 0 ,]
#dim(dat_shock)
#dat_fyke <- dat_fyke[dat_fyke$effort_FT_NET != 0 ,]
#dim(dat_fyke)
#dat_gill <- dat_gill[dat_gill$effort_GILL != 0 ,]
#dim(dat_gill)
#dat_seine <- dat_seine[dat_seine$effort_SEINE != 0 ,]
#dim(dat_seine)


# Create offset term
dat_shock$offset_shock <- log(dat_shock$effort_SHOCK)
dat_fyke$offset_fyke <- log(dat_fyke$effort_FT_NET)
dat_gill$offset_gill<-log(dat_gill$effort_GILL)
dat_seine$offset_seine<-log(dat_seine$effort_SEINE)

#################################################################
########## BUGS CODE ############################################
#################################################################
#linear poisson model with offset for effort  
sink("model.txt")
cat("
    model {
    # Likelihood: 
    for (i in 1:n){ 
    y[i] ~ dpois(lambda[i])  # Distribution 
    log(lambda[i]) <- log.lambda[i] # log-link 
    log.lambda[i] <- 1 * offset[i] + alpha + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i] + b[4] * x4[i] + 
                       b[5] * x5[i] + b[6] * x6[i] + b[7] * x7[i] + b[8] * x8[i] + b[9] * x9[i] + b[10] * x10[i] +
                      b[11] * x11[i] + b[12] * x12[i]  + b[13] * x13[i] + b[14] * x14[i] 
    } 
    
    
    # Priors
    alpha ~ dnorm(0, 0.0001)
    
    # Level-1 predictors #  horseshoe prior 
    for(k in 1:14){
    b[k] ~ dnorm(0,0.0001)
    }
    
    } # end model
    ",fill = TRUE)
sink()

### set parameters 
# Initial values
inits <- function (){
  list (alpha = rnorm(1,0,1) )
}


# Parameters monitored
parameters <- c("alpha", "b")


# MCMC settings
ni <- 10000 # number of iterations 
nt <- 1 #number to thin
nb <- 5000 #number to burn  
nc <- 3 #number of chains


#### model nets and shock separately #### 

#*shock####

# Load data
data <- list(y = dat_shock$fish_count_SHOCK, n = dim(dat_shock)[1], 
             x1=dat_shock$z_secchi, x2=dat_shock$z_lake_area, x3=dat_shock$z_perim_km, x4=dat_shock$z_order, x5=dat_shock$z_ws_area_km2, x6=dat_shock$z_ws_urban, x7=dat_shock$z_ws_forest, x8=dat_shock$z_ws_agriculture, x9=dat_shock$z_ws_shrub, x10=dat_shock$z_ws_wetland, 
             x11=dat_shock$z_ws_slope, x12=dat_shock$z_ws_elevation, x13=dat_shock$z_dd_mean, x14=dat_shock$z_houses,
             offset = dat_shock$offset_shock
) 

#initial values, parameters to monitor, and MCMC parameters are the same as above 
# Set timer 
start.time = Sys.time()         

# Call JAGS from R and run model 
outShock <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
                 n.thin = nt, n.iter = ni, n.burnin = nb)

# # Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

# Summarize posteriors
print(outShock, dig = 3)

# traceplot(outShock) ???

#*fyke and trap nets####
# Number of FMUs might be different. set these parameters again 
J <- length(unique(dat_fyke$FMU_Code))

# Create group index, must go from 1,..J
dat_fyke$FMU_Code <- droplevels(as.factor(dat_fyke$FMU_Code))
group <- as.numeric(as.factor(dat_fyke$FMU_Code))


# Load data
data <- list(y = dat_fyke$fish_count_FT_NET, group = group, n = dim(dat_fyke)[1], J = J,
             x1=dat_fyke$z_secchi, x2=dat_fyke$z_lake_area, x3=dat_fyke$z_perim_km, x4=dat_fyke$z_order, x5=dat_fyke$z_ws_area_km2, x6=dat_fyke$z_ws_urban, x7=dat_fyke$z_ws_forest, x8=dat_fyke$z_ws_agriculture, x9=dat_fyke$z_ws_shrub, x10=dat_fyke$z_ws_wetland, 
             x11=dat_fyke$z_ws_slope, x12=dat_fyke$z_ws_elevation, x13=dat_fyke$z_dd_mean, x14=dat_fyke$z_houses,
             offset = dat_fyke$offset_fyke
) 

#initial values, parameters to monitor, and MCMC parameters are the same as above 
# Set timer 
start.time = Sys.time()         

# Call JAGS from R and run model 
outNET <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
               n.thin = nt, n.iter = ni, n.burnin = nb)

# # Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

# Summarize posteriors
print(outNET, dig = 3)

#*gill nets####
# Number of FMUs might be different. set these parameters again 
J <- length(unique(dat_gill$FMU_Code))

# Create group index, must go from 1,..J
dat_gill$FMU_Code <- droplevels(as.factor(dat_gill$FMU_Code))
group <- as.numeric(as.factor(dat_gill$FMU_Code))


# Load data
data <- list(y = dat_gill$fish_count_GILL, group = group, n = dim(dat_gill)[1], J = J,
             x1=dat_gill$z_secchi, x2=dat_gill$z_lake_area, x3=dat_gill$z_perim_km, x4=dat_gill$z_order, x5=dat_gill$z_ws_area_km2, x6=dat_gill$z_ws_urban, x7=dat_gill$z_ws_forest, x8=dat_gill$z_ws_agriculture, x9=dat_gill$z_ws_shrub, x10=dat_gill$z_ws_wetland, 
             x11=dat_gill$z_ws_slope, x12=dat_gill$z_ws_elevation, x13=dat_gill$z_dd_mean, x14=dat_gill$z_houses,
             offset = dat_gill$offset_gill
) 

#initial values, parameters to monitor, and MCMC parameters are the same as above 
# Set timer 
start.time = Sys.time()         

# Call JAGS from R and run model 
outGILL <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
                n.thin = nt, n.iter = ni, n.burnin = nb)

# # Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

# Summarize posteriors
print(outGILL, dig = 3)

#*seine ####
# Number of FMUs might be different. set these parameters again 
J <- length(unique(dat_seine$FMU_Code))

# Create group index, must go from 1,..J
dat_seine$FMU_Code <- droplevels(as.factor(dat_seine$FMU_Code))
group <- as.numeric(as.factor(dat_seine$FMU_Code))


# Load data
data <- list(y = dat_seine$fish_count_SEINE, group = group, n = dim(dat_seine)[1], J = J,
             x1=dat_seine$z_secchi, x2=dat_seine$z_lake_area, x3=dat_seine$z_perim_km, x4=dat_seine$z_order, x5=dat_seine$z_ws_area_km2, x6=dat_seine$z_ws_urban, x7=dat_seine$z_ws_forest, x8=dat_seine$z_ws_agriculture, x9=dat_seine$z_ws_shrub, x10=dat_seine$z_ws_wetland, 
             x11=dat_seine$z_ws_slope, x12=dat_seine$z_ws_elevation, x13=dat_seine$z_dd_mean, x14=dat_seine$z_houses,
             offset = dat_seine$offset_seine
) 

#initial values, parameters to monitor, and MCMC parameters are the same as above 
# Set timer 
start.time = Sys.time()         

# Call JAGS from R and run model 
outSEINE <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
                 n.thin = nt, n.iter = ni, n.burnin = nb)

# # Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

# Summarize posteriors
print(outSEINE, dig = 3)

########### PLOT ####################################

#* plot for shock #### 
# Posterior means and CIs for all parameters
shock_betaEsts <- matrix(NA, nrow=14,ncol=3)
for(i in 1:14){
  shock_betaEsts[i,1] <- mean(quantile(outShock$BUGSoutput$sims.list$b[,i],c(0.05,0.95)))
  shock_betaEsts[i, 2:3] <- quantile(outShock$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
shock_betaEsts<-as.data.frame(shock_betaEsts)


colnames(shock_betaEsts) <- c("estimate", "conf.low", "conf.high")
shock_betaEsts$term <- c("Secchi", "lake_area", "perim", "order", "ws_area", "ws_urban", "ws_forest", "ws_agriculture",
                         "ws_shrub", "ws_wetland", "ws_slope", "ws_elevation", "dd_mean", "shoreline_houses"
)

#add colors for sig different than 0
shock_betaEsts$color <- as.numeric(shock_betaEsts[,2] * shock_betaEsts[,3] > 0 )

# plot using the dotwhisker package 

shock_plot<-dwplot(shock_betaEsts, style = "dotwhisker", 
                   dot_args = list(aes(colour = factor(color))),
                   whisker_args = list(aes(colour = factor(color))), 
                   vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("shock") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
shock_plot
# theme(plot.title = element_text(face="bold"),
#    legend.position = c(0.007, 0.01),
#   legend.justification = c(0, 0), 
#  legend.background = element_rect(colour="grey80"),
# legend.title = element_blank()) 

#* plot for nets #### 
# Posterior means and CIs for all parameters
net_betaEsts <- matrix(NA, nrow=14,ncol=3)
for(i in 1:14){
  net_betaEsts[i,1] <- mean(quantile(outNET$BUGSoutput$sims.list$b[,i],c(0.05,0.95)))
  net_betaEsts[i, 2:3] <- quantile(outNET$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
net_betaEsts<-as.data.frame(net_betaEsts)


colnames(net_betaEsts) <- c("estimate", "conf.low", "conf.high")
net_betaEsts$term <- c("Secchi", "lake_area", "perim", "order", "ws_area", "ws_urban", "ws_forest", "ws_agriculture",
                       "ws_shrub", "ws_wetland", "ws_slope", "ws_elevation", "dd_mean", "shoreline_houses"
)

#add colors for sig different than 0
net_betaEsts$color <- as.numeric(net_betaEsts[,2] * net_betaEsts[,3] > 0 )

# plot using the dotwhisker package 

net_plot<-dwplot(net_betaEsts, style = "dotwhisker", 
                 dot_args = list(aes(colour = factor(color))),
                 whisker_args = list(aes(colour = factor(color))), 
                 vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("net") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
net_plot

#* plot for gill #### 
# Posterior means and CIs for all parameters
net_betaEsts <- matrix(NA, nrow=14,ncol=3)
for(i in 1:14){
  net_betaEsts[i,1] <- mean(quantile(outGILL$BUGSoutput$sims.list$b[,i],c(0.05,0.95)))
  net_betaEsts[i, 2:3] <- quantile(outGILL$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
net_betaEsts<-as.data.frame(net_betaEsts)


colnames(net_betaEsts) <- c("estimate", "conf.low", "conf.high")
net_betaEsts$term <- c("Secchi", "lake_area", "perim", "order", "ws_area", "ws_urban", "ws_forest", "ws_agriculture",
                       "ws_shrub", "ws_wetland", "ws_slope", "ws_elevation", "dd_mean", "shoreline_houses"
)

#add colors for sig different than 0
net_betaEsts$color <- as.numeric(net_betaEsts[,2] * net_betaEsts[,3] > 0 )

# plot using the dotwhisker package 

gill_plot<-dwplot(net_betaEsts, style = "dotwhisker", 
                  dot_args = list(aes(colour = factor(color))),
                  whisker_args = list(aes(colour = factor(color))), 
                  vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("gill") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
gill_plot

#* plot for seine #### 
# Posterior means and CIs for all parameters
net_betaEsts <- matrix(NA, nrow=14,ncol=3)
for(i in 1:14){
  net_betaEsts[i,1] <- mean(quantile(outSEINE$BUGSoutput$sims.list$b[,i],c(0.05,0.95)))
  net_betaEsts[i, 2:3] <- quantile(outSEINE$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
net_betaEsts<-as.data.frame(net_betaEsts)


colnames(net_betaEsts) <- c("estimate", "conf.low", "conf.high")
net_betaEsts$term <- c("Secchi", "lake_area", "perim", "order", "ws_area", "ws_urban", "ws_forest", "ws_agriculture",
                       "ws_shrub", "ws_wetland", "ws_slope", "ws_elevation", "dd_mean", "shoreline_houses"
)

#add colors for sig different than 0
net_betaEsts$color <- as.numeric(net_betaEsts[,2] * net_betaEsts[,3] > 0 )

# plot using the dotwhisker package 

seine_plot<-dwplot(net_betaEsts, style = "dotwhisker", 
                   dot_args = list(aes(colour = factor(color))),
                   whisker_args = list(aes(colour = factor(color))), 
                   vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("seine") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
seine_plot

Fig1<-cowplot::plot_grid(net_plot, shock_plot)
cowplot::save_plot("Fig1.png", Fig1, base_width = 8,
                   base_height = 4, dpi=300)

