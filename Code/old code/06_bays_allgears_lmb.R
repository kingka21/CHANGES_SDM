# modeling hierarchical model with varying intercepts by FMU and all gears combined
#written by Katelyn King and Henrique Giacomini
#code adapted from Tyler Wagner 
#separating catches from different gears
#removed the lake grouping 
#final drivers of Hansen et al 2017 
#lmb - degree days, lake order, and Secchi depth - 3/4 of our obs have Secchi 

#load libraries 
library(R2jags)
library(lme4)
library(MCMCpack)
library(dotwhisker) #https://cran.r-project.org/web/packages/dotwhisker/vignettes/dotwhisker-vignette.html
library(dplyr)
library(ggmcmc) # package for analyzing output 
library(dotwhisker)

#### data ####
#need to keep all the gears used in a lake even if that gear did not catch sp of interest 

#lmb_dat_for_model has count by gear and effort by gear
#join tables for LMB 
#dat<- lmb_dat_for_model %>% # 472 unique lakes, but 41 lakes don't match the drivers with the nhdid
 # left_join(dplyr::select(lake_ll, new_key, LONG_DD, LAT_DD, FMU_Code))

dat<-read.csv("Data/lmb_dat_for_model_mar17.csv") 

dd_temp<-dd_temp %>%
  group_by(IHDLKID) %>%
  slice_max(HECTARES) %>% 
  rename(nhdid = IHDLKID)
dd_temp_long<-pivot_longer(dd_temp, 
                           cols= starts_with("DD_"),
                           names_to = "year", 
                           names_prefix = "DD_",
                           values_to = "dd_year")
dd_temp_long$year<-as.integer(as.character(dd_temp_long$year))
dat<-left_join(dat, dd_temp_long, by = c("nhdid", "year"))
#dat<-read.csv("Data/lmb_model_data_feb10.csv") %>%
#  left_join(temp_do_no_dups)

## standardize and transform predictor values ##
dat$z_lake_area<-as.numeric(scale(log(dat$lake_area_m2))) #good 
dat$z_max_depth<-as.numeric(scale(log(dat$maxdepth_m))) #good 
dat$z_dd_mean<-as.numeric(scale(log(dat$dd_year))) #good
dat$z_surface_temp_mean<-as.numeric(scale(log(dat$surface_temp_mean))) #good 
dat$z_secchi<-as.numeric(scale(log(dat$secchi_m)))
dat$z_bottom_do<-as.numeric(scale(log(dat$bottom_do_mgl+ 0.001))) ##has 0s so added 0.001 
dat$z_julian<-as.numeric(scale(dat$julian)) ##standardize julian date 
dat$z_ws_agriculture<-as.numeric(scale(asin(sqrt(dat$ws_agriculture_prop)))) #ok
dat$z_ws_forest<-as.numeric(scale(asin(sqrt(dat$ws_forest_prop)))) #good
dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat$ws_wetland_prop)))) #good
dat$z_elevation<-as.numeric(scale(log(dat$ws_mean_elevation_m)))
dat$logeffort <- log(dat$effort_new)

#check out response: log of count is mostly normally dist, so Poisson seems appropriate 
#hist(dat$fish_count_new)
#hist(log(dat$fish_count_new))

dat<-dplyr::select(dat, new_key, fish_count_new, logeffort, gear2, FMU_Code, 
                   z_lake_area, z_max_depth, z_secchi, z_julian, 
                   z_bottom_do, z_surface_temp_mean,
                   z_ws_forest, z_ws_wetland, pike_pres, wae_pres) %>%
            na.omit()

n_distinct(dat$new_key) #201 lakes 
#check for collinearity 
my_data <- dat[, c(7:15)] 
PerformanceAnalytics::chart.Correlation(my_data, histogram=TRUE, pch=19)
#lake area and watershed - keep lake area 
#max depth and bottom temp - keep max depth 
#urban and forest - keep forest. 
#note that surface temp and DD are only 0.5 corr 
#down to 12 variables 
#dat<-dat%>% dplyr::select(-c(z_ws_area_km2, z_ws_urban, z_bottom_temp, z_ws_shrub))

#### intercept only model to chec ICC ####
#determine if there is significant region-scale variation 
m_region<-lmer(log(fish_count_new+1) ~ 1 + (1|FMU_Code), dat, REML=FALSE) #use maximum likelihood instead of resricted maximum likelihod
summary(m_region) # the between region variance is the sd of the grouping variable and the 'residual variance' is the within region variance  
icc_n <- as.data.frame(VarCorr(m_region))[1,4]
icc_d <- as.data.frame(VarCorr(m_region))[1,4] + 
  as.data.frame(VarCorr(m_region))[2,4]
icc_n / icc_d #variance among regions 19%  

#################################################################
########## BUGS CODE ############################################
#################################################################

#hierarchical model without the lake groups (with varying intercepts by fmu)
sink("model.txt")
cat("
    model {
    # Likelihood: 
    
    for (i in 1:n){         
    y[i] ~ dpois(lambda[i])  # Distribution Poisson  
    log(lambda[i]) <- log.lambda[i] # log-link 
    log.lambda[i] <- alpha[group[i]] + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i] + b[4] * x4[i] + 
                       b[5] * x5[i] + b[6] * x6[i] + b[7] * x7[i] + b[8] * x8[i] + b[9]*x9[i] + b[10]*x10[i] +
                      logq[gear[i]]*IND[i] + logeffort[i]
         
    } 
    
     # Level-2 of the model: group (fisheries management unit regions) 
    for(j in 1:J){
    alpha[j] ~ dnorm(mu.alpha,tau.alpha)
    }
    
    # Priors
    mu.alpha ~ dnorm(0, 0.0001)
    sigma.alpha ~ dunif(0,100)
    
    #derived quantities 
    tau.alpha <- pow(sigma.alpha,-2)
    
    # priors for predictors 
    for(k in 1:10){
    b[k] ~ dnorm(0,0.0000001)
    }
    
     # priors for gear-specific log-catchabilities 
    for (k in 1:ngears) { 
    logq[k] ~ dnorm(0,0.0000001)   #non-informative
    } 

    } # end model
    ",fill = TRUE)
sink()


### set parameters ####
# Initial values assign for random variables (FMU there are 8 regions) 
inits <- function (){
  list (alpha = rnorm(8))
}


# Parameters monitored
parameters <- c("alpha", "mu.alpha","tau.alpha", "b", 'logq')


# MCMC settings
ni <- 10000 # number of iterations 
nt <- 1 #number to thin
nb <- 5000 #number to burn  
nc <- 3 #number of chains


#### set up data #### 


# Set the number of FMUs 
J <- length(unique(dat$FMU_Code))

# Create group index, must go from 1,..J
dat$FMU_Code <- droplevels(as.factor(dat$FMU_Code))
dat$group <- as.numeric(as.factor(dat$FMU_Code))

#set number of different sampling gears used and set a gear index 
ngears<-length(unique(dat$gear2))
gear <- as.numeric(as.factor(dat$gear2))

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
dat$IND<-ifelse(dat$gear2 == "FT_NET", 0, 1) #FT_NET is the reference gear

# Load data
data <- list(y = dat$fish_count_new, group = dat$group, gear=gear, n = dim(dat)[1], J = J, ngears = ngears,
             x1=dat$z_secchi, x2=dat$z_lake_area, x3=dat$z_surface_temp_mean, x4=dat$z_max_depth, x5=dat$z_bottom_do,
             x6=dat$z_ws_forest,x7=dat$z_ws_wetland, x8=dat$z_julian, x9=dat$wae_pres, x10=dat$pike_pres,
             logeffort=dat$logeffort, IND=dat$IND
) 

### run the model using JAGS and R ####
# Set timer 
start.time = Sys.time()         

# Call JAGS from R and run model 
output_2level<- jags(data, inits, parameters, "model.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb)


# # Calculate computation time (takes about 12 mins!) 
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

# Summarize posteriors
print(output_2level, dig = 3)
#save output 
#jag.sum<-output$BUGSoutput$summary
#write.table(x=jag.sum,file="out.txt",sep="\t")

### check convergence ####
#in this case, logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
#check convergence with Brooks-Gelman-Rubin statistic
which(output_2level$BUGSoutput$summary[, c("Rhat")] > 1.1)

which(output_2level$BUGSoutput$summary[, c("n.eff")] < 0.1*15000) #times number of retained draws 
which(output_2level$BUGSoutput$summary[, c("n.eff")] < 500)

#look at trace plots 
#R2jags::traceplot(output)
#saveRDS(output, "Data/output/model_output_lmb.rds") 

#read in output 
#output <- readRDS("Data/output/model_output.rds") 
#density plots
# use as.mcmmc to convert rjags object into mcmc.list - plots in coda
jagsfit.mcmc <- as.mcmc(output_2level)
#require(lattice)
#densityplot(jagsfit.mcmc)

#the ggs function transforms the mcmc output into a longformat tibble, that we can use to make different types of plots.
model1tranformed <- ggs(jagsfit.mcmc) 
summary(model1tranformed$Parameter)

#traceplots (caterpillars)
ggs_traceplot(model1tranformed, family = "alpha")
ggs_traceplot(model1tranformed, family = "b")
ggs_traceplot(model1tranformed, family = "logq")

# posterior density plots 
#ggplot(filter(model1tranformed,Parameter == "b[1]", Iteration > 1000),aes(x = value))+
 # geom_density(fill  = "yellow", alpha = .5)+
  #geom_vline(xintercept = 0, col  = "red", size = 1)+ 
  #theme_light() +
  #labs(title = "Posterior Density of b1 secchi")


#### effect plots #### 


#betas
bEst <- matrix(NA, nrow=10,ncol=3)
for(i in 1:10){ #parameters
  bEst[i,1] <- mean(output_2level$BUGSoutput$sims.list$b[,i])
  bEst[i,2:3] <- quantile(output_2level$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
bEst<-as.data.frame(bEst)
bEst$variable<-c("Secchi", "lake area", "mean_surface_temp", "max_depth", "bottom_do", 
                 "ws forest", "ws wetland",  "julian day", "pres_wae", "pres_pike")

#logqs
##qEst <- matrix(NA, nrow=4,ncol=3)
#for(i in 1:4){ #gear catchability without the reference gear 
 # qEst[i,1] <- mean(output$BUGSoutput$sims.list$logq[,i])
  #qEst[i,2:3] <- quantile(output$BUGSoutput$sims.list$logq[,i],c(0.05,0.95))
#}
#qEst<-as.data.frame(qEst)
#qEst$variable<-c("fyke_ref", "gill", "seine", "shock")
#qEst<- filter(qEst, variable != "fyke_ref") #remove reference

#combine datasets 
#EstsLake <- gtools::smartbind(bEst,qEst)
#colnames(EstsLake)<- c("estimate", "conf.low", "conf.high", "term")
colnames(bEst)<- c("estimate", "conf.low", "conf.high", "term")

#add colors for sig different than 0
#EstsLake$color <- as.numeric(EstsLake[,2] * EstsLake[,3] > 0 )
bEst$color <- as.numeric(bEst[,2] * bEst[,3] > 0 )
bEst$Parameter <-c("b[1]", "b[2]", "b[3]", "b[4]","b[5]","b[6]","b[7]","b[8]", "b[9]" , "b[10]")

# plot using the dotwhisker package 
dwplot(bEst, style = "dotwhisker", 
                  dot_args = list(aes(colour = factor(color))),
                  whisker_args = list(aes(colour = factor(color))), 
                  vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("Michigan LMB") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none") 

#another plot with the distribution 
model1tranformed %>%
  filter(grepl("b",Parameter))%>%
  left_join(bEst) %>% #get parameter names and colors 
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.05, .95)) + 
  theme_bw() 

### MODEL FIT CODE ####
#this has a posterior dist (30000) for coef for every param; e.g. 214 lakes, 10 betas, 4 logqs, 8 mu.alphas, etc. 
coefs <- output_2level$BUGSoutput$sims.matrix[,1:25] 

### 
# Number of desired MCMC samples used for prediction (use a subset of all samples)
nsim <- 3000
# Chain length from analysis
chainLength <- output_2level$BUGSoutput$n.sims
# Select thinned steps in chain for posterior predictions to ensure we take values from length of posterior
ID = seq( 1 , chainLength , floor(chainLength/nsim) )

# Container for predicted values

predictions_sim<- array(NA, c(nsim,length(dat[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predictions_sim)

for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(dat[[1]])){ #loop over covariate data (obs)
    predictions_sim[i,t] <- exp(coefs[ID[i],paste0('alpha[',data[['group']][t],']')] + coefs[ID[i],'b[1]']*data[['x1']][t]  + coefs[ID[i],'b[2]']*data[['x2']][t]  + coefs[ID[i],'b[3]']*data[['x3']][t]  + 
                                  coefs[ID[i],'b[4]']*data[['x4']][t]  + coefs[ID[i],'b[5]']*data[['x5']][t] + coefs[ID[i],'b[6]' ]*data[['x6']][t] + coefs[ID[i],'b[7]']*data[['x7']][t] +
                                  coefs[ID[i],'b[8]']*data[['x8']][t] + coefs[ID[i],'b[9]']*data[['x9']][t] + coefs[ID[i],'b[10]']*data[['x10']][t] + 
                                  coefs[ID[i], paste0('logq[',data[['gear']][t],']')]*data[['IND']][t] + data[['logeffort']][t]) # exp the predictions because we used a log-link in the Poisson
  }
}


#container for random values 
exp_catch <- array(NA, c(nsim,length(dat[[1]])) )

for(i in 1:nrow(predictions_sim)){ # loop over rows (obs)
  for(t in 1:ncol(predictions_sim) ){ #loop over columns(sim)
    exp_catch[i,t] <- rpois(1, predictions_sim[i,t])
  }
}

#From this distribution, you can extract the median and the 2.5% and 97.5% percentiles (95% credible interval),
col.med <- apply(exp_catch, 2, median )
col.means <- apply(exp_catch, 2, mean )
# 95% CIs for fitted values
upperCI.Group <- apply(exp_catch, 2, quantile, probs=c(0.975) )
lowerCI.Group <- apply(exp_catch, 2, quantile, probs=c(0.025) )

#prep data for plotting
med_data=data.frame(col.med, col.means, upperCI.Group, lowerCI.Group) 
med_data$row <- as.numeric(row.names(med_data))

dat$site <- as.numeric(as.factor(dat$new_key))
obs_catch<-dplyr::select(dat, site, gear2, fish_count_new) %>% 
  mutate(row = row_number())
plot_data<-left_join(med_data, obs_catch) 

plot( log(med_data$col.med + 1), log(obs_catch$fish_count_new+1))
plot( med_data$col.means, obs_catch$fish_count_new)

#which can be compared to the observed catch for that combination of lake and gear. 
#One suggestion for visualizing the model fit is to generate a separate graph for each gear, 
#plot the median predicted catches (+- the 95% ci) of different lakes on the x-axis versus the observed catches of lakes on the y-axis.
pred_plot<-ggplot()+
  geom_pointrange(data=plot_data, aes(x=col.med, y=fish_count_new, xmax=upperCI.Group, xmin=lowerCI.Group))+ 
  xlab('predicted catch')+
  ylab('observed catch')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)

pred_plot+facet_wrap(~ gear2, ncol=2, scales = 'free') #allow scales to vary 

#plot log scale 
pred_plot_log<-ggplot()+
  geom_pointrange(data=plot_data, aes(x=log(col.med+1), y=log(fish_count_new + 1), xmax=log(upperCI.Group +1), xmin=log(lowerCI.Group +1) ) )+ 
  xlab('predicted catch log')+
  ylab('observed catch log')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)

pred_plot_log+facet_wrap(~ gear2, ncol=2, scales = 'free') #allow scales to vary 

### Try to plot model fits  #### 
library(bayesplot)
yrep<-exp_catch
y<-dat$fish_count_new

color_scheme_set("brightblue")
ppc_dens_overlay(y, yrep[1:50, ])
ppc_dens_overlay(y, yrep[1:50, ]) + xlim(0, 50)
ppc_hist(y, yrep[1:5, ])
#group 
ppc_dens_overlay_grouped(y, yrep[1:50, ], group = dat$gear2) + xlim(0, 20)
#ggsave("Ag_Effect_on_intercepts.pdf", height = 12, width = 12, units="in")

### marginal effects ### 
BayesPostEst::mcmcMargEff(mod=output, 
            main='b[3]',
            int='b[2]', 
            plot=TRUE)
### EXAMPLE FROM PREDICTIONS WHEN READY TO HINDCAST #### 

for (i in 1:N) { #contemporary data
  y[i] ~ dnorm(beta0 + beta1 * x[i], sigma) #e.g. x is DD 
}
for (j in 1:Nnew) { # historical data 
  ynew[j] ~ dnorm(beta0 + beta1 * xnew[j], sigma) #e.g. xnew is hindcasted DD 
}

#where y, x and xnew are data vectors and ynew is a variable for storing predictions. 
#What you get is a distribution of values that are plausible given your estimated model. 
#Since the model is probabilistic, the prediction is also probabilistic, i.e. we get the whole distribution of possible ynew values.
#For point-values take the average of ynew

#then compare ynew to actual abundance values 
