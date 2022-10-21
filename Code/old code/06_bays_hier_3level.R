# modeling 3 level hierarchical model with varying intercepts and all gears combined
#written by Katelyn King and Henrique Giacomini
#code adapted from Tyler Wagner 
#separating catches from different gears
#add another level in the model (the site, with multiple catch values, as a random factor). 

#load libraries 
library(R2jags)
library(lme4)
library(MCMCpack)
library(dotwhisker) #https://cran.r-project.org/web/packages/dotwhisker/vignettes/dotwhisker-vignette.html
library(dplyr)
library(ggmcmc) # package for analyzing output 


#### data ####
#need to keep all the gears used in a lake even if that gear did not catch sp of interest 

#lmb_dat_for_model has count by gear and effort by gear
dat<-read.csv("Data/lmb_dat_for_model_mar17.csv") 
dd_temp<-read.csv("Data/MI_data/lake_degree_days_year.csv")
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

## standardize and transform predictor values ##
#dat$z_order<-as.numeric(scale(log(dat$lake_order+0.001))) # min is 0 because these are the isolated
#dat$z_ws_urban<-as.numeric(scale(asin(sqrt(dat$ws_urban_prop)))) #good
#dat$z_ws_agriculture<-as.numeric(scale(asin(sqrt(dat$ws_agriculture_prop)))) #ok
#dat$z_ws_elevation<-as.numeric(scale(log(dat$ws_mean_elevation_m))) #good 
#dat$z_houses<-as.numeric(scale(log(dat$houses_km+ 0.001))) #ok remove houses for now 
dat$z_lake_area<-as.numeric(scale(log(dat$lake_area_m2))) #good 
dat$z_max_depth<-as.numeric(scale(log(dat$maxdepth_m))) #good 
dat$z_dd_mean<-as.numeric(scale(log(dat$dd_year))) #good
dat$z_surface_temp_mean<-as.numeric(scale(log(dat$surface_temp_mean))) #good 
dat$z_secchi<-as.numeric(scale(log(dat$secchi_m)))
dat$z_bottom_do<-as.numeric(scale(log(dat$bottom_do_mgl+ 0.001))) ##has 0s so added 0.001 
dat$z_julian<-as.numeric(scale(dat$julian)) ##standardize julian date 
dat$z_ws_forest<-as.numeric(scale(asin(sqrt(dat$ws_forest_prop)))) #good
dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat$ws_wetland_prop)))) #good

dat$logeffort <- log(dat$effort_new)

dat<-dplyr::select(dat, new_key, fish_count_new, logeffort, gear2, FMU_Code, 
                   z_lake_area, z_max_depth, z_secchi, z_julian, 
                   z_bottom_do, z_dd_mean,
                   z_ws_forest, z_ws_wetland, pike_pres, wae_pres) %>%
  na.omit()
#put Inf values to 0 
#is.na(dat)<-sapply(dat, is.infinite)
#dat[is.na(dat)]<-0

#check out response: log of count is mostly normally dist, so Poisson seems appropriate 
hist(dat$fish_count_new)
hist(log(dat$fish_count_new))

#################################################################
########## BUGS CODE ############################################
#################################################################
#hierarchical model with varying intercepts 
sink("model.txt")
cat("
    model {
    # Likelihood: 
    for (i in 1:n){         
    y[i] ~ dpois(lambda[i])  # Distribution Poisson  
    log(lambda[i]) <- log.lambda[i] # log-link 
    log.lambda[i] <- alpha[site[i]] + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i] + b[4] * x4[i] + 
                       b[5] * x5[i] + b[6] * x6[i] + b[7] * x7[i] + b[8] * x8[i] + b[9] * x9[i] + b[10] * x10[i] +
                      logq[gear[i]]*IND[i] + logeffort[i]
         
    } 
    
    # Level-2 of the model: site 
    for(s in 1:nsites){
    alpha[s] ~ dnorm(mu.alpha[group[s]],tau.alpha)
    }
    
     # Level-3 of the model: group (fisheries management unit regions) 
    for(j in 1:J){
    mu.alpha[j] ~ dnorm(mu.alpha2,tau.alpha2)
    }
    
    # Priors
    mu.alpha2 ~ dnorm(0, 0.0001)
    sigma.alpha ~ dunif(0,100)
    sigma.alpha2 ~ dunif(0,100)
    
    #derived quantities 
    tau.alpha <- pow(sigma.alpha,-2)
    tau.alpha2 <- pow(sigma.alpha2,-2)
    
    # priors for predictors 
    for(k in 1:10){
    b[k] ~ dnorm(0,0.0000001)
    }
    
     # priors for gear-specific log-catchabilities 
    for (k in 1:ngears) { 
    logq[k] ~ dnorm(0,  0.0000001)   #non-informative
    } 

                      
    } # end model
    ",fill = TRUE)
sink()

# if you want to predict within jags add this inside the model script 
#for (i in 1:n){         
 # ysim[i] ~ dpois(lambdasim[i])  # Distribution Poisson  
#  log(lambdasim[i]) <- log.lambdasim[i] # log-link 
 # log.lambdasim[i] <- alphasim[site[i]] + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i] + b[4] * x4[i] + 
  #  b[5] * x5[i] + b[6] * x6[i] + b[7] * x7[i] + b[8] * x8[i] + b[9] * x9[i] + b[10] * x10[i] +
   # logq[gear[i]]*IND[i] + logeffort[i]
#} 

# Level-2 of the model: site 
#for(s in 1:nsites){
 # alphasim[s] ~ dnorm(mu.alphasim[group[s]],tau.alpha)
#}

# Level-3 of the model: group (fisheries management unit regions) 
#for(j in 1:J){
 # mu.alphasim[j] ~ dnorm(mu.alpha2,tau.alpha2)
#}

### set parameters ####
# Initial values
inits <- function (){
  list (mu.alpha = rnorm(8))
}


# Parameters monitored
parameters <- c("alpha", "mu.alpha","tau.alpha", "mu.alpha2","tau.alpha2","b", 'logq', 'ysim')


# MCMC settings
ni <- 100000 # number of iterations 
nt <- 10 #number to thin
nb <- 30000 #number to burn  
nc <- 3 #number of chains


#### set up data #### 

#set the number of sites and create site index 
nsites<-length(unique(dat$new_key))
site <- as.numeric(as.factor(dat$new_key))

# Set the number of FMUs 
J <- length(unique(dat$FMU_Code))

# Create group index, must go from 1,..J
dat$FMU_Code <- droplevels(as.factor(dat$FMU_Code))
dat$group <- as.numeric(as.factor(dat$FMU_Code))
group <- doBy::summaryBy(group ~ new_key, data=dat, FUN=mean) #need region to be indexed by site

#set number of different sampling gears used and set a gear index 
ngears<-length(unique(dat$gear2))
gear <- as.numeric(as.factor(dat$gear2))

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
dat$IND<-ifelse(dat$gear2 == "FT_NET", 0, 1) #FT_NET is the reference gear

# Load data
data <- list(y = dat$fish_count_new, group = group$group.mean, gear=gear, n = dim(dat)[1], J = J, ngears = ngears,
             x1=dat$z_secchi, x2=dat$z_lake_area, x3=dat$z_dd_mean, x4=dat$z_max_depth, x5=dat$z_bottom_do,
             x6=dat$z_ws_forest,x7=dat$z_ws_wetland, x8=dat$z_julian, x9=dat$wae_pres, x10=dat$pike_pres,
             logeffort=dat$logeffort, IND=dat$IND, nsites = nsites, site =site
) 

### run the model using JAGS and R ####
# Set timer 
start.time = Sys.time()         

# Call JAGS from R and run model 
output<- jags(data, inits, parameters, "model.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb)


# # Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

#save output 
#saveRDS(output, "Data/output/model_output_lmb_3level.rds") 

#jag.sum<-output$BUGSoutput$summary
#write.table(x=jag.sum,file="out.txt",sep="\t")

# Summarize posteriors
print(output, dig = 3)
#in this case, logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
#check convergence with Brooks-Gelman-Rubin statistic
which(output$BUGSoutput$summary[, c("Rhat")] > 1.1)
which(output$BUGSoutput$summary[, c("n.eff")] < 0.1*30000) #times number of retained draws 
which(output$BUGSoutput$summary[, c("n.eff")] < 500)
jagsfit.mcmc <- as.mcmc(output)
model1tranformed <- ggs(jagsfit.mcmc) 
ggs_traceplot(model1tranformed, family = "b")
#check stationarity Heidelberg-Welch convergence 
heidel.diag(jagsfit.mcmc)

#### convergence checks #### 
output1<-readRDS("Data/output/model_output_lmb_3level.rds")
#check convergence with Brooks-Gelman-Rubin statistic
#https://stats.stackexchange.com/questions/418142/low-effective-sample-size-but-good-r-hat-is-this-a-problem
which(output1$BUGSoutput$summary[, c("Rhat")] > 1.1)
#effective sample size check
which(output1$BUGSoutput$summary[, c("n.eff")] < 0.1*30000) #times number of retained draws 

#look at trace plots 
# use as.mcmmc to convert rjags object into mcmc.list - plots in coda
jagsfit.mcmc <- as.mcmc(output1)
require(lattice)
densityplot(jagsfit.mcmc)

#the ggs function transforms the mcmc output into a longformat tibble, that we can use to make different types of plots.
model1tranformed <- ggs(jagsfit.mcmc) 
summary(model1tranformed$Parameter)

#traceplots (caterpillars)
ggs_traceplot(model1tranformed, family = "alpha") #note there are 214 alphas
ggs_traceplot(model1tranformed, family = "b")
ggs_traceplot(model1tranformed, family = "logq")

#### try to update with more chains #### 
# if the model does not converge, update it!
recompile(output1)
jagsfit.upd <- update(output1, n.iter=50000)
print(jagsfit.upd, dig=3)
print(jagsfit.upd, intervals=c(0.025, 0.5, 0.975))

#check again
which(output1$BUGSoutput$summary[, c("Rhat")] > 1.1)
which(output1$BUGSoutput$summary[, c("n.eff")] < 0.1*30000) #times number of retained draws 
jagsfit.mcmc <- as.mcmc(output1)
model1tranformed <- ggs(jagsfit.mcmc) 
ggs_traceplot(model1tranformed, family = "b")

#### effect plots #### 
output1<-readRDS("Data/output/model_output_lmb_3level.rds")
library(dotwhisker)

#betas
bEst <- matrix(NA, nrow=10,ncol=3)
for(i in 1:10){ #parameters
  bEst[i,1] <- mean(output$BUGSoutput$sims.list$b[,i])
  bEst[i,2:3] <- quantile(output$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
bEst<-as.data.frame(bEst)
bEst$variable<-c("Secchi", "lake area", "dd_year", "max_depth", "bottom_do", 
                 "ws forest", "ws wetland",  "julian day", "pres_wae", "pres_pike" )

colnames(bEst)<- c("estimate", "conf.low", "conf.high", "term")
bEst$Parameter <-c("b[1]", "b[2]", "b[3]", "b[4]","b[5]","b[6]","b[7]","b[8]", "b[9]" , "b[10]")


#add colors for sig different than 0
bEst$color <- as.numeric(bEst[,2] * bEst[,3] > 0 )

# plot using the dotwhisker package 
lake_plot<-dwplot(bEst, style = "dotwhisker", 
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

#PLOT WITH DISTRIBUTIONS 
jagsfit.mcmc <- as.mcmc(output1)
model1tranformed <- ggs(jagsfit.mcmc) 

model1tranformed %>%
  filter(grepl("b",Parameter))%>%
  left_join(bEst) %>% #get parameter names and colors 
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.05, .95)) + 
  theme_bw() 


#### MODEL FIT #### 
# Simulating data from the posterior predictive distribution using the observed predictors is useful for checking the fit of the model.
#obtain our samples from the posterior
#this has a posterior dist (30000) for coef for every param; e.g. 214 lakes, 10 betas, 4 logqs, 8 mu.alphas, etc. 
coefs <- output1$BUGSoutput$sims.matrix[,1:240] 

### 
# Number of desired MCMC samples used for prediction (use a subset of all samples)
nsim <- 3000
# Chain length from analysis
chainLength <- output1$BUGSoutput$n.sims
# Select thinned steps in chain for posterior predictions to ensure we take values from length of posterior
ID = seq( 1 , chainLength , floor(chainLength/nsim) )

# Container for predicted values
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
#predictions_sim<- array(NA, c(nsim,length(dat[[1]]),nsites)) # 3 dimensions - simulations, length of data, groups 
#dim(predictions_sim)

#for(s in 1:nsites){ # loop over lakes 
 # for(i in 1:nsim ){  #loop over sims 
  #  for(t in 1:length(dat[[1]])){ #loop over covariate data (obs)
   #   predictions_sim[i,t,s] <- exp(coefs[ID[i],paste0('alpha[',s,']')] + coefs[ID[i],'b[1]']*data[['x1']][t]  + coefs[ID[i],'b[2]']*data[['x2']][t]  + coefs[ID[i],'b[3]']*data[['x3']][t]  + 
    #                               coefs[ID[i],'b[4]']*data[['x4']][t]  + coefs[ID[i],'b[5]']*data[['x5']][t] + coefs[ID[i],'b[6]' ]*data[['x6']][t] + coefs[ID[i],'b[7]']*data[['x7']][t] +
     #                                coefs[ID[i],'b[8]']*data[['x8']][t] + coefs[ID[i],'b[9]']*data[['x9']][t] + coefs[ID[i],'b[10]']*data[['x10']][t] + 
      #                             coefs[ID[i], paste0('logq[',data[['gear']][t],']')]*data[['IND']][t] + data[['logeffort']][t]) # exp the predictions because we used a log-link in the Poisson
#    }
#  }
#}

predictions_sim<- array(NA, c(nsim,length(dat[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predictions_sim)

for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(dat[[1]])){ #loop over covariate data (obs)
   predictions_sim[i,t] <- exp(coefs[ID[i],paste0('alpha[',data[['site']][t],']')] + coefs[ID[i],'b[1]']*data[['x1']][t]  + coefs[ID[i],'b[2]']*data[['x2']][t]  + coefs[ID[i],'b[3]']*data[['x3']][t]  + 
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

#ggsave("Ag_Effect_on_intercepts.pdf", height = 12, width = 12, units="in")


### Try to plot model fits  #### 
library(bayesplot)
#yrep<-as.matrix(output$BUGSoutput$sims.list$ysim ) #from the predictions within JAGS
yrep<-predictions_sim
y<-dat$fish_count_new

color_scheme_set("brightblue")
ppc_dens_overlay(y, yrep[1:50, ])
ppc_dens_overlay(y, yrep[1:50, ]) + xlim(0, 50)
ppc_hist(y, yrep[1:5, ])
#group 
ppc_dens_overlay_grouped(y, yrep[1:50, ], group = dat$gear2) + xlim(0, 20)
#look at how well it predicts the 0 obs 
prop_zero <- function(x) mean(x == 0)
prop_zero(y) 
ppc_stat(y, yrep, stat = "prop_zero", binwidth = 0.005)
