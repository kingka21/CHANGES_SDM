## Negative Binomial dist including varying gear  #### 
# betas vary by gear, but fixed 
#load libraries 
library(R2jags)
library(lme4)
library(MCMCpack)
library(dotwhisker) #https://cran.r-project.org/web/packages/dotwhisker/vignettes/dotwhisker-vignette.html
library(dplyr)
library(ggmcmc) # package for analyzing output 
library(tidybayes)

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
                   z_bottom_do, z_surface_temp_mean,
                   z_ws_forest, z_ws_wetland, pike_pres, wae_pres) %>%
  na.omit()

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
    y[i] ~ dnegbin(p[i],r)
    p[i] <- r/(r + lambda[i])
    log(lambda[i]) <- log.lambda[i] # log-link 
    log.lambda[i] <- alpha[site[i]] + b1[gear[i]] * x1[i] + b2[gear[i]] * x2[i] + b3[gear[i]] * x3[i] + b4[gear[i]] * x4[i] + 
                       b5[gear[i]]*x5[i] + b6[gear[i]]*x6[i] + b7[gear[i]]*x7[i] + b8[gear[i]]*x8[i] + b9[gear[i]]*x9[i] + b10[gear[i]]*x10[i] +
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
    
    # varying betas by gear 
    for(j in 1:ngears){
    b1[j] ~ dnorm(0,  0.0000001)
    b2[j] ~ dnorm(0,  0.0000001)
    b3[j] ~ dnorm(0,  0.0000001)
    b4[j] ~ dnorm(0,  0.0000001)
    b5[j] ~ dnorm(0,  0.0000001)
    b6[j] ~ dnorm(0,  0.0000001)
    b7[j] ~ dnorm(0,  0.0000001)
    b8[j] ~ dnorm(0,  0.0000001)
    b9[j] ~ dnorm(0,  0.0000001)
    b10[j] ~ dnorm(0,  0.0000001)
    }
    
    # Priors
    mu.alpha2 ~ dnorm(0, 0.0001)
    sigma.alpha ~ dunif(0,100)
    sigma.alpha2 ~ dunif(0,100)
    r ~ dunif(0,100)
    
    #derived quantities 
    tau.alpha <- pow(sigma.alpha,-2)
    tau.alpha2 <- pow(sigma.alpha2,-2)
    
    
    # priors for gear-specific log-catchabilities 
    for (k in 1:ngears) { 
    logq[k] ~ dnorm(0,  0.0000001)   #non-informative
    } 
    
                      
    } # end model
    ",fill = TRUE)
sink()


### set parameters ####

# Initial values #8 is the number of regions, 
inits <- function (){
  list (mu.alpha = rnorm(8), b1=rnorm(4), b2=rnorm(4), b3=rnorm(4), b4=rnorm(4), 
        b5=rnorm(4), b6=rnorm(4), b7=rnorm(4), b8=rnorm(4), b9=rnorm(4), b10=rnorm(4)
  )
}


# Parameters monitored
parameters <- c("alpha", "mu.alpha","tau.alpha", "mu.alpha2","tau.alpha2",'b1', "b2", "b3", "b4", "b5", 
                "b6", "b7", "b8", "b9", "b10",
                'logq', 'r')


# MCMC settings 
ni <- 130000 # number of iterations 130000
nt <- 20 #number to thin 20
nb <- 30000 #number to burn   30000
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
             x1=dat$z_secchi, x2=dat$z_lake_area, x3=dat$z_surface_temp_mean, x4=dat$z_max_depth, x5=dat$z_bottom_do,
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
#saveRDS(output, "Data/output/output_lmb_nbinom_surface_varslopes4.rds") 
output<-readRDS("Data/output/output_lmb_nbinom_surface_varslopes4.rds")
# Summarize posteriors
print(output, dig = 3)
#in this case, logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
#check convergence with Brooks-Gelman-Rubin statistic
which(output$BUGSoutput$summary[, c("Rhat")] > 1.1)
which(output$BUGSoutput$summary[, c("n.eff")] < 0.1*15000) #times number of retained draws 
which(output$BUGSoutput$summary[, c("n.eff")] < 1000)
jagsfit.mcmc <- as.mcmc(output)
model1tranformed <- ggs(jagsfit.mcmc) 

ggs_traceplot(model1tranformed, family = "b9")
ggs_traceplot(model1tranformed, family = "r")
ggs_traceplot(model1tranformed, family = "logq")

#### effect plots #### 
library(dotwhisker)

#betas
#1 is FT, 2 is gill, 3 is seine, 4 is shock 

EstsSecchi <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  EstsSecchi[i,1] <- mean(output$BUGSoutput$sims.list$b1[,i])
  EstsSecchi[i,2:3] <- quantile(output$BUGSoutput$sims.list$b1[,i],c(0.05,0.95))
}
EstsSecchi<-as.data.frame(EstsSecchi)
colnames(EstsSecchi)<- c("estimate", "conf.low", "conf.high")
EstsSecchi$gear<-c("Fyke", "Gill", "Seine", "Shock")
EstsSecchi$term<-'Secchi'

Estsarea <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  Estsarea[i,1] <- mean(output$BUGSoutput$sims.list$b2[,i])
  Estsarea[i,2:3] <- quantile(output$BUGSoutput$sims.list$b2[,i],c(0.05,0.95))
}
Estsarea<-as.data.frame(Estsarea)
colnames(Estsarea)<- c("estimate", "conf.low", "conf.high")
Estsarea$gear<-c("Fyke", "Gill", "Seine", "Shock")
Estsarea$term<-'lake_area'

Eststemp <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  Eststemp[i,1] <- mean(output$BUGSoutput$sims.list$b3[,i])
  Eststemp[i,2:3] <- quantile(output$BUGSoutput$sims.list$b3[,i],c(0.05,0.95))
}
Eststemp<-as.data.frame(Eststemp)
colnames(Eststemp)<- c("estimate", "conf.low", "conf.high")
Eststemp$gear<-c("Fyke", "Gill", "Seine", "Shock")
Eststemp$term<-"surface_temp"

Estsdepth <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  Estsdepth[i,1] <- mean(output$BUGSoutput$sims.list$b4[,i])
  Estsdepth[i,2:3] <- quantile(output$BUGSoutput$sims.list$b4[,i],c(0.05,0.95))
}
Estsdepth<-as.data.frame(Estsdepth)
colnames(Estsdepth)<- c("estimate", "conf.low", "conf.high")
Estsdepth$gear<-c("Fyke", "Gill", "Seine", "Shock")
Estsdepth$term<-"max_depth"

EstsDO <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  EstsDO[i,1] <- mean(output$BUGSoutput$sims.list$b5[,i])
  EstsDO[i,2:3] <- quantile(output$BUGSoutput$sims.list$b5[,i],c(0.05,0.95))
}
EstsDO<-as.data.frame(EstsDO)
colnames(EstsDO)<- c("estimate", "conf.low", "conf.high")
EstsDO$gear<-c("Fyke", "Gill", "Seine", "Shock")
EstsDO$term<-"bottom_DO"

EstsFor <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  EstsFor[i,1] <- mean(output$BUGSoutput$sims.list$b6[,i])
  EstsFor[i,2:3] <- quantile(output$BUGSoutput$sims.list$b6[,i],c(0.05,0.95))
}
EstsFor<-as.data.frame(EstsFor)
colnames(EstsFor)<- c("estimate", "conf.low", "conf.high")
EstsFor$gear<-c("Fyke", "Gill", "Seine", "Shock")
EstsFor$term<-'ws_forest'

Estswet <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  Estswet[i,1] <- mean(output$BUGSoutput$sims.list$b7[,i])
  Estswet[i,2:3] <- quantile(output$BUGSoutput$sims.list$b7[,i],c(0.05,0.95))
}
Estswet<-as.data.frame(Estswet)
colnames(Estswet)<- c("estimate", "conf.low", "conf.high")
Estswet$gear<-c("Fyke", "Gill", "Seine", "Shock")
Estswet$term<-'ws_wetland'

Estsday <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  Estsday[i,1] <- mean(output$BUGSoutput$sims.list$b8[,i])
  Estsday[i,2:3] <- quantile(output$BUGSoutput$sims.list$b8[,i],c(0.05,0.95))
}
Estsday<-as.data.frame(Estsday)
colnames(Estsday)<- c("estimate", "conf.low", "conf.high")
Estsday$gear<-c("Fyke", "Gill", "Seine", "Shock")
Estsday$term<-"julian day"

Estswae <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  Estswae[i,1] <- mean(output$BUGSoutput$sims.list$b9[,i])
  Estswae[i,2:3] <- quantile(output$BUGSoutput$sims.list$b9[,i],c(0.05,0.95))
}
Estswae<-as.data.frame(Estswae)
colnames(Estswae)<- c("estimate", "conf.low", "conf.high")
Estswae$gear<-c("Fyke", "Gill", "Seine", "Shock")
Estswae$term<-'pres_wae'

Estspike <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  Estspike[i,1] <- mean(output$BUGSoutput$sims.list$b10[,i])
  Estspike[i,2:3] <- quantile(output$BUGSoutput$sims.list$b10[,i],c(0.05,0.95))
}
Estspike<-as.data.frame(Estspike)
colnames(Estspike)<- c("estimate", "conf.low", "conf.high")
Estspike$gear<-c("Fyke", "Gill", "Seine", "Shock")
Estspike$term<-"pres_pike"

all_est<- gtools::smartbind(EstsSecchi, Estsarea, Eststemp, Estsdepth, EstsDO, EstsFor, Estswet, Estsday, Estswae, Estspike)
#add colors for sig different than 0
all_est$color <- as.numeric(all_est[,2] * all_est[,3] > 0 )
all_est$Parameter <-c("b1[1]", "b1[2]", "b1[3]", "b1[4]",
                      "b2[1]", "b2[2]", "b2[3]", "b2[4]",
                      "b3[1]", "b3[2]", "b3[3]", "b3[4]",
                      "b4[1]", "b4[2]", "b4[3]", "b4[4]",
                      "b5[1]", "b5[2]", "b5[3]", "b5[4]",
                      "b6[1]", "b6[2]", "b6[3]", "b6[4]",
                      "b7[1]", "b7[2]", "b7[3]", "b7[4]",
                      "b8[1]", "b8[2]", "b8[3]", "b8[4]",
                      "b9[1]", "b9[2]", "b9[3]", "b9[4]",
                      "b10[1]", "b10[2]", "b10[3]", "b10[4]")

#split into the gears for graphing 
fyke<-filter(all_est, gear=="Fyke")
gill<-filter(all_est, gear=="Gill")
seine<-filter(all_est, gear=="Seine")
shock<-filter(all_est, gear=="Shock")



# plot using the dotwhisker package 
fyke_plot<-dwplot(fyke, style = "dotwhisker", 
                  dot_args = list(aes(colour = factor(color))),
                  whisker_args = list(aes(colour = factor(color))), 
                  vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("fyke") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
fyke_plot

gill_plot<-dwplot(gill, style = "dotwhisker", 
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

seine_plot<-dwplot(seine, style = "dotwhisker", 
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

shock_plot<-dwplot(shock, style = "dotwhisker", 
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

#####PLOT WITH DISTRIBUTIONS ####
#1 is FT, 2 is gill, 3 is seine, 4 is shock 
model1tranformed %>%
  filter(grepl("b.*[1]",Parameter))%>% #use .* as a wild card to capture all 10 params
  left_join(fyke) %>% #get parameter names and colors 
  na.omit(term) %>%
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.05, .95)) + 
  theme_bw() 

model1tranformed %>%
  filter(grepl("b.*[2]",Parameter))%>% #use .* as a wild card to capture all 10 params
  left_join(gill) %>% #get parameter names and colors 
  na.omit(term) %>%
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.05, .95)) + 
  theme_bw() 

model1tranformed %>%
  filter(grepl("b.*[3]",Parameter))%>% #use .* as a wild card to capture all 10 params
  left_join(seine) %>% #get parameter names and colors 
  na.omit(term) %>%
  ggplot(aes(x = value, y = term, fill= factor(color)),  xmin=-5, xmax=5) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.10, .90)) + 
  theme_bw() 

model1tranformed %>%
  filter(grepl("b.*[4]",Parameter))%>% #use .* as a wild card to capture all 10 params
  left_join(shock) %>% #get parameter names and colors 
  na.omit(term) %>%
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.10, .90)) + 
  theme_bw() 

#### MODEL FIT #### 
# Simulating data from the posterior predictive distribution using the observed predictors is useful for checking the fit of the model.
#obtain our samples from the posterior
#this has a posterior dist (30000) for coef for every param; e.g. 214 lakes, 10x4 betas, 4 logqs, 8 mu.alphas, etc. 
coefs <- output$BUGSoutput$sims.matrix[,1:271] 

### 
# Number of desired MCMC samples used for prediction (use a subset of all samples)
nsim <- 3000
# Chain length from analysis
chainLength <- output$BUGSoutput$n.sims
# Select thinned steps in chain for posterior predictions to ensure we take values from length of posterior
ID = seq( 1 , chainLength , floor(chainLength/nsim) )

# Container for predicted values
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
predictions_sim<- array(NA, c(nsim,length(dat[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predictions_sim)

for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(dat[[1]])){ #loop over covariate data (obs)
    predictions_sim[i,t] <- exp(coefs[ID[i],paste0('alpha[',data[['site']][t],']')] + coefs[ID[i], paste0('b1[',data[['gear']][t],']')]*data[['x1']][t]  + coefs[ID[i],paste0('b2[',data[['gear']][t],']')]*data[['x2']][t]  + coefs[ID[i],paste0('b3[',data[['gear']][t],']')]*data[['x3']][t]  + 
                                  coefs[ID[i],paste0('b4[',data[['gear']][t],']')]*data[['x4']][t]  + coefs[ID[i],paste0('b5[',data[['gear']][t],']')]*data[['x5']][t] + coefs[ID[i],paste0('b6[',data[['gear']][t],']')]*data[['x6']][t] + coefs[ID[i],paste0('b7[',data[['gear']][t],']')]*data[['x7']][t] +
                                  coefs[ID[i],paste0('b8[',data[['gear']][t],']')]*data[['x8']][t] + coefs[ID[i],paste0('b9[',data[['gear']][t],']')]*data[['x9']][t] + coefs[ID[i],paste0('b10[',data[['gear']][t],']')]*data[['x10']][t] + 
                                  coefs[ID[i], paste0('logq[',data[['gear']][t],']')]*data[['IND']][t] + data[['logeffort']][t]) # exp the predictions because we used a log-link in the negbi
    
  }
}


#container for p values 
#parametrization is by the dispersion parameter, where prob = size/(size+mu).
#p is the success parameter and r is the dispersion parameter.
p <- array(NA, c(nsim,length(dat[[1]])) )
exp_catch <- array(NA, c(nsim,length(dat[[1]])) )

for(i in 1:nrow(predictions_sim)){ # loop over rows (sims)
  for(t in 1:ncol(predictions_sim) ){ #loop over columns(obs)
    p[i,t] <- coefs[ID[i],'r']/(coefs[ID[i],'r'] + predictions_sim[i,t])
    exp_catch[i,t] <- rnbinom(1, prob=p[i,t], size= coefs[ID[i],'r'])
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

plot( med_data$col.med , obs_catch$fish_count_new)
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

###  plot model fits  #### 
library(bayesplot)
#yrep<-as.matrix(output$BUGSoutput$sims.list$ysim ) #from the predictions within JAGS
yrep<-exp_catch
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

#*residuals for contemporary model  
#Bayes-R2 using residuals
e <- -1 * sweep(yrep, 2, y)
var_ypred <- apply(yrep, 1, var)
var_e <- apply(e, 1, var)
bayesR2res<- var_ypred / (var_ypred + var_e)
round(median(bayesR2res), 2) 

#residual plot
#https://cran.r-project.org/web/packages/tidybayes/vignettes/tidybayes-residuals.html
plot_data$res.med <- apply(e, 2, median )
res_plot<-plot_data %>%
  ggplot(aes(x = col.med, y = res.med)) +
  stat_pointinterval()
res_plot
res_plot+facet_wrap(~ gear2, ncol=2, scales = 'free')

plot_data  %>%
  ggplot(aes(sample = res.med)) +
  geom_qq() +
  geom_qq_line()

#### PREDICTION ##### 
hist_dat<-read.csv("/Users/katelynking/Desktop/historical_lmb.csv")

#fix julian day 
hist_dat$julian<-ifelse(hist_dat$new_key == '65-63', 181, hist_dat$julian)


#set up data 
hist_dat$logeffort<-log(hist_dat$effort_sum)
## standardize and transform predictor values ##
hist_dat$z_lake_area<-as.numeric(scale(log(hist_dat$lake_area_m2))) #good 
hist_dat$z_max_depth<-as.numeric(scale(log(hist_dat$max_depth_m))) #good 
hist_dat$z_dd_mean<-as.numeric(scale(log(hist_dat$dd_year))) #good
hist_dat$z_surface_temp_mean<-as.numeric(scale(log(hist_dat$surface_temp_year))) #good 
hist_dat$z_secchi<-as.numeric(scale(log(hist_dat$secchi_m)))
hist_dat$z_bottom_do<-as.numeric(scale(log(hist_dat$bottom_do + 0.001))) ##has 0s so added 0.001 
hist_dat$z_julian<-as.numeric(scale(hist_dat$julian)) ##standardize julian date 
hist_dat$z_ws_forest<-as.numeric(scale(asin(sqrt(hist_dat$ws_forest_prop)))) #good
hist_dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(hist_dat$ws_wetland_prop)))) #good

# need group number to match initial data 
#MCM =1 Central Lake Michigan
#MER = 2 Lake Erie
# MES = 3 Eastern Lake Superior
#MNH = 4 Northern Lake Huron
# MNM = 5 Northern Lake Michigan
#MSH = 6 Southern Lake Huron
#MSM = 7 Southern Lake Michigan
#MWS= 8 Western Lake Superior
hist_dat<-hist_dat %>% 
  mutate(
    group = case_when( MGMT_UNIT == "Central Lake Michigan" ~ 1,
                       MGMT_UNIT == "Lake Erie" ~ 2,
                       MGMT_UNIT == "Eastern Lake Superior" ~ 3, 
                       MGMT_UNIT == "Northern Lake Huron" ~ 4,
                       MGMT_UNIT == "Northern Lake Michigan" ~ 5,
                       MGMT_UNIT== "Southern Lake Huron" ~ 6,
                       MGMT_UNIT== "Southern Lake Michigan" ~ 7,
                       MGMT_UNIT== "Western Lake Superior" ~ 8)
  )

#need gear number to match original data (no shock) 
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
hist_dat<-hist_dat %>% 
  mutate(
    gear_index = case_when( gear == "FT_net" ~ 1,
                            gear == "gill" ~ 2,
                            gear == "seine" ~ 3)
  )

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
hist_dat$IND<-ifelse(hist_dat$gear == "FT_net", 0, 1) #FT_NET is the reference gear

#wae and pike pres/abs 
hist_dat$pres_wae<- ifelse(hist_dat$walleye_sum >= 1, 1, 0)
hist_dat$pres_pike<- ifelse(hist_dat$pike_sum >= 1, 1, 0)

# Container for predicted values
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
predict_hist<- array(NA, c(nsim,length(hist_dat[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predict_hist)

#use regional mu.alphas instead of site specific, since dif sites (this takes like 5 mins)
for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(hist_dat[[1]])){ #loop over covariate data (obs)
    predict_hist[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',hist_dat[['group']][t],']')] + coefs[ID[i],paste0('b1[',hist_dat[['gear_index']][t],']')]*hist_dat[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',hist_dat[['gear_index']][t],']')]*hist_dat[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',hist_dat[['gear_index']][t],']')]*hist_dat[['z_surface_temp_mean']][t]  + 
                               coefs[ID[i],paste0('b4[',hist_dat[['gear_index']][t],']')]*hist_dat[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',hist_dat[['gear_index']][t],']')]*hist_dat[['z_bottom_do']][t] + coefs[ID[i],paste0('b6[',hist_dat[['gear_index']][t],']')]*hist_dat[['z_ws_forest']][t] + coefs[ID[i],paste0('b7[',hist_dat[['gear_index']][t],']')]*hist_dat[['z_ws_wetland']][t] +
                               coefs[ID[i],paste0('b8[',hist_dat[['gear_index']][t],']')]*hist_dat[['z_julian']][t] + coefs[ID[i],paste0('b9[',hist_dat[['gear_index']][t],']')]*hist_dat[['pres_wae']][t] + coefs[ID[i],paste0('b10[',hist_dat[['gear_index']][t],']')]*hist_dat[['pres_pike']][t] + 
                               coefs[ID[i], paste0('logq[',hist_dat[['gear_index']][t],']')]*hist_dat[['IND']][t] + hist_dat[['logeffort']][t]) # exp the predictions because we used a log-link in the Poisson
    
  }
}


#container for p values 
#parametrization is by the dispersion parameter, where prob = size/(size+mu).
#p is the success parameter and r is the dispersion parameter.
p_hist <- array(NA, c(nsim,length(hist_dat[[1]])) )
exp_catch_hist <- array(NA, c(nsim,length(hist_dat[[1]])) )

for(i in 1:nrow(predict_hist)){ # loop over rows (sims)
  for(t in 1:ncol(predict_hist) ){ #loop over columns(obs)
    p_hist[i,t] <- coefs[ID[i],'r']/(coefs[ID[i],'r'] + predict_hist[i,t])
    exp_catch_hist[i,t] <- rnbinom(1, prob=p_hist[i,t], size= coefs[ID[i],'r'])
  }
}

#From this distribution, you can extract the median and the 2.5% and 97.5% percentiles (95% credible interval),
hist.med <- apply(exp_catch_hist, 2, median )
hist.means <- apply(exp_catch_hist, 2, mean )
# 95% CIs for fitted values
upperPI <- apply(exp_catch_hist, 2, quantile, probs=c(0.975) )
lowerPI <- apply(exp_catch_hist, 2, quantile, probs=c(0.025) )

#prep data for plotting
med_data=data.frame(hist.med, hist.means, upperPI, lowerPI) 
med_data$row <- as.numeric(row.names(med_data))

hist_dat$site <- as.numeric(as.factor(hist_dat$new_key))
hist_catch<-dplyr::select(hist_dat, site, gear, largemouthbass_sum) %>% 
  mutate(row = row_number())
plot_data<-left_join(med_data, hist_catch) 

plot( med_data$hist.med , hist_catch$largemouthbass_sum)
plot( med_data$hist.means, hist_catch$largemouthbass_sum)

#which can be compared to the observed catch for that combination of lake and gear. 
#One suggestion for visualizing the model fit is to generate a separate graph for each gear, 
#plot the median predicted catches (+- the 95% ci) of different lakes on the x-axis versus the observed catches of lakes on the y-axis.
pred_plot<-ggplot()+
  geom_pointrange(data=plot_data, aes(x=hist.med, y=largemouthbass_sum, xmax=upperPI, xmin=lowerPI))+ 
  xlab('predicted catch')+
  ylab('historical catch')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)

pred_plot+facet_wrap(~ gear, ncol=2, scales = 'free') #allow scales to vary 

#plot log scale 
pred_plot_log<-ggplot()+
  geom_pointrange(data=plot_data, aes(x=log(hist.med+1), y=log(largemouthbass_sum + 1), xmax=log(upperPI +1), xmin=log(lowerPI +1) ) )+ 
  xlab('predicted catch log')+
  ylab('historical catch log')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)

pred_plot_log+facet_wrap(~ gear, ncol=2, scales = 'free') #allow scales to vary 

#means log scale 
hist_plot_log<-ggplot()+
  geom_pointrange(data=plot_data, aes(x=log(hist.means+1), y=log(largemouthbass_sum + 1), xmax=log(upperPI +1), xmin=log(lowerPI +1) ) )+ 
  xlab('predicted catch log')+
  ylab('historical catch log')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)

hist_plot_log+facet_wrap(~ gear, ncol=2, scales = 'free') #allow scales to vary 

#### RESIDUALS #### 
#Bayes-R2 using residuals
#https://avehtari.github.io/bayes_R2/bayes_R2.html#2_Functions_for_Bayesian_R-squared_for_stan_glm_models
ypred<-exp_catch_hist
y<-plot_data$largemouthbass_sum
e <- -1 * sweep(ypred, 2, y) #residual distribution 
var_ypred <- apply(ypred, 1, var) # variance of modelled predictive means
var_e <- apply(e, 1, var) #variance of modeled residuals 
bayesR2res<- var_ypred / (var_ypred + var_e) #R^2
round(median(bayesR2res), 2) #0.34

#plot the R2 posterior and median
mcmc_hist(data.frame(bayesR2res), binwidth=0.02)+xlim(c(0,1))+
  xlab(('R^2'))+
  geom_vline(xintercept=median(bayesR2res))+
  ggtitle('Bayesian R squared posterior and median')

#plot fits 
color_scheme_set("brightblue")
ppc_dens_overlay(y, ypred[1:50, ])
ppc_dens_overlay(y, ypred[1:50, ]) + xlim(0, 50)
ppc_hist(y, ypred[1:5, ])
#group 
ppc_dens_overlay_grouped(y, ypred[1:50, ], group = hist_dat$gear) + xlim(0, 20)

#residual plot
#https://cran.r-project.org/web/packages/tidybayes/vignettes/tidybayes-residuals.html
plot_data$res.med <- apply(e, 2, median )
res_plot<-plot_data %>%
  ggplot(aes(x = hist.med, y = res.med)) +
  stat_pointinterval()

res_plot+facet_wrap(~ gear, ncol=2, scales = 'free')

plot_data  %>%
  ggplot(aes(sample = res.med)) +
  geom_qq() +
  geom_qq_line()
