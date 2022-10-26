## Negative Binomial - remove 4 in-situ variables #### 

#load libraries 
library(R2jags)
library(lme4)
library(MCMCpack)
library(dotwhisker) #https://cran.r-project.org/web/packages/dotwhisker/vignettes/dotwhisker-vignette.html
library(dplyr)
library(ggmcmc) # package for analyzing output 
library(tidybayes)
library(splitTools)
#### data ####
#lmb_dat_for_model has count by gear and effort by gear
lmb_data<-read.csv("Data/lmb_dat_for_model_aug30.csv") 

#surface temp by year 
surface_temp<-read.csv("Data/MI_data/lake_surface_temp.csv")%>%
  group_by(IHDLKID) %>%
  slice_max(HECTARES) %>% 
  rename(nhdid = IHDLKID)%>% 
  pivot_longer(
    cols= starts_with("TAVE_"),
    names_to = "year", 
    names_prefix = "TAVE_",
    values_to = "surf_temp_year") %>% 
  mutate(year = as.integer(as.character(year)))

dd_temp<-read.csv("Data/MI_data/lake_degree_days_year.csv")%>%
  group_by(IHDLKID) %>%
  slice_max(HECTARES) %>% 
  rename(nhdid = IHDLKID)%>% 
  pivot_longer(
    cols= starts_with("DD_"),
    names_to = "year", 
    names_prefix = "DD_",
    values_to = "dd_year")%>% 
  mutate(year = as.integer(as.character(year)))

dat<-left_join(lmb_data, dd_temp, by = c("nhdid", "year")) %>% 
  left_join(surface_temp, by = c("nhdid", "year"))

## standardize and transform predictor values ##
dat$z_lake_area<-as.numeric(scale(log(dat$lake_area_m2))) #good 
dat$z_max_depth<-as.numeric(scale(log(dat$maxdepth_m))) #good 
dat$z_surface_temp_year<-as.numeric(scale(log(dat$surf_temp_year))) #good 
dat$z_secchi<-as.numeric(scale(log(dat$secchi_m)))
dat$z_bottom_do<-as.numeric(scale(log(dat$bottom_do_mgl+ 0.001))) ##has 0s so added 0.001 
dat$z_doy<-as.numeric(scale(dat$day_of_year)) ##standardize julian date 
dat$z_ws_forest<-as.numeric(scale(asin(sqrt(dat$ws_forest_prop)))) #good
dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat$ws_wetland_prop)))) #good

dat$logeffort <- log(dat$effort_new)

dat<-dplyr::select(dat, new_key, fish_count_new, logeffort, gear2, FMU_Code, 
                   z_lake_area, z_max_depth, z_secchi, z_doy, 
                  z_surface_temp_year,
                   z_ws_forest, z_ws_wetland) %>%
  na.omit()
#note 765 obs and 214 lakes with in-situ, without 1255 obs and 361 lakes  
n_distinct(dat$new_key) 

#check out response: log of count is mostly normally dist
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
    log.lambda[i] <- alpha[site[i]] + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i] + b[4] * x4[i] + 
                       b[5] * x5[i] + b[6] * x6[i] + b[7] * x7[i] +
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
    r ~ dunif(0,100)
    
    #derived quantities 
    tau.alpha <- pow(sigma.alpha,-2)
    tau.alpha2 <- pow(sigma.alpha2,-2)
    
    # priors for predictors 
    for(k in 1:7){
    b[k] ~ dnorm(0,0.0000001)
    }
    
     # priors for gear-specific log-catchabilities 
    for (k in 1:ngears) { 
    logq[k] ~ dnorm(0,  0.0000001)   #non-informative
    } 

                      
    } # end model
    ",fill = TRUE)
sink()


### set parameters ####
# Initial values
inits <- function (){
  list (mu.alpha = rnorm(8))
}


# Parameters monitored
parameters <- c("alpha", "mu.alpha","tau.alpha", "mu.alpha2","tau.alpha2","b", 'logq', 'r')


# MCMC settings
ni <- 100000 # number of iterations 
nt <- 10 #number to thin
nb <- 30000 #number to burn  
nc <- 3 #number of chains

#### split data into training and validation sets for cross-fold validation #### 
part1<-partition(dat$new_key, p=c(train=0.8, test = 0.2), seed = 1,  type =c( "grouped")) 
train1<-dat[part1$train,] 
test1<-dat[part1$test,]

part2<-partition(dat$new_key, p=c(train=0.8, test = 0.2), seed = 2, type =c( "grouped"))
train2<-dat[part2$train,]
test2<-dat[part2$test,]

part3<-partition(dat$new_key, p=c(train=0.8, test = 0.2), seed = 3, type =c( "grouped"))
train3<-dat[part3$train,]
test3<-dat[part3$test,]

part4<-partition(dat$new_key, p=c(train=0.8, test = 0.2), seed = 4, type =c( "grouped"))
train4<-dat[part4$train,]
test4<-dat[part4$test,]

part5<-partition(dat$new_key, p=c(train=0.8, test = 0.2), seed = 5, type =c( "grouped"))
train5<-dat[part5$train,]
test5<-dat[part5$test,]

#### set up data #### 

#set the number of sites and create site index 
nsites<-length(unique(train5$new_key))
site <- as.numeric(as.factor(train5$new_key))

# Set the number of FMUs 
J <- length(unique(train5$FMU_Code))

# Create group index, must go from 1,..J
train5$FMU_Code <- droplevels(as.factor(train5$FMU_Code))
train5$group <- as.numeric(as.factor(train5$FMU_Code))
group <- doBy::summaryBy(group ~ new_key, data=train5, FUN=mean) #need region to be indexed by site

#set number of different sampling gears used and set a gear index 
ngears<-length(unique(train5$gear2))
gear <- as.numeric(as.factor(train5$gear2))

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
train5$IND<-ifelse(train5$gear2 == "FT_NET", 0, 1) #FT_NET is the reference gear

# Load data
data <- list(y = train5$fish_count_new, group = group$group.mean, gear=gear, n = dim(train5)[1], J = J, ngears = ngears,
             x1=train5$z_secchi, x2=train5$z_lake_area, x3=train5$z_surface_temp_year, x4=train5$z_max_depth,
             x5=train5$z_ws_forest,x6=train5$z_ws_wetland, x7=train5$z_doy,
             logeffort=train5$logeffort, IND=train5$IND, nsites = nsites, site =site
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
#saveRDS(output, "Data/output/output_model1_secchi_fold5_lakes.rds") 
output<-readRDS("Data/output/output_model1_secchi_fold3_lakes.rds")

#### Summarize posteriors ####
print(output, dig = 3)
#in this case, logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
#check convergence with Brooks-Gelman-Rubin statistic
which(output$BUGSoutput$summary[, c("Rhat")] > 1.1)
which(output$BUGSoutput$summary[, c("n.eff")] < 0.1*21000) #times number of retained draws 
which(output$BUGSoutput$summary[, c("n.eff")] < 500)
jagsfit.mcmc <- as.mcmc(output)
model1tranformed <- ggs(jagsfit.mcmc) 
ggs_traceplot(model1tranformed, family = "b")
ggs_traceplot(model1tranformed, family = "r")

#### effect plots #### 
library(dotwhisker)

#betas
bEst <- matrix(NA, nrow=7,ncol=3)
for(i in 1:7){ #parameters
  bEst[i,1] <- mean(output$BUGSoutput$sims.list$b[,i])
  bEst[i,2:3] <- quantile(output$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
bEst<-as.data.frame(bEst)
bEst$variable<-c("Secchi depth",  "lake area", "surface temp", "max depth",  
                 "ws forest", "ws wetland", "day-of-year")

colnames(bEst)<- c("estimate", "conf.low", "conf.high", "term")
bEst$Parameter <-c("b[1]", "b[2]", "b[3]", "b[4]","b[5]","b[6]", "b[7]")


#add colors for sig different than 0
bEst$color <- as.numeric(bEst[,2] * bEst[,3] > 0 )

# plot using the dotwhisker package 
lake_plot<-dwplot(bEst, style = "dotwhisker", 
                  dot_args = list(aes(colour = factor(color))),
                  whisker_args = list(aes(colour = factor(color))), 
                  vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Variable") +
  ggtitle("lakes") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
lake_plot

#PLOT WITH DISTRIBUTIONS 
model_1_covariates<-model1tranformed %>%
  filter(grepl("b",Parameter))%>%
  left_join(bEst) %>% #get parameter names and colors 
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.05, .95)) + 
  theme_bw() + 
  theme(legend.position = "none") +
  ylab("predictor")

ggsave(plot=model_1_covariates, 
       device = "png", 
       filename = "figures/model_1_covariates.png", 
       dpi = 600, height = 8, width = 5, units = "in")

#### MODEL FIT #### 
# Simulating data from the posterior predictive distribution using the observed predictors is useful for checking the fit of the model.
#obtain our samples from the posterior
#this has a posterior dist (21000) for coef for every param; e.g. 255 lakes, 7 betas, 4 logqs, 8 mu.alphas, etc. 
coefs <- output$BUGSoutput$sims.matrix[,1:228] 

### 
# Number of desired MCMC samples used for prediction (use a subset of all samples)
nsim <- 3000
# Chain length from analysis
chainLength <- output$BUGSoutput$n.sims
# Select thinned steps in chain for posterior predictions to ensure we take values from length of posterior
ID = seq( 1 , chainLength , floor(chainLength/nsim) )

#need to add management unit index and gear index
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
test1<-test1 %>% 
  mutate(group = case_when( FMU_Code == "MI-MCM" ~ 1, 
                            FMU_Code == "MI-MER" ~ 2, 
                            FMU_Code == "MI-MES" ~ 3, 
                            FMU_Code == "MI-MNH" ~ 4,  
                            FMU_Code == "MI-MNM" ~ 5, 
                            FMU_Code== "MI-MSH" ~ 6,
                            FMU_Code== "MI-MSM" ~ 7, 
                            FMU_Code== "MI-MWS" ~ 8)
  )%>% 
  mutate(gear = case_when(gear2 == "FT_NET" ~ 1,
                          gear2 == "GILL" ~ 2,
                          gear2 == "SEINE" ~ 3, 
                          gear2 == "SHOCK" ~ 4)
  )

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
test1$IND<-ifelse(test1$gear2 == "FT_NET", 0, 1) #FT_NET is the reference gear

# Container for predicted values
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
predictions_sim<- array(NA, c(nsim,length(test1[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predictions_sim)

for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(test1[[1]])){ #loop over covariate data (obs)
    predictions_sim[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',test1[['group']][t],']')] + coefs[ID[i],'b[1]']*test1[['z_secchi']][t]  + coefs[ID[i],'b[2]']*test1[['z_lake_area']][t]  + coefs[ID[i],'b[3]']*test1[['z_surface_temp_year']][t]  + 
                                  coefs[ID[i],'b[4]']*test1[['z_max_depth']][t]  + coefs[ID[i],'b[5]']*test1[['z_ws_forest']][t] + coefs[ID[i],'b[6]' ]*test1[['z_ws_wetland']][t] + coefs[ID[i],'b[7]' ]*test1[['z_doy']][t] + 
                                  coefs[ID[i], paste0('logq[',test1[['gear']][t],']')]*test1[['IND']][t] + test1[['logeffort']][t]) # exp the predictions because we used a log-link in the Poisson
    
  }
}

#container for p values 
#parametrization is by the dispersion parameter, where prob = size/(size+mu).
#p is the success parameter and r is the dispersion parameter.
p <- array(NA, c(nsim,length(test1[[1]])) )
exp_catch <- array(NA, c(nsim,length(test1[[1]])) )

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

test1$site <- as.numeric(as.factor(test1$new_key))
obs_catch<-dplyr::select(test1, site, new_key, gear2, fish_count_new) %>% 
  mutate(row = row_number())
plot_data<-left_join(med_data, obs_catch) %>%
  mutate(gear = case_when(       #rename gear types for publication
    gear2 == "FT_NET" ~ 'fyke net', 
    gear2 == "GILL" ~ 'gillnet', 
    gear2 == "SEINE" ~ 'seine', 
    gear2 == "SHOCK" ~ 'shock')) 

plot( med_data$col.med , obs_catch$fish_count_new)
plot( med_data$col.means, obs_catch$fish_count_new)

#which can be compared to the observed catch for that combination of lake and gear. 
#One suggestion for visualizing the model fit is to generate a separate graph for each gear, 
#plot the mean predicted catches (+- the 95% ci) of different lakes on the x-axis versus the observed catches of lakes on the y-axis.
pred_plot<-ggplot()+
  geom_pointrange(data=plot_data, aes(x=col.means, y=fish_count_new, xmax=upperCI.Group, xmin=lowerCI.Group))+ 
  xlab('predicted catch')+
  ylab('observed catch')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)

pred_plot+facet_wrap(~ gear, ncol=2, scales = 'free') #allow scales to vary 

#plot log scale 
pred_plot_log<-ggplot(data=plot_data, aes(x=log(col.means +1), y=log(fish_count_new +1), xmax=log(upperCI.Group +1), xmin=log(lowerCI.Group +1) ) )+
  geom_pointrange( color="grey" )+ 
  geom_point(color="black")+
  xlab('predicted catch (log)')+
  ylab('observed catch (log)')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)

model_1_pred_obs<-pred_plot_log+facet_wrap(~ gear, ncol=2, scales = 'free') #allow scales to vary 

ggsave(plot=model_1_pred_obs, 
       device = "png", 
       filename = "figures/model_1_pred_obs.png", 
       dpi = 600, height = 8, width = 12, units = "in")
###  plot model fits  #### 
library(bayesplot)
#yrep<-as.matrix(output$BUGSoutput$sims.list$ysim ) #from the predictions within JAGS
yrep<-exp_catch
y<-test1$fish_count_new

color_scheme_set("brightblue")
ppc_dens_overlay(y, yrep[1:50, ])
ppc_dens_overlay(y, yrep[1:50, ]) + xlim(0, 50)
#ppc_hist(y, yrep[1:5, ])
#group 
model_1_dens<-ppc_dens_overlay_grouped(y, yrep[1:50, ], group = test5$gear2) + xlim(0, 20)

ggsave(plot=model_1_dens, 
       device = "png", 
       filename = "figures/model_1_dens.png", 
       dpi = 600, height = 8, width = 12, units = "in")

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

# RMSE 
sqrt(mean((plot_data$col.means-plot_data$fish_count_new)^2))

#residual plot
#https://cran.r-project.org/web/packages/tidybayes/vignettes/tidybayes-residuals.html
plot_data$res.med <- apply(e, 2, median )
res_plot<-plot_data %>%
  ggplot(aes(x = col.med, y = res.med)) +
  stat_pointinterval()
res_plot
res_plot+facet_wrap(~ gear, ncol=2, scales = 'free')
#plot log scale 
res_plot_log<-plot_data %>%
  ggplot(aes(x = log(col.med+1), y = log(res.med+1))) +
  stat_pointinterval()
res_plot_log

res_plot_log+facet_wrap(~ gear, ncol=2, scales = 'free')

#qq line plot 
plot_data  %>%
  ggplot(aes(sample = res.med)) +
  geom_qq() +
  geom_qq_line()

#### PREDICTION to Historical data ##### 
hist_dat<-read.csv("Data/historical_lmb2_with_secchi.csv")%>% 
  drop_na(secchi_combo) #only use samples with secchi from contemp data 

#set up data 
hist_dat$logeffort<-log(hist_dat$effort_sum)
## standardize and transform predictor values ##
hist_dat$z_lake_area<-as.numeric(scale(log(hist_dat$lake_area_m2))) #good 
hist_dat$z_max_depth<-as.numeric(scale(log(hist_dat$max_depth_m))) #good 
hist_dat$z_dd_mean<-as.numeric(scale(log(hist_dat$dd_year))) #good
hist_dat$z_surface_temp_year<-as.numeric(scale(log(hist_dat$surface_temp_year))) #good 
hist_dat$z_secchi<-as.numeric(scale(log(hist_dat$secchi_combo)))
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

#### *multiple test sets for each fold ####
#output<-readRDS("Data/output/output_model1_secchi_lmb_fold1.rds")
#coefs <- output$BUGSoutput$sims.matrix[,1:277] 

# Number of desired MCMC samples used for prediction (use a subset of all samples)
nsim <- 3000
# Chain length from analysis
chainLength <- output$BUGSoutput$n.sims
# Select thinned steps in chain for posterior predictions to ensure we take values from length of posterior
ID = seq( 1 , chainLength , floor(chainLength/nsim) )

#pull out about 50 lakes to match contemp test dataset 
hist1<-partition(hist_dat$new_key, p=c(train=0.4,test=0.6), seed = 1, type =c("grouped"))
hist_test1<-hist_dat[hist1$test,]
hist2<-partition(hist_dat$new_key, p=c(train=0.4,test=0.6), seed = 2, type =c("grouped"))
hist_test2<-hist_dat[hist2$test,]
hist3<-partition(hist_dat$new_key, p=c(train=0.4,test=0.6), seed = 3, type =c("grouped"))
hist_test3<-hist_dat[hist3$test,]
hist4<-partition(hist_dat$new_key, p=c(train=0.4,test=0.6), seed = 4, type =c("grouped"))
hist_test4<-hist_dat[hist4$test,]
hist5<-partition(hist_dat$new_key, p=c(train=0.4,test=0.6), seed = 5, type =c("grouped"))
hist_test5<-hist_dat[hist5$test,]

hist_test<-hist_test5# just change this input each time 
n_distinct(hist_test$new_key)

#* start prediction #### 
# Container for predicted values
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
predict_hist<- array(NA, c(nsim,length(hist_test[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predict_hist)

#use regional mu.alphas instead of site specific, since dif sites 
for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(hist_test[[1]])){ #loop over covariate data (obs)
    predict_hist[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',hist_test[['group']][t],']')] + coefs[ID[i],'b[1]']*hist_test[['z_secchi']][t]  + coefs[ID[i],'b[2]']*hist_test[['z_lake_area']][t]  + coefs[ID[i],'b[3]']*hist_test[['z_surface_temp_year']][t]  + 
                               coefs[ID[i],'b[4]']*hist_test[['z_max_depth']][t]  + coefs[ID[i],'b[5]' ]*hist_test[['z_ws_forest']][t] + coefs[ID[i],'b[6]']*hist_test[['z_ws_wetland']][t] + coefs[ID[i],'b[7]']*hist_test[['z_julian']][t] +
                               coefs[ID[i], paste0('logq[',hist_test[['gear_index']][t],']')]*hist_test[['IND']][t] + hist_test[['logeffort']][t]) # exp the predictions because we used a log-link in the Poisson
    
  }
}


#container for p values 
#parametrization is by the dispersion parameter, where prob = size/(size+mu).
#p is the success parameter and r is the dispersion parameter.
p_hist <- array(NA, c(nsim,length(hist_test[[1]])) )
exp_catch_hist <- array(NA, c(nsim,length(hist_test[[1]])) )

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

hist_test$site <- as.numeric(as.factor(hist_test$new_key))
hist_catch<-dplyr::select(hist_test, site, gear, largemouthbass_sum) %>% 
  mutate(row = row_number())
plot_data<-left_join(med_data, hist_catch) 

plot( med_data$hist.med , hist_catch$largemouthbass_sum)
plot( med_data$hist.means, hist_catch$largemouthbass_sum)

#remove fyke before plotting 
plot_data_no_fyke=filter(plot_data, gear != 'FT_net')
#which can be compared to the observed catch for that combination of lake and gear. 
#One suggestion for visualizing the model fit is to generate a separate graph for each gear, 
#plot the mean or median predicted catches (+- the 95% ci) of different lakes on the x-axis versus the observed catches of lakes on the y-axis.
pred_plot<-ggplot()+
  geom_pointrange(data=plot_data_no_fyke, aes(x=hist.means, y=largemouthbass_sum, xmax=upperPI, xmin=lowerPI))+ 
  xlab('predicted catch')+
  ylab('historical catch')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)

pred_plot+facet_wrap(~ gear, ncol=2, scales = 'free') #allow scales to vary 

#plot log scale 
pred_plot_log<-ggplot(data=plot_data_no_fyke, aes(x=log(hist.means+1), y=log(largemouthbass_sum + 1), xmax=log(upperPI +1), xmin=log(lowerPI +1) ) )+
  geom_pointrange( color="grey" )+ 
  geom_point(color="black")+
  xlab('predicted catch (log)')+
  ylab('historical catch (log)')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)

model_1_hist_pred<-pred_plot_log+facet_wrap(~ gear, ncol=2, scales = 'free') #allow scales to vary 

ggsave(plot=model_1_hist_pred, 
       device = "png", 
       filename = "figures/model_1_hist_pred.png", 
       dpi = 600, height = 8, width = 12, units = "in")

#### RESIDUALS #### 
#Bayes-R2 using residuals
#https://avehtari.github.io/bayes_R2/bayes_R2.html#2_Functions_for_Bayesian_R-squared_for_stan_glm_models
ypred<-exp_catch_hist
y<-plot_data$largemouthbass_sum
e <- -1 * sweep(ypred, 2, y) #residual distribution 
var_ypred <- apply(ypred, 1, var) # variance of modelled predictive means
var_e <- apply(e, 1, var) #variance of modeled residuals 
bayesR2res<- var_ypred / (var_ypred + var_e) #R^2
round(median(bayesR2res), 2) 

sqrt(mean((plot_data$hist.means-plot_data$largemouthbass_sum)^2))

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

#### relative abundance hist ####
#extract just lake data (need a single row for each lake) but also need pres/abs of walleye and pike in a lake, not by gear
hist_just_lake<-hist_dat %>% 
  distinct(new_key, .keep_all = TRUE) #83 lakes 

# Container for predicted values
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
predict_hist<- array(NA, c(nsim,length(hist_just_lake[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predict_hist)

#use regional mu.alphas instead of site specific, since dif sites, and environment, no gear data 
for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(hist_just_lake[[1]])){ #loop over covariate data (obs)
    predict_hist[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',hist_just_lake[['group']][t],']')] + coefs[ID[i],'b[1]']*hist_just_lake[['z_secchi']][t]  + coefs[ID[i],'b[2]']*hist_just_lake[['z_lake_area']][t]  + coefs[ID[i],'b[3]']*hist_just_lake[['z_surface_temp_year']][t]  + 
                               coefs[ID[i],'b[4]']*hist_just_lake[['z_max_depth']][t]  + coefs[ID[i],'b[5]' ]*hist_just_lake[['z_ws_forest']][t] + coefs[ID[i],'b[6]']*hist_just_lake[['z_ws_wetland']][t] +
                               coefs[ID[i],'b[7]']*hist_just_lake[['z_julian']][t] ) 
    
    
  }
}


#container for p values 
#parametrization is by the dispersion parameter, where prob = size/(size+mu).
#p is the success parameter and r is the dispersion parameter.
p_hist <- array(NA, c(nsim,length(hist_just_lake[[1]])) )
exp_abund_hist <- array(NA, c(nsim,length(hist_just_lake[[1]])) )

for(i in 1:nrow(predict_hist)){ # loop over rows (sims)
  for(t in 1:ncol(predict_hist) ){ #loop over columns(obs)
    p_hist[i,t] <- coefs[ID[i],'r']/(coefs[ID[i],'r'] + predict_hist[i,t])
    exp_abund_hist[i,t] <- rnbinom(1, prob=p_hist[i,t], size= coefs[ID[i],'r'])
  }
}

#From this distribution, you can extract the median and the 2.5% and 97.5% percentiles (95% credible interval),
hist.med <- apply(exp_abund_hist, 2, median )
hist.means <- apply(exp_abund_hist, 2, mean )

#prep data for plotting
abund_data=data.frame(hist.med, hist.means) 
abund_data$row <- as.numeric(row.names(abund_data))

hist_just_lake$site <- as.numeric(as.factor(hist_just_lake$new_key))
hist_ids<-dplyr::select(hist_just_lake, site, new_key) %>% 
  mutate(row = row_number())
map_data<-left_join(abund_data, hist_ids) 

#### relative abundance contemporary ####
#extract just lake data (need a single row for each lake) but also need pres/abs of walleye and pike in a lake, not by gear
cont_just_lake<-test1 %>% 
  distinct(new_key, .keep_all = TRUE) #135 test 1 subset 

# Container for predicted values
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
predict_cont<- array(NA, c(nsim,length(cont_just_lake[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predict_cont)

#use lake/site specific alphas and environment, no gear data 
for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(cont_just_lake[[1]])){ #loop over covariate data (obs)
    predict_cont[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',cont_just_lake[['group']][t],']')] + coefs[ID[i],'b[1]']*cont_just_lake[['z_secchi']][t]  + coefs[ID[i],'b[2]']*cont_just_lake[['z_lake_area']][t]  + coefs[ID[i],'b[3]']*cont_just_lake[['z_surface_temp_year']][t]  + 
                               coefs[ID[i],'b[4]']*cont_just_lake[['z_max_depth']][t]  + coefs[ID[i],'b[5]' ]*cont_just_lake[['z_ws_forest']][t] + coefs[ID[i],'b[6]']*cont_just_lake[['z_ws_wetland']][t] +
                               coefs[ID[i],'b[7]']*cont_just_lake[['z_doy']][t] ) 
    
    
  }
}


#container for p values 
#parametrization is by the dispersion parameter, where prob = size/(size+mu).
#p is the success parameter and r is the dispersion parameter.
p_cont <- array(NA, c(nsim,length(cont_just_lake[[1]])) )
exp_abund_cont <- array(NA, c(nsim,length(cont_just_lake[[1]])) )

for(i in 1:nrow(predict_cont)){ # loop over rows (sims)
  for(t in 1:ncol(predict_cont) ){ #loop over columns(obs)
    p_cont[i,t] <- coefs[ID[i],'r']/(coefs[ID[i],'r'] + predict_cont[i,t])
    exp_abund_cont[i,t] <- rnbinom(1, prob=p_cont[i,t], size= coefs[ID[i],'r'])
  }
}

#From this distribution, you can extract the median and the 2.5% and 97.5% percentiles (95% credible interval),
cont.med.abund <- apply(exp_abund_cont, 2, median )
cont.means.abund <- apply(exp_abund_cont, 2, mean )

#prep data for plotting
cont_abund_data=data.frame(cont.med.abund, cont.means.abund) 
cont_abund_data$row <- as.numeric(row.names(cont_abund_data))

cont_just_lake$site <- as.numeric(as.factor(cont_just_lake$new_key))
cont_ids<-dplyr::select(cont_just_lake, site, new_key) %>% 
  mutate(row = row_number())
cont_map_data<-left_join(cont_abund_data, cont_ids) 



#*map predicted abundance ####
library(ggmap)
library(stringr)
library(devtools)
library(mapdata)
library(dplyr)

#get the lat/long points - need lat/lon from Humphries table 
points<-read.csv("Data/MI_data/Humphries_table.csv") %>%
  dplyr::select( New_Key, LAT_DD, LONG_DD) %>%
  rename(new_key = New_Key) 

hist_abund_map_dat<-left_join(map_data, points)
cont_abund_map_dat<-left_join(cont_map_data, points)

summary(hist_abund_map_dat$hist.means)
summary(cont_abund_map_dat$cont.means.abund)

#get MI basemap 
MI_basemap<-map_data("state") %>%
  subset(region %in% c("michigan")) # select michigan 
p<-ggplot(data = MI_basemap) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + #this fill MI white and outline is black
  coord_fixed(1.3) 
#map of lmb abund historical 
hist_map<-p+ geom_point(data=hist_abund_map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(hist.means)), alpha=0.7) + 
  scale_colour_gradientn(colours = c("tan", "green","skyblue", "midnightblue"),
                         )+ 
  labs(color="lmb relative density", tag = "a") + #changes the labels on the legend 
  theme(legend.position = "bottom", 
        plot.tag = element_text(), 
        plot.tag.position = c(0.1,0.95)) +
  ylab(NULL) + xlab(NULL)

hist_map
#log
p+ geom_point(data=hist_abund_map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(log(hist.means + 1)))) + 
  scale_colour_gradient(low = "yellow", high = "darkblue")  + 
  labs(color="lmb log abund")  #changes the labels on the legend 

#map of lmb abund contemporary 
cont_map<-p+ geom_point(data=cont_abund_map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(cont.means.abund)), alpha=0.7) + 
  scale_colour_gradientn(colours = c("tan", "green","skyblue", "midnightblue"),
                         )+ 
  labs(color="lmb relative density", tag= "b")  +  #changes the labels on the legend
  theme(legend.position = "bottom", 
        plot.tag = element_text(), 
        plot.tag.position = c(0.1,0.95)) + 
  ylab(NULL) + xlab(NULL)

cont_map
#log
p+ geom_point(data=cont_abund_map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(log(cont.means.abund + 1)))) + 
  scale_colour_gradient(low = "yellow", high = "darkblue")  + 
  labs(color="lmb log abund") +  #changes the labels on the legend
  theme(legend.position = "bottom") + ylab(NULL) + xlab(NULL)

#* save as 2-panel plot 
map_model1<-cowplot::plot_grid(hist_map, cont_map)

ggsave(plot=map_model1, 
       device = "png", 
       filename = "figures/model_1_map_abund.png", 
       dpi = 600, height = 8, width = 12, units = "in",
       bg="#ffffff") #sets background to white 

#other stats 
mean(output$BUGSoutput$sims.list$b[,1] > 0) # Posterior probability that Secchi effect is positive 
quantile(output$BUGSoutput$sims.list$b[,1],c(0.05,0.95))
mean(output$BUGSoutput$sims.list$b[,1])

#save predicted abundance in each lake 
write.csv(hist_abund_map_dat, "Data/output/model1_hist_abund_pred.csv", row.names = FALSE)
write.csv(cont_abund_map_dat, "Data/output/model1_cont_abund_pred.csv", row.names = FALSE)

model2_hist_est<-read.csv("Data/output/model2_hist_abund_pred.csv")
model2_cont_est<-read.csv("Data/output/model2_cont_abund_pred.csv")
