## Negative Binomial gear-general model #### 

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
#lmb_dat has count by gear and effort by gear
lmb_data<-read.csv("Data/contemp_lmb_dat.csv") 

#surface temp by year 
surface_temp<-read.csv("Data/lake_surface_temp.csv")%>%
  group_by(IHDLKID) %>%
  slice_max(HECTARES) %>% 
  rename(nhdid = IHDLKID)%>% 
  pivot_longer(
    cols= starts_with("TAVE_"),
    names_to = "year", 
    names_prefix = "TAVE_",
    values_to = "surf_temp_year") %>% 
  mutate(year = as.integer(as.character(year)))


dat<-left_join(lmb_data, surface_temp, by = c("nhdid", "year")) 

## standardize and transform predictor values ##
dat$z_lake_area<-as.numeric(scale(log(dat$lake_area_m2))) 
dat$z_max_depth<-as.numeric(scale(log(dat$maxdepth_m))) 
dat$z_surface_temp_year<-as.numeric(scale(log(dat$surf_temp_year))) 
dat$z_secchi<-as.numeric(scale(log(dat$secchi_m)))
dat$z_bottom_do<-as.numeric(scale(log(dat$bottom_do_mgl+ 0.001))) ##has 0s so added 0.001 
dat$z_doy<-as.numeric(scale(dat$day_of_year)) ##standardize doy date 
dat$z_ws_forest<-as.numeric(scale(asin(sqrt(dat$ws_forest_prop)))) 
dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat$ws_wetland_prop)))) 

dat$logeffort <- log(dat$effort_new)

dat<-dplyr::select(dat, new_key, fish_count_new, logeffort, gear2, FMU_Code, 
                   z_lake_area, z_max_depth, z_secchi, z_doy, 
                  z_surface_temp_year,
                   z_ws_forest, z_ws_wetland) %>%
  na.omit()

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
output<-readRDS("Data/output/output_model1_secchi_fold5_lakes.rds")

#### Summarize posteriors ####
print(output, dig = 3)
#in this case, logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
#check convergence with Brooks-Gelman-Rubin statistic and traceplots 
which(output$BUGSoutput$summary[, c("Rhat")] > 1.1)
which(output$BUGSoutput$summary[, c("n.eff")] < 500)
jagsfit.mcmc <- as.mcmc(output)
model1tranformed <- ggs(jagsfit.mcmc) 
ggs_traceplot(model1tranformed, family = "b")
ggs_traceplot(model1tranformed, family = "r")

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

cont_test<-test5
#need to add management unit index and gear index
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
cont_test<-cont_test %>% 
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
  ) %>% 
  mutate(IND = ifelse(gear2 == "FT_NET", 0, 1))#create an indicator variable “IND”, with value 0 for every sample using the reference gear (FT_NET)


# predict catch values for each gear/lake 
predictions_sim<- array(NA, c(nsim,length(cont_test[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predictions_sim)

for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(cont_test[[1]])){ #loop over covariate data (obs)
    predictions_sim[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',cont_test[['group']][t],']')] + coefs[ID[i],'b[1]']*cont_test[['z_secchi']][t]  + coefs[ID[i],'b[2]']*cont_test[['z_lake_area']][t]  + coefs[ID[i],'b[3]']*cont_test[['z_surface_temp_year']][t]  + 
                                  coefs[ID[i],'b[4]']*cont_test[['z_max_depth']][t]  + coefs[ID[i],'b[5]']*cont_test[['z_ws_forest']][t] + coefs[ID[i],'b[6]' ]*cont_test[['z_ws_wetland']][t] + coefs[ID[i],'b[7]' ]*cont_test[['z_doy']][t] + 
                                  coefs[ID[i], paste0('logq[',cont_test[['gear']][t],']')]*cont_test[['IND']][t] + cont_test[['logeffort']][t]) # exp the predictions because we used a log-link in the Poisson
    
  }
}

#*predicted catch distribution #### 
#extract the median and the 2.5% and 97.5% percentiles (95% credible interval),
col.med <- apply(predictions_sim, 2, median )
col.means <- apply(predictions_sim, 2, mean )
# 95% CIs for fitted values
upperCI.Group <- apply(predictions_sim, 2, quantile, probs=c(0.975) )
lowerCI.Group <- apply(predictions_sim, 2, quantile, probs=c(0.025) )

#prep data for plotting
med_data=data.frame(col.med, col.means, upperCI.Group, lowerCI.Group) 
med_data$row <- as.numeric(row.names(med_data))

cont_test$site <- as.numeric(as.factor(cont_test$new_key))
obs_catch<-dplyr::select(cont_test, site, new_key, gear2, fish_count_new) %>% 
  mutate(row = row_number())
plot_data<-left_join(med_data, obs_catch) %>%
  mutate(gear = case_when(       #rename gear types for publication
    gear2 == "FT_NET" ~ 'fyke net', 
    gear2 == "GILL" ~ 'gillnet', 
    gear2 == "SEINE" ~ 'seine', 
    gear2 == "SHOCK" ~ 'shock')) 

plot( med_data$col.med , obs_catch$fish_count_new)
plot( med_data$col.means, obs_catch$fish_count_new)

#*get deviance and residual deviance #### 
#p is the success parameter (where prob = size/(size+mu)) and r is the dispersion parameter.
#calculate the log-likelihood (L) for each fold (each observation gets a prob and then sum the obs)
p <- array(NA, c(nsim,length(cont_test[[1]])) ) #container for p values 
likeli_cont <- array(NA, c(nsim,length(cont_test[[1]])) )

for(i in 1:nrow(predictions_sim)){ # loop over rows (sims)
  for(t in 1:ncol(predictions_sim) ){ #loop over columns(obs)
    p[i,t] <- coefs[ID[i],'r']/(coefs[ID[i],'r'] + predictions_sim[i,t])
    likeli_cont[i,t] = -2*(log(dnbinom(cont_test[['fish_count_new']][t], mu=p[i,t], size= coefs[ID[i],'r']))) #deviance for each obs
  }
}

#*#deviance 
#need to sum across columns (observations) to get a dist of the deviance with 3000 values and then take mean
likeli_cont[is.infinite(likeli_cont)] <- NA
deviance<-rowSums(likeli_cont, na.rm = TRUE)
mean(deviance, na.rm = TRUE)

#* deviance residuals 
dev_resid_cont <- data.frame(matrix(ncol=1, nrow =c(length(cont_test[[1]]))) ) # will hold all of the 5 test sets for a fold 
colnames(dev_resid_cont) <- c("deviance_res")

#From the posterior distribution of deviance get the mean for each observation 
plot_data$dev_means <- apply(likeli_cont, 2, mean,  na.rm = TRUE )

#add deviance residuals
for(i in 1:nrow(plot_data)){ # loop over rows
  dev_resid_cont[i,1] = sign(plot_data[['fish_count_new']][i] - plot_data[['col.means']][i]) * sqrt(plot_data[['dev_means']][i])
}

#write residuals for each fold
#write.csv(dev_resid_cont, "Data/output/model1_cont_res_dev_fold5.csv", row.names = FALSE)


#*compare the observed catch for that combination of lake and gear. ####
# visualize the model fit: generate a separate graph for each gear 
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
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed')

model_1_pred_obs<-pred_plot_log+facet_wrap(~ gear, ncol=2, scales = 'free') #allow scales to vary 
model_1_pred_obs

ggsave(plot=model_1_pred_obs, 
       device = "png", 
       filename = "figures/fig2_model_1_pred_obs.png", 
       dpi = 600, height = 8, width = 12, units = "in")



#*residuals for contemporary model  

#Bayes-R2 using residuals
e <- -1 * sweep(predictions_sim, 2, y)
var_ypred <- apply(predictions_sim, 1, var)
var_e <- apply(e, 1, var)
bayesR2res<- var_ypred / (var_ypred + var_e)
round(median(bayesR2res), 2) 

# RMSE 
sqrt(mean((plot_data$col.means-plot_data$fish_count_new)^2))


#### PREDICTION to Historical data ##### 
hist_dat<-read.csv("Data/historical_lmb2_with_secchi.csv")%>% 
  drop_na(secchi_combo) #only use samples with secchi from contemp data 

#set up data 
hist_dat$logeffort<-log(hist_dat$effort_sum)
## standardize and transform predictor values ##
hist_dat$z_lake_area<-as.numeric(scale(log(hist_dat$lake_area_m2))) #good 
hist_dat$z_max_depth<-as.numeric(scale(log(hist_dat$max_depth_m))) #good 
hist_dat$z_surface_temp_year<-as.numeric(scale(log(hist_dat$surface_temp_year))) #good 
hist_dat$z_secchi<-as.numeric(scale(log(hist_dat$secchi_combo)))
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
  ) %>% 
  mutate( #need gear number to match original data (no shock) 
    gear_index = case_when( gear == "FT_net" ~ 1,
                            gear == "gill" ~ 2,
                            gear == "seine" ~ 3)
  ) %>% 
  mutate(IND = ifelse(gear == "FT_net", 0, 1)) #FT_NET is the reference gear

#### *multiple test sets for each fold ####
output<-readRDS("Data/output/output_model1_secchi_fold5_lakes.rds")
coefs <- output$BUGSoutput$sims.matrix[,1:228] 

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

#* start prediction #### 
hist_test<-hist_test5# just change this input each time 
n_distinct(hist_test$new_key)
# Container for predicted values
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

#*get residual deviance #### 
#p is the success parameter (prob = size/(size+mu) and r is the dispersion parameter.
p_hist <- array(NA, c(nsim,length(hist_test[[1]])) )
likeli_hist <- array(NA, c(nsim,length(hist_test[[1]])) )

for(i in 1:nrow(predict_hist)){ # loop over rows (sims)
  for(t in 1:ncol(predict_hist) ){ #loop over columns(obs)
    p_hist[i,t] <- coefs[ID[i],'r']/(coefs[ID[i],'r'] + predict_hist[i,t])
    likeli_hist[i,t] <- -2*(log(dnbinom(hist_test[['largemouthbass_sum']][t], mu=p_hist[i,t], size= coefs[ID[i],'r']))) #deviance 
  }
}

#*From the posterior distribution of expected catch prediction extract the mean ####
#and the 2.5% and 97.5% percentiles (95% credible interval),
hist.med <- apply(predict_hist, 2, median )
hist.means <- apply(predict_hist, 2, mean )
# 95% CIs for fitted values
upperPI <- apply(predict_hist, 2, quantile, probs=c(0.975) )
lowerPI <- apply(predict_hist, 2, quantile, probs=c(0.025) )

#prep data for plotting
med_data=data.frame(hist.med, hist.means, upperPI, lowerPI) 
med_data$row <- as.numeric(row.names(med_data))

hist_test$site <- as.numeric(as.factor(hist_test$new_key))
hist_catch<-dplyr::select(hist_test, site, gear, largemouthbass_sum) %>% 
  mutate(row = row_number())
plot_data<-left_join(med_data, hist_catch) 

#plot( med_data$hist.med , hist_catch$largemouthbass_sum)
#plot( med_data$hist.means, hist_catch$largemouthbass_sum)

#*residual deviance 
#From the posterior distribution of deviance get the mean for each observation 
plot_data$dev_means <- apply(likeli_hist, 2, mean )

#dev_resid <- data.frame(matrix(ncol=5, nrow =79) ) # will hold all of the 5 test sets for a fold 
#colnames(dev_resid) <- c("test1", "test2", "test3", "test4", "test5")

#add residual deviance of each test to a new column 
for(i in 1:nrow(hist_test)){ # loop over rows
  dev_resid[i,5] = sign(plot_data[['largemouthbass_sum']][i] - plot_data[['hist.means']][i]) * sqrt(plot_data[['dev_means']][i])
}

#write.csv(dev_resid, "Data/output/model1_hist_res_dev_fold5.csv", row.names = FALSE)

#remove fyke before plotting 
plot_data_no_fyke=filter(plot_data, gear != 'FT_net')
#which can be compared to the observed catch for that combination of lake and gear. 
#visualize the model fit: generate a separate graph for each gear 
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


#### historical Bayes-R2 using residuals #### 
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

