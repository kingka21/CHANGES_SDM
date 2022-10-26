#### generate observed vs predicted plot with contemp and historical data #### 

#use model 2 - better model 

#####read in model output ####
output<-readRDS("Data/output/output_model2_secchi_fold5_lakes.rds")

##### read in and partition data ####
#lmb_dat_for_model has count by gear and effort by gear
dat<-read.csv("Data/lmb_dat_for_model_aug30.csv") 
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

dat<-left_join(dat, dd_temp, by = c("nhdid", "year")) %>% 
  left_join(surface_temp, by = c("nhdid", "year"))

## standardize and transform predictor values ##
#dat$z_order<-as.numeric(scale(log(dat$lake_order+0.001))) # min is 0 because these are the isolated
dat$z_lake_area<-as.numeric(scale(log(dat$lake_area_m2))) #good 
dat$z_max_depth<-as.numeric(scale(log(dat$maxdepth_m))) #good 
dat$z_dd_year<-as.numeric(scale(log(dat$dd_year))) #good
dat$z_surface_temp_year<-as.numeric(scale(log(dat$surf_temp_year))) #good 
dat$z_secchi<-as.numeric(scale(log(dat$secchi_m)))
dat$z_bottom_do<-as.numeric(scale(log(dat$bottom_do_mgl+ 0.001))) ##has 0s so added 0.001 
dat$z_doy<-as.numeric(scale(dat$day_of_year)) ##standardize julian date 
dat$z_ws_forest<-as.numeric(scale(asin(sqrt(dat$ws_forest_prop)))) #good
dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat$ws_wetland_prop)))) #good

dat$logeffort <- log(dat$effort_new)

#not using wae/pike or DO but use all available Secchi! 
dat<-dplyr::select(dat, new_key, fish_count_new, logeffort, gear2, FMU_Code, 
                   z_lake_area, z_max_depth, z_secchi, z_doy, 
                   z_surface_temp_year,
                   z_ws_forest, z_ws_wetland) %>%
  na.omit()

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

#### contemporary fit #### 
# Simulating data from the posterior predictive distribution using the observed predictors 
#obtain our samples from the posterior
coefs <- output$BUGSoutput$sims.matrix[,1:256] 
# Number of desired MCMC samples used for prediction (use a subset of all samples)
nsim <- 3000
# Chain length from analysis
chainLength <- output$BUGSoutput$n.sims
# Select thinned steps in chain for posterior predictions to ensure we take values from length of posterior
ID = seq( 1 , chainLength , floor(chainLength/nsim) )

cont_test<-test5# just change this input each time 
n_distinct(cont_test$new_key)

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
  )

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
cont_test$IND<-ifelse(cont_test$gear2 == "FT_NET", 0, 1) #FT_NET is the reference gear

# Container for predicted values using the test dataset 5
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
#x1=z_secchi, x2=z_lake_area, x3=z_surface_temp_year, x4=z_max_depth,
#x5=z_ws_forest,x6=z_ws_wetland, x7=z_doy,
predictions_sim<- array(NA, c(nsim,length(cont_test[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predictions_sim)

for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(cont_test[[1]])){ #loop over covariate data (obs) #use regional mu.alphas instead of site specific
    predictions_sim[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',cont_test[['group']][t],']')] + coefs[ID[i], paste0('b1[',cont_test[['gear']][t],']')]*cont_test[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',cont_test[['gear']][t],']')]*cont_test[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',cont_test[['gear']][t],']')]*cont_test[['z_surface_temp_year']][t]  + 
                                  coefs[ID[i],paste0('b4[',cont_test[['gear']][t],']')]*cont_test[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',cont_test[['gear']][t],']')]*cont_test[['z_ws_forest']][t] + coefs[ID[i],paste0('b6[',cont_test[['gear']][t],']')]*cont_test[['z_ws_wetland']][t] + coefs[ID[i],paste0('b7[',cont_test[['gear']][t],']')]*cont_test[['z_doy']][t] +
                                  coefs[ID[i], paste0('logq[',cont_test[['gear']][t],']')]*cont_test[['IND']][t] + cont_test[['logeffort']][t]) # exp the predictions because we used a log-link in the negbi
    
  }
}


#container for p values 
#parametrization is by the dispersion parameter, where prob = size/(size+mu).
#p is the success parameter and r is the dispersion parameter.
p <- array(NA, c(nsim,length(cont_test[[1]])) )
exp_catch <- array(NA, c(nsim,length(cont_test[[1]])) )

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
upperPI <- apply(exp_catch, 2, quantile, probs=c(0.975) )
lowerPI <- apply(exp_catch, 2, quantile, probs=c(0.025) )

#prep data for plotting
med_data=data.frame(col.med, col.means, upperPI, lowerPI) 
med_data$row <- as.numeric(row.names(med_data))

cont_test$site <- as.numeric(as.factor(cont_test$new_key))
obs_catch<-dplyr::select(cont_test, new_key, site, gear2, fish_count_new) %>% 
  mutate(row = row_number())
cont_plot_data<-left_join(med_data, obs_catch) %>% 
  mutate(gear = case_when(       #rename gear types for publication
    gear2 == "FT_NET" ~ 'fyke net', 
    gear2 == "GILL" ~ 'gillnet', 
    gear2 == "SEINE" ~ 'seine', 
    gear2 == "SHOCK" ~ 'shock')) %>% 
  mutate(type = "contemporary") %>% 
  select(-c(gear2))

#### historical fit #### 
hist_dat<-read.csv("Data/historical_lmb2_with_secchi.csv") %>% 
  drop_na(secchi_combo) #only use samples with secchi from contemp data 

#set up data 
hist_dat$logeffort<-log(hist_dat$effort_sum)
## standardize and transform predictor values ##
hist_dat$z_lake_area<-as.numeric(scale(log(hist_dat$lake_area_m2))) #good 
hist_dat$z_max_depth<-as.numeric(scale(log(hist_dat$max_depth_m))) #good 
hist_dat$z_dd_year-as.numeric(scale(log(hist_dat$dd_year))) #good
hist_dat$z_surface_temp_year<-as.numeric(scale(log(hist_dat$surface_temp_year))) #good 
hist_dat$z_secchi<-as.numeric(scale(log(hist_dat$secchi_combo)))
hist_dat$z_bottom_do<-as.numeric(scale(log(hist_dat$bottom_do + 0.001))) ##has 0s so added 0.001 
hist_dat$z_julian<-as.numeric(scale(hist_dat$julian)) ##standardize julian date 
hist_dat$z_ws_forest<-as.numeric(scale(asin(sqrt(hist_dat$ws_forest_prop)))) #good
hist_dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(hist_dat$ws_wetland_prop)))) #good

# need group number to match initial data 
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
  mutate(
    gear_index = case_when( gear == "FT_net" ~ 1,
                            gear == "gill" ~ 2,
                            gear == "seine" ~ 3)
  )

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
hist_dat$IND<-ifelse(hist_dat$gear == "FT_net", 0, 1) #FT_NET is the reference gear

#pull out about 50 lakes to match contemporary test dataset 
hist1<-partition(hist_dat$new_key, p=c(train=0.4,test=0.6), seed = 1, type =c("grouped"))
hist_test1<-hist_dat[hist1$test,]

hist_test<-hist_test1# just change this input each time 
n_distinct(hist_test$new_key)

# Container for predicted values
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
predict_hist<- array(NA, c(nsim,length(hist_test[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predict_hist)

#use regional mu.alphas instead of site specific, since dif sites (this takes like 5 mins)
for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(hist_test[[1]])){ #loop over covariate data (obs)
    predict_hist[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',hist_test[['group']][t],']')] + coefs[ID[i],paste0('b1[',hist_test[['gear_index']][t],']')]*hist_test[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',hist_test[['gear_index']][t],']')]*hist_test[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',hist_test[['gear_index']][t],']')]*hist_test[['z_surface_temp_year']][t]  + 
                               coefs[ID[i],paste0('b4[',hist_test[['gear_index']][t],']')]*hist_test[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',hist_test[['gear_index']][t],']')]*hist_test[['z_ws_forest']][t] + coefs[ID[i],paste0('b6[',hist_test[['gear_index']][t],']')]*hist_test[['z_ws_wetland']][t] +
                               coefs[ID[i],paste0('b7[',hist_test[['gear_index']][t],']')]*hist_test[['z_julian']][t] + 
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
hist_catch<-dplyr::select(hist_test, new_key, site, gear, largemouthbass_sum) %>% 
  mutate(row = row_number())
hist_plot_data<-left_join(med_data, hist_catch) %>% 
  mutate(gear = case_when(       #rename gear types for publication
    gear == "FT_net" ~ 'fyke net', 
    gear == "gill" ~ 'gillnet', 
    gear == "seine" ~ 'seine')) %>% 
  rename(col.med=hist.med, col.means=hist.means, fish_count_new=largemouthbass_sum ) %>% 
  mutate(type = "historical")


#* plot #### 
#combine cont and hist data to plot together 
both_data<-rbind(cont_plot_data, hist_plot_data) %>% 
  filter(gear != "fyke net" & gear != "shock")

#plot log scale 
pred_plot_log<-ggplot(data=both_data, aes(x=log(col.means+1), y=log(fish_count_new + 1), xmax=log(upperPI +1), xmin=log(lowerPI +1), colour = type, alpha=I(0.8)) )+
  geom_pointrange( ) + 
  scale_colour_manual(values = c("darkblue", "lightsalmon"), 
                      name= '') +
  xlab('predicted catch (log)')+
  ylab('observed catch (log)')+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  theme(legend.position = "bottom", legend.text=element_text(size=14)) +
  geom_abline(intercept = 0, slope = 1)


model2_hist_cont_pred_obs<-pred_plot_log+facet_wrap(~ gear, ncol=2, scales = 'free') #allow scales to vary 


ggsave(plot=model2_hist_cont_pred_obs, 
       device = "png", 
       filename = "figures/model2_hist_cont_pred_obs.png", 
       dpi = 600, height = 8, width = 12, units = "in",
       bg="#ffffff") #sets background to white 
