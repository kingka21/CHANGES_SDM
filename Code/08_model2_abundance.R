#### estimation of relative abundance
#use model 2 - better model 

library(dplyr)
library(ggplot2)
library(cowplot)

#### relative abundance historical ####
#*read in model output and data ####
output<-readRDS("Data/output/output_model2_secchi_fold5_lakes.rds")

#set up simulation parameters 
coefs <- output$BUGSoutput$sims.matrix[,1:256] #obtain our samples from the posterior
nsim <- 3000 # Number of desired MCMC samples used for prediction (use a subset of all samples)
chainLength <- output$BUGSoutput$n.sims # Chain length from analysis
ID = seq( 1 , chainLength , floor(chainLength/nsim) ) #take values from length of posterior

#*read in data 
hist_dat<-read.csv("Data/historical_lmb2_with_secchi.csv") %>% 
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
  ) %>% 
  mutate(IND = ifelse(gear == "FT_net", 0, 1)) #create an indicator variable “IND”, with value 0 for reference gear and value 1 otherwise

#pull out about 50 lakes to match contemporary test dataset 
hist1<-partition(hist_dat$new_key, p=c(train=0.4,test=0.6), seed = 1, type =c("grouped"))
hist_test1<-hist_dat[hist1$test,]

#*abundance 
#index of density that is the sum of expected catches of all gears, using a common effort (E=1)
#extract just lake data (need a single row for each lake) 
hist_just_lake<-hist_test1 %>% 
  distinct(new_key, .keep_all = TRUE) %>% #50 lakes 
  select(-c(gear, effort_sum, logeffort, gear_index)) %>% #remove gear and effort info - will use common effort and both gear
  mutate(effort = 1) 

dup<-hist_just_lake[rep(c(1:50),2),] # duplicate all rows 

hist_just_lake<-dup %>% 
  group_by(new_key) %>% 
  mutate(gear = if_else(row_number() %% 2 ==1, "seine", "gill")) %>% #add just seine and gill 
  ungroup %>% 
  mutate(gear_index = ifelse(gear == "gill", 2, 3))

#* predict density histrocial ####
predict_hist<- array(NA, c(nsim,length(hist_just_lake[[1]]))) 
dim(predict_hist)

#use regional mu.alphas instead of site specific, since dif sites, use estimated q and effort = 1
for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(hist_just_lake[[1]])){ #loop over covariate data (obs)
    predict_hist[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',hist_just_lake[['group']][t],']')] + coefs[ID[i],paste0('b1[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_surface_temp_year']][t]  + 
                               coefs[ID[i],paste0('b4[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_ws_forest']][t] + coefs[ID[i],paste0('b6[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_ws_wetland']][t] +
                               coefs[ID[i],paste0('b7[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_julian']][t] + 
                               coefs[ID[i], paste0('logq[',hist_just_lake[['gear_index']][t],']')] + hist_just_lake[['effort']][t]) # exp the predictions because we used a log-link in the Poisson
    
  }
}

#From this distribution, you can extract the median/mean and the 2.5% and 97.5% percentiles (95% credible interval),
hist.med <- apply(predict_hist, 2, median )
hist.means <- apply(predict_hist, 2, mean )

#prep data for plotting
abund_data=data.frame(hist.med, hist.means) 
abund_data$row <- as.numeric(row.names(abund_data))
abund_data$log_mean<-log(abund_data$hist.means)

hist_just_lake$site <- as.numeric(as.factor(hist_just_lake$new_key))
hist_ids<-dplyr::select(hist_just_lake, site, new_key) %>% 
  mutate(row = row_number())
map_data<-left_join(abund_data, hist_ids)

#add up the densities from both gears 
hist_sum_dens<-aggregate(map_data$hist.means, by=list(new_key=map_data$new_key), FUN=sum) %>% 
  rename(sum_dens = x)

#### relative abundance contemporary ####
#*read in data ####
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

dat<-left_join(dat, surface_temp, by = c("nhdid", "year"))

## standardize and transform predictor values ##
dat$z_lake_area<-as.numeric(scale(log(dat$lake_area_m2))) #good 
dat$z_max_depth<-as.numeric(scale(log(dat$maxdepth_m))) #good 
dat$z_surface_temp_year<-as.numeric(scale(log(dat$surf_temp_year))) #good 
dat$z_secchi<-as.numeric(scale(log(dat$secchi_m)))
dat$z_doy<-as.numeric(scale(dat$day_of_year)) ##standardize julian date 
dat$z_ws_forest<-as.numeric(scale(asin(sqrt(dat$ws_forest_prop)))) #good
dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat$ws_wetland_prop)))) #good

dat$logeffort <- log(dat$effort_new)

dat<-dplyr::select(dat, new_key, fish_count_new, logeffort, gear2, FMU_Code, 
                   z_lake_area, z_max_depth, z_secchi, z_doy, 
                   z_surface_temp_year,
                   z_ws_forest, z_ws_wetland) %>%
  na.omit()

part5<-partition(dat$new_key, p=c(train=0.8, test = 0.2), seed = 5, type =c( "grouped"))
train5<-dat[part5$train,]
test5<-dat[part5$test,]

#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
cont_test<-test5 %>% 
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
  mutate(IND = ifelse(gear2 == "FT_NET", 0, 1)) #create an indicator variable “IND”

#extract just lake data (need a single row for each lake) 
cont_just_lake<-cont_test %>% 
  distinct(new_key, .keep_all = TRUE) %>% #51 lakes 
  select(-c(gear2, gear, logeffort, IND)) %>% 
  mutate(effort = 1) 

dup<-cont_just_lake[rep(c(1:51),2),] # duplicate all rows 

cont_just_lake<-dup %>% 
  group_by(new_key) %>% 
  mutate(gear = if_else(row_number() %% 2 ==1, "seine", "gill")) %>% #add just seine and gill 
  ungroup %>% 
  mutate(gear_index = ifelse(gear == "gill", 2, 3))


#* predict contemporary ####
predict_cont<- array(NA, c(nsim,length(cont_just_lake[[1]]))) 
dim(predict_cont)

#use lake/site specific alphas and environment, no gear data 
for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(cont_just_lake[[1]])){ #loop over covariate data (obs)
    predict_cont[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',cont_just_lake[['group']][t],']')] + coefs[ID[i],paste0('b1[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_surface_temp_year']][t]  + 
                               coefs[ID[i],paste0('b4[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_ws_forest']][t] + coefs[ID[i],paste0('b6[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_ws_wetland']][t] +
                               coefs[ID[i],paste0('b7[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_doy']][t] + 
                               coefs[ID[i], paste0('logq[',cont_just_lake[['gear_index']][t],']')] + cont_just_lake[['effort']][t]) # exp the predictions because we used a log-link in the Poisson
    
  }
}


#From this distribution, you can extract the median and the 2.5% and 97.5% percentiles (95% credible interval),
cont.med.abund <- apply(predict_cont, 2, median )
cont.means.abund <- apply(predict_cont, 2, mean )

#prep data for plotting
cont_abund_data=data.frame(cont.med.abund, cont.means.abund) 
cont_abund_data$row <- as.numeric(row.names(cont_abund_data))
cont_abund_data$log_mean<-log(cont_abund_data$cont.means.abund)

cont_just_lake$site <- as.numeric(as.factor(cont_just_lake$new_key))
cont_ids<-dplyr::select(cont_just_lake, site, new_key) %>% 
  mutate(row = row_number())
cont_map_data<-left_join(cont_abund_data, cont_ids) 

#add up the densities from both gears 
cont_sum_dens<-aggregate(cont_map_data$cont.means.abund, by=list(new_key=cont_map_data$new_key), FUN=sum) %>% 
  rename(sum_dens = x)


##### map predicted abundance ####
library(ggmap)
library(stringr)
library(devtools)
library(mapdata)
library(dplyr)

#get the lat/long points - need lat/lon from Humphries table 
points<-read.csv("Data/MI_data/Humphries_table.csv") %>%
  dplyr::select( New_Key, LAT_DD, LONG_DD) %>%
  rename(new_key = New_Key) 

hist_abund_map_dat<-left_join(hist_sum_dens, points)
cont_abund_map_dat<-left_join(cont_sum_dens, points)

summary(log(hist_abund_map_dat$sum_dens))
summary(log(cont_abund_map_dat$sum_dens ))
summary(hist_abund_map_dat$sum_dens)
summary(cont_abund_map_dat$sum_dens)

#get MI basemap 
MI_basemap<-map_data("state") %>%
  subset(region %in% c("michigan")) # select michigan 
p<-ggplot(data = MI_basemap) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + #this fill MI white and outline is black
  coord_fixed(1.3) 
#map of lmb abund historical 
#log
hist_map<-p+ geom_point(data=hist_abund_map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(log(sum_dens))), size = 3) + 
  scale_color_viridis(direction = -1, limits = c(-3.5,3) )  + #
  labs(color="log relative density", tag = "a")  +#changes the labels on the legend 
  theme(legend.position = "bottom", 
        plot.tag = element_text(), 
        plot.tag.position = c(0.1,0.95)) + 
  ylab(NULL) + xlab(NULL)
hist_map

#map of lmb abund contemporary 
#log
cont_map<-p+ geom_point(data=cont_abund_map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(log(sum_dens))), size = 3) + 
  scale_color_viridis(direction = -1, limits = c(-3.5,3) )  + #
  labs(color="log relative density", tag = "b") +  #changes the labels on the legend
  theme(legend.position = "bottom", 
        plot.tag = element_text(), 
        plot.tag.position = c(0.1,0.95)) + 
  ylab(NULL) + xlab(NULL)

cont_map
#* save as 2-panel plot 
map_model2<-cowplot::plot_grid(hist_map, cont_map)
map_model2

ggsave(plot=map_model2, 
       device = "png", 
       filename = "figures/model2_map_abund.png", 
       dpi = 600, height = 8, width = 12, units = "in",
       bg="#ffffff") #sets background to white 

#save predicted abundance in each lake 
#write.csv(hist_abund_map_dat, "Data/output/model2_hist_abund_pred.csv", row.names = FALSE)
#write.csv(cont_abund_map_dat, "Data/output/model2_cont_abund_pred.csv", row.names = FALSE)

### HISTOGRAMS ####

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

hist_1<-hist(cont_abund_map_dat$sum_dens)
hist_2<-hist(hist_abund_map_dat$sum_dens)
plot(hist_1, col = c1, xlim=c(0,10), xlab="relative density")
plot(hist_2, col = c2, add=T)

hist_3<-hist(log(cont_abund_map_dat$sum_dens))
hist_4<-hist(log(hist_abund_map_dat$sum_dens), breaks=5)
plot(hist_3, col = c1, xlab="relative density (log)")
plot(hist_4, col = c2, add=T)

#ggplot 
cont_histogram<-cont_abund_map_dat %>% 
  rename(mean= sum_dens) %>% 
  select(new_key, mean) %>% 
  mutate(time = "contemporary")

hist_histogram<-hist_abund_map_dat %>% 
  rename(mean= sum_dens) %>% 
  select(new_key, mean) %>% 
  mutate(time = "historical")

abund_histogram<-rbind(cont_histogram, hist_histogram)


abund_histogram_fig7<-ggplot(abund_histogram, aes(x=log(mean), fill=time)) +
  geom_histogram(alpha=0.25, position="identity", aes(y = ..density..), color="black", bins = 10) +
  labs(x="relative density (log)", y = "frequency") + 
  scale_fill_manual(values=c("lightblue", "lightsalmon")) +
  guides(fill=guide_legend(title='')) + 
  theme(legend.position = c(0.9,0.8), legend.text=element_text(size=14))

ggsave(plot=abund_histogram_fig7, 
       device = "png", 
       filename = "figures/abund_histogram_fig7.png", 
       dpi = 600, height = 8, width = 12, units = "in",
       bg="#ffffff") #sets background to white 

### GAM latitude plot  #### 
cont_dat<-read.csv("Data/output/model2_cont_abund_pred.csv") %>% 
  mutate(type='contemporary')
hist_dat<-read.csv("Data/output/model2_hist_abund_pred.csv") %>% 
  mutate(type='historical')

library(mgcv)


cont_dat<-cont_just_lake %>% 
  distinct(new_key, .keep_all = TRUE) %>% 
  left_join(cont_abund_map_dat) %>% 
  select(new_key, z_surface_temp_year, sum_dens, LAT_DD, LONG_DD) %>% 
  mutate(type='contemporary')
hist_dat<-hist_just_lake %>% 
  distinct(new_key, .keep_all = TRUE)%>% 
  left_join(hist_abund_map_dat) %>% 
  select(new_key, z_surface_temp_year, sum_dens, LAT_DD, LONG_DD)%>% 
  mutate(type='historical')

cont_gam<-mgcv::gam(log(sum_dens) ~ s(LAT_DD, bs="cs"), data = cont_dat)
summary(cont_gam)

#pull out data
cont_plot_dat<-plot(cont_gam)

#turn into dataframe
dat2<-data.frame(x = cont_plot_dat[[1]]$x,
                 y = cont_plot_dat[[1]]$fit,
                 ymin = cont_plot_dat[[1]]$fit - cont_plot_dat[[1]]$se,
                 ymax = cont_plot_dat[[1]]$fit + cont_plot_dat[[1]]$se)

dat2_obs <- data.frame(x = cont_plot_dat[[1]]$raw)

# plot cont gam
p2<-ggplot(dat2, aes(x = x, y = y)) +
  geom_line() +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = 'gray', alpha = 0.3) +
  labs(x = "latitude", y = "density (log)") +
  geom_point(data = dat1_obs, aes(x = x, y = -Inf), shape = '|', size = 6, alpha = 0.10)
p2


hist_gam<-mgcv::gam(log(sum_dens) ~ s(LAT_DD), data = hist_dat)
summary(hist_gam)

AIC(cont_gam, hist_gam)
anova(cont_gam, hist_gam)

#pull out data
hist_plot_dat<-plot(hist_gam)

#turn into dataframe
dat1<-data.frame(x = hist_plot_dat[[1]]$x,
                 y = hist_plot_dat[[1]]$fit,
                 ymin = hist_plot_dat[[1]]$fit - hist_plot_dat[[1]]$se,
                 ymax = hist_plot_dat[[1]]$fit + hist_plot_dat[[1]]$se)

dat1_obs <- data.frame(x = hist_plot_dat[[1]]$raw)

# plot gam1
p1<-ggplot(dat1, aes(x = x, y = y)) +
  geom_line() +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), color = 'gray', alpha = 0.3) +
  labs(x = "latitude", y = "density (log)") +
  geom_point(data = dat1_obs, aes(x = x, y = -Inf), shape = '|', size = 6, alpha = 0.10)
p1



#other plot where the two datasets are combined 
plot_dat<-rbind(cont_dat, hist_dat)
gam1 <- gam(log(sum_dens) ~ type + s(LAT_DD), data = plot_dat, method = "ML")
summary(gam1) #looks like the time period is not really different 

gam_plot<-ggplot(data=plot_dat, aes(x=LAT_DD, y=log(sum_dens), color=type)) + 
  geom_point()  +
  scale_color_manual(values=c("darkblue", "lightsalmon")) + 
  geom_smooth(method='gam') +
  guides(color=guide_legend(title='')) + 
  xlab("latitude") + ylab("density (log)") + 
  theme(legend.position = c(0.9,0.8), legend.text=element_text(size=14), 
        axis.title = element_text(size=16), axis.text = element_text(size=14))
gam_plot
ggsave("figures/gam_plot.png", gam_plot)
