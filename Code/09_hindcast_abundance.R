#### Hindacst to same contemp lakes to compare contemporary and  historical densities ####
library(dplyr)
library(tidyr)
#### data ####
#surface temp average across years 
surface_temp<-read.csv("Data/MI_data/lake_surface_temp.csv") %>%
  group_by(IHDLKID) %>%
  slice_max(HECTARES) %>% 
  rename(nhdid = IHDLKID)%>% 
  mutate(temp_mean_contemp = mean(TAVE_2003:TAVE_2019), 
         temp_mean_hist = mean(TAVE_1936:TAVE_1964)) %>% 
  select(nhdid, temp_mean_contemp, temp_mean_hist)
  
  
#get temperature for each lake from the beginning of historical years (1936)
surface_temp_hist<-read.csv("Data/MI_data/lake_surface_temp.csv")%>%
  group_by(IHDLKID) %>%
  slice_max(HECTARES) %>% 
  rename(nhdid = IHDLKID) %>% 
  mutate(hist_temp = TAVE_1936) %>% 
  select(nhdid, hist_temp)
  
#merge data 
dat<-read.csv("Data/lmb_dat_for_model_aug30.csv") %>%
  left_join(surface_temp) %>% 
  mutate(z_lake_area = as.numeric(scale(log(lake_area_m2))), #transform/standardize
         z_max_depth=as.numeric(scale(log(maxdepth_m))),  
         z_temp_mean_contemp=as.numeric(scale(log(temp_mean_contemp))),  
         z_temp_mean_hist=as.numeric(scale(log(temp_mean_hist)))  ,
         z_secchi=as.numeric(scale(log(secchi_m))),
         z_doy=as.numeric(scale(day_of_year)) ,
         z_ws_forest=as.numeric(scale(asin(sqrt(ws_forest_prop)))) ,
         z_ws_wetland=as.numeric(scale(asin(sqrt(ws_wetland_prop)))) ,
         logeffort = log(effort_new)) %>% 
  select(new_key, fish_count_new, logeffort, gear2, FMU_Code, 
         z_lake_area, z_max_depth, z_secchi, z_doy, 
         z_temp_mean_contemp, 
         z_ws_forest, z_ws_wetland, z_temp_mean_hist, temp_mean_contemp, temp_mean_hist) %>%
  na.omit() %>% 
  mutate(group = as.numeric(as.factor(FMU_Code)),  #group
         gear = as.numeric(as.factor(gear2))  # gear 
         )


#read in model output 
output<-readRDS("Data/output/output_model2_secchi_fold5_lakes.rds")

#*density estimate for contemp and hindcast to contemp lakes ####
#index of density that is the sum of expected catches of all gears, using a common effort (E=1)
#extract just lake data (need a single row for each lake) 
just_lake<-dat %>% 
  distinct(new_key, .keep_all = TRUE) %>% #255 lakes 
  select(-c(gear2, logeffort, gear)) %>% #remove gear and effort info - will use common effort and all gear
  mutate(effort = 1) 

dup<-just_lake[rep(c(1:255),4),] # duplicate all rows 4x

dup_just_lake<-dup %>% 
  group_by(new_key) %>% 
  mutate(gear = case_when(row_number() %% 4 ==1 ~ "seine",
                          row_number() %% 4 ==2 ~ "gill",
                          row_number() %% 4 ==3 ~ "shock", 
                          TRUE ~ "fyke"
                           ) ) %>% 
  ungroup %>% 
  mutate(gear_index = case_when(gear == "fyke" ~ 1, 
                                gear == "gill" ~ 2,
                                gear == "seine" ~ 3,
                                gear == "shock" ~ 4,)) %>%
  mutate(IND = ifelse(gear == "fyke", 0, 1), #create an indicator variable
         doy_median = -0.1325897) #set DOY to the z_median 
  
#### set up simulation parameters #### 
coefs <- output$BUGSoutput$sims.matrix[,1:256] 
nsim <- 3000
chainLength <- output$BUGSoutput$n.sims
ID = seq( 1 , chainLength , floor(chainLength/nsim) )

#*current density estimates #### 
predictions_sim<- array(NA, c(nsim,length(dup_just_lake[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predictions_sim)

for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(dup_just_lake[[1]])){ #loop over covariate data (obs) #use regional mu.alphas instead of site specific
    predictions_sim[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',dup_just_lake[['group']][t],']')] + coefs[ID[i], paste0('b1[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_temp_mean_contemp']][t]  + 
                                  coefs[ID[i],paste0('b4[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_ws_forest']][t] + coefs[ID[i],paste0('b6[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_ws_wetland']][t] + coefs[ID[i],paste0('b7[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['doy_median']][t] +
                                  coefs[ID[i], paste0('logq[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['IND']][t] + dup_just_lake[['effort']][t]) # exp the predictions because we used a log-link in the negbi
    
  }
}

#From this distribution, you can extract the median/mean and the 2.5% and 97.5% percentiles (95% credible interval),
med <- apply(predictions_sim, 2, median )
means <- apply(predictions_sim, 2, mean )

#prep data for plotting
abund_data=data.frame(med, means) 
abund_data$row <- as.numeric(row.names(abund_data))
abund_data$log_mean<-log(abund_data$means)

dup_just_lake$site <- as.numeric(as.factor(dup_just_lake$new_key))
lake_ids<-dplyr::select(dup_just_lake, site, new_key) %>% 
  mutate(row = row_number())
cont_abund_data<-left_join(abund_data, lake_ids)

#add up the densities from all gears 
cont_sum_dens<-aggregate(cont_abund_data$means, by=list(new_key=cont_abund_data$new_key), FUN=sum) %>% 
  rename(sum_dens = x)

#*hindcast density estimates #### 
#*#same as above just change the temperature to the historical temp data 
hindcast<- array(NA, c(nsim,length(dup_just_lake[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(hindcast)

for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(dup_just_lake[[1]])){ #loop over covariate data (obs) #use regional mu.alphas instead of site specific
    hindcast[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',dup_just_lake[['group']][t],']')] + coefs[ID[i], paste0('b1[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_temp_mean_hist']][t]  + 
                                  coefs[ID[i],paste0('b4[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_ws_forest']][t] + coefs[ID[i],paste0('b6[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_ws_wetland']][t] + coefs[ID[i],paste0('b7[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['doy_median']][t] +
                                  coefs[ID[i], paste0('logq[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['IND']][t] + dup_just_lake[['effort']][t]) # exp the predictions because we used a log-link in the negbi
    
  }
}

#From this distribution, you can extract the median/mean and the 2.5% and 97.5% percentiles (95% credible interval),
med <- apply(hindcast, 2, median )
means <- apply(hindcast, 2, mean )

#prep data for plotting
hind_abund_data=data.frame(med, means) 
hind_abund_data$row <- as.numeric(row.names(hind_abund_data))
hind_abund_data$log_mean<-log(hind_abund_data$means)

dup_just_lake$site <- as.numeric(as.factor(dup_just_lake$new_key))
lake_ids<-dplyr::select(dup_just_lake, site, new_key) %>% 
  mutate(row = row_number())
hindcast_data<-left_join(hind_abund_data, lake_ids)

#add up the densities from all gears 
hind_sum_dens<-aggregate(hindcast_data$means, by=list(new_key=hindcast_data$new_key), FUN=sum) %>% 
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

hist_abund_map_dat<-left_join(hind_sum_dens, points)
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
library(viridis)
hist_map<-p+ geom_point(data=hist_abund_map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(log(sum_dens))), size = 3) + 
  scale_color_viridis(direction = -1, limits = c(-1.2,4.2) )  + #
  labs(color="log relative density", tag = "a")  +#changes the labels on the legend 
  theme(legend.position = "bottom", 
        plot.tag = element_text(), 
        plot.tag.position = c(0.1,0.95)) + 
  ylab(NULL) + xlab(NULL)
hist_map

#map of lmb abund contemporary 
#log
cont_map<-p+ geom_point(data=cont_abund_map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(log(sum_dens))), size = 3) + 
  scale_color_viridis(direction = -1, limits = c(-1.2,4.2) )  + #
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

#### look at changes in relative abundance #### 
abund_est<-left_join(cont_abund_map_dat, hist_abund_map_dat, by = c("new_key")) %>% 
  rename(dens_cont= sum_dens.x, dens_hist = sum_dens.y) %>% 
  mutate(change = dens_cont - dens_hist, 
         color = ifelse(change>0, "increase", "decrease")) 

plot(abund_est$change)

p+ geom_point(data=abund_est, aes(x = LONG_DD.x, y = LAT_DD.x, colour = color)) + 
  scale_color_manual(values= c("increase" = "darkgreen", 
                               "decrease" =  'coral3') ) +
                labs(color="density change") +  #changes the labels on the legend
  theme(legend.position = "bottom") + 
  ylab(NULL) + xlab(NULL) + 
  theme_void()


p+ geom_point(data=abund_est, aes(x = LONG_DD.x, y = LAT_DD.x, colour = change), size = 3) + 
  scale_color_gradient(low = "coral3", high = "darkgreen")  + #
  labs(color="change in density") +  #changes the labels on the legend
  theme(legend.position = "bottom") + 
  ylab(NULL) + xlab(NULL)

plot(abund_est$change)
hist(abund_est$change)
plot(abund_est$LAT_DD.x,abund_est$change)


### HISTOGRAMS ####
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


ggplot(abund_histogram, aes(x=log(mean), fill=time)) +
  geom_histogram(alpha=0.5, position="identity", aes(y = ..density..), color="black", bins = 10) +
  labs(x="relative density (log)", y = "frequency") + 
  scale_fill_manual(values=c("lightblue", "lightsalmon")) +
  guides(fill=guide_legend(title='')) + 
  theme(legend.position = c(0.9,0.8), legend.text=element_text(size=14), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  theme_bw() 


#### changes in water temp of the lakes #### 
just_lake<-dat %>% 
  distinct(new_key, .keep_all = TRUE) %>% #255 lakes 
  mutate(temp_change = temp_mean_contemp - temp_mean_hist) %>% #pos means got warmer 
left_join(abund_est, by="new_key") #pos change means increase in abund

summary(just_lake$temp_change)

plot(just_lake$temp_change, just_lake$change, 
     xlab = "temperature change (degC)", ylab="projected change in relative density")
abline(v = 0, h = 0, lty=2)

hist(just_lake$change,  xlab ="projected change in relative density", main = "")

#*ggplot
ggplot(just_lake, aes(x=change)) +
  geom_histogram() +
  labs(x="estimated change in relative density", y = "frequency") + 
  theme_bw() 
 
ggplot(just_lake, aes(x=temp_change, y = change)) + 
  geom_point() + 
  scale_color_gradient(low = "coral3", high = "darkgreen")  +
  geom_hline(yintercept = 0, linetype = 'dashed')  + 
  geom_vline(xintercept = 0, linetype = 'dashed') +
  labs(x="temperature change (degC)", y = "projected change in relative density") +
  theme_bw()

ggsave(plot=fig5, 
       device = "png", 
       filename = "figures/abund_histogram_fig7.png", 
       dpi = 600, height = 8, width = 12, units = "in",
       bg="#ffffff") #sets background to white 

#map outliers 
#new key 12-21
new_1221<-filter(cont_abund_map_dat, new_key == "12-21")
p+ geom_point(data=new_1221, aes(x = LONG_DD, y = LAT_DD)) 

#new_key 41-464
new_41464<-filter(cont_abund_map_dat, new_key == "41-464")
p+ geom_point(data=new_41464, aes(x = LONG_DD, y = LAT_DD)) 

#new_key 14-154
new_14154<-filter(cont_abund_map_dat, new_key == "14-154")
p+ geom_point(data=new_14154, aes(x = LONG_DD, y = LAT_DD)) 
