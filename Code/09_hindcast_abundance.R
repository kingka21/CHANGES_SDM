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

################################################ 
#### take the mean of differences #### 
########################################################
#get difference first, then average the difference 
predictions_sim<- array(NA, c(nsim,length(dup_just_lake[[1]]))) # 2 dimensions - length of data as columns and sims as rows 
dim(predictions_sim)
hindcast<- array(NA, c(nsim,length(dup_just_lake[[1]]))) 
dim(hindcast)
changes <- array(NA, c(nsim,length(dup_just_lake[[1]])))
dim(changes)

for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(dup_just_lake[[1]])){ #loop over covariate data (obs) #use regional mu.alphas instead of site specific
    predictions_sim[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',dup_just_lake[['group']][t],']')] + coefs[ID[i], paste0('b1[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_temp_mean_contemp']][t]  + 
                                  coefs[ID[i],paste0('b4[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_ws_forest']][t] + coefs[ID[i],paste0('b6[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_ws_wetland']][t] + coefs[ID[i],paste0('b7[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['doy_median']][t] +
                                  coefs[ID[i], paste0('logq[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['IND']][t] + dup_just_lake[['effort']][t]) # exp the predictions because we used a log-link in the negbi
    
    hindcast[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',dup_just_lake[['group']][t],']')] + coefs[ID[i], paste0('b1[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_temp_mean_hist']][t]  + 
                           coefs[ID[i],paste0('b4[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_ws_forest']][t] + coefs[ID[i],paste0('b6[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['z_ws_wetland']][t] + coefs[ID[i],paste0('b7[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['doy_median']][t] +
                           coefs[ID[i], paste0('logq[',dup_just_lake[['gear_index']][t],']')]*dup_just_lake[['IND']][t] + dup_just_lake[['effort']][t]) # exp the predictions because we used a log-link in the negbi
    changes[i,t] <- predictions_sim[i,t] - hindcast[i,t]

  }
}


#From this distribution, you can extract the median/mean and the 2.5% and 97.5% percentiles (95% credible interval),
med <- apply(changes, 2, median ) #2 manipulation is performed on columns 
means <- apply(changes, 2, mean )

#prep data for plotting
abund_data=data.frame(med, means) 
abund_data$row <- as.numeric(row.names(abund_data))

dup_just_lake$site <- as.numeric(as.factor(dup_just_lake$new_key))
lake_ids<-dplyr::select(dup_just_lake, site, new_key) %>% 
  mutate(row = row_number())
change_abund_data<-left_join(abund_data, lake_ids)

#add up the change in densities from all gears 
sum_change_dens<-aggregate(change_abund_data$means, by=list(new_key=change_abund_data$new_key), FUN=sum) %>% 
  rename(dens_change = x)

#*get means from contemporary 
#From this distribution, you can extract the median/mean and the 2.5% and 97.5% percentiles (95% credible interval),
med <- apply(predictions_sim, 2, median ) #2 manipulation is performed on columns 
mean_cont_dens <- apply(predictions_sim, 2, mean )

#prep data for plotting
dens_contemp=data.frame(mean_cont_dens) 
dens_contemp$row <- as.numeric(row.names(dens_contemp))

dup_just_lake$site <- as.numeric(as.factor(dup_just_lake$new_key))
lake_ids<-dplyr::select(dup_just_lake, site, new_key) %>% 
  mutate(row = row_number())
cont_abund_data<-left_join(dens_contemp, lake_ids)

#add dens across gears
cont_sum_dens<-aggregate(cont_abund_data$mean_cont_dens, by=list(new_key=cont_abund_data$new_key), FUN=sum) %>% 
  rename(sum_dens_cont = x)

#* get hindcast density estimates
med <- apply(hindcast, 2, median )
mean_hist_dens <- apply(hindcast, 2, mean )

#prep data for plotting
hind_abund_data=data.frame(mean_hist_dens) 
hind_abund_data$row <- as.numeric(row.names(hind_abund_data))

hindcast_data<-left_join(hind_abund_data, lake_ids)

#add up the densities from all gears 
hind_sum_dens<-aggregate(hindcast_data$mean_hist_dens, by=list(new_key=hindcast_data$new_key), FUN=sum) %>% 
  rename(sum_dens_hist = x)

#* join everything together and add temp change 
just_lake<-dup_just_lake %>% 
  distinct(new_key, .keep_all = TRUE) %>% #255 lakes 
  mutate(z_temp_change = z_temp_mean_contemp - z_temp_mean_hist) %>% #pos means got warmer 
  left_join(sum_change_dens, by="new_key") %>% #pos change means increase in abund
  left_join(hind_sum_dens) %>% 
  left_join(cont_sum_dens)
  
#save data 
#write.csv(just_lake, "Data/output/density_temp_changes.csv", row.names=FALSE)

##### map and graph of predicted abundance Fig 5 ####
library(ggmap)
library(stringr)
library(devtools)
library(mapdata)
library(dplyr)

#get the lat/long points - need lat/lon from Humphries table 
points<-read.csv("Data/MI_data/Humphries_table.csv") %>%
  dplyr::select( New_Key, LAT_DD, LONG_DD) %>%
  rename(new_key = New_Key) 

map_dat<-left_join(just_lake, points) 
summary(log(map_dat$sum_dens_hist))
summary(log(map_dat$sum_dens_cont ))
summary(map_dat$sum_dens_hist)
summary(map_dat$sum_dens_cont )
summary(map_dat$dens_change )
summary(map_dat$z_temp_change)

#get MI basemap 
MI_basemap<-map_data("state") %>%
  subset(region %in% c("michigan")) # select michigan 
p<-ggplot(data = MI_basemap) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + #this fill MI white and outline is black
  coord_fixed(1.3) 

#map of lmb abund historical 
#log
library(viridis)
hist_map<-p+ geom_point(data=map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(log(sum_dens_hist))), size = 3) + 
  scale_color_viridis(direction = -1, limits = c(-1.2,4.2) )  + #
  labs(color="log relative density", tag = "a")  +#changes the labels on the legend 
  theme_bw() + 
  theme( legend.position = c(0.25,0.25),
        plot.tag = element_text(), 
        plot.tag.position = c(0.1,0.95)) + 
  ylab(NULL) + xlab(NULL)
hist_map

#map of lmb abund contemporary 
#log
cont_map<-p+ geom_point(data=map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(log(sum_dens_cont))), size = 3) + 
  scale_color_viridis(direction = -1, limits = c(-1.2,4.2) )  + #
  labs(color="log relative density", tag = "b") +  #changes the labels on the legend
  theme(legend.position = "bottom", 
        plot.tag = element_text(), 
        plot.tag.position = c(0.1,0.95)) + 
  ylab(NULL) + xlab(NULL)

cont_map

#map of change in density from hist to contemp 
den_change_map<-p+ geom_point(data=map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(dens_change)), size = 3) +  #
  scale_color_gradientn(
    colours=colorRampPalette((RColorBrewer::brewer.pal(11,"PRGn")))(255),
    values = c(1.0, (0-min(map_dat$dens_change))/(max(map_dat$dens_change)-min(map_dat$dens_change)),0)
  ) + 
  labs(color="density change", tag = "b") + #changes the labels on the legend 
  theme_bw() + 
  theme(legend.position = c(0.25,0.25),
        plot.tag = element_text(), 
        plot.tag.position = c(0.1,0.95)) + 
  ylab(NULL) + xlab(NULL) 

temp_change_map<-p+ geom_point(data=map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(z_temp_change)), size = 3) +  #
  scale_color_gradientn(
    colours=colorRampPalette((RColorBrewer::brewer.pal(11,"RdBu")))(255),
    values = c(1.0, (0-min(map_dat$z_temp_change))/(max(map_dat$z_temp_change)-min(map_dat$z_temp_change)),0)
  ) + 
  labs(color="temperature change", tag = "c") + #changes the labels on the legend 
  theme_bw() + 
  theme( legend.position = c(0.25,0.25),
        plot.tag = element_text(), 
        plot.tag.position = c(0.1,0.95)) + 
  ylab(NULL) + xlab(NULL) 

temp_dens_plot<-ggplot(map_dat, aes(x=z_temp_change, y = dens_change, color = z_lake_area )) + 
  geom_point() + 
  scale_color_viridis(direction = -1 ) +
  #scale_color_gradient(low = "coral3", high = "darkgreen")  +
  geom_hline(yintercept = 0, linetype = 'dashed')  + 
  geom_vline(xintercept = 0, linetype = 'dashed') +
  labs(x="temperature change (degC)", y = "change in estimated relative density", color = "lake area (scaled)", tag = "d") +
  theme_bw() + 
  theme( legend.position = c(0.80,0.25),
        plot.tag = element_text(), 
        plot.tag.position = c(0.13,0.95))

#to to align plots 
right_plots<-cowplot::plot_grid(den_change_map, temp_dens_plot, ncol=1, align="v", axis = "l")
left_plots<-cowplot::plot_grid(hist_map, temp_change_map, ncol=1)
figure5<-cowplot::plot_grid(left_plots, right_plots)

ggsave(plot=figure5, 
       device = "png", 
       filename = "figures/figure5.png", 
       dpi = 600, height = 8, width = 8, units = "in",
       bg="#ffffff") #sets background to white 


### HISTOGRAM

ggplot(just_lake, aes(x=dens_change)) +
  geom_histogram() +
  labs(x="estimated change in relative density", y = "frequency") + 
  theme_bw() 
 
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
