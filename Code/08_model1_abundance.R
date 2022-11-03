#### estimate relative abundance using model 1 

#### estimate relative abundance historical ####

#*read in model output ####
output<-readRDS("Data/output/output_model1_secchi_fold5_lakes.rds")


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


#save predicted abundance in each lake 
write.csv(hist_abund_map_dat, "Data/output/model1_hist_abund_pred.csv", row.names = FALSE)
write.csv(cont_abund_map_dat, "Data/output/model1_cont_abund_pred.csv", row.names = FALSE)

model2_hist_est<-read.csv("Data/output/model2_hist_abund_pred.csv")
model2_cont_est<-read.csv("Data/output/model2_cont_abund_pred.csv")