### Map figures #### 
#devtools::install_github("dkahle/ggmap")
library(ggmap)
library(stringr)
library(devtools)
library(mapdata)
library(dplyr)
library(ggpubr)


#get the lat/long points 
lmb_points<-read.csv("Data/lmb_dat_for_model_aug30.csv") %>%
  dplyr::select( new_key, LAT_DD, LONG_DD, fish_count_new, dd_mean, surface_temp_mean) %>%
  distinct(new_key, .keep_all = TRUE)

# get historical data 
hist_dat<-read.csv("/Users/katelynking/Desktop/historical_lmb.csv")

points<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/IFR_Lake_Points.csv")%>%
  rename(new_key = NEW_KEY) %>% 
  select(new_key, LONG_DD, LAT_DD)

for_map<-hist_dat %>% 
  distinct(new_key, .keep_all = TRUE) %>% 
  left_join(points) %>% 
  select(new_key, LONG_DD, LAT_DD)%>% 
  mutate(data = 'historical')

#get MI basemap 
MI_basemap<-map_data("state",  region = c("michigan"))  # select michigan 
map<-ggplot(data = MI_basemap) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + #this fill MI white and outline is black
  coord_fixed(1.3) 


######### map of cont and hist data  ####
#is 201 lakes when you have to have environ variables that match 
lmb_model<-lmb_points%>%
            na.omit() %>% 
  select(new_key, LONG_DD, LAT_DD) %>% 
  mutate(data = 'contemporary')
#combine data for map 
map_data<-rbind(lmb_model, for_map)

map_of_data<-map + 
  geom_point(data=map_data, aes(x = LONG_DD, y = LAT_DD, color=c(data) )) +
  theme_bw() +  
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank()) +
  clean_theme() +   
  scale_colour_manual(values = c("darkblue", "lightsalmon"), 
                      name= '') +
  theme(legend.position = c(0.20, 0.25), 
        legend.text = element_text(size=12)) 

  

ggsave("/Users/katelynking/Desktop/map_of_data.png", map_of_data)


####map of lmb abund #### 
map+ geom_point(data=lmb_points, aes(x = LONG_DD, y = LAT_DD, colour = c(log(fish_count_new + 1)))) + 
                scale_colour_gradient(low = "yellow", high = "darkblue")  + 
                labs(color="lmb log abund")  #changes the labels on the legend 

#map of temps  
p+ geom_point(data=lmb_points, aes(x = LONG_DD, y = LAT_DD, colour = c(dd_mean))) + 
  scale_colour_gradient(low = "blue", high = "red")  + 
  labs(color="gdd")  #changes the labels on the legend 

