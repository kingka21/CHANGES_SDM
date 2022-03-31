### Map figures #### 
#devtools::install_github("dkahle/ggmap")
library(ggmap)
library(stringr)
library(devtools)
library(mapdata)
library(dplyr)

#get the lat/long points 
lmb_points<-read.csv("Data/lmb_dat_for_model_mar14.csv") %>%
  dplyr::select( new_key, LAT_DD, LONG_DD, fish_count_new, dd_mean, surface_temp_mean) %>%
  distinct(new_key, .keep_all = TRUE)

#get MI basemap 
MI_basemap<-map_data("state") %>%
        subset(region %in% c("michigan")) # select michigan 
p<-ggplot(data = MI_basemap) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + #this fill MI white and outline is black
  coord_fixed(1.3) 
#map of lmb abund 
p+ geom_point(data=lmb_points, aes(x = LONG_DD, y = LAT_DD, colour = c(log(fish_count_new + 1)))) + 
                scale_colour_gradient(low = "yellow", high = "darkblue")  + 
                labs(color="lmb log abund")  #changes the labels on the legend 

#map of temps  
p+ geom_point(data=lmb_points, aes(x = LONG_DD, y = LAT_DD, colour = c(dd_mean))) + 
  scale_colour_gradient(low = "blue", high = "red")  + 
  labs(color="gdd")  #changes the labels on the legend 

