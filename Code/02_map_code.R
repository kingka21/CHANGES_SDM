### Map figures #### 
#devtools::install_github("dkahle/ggmap")
library(ggmap)
library(stringr)
library(devtools)
library(mapdata)
library(dplyr)
library(ggpubr)
library(sf)
library(ggplot2)
library(cowplot)
library(spData)

#### For Figure 1 sampling locations #### 
#get the lat/long points 
cont_points<-read.csv("Data/contemp_lmb_dat.csv") %>%
  dplyr::select( new_key, LAT_DD, LONG_DD, fish_count_new, lake_area_m2, maxdepth_m, surface_temp_mean, secchi_m, day_of_year, ws_forest_prop, ws_wetland_prop) %>%
  distinct(new_key, .keep_all = TRUE) %>% 
  na.omit()%>% 
  mutate(data = 'contemporary') %>%
  dplyr::select(new_key, LONG_DD, LAT_DD, data)

# get historical data 
hist_dat<-read.csv("Data/historical_lmb2_with_secchi.csv") %>%
  na.omit() %>% 
  distinct(new_key, .keep_all = TRUE)

points<-read.csv("Data/Humphries_table.csv")%>%
  rename(new_key = New_Key) %>% 
  dplyr::select(new_key, LONG_DD, LAT_DD)

hist_points<-hist_dat  %>% 
  left_join(points) %>% 
  dplyr::select(new_key, LONG_DD, LAT_DD)%>% 
  mutate(data = 'historical')



#shapefile with the FMU boundaries 
library(sf)
fmu_shape<- st_read(dsn = "Data/MDNR_FMUs", layer = "Michigan_DNR_Fisheries_Management_Units") 
raster::crs(fmu_shape)

#combine data for map 
map_data<-rbind(cont_points, hist_points)

map_of_data<-ggplot() +
  geom_sf(data=fmu_shape, aes(), fill = "white") + 
  geom_point(data=map_data, aes(x = LONG_DD, y = LAT_DD, color=c(data) )) +
  theme_bw() +  
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank()) +
  #clean_theme() + 
  theme_void() +
  scale_colour_manual(values = c("darkblue", "lightsalmon"), 
                      name= '') +
  theme(legend.position = c(0.20, 0.25), 
        legend.text = element_text(size=12)) 

map_of_data


# get US inset 
data("us_states", package = "spData")
us_states = st_transform(us_states, crs = "+proj=longlat +datum=WGS84 +no_defs") #match the FMU shapefile
MI_box = st_as_sfc(st_bbox(fmu_shape))

us_map<-ggplot() + 
  geom_sf(data = us_states, fill = "white") + 
  geom_sf(data = MI_box, fill = NA, color = "red", size = 1.2) +
  theme_void()

us_map

inset_map = ggdraw() +
  draw_plot(map_of_data) +
  draw_plot(us_map, x = 0.05, y = 0.35, width = 0.3, height = 0.3)

inset_map

ggsave( inset_map, 
          device = "png", 
          filename = "figures/fig1_map_of_data.png", 
          dpi = 600, height = 8, width = 8, units = "in",
          bg="#ffffff") #sets background to white 

