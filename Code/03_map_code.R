### Map figures #### 
#devtools::install_github("dkahle/ggmap")
library(ggmap)
library(stringr)
library(devtools)
library(mapdata)
library(dplyr)
library(ggpubr)

#### For Figure 1 sampling locations #### 
#get the lat/long points 
cont_points<-read.csv("Data/lmb_dat_for_model_aug30.csv") %>%
  dplyr::select( new_key, LAT_DD, LONG_DD, fish_count_new, lake_area_m2, maxdepth_m, surface_temp_mean, secchi_m, day_of_year, ws_forest_prop, ws_wetland_prop) %>%
  distinct(new_key, .keep_all = TRUE) %>% 
  na.omit()%>% 
  mutate(data = 'contemporary') %>%
  select(new_key, LONG_DD, LAT_DD, data)

# get historical data 
hist_dat<-read.csv("Data/historical_lmb2_with_secchi.csv") %>%
  na.omit() %>% 
  distinct(new_key, .keep_all = TRUE)

points<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/IFR_Lake_Points.csv")%>%
  rename(new_key = NEW_KEY) %>% 
  select(new_key, LONG_DD, LAT_DD)

hist_points<-hist_dat  %>% 
  left_join(points) %>% 
  select(new_key, LONG_DD, LAT_DD)%>% 
  mutate(data = 'historical')



#shapefile with the FMU boundaries 
library(sf)
fmu_shape<- st_read(dsn = "/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/MDNR_FMUs", layer = "Michigan_DNR_Fisheries_Management_Units") 

#combine data for map 
map_data<-rbind(cont_points, hist_points)

ggplot() +
  geom_sf(data=fmu_shape, aes(), fill = "white") + 
  geom_point(data=map_data, aes(x = LONG_DD, y = LAT_DD, color=c(data) )) +
  theme_bw() +  
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank()) +
  clean_theme() +   
  scale_colour_manual(values = c("darkblue", "lightsalmon"), 
                      name= '') +
  theme(legend.position = c(0.20, 0.25), 
        legend.text = element_text(size=12)) 

map_of_data

ggsave("/Users/katelynking/Desktop/map_of_data.png", map_of_data)


####map of lmb abund #### 
model1_cont_preds<-read.csv("Data/output/model1_cont_abund_pred.csv")
model1_hist_preds<-read.csv("Data/output/model1_hist_abund_pred.csv")
model2_cont_preds<-read.csv("Data/output/model2_cont_abund_pred.csv")
model2_hist_preds<-read.csv("Data/output/model2_hist_abund_pred.csv")


#get MI basemap 
MI_basemap<-map_data("state",  region = c("michigan"))  # select michigan 
map<-ggplot(data = MI_basemap) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + #this fill MI white and outline is black
  coord_fixed(1.3) 

#model 1 hist
mod1_hist<-map+ geom_point(data=model1_hist_preds, aes(x = LONG_DD, y = LAT_DD, colour = c(hist.means))) + 
  scale_colour_gradient(low = "yellow", high = "darkblue")  + 
  labs(color="relative density", tag = "a") + #changes the labels on the legend and add panel letter
  ylab(NULL) + xlab(NULL) + 
  theme_bw() + theme(legend.position = "bottom", 
                     plot.tag=element_text(), 
                     plot.tag.position = c(0.1,0.95)) 

#model1 contemp
mod1_cont<-map+ geom_point(data=model1_cont_preds, aes(x = LONG_DD, y = LAT_DD, colour = c(cont.means.abund))) + 
                scale_colour_gradient(low = "yellow", high = "darkblue")  + 
                labs(color="relative density", tag = "b") + #changes the labels on the legend and add panel letter
                  ylab(NULL) + xlab(NULL) + 
               theme_bw() + theme(legend.position = "bottom", 
                     plot.tag=element_text(), 
                     plot.tag.position = c(0.1,0.95)) 

#model 2 hist
mod2_hist<-map+ geom_point(data=model2_hist_preds, aes(x = LONG_DD, y = LAT_DD, colour = c(hist.means))) + 
  scale_colour_gradient(low = "yellow", high = "darkblue")  + 
  labs(color="relative density", tag = "c") + #changes the labels on the legend and add panel letter
  theme_bw() + theme(legend.position = "bottom", 
                     plot.tag=element_text(), 
                     plot.tag.position = c(0.1,0.95)) + 
  ylab(NULL) + xlab(NULL)

#model 2 contemp
mod2_cont<-map+ geom_point(data=model2_cont_preds, aes(x = LONG_DD, y = LAT_DD, colour = c(cont.means.abund))) + 
  scale_colour_gradient(low = "yellow", high = "darkblue", limits=c(0,15))  + 
  labs(color="relative density", tag = "d") + #changes the labels on the legend and add panel letter
  theme_bw() + theme(legend.position = "bottom", 
                     plot.tag=element_text(), 
                     plot.tag.position = c(0.1,0.95)) +
  ylab(NULL) + xlab(NULL)

mapfig7<-cowplot::plot_grid(mod1_hist, mod1_cont, mod2_hist, mod2_cont)


ggsave(plot=mapfig7, 
       device="png", 
       filename = "figures/mapfig7.png", 
       dpi=600, height = 8, width=8, units ="in", 
       bg="#ffffff")


### HISTOGRAMS ####
model2_cont_preds<-read.csv("Data/output/model2_cont_abund_pred.csv")
model2_hist_preds<-read.csv("Data/output/model2_hist_abund_pred.csv")

c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")


hist_1<-hist(model1_cont_preds$sum_dens)
hist_2<-hist(model1_hist_preds$sum_dens)
plot(hist_1, col = c1, xlim=c(0,10), xlab="relative density")
plot(hist_2, col = c2, add=T)

hist_3<-hist(log(model2_cont_preds$sum_dens))
hist_4<-hist(log(model2_hist_preds$sum_dens), breaks=5)
plot(hist_3, col = c1, xlab="relative density (log)")
plot(hist_4, col = c2, add=T)

#ggplot 
model1_cont_hist<-model1_cont_preds %>% 
  rename(mean= cont.means.abund) %>% 
  select(new_key, mean) %>% 
  mutate(time = "contemporary")

model1_hist_hist<-model1_hist_preds %>% 
  rename(mean= hist.means) %>% 
  select(new_key, mean) %>% 
  mutate(time = "historical")

model_1<-rbind(model1_cont_hist, model1_hist_hist)

plot_multi_histogram <- function(df, feature, label_column) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.25, position="identity", aes(y = ..density..), color="black", bins = 10) +
    labs(x=feature, y = "Density")
  plt + guides(fill=guide_legend(title=label_column))
}

plot_multi_histogram(model_1, 'mean', 'time')

#### boxplot ####

