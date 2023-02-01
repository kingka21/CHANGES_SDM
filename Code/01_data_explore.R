# code written by Katelyn King 
# created: April 2021 
#Data exploration 

library(PerformanceAnalytics)

#### look at correlation among driver variables #### 
driver_vars<-read.csv("Data/contemp_lmb_dat.csv") %>% 
  distinct(new_key, .keep_all = TRUE)

driver_nonas<-na.omit(driver_vars) # find out how many rows have all of the variables #371 

my_data <- driver_vars[, c(8,10,13,14,25,27:37,46)] 
PerformanceAnalytics::chart.Correlation(my_data, histogram=TRUE, pch=19)
 
#pick dd or surface temp 
my_data <- driver_vars[, c(8,10,25,29,32,35,42)] #surface temp year 
PerformanceAnalytics::chart.Correlation(my_data, histogram=TRUE, pch=19)
my_data <- driver_vars[, c(8,10,25,29,32,35,36)] #DD mean
PerformanceAnalytics::chart.Correlation(my_data, histogram=TRUE, pch=19)

#### plot drivers vs. response variables #### 
my_data <- driver_vars[, c(5,8,10,25,29,32,35,42)] #surface temp year 
par(mfrow = c(2,4)) #15 drivers 
#ggplot 
x_axis_vars = names(my_data[,])
for (i in 1:length(x_axis_vars)) {
  x_vars = x_axis_vars[i]
  print(ggplot(my_data, aes(x=.data[[x_vars]], y= fish_count_new)) + 
    geom_point() + 
      geom_smooth(method="lm") + 
      xlab(paste(x_vars))
    )
  #Sys.sleep(2) #this makes it shuffle through 
 }


#### lake area distribution across gears #### 
count_data<-read.csv("Data/contemp_lmb_dat.csv")

ggplot(count_data, aes(x=log(lake_area_ha), color=gear2, linetype=gear2)) +
  geom_density() + 
  scale_linetype_manual(values=c("twodash", "dotted", "solid", "longdash"))+
  scale_color_manual(values=c('#000000','#E69F00', "#009E73","#56B4E9" )) + 
  theme_bw()

#### temporal variation in gear use #### 
ggplot(count_data, aes(x=year, color=gear2, linetype=gear2)) +
  geom_density() + 
  scale_linetype_manual(values=c("twodash", "dotted", "solid", "longdash"))+
  scale_color_manual(values=c('#000000','#E69F00', "#009E73","#56B4E9" )) + 
  theme_bw()

#### map gear use to look for spatial autocorrelation #### 
#pull Michigan map from the map package and plot points just to examine
library(ggmap)
library(mapdata)
#create color palette 
colfunc = colorRampPalette(c("cyan", "red")) 
ncols = 208 
cols = colfunc(208)
par(mfrow=c(1,2)) 


MI_map<-map_data(map = "state", region = c("michigan"))  #pull out the usa map
map<-ggplot(data = MI_map) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + #this fill MI white and outline is black
  coord_fixed(1.3) 

#map LMB shock
dat_shock<-filter(count_data, gear2 == "SHOCK")
map + 
    geom_point(data=dat_shock, aes(x = LONG_DD, y = LAT_DD, colour = c(log(fish_count_new)))) + 
    scale_colour_gradient(low = "cyan", high = "red") + 
    labs(color="shock catch")  #changes the labels on the legend 

#map LMB net 
dat_seine<-filter(count_data, gear2 == "SEINE")
map + 
  geom_point(data=dat_seine, aes(x = LONG_DD, y = LAT_DD, colour = c(log(fish_count_new)))) + 
  scale_colour_gradient(low = "cyan", high = "red") + 
  labs(color="seine catch")  #changes the labels on the legend 

dat_gill<-filter(count_data, gear2 == "GILL")
map + 
  geom_point(data=dat_gill, aes(x = LONG_DD, y = LAT_DD, colour = c(log(fish_count_new)))) + 
  scale_colour_gradient(low = "cyan", high = "red") + 
  labs(color="gill catch")  #changes the labels on the legend 

dat_fyke<-filter(count_data, gear2 == "FT_NET")
map + 
  geom_point(data=dat_fyke, aes(x = LONG_DD, y = LAT_DD, colour = c(log(fish_count_new)))) + 
  scale_colour_gradient(low = "cyan", high = "red") + 
  labs(color="fyke catch")  #changes the labels on the legend 


#boxplot
boxplot(count_data$fish_count_new ~ count_data$gear2)


