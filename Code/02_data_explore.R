# code written by Katelyn King 
# created: April 2021 
#Data exploration 

library(PerformanceAnalytics)

driver_vars<-read.csv("Data/driver_varibles.csv")

driver_nonas<-na.omit(driver_vars) # find out how many rows have all of the variables #371 
#need to figure out what to do with the -99 (e.g. does 0 make more sense)
my_data <- driver_vars[, c(5,6,13:28)] 
PerformanceAnalytics::chart.Correlation(my_data, histogram=TRUE, pch=19)
# use houses instead of docks... don't think canada has the dock data and they are highly correlated  
#pick dd or surface temp 
my_data <- driver_vars[, c(5,6,13:21, 28,29)] 


##### map variables to look for spatial autocorrelation #### 
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
map + 
    geom_point(data=dat_shock, aes(x = LONG_DD, y = LAT_DD, colour = c(log(SHOCK)))) + 
    scale_colour_gradient(low = "cyan", high = "red") + 
    labs(color="shock cpue")  #changes the labels on the legend 

#map LMB net 
map + 
  geom_point(data=dat_net, aes(x = LONG_DD, y = LAT_DD, colour = c(log(NET)))) + 
  scale_colour_gradient(low = "cyan", high = "red") + 
  labs(color="net cpue")  #changes the labels on the legend 

#example plotting the predictor (x) and the response (Y) to look for spatial pattern
for (i in 1:2){ 
  if (i==1) value = log(dat_net$NET) 
  if (i==2) value = dat_shock$SHOCK
  value = value-min(value) 
  value = 1+(ncols-1)*value/max(value) 
  plot(dat_pres$LONG_DD,dat_pres$LAT_DD,col=cols[value],pch=16,main=c("net","shock")[i], asp=1)
}

#example of predictors 
for (i in 1:2){ 
  if (i==1) value = dat_net$dd_mean
  if (i==2) value = dat_net$lake_order
  value = value-min(value) 
  value = 1+(ncols-1)*value/max(value) 
  plot(dat_pres$LONG_DD,dat_pres$LAT_DD,col=cols[value],pch=16,main=c("DD","order")[i], asp=1)
}

#boxplot
boxplot(dat_shock$SHOCK)
boxplot(dat_net$NET)
boxplot(log(dat_shock$SHOCK), ylab="logshock")
boxplot(log(dat_net$NET), ylab="lognet")
