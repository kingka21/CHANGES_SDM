# code written by Katelyn King 
# created: April 2021 
#Data exploration 

library(PerformanceAnalytics)

#### look at correlation among driver variables #### 
driver_vars<-read.csv("Data/driver_varibles.csv")

driver_nonas<-na.omit(driver_vars) # find out how many rows have all of the variables #371 

my_data <- driver_vars[, c(5,6,13:28)] 
PerformanceAnalytics::chart.Correlation(my_data, histogram=TRUE, pch=19)
# use houses instead of docks... don't think Canada has the dock data and they are highly correlated  
#pick dd or surface temp 
my_data <- driver_vars[, c(5,6,13:21, 28,29)] 

#### plot drivers vs. response variables #### 
par(mfrow = c(5, 3)) #15 drivers 
#ggplot 
x_axis_vars = names(dat_shock[, 51:66])
for (i in 1:length(x_axis_vars)) {
  x_vars = x_axis_vars[i]
  print(ggplot(dat_shock, aes(x=.data[[x_vars]], y= SHOCK)) + 
    geom_point() + 
      geom_smooth(method="lm") + 
      xlab(paste(x_vars))
    )
  #Sys.sleep(2) #this makes it shuffle through 
 }

#base R plot 
x_axis_vars = names(dat_shock[, 51:66])

for(i in 1:length(x_axis_vars)){
  curr_x_axis_var = x_axis_vars[i] #set the current x axis value 
  plot(y = dat_shock$SHOCK,
       x = dat_shock[[curr_x_axis_var]],
       xlab = paste(curr_x_axis_var))
  Sys.sleep(2) #this makes it shuffle through
}

####plot catch data by gear #### 
lmb<-filter(count_data, species == 'LMB')
ggplot(lmb, aes(x=log(fish_count), color=gear2, linetype=gear2)) +
  geom_density() + 
  scale_linetype_manual(values=c("twodash", "dotted", "solid", "longdash"))+
  scale_color_manual(values=c('#000000','#E69F00', "#009E73","#56B4E9" )) + 
  theme_bw()  + 
  ggtitle("LMB")


wae<-filter(count_data, species == 'WAE')
ggplot(wae, aes(x=log(fish_count), color=gear2, linetype=gear2)) +
  geom_density() + 
  scale_linetype_manual(values=c("twodash", "dotted", "solid", "longdash"))+
  scale_color_manual(values=c('#000000','#E69F00', "#009E73","#56B4E9" )) + 
  theme_bw()+ 
  ggtitle("WAE")

cis<-filter(count_data, species == 'CIS')
ggplot(cis, aes(x=log(fish_count), color=gear2, linetype=gear2)) +
  geom_density() + 
  scale_linetype_manual(values=c("twodash", "dotted", "solid", "longdash"))+
  scale_color_manual(values=c('#000000','#E69F00', "#009E73","#56B4E9" )) + 
  theme_bw()+ 
  ggtitle("CIS")

#### area distribution across gears #### 
#join tables for LMB 
new<-count_data %>%
  group_by(new_key, gear2) %>%
  summarise(sum = sum(fish_count))

dat<-left_join(new, drivers ) %>% 
  select(new_key, gear2, lake_area_ac)

ggplot(dat, aes(x=lake_area_ac, color=gear2, linetype=gear2)) +
  geom_density() + 
  scale_linetype_manual(values=c("twodash", "dotted", "solid", "longdash"))+
  scale_color_manual(values=c('#000000','#E69F00', "#009E73","#56B4E9" )) + 
  theme_bw()

#### temporal variation in gear use #### 
new<-count_data %>%
  group_by(new_key, gear2) %>%
  summarise(sum = sum(fish_count))

MI_data_no_dups
MI_data<- MI_data_no_dups[!duplicated(paste(MI_data_no_dups$new_key)),] %>%
  select(new_key, year)

dat<-left_join(new, MI_data, by=c("new_key" )) %>% 
  select(new_key, gear2, year)

ggplot(dat, aes(x=year, color=gear2, linetype=gear2)) +
  geom_density() + 
  scale_linetype_manual(values=c("twodash", "dotted", "solid", "longdash"))+
  scale_color_manual(values=c('#000000','#E69F00', "#009E73","#56B4E9" )) + 
  theme_bw()
#### map variables to look for spatial autocorrelation #### 
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

####map of sp obs across the state (just presence) ####
#note these are obs before removing ones that dont have predictor variables. so numbers go down 
# lmb 
map + 
  geom_point(data=lmb_count_dat, aes(x = LONG_DD, y = LAT_DD)) +
  ggtitle("LMB n=378")

#wae
map + 
  geom_point(data=wae_count_dat, aes(x = LONG_DD, y = LAT_DD)) +
  ggtitle("WAE n=226")
#cis
map + 
  geom_point(data=cis_count_dat, aes(x = LONG_DD, y = LAT_DD)) +
  ggtitle("CIS n=32")
