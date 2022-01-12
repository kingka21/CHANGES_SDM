#### Ontario Lake and Watershed Derived Variables #### 

#### load libraries ####
library(sp)
library(sf)
library(dplyr)
library(raster)
library(ggplot2)
library(exactextractr)

#### summarize raster data relative to WS  #### 
#http://www.wvview.org/spatial_analytics/Raster_Analysis/_site/index.html

#read in raster data using the raster package and shp data using the sf package 
lc <- raster("/Users/katelynking/Desktop/UofM/CHANGES/ON_data/Provincial_Land_Cover_1996/Provincial_Land_Cover_1996/plc1996_28/") #LULC data 
elev <- raster("/Users/katelynking/Desktop/UofM/CHANGES/ON_data/PDEM_data/PDEM_North.tif") #elevation data for north province 
elev2<-raster("/Users/katelynking/Desktop/UofM/CHANGES/ON_data/PDEM_data/PDEM_South.tif") #elevation data for the south province
ws1<- st_read(dsn = "/Users/katelynking/Desktop/UofM/CHANGES/ON_data/BsM_lake_catchments", layer = "Cy1_BsM_lakecatchments")
ws2<- st_read(dsn = "/Users/katelynking/Desktop/UofM/CHANGES/ON_data/BsM_lake_catchments", layer = "Cy2_BsM_lakecatchments")

#* elevation #### 
#extract mean elevation for each watershed North and South
ws_elev_mean <- st_as_sf(extract(elev, ws1, fun=mean, sp=TRUE)) #north elev
ws1_elev_mean <- st_set_geometry(ws_elev_mean, NULL)
write.csv(ws1_elev_mean, 'Data/ON_data/ws_elev_mean_cy1.csv', row.names = FALSE)

ws_elev_mean_south <- st_as_sf(extract(elev2, ws1, fun=mean, sp=TRUE)) #south elev
ws1_elev_mean_south <- st_set_geometry(ws_elev_mean_south, NULL)
write.csv(ws1_elev_mean_south, 'Data/ON_data/ws_elev_mean_cy1_south.csv', row.names = FALSE)

#watersheds from cycle 2 
ws2_elev_mean <- st_as_sf(extract(elev, ws2, fun=mean, sp=TRUE)) #takes like 30 min to run 
ws2_elev_mean <- st_set_geometry(ws2_elev_mean, NULL)
write.csv(ws2_elev_mean, 'Data/ON_data/ws_elev_mean_cy2.csv', row.names = FALSE) 

ws2_elev_mean_south <- st_as_sf(extract(elev2, ws2, fun=mean, sp=TRUE)) #south elev  
ws2_elev_mean_south <- st_set_geometry(ws2_elev_mean_south, NULL)
write.csv(ws2_elev_mean_south, 'Data/ON_data/ws_elev_mean_cy2_south.csv', row.names = FALSE) 

elev_cy1_s<-read.csv('Data/ON_data/ws_elev_mean_cy1_south.csv') %>% 
            dplyr::select(My_intwbdy, My_intwb_1, Shape_Area, PDEM_South) %>% 
  rename(ws_area=Shape_Area, ws_mean_elev=PDEM_South)
elev_cy1_n<-read.csv('Data/ON_data/ws_elev_mean_cy1.csv')%>% 
  dplyr::select(My_intwbdy, My_intwb_1, Shape_Area, PDEM_North)%>% 
  rename(ws_area=Shape_Area, ws_mean_elev=PDEM_North)
elev_cy2_s<-read.csv('Data/ON_data/ws_elev_mean_cy2_south.csv')%>% 
  dplyr::select(My_intwbdy, My_intwb_1, Shape_Area, PDEM_South)%>% 
  rename(ws_area=Shape_Area, ws_mean_elev=PDEM_South)
elev_cy2_n<-read.csv('Data/ON_data/ws_elev_mean_cy2.csv') %>% 
  dplyr::select(My_intwbdy, My_intwb_1, Shape_Area, PDEM_North)%>% 
  rename(ws_area=Shape_Area, ws_mean_elev=PDEM_North)

#join all 
ws_elev<-rbind(elev_cy1_s, elev_cy1_n, elev_cy2_s, elev_cy2_n) 

ws_elev$ws_id<-ifelse(ws_elev$My_intwbdy == 0, ws_elev$My_intwb_1, ws_elev$My_intwbdy)
ws_elev<-dplyr::select(ws_elev, ws_id, ws_area, ws_mean_elev) %>%
            na.omit(ws_mean_elev)

#look at duplicates : many duplicates because of sout, north, and cycles overlap
#duplicates<-ws_elev %>% group_by(ws_id) %>% filter(n() > 1)
elev_ws_cy1and2<- ws_elev[!duplicated(paste(ws_elev$ws_id)),] #1727 unique ws 
write.csv(elev_ws_cy1and2, "Data/ON_data/elev_ws_cy1and2.csv", row.names = FALSE)

#* LULC #####################
#watersheds have two columns with ws ids - merge these into one column 
ws1$ws_id<-ifelse(ws1$My_intwbdy == 0, ws1$My_intwb_1, ws1$My_intwbdy )
ws1<-dplyr::select(ws1, ws_id)
ws2$ws_id<-ifelse(ws2$My_intwbdy == 0, ws2$My_intwb_1, ws2$My_intwbdy )
ws2<-dplyr::select(ws2, ws_id)

#start by cropping the LULC data to match the boundaries of the watershed. 
landcov_ws <- crop(lc, ws1)

#exact_extract function is really fast! 
ex <- exact_extract(landcov_ws,ws1,include_cols="ws_id") 
ex <- do.call(rbind,ex)
lulc_cy1 <- ex %>% 
  group_by(ws_id) %>% 
  summarise(
  landuse = sort(unique(value)),
  count = table(value),
  freq = round(table(value) / length(value), 3)
)

#pull into a data frame 
proportion<-as.data.frame(lulc_cy1$freq)
proportion$ws_id<-lulc_cy1$ws_id

#join with the class names and spread the data 
lulc_class<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/ON_data/Provincial_Land_Cover_1996/lulc_class.csv", stringsAsFactors = TRUE)
lulc_class$landuse<-as.factor(lulc_class$landuse)
cy1_lulc<-left_join(proportion, lulc_class, by=c('Var1' = 'landuse'))

lc1_wide <- tidyr::pivot_wider(data= cy1_lulc, 
             id_cols = c(ws_id),
             names_from = name, 
             values_from = Freq, 
             values_fill = list(Freq = 0)) #replace NAs with 0)

#second cycle set of watersheds 
landcov_ws2 <- crop(lc, ws2)

#exact_extract function 
ex2 <- exact_extract(landcov_ws2, ws2,include_cols="ws_id") 
ex2 <- do.call(rbind,ex2)
lulc_cy2 <- ex2 %>% group_by(ws_id) %>% summarise(
  landuse = sort(unique(value)),
  count = table(value),
  freq = round(table(value) / length(value), 3)
)
#pull into a data frame 
proportion2<-as.data.frame(lulc_cy2$freq)
proportion2$ws_id<-lulc_cy2$ws_id

cy2_lulc<-left_join(proportion2, lulc_class, by=c('Var1' = 'landuse'))

lc2_wide <- tidyr::pivot_wider(data= cy2_lulc, 
                               id_cols = c(ws_id),
                               names_from = name, 
                               values_from = Freq, 
                               values_fill = list(Freq = 0)) #replace NAs with 0)

#add missing columns so we can bind Coastal_mudflats, Intertidal_marsh, Supertidal_marsh, Tundra_heath

lc2_wide$Coastal_mudflats<-0
lc2_wide$Intertidal_marsh<-0
lc2_wide$Supertidal_marsh<-0
lc2_wide$Tundra_heath<-0

#bind the two watershed datasets 
land_usecover_cy12<-rbind(lc1_wide, lc2_wide) #783 ws and 728 ws 

#look at duplicates 
#duplicates<-land_usecover_cy12 %>% group_by(ws_id) %>% filter(n() > 1)
#all other duplicates look like they match exactly ~500 lakes were sampled in both cycles 
land_usecover_cy12<- land_usecover_cy12[!duplicated(paste(land_usecover_cy12$ws_id)),] #1757 unique ws 

write.csv(land_usecover_cy12, "Data/ON_data/lulc_ws_cycle1and2.csv", row.names = FALSE)
land_usecover_cy12<-read.csv("Data/ON_data/lulc_ws_cycle1and2.csv") 

#other method 
#https://bookdown.org/mcwimberly/gdswr-book/combining-raster-and-vector-data-2-land-cover-data.html
#landcov_ws <- mask(landcov_ws, ws1) #this takes so long I never got it to finish 
#estimate the relative areas of different classes within each polygon
#ws_raster <- rasterize(ws1,landcov_ws, field = "ID") #rasterize the watersheds
#xtab <- crosstab(ws_raster, landcov_ws) #function to calculate the number of cells of each lc in each watershed polygon
#lcsum <- data.frame(as.matrix(xtab)) #create table of sums 


#### assign each lake to a watershed to extract data #### 

#ON lakes from cylcle 1 and 2 shapefiles 
lakes_shape<- st_read(dsn = "/Users/katelynking/Desktop/UofM/CHANGES/ON_data/BsM_lake_catchments", layer = "Cy1_BsM_lakes") %>%
  dplyr::select(WbyLID)
lakes_shape2<- st_read(dsn = "/Users/katelynking/Desktop/UofM/CHANGES/ON_data/BsM_lake_catchments", layer = "Cy2_BsM_lakes")%>%
  dplyr::select(WbyLID)
#crs(lakes_shape) 

#transform objects to spatial for over to work 
ws1_sp<-as(ws1, Class="Spatial")
ws2_sp<-as(ws2, Class="Spatial")
lake_shape_sp<-as(lakes_shape, Class="Spatial")
lake_shape_sp2<-as(lakes_shape2, Class="Spatial")

prjnew<-crs(ws1_sp) #get the crs from the ws file

#over function: at the spatial locations of object x retrieves the indexes or attributes from spatial object y
#over function only works if they are exact same CRS 
lake_shape_sp<- spTransform(lake_shape_sp, prjnew) #project to the exact same 
lake_shape_sp2<- spTransform(lake_shape_sp2, prjnew) #project to the exact same 

ws_lakeshp<-over(lake_shape_sp, ws1_sp)
ws_lakeshp2<-over(lake_shape_sp2, ws2_sp)

# Add data back (all match but one lake!) 
ws_lakeshp$lakeid <- lake_shape_sp$WbyLID
ws_lakeshp2$lakeid <- lake_shape_sp2$WbyLID

ws_match<-rbind(ws_lakeshp,ws_lakeshp2 )
ws_match<- ws_match[!duplicated(paste(ws_match$lakeid)),] #953 unique lakes 
write.csv(ws_match, "Data/ON_data/ws_lake_ids.csv", row.names = FALSE)

############################ #### #### #### #### #### #### #### #### #### #### #### 
#### base maps for the US and Canada from the raster package - pulls from the GADM site ####
#will download gdb to your directory 
#select states and provinces of interest 
states    <- c('Michigan')
provinces <- c("Ontario")

us <- getData("GADM",country="USA",level=1)
canada <- getData("GADM",country="CAN",level=1)

us.states <- us[us$NAME_1 %in% states,]
ca.provinces <- canada[canada$NAME_1 %in% provinces,]

ggplot(us.states,aes(x=long,y=lat,group=group))+
  geom_path()+
  geom_path(data=ca.provinces)+
  geom_path(data=lakes_shape)+
  coord_map() 

#get lake points 
lakes<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/ON_data/BsM_data_2021/lake_info_BsM_2021.csv") %>%
  dplyr::select(BsM_Cycle, WbyLID, WbyLat, WbyLong) %>%
  rename(lakeid=WbyLID, Lat=WbyLat, Lon=WbyLong) %>%
  filter(BsM_Cycle == 1 | BsM_Cycle == 2) #select only 1 or 2, this is what we have WS for 
# some lakes are in multiple cycles and some have different days of netting in the same cycle  
#just need lat/lon for lake - so select only one row of each lake n = 1147
lakes<- lakes[!duplicated(paste(lakes$lakeid)),] 

# project points  
lakes_sf = st_as_sf(lakes, coords= c("Lon" , "Lat"), crs = "+proj=longlat +datum=NAD83 +ellps=GRS80") %>% 
  dplyr::select(lakeid)

#plot ws and points to see if they match 
ggplot() + 
  geom_sf(data = ws1) +
  geom_sf(data = lakes_sf)
