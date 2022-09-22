## code written by Katelyn King 
# created: April 2021 
 
##### matching data to LAGOS lakes
#### load libraries ####
library(rgdal)
library(sp)
library(dplyr)
library(raster)
library(ggplot2)

#* match fish data to lagos lakes using lake link#### 
MIpoints<-dplyr::select(MI_data, new_key, ihdlkid, long_dd, lat_dd) %>% 
  distinct(new_key) #select only one row of each lake
MIpoints[is.na(MIpoints)] <- 0
MIpoints$ihdlkid<-as.factor(MIpoints$ihdlkid)
#MI_lagos_link<-left_join(MIpoints, lake_link, by = c('ihdlkid' = 'lake_nhdid'))
#MI_lagos_link<- MI_lagos_link[!duplicated(paste(MI_lagos_link$new_key)),] #27 lakes don't match 

#* match fish data to lagos lakes using shapefile of lakes #### 
#read in polygons 
#4 ha LAGOS-NE
lakes_shape<- readOGR(dsn = "/Users/katelynking/Desktop/Cont Limno/LAGOS_NE_All_Lakes_4ha", layer = "LAGOS_NE_All_Lakes_4ha")

#>1 ha LAGOS-US 
#lakes_shape<- readOGR(dsn = "/Users/katelynking/Desktop/Cont Limno/LAGOS_US/LAGOS_US_polys", layer = "LAGOS_US_All_lakes_1ha_v0.5")

#change to an sf object   
#sf_US<-sf::st_as_sf(lakes_shape)

#select out Michigan polys only 
#MI_poly<-dplyr::filter(sf_US, STATE == "MI")  

#convert it back to a spatial polygons data frame
#MI_poly.sp<-as(MI_poly, "Spatial")
#crs(MI_poly.sp) #shows me the projection of the shapefiles so that I can project the same to the points 
#writeOGR(MI_poly.sp, ".", "MI_poly_1ha", 
#       driver = "ESRI Shapefile")

#>1 ha LAGOS-US Michigan only 
#MI_poly.sp<- readOGR(dsn = "/MI_poly_shapefile", layer = "MI_poly_1ha")

#get the lat long data for the Michigan lakes with fish data 
x <- MIpoints$long_dd
y <- MIpoints$lat_dd

#make a dataframe of the coordinates and project
lake.ll <- SpatialPoints(data.frame(x=x,y=y), proj4string=CRS("+proj=longlat +datum=NAD83"))   

prjnew <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80
+towgs84=0,0,0,0,0,0,0 +units=m +no_defs ")

MI.aea <- spTransform(lake.ll, prjnew)
plot(MI.aea)

#the over function, where x = "SpatialPoints", y = "SpatialPolygonsDataFrame" from which the geometries are queried
#join with LAGOS-US 1ha shapefile 
#MI_fish<-sp::over(MI.aea, MI_poly.sp) %>%
# dplyr::select(lagoslakei, STATE, Hectares) %>%
#rename(lagoslakeid = lagoslakei, state = STATE, area_ha = Hectares)
#MI_fish$new_key<-paste(MIpoints$new_key)  #16 lakes don't match 

#join with the 4ha shapefile from LAGOS_NE
MI_LAGOS<-sp::over(MI.aea, lakes_shape) %>%
  dplyr::select(lagoslakei, STATE, COUNTY_Nam, Lake_Area_) %>%
  rename(lagoslakeid = lagoslakei, state = STATE, county= COUNTY_Nam, area_ha = Lake_Area_)
MI_LAGOS$new_key<-as.factor(paste(MIpoints$new_key))  #10 lakes don't match 

#match fish data to the lagoslakeids 
MI_data<-left_join(MI_data, MI_LAGOS, by = "new_key")

#### investigate lakes that don't have nhdids and dont match polys 
#no match with poly or nhd join new_key: 22-123, 29-59, 29-101, 3-297,  52-839, 79_44,  81-194, 82-235, 

#MI_invest<- MI_data[!duplicated(paste(MI_data$new_key)),]
#MI_invest<-MI_invest[is.na(MI_invest$lagoslakeid),]
#get the lat long to map data that dont match 
#x <- MI_invest$long_dd
#y <- MI_invest$lat_dd

#make a dataframe of the coordinates and project
#coords<-data.frame(x=x,y=y)
#lake.ll <- SpatialPoints(coords, proj4string=CRS("+proj=longlat +datum=NAD83"))   
#mapview::mapview(lake.ll) + mapview::mapview(lakes_shape)

#assign lagoslakeids from investigation 
#63-2 = lagoslakeid 1997
#48-485 = 75909
MI_data_lagos <- MI_data %>%
  mutate(lagoslakeid2 = case_when(new_key =="63-2" ~ "1997", #if ~ then 
                                  new_key =="48-485" ~ "75909", 
                                  TRUE ~ as.character(MI_data$lagoslakeid) #else use the lagoslakeid that exists 
  )) 

#### add lagos data #### 
library(LAGOSNE)
#load the data and call it lg
lg <- lagosne_load()

IWS_LULC <- lg$iws.lulc #this pulls out the watershed land use metrics
IWS_LULC <- data.frame(lagoslakeid=IWS_LULC$lagoslakeid,
                       developed_open2006_pct=IWS_LULC$iws_nlcd2006_pct_21,
                       developed_low2006_pct=IWS_LULC$iws_nlcd2006_pct_22,
                       developed_med2006_pct=IWS_LULC$iws_nlcd2006_pct_23,
                       developed_high2006_pct=IWS_LULC$iws_nlcd2006_pct_24,
                       deciduous2006_pct=IWS_LULC$iws_nlcd2006_pct_41,
                       evergreen2006_pct=IWS_LULC$iws_nlcd2006_pct_42,
                       mixedforest2006_pct=IWS_LULC$iws_nlcd2006_pct_43,
                       scrubshrub2006_pct=IWS_LULC$iws_nlcd2006_pct_52,
                       grasslands_herbaceous2006_pct=IWS_LULC$iws_nlcd2006_pct_71,
                       pasture_hay2006_pct=IWS_LULC$iws_nlcd2006_pct_81,
                       rowcrops2006_pct=IWS_LULC$iws_nlcd2006_pct_82,
                       roaddensity_mperha=IWS_LULC$iws_roaddensity_density_mperha,
                       wetland_woody2006_pct=IWS_LULC$iws_nlcd2006_pct_90,
                       wetland_emergent2006_pct=IWS_LULC$iws_nlcd2006_pct_95,
                       damdens_ws = IWS_LULC$iws_damdensity_pointspersqkm)

#add up all the developed categories and call it urban
IWS_LULC$urban_ws <- IWS_LULC$developed_open2006_pct +
  IWS_LULC$developed_low2006_pct + IWS_LULC$developed_med2006_pct +
  IWS_LULC$developed_high2006_pct

IWS_LULC$forest_ws <- IWS_LULC$deciduous2006_pct +
  IWS_LULC$evergreen2006_pct + IWS_LULC$mixedforest2006_pct

IWS_LULC$agriculture_ws <- IWS_LULC$pasture_hay2006_pct +
  IWS_LULC$rowcrops2006_pct

IWS_LULC$wetland_ws <- IWS_LULC$wetland_woody2006_pct +
  IWS_LULC$wetland_emergent2006_pct

IWS_LULC$shrub_ws<-IWS_LULC$scrubshrub2006_pct +
  IWS_LULC$grasslands_herbaceous2006_pct

#change road density to km/sqkm
IWS_LULC$roaddens_ws<-IWS_LULC$roaddensity_mperha / 10

#select just the columns I want
IWS_table<-dplyr::select(IWS_LULC, lagoslakeid, damdens_ws, urban_ws, forest_ws, agriculture_ws, wetland_ws, roaddens_ws, shrub_ws)

#same thing for the local buffer around lakes (100 m)
Buff100_LULC <- lg$buffer100m.lulc
Buff100_LULC <- data.frame(lagoslakeid=Buff100_LULC$lagoslakeid,
                           buff_developed_open2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_21,
                           buff_developed_low2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_22,
                           buff_developed_med2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_23,
                           buff_developed_high2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_24,
                           buff_deciduous2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_41,
                           buff_evergreen2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_42,
                           buff_mixedforest2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_43,
                           buff_scrubshrub2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_52,
                           buff_grasslands_herbaceous2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_71,
                           buff_pasture_hay2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_81,
                           buff_rowcrops2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_82,
                           buff_roaddensity_mperha=Buff100_LULC$buffer100m_roaddensity_density_mperha,
                           buff_wetland_woody2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_90,
                           buff_wetland_emergent2006_pct=Buff100_LULC$buffer100m_nlcd2006_pct_95)

Buff100_LULC$urban_buf <- Buff100_LULC$buff_developed_open2006_pct +
  Buff100_LULC$buff_developed_low2006_pct +
  Buff100_LULC$buff_developed_med2006_pct +
  Buff100_LULC$buff_developed_high2006_pct

Buff100_LULC$forest_buf <- Buff100_LULC$buff_deciduous2006_pct +
  Buff100_LULC$buff_evergreen2006_pct +
  Buff100_LULC$buff_mixedforest2006_pct

Buff100_LULC$agriculture_buf <-Buff100_LULC$buff_pasture_hay2006_pct +
  Buff100_LULC$buff_rowcrops2006_pct

Buff100_LULC$wetland_buf<-Buff100_LULC$buff_wetland_woody2006_pct +
  Buff100_LULC$buff_wetland_emergent2006_pct

#change road density to km/sqkm
Buff100_LULC$roaddens_buf<-Buff100_LULC$buff_roaddensity_mperha / 10

buff100<-dplyr::select(Buff100_LULC, lagoslakeid, urban_buf, forest_buf, agriculture_buf, wetland_buf, roaddens_buf)


# watershed area
lake_ws_area <- data.frame(lagoslakeid=lg$iws$lagoslakeid,
                           iwsarea_ha=lg$iws$iws_ha)
lake_ws_area$area_ws <- lake_ws_area$iwsarea_ha / 100

area_table<-dplyr::select(lake_ws_area, lagoslakeid, area_ws )

#locus table has lat/lon, lake area, lake perimiter, lake elevation, all of the zones that lake is a part of
lake_info <- data.frame(lagoslakeid=lg$locus$lagoslakeid,
                        hu8_id=lg$locus$hu8_zoneid,
                        hu6_id=lg$locus$hu6_zoneid,
                        hu4_id=lg$locus$hu4_zoneid,
                        elev_ws=lg$locus$elevation_m)

#merge all ecological context data by lagoslakeid should be #51,065 observations of lakes > 4ha
all_vars <- left_join( area_table, lake_info, by = "lagoslakeid") %>%
  left_join(buff100, by = "lagoslakeid") %>%
  left_join(IWS_table, by = "lagoslakeid")
all_vars$lagoslakeid<-as.character(all_vars$lagoslakeid)

### join predictors and fish data 
MI_data_with_pred<-left_join(MI_data_LMB, all_vars)


#### link to historical lakes #### 
# to see if there is more contemp Secchi data 
#read in lat/lon
points<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/IFR_Lake_Points.csv")%>%
  rename(new_key = NEW_KEY) %>% 
  dplyr::select(new_key, LONG_DD, LAT_DD)
#read in historical data
hist_dat<-read.csv("Data/historical_lmb2.csv") %>% 
  dplyr::select(new_key) %>% 
  distinct(new_key) %>% 
  left_join(points)
#read in SNT secchi data 
snt_secchi<-read.csv("Data/MI_data/snt_secchi.csv") %>%
  dplyr:: select(Survey_Number, New_Key, SECCHI_FT)%>% 
  rename(survey_number=Survey_Number, new_key=New_Key) %>% #good new_key 
  mutate(secchi_m = SECCHI_FT/3.28) #convert to meters
#lagos shapefile 
lakes_shape<- readOGR(dsn = "/Users/katelynking/Desktop/LAGOS_US/LAGOS_NE_All_Lakes_4ha", layer = "LAGOS_NE_All_Lakes_4ha")
crs(lakes_shape)
#get the lat long data for the Michigan lakes with historical data 
x <- hist_dat$LONG_DD
y <- hist_dat$LAT_DD

#make a dataframe of the coordinates and project
lake.ll <- SpatialPoints(data.frame(x=x,y=y), proj4string=CRS("+proj=longlat +datum=NAD83"))   

prjnew <- CRS("+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs ")

MI.aea <- spTransform(lake.ll, prjnew)
plot(MI.aea)

#join with the 4ha shapefile from LAGOS_NE
MI_LAGOS<-sp::over(MI.aea, lakes_shape) %>%
  dplyr::select(lagoslakei, STATE, COUNTY_Nam, Lake_Area_) %>%
  rename(lagoslakeid = lagoslakei, state = STATE, county= COUNTY_Nam, area_ha = Lake_Area_) 
MI_LAGOS$new_key<-paste(hist_dat$new_key)  #only 6 lakes don't match 

####get LAGOS Secchi data #### 
#from NE - dataset is version 1.087.3. downloaded Aug 2022
library(LAGOSNE)
#lagosne_get(dest_folder = lagos_path()) run one time to get files on computer
lagos <- lagosne_load()
names(lagos)
query_lagos_names("secchi")
epi_nutr <- lagos$epi_nutr
secchi <- epi_nutr %>% 
  dplyr::select(lagoslakeid, sampledate, secchi) %>% 
  tidyr::drop_na(secchi)

# time period: >=2002, using samples collected between Jun 15 and Sep 15 (stratification period)
# define time parameters
first_year <- 2002
first_day  <- '0615' #i.e., '0615' for Jun 15
last_day   <- '0915'

# convert sampledate to date format and create "monthday" column to filter by day
secchi$sampledate   <- as.Date(secchi$sampledate, format="%m/%d/%Y")
secchi$monthday     <- format(secchi$sampledate, format="%m%d")
secchi$sampleyear     <- format(secchi$sampledate, format="%Y")

secchi_subset <- secchi %>% 
  filter(monthday >= first_day & monthday <= last_day) %>% 
  filter(sampleyear >= first_year) %>% 
  group_by(lagoslakeid) %>% 
  summarise(mean_secchi = round(mean(secchi), 3)) %>% 
  mutate(lagoslakeid = as.character(as.integer(lagoslakeid))
         )

hist_lakes_lagos_secchi<-left_join(MI_LAGOS, secchi_subset) %>% 
  left_join(snt_secchi) %>% 
  mutate(secchi_combo= ifelse(is.na(secchi_m), mean_secchi, secchi_m), #use contemp data, if na, use lagos 
         secchi_combo = round(secchi_combo,3)) %>% 
  left_join(hist_dat) %>% 
  tidyr::drop_na(secchi_combo) %>% 
  dplyr::select(new_key, secchi_combo)

#add data to the historical dataset 
historical_lmb2<-read.csv("Data/historical_lmb2.csv")
hist_dat_with_secchi<-left_join(historical_lmb2, hist_lakes_lagos_secchi)

#save the data 
write.csv(hist_dat_with_secchi, "Data/historical_lmb2_with_secchi.csv", row.names = FALSE)

#get MI basemap 
MI_basemap<-map_data("state",  region = c("michigan"))  # select michigan 
map<-ggplot(data = MI_basemap) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + #this fill MI white and outline is black
  coord_fixed(1.3) 

map + 
  geom_point(data=hist_lakes_lagos_secchi, aes(x = LONG_DD, y = LAT_DD, colour=c(secchi_combo) )) +
  theme_bw() +  
  theme( panel.grid.major = element_blank(),
         panel.grid.minor = element_blank()) +
  scale_colour_gradient(low = "darkgreen", high="blue") +
  theme(legend.position = c(0.20, 0.25), 
        legend.text = element_text(size=12)) 
