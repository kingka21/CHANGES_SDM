#### Lake matching between Zooniverse cards and MI lakes ####
#Public Land Survey System (PLSS) uses  Township Range and Section (TRS)
#### load libraries ####
library(rgdal)
library(sp)
library(sf)
library(dplyr)
library(raster)
library(ggplot2)
library(mapview)
library(rgeos)

#read in lagos lake link and see if that match works
lagos_lake_link<-read.csv("/Users/katelynking/Desktop/Cont Limno/lake_link_Jan2021.csv") %>%
          dplyr::select(lagoslakeid, lake_nhdid) 
lagos_lakes<- lagos_lake_link[!duplicated(paste(lagos_lake_link$lake_nhdid)),]
lake_link<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/new_key_comid_translate.csv")
join<-left_join(lake_link, lagos_lakes, by=c('COMID' = "lake_nhdid"))
join<-join[!is.na(join$lagoslakei),] #only 7751 - not much more than the join with the lake polys 

#read in polygons 
#MI assigned new_key to lake >5 acres (2ha) and is ~9,100 lakes 
##>1ha lakes LAGOS lakes 15,886
lakes_shape<- readOGR(dsn = "MI_poly_shapefile/", layer = "MI_poly_1ha")
crs(lakes_shape) #shows me the projection (aea projection)

#try to match new_key by using the link table
lake_link<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/NEW_KEY_NHD_link.csv") #Jim
lake_link<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/new_key_comid_translate.csv") #Kevin
lake_link<-filter(lake_link, NEW_KEY != "") #9116 obs with new key
lake_link$COMID<-as.factor(lake_link$COMID)
poly_data<-sf::st_as_sf(lakes_shape) #get data table from lake polygons 
lake_join<-merge(poly_data, lake_link, by.x='Permanent_', by.y= 'COMID' )
lake_join<-lake_join[!is.na(lake_join$lagoslakei),] #7,969 lakes matched based on NHD ID using Jim's table but only 7730 match from Kevin's 
MI_poly.sp<-as(lake_join, "Spatial") # put polygon back to spatial 

#TRS - Township Range and Section
#pull gis layer - had to do this on a PC because GDAL is not supported by MAC 
#devtools::install_github("yonghah/esri2sf")
#library("esri2sf")
#url <- "https://arcgis-server.lsa.umich.edu/arcgis/rest/services/IFR/mi_boundaries/MapServer/3"
#df <- esri2sf(url)
#plot(df)
#sf::st_write(df, "TRS.shp")
#read in layer 
TRS_shape<- readOGR(dsn = "/Users/katelynking/Desktop/UofM/CHANGES/lake_matching/TRS_shape", layer = "TRS")
crs(TRS_shape) #shows me the projection 
plot(TRS_shape)
TRS.aea <- spTransform(TRS_shape, crs(MI_poly.sp)) #transform to aea projection 

#overlay the lakes on TRS to see matches 
overlay<-over(MI_poly.sp, TRS.aea ) # keeps only attributes of TRS 
overlay$nhd_comid<-paste(MI_poly.sp@data$Permanent_) #add back NHD ids
overlay$lagoslakeid<-paste(MI_poly.sp@data$lagoslakei )  # add back lagos ids
overlay$NEW_KEY<-paste(MI_poly.sp@data$NEW_KEY) #add in new key
overlay$Hectares<-paste(MI_poly.sp@data$Hectares)

#try matching with example mismatches from Zooniverse 
no_match<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/lake_matching/no_match.csv", header=TRUE)
test<-left_join(no_match, overlay, by=c("TRS" = "TWNRNGS"))

humph<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/lake_matching/HUMPHRIES_COUNTY_LATLONG_IFR_LAKES.csv", header=TRUE)

#### try adding a buffer around lakes ####
#add a buffer around lakes using rgeos
lakes_buf<-gBuffer(MI_poly, byid=TRUE, width=50) #must be sp object, byid buffers all geometries, 



