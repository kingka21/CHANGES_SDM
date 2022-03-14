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
library(tidyr)

#read in lagos lake link and see if that match works
lagos_lake_link<-read.csv("/Users/katelynking/Desktop/LAGOS_US/lake_link.csv") %>%
          dplyr::select(lagoslakeid, lake_nhdid) 
lagos_lakes<- lagos_lake_link[!duplicated(paste(lagos_lake_link$lake_nhdid)),]
lake_link<-read.csv("Data/MI_data/new_key_comid_translate.csv")
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

########################################
### try matching to Winslow lakes ###########
##################################################
winslow<-read.table("/Users/katelynking/Desktop/Winslow data/NLDAS_thermal_metrics.tsv", sep = '\t', header = TRUE) %>%
  select(year, site_id, gdd_wtr_0c, gdd_wtr_5c, gdd_wtr_10c)
dd_temp<-read.csv("Data/MI_data/lake_degree_days_year.csv") %>% 
  rename(dnr_nhdid=IHDLKID) %>% 
  select(-c(FETCH_M, ZMAX_M, ZAVE_M, D)) %>%
  pivot_longer(cols = starts_with("DD"),
               names_to = "year",
               names_prefix = "DD_",
               values_to = "ggd_dnr")
dd_temp$dnr_nhdid<-as.character(as.integer(dd_temp$dnr_nhdid))
dd_temp$year<-as.integer(as.character(dd_temp$year))

#dnr_winslow<-left_join(lake_link, winslow, by=c("nhdid" = "site_id")) #nothing matches :( 

winslow_lagos<-read.csv("/Users/katelynking/Desktop/Winslow data/Winslow_LAGOS_Xwalk.csv")
winslow_lagos$winslow_nhdid<- gsub("nhd_", "", winslow_lagos$site_id)
lagos_lakes<- lagos_lake_link[!duplicated(paste(lagos_lake_link$lake_nhdid)),]

lagos_cross<-left_join(winslow_lagos, lagos_lakes) %>% 
  select(lagoslakeid, winslow_nhdid, lake_nhdid) #here nhdid is from Winslow and lake_nhdid is from lagos crosswalk

lake_link$dnr_nhdid<-as.character(as.integer(lake_link$nhdid))
winslow$site_id<-as.character(as.integer(winslow$site_id))

dnr_winslow<-left_join(lake_link, lagos_cross, by=c("dnr_nhdid"="lake_nhdid")) %>% #1706 to match 
        left_join(winslow, by=c("winslow_nhdid" = "site_id")) %>%
        drop_na(winslow_nhdid) %>% 
  left_join(dd_temp, by=c("dnr_nhdid","year")) %>% 
  group_by(dnr_nhdid, year) %>% 
  filter(n() == 1) #get rid of duplicates 
plot(dnr_winslow$gdd_wtr_0c, dnr_winslow$ggd_dnr)
boxplot(dnr_winslow$gdd_wtr_0c, dnr_winslow$ggd_dnr)
t.test(dnr_winslow$gdd_wtr_0c, dnr_winslow$ggd_dnr)

#look at duplicates in the temp variables
dups<-dnr_winslow %>% 
      group_by(dnr_nhdid, year) %>% 
      filter(n() > 1)
