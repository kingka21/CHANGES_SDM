# code written by Katelyn King 
# created: April 2021 
# Data wrangling 

#### load libraries ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

##### data merging and cleaning ####
#*driver variables ####
lake_vars<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/lake_attributes.csv") %>%
  dplyr::select(IHDLKID,  ACRES_IHD, PERIM_KM, STRAHLER, LINK, UP_LAKES, UP_LAKE_ACRES, DOWN_LAKES, DOWN_LAKE_ACRES, WL_KM2 ) %>%
  rename(ihdlkid = IHDLKID, lake_area_ac=ACRES_IHD, lake_perim_km=PERIM_KM, lake_order=STRAHLER, lake_link_num = LINK, uplakes_n=UP_LAKES, uplakes_area_ac=UP_LAKE_ACRES, downlakes_n=DOWN_LAKES, downlakes_area_ac=DOWN_LAKE_ACRES, ws_area_km2=WL_KM2 )
lake_link<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/new_key_comid_translate.csv") %>% #the lake variables table does not have many new_key ids - so link to by this table from KEvin 
            dplyr::select(COMID, NEW_KEY) %>%
            rename(new_key=NEW_KEY)
lake_depth<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/lake_depth_2021.csv") %>%
  dplyr::select(new_key, MaxDepth_ft, MeanDepth_ft)
lake_join<-left_join(lake_vars, lake_link, by=c('ihdlkid' = "COMID")) %>% #some comids have multiple new_keys, keep duplicates 
            left_join(lake_depth)
ws_vars<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/local_catchment_attributes.csv") %>%
  dplyr::select(IHDID, WL_NLCD21, WL_NLCD22, WL_NLCD23, WL_NLCD24, WL_NLCD41, WL_NLCD42, WL_NLCD43, WL_NLCD52, WL_NLCD71, WL_NLCD81, WL_NLCD82, WL_NLCD90, WL_NLCD95, WL_SLOPE, WL_ELEVMEAN ) %>%
  rename(ihdlkid = IHDID,  ws_slope_deg=WL_SLOPE, ws_mean_elevation_m=WL_ELEVMEAN)
#add up all the developed categories and call it urban
ws_vars$ws_urban_prop <- ws_vars$WL_NLCD21 + ws_vars$WL_NLCD22 + ws_vars$WL_NLCD23 + ws_vars$WL_NLCD24
ws_vars$ws_forest_prop <- ws_vars$WL_NLCD41 + ws_vars$WL_NLCD42 + ws_vars$WL_NLCD43
ws_vars$ws_agriculture_prop <- ws_vars$WL_NLCD81 + ws_vars$WL_NLCD82
ws_vars$ws_wetland_prop <- ws_vars$WL_NLCD90 + ws_vars$WL_NLCD95
ws_vars$ws_shrub_prop<-ws_vars$WL_NLCD52 + ws_vars$WL_NLCD71
ws_vars2<-dplyr::select(ws_vars, ihdlkid, ws_urban_prop, ws_forest_prop, ws_agriculture_prop, ws_shrub_prop, ws_wetland_prop, ws_slope_deg, ws_mean_elevation_m)
secchi<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/snt_secchi.csv") %>%
          dplyr:: select(Survey_Number, New_Key, SECCHI_FT)%>% 
          rename(survey_number=Survey_Number, new_key=New_Key)
shore_vars<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/snt_shoreline.csv") %>% 
          rename(survey_number=Survey_Number, new_key=NEW_KEY)
surface_temp<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/lake_surface_temp.csv") %>%
          dplyr:: select(IHDLKID, TAVE_2002:TAVE_2019) 
surface_temp$surface_temp_mean<-rowMeans(surface_temp[,2:19])# take the means across the current years (2002-2019)
surface_mean<-dplyr::select(surface_temp, IHDLKID, surface_temp_mean) %>% 
                rename(ihdlkid=IHDLKID)
dd_temp<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/lake_degree_days_year.csv") %>%
          dplyr:: select(IHDLKID, DD_2002:DD_2019)
dd_temp$dd_mean<-rowMeans(dd_temp[,2:19])# take the means across the current years (2002-2019)
dd_mean<-dplyr::select(dd_temp, IHDLKID, dd_mean)%>% 
          rename(ihdlkid=IHDLKID)


#look at duplicates in the temp variables #randomly remove? remove longer or shorter? 
dups<-dd_mean %>% 
  group_by(ihdlkid) %>% 
  filter(n() > 1)
dups<-surface_mean %>% 
  group_by(ihdlkid) %>% 
  filter(n() > 1)

surface_temp_mean<- surface_mean[!duplicated(paste(surface_mean$ihdlkid)),] #select only one row of each lake (randomly remove)
degdays_mean<- dd_mean[!duplicated(paste(dd_mean$ihdlkid)),] #select only one row of each lake

#join driver variables 
driver_vars<-left_join(lake_join, ws_vars2) %>%
  left_join(surface_temp_mean) %>%
  left_join(degdays_mean) %>%
  left_join(shore_vars, by=("new_key")) %>% # will be more obs because multiple surveys 
  left_join(secchi, by=c("new_key"))
#put all column names to lowercase for consistency  
names(driver_vars)<-tolower(names(driver_vars))  
#re-order columns 
driver_vars <- driver_vars[, c(1, 11, 23, 30, 2, 3, 4, 5, 6, 7,8,9,10,12,13,14,15,16,17,18,19,20,21,22,24,25,26,27,28,29,31)]

#duplicates come from multiple surveys e.g. shoreline and secchi 
drivers<- driver_vars[!duplicated(paste(driver_vars$new_key)),] #select only one row of each lake - need to do this by new_key because some nhdids have multiple new_keys

#set lake order to 0 if it is an isolated lake (currently -99)
drivers <- drivers %>%
  mutate(
    lake_order = case_when(
      lake_order == -99 ~ as.integer(0),   #if ~ then 
      TRUE      ~ drivers$lake_order))  ### the else part 

#convert secchi to meters and lake area to m2 and depth to m
drivers$secchi_m<-drivers$secchi_ft/3.28
drivers$lake_area_m2<-drivers$lake_area_ac/0.00024711 #change from acre to m2
drivers$maxdepth_m<-drivers$maxdepth_ft/3.28
drivers$meandepth_ft<-drivers$meandepth_ft/3.28

#make the variable geometry ratio - Area (m2) ^ .25 / max depth (m)
drivers$geom_ratio<-(drivers$lake_area_m2^0.25) / drivers$maxdepth_m

write.csv(drivers, "Data/driver_varibles.csv", row.names = FALSE)


############################
##### catch data #### 
############################
lake_info<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/snt_lake_info.csv")%>%
  rename(new_key=New_Key, survey_number=Survey_Number)
catch<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/snt_catch_data.csv") %>%
  dplyr:: select(NEW_KEY, Survey_Number, GEAR, SPECIES, FISH_COUNT)%>%
  rename(new_key=NEW_KEY, survey_number=Survey_Number)
effort<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/snt_effort_data.csv") %>%
  dplyr::select(NEW_KEY, Survey_Number, GEAR, EFFORT_MEASURE, EFFORT, SAMPLE_START_DATE, SAMPLE_END_DATE) %>%
  rename(survey_number=Survey_Number, new_key=NEW_KEY)

#join raw data tables 
catch_data<-left_join(lake_info, effort, by=c("new_key","survey_number")) %>%#join by survey number
  full_join(catch, by = c("new_key", "survey_number", "GEAR")) %>%  # join by survey number and gear to get effort #full join keeps all of the gears/efforts used 
  drop_na('GEAR') #drop lakes where we don't have sample data 

catch_data$SPECIES<-ifelse(is.na(catch_data$SPECIES), 'NONE', catch_data$SPECIES) #keep effort even when no sp caught 
catch_data$FISH_COUNT<-ifelse(is.na(catch_data$FISH_COUNT), 0, catch_data$FISH_COUNT)

#put all column names to lowercase for consistency  
names(catch_data)<-tolower(names(catch_data))

#re-order columns 
#catch_data <- catch_data[, c(1, 2, 19, 6, 7, 8, 14, 11, 12, 13, 17, 18, 9, 15, 16, 22, 23, 4, 5, 3, 20, 21, 10)]

#convert date to julian (year) day 
catch_data$date<- as.Date(catch_data$sample_start_date, "%m/%d/%y") # convert "date" from chr to a Date class and specify current date format
catch_data$julian <- yday(catch_data$date)  

#MI_data$macth<-MI_data$survey_year == MI_data$year #test years, all TRUE, remove SURVEY_YEAR  
#MI_data$macth<-MI_data$new_key.x == MI_data$new_key.y #test to see if the NEW_KEY matches, all TRUE, so can drop NEW_KEY from effort dataset

#*look for and remove duplicate surveys on lakes ####
# distinct lake/survey combos n=505
#unique_LS<-catch_data %>% 
 # distinct(new_key, survey_number, .keep_all = T)

# look at non-distinct: only 33 lakes with multiple surveys
#lakes_multi_survey<-unique_LS %>% 
 # group_by(new_key) %>% 
  #filter(n() > 1)

#choose the most recent survey 
MI_data_no_dups<-catch_data %>%
  group_by(new_key) %>%
  slice_max(year)

#problem is there are still 3 lakes with two basins, so two surveys but one new_key 
  #unique_LS<-MI_data_no_dups %>% 
    # distinct(new_key, survey_number, .keep_all = T)

# non-distinct only 3 lakes with multiple surveys
  #lakes_multi_survey<-unique_LS %>% 
    #group_by(new_key) %>% 
      #filter(n() > 1)

#combine gear types for effort 
MI_data_no_dups$gear2<-sapply( MI_data_no_dups$gear, function(x) {
  if(x == "LMFYKE") {'FT_NET'}
  else {
    if(x == 'SMFYKE') {'FT_NET'}
    else {
      if(x == 'GLGNET') {'GILL'}
      else {
        if(x== 'IGNET') {'GILL'}
        else {
          if(x== 'SEINE') {'SEINE'}
          else {
            if(x== 'TRAPNET') {'FT_NET'}
            else {
              if(x== 'BOOMSHK') {'SHOCK'}
                    }}}}}}})

# group by lake, species, and gear then sum counts and effort 
count_data<-dplyr::select(MI_data_no_dups, new_key, species, fish_count, gear2, effort) %>%
    group_by(new_key, species, gear2) %>%
    summarise_if(is.numeric, sum)

#calculate CPUE ->abundance divided by effort 
count_data$cpue<-round(count_data$fish_count/count_data$effort, 3) 

#need to "spread" the table so that lakes with missing sp can be included as 0s in the dist model 
#rearrange table to have species along the top, will aggregate by the lake name column (#472 lakes) and gear 
lake_sp_wide<-tidyr::pivot_wider(data= count_data, 
                          id_cols = c(new_key, gear2),
                           names_from = species, 
                           values_from = cpue, 
                           values_fill = list(cpue = 0)) #replace NAs with 0
                         #  values_fn = list(cpue = sum)) #adds up all of the catch data

#### LARGE MOUTH BASS DATA #### 
#select out LMB # this keeps effort for absences
MI_data_LMB<- select(lake_sp_wide, new_key, gear2, LMB) 

names(MI_data_LMB)<-tolower(names(MI_data_LMB))


lmb_format<-tidyr::pivot_wider(data= MI_data_LMB, 
                                  id_cols = c(new_key),
                                  names_from = gear2, 
                                  values_from = lmb, 
                                  values_fill = list(lmb = 0))  #replace NAs with 0
                              


#use catch data not cpue #this does not have lakes with no LMB 
count_lmb<-filter(count_data, species=='LMB')
count_lmb_format<-tidyr::pivot_wider(data= count_lmb, 
                               id_cols = c(new_key),
                               names_from = gear2, 
                               values_from = c(fish_count, effort),
                               values_fill = list(fish_count = 0, effort = 0)) #replace NAs with 0

#### WALLEYE DATA #### 
#select out WAE # this keeps effort for absences
MI_data_WAE<- select(lake_sp_cpue, new_key, gear2, WAE) 

names(MI_data_WAE)<-tolower(names(MI_data_WAE))


wae_format<-tidyr::pivot_wider(data= MI_data_WAE, 
                               id_cols = c(new_key),
                               names_from = gear2, 
                               values_from = wae, 
                               values_fill = list(wae = 0))  #replace NAs with 0
          


#use catch data not cpue 
count_wae<-filter(count_data, species=='WAE')
count_wae_format<-tidyr::pivot_wider(data= count_wae, 
                                     id_cols = c(new_key),
                                     names_from = gear2, 
                                     values_from = c(fish_count, effort),
                                     values_fill = list(fish_count = 0, effort = 0)) #replace NAs with 0

#### CISCO DATA #### 
#select out WAE # this keeps effort for absences
MI_data_CIS<- select(lake_sp_cpue, new_key, gear2, CIS) 

names(MI_data_CIS)<-tolower(names(MI_data_CIS))


cis_format<-tidyr::pivot_wider(data= MI_data_CIS, 
                               id_cols = c(new_key),
                               names_from = gear2, 
                               values_from = cis, 
                               values_fill = list(cis = 0))  #replace NAs with 0
       


#use catch data not cpue 
count_cis<-filter(count_data, species=='CIS')
count_cis_format<-tidyr::pivot_wider(data= count_cis, 
                                     id_cols = c(new_key),
                                     names_from = gear2, 
                                     values_from = c(fish_count, effort),
                                     values_fill = list(fish_count = 0, effort = 0)) #replace NAs with 0

#### link catch data to drivers ####
drivers<-read.csv( "Data/driver_varibles.csv")

#pull out lat/lon for spatial autocorrelation testing, mapping, modeling 
lake_ll<- lake_info[!duplicated(paste(lake_info$new_key)),] 

#join tables for LMB 
dat<-left_join(lmb_format, drivers ) %>% # 472 unique lakes, but 41 lakes don't match the drivers with the nhdid
      left_join(dplyr::select(lake_ll, new_key, LONG_DD, LAT_DD, FMU_Code))

lmb_count_dat<-left_join(count_lmb_format, dat)  %>% #count data don't have lakes with 0 obs so less rows
      rename(SHOCK_cpue = SHOCK, FT_NET_cpue = FT_NET, GILL_cpue = GILL , SEINE_cpue = SEINE)

#join tables for WAE
dat<-left_join(wae_format, drivers ) %>% # 472 unique lakes, but 41 lakes don't match the drivers with the nhdid
  left_join(dplyr::select(lake_ll, new_key, LONG_DD, LAT_DD, FMU_Code))

wae_count_dat<-left_join(count_wae_format, dat) #count data don't have lakes with 0 obs so less rows

#join tables for CIS 
dat<-left_join(cis_format, drivers ) %>% # 472 unique lakes, but 41 lakes don't match the drivers with the nhdid
  left_join(dplyr::select(lake_ll, new_key, LONG_DD, LAT_DD, FMU_Code))

cis_count_dat<-left_join(count_cis_format, dat) #count data don't have lakes with 0 obs so less rows

#final drivers of Hansen et al 2017 
#walleye- degree days, lake area, conductivity, and shoreline complexity 
#lmb - degree days, lake order, and Secchi depth - 3/4 of our obs have secchi - could pull secchi from LAGOS 

##### quick check across gear types ####
hist(log(MI_data_LMB$lmb))
anova(lm(log(lmb + 0.001) ~ gear2, data = MI_data_LMB, na.action = na.omit ))
library(agricolae) #LSD test
LSD<-LSD.test(log(MI_data_LMB$lmb + 0.001), MI_data_LMB$gear2, 814, 8.9171, alpha=0.05)  #specify the DF and MSE of the residuals
LSD$groups

#try without the 0s
dat_no0<-filter(MI_data_LMB, lmb >0) 
hist(log(dat_no0$lmb))
anova(lm(log(lmb) ~ gear2, data = dat_no0, na.action = na.omit ))
library(agricolae) #LSD test
LSD<-LSD.test(log(dat_no0$lmb), dat_no0$gear2, 645, 1.559, alpha=0.05)  #specify the DF and MSE of the residuals
LSD$groups
