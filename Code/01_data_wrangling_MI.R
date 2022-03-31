##### data merging and cleaning ####
# code written by Katelyn King 
# created: April 2021 
# Data wrangling 

#### load libraries ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(rLakeAnalyzer)

#### morphometry and watershed driver variables ####
lake_vars<-read.csv("Data/MI_data/lake_attributes.csv") %>%
  dplyr::select(IHDLKID,  ACRES_IHD, PERIM_KM, STRAHLER, LINK, UP_LAKES, UP_LAKE_ACRES, DOWN_LAKES, DOWN_LAKE_ACRES, WL_KM2 ) %>%
  rename(nhdid = IHDLKID, lake_area_ac=ACRES_IHD, lake_perim_km=PERIM_KM, lake_order=STRAHLER, lake_link_num = LINK, uplakes_n=UP_LAKES, uplakes_area_ac=UP_LAKE_ACRES, downlakes_n=DOWN_LAKES, downlakes_area_ac=DOWN_LAKE_ACRES, ws_area_km2=WL_KM2 ) %>%
  mutate(lake_area_km2=round(lake_area_ac/247,2), #change from acre to km2
         lake_area_m2= lake_area_ac/0.00024711, #change from acre to m2
         lake_order = case_when( #set lake order to 0 if it is an isolated lake (currently -99)
           lake_order == -99 ~ as.integer(0),   #if ~ then 
           TRUE      ~ lake_order))  ### the else part   

lake_link<-read.csv("Data/MI_data/new_key_comid_translate.csv") #lake var table does not have many new_key- so link to by this table from KEvin (includes ~25 lakes that were updated with nhdids)
            
lake_depth<-read.csv("Data/MI_data/lake_depth.csv") %>%
  dplyr::select(New_Key, MaxDepth_ft, MeanDepth_ft) %>% 
  rename(new_key=New_Key) %>%
  mutate(maxdepth_m = round(MaxDepth_ft/3.28, 2), 
         meandepth_m= round(MeanDepth_ft/3.28 ,2))      #change from feet to meters

lake_join<-left_join(lake_vars, lake_link) %>% #some nhdids have multiple new_keys, keep duplicates 
            left_join(lake_depth) %>%
          mutate(geom_ratio = lake_area_m2^0.25 / maxdepth_m)  #make the variable geometry ratio - Area (m2) ^ .25 / max depth (m)

#watershed variables 
ws_vars<-read.csv("Data/MI_data/local_catchment_attributes.csv") %>%
  dplyr::select(IHDID, WL_NLCD21, WL_NLCD22, WL_NLCD23, WL_NLCD24, WL_NLCD41, WL_NLCD42, WL_NLCD43, WL_NLCD52, WL_NLCD71, WL_NLCD81, WL_NLCD82, WL_NLCD90, WL_NLCD95, WL_SLOPE, WL_ELEVMEAN ) %>%
  rename(nhdid = IHDID,  ws_slope_deg=WL_SLOPE, ws_mean_elevation_m=WL_ELEVMEAN) %>% 
  mutate(ws_urban_prop = WL_NLCD21 + WL_NLCD22 + WL_NLCD23 + WL_NLCD24, #add up categories 
         ws_forest_prop= WL_NLCD41 + WL_NLCD42 + WL_NLCD43,
         ws_agriculture_prop=WL_NLCD81 + WL_NLCD82,
         ws_wetland_prop= WL_NLCD90 + WL_NLCD95,
         ws_shrub_prop=WL_NLCD52 + WL_NLCD71) %>%
  dplyr::select(nhdid, ws_urban_prop, ws_forest_prop, ws_agriculture_prop, ws_shrub_prop, ws_wetland_prop, ws_slope_deg, ws_mean_elevation_m)


#### modeled temperature variables ####
surface_temp<-read.csv("Data/MI_data/lake_surface_temp.csv") %>%
          dplyr:: select(IHDLKID, HECTARES, TAVE_2002:TAVE_2019) #keep the area, this will help with matching to lakes in DNR, some NHD lakes are split into multiple DNR lakes
surface_temp$surface_temp_mean<-rowMeans(surface_temp[,3:20])# take the means across the current years (2002-2019)
surface_mean<-dplyr::select(surface_temp, IHDLKID, HECTARES, surface_temp_mean) %>% 
                rename(nhdid=IHDLKID)

dd_temp<-read.csv("Data/MI_data/lake_degree_days_year.csv") %>%
          dplyr:: select(IHDLKID, HECTARES, DD_2002:DD_2019)
dd_temp$dd_mean<-rowMeans(dd_temp[,3:20])# take the means across the current years (2002-2019)
dd_mean<-dplyr::select(dd_temp, IHDLKID, HECTARES, dd_mean)%>% 
          rename(nhdid=IHDLKID)

#combine temperature by id and area # no dups if you include different lake sizes! 
temp_dat<-left_join(dd_mean, surface_mean, by=c("nhdid", "HECTARES")) %>% 
  rename(lake_area_ha = HECTARES) %>% 
  group_by(nhdid) %>%
  slice_max(lake_area_ha) #use the larger area for dups, usually separated basins

#look at duplicates in the temp variables, multiple temp measures on diff polys based on area (selected largest above)
#dups<-temp_dat %>% 
 # group_by(nhdid, lake_area_ha) %>% 
  #filter(n() > 1)

#join driver variables 
driver_vars<-left_join(lake_join, ws_vars) %>%
  left_join(temp_dat, by=c('nhdid')) 
#put all column names to lowercase for consistency  
names(driver_vars)<-tolower(names(driver_vars))  
#re-order columns 
#driver_vars <- driver_vars[, c(1, 13, 2, 3, 4, 5, 6, 7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28)]

#no dups! 
#dups<-driver_vars %>% 
 # drop_na(new_key) %>% 
  #group_by(new_key) %>% 
  #filter(n() > 1)

#write.csv(driver_vars, "Data/driver_varibles.csv", row.names = FALSE)

#### secchi and shoreline driver variables ####
#there are duplicates of these variables because of multiple surveys 
secchi<-read.csv("Data/MI_data/snt_secchi.csv") %>%
  dplyr:: select(Survey_Number, New_Key, SECCHI_FT)%>% 
  rename(survey_number=Survey_Number, new_key=New_Key) %>% #good new_key 
  mutate(secchi_m = SECCHI_FT/3.28) #convert to meters

shore_vars<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/snt_shoreline.csv") %>% 
  rename(survey_number=Survey_Number, new_key=NEW_KEY) #good new_key

#### DO and TEMP PROFILE ####
temp_do<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/snt_tempdo_2002_2020.csv")
names(temp_do)<-tolower(names(temp_do))

#look for cases when depth is repeated (don't want this for profile calcs)
unique<-temp_do %>% 
  group_by(survey_number, survey_effort_key, basin_number, depth_new_ft) %>%
  mutate(n = n()) %>%
  filter(n > 1)

#remove some obs where there are different dates than the other measures 
temp_do<-temp_do[!(temp_do$survey_number =="3792" & temp_do$reading_datetime =="8/30/02"),]
temp_do<-temp_do[!(temp_do$survey_number =="5473" & temp_do$reading_datetime =="8/19/05"),]
#remove one of the two measures at the same depth - these are within a tenth of a degree, so removing one doesn't change results 
temp_do<-temp_do[!(temp_do$survey_number =="5509" & temp_do$depth_new_ft =="0" & temp_do$temp_new_f == "76.28"),]
temp_do<-temp_do[!(temp_do$survey_number =="5999" & temp_do$basin_number =="1" & temp_do$depth_new_ft =="27" & temp_do$temp_new_f == "50.04"),]
temp_do<-temp_do[!(temp_do$survey_number =="6890" & temp_do$depth_new_ft =="0" & temp_do$temp_new_f == "72.46"),]
temp_do<-temp_do[!(temp_do$survey_number =="11982" & temp_do$depth_new_ft =="56" & temp_do$temp_new_f == "40.6"),]

#convert ft to m and F to C 
temp_do$depth_m<-round(temp_do$depth_new_ft/3.281,3)
temp_do$temp_c<-round((temp_do$temp_new_f-32)*(5/9), 3)

#calc thermo depths using the Rlakeanalyzer package : need unique depths  
thermo_depth<-temp_do %>% 
              group_by(survey_number, survey_effort_key, basin_number) %>%
              summarise(thermo_depth_m = thermo.depth(wtr=temp_c, depths=depth_m, seasonal = FALSE)) #location of the thermocline from a temperature profile

#plot one of the lake surveys to see if the calculated depth matches visual 
test<-filter(temp_do, survey_number == 3415 & survey_effort_key == 56)
ggplot(test, aes(x = temp_c, y = depth_m, colour = temp_c)) +
  geom_point() + 
  scale_y_reverse() +
  scale_colour_gradient2(
    midpoint = 15, 
    high = ("red"), #scales::muted("red") will give you a muted color! 
    low = ("blue") #scales::muted("blue")
  ) + 
  theme_bw()

#pull out the top/surface value for each combo 
top<-temp_do %>% 
  group_by(survey_number, survey_effort_key, basin_number) %>%
  slice_head()  %>% #oull first value
  dplyr::select(survey_number, survey_effort_key, basin_number, reading_datetime, depth_m, temp_c, oxygen_new_mgl)%>% 
  rename(surface_depth_m  = depth_m, surface_temp_c= temp_c, surface_do_mgl = oxygen_new_mgl)


#pull out the bottom value for each combo 
bottom<-temp_do %>% 
  group_by(survey_number, survey_effort_key, basin_number) %>%
  slice_tail() %>% #pulls last value
  dplyr::select(survey_number, survey_effort_key, basin_number, reading_datetime, depth_m, temp_c, oxygen_new_mgl) %>% 
  rename(bottom_depth_m  = depth_m, bottom_temp_c= temp_c, bottom_do_mgl = oxygen_new_mgl)

#join new temp/do measures 
temp_do_measures<-left_join(top, bottom) %>%
  left_join(thermo_depth) %>%
  group_by(survey_number) %>%
  slice_max(bottom_depth_m)

### for all DO/temp measures, use the deepest basin in a lake. 
temp_do<-read.csv("Data/MI_data/temp_do_measures_MI.csv")
temp_do_no_dups<-temp_do %>%
  group_by(survey_number) %>%
  slice_max(bottom_depth_m)
## SAVE  AGAIN WITH NO DUPS BUT ADD TO FUNCTION ABOVE ON NEW COMP
#write.csv(temp_do_measures, "Data/temp_do_measures_MI.csv", row.names = FALSE)

#### TDO3 #####
#remove NAs to interpolate temperature at 3mg/L DO 
temp_do_noNA<-temp_do %>%
  drop_na('oxygen_new_mgl')
#surveys had 0 for all DO measures - check with Kevin to see if this is a mistake 
temp_do_noNA<-temp_do_noNA[!(temp_do_noNA$survey_number =="11294"),]
temp_do_noNA<-temp_do_noNA[!(temp_do_noNA$survey_number =="11381"),]
temp_do_noNA<-temp_do_noNA[!(temp_do_noNA$survey_number =="11384"),]

#use the approx function to interpolate values at 3mgl DO
data2<-temp_do_noNA %>%
    group_by(survey_number, survey_effort_key, basin_number)%>% 
    summarise(TDO3 = approx(oxygen_new_mgl,temp_c, xout=3)) #approx 
data2$ids<-rep(c('x', 'y'),times=468) #add an id column so can pivot 
TDO3 <- tidyr::pivot_wider(data = data2, 
                              id_cols = c(survey_number, survey_effort_key, basin_number),
                              names_from = ids, 
                              values_from = TDO3) %>%
                dplyr::select(-c(x)) %>% 
                rename(TDO3=y)

#the approx function outputs as a list, so need to unlist: 
TDO3[sapply(TDO3, is.list)] <- apply(TDO3[sapply(TDO3, is.list)], 
                                             1, function(x) 
                                               paste(unlist(x), 
                                                     sep=", ", collapse=", "))
#some values are NA: because the values don't go down to 3mgl and no extrapolation (called interpolation)
# maybe if I use lm and predict then I would get it - but maybe dont want to do this because NA would mean that DO doesn't get down to lethal in that lake 

#create a plot to check 
test<-filter(temp_do_noNA, survey_number == 3451 & survey_effort_key == 1)
plot (test$oxygen_new_mgl, test$temp_c, type="o")
abline(v = 3, col = "#ff0000") # x value 
abline(h = 23.9, col = 'blue') # y value 

#
#write.csv(TDO3, "Data/TDO3_MI.csv", row.names = FALSE)


############################
##### catch data #### 
############################
lake_info<-read.csv("Data/MI_data/snt_lake_info.csv")

catch<-read.csv("Data/MI_data/snt_catch_data.csv") %>%
  dplyr:: select(NEW_KEY, Survey_Number, GEAR, SPECIES, FISH_COUNT)%>%
  rename(new_key=NEW_KEY, survey_number=Survey_Number)

effort<-read.csv( "Data/MI_data/snt_effort_data.csv") %>%
  dplyr::select(NEW_KEY, Survey_Number, GEAR, EFFORT_MEASURE, EFFORT, SAMPLE_START_DATE, SAMPLE_END_DATE) %>%
  rename(survey_number=Survey_Number, new_key=NEW_KEY)

#join raw data tables 
catch_data<-left_join(lake_info, effort, by=c("new_key", "survey_number")) %>% #join by survey number
  full_join(catch, by = c("new_key", "survey_number", "GEAR")) %>%  # join by survey number and gear to get effort #full join keeps all of the gears/efforts used 
  drop_na('GEAR')  #drop lakes where we don't have sample data 
 
catch_data$SPECIES<-ifelse(is.na(catch_data$SPECIES), 'NONE', catch_data$SPECIES) #keep effort even when no sp caught 
catch_data$FISH_COUNT<-ifelse(is.na(catch_data$FISH_COUNT), 0, catch_data$FISH_COUNT)

#put all column names to lowercase for consistency  
names(catch_data)<-tolower(names(catch_data))

#re-order columns 
#catch_data <- catch_data[, c(1, 2, 19, 6, 7, 8, 14, 11, 12, 13, 17, 18, 9, 15, 16, 22, 23, 4, 5, 3, 20, 21, 10)]

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

#there are still 3 lakes with two basins, so two surveys but one new_key 
 # unique_LS<-MI_data_no_dups %>% 
    # distinct(new_key, survey_number, .keep_all = T)

# non-distinct only 3 lakes with multiple surveys, will keep these 
  #lakes_multi_survey<-unique_LS %>% 
    #group_by(new_key) %>% 
      #filter(n() > 1)


# group by lake, species, and gear then sum counts and effort 
count_data<-dplyr::select(MI_data_no_dups, new_key, fmu_code, species, fish_count, gear, effort) %>%
    group_by(new_key, species, gear) %>%
    summarise_if(is.numeric, sum) %>% 
  ungroup()


#### LARGE MOUTH BASS DATA #### 

#use catch data #this does not have lakes with no LMB 
count_lmb<-filter(count_data, species=='LMB')
#link to effort table to include gear and lakes that did not catch LMB or where LMB was not present 
all_effort<-left_join(effort, count_lmb, by=c('new_key' = 'new_key', 'GEAR'= 'gear')) %>%
  dplyr::select(-c(effort, species ))

#put all column names to lowercase for consistency  
names(all_effort)<-tolower(names(all_effort))

all_effort$fish_count<-ifelse(is.na(all_effort$fish_count), 0, all_effort$fish_count)

#convert date to a year and julian (year) day 
all_effort$date<- as.Date(all_effort$sample_start_date, "%m/%d/%y") # convert "date" from chr to a Date class and specify current date format
all_effort$year<-year(all_effort$date)
all_effort$julian <- yday(all_effort$date)  

#make a dataframe with the date info to join back 
sample_dates<-dplyr::select(all_effort, new_key, survey_number, date, year, julian, gear)%>%
  mutate(
  gear2 = case_when(       #combine gear types 
    gear == "LMFYKE" ~ 'FT_NET', 
    gear == "SMFYKE" ~ 'FT_NET', 
    gear == "TRAPNET" ~ 'FT_NET', 
    gear == "GLGNET" ~ 'GILL', 
    gear == "IGNET" ~ 'GILL', 
    gear == "SEINE" ~ 'SEINE', 
    gear == "BOOMSHK" ~ 'SHOCK')) %>%
  distinct(new_key, survey_number, gear2, .keep_all = TRUE) %>%
  dplyr::select(-c(gear))

#remove duplicate surveys, combine gears
lmb_new<-all_effort %>%
  group_by(new_key) %>%
  slice_max(year) %>%     #choose the most recent survey 
  mutate(
    gear2 = case_when(       #combine gear types 
    gear == "LMFYKE" ~ 'FT_NET', 
    gear == "SMFYKE" ~ 'FT_NET', 
    gear == "TRAPNET" ~ 'FT_NET', 
    gear == "GLGNET" ~ 'GILL', 
    gear == "IGNET" ~ 'GILL', 
    gear == "SEINE" ~ 'SEINE', 
    gear == "BOOMSHK" ~ 'SHOCK'))  %>%
  group_by(new_key, survey_number, gear2) %>% # group by lake, survey, and gear then sum counts and effort 
  summarize(effort_new = sum(effort),
            fish_count_new = sum(fish_count)) %>% #join 
  left_join(sample_dates)

#add predator information (NOP=northern pike and WAE = walleye)
pike<-filter(count_data, species=='NOP') %>% 
      select(new_key, gear) %>% 
        mutate(pike_pres = 1) %>%
  mutate(
    gear2 = case_when(       #combine gear types 
      gear == "LMFYKE" ~ 'FT_NET', 
      gear == "SMFYKE" ~ 'FT_NET', 
      gear == "TRAPNET" ~ 'FT_NET', 
      gear == "GLGNET" ~ 'GILL', 
      gear == "IGNET" ~ 'GILL', 
      gear == "SEINE" ~ 'SEINE', 
      gear == "BOOMSHK" ~ 'SHOCK'))  %>%
  group_by(new_key, gear2) %>% # group by lake and gear to combine gears 
  summarize(pike_pres2 = sum(pike_pres)) %>%
  mutate(pike_pres = 1) %>%
  select(-c(pike_pres2)) %>% 
  ungroup()

walleye<-filter(count_data, species=='WAE') %>% 
  select(new_key, gear) %>% 
  mutate(wae_pres = 1) %>%
  mutate(
    gear2 = case_when(       #combine gear types 
      gear == "LMFYKE" ~ 'FT_NET', 
      gear == "SMFYKE" ~ 'FT_NET', 
      gear == "TRAPNET" ~ 'FT_NET', 
      gear == "GLGNET" ~ 'GILL', 
      gear == "IGNET" ~ 'GILL', 
      gear == "SEINE" ~ 'SEINE', 
      gear == "BOOMSHK" ~ 'SHOCK'))  %>%
  group_by(new_key, gear2) %>% # group by lake and gear to combine gears 
  summarize(wae_pres2 = sum(wae_pres)) %>%
  mutate(wae_pres = 1) %>%
  select(-c(wae_pres2)) %>% 
  ungroup()


lake_ll<- lake_info[!duplicated(paste(lake_info$new_key)),] 

#join for modeling: 
lmb_dat_for_model<-left_join(lmb_new, secchi, by=c('new_key', "survey_number")) %>% 
  left_join(driver_vars, by=c('new_key')) %>%
  left_join(temp_do_no_dups) %>%
  left_join(pike, by=c('new_key', 'gear2')) %>% 
  left_join(walleye, by=c('new_key', 'gear2')) %>% 
  mutate(pike_pres = ifelse(is.na(pike_pres), 0, pike_pres),
           wae_pres= ifelse(is.na(wae_pres), 0, wae_pres)
           ) %>%
  left_join(dplyr::select(lake_ll, new_key, LONG_DD, LAT_DD, FMU_Code))
# left_join(temp_do, by = "survey_number") - remove multiple basin, need just deepest

write.csv(lmb_dat_for_model, "Data/lmb_dat_for_model_mar17.csv", row.names = FALSE)

#if you want one row per lake 
#count_lmb_format<-tidyr::pivot_wider(data= count_lmb, 
 #                                    id_cols = c(new_key),
  #                                   names_from = gear2, 
   #                                  values_from = c(fish_count, effort),
    #                                 values_fill = list(fish_count = 0, effort = 0)) #replace NAs with 0



#### WALLEYE DATA #### 
#use catch data #this does not have lakes with no WAE 
count_wae<-filter(count_data, species=='WAE')
#link to effort table to include gear and lakes that did not catch WAE or where WAE was not present 
all_effort<-left_join(effort, count_wae, by=c('new_key' = 'new_key', 'GEAR'= 'gear')) %>%
  dplyr::select(-c(effort, species ))

#put all column names to lowercase for consistency  
names(all_effort)<-tolower(names(all_effort))

all_effort$fish_count<-ifelse(is.na(all_effort$fish_count), 0, all_effort$fish_count)

#convert date to a year and julian (year) day 
all_effort$date<- as.Date(all_effort$sample_start_date, "%m/%d/%y") # convert "date" from chr to a Date class and specify current date format
all_effort$year<-year(all_effort$date)
all_effort$julian <- yday(all_effort$date)  

#make a dataframe with the date info to join back 
sample_dates<-dplyr::select(all_effort, new_key, survey_number, date, year, julian, gear)%>%
  mutate(
    gear2 = case_when(       #combine gear types 
      gear == "LMFYKE" ~ 'FT_NET', 
      gear == "SMFYKE" ~ 'FT_NET', 
      gear == "TRAPNET" ~ 'FT_NET', 
      gear == "GLGNET" ~ 'GILL', 
      gear == "IGNET" ~ 'GILL', 
      gear == "SEINE" ~ 'SEINE', 
      gear == "BOOMSHK" ~ 'SHOCK')) %>%
  distinct(new_key, survey_number, gear2, .keep_all = TRUE) %>%
  dplyr::select(-c(gear))

#remove duplicate surveys, combine gears
wae_new<-all_effort %>%
  group_by(new_key) %>%
  slice_max(year) %>%     #choose the most recent survey 
  mutate(
    gear2 = case_when(       #combine gear types 
      gear == "LMFYKE" ~ 'FT_NET', 
      gear == "SMFYKE" ~ 'FT_NET', 
      gear == "TRAPNET" ~ 'FT_NET', 
      gear == "GLGNET" ~ 'GILL', 
      gear == "IGNET" ~ 'GILL', 
      gear == "SEINE" ~ 'SEINE', 
      gear == "BOOMSHK" ~ 'SHOCK'))  %>%
  group_by(new_key, survey_number, gear2) %>% # group by lake, survey, and gear then sum counts and effort 
  summarize(effort_new = sum(effort),
            fish_count_new = sum(fish_count)) %>% #join 
  left_join(sample_dates)

#join for modeling: 
wae_dat_for_model<-left_join(wae_new, secchi, by=c('new_key', "survey_number")) %>% 
  left_join(drivers, by=c('new_key')) %>%
  left_join(temp_do_no_dups) %>%
  left_join(dplyr::select(lake_ll, new_key, LONG_DD, LAT_DD, FMU_Code))

#### CISCO DATA #### 
#select out WAE # this keeps effort for absences
#use catch data #this does not have lakes with no WAE 
count_cis<-filter(count_data, species=='CIS')
#link to effort table to include gear and lakes that did not catch WAE or where WAE was not present 
all_effort<-left_join(effort, count_cis, by=c('new_key' = 'new_key', 'GEAR'= 'gear')) %>%
  dplyr::select(-c(effort, species ))

#put all column names to lowercase for consistency  
names(all_effort)<-tolower(names(all_effort))

all_effort$fish_count<-ifelse(is.na(all_effort$fish_count), 0, all_effort$fish_count)

#convert date to a year and julian (year) day 
all_effort$date<- as.Date(all_effort$sample_start_date, "%m/%d/%y") # convert "date" from chr to a Date class and specify current date format
all_effort$year<-year(all_effort$date)
all_effort$julian <- yday(all_effort$date)  

#make a dataframe with the date info to join back 
sample_dates<-dplyr::select(all_effort, new_key, survey_number, date, year, julian, gear)%>%
  mutate(
    gear2 = case_when(       #combine gear types 
      gear == "LMFYKE" ~ 'FT_NET', 
      gear == "SMFYKE" ~ 'FT_NET', 
      gear == "TRAPNET" ~ 'FT_NET', 
      gear == "GLGNET" ~ 'GILL', 
      gear == "IGNET" ~ 'GILL', 
      gear == "SEINE" ~ 'SEINE', 
      gear == "BOOMSHK" ~ 'SHOCK')) %>%
  distinct(new_key, survey_number, gear2, .keep_all = TRUE) %>%
  dplyr::select(-c(gear))

#remove duplicate surveys, combine gears
cis_new<-all_effort %>%
  group_by(new_key) %>%
  slice_max(year) %>%     #choose the most recent survey 
  mutate(
    gear2 = case_when(       #combine gear types 
      gear == "LMFYKE" ~ 'FT_NET', 
      gear == "SMFYKE" ~ 'FT_NET', 
      gear == "TRAPNET" ~ 'FT_NET', 
      gear == "GLGNET" ~ 'GILL', 
      gear == "IGNET" ~ 'GILL', 
      gear == "SEINE" ~ 'SEINE', 
      gear == "BOOMSHK" ~ 'SHOCK'))  %>%
  group_by(new_key, survey_number, gear2) %>% # group by lake, survey, and gear then sum counts and effort 
  summarize(effort_new = sum(effort),
            fish_count_new = sum(fish_count)) %>% #join 
  left_join(sample_dates)

#join for modeling: 
cis_dat_for_model<-left_join(cis_new, secchi, by=c('new_key', "survey_number")) %>% 
  left_join(drivers, by=c('new_key')) %>%
  left_join(temp_do_no_dups) %>%
  left_join(dplyr::select(lake_ll, new_key, LONG_DD, LAT_DD, FMU_Code))


#### link catch data to drivers ####
driver_vars<-read.csv("Data/driver_varibles.csv") %>%
  distinct( nhdid, .keep_all = TRUE)
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
#lmb - degree days, lake order, and Secchi depth - 3/4 of our obs have Secchi - could pull Secchi from LAGOS 


#### if you want cpue data #### 
#calculate CPUE ->abundance divided by effort 
count_data$cpue<-round((count_data$fish_count/count_data$effort), 3) 

#need to "spread" the table so that lakes with missing sp can be included as 0s in the dist model 
#rearrange table to have species along the top, will aggregate by the lake name column (#472 lakes) and gear 
sp_wide_cpue<-tidyr::pivot_wider(data= count_data, 
                                 id_cols = c(new_key.y, gear),
                                 names_from = species, 
                                 values_from = cpue, 
                                 values_fill = list(cpue = 0)) #replace NAs with 0
#  values_fn = list(cpue = sum)) #adds up all of the catch data

sp_wide_count<-tidyr::pivot_wider(data= count_data, 
                                  id_cols = c(new_key.y, gear), #note that you lose information on effort 
                                  names_from = species, 
                                  values_from = fish_count, 
                                  values_fill = list(fish_count = 0)) #replace NAs with 0


#select out cpue LMB # this keeps cpue for absences
mi_lmb_cpue<- dplyr::select(sp_wide_cpue, new_key, gear2, LMB) 
lmb_format<-tidyr::pivot_wider(data= mi_lmb_cpue, 
                               id_cols = c(new_key),
                               names_from = gear2, 
                               values_from = LMB, 
                               values_fill = list(LMB = 0))  #replace NAs with 0

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


### match to temp data from Winslow 
winslow_lakes<-read.table(file = "/Users/katelynking/Desktop/NLDAS_thermal_metrics.tsv", sep = '\t', header = TRUE) %>%
  distinct(site_id, .keep_all = TRUE)

all_lakes<-catch_data%>% distinct(ihdlkid)
lakes_with_temps<-left_join(all_lakes, winslow_lakes, by=c('ihdlkid' = 'site_id'))
  
