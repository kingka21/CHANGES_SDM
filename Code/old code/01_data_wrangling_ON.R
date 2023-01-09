# code written by Katelyn King 
# created: Sept 2021 
# Data wrangling of Ontario Data 

#### load libraries ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

##### data cleaning and merging ####

#* catch data #### 
##filter out "non-standard" gearuse, no effortstatus = 3 (major problems with sampling), and no CatchExcludeFlag (trouble identifying fish)
#Geartype code 1 =large mesh and 2= small mesh 
indv_catch<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/ON_data/BsM_data_2021/Netting_Individual_Fish_2021.csv") %>% #this had fish listed individually
  select(BsM_ProjectName, GearTypeCode, SpecName, GearUse, EffortStatusCode,  CatchExcludeFlag) %>%
  filter(GearUse == "Standard" & EffortStatusCode != 3 & CatchExcludeFlag != "Yes" ) %>%
  group_by(BsM_ProjectName, GearTypeCode, SpecName) %>%
  count()

group_catch<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/ON_data/BsM_data_2021/Netting_fish_sizes_2021.csv") %>% #this has some in size classes, mostly small mesh  
  select(BsM_ProjectName, GearTypeCode, SpecName, GearUse, EffortStatusCode,  CatchExcludeFlag) %>%
  filter(GearUse == "Standard" & EffortStatusCode != 3 & CatchExcludeFlag != "Yes" ) %>% 
  group_by(BsM_ProjectName, GearTypeCode, SpecName) %>%
  count()

nets<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/ON_data/BsM_data_2021/Netting_2021.csv") %>% #effort 
  select(BsM_Cycle, BsM_ProjectName, GearTypeCode, GearGangTotal, GearUse, EffortStatusCode) %>%
  filter(GearUse == "Standard" & EffortStatusCode != 3 ) %>% 
  select(-c(GearUse, EffortStatusCode)) %>%
  group_by(BsM_ProjectName, GearTypeCode) %>% #sum up the total gangs per project and gear
  summarise(TotalGangLake = sum(GearGangTotal, na.rm = TRUE))

indv_effort<-left_join(indv_catch, nets, by = c('BsM_ProjectName', 'GearTypeCode'))
group_effort<-left_join(group_catch, nets, by = c('BsM_ProjectName', 'GearTypeCode'))

#*lake variables ####
lake_vars<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/ON_data/BsM_data_2021/lake_info_BsM_2021.csv") %>%
        select(BsM_Cycle, FMZ, BsM_ProjectName, WbyLID, WbyLat, WbyLong, SurfArea, DepthMax, GrowingDegreeDayMean8110, JulianDayFirstNet) %>%
  rename(lakeid=WbyLID, Lat=WbyLat, Lon=WbyLong, lake_area_ha=SurfArea) 
#look at duplicate rows # some lakes are in multiple cycles and some have different days of netting in the same cycle  
lakes_multi_survey<-lake_vars %>% 
 group_by(lakeid) %>% 
  filter(n() > 1)
#probably take the most recent cycle and then the most max julian date 
lake_vars<- lake_vars[!duplicated(paste(lake_vars$lakeid)),] #select only one row of each lake n=1238

#*water chemistry ####
chem<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/ON_data/BsM_data_2021/Water_Chemistry_2021.csv") %>%
  select(Cycle, lakeid, Secchi_m, Total_Phosphorus_ugL, pH, Conductivity_uScms)

#check for duplicates - will have to match by lake and cycle to get one obs per lake 
lakes_multi_survey<-chem %>% 
  group_by(lakeid, Cycle) %>% 
  filter(n() > 1)

#*DO and temp profile data ####
#will have to match by lake and cycle #0.5 m is the surface temp 
temp_do<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/ON_data/BsM_data_2021/WaterTemperature_DO_profiles_2021.csv") %>%
  dplyr::select(Cycle, lakeid, Date, Depth_m, Temp_degC, DO_mgL, Primary_Profile) %>%
  filter(Primary_Profile == "Yes")

names(temp_do)<-tolower(names(temp_do))  

#look for cases when depth is repeated (don't want this for profile calcs)
unique<-temp_do %>% 
  group_by(cycle, lakeid, date, depth_m) %>%
  mutate(n = n()) %>%
  filter(n > 1)

#calc thermo depths using the Rlakeanalyzer package : need unique depths  
thermo_depth<-temp_do %>% 
  group_by(cycle, lakeid, date) %>%
  summarise(thermo_depth_m = thermo.depth(wtr=temp_degc, depths=depth_m, seasonal = FALSE)) #location of the thermocline from a temperature profile

#plot one of the lake surveys to see if the calculated depth matches visual 
test<-filter(temp_do, lakeid == '15-3532-55170' & cycle == '1')
ggplot(test, aes(x = temp_degc, y = depth_m, colour = temp_degc)) +
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
  group_by(cycle, lakeid, date) %>%
  slice_head()  %>% 
  dplyr::select(cycle, lakeid, date, depth_m, temp_degc, do_mgl)%>% 
  rename(surface_depth_m  = depth_m, surface_temp_c= temp_degc, surface_do_mgl = do_mgl)


#pull out the bottom value for each combo 
bottom<-temp_do %>% 
  group_by(cycle, lakeid, date) %>%
  slice_tail() %>% 
  dplyr::select(cycle, lakeid, date, depth_m, temp_degc, do_mgl) %>% 
  rename(bottom_depth_m  = depth_m, bottom_temp_c= temp_degc, bottom_do_mgl = do_mgl)

#join new temp/do measures 
temp_do_ON<-left_join(top, bottom) %>%
  left_join(thermo_depth)

write.csv(temp_do_ON, "Data/ON_data/temp_do_measures_ON.csv", row.names = FALSE)

#* shorline data #### 
shoreline_dat<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/ON_data/shoreline_data_2017.csv") %>% #there is more in this dataset but it is old - so using variables that may not change 
                select(Waterbody_ID_AHI, Waterbody_LID, AHI_Year_End, Perimeter_km, Elevation_m, Crown_Land_pc, Resort_cnt, Cottage_cnt)
#lakes_multi_survey<-shoreline_dat %>%  # looks like the AHI considered some lakes as separate but the new survey has them as one ?
 # group_by(Waterbody_LID) %>% 
  #filter(n() > 1)
shoreline_dat<- shoreline_dat[!duplicated(paste(shoreline_dat$Waterbody_LID)),] #select only one row of each lake n=720


#can join all tables by Waterbody_LID 
ON_dat<-left_join(catch_dat, lake_vars, by = "Waterbody_LID") %>%
        left_join(shoreline_dat)

names(ON_dat)<-tolower(names(ON_dat))  
# split out month and year  
ON_dat<-ON_dat %>%
  separate(survey_year_month, c("year", "month"), sep = "-")
#LmBas is LMB
#Walle is WAE 
#LaHer is #CIS 

# effort is already calculated Catch-per-unit-effort (fish/net) CUE is calculated as area-weighted mean (based on CUE and Benthic area in each depth stratum)

#### LARGE MOUTH BASS DATA #### 
#select out LMB 
#use catch data not cpue 
### need to split out the large vs small mesh 
on_lmb_s<-filter(ON_dat, spc_label=='LmBas' & gear == "S")
on_lmb_l<-filter(ON_dat, spc_label=='LmBas' & gear == "L")

#### WALLEYE DATA #### 
on_wae<-filter(ON_dat, spc_label=='Walle')

#### CISCO DATA #### 
on_cis<-filter(ON_dat, spc_label=='LaHer' & gear == "S")
on_cis<-filter(ON_dat, spc_label=='LaHer' & gear == "L")
