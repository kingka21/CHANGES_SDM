#### size investigation ####
library(dplyr)
library(tidyr)
library(ggplot2)

lake_info<-read.csv("Data/MI_data/snt_lake_info.csv") %>% 
  select(new_key, survey_number)

effort<-read.csv( "Data/MI_data/snt_effort_data.csv") %>%
  dplyr::select(NEW_KEY, Survey_Number, GEAR,  EFFORT, SURVEY_YEAR) %>%
  rename(survey_number=Survey_Number, new_key=NEW_KEY, gear=GEAR, year = SURVEY_YEAR) %>% 
  filter(gear != "SMFYKE") %>%#remove mini fyke 
  group_by(new_key) %>%
  slice_max(year) %>%  #choose the most recent survey 
  ungroup()%>% 
  mutate(
    gear2 = case_when(       #combine gear types 
      gear == "LMFYKE" ~ 'FT_NET', 
      gear == "TRAPNET" ~ 'FT_NET', 
      gear == "GLGNET" ~ 'GILL', 
      gear == "IGNET" ~ 'GILL', 
      gear == "SEINE" ~ 'SEINE', 
      gear == "BOOMSHK" ~ 'SHOCK'))

in_grp<-read.csv("/Users/katelynking/Desktop/UofM/CHANGES/SNT_data/mi_snt_catch_inchgrp_0722.csv") %>% 
  filter(gear != "SMFYKE") %>% #remove mini fyke 
  filter(SPECIES_CODE == 'LMB') %>% 
  mutate(
    gear2 = case_when(       #combine gear types 
      gear == "LMFYKE" ~ 'FT_NET', 
      gear == "TRAPNET" ~ 'FT_NET', 
      gear == "GLGNET" ~ 'GILL', 
      gear == "IGNET" ~ 'GILL', 
      gear == "SEINE" ~ 'SEINE', 
      gear == "BOOMSHK" ~ 'SHOCK')) %>% 
  group_by(Survey_Number, Inch_Group, gear2) %>% # group by survey, inch group, and gear then sum counts and effort 
  summarize(effort_new = sum(EFFORT),
            fish_count_new = sum(Total_Number_Caught))  %>% 
  mutate(cpue= fish_count_new/effort_new,
         age= ifelse(Inch_Group <=5, 'juvenile', 'adult')) %>% ##anything < 1 year old is juvenile, from Kevin's report 5in was maax of 1 year olds  
  left_join(lake_info, by = c('Survey_Number' = 'survey_number')) %>% 
  ungroup()


#frequency plot by gear #but how to include total caught? can I extend table so each row is one fish caught? 
ggplot(in_grp, aes(x = Inch_Group)) +  geom_density(aes(color = gear2))

ggplot(in_grp, aes(x = Inch_Group, y=fish_count_new, group=gear2)) +  
  geom_line(aes(color = gear2))

ggplot(in_grp, aes(x = Inch_Group, y=cpue, group=gear2)) +  
  geom_line(aes(color = gear2))

ggplot(in_grp, aes(x = Inch_Group, y  = fish_count_new)) +
  geom_point(aes(color = gear2))

## pull out just the min and max age of a survey to compare to historical 
contemp_min<- in_grp %>% 
  group_by(Survey_Number, gear2) %>%
  slice(which.min(Inch_Group)) %>% 
  filter(gear2 == "GILL" | gear2 == "SEINE") #just gill and seine to compare 

contemp_max<- in_grp %>% 
  group_by(Survey_Number, gear2) %>%
  slice(which.max(Inch_Group)) %>% 
  filter(gear2 == "GILL" | gear2 == "SEINE") #just gill and seine to compare
 
ggplot(contemp_min, aes(x = Inch_Group)) +  geom_density(aes(color = gear2))
ggplot(contemp_max, aes(x = Inch_Group)) +  geom_density(aes(color = gear2))

#calculate length ranges by gear type of SnT data ####
contemp_min %>%                               # Summary by group using dplyr
  group_by(gear2) %>% 
  summarize(min = min(Inch_Group),
            median = median(Inch_Group),
            mean = mean(Inch_Group),
            max = max(Inch_Group), 
            mean_range = (min+max)/2)

contemp_max %>%                               # Summary by group using dplyr
  group_by(gear2) %>% 
  summarize(min = min(Inch_Group),
            median = median(Inch_Group),
            mean = mean(Inch_Group),
            max = max(Inch_Group), 
            mean_range = (min+max)/2)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### Does Juvenile abundance predict adult abundance ? - are they correlated? ####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

#include all effort so that you have 0 lmb caught even if a gear was used - need 0 juv and adult for other graphs
all_effort<-left_join(effort, in_grp, by=c('survey_number' = 'Survey_Number', 'gear2'))

#try using just count data #pivot so that you get 0s for juv and adult (this combines gears)
#there are 475 surveys comparing here 
life_stage_comp <- all_effort %>% 
  pivot_wider(id_cols = c(new_key.x, survey_number),
              names_from = age, 
              values_from = fish_count_new, 
              values_fn = sum)  

life_stage_comp[is.na(life_stage_comp)] <- 0
life_stage_comp<-life_stage_comp%>% 
  mutate(total= adult+juvenile)

cor(life_stage_comp$adult, life_stage_comp$juvenile, use="complete.obs") #0.30
cor(life_stage_comp$adult, life_stage_comp$total, use="complete.obs") #0.97
cor(life_stage_comp$juvenile, life_stage_comp$total, use="complete.obs") #0.50
plot( life_stage_comp$juvenile, life_stage_comp$adult)
plot( life_stage_comp$juvenile, life_stage_comp$total)
plot( life_stage_comp$adult, life_stage_comp$total)

#filter out lake with no lmb # 311 lakes 
life_sub<-filter(life_stage_comp, total != 0)

cor(life_sub$adult, life_sub$juvenile, use="complete.obs") #0.13
cor(life_sub$adult, life_sub$total, use="complete.obs") #0.95
cor(life_sub$juvenile, life_sub$total, use="complete.obs") #0.41
plot( life_sub$juvenile, life_sub$adult)
plot( life_sub$juvenile, life_sub$total)
plot( life_sub$adult, life_sub$total)

#compare only lakes where both seine and gill were used (split by gear if you want to use cpue )
gill_seine_combo<-all_effort %>% 
  group_by(survey_number) %>% 
  filter('GILL' %in% gear2 & 'SEINE' %in% gear2) %>% 
  filter(gear2 == "GILL" | gear2 == "SEINE")

n_distinct(gill_seine_combo$survey_number) #337 surveys use both 

gill_seine_combo <- gill_seine_combo %>% 
  pivot_wider(id_cols = survey_number,
              names_from = age, 
              values_from = fish_count_new, 
              values_fn = sum) 
  
gill_seine_combo[is.na(gill_seine_combo)] <- 0
gill_seine_combo<-gill_seine_combo%>% 
  mutate(total= adult+juvenile)

cor(gill_seine_combo$adult, gill_seine_combo$juvenile, use="complete.obs") #0.19
cor(gill_seine_combo$adult, gill_seine_combo$total, use="complete.obs") #0.68
cor(gill_seine_combo$juvenile, gill_seine_combo$total, use="complete.obs") #0.85
plot( gill_seine_combo$juvenile, gill_seine_combo$adult)
plot( gill_seine_combo$juvenile, gill_seine_combo$total)
plot( gill_seine_combo$adult, gill_seine_combo$total)

#lakes/surveys that used seine 
seine_combo<-all_effort %>% 
  group_by(survey_number) %>% 
  filter( 'SEINE' %in% gear2)

n_distinct(seine_combo$survey_number) # 361 use seine

#*how that relates to the total population using estimates of total abundance from the model 2 ####
cont_abund_pred<-read.csv("Data/output/model2_cont_abund_pred.csv") %>% 
  left_join(life_stage_comp, by = c("new_key" = "new_key.x"))

cor(cont_abund_pred$adult, cont_abund_pred$cont.med.abund, use="complete.obs") #0.32
plot( cont_abund_pred$adult, cont_abund_pred$cont.med.abund)
cor(cont_abund_pred$juvenile, cont_abund_pred$cont.med.abund, use="complete.obs") #0.21
plot( cont_abund_pred$juvenile, cont_abund_pred$cont.med.abund)

#think you can't really compare this because the count is way higher than predicted abund ( abund scaled based on gear)
#so cpue should be a better estimate 
full_lake <- all_effort %>% 
  pivot_wider(id_cols = c(new_key.x, survey_number),
              names_from = age, 
              values_from = cpue, 
              values_fn = sum)  

cont_abund_pred<-read.csv("Data/output/model2_cont_abund_pred.csv") %>% 
  left_join(full_lake, by = c("new_key" = "new_key.x"))

cor(cont_abund_pred$adult, cont_abund_pred$cont.med.abund, use="complete.obs") #0.64
plot( cont_abund_pred$adult, cont_abund_pred$cont.med.abund)
cor(cont_abund_pred$juvenile, cont_abund_pred$cont.med.abund, use="complete.obs") #0.33
plot( cont_abund_pred$juvenile, cont_abund_pred$cont.med.abund)


#old code 
#how that relates to the total population using cpue of the cont data (average of two or closer to one or the other)
cpue_age<-in_grp %>% 
  group_by(new_key, gear2, age) %>% # group by lake, and gear then sum counts and effort 
  summarize(effort_new = sum(EFFORT),
            fish_count_new = sum(Total_Number_Caught), 
            cpue_new= fish_count_new/effort_new) 
full_lake<-pivot_wider(cpue_age, 
                      id_cols = new_key,
                      names_from = age, 
                      values_from = cpue_new, 
                      values_fn = sum) %>% 
  mutate(total= sum(adult, juvenile, na.rm=TRUE))

cor(full_lake$adult, full_lake$juvenile, use="complete.obs") #0.14
cor(full_lake$adult, full_lake$total, use="complete.obs") #0.48 
cor(full_lake$juvenile, full_lake$total, use="complete.obs") #0.94 
plot( full_lake$juvenile, full_lake$adult)
plot( full_lake$juvenile, full_lake$total)
plot( full_lake$adult, full_lake$total)

### compare by gear #### 
by_gear<-pivot_wider(cpue_age, 
            id_cols = c(new_key, gear2),
            names_from = age, 
            values_from = cpue_new)
FT<-filter(by_gear, gear2 == 'FT_NET')
GILL<-filter(by_gear, gear2 == 'GILL')
SHOCK<-filter(by_gear, gear2 == 'SHOCK')
SEINE<-filter(by_gear, gear2 == 'SEINE')
ggplot(by_gear, aes( juvenile, adult)) + geom_point() + facet_grid(vars(gear2))

cor( FT$juvenile,FT$adult, use="complete.obs") # -0.05
cor( GILL$juvenile,GILL$adult, use="complete.obs") #-0.11
cor( SHOCK$juvenile,SHOCK$adult, use="complete.obs") #0.22
cor( SEINE$juvenile,SEINE$adult, use="complete.obs") #-0.24

### add length to catch data for modeling #### 
lake_info<-read.csv("Data/MI_data/snt_lake_info.csv")

#add size data 
size_dat<-select(in_grp, Inch_Group, Total_Number_Caught, Survey_Number, SPECIES_CODE, EFFORT, gear, gear2, age, new_key) %>%
  rename(species = SPECIES_CODE)
names(size_dat)<-tolower(names(size_dat))

#join raw data tables use size data in place of other catch data 
catch_data<-left_join(lake_info, size_dat, by=c("new_key", "survey_number")) %>% #join by survey number
  drop_na('gear')  #drop lakes where we don't have sample data 


#put all column names to lowercase for consistency  
names(catch_data)<-tolower(names(catch_data))

#choose the most recent survey 
MI_data_no_dups<-catch_data %>%
  group_by(new_key) %>%
  slice_max(year)

# group by lake, species, and gear then sum counts (add effort later)
count_lmb<-dplyr::select(MI_data_no_dups, new_key, survey_number, fmu_code, species, gear, age, total_number_caught) %>%
  filter( age=='adult')%>% #select adults
  group_by(new_key, survey_number,species, gear) %>%
  summarise_if(is.numeric, sum) %>% 
  ungroup() %>%
  filter(gear != "SMFYKE") #remove small mesh (mini) catches 


#link to effort table to include gear and lakes that did not catch LMB or where LMB was not present 
effort<-read.csv( "Data/MI_data/snt_effort_data.csv") %>%
  dplyr::select(NEW_KEY, Survey_Number, GEAR, EFFORT_MEASURE, EFFORT, SAMPLE_START_DATE, SAMPLE_END_DATE) %>%
  rename(survey_number=Survey_Number, new_key=NEW_KEY, gear = GEAR)

all_effort<-left_join(effort, count_lmb, by=c('new_key', "survey_number", 'gear')) %>%
  filter(gear != "SMFYKE") #remove small mesh (mini) catches 

#put all column names to lowercase for consistency  
names(all_effort)<-tolower(names(all_effort))

all_effort$total_number_caught<-ifelse(is.na(all_effort$total_number_caught), 0, all_effort$total_number_caught)

#convert date to a year and julian (year) day 
all_effort$date<- as.Date(all_effort$sample_start_date, "%m/%d/%y") # convert "date" from chr to a Date class and specify current date format
all_effort$year<-year(all_effort$date)
all_effort$julian <- yday(all_effort$date)  

#make a dataframe with the date info to join back 
sample_dates<-dplyr::select(all_effort, new_key, survey_number, date, year, julian, gear)%>%
  mutate(
    gear2 = case_when(       #combine gear types 
      gear == "LMFYKE" ~ 'FT_NET', 
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
      gear == "TRAPNET" ~ 'FT_NET', 
      gear == "GLGNET" ~ 'GILL', 
      gear == "IGNET" ~ 'GILL', 
      gear == "SEINE" ~ 'SEINE', 
      gear == "BOOMSHK" ~ 'SHOCK'))  %>%
  group_by(new_key, survey_number, gear2) %>% # group by lake, survey, and gear then sum counts and effort 
  summarize(effort_new = sum(effort),
            fish_count_new = sum(total_number_caught)) %>% #join 
  left_join(sample_dates)

#join with all the other data from data_wrangling_MI
write.csv(lmb_dat_for_model, "Data/lmb_dat_for_model_adult.csv", row.names = FALSE) #jul removed mini fyke and added pike/wae cpue

###Compare length ranges of SNT and historical data to see that populations are similar ####
historical_lengths<-read.csv("/Users/katelynking/Desktop/subset_data_lmb/lmb_fish_length.csv",na.strings = c("", "NA")) %>% 
  drop_na(lmb_length_min) %>%
  mutate(lmb_length_min = as.numeric(lmb_length_min))%>%
mutate(
  unit2 = ifelse( is.na(units) & lmb_length_min >20, 'mm', NA),  
  unit2 = ifelse( is.na(units) & lmb_length_min <=20, 'inches', unit2),
  unit2 = ifelse( is.na(unit2), units, unit2)
) %>%
  mutate(length_min = ifelse(unit2== "inches", lmb_length_min, NA),
         length_min = ifelse(unit2== "mm", lmb_length_min/25.4, length_min), 
         length_max = ifelse(unit2== "inches", lmb_length_max, NA),
         length_max = ifelse(unit2== "mm", lmb_length_max/25.4, length_max)
  ) %>% 
  mutate(inch_group_min = case_when(
    length_min < 1 ~ "0",
    length_min < 2 ~ "1",
    length_min < 3 ~ "2",
    length_min < 4 ~ "3",
    length_min < 5 ~ "4",
    length_min < 6 ~ "5",
    length_min < 7 ~ "6",
    length_min < 8 ~ "7",
    length_min < 9 ~ "8",
    length_min < 10 ~ "9",
    length_min < 11 ~ "10",
    length_min < 12 ~ "11",
    length_min < 13 ~ "12",
    length_min < 14 ~ "13",
    length_min < 15 ~ "14",
    length_min < 16 ~ "15",
    length_min < 17 ~ "16",
    length_min < 18 ~ "17",
    length_min < 19 ~ "18",
  ), 
  inch_group_max = case_when(
      length_max < 1 ~ "0",
      length_max < 2 ~ "1",
      length_max < 3 ~ "2",
      length_max < 4 ~ "3",
      length_max < 5 ~ "4",
      length_max < 6 ~ "5",
      length_max < 7 ~ "6",
      length_max < 8 ~ "7",
      length_max < 9 ~ "8",
      length_max < 10 ~ "9",
      length_max < 11 ~ "10",
      length_max < 12 ~ "11",
      length_max < 13 ~ "12",
      length_max < 14 ~ "13",
      length_max < 15 ~ "14",
      length_max < 16 ~ "15",
      length_max < 17 ~ "16",
      length_max < 18 ~ "17",
      length_max < 19 ~ "18",
      length_max < 20 ~ "19",
    )) %>% 
  mutate(
    inch_group_max = ifelse(is.na(inch_group_max), inch_group_min, inch_group_max) 
  )

ggplot(historical_lengths, aes(x = inch_group_min)) +  geom_density(aes(color = gear2))  
#write.csv(historical_lengths, "/Users/katelynking/Desktop/subset_data_lmb/lmb_fish_length_clean.csv")

#would need to go back to historical data code and include length there by subject.id before combining info from the same gear across cards. 
#this would be done in the "data2" step in kk_aggregate_for_model.R 

hist_dat<-read.csv("Data/historical_lmb2.csv")
n_distinct(hist_dat$new_key
           )
