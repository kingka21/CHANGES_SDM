# Random Forest model for Ontario Data 

#### R libraries ####
library(dplyr)
library(randomForest)
library(ggplot2)

#### load data ####
gill_dat<-select(on_lmb_l, cue_fish_gang, area_ha, perimeter_km, depth_max, secchi_su, totalphosphorus_ugl, ph_pctl, conductivity_uscms, wc_tds, dd5_8110, elevation_m, cottage_cnt)
gill_dat$houses_km<-gill_dat$cottage_cnt/gill_dat$perimeter_km
#don't have lake order or any of the ws variables,  but have in-situ data  

### try to see the top drivers of CPUE 
## Random forest procedure
#set up parameters 
ntreez <- 500
par(mfrow=c(2,2))
par(mar=c(2,2,1,1), oma=c(2,2,0,0))
set_string <-c(8, 188, 18, 999) #TRY DIFFERENT SEED SETS 188, 18, 999 top predictors don't change 

### gill net (large and small combined) #### 
for (i in 1:length(set_string)) {
  set.seed(set_string[i])
  RF_gill <- randomForest(log(cue_fish_gang) ~  area_ha+ perimeter_km+ depth_max+ secchi_su+ 
                           totalphosphorus_ugl+ ph_pctl+ conductivity_uscms+ wc_tds+ dd5_8110+ elevation_m+ houses_km,  
                           data=gill_dat, ntree=ntreez, importance=T, na.action=na.omit)
  RF_gill
  shock_imp <-randomForest::importance(RF_gill, type=1, scale=FALSE) #mean decrease in accuracy (also called permutation accuracy importance).
  imp<-as.data.frame(shock_imp) 
  imp$'Pred'   <-rownames(imp)
  imp <- structure(imp$`%IncMSE`, names = as.character(imp$Pred))
  
  dotchart(imp[order(imp)], xlab = "mean decrease in accuracy",
           main = set_string[i])
  
} 

# diagnostic for if I ran enough trees to reduce MSE  
plot(RF_gill, xlim=c(0,500))

