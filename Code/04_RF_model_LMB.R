install.packages('randomForest')
#### R libraries ####
library(dplyr)
library(randomForest)
library(ggplot2)

#### load data ####
shock_dat<-select(lmb_count_dat, SHOCK, lake_area_m2, geom_ratio, lake_order,  ws_urban_prop, ws_forest_prop, ws_agriculture_prop, ws_shrub_prop, ws_wetland_prop, ws_slope_deg, ws_mean_elevation_m, dd_mean, houses_km, secchi_m)
fyke_dat<-select(lmb_count_dat, FT_NET, lake_area_m2, geom_ratio, lake_order,  ws_urban_prop, ws_forest_prop, ws_agriculture_prop, ws_shrub_prop, ws_wetland_prop, ws_slope_deg, ws_mean_elevation_m, dd_mean, houses_km, secchi_m)
gill_dat<-select(lmb_count_dat, GILL, lake_area_m2, geom_ratio, lake_order,  ws_urban_prop, ws_forest_prop, ws_agriculture_prop, ws_shrub_prop, ws_wetland_prop, ws_slope_deg, ws_mean_elevation_m, dd_mean, houses_km, secchi_m)
seine_dat<-select(lmb_count_dat, SEINE, lake_area_m2, geom_ratio, lake_order, ws_urban_prop, ws_forest_prop, ws_agriculture_prop, ws_shrub_prop, ws_wetland_prop, ws_slope_deg, ws_mean_elevation_m, dd_mean, houses_km, secchi_m)


### try to see the top drivers of CPUE 
## Random forest procedure
#set up parameters 
ntreez <- 500
par(mfrow=c(2,2))
par(mar=c(2,2,1,1), oma=c(2,2,0,0))
set_string <-c(8, 188, 18, 999) #TRY DIFFERENT SEED SETS 188, 18, 999 top predictors don't change 

### shock #### 
noNAs<-filter(shock_dat, SHOCK > 0)
for (i in 1:length(set_string)) {
  set.seed(set_string[i])
  RF_shock <- randomForest(log(SHOCK) ~  lake_area_m2+ geom_ratio + lake_order+ 
                        ws_urban_prop+ ws_forest_prop+ ws_agriculture_prop+ ws_shrub_prop+ ws_wetland_prop+ ws_slope_deg+ ws_mean_elevation_m+ 
                        dd_mean+ houses_km+ secchi_m, 
                        data=noNAs, ntree=ntreez, importance=T, na.action=na.omit)
  RF_shock
  shock_imp <-randomForest::importance(RF_shock, type=1, scale=FALSE) #mean decrease in accuracy (also called permutation accuracy importance).
  imp<-as.data.frame(shock_imp) 
  imp$'Pred'   <-rownames(imp)
  imp <- structure(imp$`%IncMSE`, names = as.character(imp$Pred))
  
  dotchart(imp[order(imp)], xlab = "mean decrease in accuracy",
           main = set_string[i])
  
} 

# diagnostic for if I ran enough trees to reduce MSE  
plot(RF_shock, xlim=c(0,500))

### fyke #### 
noNAs<-filter(fyke_dat, FT_NET > 0)
for (i in 1:length(set_string)) {
  set.seed(set_string[i])
  RF_fyke <- randomForest(log(FT_NET) ~  lake_area_km2+ lake_perim_km+ lake_order+ ws_area_km2+ 
                             ws_urban_prop+ ws_forest_prop+ ws_agriculture_prop+ ws_shrub_prop+ ws_wetland_prop+ ws_slope_deg+ ws_mean_elevation_m+ 
                             dd_mean+ houses_km+ secchi_m, 
                           data=noNAs, ntree=ntreez, importance=T, na.action=na.omit)
  RF_fyke
  shock_imp <-randomForest::importance(RF_fyke, type=1, scale=FALSE) #mean decrease in accuracy (also called permutation accuracy importance).
  imp<-as.data.frame(shock_imp) 
  imp$'Pred'   <-rownames(imp)
  imp <- structure(imp$`%IncMSE`, names = as.character(imp$Pred))
  
  dotchart(imp[order(imp)], xlab = "mean decrease in accuracy",
           main = set_string[i])
  
} 

# diagnostic for if I ran enough trees to reduce MSE  
plot(RF_fyke, xlim=c(0,500))

### gill #### 
noNAs<-filter(gill_dat, GILL > 0 ) %>% drop_na()
for (i in 1:length(set_string)) {
  set.seed(set_string[i])
  RF_gill <- randomForest(log(GILL) ~  lake_area_m2+ geom_ratio+ lake_order+ 
                            ws_urban_prop+ ws_forest_prop+ ws_agriculture_prop+ ws_shrub_prop+ ws_wetland_prop+ ws_slope_deg+ ws_mean_elevation_m+ 
                            dd_mean+ houses_km+ secchi_m, 
                          data=noNAs, ntree=ntreez, importance=T, na.action=na.omit)
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

### seine #### 
noNAs<-filter(seine_dat, SEINE > 0)
for (i in 1:length(set_string)) {
  set.seed(set_string[i])
  RF_seine <- randomForest(log(SEINE) ~  lake_area_km2+ lake_perim_km+ lake_order+ ws_area_km2+ 
                            ws_urban_prop+ ws_forest_prop+ ws_agriculture_prop+ ws_shrub_prop+ ws_wetland_prop+ ws_slope_deg+ ws_mean_elevation_m+ 
                            dd_mean+ houses_km+ secchi_m, 
                          data=noNAs, ntree=ntreez, importance=T, na.action=na.omit)
  RF_seine
  shock_imp <-randomForest::importance(RF_seine, type=1, scale=FALSE) #mean decrease in accuracy (also called permutation accuracy importance).
  imp<-as.data.frame(shock_imp) 
  imp$'Pred'   <-rownames(imp)
  imp <- structure(imp$`%IncMSE`, names = as.character(imp$Pred))
  
  dotchart(imp[order(imp)], xlab = "mean decrease in accuracy",
           main = set_string[i])
  
} 

# diagnostic for if I ran enough trees to reduce MSE  
plot(RF_seine, xlim=c(0,500))

#### INVESTIGATING EFFECTS ####
#*Partial Dependence Plots ####
#There are a number of R packages that implement PDPs. 
# the iml package, pdp, or DALEX.
library(iml)
X <- noNAs[which(names(noNAs) != "GILL")]
predictor <- Predictor$new(RF_gill, data = X, y = noNAs$GILL)
pdp.obj = FeatureEffect$new(predictor, feature = "lake_area_m2") #look at a specific pred
pdp.obj$plot()
pdp.obj$set.feature("dd_mean") # look at another feature 
pdp.obj$plot()
pdp.obj$set.feature("geom_ratio") # look at another feature 
pdp.obj$plot()
effs <- FeatureEffects$new(predictor) #look at all of the preds 
plot(effs)

pdp.obj = FeatureEffect$new(predictor, feature = c("lake_area_m2", "houses_km")) #look at two predictors 
pdp.obj$plot()
### how strongly features interact with each other. 
#The interaction measure regards how much of the variance of $f(x)$ is explained by the interaction. 
#The measure is between 0 (no interaction) and 1 (= 100% of variance of $f(x)$ due to interactions).
interact = Interaction$new(predictor)
plot(interact)

interact = Interaction$new(predictor, feature = "lake_area_m2") #see how lake area interacts with everything else 
plot(interact)
