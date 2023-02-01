# CHANGES_SDM


This repository has the data and code that supports the manuscript currently under review at Global Change Biology: ** Katelyn B.S. King, Henrique C. Giacomini, Kevin Wehrly, Hernán López-Fernández, Andrea K. Thomer, Karen M. Alofs
 Using historical catch data to evaluate predicted changes in fish relative abundance in response to a warming climate 
**  

The manuscript investigates if models built using spatial environmental gradients can reliably predict population changes through time in response to climate warming. We use largemouth bass abundance across Michigan lakes to validate hindcast. 


**The 'Code' folder includes the following:** \
**01_data_explore** includes exploring driver variables \
**02_map_code** includes code for mapping sample sites
**03_model1** includes the methods for running model 1 and calculating deviance,  residual deviance, and Bayes R2 \
**03_model2** includes the methods for running model 2 and calculating deviance,  residual deviance, and Bayes R2 \
**04_deviance_resid** includes the methods for comparing the residuals from the two time periods \
**05_effects_plots_model1** includes the methods for creating effects plots for model1 \
**05_effects_plots_model2** includes the methods for creating effects plots for model2 \
**06_obs_vs_pred_Fig6** includes the methods for plotting the observed vs predicted values from the two time periods \
**07_hindcast_abundance** includes the methods for hindcasting abundace from model 2 \


**The 'Data' folder includes the following:** \
**contemp_lmb_dat.csv** which are the contemporary lake attributes and catch data used in the models. See the Metadata.xlsx file which describes the variables. \
**historical_lmb2_with_secchi.csv** which are the historical lake attributes and catch data used to validate the models. See the Metadata.xlsx file which describes the variables. 
**Humphries_table.csv** which is the authority file for lake ids and latitude/longitude for Michigan lakes. See the Metadata.xlsx file which describes the variables. \
**lake_surface_temp.csv** which are the modeled yearly lake temperatures from James Breck Michigan Department of Natural Resources. See the Metadata.xlsx file which describes the variables. \
**MDNR_FMUs** folder contains the shapefiles for the Fisheries Management Units \
**output** folder contains the density estimates, temperature changes, and the residual deviance for all of the folds from the model for the contemporary and historical catch data. \

**The 'figures' folder includes the figures made in R for the manuscript**
