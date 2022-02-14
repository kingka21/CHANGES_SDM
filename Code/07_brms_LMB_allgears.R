### Varying intercepts 3-level model 
#Bayesian approach - BRMS package 
install.packages('brms')
install.packages('rstan')
install.packages('ggmcmc') #package for analyzing MCMC output 
install.packages('tidybayes')
library(brms)
library(rstan)
library(ggmcmc)
library(tidybayes)
library(dplyr)

#test run to see if brms works 
BRM1 <- brm(weight ~ Diet, data = ChickWeight) 
summary(BRM1)

#### data ####
#need to keep all the gears used in a lake even if that gear did not catch sp of interest 
#need to add a gear with effor 0 and catch 0 to all lakes 

#count_lmb has count, effort by gear
#join tables for LMB 
#dat<-lmb_dat_for_model %>% # 472 unique lakes, but 41 lakes don't match the drivers with the nhdid
  #left_join(dplyr::select(lake_ll, new_key, LONG_DD, LAT_DD, FMU_Code))
#save data for now so that I can update model easily while still in the trial and error phase 
#write.csv(dat, 'Data/lmb_model_data_feb10.csv', row.names = FALSE)
dat<-read.csv('Data/lmb_model_data_feb10.csv')

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
dat$IND<-ifelse(dat$gear2 == "FT_NET", 0, 1) #FT_NET is the reference gear


## standardize and transform predictor values ##
# Note that scaling numeric predictors benefits here and makes specifying the prior easier as well.
dat$z_lake_area<-as.numeric(scale(log(dat$lake_area_m2))) #good 
dat$z_dd_mean<-as.numeric(scale(log(dat$dd_mean))) #good
dat$z_secchi<-as.numeric(scale(log(dat$secchi_m)))
dat$z_order<-as.numeric(scale(log(dat$lake_order+0.001))) # min is 0 because these are the isolated
dat$z_perim_km<-as.numeric(scale(log(dat$lake_perim_km))) 
dat$z_ws_area_km2<-as.numeric(scale(log(dat$ws_area_km2)))  
dat$z_ws_urban<-as.numeric(scale(asin(sqrt(dat$ws_urban_prop)))) #good
dat$z_ws_agriculture<-as.numeric(scale(asin(sqrt(dat$ws_agriculture_prop)))) #ok
dat$z_ws_forest<-as.numeric(scale(asin(sqrt(dat$ws_forest_prop)))) #good
dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat$ws_wetland_prop)))) #good
dat$z_ws_shrub<-as.numeric(scale(asin(sqrt(dat$ws_shrub_prop)))) #ok
dat$z_ws_slope<-as.numeric(scale(log(dat$ws_slope_deg + 0.001))) #has 0s so added 0.001 
dat$z_ws_elevation<-as.numeric(scale(log(dat$ws_mean_elevation_m))) #good 
dat$z_houses<-as.numeric(scale(log(dat$houses_km+ 0.001))) #ok 

dat$logeffort <- log(dat$effort_new)

#### plot some data #####
ggplot(data  = dat,
       aes(x = z_dd_mean,
           y = log(fish_count_new)))+
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter")+# to add some random noise for plotting purposes
  geom_smooth(method = lm,
              se     = FALSE, 
              col    = "black",
              size   = .5, 
              alpha  = .8)+ # to add regression line
  theme_minimal()+
  labs(title = "fish count vs. DD temp ")

#### plot by group #####
ggplot(data  = dat,
       aes(x = z_dd_mean,
           y = log(fish_count_new), 
           col = gear2, 
           group = gear2))+ #add colors for gear
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter")+# to add some random noise for plotting purposes
  geom_smooth(method = lm,
              se     = FALSE, 
              size   = .5, 
              alpha  = .8)+ # to add regression line
  theme_minimal()+
  labs(title = "fish count vs. DD temp ")

#### intercept only model #### 
interceptonlymodel <- brm(fish_count_new ~ 1 + (1|FMU_Code),  
                          data = dat, 
                          warmup = 1000, iter = 3000, 
                          cores = 2, chains = 2, 
                          seed = 123
                          ) 
summary(interceptonlymodel)
hyp <- "sd_FMU_Code__Intercept^2 / (sd_FMU_Code__Intercept^2 + sigma^2) = 0"
hypothesis(interceptonlymodel, hyp, class = NULL)

#### first level predictors #### 
#note that for this model you will have to ignore the convergence checks for anything with IND because it is neutral 
model_1<-brm(fish_count_new ~ 1 + z_lake_area + z_dd_mean + z_secchi + julian + maxdepth_m + gear2*IND + offset(logeffort), #+ #offset is not a parameter to be estimated - its value is added directly onto the right side of the formula  #a “1” in the formula the function indicates the intercept.
               #(1 | FMU_Code),  # the lake level effect).#varying intercept if there is a 1 in front of the line 
             data = dat,
             family = poisson(link = "log"), #family argument to specify distribution of the response 
             chains = 2,cores = 2, 
             warmup = 1000, iter = 3000)
summary(model_1)
#posterior checks 
plot(model_1) #can use the brms default to plot or the code below 

#plot conditional effects of each pred
conditional_effects(model_1)
get_variables(model_1) #see all of the estimated parameters 

model1tranformed <- ggs(model_1) # the ggs function transforms the brms output into a longformat tibble, that we can use to make different types of plots.

#caterpillar plots 
ggplot(filter(model1tranformed, Parameter %in% c("b_Intercept", "b_z_lake_area", "b_z_dd_mean")),
       aes(x   = Iteration,
           y   = value, 
           col = as.factor(Chain)))+
  geom_line() +
  geom_vline(xintercept = 1000)+
  facet_grid(Parameter ~ . ,
             scale  = 'free_y',
             switch = 'y')+
  labs(title = "Caterpillar Plots", 
       col   = "Chains")

# posterior density plots 
ggplot(filter(model1tranformed,Parameter == "b_Intercept", Iteration > 1000),aes(x = value))+
  geom_density(fill  = "yellow", alpha = .5)+
  geom_vline(xintercept = 0, col  = "red", size = 1)+ 
  geom_vline(xintercept = summary(model_1)$fixed[1,3], #lower CI 
             col = "blue",
             linetype = 2) +
  geom_vline(xintercept = summary(model_1)$fixed[1,4], #upper CI 
             col = "blue",
             linetype = 2) +
  theme_light() +
  labs(title = "Posterior Density of Intercept")

#lake area
ggplot(filter(model1tranformed,Parameter == "b_z_lake_area", Iteration > 1000), aes(x = value))+
  geom_density(fill  = "orange", alpha = .5)+
  geom_vline(xintercept = 0,  col  = "red", size = 1)+ 
  geom_vline(xintercept = summary(model_1)$fixed[2,3], col = "blue",  linetype = 2) + #lower CI 
  geom_vline(xintercept = summary(model_1)$fixed[2,4], col = "blue",linetype = 2) +  #upper CI 
  theme_light() +
  labs(title = "Posterior Density of area")

#temp
ggplot(filter(model1tranformed, Parameter == "b_z_dd_mean", Iteration > 1000), aes(x = value))+
  geom_density(fill  = "orange",   alpha = .5)+
  geom_vline(xintercept = 0, col  = "red",  size = 1)+ 
  geom_vline(xintercept = summary(model_1)$fixed[3,3], col = "blue",  linetype = 2) + #lower CI
  geom_vline(xintercept = summary(model_1)$fixed[3,4],  col = "blue", linetype = 2) + #upper CI 
  theme_light() +
  labs(title = "Posterior Density of temp")

#### add 2nd level predictors #### 
#just add the predictors e.g. a hu4 variable 
#fish_count_new ~ z_lake_area + z_dd_mean + ag_hu4 + offset(logeffort) + (1 | gear2)

#calculate explained variance at each level if you have predictors at a different level 
#use the intercept only model and residual variance at the group-level (level 2) is 0.85^2 - 2level model 0.55^2 / intercept only model 0.85^2 = 0.58 
(0.85^2 -  0.55^2) / 0.85^2
#residual variance family level (individual level 1): intercept only model 1.11^2 - 2level model 0.77^2 / intercept only model 1.11^2 = 0.52
(1.11^2 - 0.77^2) /  1.11^2

#### random slopes  #### 
model2<-brm(fish_count_new ~ 1 + z_lake_area + z_dd_mean + offset(logeffort) + #offset is not a parameter to be estimated - its value is added directly onto the right side of the formula  #a “1” in the formula the function indicates the intercept.
               (1 + z_lake_area + z_dd_mean | gear2),  # varying slopes for both predictor variables 
             data = dat,
             family = poisson(link = "log"), #family argument to specify distribution of the response 
             chains = 2,cores = 2, 
             warmup = 1000, iter = 3000)
summary(model2)
get_variables(model2) #see all of the estimated parameters 

#plot varying slopes and intercepts


#nice summary table of gears 
model2 %>%
  spread_draws(r_gear2[gear,term]) %>%
  summarise_draws() #

model2 %>%
  spread_draws(r_gear2[gear,term]) %>%
  mean_qi() #get the mean and 95% CI for each variable #can also do median_qi 

#### posterior means and predictions ####
#plot posterior means 
dat %>%
  modelr::data_grid(gear2, z_lake_area, z_dd_mean, logeffort) %>% #think you have to include all variables in the model 
  add_epred_draws(model2) %>%
  ggplot(aes(x = .epred, y = gear2)) +
  stat_halfeye(.width = c(.95)) #specify CI usiung the width command 

#plot posterior prediction distributions 
dat %>%
  modelr::data_grid(gear2, z_lake_area, z_dd_mean, logeffort) %>%
  add_predicted_draws(model2) %>%
  ggplot(aes(x = .prediction, y = gear2)) +
  ggdist::stat_slab()

#this timed out - not enough memory 
dat %>%
  group_by(gear2) %>%
  modelr::data_grid(gear2, z_lake_area, z_dd_mean, logeffort) %>%
  add_epred_draws(model2) %>%
  ggplot(aes(x = z_lake_area, y = fish_count_new, color = ordered(gear2))) +
  stat_lineribbon(aes(y = .epred)) +
  geom_point(data = dat) +
  scale_fill_brewer(palette = "Greys") +
  scale_color_brewer(palette = "Set2")

#### cross-scale interactions #### 
#texp is a class level  (level-2) variable and extrav is an individual level (level-1) variable 
model5 <- brm(popular ~ 1 + sex + extrav + texp + extrav:texp + #the : between them is the cross-scale interaction
                (1 + extrav|class), #still include the random slopes of extrav with class 
              data  = popular2data, warmup = 1000,
              iter  = 3000, chains = 2, 
              seed  = 123, control = list(adapt_delta = 0.97),
              cores = 2) # 

#plot varying slopes and intercepts  of the CSI 
ggplot(data = popular2data, 
       aes(x   = extrav,
           y   = popular,
           col = as.factor(texp)))+
  viridis::scale_color_viridis(discrete = TRUE)+
  geom_point(size     = .7,
             alpha    = .8,
             position = "jitter")+
  geom_smooth(method = lm,
              se     = FALSE, 
              size   = 2,
              alpha  = .8)+
  theme_minimal()+
  labs(title    = "Linear Relationship between Different Years of Teacher Experience and Extrav")

#### setting up the model in BRMS #### 
for (i in 1:n){         
  y[i] ~ dpois(lambda[i])  # Distribution Poisson  
  log(lambda[i]) <- log.lambda[i] # log-link 
  log.lambda[i] <- alpha[site[i]] + b[1] * x1[i] + b[2] * x2[i] + logq[gear[i]]*IND[i] + logeffort[i]
  
} 

# Level-2 of the model: site / lake 
for(s in 1:nsites){
  alpha[s] ~ dnorm(mu.alpha[group[s]],tau.alpha)
}

# Level-3 of the model: group (fisheries management unit regions) 
for(j in 1:J){
  mu.alpha[j] ~ dnorm(mu.alpha2,tau.alpha2)
}


# bays package brms 
#bayes_R2 function to get R2
#( <varying parameter(s)> | <grouping variable(s)> )
#Varying intercepts are regularized by estimating how diverse the clusters are while estimating the features of each cluster. 
#A major benefit of using varying effects estimates is that they provide more accurate estimates of the individual cluster intercepts. 
#On average, the varying effects actually provide a better estimate of the individual cluster means. The reason that the varying intercepts provide better estimates is that they do a better job of trading off underfitting and overfitting.
model_1<-brm(fish_count_new ~ 1 + z_lake_area + z_dd_mean + gear2 + offset(logeffort) + #offset is not a parameter to be estimated - its value is added directly onto the right side of the formula  
              (1 | gear2),  # the lake level effect).#varying intercept if there is a 1 in front of the line 
            data = dat,
            family = poisson(link = "log"), #family argument to specify distribution of the response 
            chains = 2,
            cores = 2, 
            iter = 4000)

#old code with 3 levels - not sure where I found this 
model_1<-brm(fish_count_new ~ z_lake_area + z_dd_mean + gear2 + offset(logeffort) + #offset is not a parameter to be estimated - its value is added directly onto the right side of the formula  
               (1 | FMU_Code:new_key) + #interaction between region and lake (which is the lake level effect).#varying intercept if there is a 1 in front of the line 
               (1 | FMU_Code), #main effect of region
             data = dat,
             family = poisson(link = "log"), #family argument to specify distribution of the response 
             prior = prior(normal(0, 1), class = b,
                           cauchy(0, 5), class = sd),
             chains = 2,
             cores = 2, 
             iter = 2000)
# example using an ifelse statement 
for (x in 1:N)
  if (x < t)
    y[i] ~ normal(a + b * m, sigma1)
else if (x >= t)
  y[i] ~ normal(p + exp(b * r), sigma2)

xlt = x < t
bform <- bf(
  xlt * (a + b * m) + (1 - xlt) * (p + exp(b * r)),
  a + b + p ~ 1,
  sigma ~ 0 + xlt
  nl = TRUE
)

xlt = gear < 4
bform <- bf(
  xlt * (z_lake_area + z_dd_mean + logq*gear + offset(logeffort)) + (1 - xlt) * (z_lake_area + z_dd_mean + offset(logeffort)),
  a + b + p ~ 1,
  nl = TRUE
)

#what I want to do 
if (gear<4) { # if gear is not the reference gear (1,2,3)
  log.lambda[i] <- alpha[site[i]] + b[1]*x1[i] + b[2]*x2[i] + b[3]*x3[i] + logq[gear[i]] + logeffort[i]  #q and effort are specific to sample i
  #logq[gear[i]] is the log(q) for the gear
}
else (gear=4){ # if gear is the reference gear
  log.lambda[i] <- alpha[site[i]] + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i] + logeffort[i]  #logq catchability is set to zero
}


#PRIORS# 
#fixed effect regression coefficients, normal and student t would be the most common prior distributions
#the brms defaults to a uniform/improper prior, which is a poor choice. 
#You will want to set this for your models. Note that scaling numeric predictors benefits here and makes specifying the prior easier as well.
#this prior is set on an intercept that results when internally centering all population-level predictors around zero to improve sampling efficiency.



##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### check the trace plots ####
##### ##### ##### ##### ##### ##### ##### ##### ##### 
library(bayesplot)

color_scheme_set("orange")

post <- posterior_samples(model_1, add_chain = T)

post %>% 
  mcmc_trace(pars = vars(-iter, -lp__),
             facet_args = list(ncol = 4), 
             size = .15) +
  theme(legend.position = "none")

#### SUMMARIES #### 
print(model_1) #hierarchical  group summaries shown
posterior_summary(model_1) %>% round(digits = 2) #all parameter summaries shown
#coefficient plots similar to mine 
mcmc_plot(model_1, pars = c("^r_", "^b_", "^sd_")) +
  theme(axis.text.y = element_text(hjust = 0))

#check out variation among obs within each of the grouping variables (region and lake), 
#could remove one grouping variable if the within group variation is small, then compare models using the WAIC  
post %>%
  pivot_longer(starts_with("sd")) %>% 
  ggplot(aes(x = value, fill = name)) +
  geom_density(size = 0, alpha = 3/4, adjust = 2/3, show.legend = F) +
  annotate(geom = "text", x = 0.67, y = 2, label = "block", color = "orange4") +
  annotate(geom = "text", x = 2.725, y = 0.5, label = "actor", color = "orange1") +
 scale_fill_manual(values = str_c("orange", c(1, 4))) +
  scale_y_continuous(NULL, breaks = NULL) +
  ggtitle(expression(sigma["<group>"])) +
  coord_cartesian(xlim = c(0, 4))

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
#### compare WAIC for two models #### 
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
model_1 <- add_criterion(model_1, "waic")
model_2 <- add_criterion(model_2, "waic")

w <- loo_compare(model_1, model_2, criterion = "waic")

print(w, simplify = F) #look at p_waic for results - the lower the better 
