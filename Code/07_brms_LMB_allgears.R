### Varying intercepts 3-level model 
#Bayesian approach - BRMS package 
install.packages('brms')
install.packages('rstan')
library(brms)
library(rstan)

#test run to see if brms works 
BRM1 <- brm(weight ~ Diet, data = ChickWeight) 
summary(BRM1)

#### data ####
#need to keep all the gears used in a lake even if that gear did not catch sp of interest 

#count_lmb has count, effort by gear
#join tables for LMB 
dat<-left_join(count_lmb, drivers ) %>% # 472 unique lakes, but 41 lakes don't match the drivers with the nhdid
  left_join(dplyr::select(lake_ll, new_key, LONG_DD, LAT_DD, FMU_Code))

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

dat$logeffort <- log(dat$effort)


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
model_1<-brm(fish_count ~ z_lake_area + z_dd_mean + gear + offset(logeffort) + #offset is not a parameter to be estimated - its value is added directly onto the right side of the formula  
              (1 | region:lake) + #interaction between region and lake (which is the lake level effect).#varying intercept if there is a 1 infront of the line 
              (1 | region), #main effect of region
            data = dat,
            family = poisson(link = "log"), #family agrgument to specify distribution of the response 
            prior = prior(normal(0, 1), class = b,
                          cauchy(0, 5), class = sd),
            chains = 2,
            cores = 2, 
            iter = 2000)


#using an ifelse statement 
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

if (gear<4) { # if gear is not the reference gear (1,2,3)
  log.lambda[i] <- alpha[site[i]] + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i] + logq[gear[i]] + logeffort[i]  #q and effort are specific to sample i
  #logq[gear[i]] is the log(q) for the gear
}
else (gear=4){ # if gear is the reference gear
  log.lambda[i] <- alpha[site[i]] + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i] + logeffort[i]  #logq is set to zero
}
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
bayesplot::mcmc_plot(model_1, pars = c("^r_", "^b_", "^sd_")) +
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
