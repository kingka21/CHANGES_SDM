# modeling 3 level hierarchical model with varying intercepts and all gears combined
#written by Katelyn King and Henrique Giacomini
#code adapted from Tyler Wagner 
#separating catches from different gears
#add another level in the model (the site, with multiple catch values, as a random factor). 

#load libraries 
library(R2jags)
library(lme4)
library(MCMCpack)
library(dotwhisker) #https://cran.r-project.org/web/packages/dotwhisker/vignettes/dotwhisker-vignette.html
library(dplyr)

#### data ####
#need to keep all the gears used in a lake even if that gear did not catch sp of interest 

#lmb_dat_for_model has count by gear and effort by gear
#join tables for LMB 
dat<- lmb_dat_for_model %>% # 472 unique lakes, but 41 lakes don't match the drivers with the nhdid
  left_join(dplyr::select(lake_ll, new_key, LONG_DD, LAT_DD, FMU_Code))
dat<-read.csv('Data/lmb_model_data_feb10.csv')
## standardize and transform predictor values ##
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
#dat$z_houses<-as.numeric(scale(log(dat$houses_km+ 0.001))) #ok remove houses for now 

dat$logeffort <- log(dat$effort_new)

dat<-dplyr::select(dat, new_key, fish_count_new, logeffort, gear2, FMU_Code, z_lake_area, z_dd_mean, z_secchi, z_perim_km, z_order, 
                   z_ws_area_km2, z_ws_urban, z_ws_forest, z_ws_agriculture, z_ws_shrub, z_ws_wetland, z_ws_slope, z_ws_elevation) %>%
  na.omit()

#put Inf values to 0 
#is.na(dat)<-sapply(dat, is.infinite)
#dat[is.na(dat)]<-0

#check out response: log of count is mostly normally dist, so Poisson seems appropriate 
hist(dat$fish_count_new)
hist(log(dat$fish_count_new))

#################################################################
########## BUGS CODE ############################################
#################################################################
#hierarchical model with varying intercepts 
sink("model.txt")
cat("
    model {
    # Likelihood: 
    
    for (i in 1:n){         
    y[i] ~ dpois(lambda[i])  # Distribution Poisson  
    log(lambda[i]) <- log.lambda[i] # log-link 
    log.lambda[i] <- alpha[site[i]] + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i] + b[4] * x4[i] + 
                       b[5] * x5[i] + b[6] * x6[i] + b[7] * x7[i] + b[8] * x8[i] + b[9] * x9[i] + b[10] * x10[i] +
                      b[11] * x11[i] + b[12] * x12[i]  + b[13] * x13[i] + 
                      logq[gear[i]]*IND[i] + logeffort[i]
         
    } 
    
    # Level-2 of the model: site 
    for(s in 1:nsites){
    alpha[s] ~ dnorm(mu.alpha[group[s]],tau.alpha)
    }
    
     # Level-3 of the model: group (fisheries management unit regions) 
    for(j in 1:J){
    mu.alpha[j] ~ dnorm(mu.alpha2,tau.alpha2)
    }
    
    # Priors
    mu.alpha2 ~ dnorm(0, 0.0001)
    sigma.alpha ~ dunif(0,100)
    sigma.alpha2 ~ dunif(0,100)
    
    #derived quantities 
    tau.alpha <- pow(sigma.alpha,-2)
    tau.alpha2 <- pow(sigma.alpha2,-2)
    
    # priors for predictors 
    for(k in 1:13){
    b[k] ~ dnorm(0,3)
    }
    
     # priors for gear-specific log-catchabilities 
    for (k in 1:ngears) { 
    logq[k] ~ dnorm(0,  0.0000001)   #non-informative
    } 

    } # end model
    ",fill = TRUE)
sink()


### set parameters ####
# Initial values
inits <- function (){
  list (mu.alpha = rnorm(8))
}


# Parameters monitored
parameters <- c("alpha", "mu.alpha","tau.alpha", "mu.alpha2","tau.alpha2","b", 'logq')


# MCMC settings
ni <- 10000 # number of iterations 
nt <- 1 #number to thin
nb <- 5000 #number to burn  
nc <- 3 #number of chains


#### set up data #### 

#set the number of sites and create site index 
nsites<-length(unique(dat$new_key))
site <- as.numeric(as.factor(dat$new_key))

# Set the number of FMUs 
J <- length(unique(dat$FMU_Code))

# Create group index, must go from 1,..J
dat$FMU_Code <- droplevels(as.factor(dat$FMU_Code))
dat$group <- as.numeric(as.factor(dat$FMU_Code))
group <- doBy::summaryBy(group ~ new_key, data=dat, FUN=mean) #need region to be indexed by site

#set number of different sampling gears used and set a gear index 
ngears<-length(unique(dat$gear2))
gear <- as.numeric(as.factor(dat$gear2))

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
dat$IND<-ifelse(dat$gear2 == "SHOCK", 0, 1) #FT_NET is the reference gear

# Load data
data <- list(y = dat$fish_count_new, group = group$group.mean, gear=gear, n = dim(dat)[1], J = J, ngears = ngears,
             x1=dat$z_secchi, x2=dat$z_lake_area, x3=dat$z_dd_mean, x4=dat$z_perim_km, x5=dat$z_order, x6=dat$z_ws_area_km2, 
             x7=dat$z_ws_urban, x8=dat$z_ws_forest, x9=dat$z_ws_agriculture, x10=dat$z_ws_shrub, x11=dat$z_ws_wetland, 
             x12=dat$z_ws_slope, x13=dat$z_ws_elevation,
             logeffort=dat$logeffort, IND=dat$IND, nsites = nsites, site =site
) 

### run the model using JAGS and R ####
# Set timer 
start.time = Sys.time()         

# Call JAGS from R and run model 
output<- jags(data, inits, parameters, "model.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb)


# # Calculate computation time (takes about 12 mins!) 
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

# Summarize posteriors
print(output, dig = 3)
#save output 
jag.sum<-output$BUGSoutput$summary
write.table(x=jag.sum,file="out.txt",sep="\t")

#in this case, logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
#check convergence with Brooks-Gelman-Rubin statistic
which(output$BUGSoutput$summary[, c("Rhat")] > 1.1)

#### plots #### 
library(dotwhisker)

#betas
bEst <- matrix(NA, nrow=13,ncol=3)
for(i in 1:13){ #parameters
  bEst[i,1] <- mean(output$BUGSoutput$sims.list$b[,i])
  bEst[i,2:3] <- quantile(output$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
bEst<-as.data.frame(bEst)
bEst$variable<-c("Secchi", "lake area", "degree days", "perim", "order", "ws area", "ws urban", "ws forest", "ws ag",  "ws shrub",  "ws wetland",
                 "ws slope", "elevation" )

#logqs
qEst <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gear catchability without the reference gear 
  qEst[i,1] <- mean(output$BUGSoutput$sims.list$logq[,i])
  qEst[i,2:3] <- quantile(output$BUGSoutput$sims.list$logq[,i],c(0.05,0.95))
}
qEst<-as.data.frame(qEst)
qEst$variable<-c("fyke_ref", "gill", "seine", "shock")
qEst<- filter(qEst, variable != "fyke_ref") #remove reference

#combine datasets 
EstsLake <- gtools::smartbind(bEst,qEst)
colnames(EstsLake)<- c("estimate", "conf.low", "conf.high", "term")

#add colors for sig different than 0
EstsLake$color <- as.numeric(EstsLake[,2] * EstsLake[,3] > 0 )

# plot using the dotwhisker package 
lake_plot<-dwplot(EstsLake, style = "dotwhisker", 
                  dot_args = list(aes(colour = factor(color))),
                  whisker_args = list(aes(colour = factor(color))), 
                  vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("lakes") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
lake_plot
