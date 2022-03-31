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
library(ggmcmc) # package for analyzing output 


#### data ####
#need to keep all the gears used in a lake even if that gear did not catch sp of interest 

#lmb_dat_for_model has count by gear and effort by gear
dat<-read.csv("Data/lmb_dat_for_model_mar17.csv") 

dd_temp<-dd_temp %>%
  group_by(IHDLKID) %>%
  slice_max(HECTARES) %>% 
  rename(nhdid = IHDLKID)
dd_temp_long<-pivot_longer(dd_temp, 
                           cols= starts_with("DD_"),
                           names_to = "year", 
                           names_prefix = "DD_",
                           values_to = "dd_year")
dd_temp_long$year<-as.integer(as.character(dd_temp_long$year))
dat<-left_join(dat, dd_temp_long, by = c("nhdid", "year"))

## standardize and transform predictor values ##
#dat$z_order<-as.numeric(scale(log(dat$lake_order+0.001))) # min is 0 because these are the isolated
#dat$z_ws_urban<-as.numeric(scale(asin(sqrt(dat$ws_urban_prop)))) #good
#dat$z_ws_agriculture<-as.numeric(scale(asin(sqrt(dat$ws_agriculture_prop)))) #ok
#dat$z_ws_elevation<-as.numeric(scale(log(dat$ws_mean_elevation_m))) #good 
#dat$z_houses<-as.numeric(scale(log(dat$houses_km+ 0.001))) #ok remove houses for now 
dat$z_lake_area<-as.numeric(scale(log(dat$lake_area_m2))) #good 
dat$z_max_depth<-as.numeric(scale(log(dat$maxdepth_m))) #good 
dat$z_dd_mean<-as.numeric(scale(log(dat$dd_year))) #good
dat$z_surface_temp_mean<-as.numeric(scale(log(dat$surface_temp_mean))) #good 
dat$z_secchi<-as.numeric(scale(log(dat$secchi_m)))
dat$z_bottom_do<-as.numeric(scale(log(dat$bottom_do_mgl+ 0.001))) ##has 0s so added 0.001 
dat$z_julian<-as.numeric(scale(dat$julian)) ##standardize julian date 
dat$z_ws_forest<-as.numeric(scale(asin(sqrt(dat$ws_forest_prop)))) #good
dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat$ws_wetland_prop)))) #good

dat$logeffort <- log(dat$effort_new)

dat<-dplyr::select(dat, new_key, fish_count_new, logeffort, gear2, FMU_Code, 
                   z_lake_area, z_max_depth, z_secchi, z_julian, 
                   z_bottom_do, z_dd_mean,
                   z_ws_forest, z_ws_wetland, pike_pres, wae_pres) %>%
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
    for(k in 1:10){
    b[k] ~ dnorm(0,0.0000001)
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
ni <- 80000 # number of iterations 
nt <- 3 #number to thin
nb <- 30000 #number to burn  
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
dat$IND<-ifelse(dat$gear2 == "FT_NET", 0, 1) #FT_NET is the reference gear

# Load data
data <- list(y = dat$fish_count_new, group = group$group.mean, gear=gear, n = dim(dat)[1], J = J, ngears = ngears,
             x1=dat$z_secchi, x2=dat$z_lake_area, x3=dat$z_dd_mean, x4=dat$z_max_depth, x5=dat$z_bottom_do,
             x6=dat$z_ws_forest,x7=dat$z_ws_wetland, x8=dat$z_julian, x9=dat$wae_pres, x10=dat$pike_pres,
             logeffort=dat$logeffort, IND=dat$IND, nsites = nsites, site =site
) 

### run the model using JAGS and R ####
# Set timer 
start.time = Sys.time()         

# Call JAGS from R and run model 
output<- jags(data, inits, parameters, "model.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb)


# # Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 

# Summarize posteriors
print(output, dig = 3)
#in this case, logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
#check convergence with Brooks-Gelman-Rubin statistic
which(output$BUGSoutput$summary[, c("Rhat")] > 1.1)


#save output 
saveRDS(output, "Data/output/model_output_lmb_3level.rds") 

jag.sum<-output$BUGSoutput$summary
write.table(x=jag.sum,file="out.txt",sep="\t")

#### plots #### 
library(dotwhisker)

#betas
bEst <- matrix(NA, nrow=10,ncol=3)
for(i in 1:10){ #parameters
  bEst[i,1] <- mean(output$BUGSoutput$sims.list$b[,i])
  bEst[i,2:3] <- quantile(output$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
bEst<-as.data.frame(bEst)
bEst$variable<-c("Secchi", "lake area", "dd_year", "max_depth", "bottom_do", 
                 "ws forest", "ws wetland",  "julian day", "pres_wae", "pres_pike" )

colnames(bEst)<- c("estimate", "conf.low", "conf.high", "term")
bEst$Parameter <-c("b[1]", "b[2]", "b[3]", "b[4]","b[5]","b[6]","b[7]","b[8]", "b[9]" , "b[10]")


#add colors for sig different than 0
bEst$color <- as.numeric(bEst[,2] * bEst[,3] > 0 )

# plot using the dotwhisker package 
lake_plot<-dwplot(bEst, style = "dotwhisker", 
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

#PLOT WITH DISTRIBUTIONS 
jagsfit.mcmc <- as.mcmc(output2)
model1tranformed <- ggs(jagsfit.mcmc) 

model1tranformed %>%
  filter(grepl("b",Parameter))%>%
  left_join(bEst) %>% #get parameter names and colors 
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.05, .95)) + 
  theme_bw() 
