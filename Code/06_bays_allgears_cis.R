# CISCO SDM 
# modeling hierarchical model with varying intercepts by FMU and all gears combined 
#written by Katelyn King and Henrique Giacomini
#code adapted from Tyler Wagner 
#separating catches from different gears
#removed the lake grouping 

#load libraries 
library(R2jags)
library(lme4)
library(MCMCpack)
library(dotwhisker) #https://cran.r-project.org/web/packages/dotwhisker/vignettes/dotwhisker-vignette.html
library(dplyr)
library(ggmcmc) # package for analyzing output 

#### data ####
TDO3<-read.csv("Data/MI_data/TDO3_MI.csv")
dat <-cis_dat_for_model %>% 
  left_join(TDO3)

## standardize and transform predictor values ##
dat$z_lake_area<-as.numeric(scale(log(dat$lake_area_m2))) #good 
dat$z_geom<-as.numeric(scale(log(dat$geom_ratio))) #good 
dat$z_dd_mean<-as.numeric(scale(log(dat$dd_mean))) #good
dat$z_surface_temp<-as.numeric(scale(log(dat$surface_temp_c))) #good 
dat$z_secchi<-as.numeric(scale(log(dat$secchi_m)))
dat$z_bottom_do<-as.numeric(scale(log(dat$bottom_do_mgl+ 0.001))) ##has 0s so added 0.001 
dat$z_julian<-as.numeric(scale(dat$julian)) ##standardize julian date 
dat$z_ws_agriculture<-as.numeric(scale(asin(sqrt(dat$ws_agriculture_prop)))) #ok
dat$z_ws_forest<-as.numeric(scale(asin(sqrt(dat$ws_forest_prop)))) #good
dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat$ws_wetland_prop)))) #good
dat$z_tdo3<-as.numeric(scale(log(dat$TDO3)))

dat$logeffort <- log(dat$effort_new)

#check out response: log of count is mostly normally dist, so Poisson seems appropriate 
hist(dat$fish_count_new)
hist(log(dat$fish_count_new))

dat<-dplyr::select(dat, new_key, fish_count_new, logeffort, gear2, FMU_Code, 
                   z_lake_area, z_geom, z_secchi, z_julian, 
                   z_dd_mean, z_bottom_do, z_tdo3,
                   z_ws_forest, z_ws_wetland) %>%
  na.omit()

n_distinct(dat$new_key) #157 lakes (getting rid of depth gives more lakes)
#check for collinearity 
my_data <- dat[, c(7:15)] 
PerformanceAnalytics::chart.Correlation(my_data, histogram=TRUE, pch=19)
#lake area and watershed - keep lake area 
#max depth and bottom temp - keep max depth 
#urban and forest - keep forest. 
#note that surface temp and DD are only 0.5 corr 
#down to 12 variables 
#dat<-dat%>% dplyr::select(-c(z_ws_area_km2, z_ws_urban, z_bottom_temp, z_ws_shrub))

#### intercept only model to chec ICC ####
#determine if there is significant region-scale variation 
m_region<-lmer(log(fish_count_new+1) ~ 1 + (1|FMU_Code), dat, REML=FALSE) #use maximum likelihood instead of resricted maximum likelihod
summary(m_region) # the between region variance is the sd of the grouping variable and the 'residual variance' is the within region variance  
icc_n <- as.data.frame(VarCorr(m_region))[1,4]
icc_d <- as.data.frame(VarCorr(m_region))[1,4] + 
  as.data.frame(VarCorr(m_region))[2,4]
icc_n / icc_d #variance among regions 0%  

#################################################################
########## BUGS CODE ############################################
#################################################################

#model - not hierarchical 
sink("model.txt")
cat("
    model {
    # Likelihood: 
    
    for (i in 1:n){         
    y[i] ~ dpois(lambda[i])  # Distribution Poisson  
    log(lambda[i]) <- log.lambda[i] # log-link 
    log.lambda[i] <- alpha + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i] + b[4] * x4[i] + 
                       b[5] * x5[i] + b[6] * x6[i] + b[7] * x7[i] + b[8] * x8[i] + b[9] * x9[i] + 
                      logq[gear[i]]*IND[i] + logeffort[i]
         
    } 
    
   # Priors
    alpha ~ dnorm(0, 0.0001)
    
    # priors for predictors 
    for(k in 1:9){
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
# Initial values assign for random variables (FMU there are 8 regions) 
inits <- function (){
  list (alpha = rnorm(1))
}


# Parameters monitored
parameters <- c("alpha", "b", 'logq')


# MCMC settings
ni <- 10000 # number of iterations 
nt <- 1 #number to thin
nb <- 5000 #number to burn  
nc <- 3 #number of chains


#### set up data #### 

#set number of different sampling gears used and set a gear index 
ngears<-length(unique(dat$gear2))
gear <- as.numeric(as.factor(dat$gear2))

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
dat$IND<-ifelse(dat$gear2 == "FT_NET", 0, 1) #FT_NET is the reference gear

# Load data
data <- list(y = dat$fish_count_new, gear=gear, n = dim(dat)[1], ngears = ngears,
             x1=dat$z_secchi, x2=dat$z_lake_area, x3=dat$z_dd_mean, x4=dat$z_geom, x5=dat$z_julian, x6=dat$z_bottom_do,
             x7=dat$z_ws_forest, x8=dat$z_ws_wetland, x9=dat$z_tdo3,
             logeffort=dat$logeffort, IND=dat$IND
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
#jag.sum<-output$BUGSoutput$summary
#write.table(x=jag.sum,file="out.txt",sep="\t")

### check convergence ####
#in this case, logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
#check convergence with Brooks-Gelman-Rubin statistic
which(output$BUGSoutput$summary[, c("Rhat")] > 1.1)
#look at trace plots 
#R2jags::traceplot(output)

saveRDS(output, "Data/output/model_output_cis.rds") 

#read in output 
#output <- readRDS("Data/output/model_output.rds") 
#density plots
# use as.mcmmc to convert rjags object into mcmc.list - plots in coda
jagsfit.mcmc <- as.mcmc(output)
require(lattice)
densityplot(jagsfit.mcmc)

#the ggs function transforms the mcmc output into a longformat tibble, that we can use to make different types of plots.
model1tranformed <- ggs(jagsfit.mcmc) 
summary(model1tranformed$Parameter)

#traceplots (caterpillars)
ggs_traceplot(model1tranformed, family = "alpha")
ggs_traceplot(model1tranformed, family = "b")
ggs_traceplot(model1tranformed, family = "logq")

# posterior density plots 
ggplot(filter(model1tranformed,Parameter == "b[1]", Iteration > 1000),aes(x = value))+
  geom_density(fill  = "yellow", alpha = .5)+
  geom_vline(xintercept = 0, col  = "red", size = 1)+ 
  theme_light() +
  labs(title = "Posterior Density of b1 secchi")


#### effect plots #### 
library(dotwhisker)

#betas
bEst <- matrix(NA, nrow=8,ncol=3)
for(i in 1:8){ #parameters
  bEst[i,1] <- mean(output$BUGSoutput$sims.list$b[,i])
  bEst[i,2:3] <- quantile(output$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
bEst<-as.data.frame(bEst)
bEst$variable<-c("Secchi", "lake area", "surface temp", "shoreline complexity", "julian day", "bottom_do", 
                 "ws agriculture", "ws wetland")


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
colnames(bEst)<- c("estimate", "conf.low", "conf.high", "term")

#add colors for sig different than 0
EstsLake$color <- as.numeric(EstsLake[,2] * EstsLake[,3] > 0 )
bEst$color <- as.numeric(bEst[,2] * bEst[,3] > 0 )
bEst$Parameter <-c("b[1]", "b[2]", "b[3]", "b[4]","b[5]","b[6]","b[7]","b[8]")

# plot using the dotwhisker package 
dwplot(bEst, style = "dotwhisker", 
       dot_args = list(aes(colour = factor(color))),
       whisker_args = list(aes(colour = factor(color))), 
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("Michigan LMB") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none") 

#another plot with the distribution 
model1tranformed %>%
  filter(grepl("b",Parameter))%>%
  left_join(bEst) %>% #get parameter names and colors 
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  stat_halfeye(.width = c(.05, .95)) + 
  theme_bw() 


### EXAMPLE FROM PREDICTIONS WHEN READY TO HINDCAST #### 

for (i in 1:N) { #contemporary data
  y[i] ~ dnorm(beta0 + beta1 * x[i], sigma) #e.g. x is DD 
}
for (j in 1:Nnew) { # historical data 
  ynew[j] ~ dnorm(beta0 + beta1 * xnew[j], sigma) #e.g. xnew is hindcasted DD 
}

#where y, x and xnew are data vectors and ynew is a variable for storing predictions. 
#What you get is a distribution of values that are plausible given your estimated model. 
#Since the model is probabilistic, the prediction is also probabilistic, i.e. we get the whole distribution of possible ynew values.
#For point-values take the average of ynew

#then compare ynew to actual abundance values 
