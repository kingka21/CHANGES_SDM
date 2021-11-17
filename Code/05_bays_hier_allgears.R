# modeling hierarchical model with varying intercepts and all gears combined
#written by Katelyn King and Henrique Giacomini
#code adapted from Tyler Wagner 
#using total catch as response, use Poisson distribution 
# there is the option of separating catches from different gears, would add another level in the model (the site, with multiple catch values, as a random factor). 


#load libraries 
library(R2jags)
library(lme4)
library(MCMCpack)
library(dotwhisker) #https://cran.r-project.org/web/packages/dotwhisker/vignettes/dotwhisker-vignette.html
library(dplyr)

#### data ####
#lmb_count_dat has count, effort, and cpue
#add in the cpue ratios
#the ratio as a fixed value
lmb_count_dat$gill_shock_ratio<-lmb_count_dat$GILL_cpue/lmb_count_dat$SHOCK_cpue
lmb_count_dat$gill_ft_ratio<-lmb_count_dat$GILL_cpue/lmb_count_dat$FT_NET_cpue
lmb_count_dat$gill_seine_ratio<-lmb_count_dat$GILL_cpue/lmb_count_dat$SEINE_cpue
lmb_count_dat$shock_ft_ratio<-lmb_count_dat$SHOCK_cpue/lmb_count_dat$FT_NET_cpue
lmb_count_dat$shock_seine_ratio<-lmb_count_dat$SHOCK_cpue/lmb_count_dat$SEINE_cpue
lmb_count_dat$ft_seine_ratio<-lmb_count_dat$FT_NET_cpue/lmb_count_dat$SEINE_cpue

#add up all of the count data into one column as the response 
lmb_count_dat$count<-lmb_count_dat$fish_count_FT_NET + lmb_count_dat$fish_count_SHOCK + lmb_count_dat$fish_count_GILL + lmb_count_dat$fish_count_SEINE

## standardize and transform predictor values ##
lmb_count_dat$z_lake_area<-as.numeric(scale(log(lmb_count_dat$lake_area_m2))) #good 
lmb_count_dat$z_dd_mean<-as.numeric(scale(log(lmb_count_dat$dd_mean))) #good
lmb_count_dat$z_secchi<-as.numeric(scale(log(lmb_count_dat$secchi_m)))

dat<-dplyr::select(lmb_count_dat, new_key, count, FMU_Code, z_lake_area, z_dd_mean, z_secchi, gill_shock_ratio, gill_ft_ratio, gill_seine_ratio, shock_ft_ratio, shock_seine_ratio) %>%
  na.omit()

#put Inf values to 0 
is.na(dat)<-sapply(dat, is.infinite)
dat[is.na(dat)]<-0

#check out response: log of count is actually normally dist. 
hist(dat$count)
hist(log(dat$count))

#################################################################
########## BUGS CODE ############################################
#################################################################
#hierarchical model with varying intercepts 
sink("model.txt")
cat("
    model {
    # Likelihood: 
    ngears<-4 #total number of different sampling gears used
    for (i in 1:n){ 
    y[i] ~ dpois(lambda[i])  # Distribution Poisson  
    log(lambda[i]) <- log.lambda[i] # log-link 
    log.lambda[i] <- alpha[group[i]] + b[1] * x1[i] + b[2] * x2[i] + b[3] * x3[i] + QE #QE is the combined effect of all the gears
    QE<-log(sum(c(q,1)*effort[i,]))  #the vector with all catchability values is c(q,1), which adds the reference catchability =1 for the last gear type
          #q is a vector with catchabilities (to be estimated) of all gears except the last type (for which catchability is fixed at 1).
          #effort is a matrix with sampling effort per site (rows) and gear type (columns)  
           
    } 
    
    # Level-2 of the model 
    for(j in 1:J){
    alpha[j] ~ dnorm(mu.alpha,tau.alpha)
    }
    
    # Priors
    mu.alpha ~ dnorm(0, 0.0001)
    sigma.alpha ~ dunif(0,100)
    tau.alpha <- pow(sigma.alpha,-2)
    
    # priors for predictors 
    for(k in 1:3){
    b[k] ~ dnorm(0,3)
    }
    
     # priors for gear-specific catchabilities 
    for(k in 1:(ngears-1)){
      q[k] ~ dnorm(mu.q[k],tau.q[k]) #elements of vectors muq and tauq have to be specified below
    }
    mu.q<-c(1,1,1,1) #expected values of catchabilities relative to last gear type (1 means equal catchability), from independent studies 
    tau.q<-c(0,0,0,0) #precision of priors for relative catchabilities, to be filled with actual values
   #ratios would be useful to estimate prior distributions for those catchabilities, relative to a reference gear.

    } # end model
    ",fill = TRUE)
sink()

### set parameters 
# Initial values
inits <- function (){
  list (mu.alpha = rnorm(1), sigma.alpha=runif(1) )
}


# Parameters monitored
parameters <- c("mu.alpha", "sigma.alpha","b")


# MCMC settings
ni <- 10000 # number of iterations 
nt <- 1 #number to thin
nb <- 5000 #number to burn  
nc <- 3 #number of chains


#### set up data #### 

# Set the number of FMUs 
J <- length(unique(dat$FMU_Code))

# Create group index, must go from 1,..J
dat$FMU_Code <- droplevels(as.factor(dat$FMU_Code))
group <- as.numeric(as.factor(dat$FMU_Code))

# Load data
data <- list(y = dat$count, group = group, n = dim(dat)[1], J = J,
             x1=dat$z_secchi, x2=dat$z_lake_area, x3=dat$z_dd_mean, 
             x4=dat$gill_shock_ratio, x5=dat$gill_ft_ratio, x6=dat$gill_seine_ratio, x7=dat$shock_ft_ratio, x8=dat$shock_seine_ratio
) 


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

#### Posterior means and CIs for all parameters ####
betaEsts <- matrix(NA, nrow=8,ncol=3)
for(i in 1:8){
  betaEsts[i,1] <- mean(quantile(output$BUGSoutput$sims.list$b[,i],c(0.05,0.95)))
  betaEsts[i, 2:3] <- quantile(output$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
betaEsts<-as.data.frame(betaEsts)


colnames(betaEsts) <- c("estimate", "conf.low", "conf.high")
betaEsts$term <- c("Secchi", "lake_area", "dd_mean", 
                         "gill_shock_ratio", "gill_ft_ratio", "gill_seine_ratio", "shock_ft_ratio", "shock_seine_ratio")

#add colors for sig different than 0
betaEsts$color <- as.numeric(betaEsts[,2] * betaEsts[,3] > 0 )

# plot using the dotwhisker package 
beta_plot<-dwplot(betaEsts, style = "dotwhisker", 
                   dot_args = list(aes(colour = factor(color))),
                   whisker_args = list(aes(colour = factor(color))), 
                   vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
beta_plot


