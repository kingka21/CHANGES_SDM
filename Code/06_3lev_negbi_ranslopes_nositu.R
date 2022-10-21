## Negative Binomial dist including random effect of gear #### 
#no correlations or dependencies of betas across gears 
#load libraries 
library(R2jags)
library(lme4)
library(MCMCpack)
library(dotwhisker) #https://cran.r-project.org/web/packages/dotwhisker/vignettes/dotwhisker-vignette.html
library(dplyr)
library(ggmcmc) # package for analyzing output 
library(tidybayes)
library(splitTools)

#### data ####

#lmb_dat_for_model has count by gear and effort by gear
dat<-read.csv("Data/lmb_dat_for_model_aug30.csv") 
surface_temp<-read.csv("Data/MI_data/lake_surface_temp.csv")%>%
  group_by(IHDLKID) %>%
  slice_max(HECTARES) %>% 
  rename(nhdid = IHDLKID)%>% 
  pivot_longer(
    cols= starts_with("TAVE_"),
    names_to = "year", 
    names_prefix = "TAVE_",
    values_to = "surf_temp_year") %>% 
  mutate(year = as.integer(as.character(year)))

dd_temp<-read.csv("Data/MI_data/lake_degree_days_year.csv")%>%
  group_by(IHDLKID) %>%
  slice_max(HECTARES) %>% 
  rename(nhdid = IHDLKID)%>% 
  pivot_longer(
    cols= starts_with("DD_"),
    names_to = "year", 
    names_prefix = "DD_",
    values_to = "dd_year")%>% 
  mutate(year = as.integer(as.character(year)))

dat<-left_join(dat, dd_temp, by = c("nhdid", "year")) %>% 
  left_join(surface_temp, by = c("nhdid", "year"))

## standardize and transform predictor values ##
#dat$z_order<-as.numeric(scale(log(dat$lake_order+0.001))) # min is 0 because these are the isolated
dat$z_lake_area<-as.numeric(scale(log(dat$lake_area_m2))) #good 
dat$z_max_depth<-as.numeric(scale(log(dat$maxdepth_m))) #good 
dat$z_dd_year<-as.numeric(scale(log(dat$dd_year))) #good
dat$z_surface_temp_year<-as.numeric(scale(log(dat$surf_temp_year))) #good 
dat$z_secchi<-as.numeric(scale(log(dat$secchi_m)))
dat$z_bottom_do<-as.numeric(scale(log(dat$bottom_do_mgl+ 0.001))) ##has 0s so added 0.001 
dat$z_doy<-as.numeric(scale(dat$day_of_year)) ##standardize julian date 
dat$z_ws_forest<-as.numeric(scale(asin(sqrt(dat$ws_forest_prop)))) #good
dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(dat$ws_wetland_prop)))) #good

dat$logeffort <- log(dat$effort_new)

#not using wae/pike or DO but use all available Secchi! 
dat<-dplyr::select(dat, new_key, fish_count_new, logeffort, gear2, FMU_Code, 
                   z_lake_area, z_max_depth, z_secchi, z_doy, 
                  z_surface_temp_year,
                   z_ws_forest, z_ws_wetland) %>%
  na.omit()

n_distinct(dat$new_key)
#check out response: log of count is mostly normally dist, but begbinomial for dispersion 
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
    y[i] ~ dnegbin(p[i],r)
    p[i] <- r/(r + lambda[i])
    log(lambda[i]) <- log.lambda[i] # log-link 
    log.lambda[i] <- alpha[site[i]] + b1[gear[i]] * x1[i] + b2[gear[i]] * x2[i] + b3[gear[i]] * x3[i] + b4[gear[i]] * x4[i] + 
                       b5[gear[i]]*x5[i] + b6[gear[i]]*x6[i] + b7[gear[i]]*x7[i] +
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
    
    # varying betas by gear 
    for(j in 1:ngears){
    b1[j] ~ dnorm(mu.b1, tau.b1)
    b2[j] ~ dnorm(mu.b2, tau.b2)
    b3[j] ~ dnorm(mu.b3, tau.b3)
    b4[j] ~ dnorm(mu.b4, tau.b4)
    b5[j] ~ dnorm(mu.b5, tau.b5)
    b6[j] ~ dnorm(mu.b6, tau.b6)
    b7[j] ~ dnorm(mu.b7, tau.b7)
    }
    
    # Priors
    mu.alpha2 ~ dnorm(0, 0.0001)
    sigma.alpha ~ dunif(0,100)
    sigma.alpha2 ~ dunif(0,100)
    r ~ dunif(0,100)
    
    #derived quantities 
    tau.alpha <- pow(sigma.alpha,-2)
    tau.alpha2 <- pow(sigma.alpha2,-2)
    
    # priors for predictors 
     
    mu.b1 ~ dnorm(0,  0.0000001)   #non-informative
    mu.b2 ~ dnorm(0,  0.0000001)
    mu.b3 ~ dnorm(0,  0.0000001)
    mu.b4 ~ dnorm(0,  0.0000001)
    mu.b5 ~ dnorm(0,  0.0000001)
    mu.b6 ~ dnorm(0,  0.0000001)
    mu.b7 ~ dnorm(0,  0.0000001)
  

    
    sigma.b1 ~ dunif(0,100)
    sigma.b2 ~ dunif(0,100)
    sigma.b3 ~ dunif(0,100)
    sigma.b4 ~ dunif(0,100)
    sigma.b5 ~ dunif(0,100)
    sigma.b6 ~ dunif(0,100)
    sigma.b7 ~ dunif(0,100)

  tau.b1 <- pow(sigma.b1,-2)
  tau.b2 <- pow(sigma.b2,-2)
  tau.b3 <- pow(sigma.b3,-2)
  tau.b4 <- pow(sigma.b4,-2)
  tau.b5 <- pow(sigma.b5,-2)
  tau.b6 <- pow(sigma.b6,-2)
  tau.b7 <- pow(sigma.b6,-2)

    # priors for gear-specific log-catchabilities 
    for (k in 1:ngears) { 
    logq[k] ~ dnorm(0,  0.0000001)   #non-informative
    } 
    
                      
    } # end model
    ",fill = TRUE)
sink()


### set parameters ####

# Initial values #8 is the number of regions, 4 is num of gears 
inits <- function (){
  list (mu.alpha = rnorm(8), b1=rnorm(4), b2=rnorm(4), b3=rnorm(4), b4=rnorm(4), 
        b5=rnorm(4), b6=rnorm(4), b7=rnorm(4)
  )
}


# Parameters monitored
parameters <- c("alpha", "mu.alpha","tau.alpha", "mu.alpha2","tau.alpha2",'b1', "b2", "b3", "b4", "b5", 
                "b6", "b7", "mu.b1", "mu.b2", "mu.b3", "mu.b4", "mu.b5", "mu.b6", "mu.b7",
                'logq', 'r')


# MCMC settings 
ni <- 100000 # number of iterations 
nt <- 10 #number to thin 20
nb <- 30000 #number to burn   30000
nc <- 3 #number of chains

#### split data into training and validation sets for cross-fold validation #### 
part1<-partition(dat$new_key, p=c(train=0.8, test = 0.2), seed = 1,  type =c( "grouped")) 
train1<-dat[part1$train,] 
test1<-dat[part1$test,]

part2<-partition(dat$new_key, p=c(train=0.8, test = 0.2), seed = 2, type =c( "grouped"))
train2<-dat[part2$train,]
test2<-dat[part2$test,]

part3<-partition(dat$new_key, p=c(train=0.8, test = 0.2), seed = 3, type =c( "grouped"))
train3<-dat[part3$train,]
test3<-dat[part3$test,]

part4<-partition(dat$new_key, p=c(train=0.8, test = 0.2), seed = 4, type =c( "grouped"))
train4<-dat[part4$train,]
test4<-dat[part4$test,]

part5<-partition(dat$new_key, p=c(train=0.8, test = 0.2), seed = 5, type =c( "grouped"))
train5<-dat[part5$train,]
test5<-dat[part5$test,]

#### set up data #### 
train_cont<-train5 #just change this every time 

#set the number of sites and create site index 
nsites<-length(unique(train_cont$new_key))
site <- as.numeric(as.factor(train_cont$new_key))

# Set the number of FMUs 
J <- length(unique(train_cont$FMU_Code))

# Create group index, must go from 1,..J
train_cont$FMU_Code <- droplevels(as.factor(train_cont$FMU_Code))
train_cont$group <- as.numeric(as.factor(train_cont$FMU_Code))
group <- doBy::summaryBy(group ~ new_key, data=train_cont, FUN=mean) #need region to be indexed by site

#set number of different sampling gears used and set a gear index 
ngears<-length(unique(train_cont$gear2))
gear <- as.numeric(as.factor(train_cont$gear2))

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
train_cont$IND<-ifelse(train_cont$gear2 == "FT_NET", 0, 1) #FT_NET is the reference gear

# Load data
data <- list(y = train_cont$fish_count_new, group = group$group.mean, gear=gear, n = dim(train_cont)[1], J = J, ngears = ngears,
             x1=train_cont$z_secchi, x2=train_cont$z_lake_area, x3=train_cont$z_surface_temp_year, x4=train_cont$z_max_depth,
             x5=train_cont$z_ws_forest,x6=train_cont$z_ws_wetland, x7=train_cont$z_doy,
             logeffort=train_cont$logeffort, IND=train_cont$IND, nsites = nsites, site =site
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

#save output 
#saveRDS(output, "Data/output/output_model2_secchi_fold5_lakes.rds") 
output<-readRDS("Data/output/output_model2_secchi_fold3_lakes.rds")


#* Summarize posteriors ####
print(output, dig = 3)
#in this case, logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
#check convergence with Brooks-Gelman-Rubin statistic
which(output$BUGSoutput$summary[, c("Rhat")] > 1.1)
which(output$BUGSoutput$summary[, c("n.eff")] < 0.1*21000) #times number of retained draws 
which(output$BUGSoutput$summary[, c("n.eff")] < 500)
jagsfit.mcmc <- as.mcmc(output)
model1tranformed <- ggs(jagsfit.mcmc) 

ggs_traceplot(model1tranformed, family = "b7")
ggs_traceplot(model1tranformed, family = "b6")
ggs_traceplot(model1tranformed, family = "b5")
ggs_traceplot(model1tranformed, family = "b4")
ggs_traceplot(model1tranformed, family = "b3")
ggs_traceplot(model1tranformed, family = "b2")
ggs_traceplot(model1tranformed, family = "b1")

ggs_traceplot(model1tranformed, family = "r")
ggs_traceplot(model1tranformed, family = "logq")
ggs_traceplot(model1tranformed, family = "mu.b")

#### effect plots #### 
library(dotwhisker)

#betas
#1 is FT, 2 is gill, 3 is seine, 4 is shock 
#x1=z_secchi, x2=z_lake_area, x3=z_surface_temp_year, x4=z_max_depth,
#x5=z_ws_forest,x6=z_ws_wetland, x7=z_doy,
EstsSecchi <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  EstsSecchi[i,1] <- mean(output$BUGSoutput$sims.list$b1[,i])
  EstsSecchi[i,2:3] <- quantile(output$BUGSoutput$sims.list$b1[,i],c(0.05,0.95))
}
EstsSecchi<-as.data.frame(EstsSecchi)
colnames(EstsSecchi)<- c("estimate", "conf.low", "conf.high")
EstsSecchi$gear<-c("Fyke", "Gill", "Seine", "Shock")
EstsSecchi$term<-'Secchi'

Estsarea <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  Estsarea[i,1] <- mean(output$BUGSoutput$sims.list$b2[,i])
  Estsarea[i,2:3] <- quantile(output$BUGSoutput$sims.list$b2[,i],c(0.05,0.95))
}
Estsarea<-as.data.frame(Estsarea)
colnames(Estsarea)<- c("estimate", "conf.low", "conf.high")
Estsarea$gear<-c("Fyke", "Gill", "Seine", "Shock")
Estsarea$term<-'lake_area'

Eststemp <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  Eststemp[i,1] <- mean(output$BUGSoutput$sims.list$b3[,i])
  Eststemp[i,2:3] <- quantile(output$BUGSoutput$sims.list$b3[,i],c(0.05,0.95))
}
Eststemp<-as.data.frame(Eststemp)
colnames(Eststemp)<- c("estimate", "conf.low", "conf.high")
Eststemp$gear<-c("Fyke", "Gill", "Seine", "Shock")
Eststemp$term<-"surface_temp"

Estsdepth <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  Estsdepth[i,1] <- mean(output$BUGSoutput$sims.list$b4[,i])
  Estsdepth[i,2:3] <- quantile(output$BUGSoutput$sims.list$b4[,i],c(0.05,0.95))
}
Estsdepth<-as.data.frame(Estsdepth)
colnames(Estsdepth)<- c("estimate", "conf.low", "conf.high")
Estsdepth$gear<-c("Fyke", "Gill", "Seine", "Shock")
Estsdepth$term<-"max_depth"

EstsFor <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  EstsFor[i,1] <- mean(output$BUGSoutput$sims.list$b5[,i])
  EstsFor[i,2:3] <- quantile(output$BUGSoutput$sims.list$b5[,i],c(0.05,0.95))
}
EstsFor<-as.data.frame(EstsFor)
colnames(EstsFor)<- c("estimate", "conf.low", "conf.high")
EstsFor$gear<-c("Fyke", "Gill", "Seine", "Shock")
EstsFor$term<-'ws_forest'

Estswet <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  Estswet[i,1] <- mean(output$BUGSoutput$sims.list$b6[,i])
  Estswet[i,2:3] <- quantile(output$BUGSoutput$sims.list$b6[,i],c(0.05,0.95))
}
Estswet<-as.data.frame(Estswet)
colnames(Estswet)<- c("estimate", "conf.low", "conf.high")
Estswet$gear<-c("Fyke", "Gill", "Seine", "Shock")
Estswet$term<-'ws_wetland'

EstsDOY <- matrix(NA, nrow=4,ncol=3)
for(i in 1:4){ #gears
  EstsDOY[i,1] <- mean(output$BUGSoutput$sims.list$b7[,i])
  EstsDOY[i,2:3] <- quantile(output$BUGSoutput$sims.list$b7[,i],c(0.05,0.95))
}
EstsDOY<-as.data.frame(EstsDOY)
colnames(EstsDOY)<- c("estimate", "conf.low", "conf.high")
EstsDOY$gear<-c("Fyke", "Gill", "Seine", "Shock")
EstsDOY$term<-'day-of-year'



all_est<- gtools::smartbind(EstsSecchi, Estsarea, Eststemp, Estsdepth, EstsFor, Estswet, EstsDOY)
#add colors for sig different than 0
all_est$color <- as.numeric(all_est[,2] * all_est[,3] > 0 )
all_est$Parameter <-c("b1[1]", "b1[2]", "b1[3]", "b1[4]",
                      "b2[1]", "b2[2]", "b2[3]", "b2[4]",
                      "b3[1]", "b3[2]", "b3[3]", "b3[4]",
                      "b4[1]", "b4[2]", "b4[3]", "b4[4]",
                      "b5[1]", "b5[2]", "b5[3]", "b5[4]",
                      "b6[1]", "b6[2]", "b6[3]", "b6[4]", 
                      "b7[1]", "b7[2]", "b7[3]", "b7[4]")

#split into the gears for graphing 
fyke<-filter(all_est, gear=="Fyke")
gill<-filter(all_est, gear=="Gill")
seine<-filter(all_est, gear=="Seine")
shock<-filter(all_est, gear=="Shock")



# plot using the dotwhisker package 
fyke_plot<-dwplot(fyke, style = "dotwhisker", 
                  dot_args = list(aes(colour = factor(color))),
                  whisker_args = list(aes(colour = factor(color))), 
                  vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("fyke") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
fyke_plot

gill_plot<-dwplot(gill, style = "dotwhisker", 
                  dot_args = list(aes(colour = factor(color))),
                  whisker_args = list(aes(colour = factor(color))), 
                  vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("gill") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
gill_plot

seine_plot<-dwplot(seine, style = "dotwhisker", 
                   dot_args = list(aes(colour = factor(color))),
                   whisker_args = list(aes(colour = factor(color))), 
                   vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("seine") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
seine_plot

shock_plot<-dwplot(shock, style = "dotwhisker", 
                   dot_args = list(aes(colour = factor(color))),
                   whisker_args = list(aes(colour = factor(color))), 
                   vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Covariate") +
  ggtitle("shock") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
shock_plot

#####PLOT WITH DISTRIBUTIONS ####
#1 is FT, 2 is gill, 3 is seine, 4 is shock 
fyke_dist<-model1tranformed %>%
  filter(grepl("b.*[1]",Parameter))%>% #use .* as a wild card to capture all 10 params
  left_join(fyke) %>% #get parameter names and colors 
  na.omit(term) %>%
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.05, .95)) + 
  theme_bw() +
  ggtitle("fyke net") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none") +
  ylab("predictor")

gill_dist<-model1tranformed %>%
  filter(grepl("b.*[2]",Parameter))%>% #use .* as a wild card to capture all 10 params
  left_join(gill) %>% #get parameter names and colors 
  na.omit(term) %>%
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.05, .95)) + 
  theme_bw() +
  ggtitle("gill net") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none") +
  ylab("predictor")


seine_dist<-model1tranformed %>%
  filter(grepl("b.*[3]",Parameter))%>% #use .* as a wild card to capture all 10 params
  left_join(seine) %>% #get parameter names and colors 
  na.omit(term) %>%
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.10, .90)) + 
  theme_bw() +
  ggtitle("seine") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none") +
  ylab("predictor")


shock_dist<-model1tranformed %>%
  filter(grepl("b.*[4]",Parameter))%>% #use .* as a wild card to capture all 10 params
  left_join(shock) %>% #get parameter names and colors 
  na.omit(term) %>%
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.10, .90)) + 
  theme_bw() +
  ggtitle("shock") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none") +
  ylab("predictor")

#* save as 4-panel plot 
model2_covariates<-cowplot::plot_grid(fyke_dist, gill_dist, seine_dist, shock_dist, labels=c("a", "b", "c", "d"), ncol=2)

ggsave(plot=model2_covariates, 
       device = "png", 
       filename = "figures/model2_covariates.png", 
       dpi = 600, height = 8, width = 12, units = "in",
       bg="#ffffff") #sets background to white 


#### MODEL FIT #### 
# Simulating data from the posterior predictive distribution using the observed predictors is useful for checking the fit of the model.
#obtain our samples from the posterior
#this has a posterior dist (30000) for coef for every param; e.g. n lakes, 7x4 betas, 4 logqs, 8 mu.alphas, etc. 
coefs <- output$BUGSoutput$sims.matrix[,1:256] 

### 
# Number of desired MCMC samples used for prediction (use a subset of all samples)
nsim <- 3000
# Chain length from analysis
chainLength <- output$BUGSoutput$n.sims
# Select thinned steps in chain for posterior predictions to ensure we take values from length of posterior
ID = seq( 1 , chainLength , floor(chainLength/nsim) )

#* predict to test dataset ####
cont_test<-test3# just change this input each time 
n_distinct(cont_test$new_key)

#need to add management unit index and gear index
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
cont_test<-cont_test %>% 
  mutate(group = case_when( FMU_Code == "MI-MCM" ~ 1, 
                       FMU_Code == "MI-MER" ~ 2, 
                       FMU_Code == "MI-MES" ~ 3, 
                       FMU_Code == "MI-MNH" ~ 4,  
                       FMU_Code == "MI-MNM" ~ 5, 
                       FMU_Code== "MI-MSH" ~ 6,
                       FMU_Code== "MI-MSM" ~ 7, 
                       FMU_Code== "MI-MWS" ~ 8)
  )%>% 
  mutate(gear = case_when(gear2 == "FT_NET" ~ 1,
                            gear2 == "GILL" ~ 2,
                            gear2 == "SEINE" ~ 3, 
                            gear2 == "SHOCK" ~ 4)
  )

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
cont_test$IND<-ifelse(cont_test$gear2 == "FT_NET", 0, 1) #FT_NET is the reference gear

# Container for predicted values using the test dataset 5
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
#x1=z_secchi, x2=z_lake_area, x3=z_surface_temp_year, x4=z_max_depth,
#x5=z_ws_forest,x6=z_ws_wetland, x7=z_doy,
predictions_sim<- array(NA, c(nsim,length(cont_test[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predictions_sim)

for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(cont_test[[1]])){ #loop over covariate data (obs) #use regional mu.alphas instead of site specific
    predictions_sim[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',cont_test[['group']][t],']')] + coefs[ID[i], paste0('b1[',cont_test[['gear']][t],']')]*cont_test[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',cont_test[['gear']][t],']')]*cont_test[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',cont_test[['gear']][t],']')]*cont_test[['z_surface_temp_year']][t]  + 
                                  coefs[ID[i],paste0('b4[',cont_test[['gear']][t],']')]*cont_test[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',cont_test[['gear']][t],']')]*cont_test[['z_ws_forest']][t] + coefs[ID[i],paste0('b6[',cont_test[['gear']][t],']')]*cont_test[['z_ws_wetland']][t] + coefs[ID[i],paste0('b7[',cont_test[['gear']][t],']')]*cont_test[['z_doy']][t] +
                                  coefs[ID[i], paste0('logq[',cont_test[['gear']][t],']')]*cont_test[['IND']][t] + cont_test[['logeffort']][t]) # exp the predictions because we used a log-link in the negbi
    
  }
}


#container for p values 
#parametrization is by the dispersion parameter, where prob = size/(size+mu).
#p is the success parameter and r is the dispersion parameter.
p <- array(NA, c(nsim,length(cont_test[[1]])) )
exp_catch <- array(NA, c(nsim,length(cont_test[[1]])) )

for(i in 1:nrow(predictions_sim)){ # loop over rows (sims)
  for(t in 1:ncol(predictions_sim) ){ #loop over columns(obs)
    p[i,t] <- coefs[ID[i],'r']/(coefs[ID[i],'r'] + predictions_sim[i,t])
    exp_catch[i,t] <- rnbinom(1, prob=p[i,t], size= coefs[ID[i],'r'])
  }
}

#From this distribution, you can extract the median and the 2.5% and 97.5% percentiles (95% credible interval),
col.med <- apply(exp_catch, 2, median )
col.means <- apply(exp_catch, 2, mean )
# 95% CIs for fitted values
upperCI.Group <- apply(exp_catch, 2, quantile, probs=c(0.975) )
lowerCI.Group <- apply(exp_catch, 2, quantile, probs=c(0.025) )

#prep data for plotting
med_data=data.frame(col.med, col.means, upperCI.Group, lowerCI.Group) 
med_data$row <- as.numeric(row.names(med_data))

cont_test$site <- as.numeric(as.factor(cont_test$new_key))
obs_catch<-dplyr::select(cont_test, new_key, site, gear2, fish_count_new) %>% 
  mutate(row = row_number())
plot_data<-left_join(med_data, obs_catch) %>% 
  mutate(gear = case_when(       #rename gear types for publication
    gear2 == "FT_NET" ~ 'fyke net', 
    gear2 == "GILL" ~ 'gillnet', 
    gear2 == "SEINE" ~ 'seine', 
    gear2 == "SHOCK" ~ 'shock'))

plot( med_data$col.med , obs_catch$fish_count_new)
plot( med_data$col.means, obs_catch$fish_count_new)

#which can be compared to the observed catch for that combination of lake and gear. 
#One suggestion for visualizing the model fit is to generate a separate graph for each gear, 
#plot the median predicted catches (+- the 95% ci) of different lakes on the x-axis versus the observed catches of lakes on the y-axis.
pred_plot<-ggplot()+
  geom_pointrange(data=plot_data, aes(x=col.means, y=fish_count_new, xmax=upperCI.Group, xmin=lowerCI.Group))+ 
  xlab('predicted catch')+
  ylab('observed catch')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)

pred_plot+facet_wrap(~ gear2, ncol=2, scales = 'free') #allow scales to vary 

#plot log scale 
pred_plot_log<-ggplot(data=plot_data, aes(x=log(col.means+1), y=log(fish_count_new + 1), xmax=log(upperCI.Group +1), xmin=log(lowerCI.Group +1) ) )+
  geom_pointrange( color="grey" )+ 
  geom_point(color="black")+ 
  xlab('predicted catch (log)')+
  ylab('observed catch (log)')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)
pred_plot_log

model2_pred_obs<-pred_plot_log+facet_wrap(~ gear, ncol=2, scales = 'free') #allow scales to vary 


ggsave(plot=model2_pred_obs, 
       device = "png", 
       filename = "figures/model2_pred_obs.png", 
       dpi = 600, height = 8, width = 12, units = "in",
       bg="#ffffff") #sets background to white 

###  plot model fits  #### 
library(bayesplot)
#yrep<-as.matrix(output$BUGSoutput$sims.list$ysim ) #from the predictions within JAGS
yrep<-exp_catch
y<-cont_test$fish_count_new

color_scheme_set("brightblue")
ppc_dens_overlay(y, yrep[1:50, ])
ppc_dens_overlay(y, yrep[1:50, ]) + xlim(0, 50)
#ppc_hist(y, yrep[1:5, ])
#group 
model2_dens<-ppc_dens_overlay_grouped(y, yrep[1:50, ], group = test5$gear2) + xlim(0, 20)

ggsave(plot=model2_dens, 
       device = "png", 
       filename = "figures/model2_dens.png", 
       dpi = 600, height = 8, width = 12, units = "in",
       bg="#ffffff") 

#look at how well it predicts the 0 obs 
prop_zero <- function(x) mean(x == 0)
prop_zero(y) 
ppc_stat(y, yrep, stat = "prop_zero", binwidth = 0.005)

#*residuals for contemporary model  
yrep<-exp_catch
y<-cont_test$fish_count_new
#Bayes-R2 using residuals
e <- -1 * sweep(yrep, 2, y)
var_ypred <- apply(yrep, 1, var)
var_e <- apply(e, 1, var)
bayesR2res<- var_ypred / (var_ypred + var_e)
round(median(bayesR2res), 2) 

#RMSE 
sqrt(mean((plot_data$col.means-plot_data$fish_count_new)^2))

#residual plot
#https://cran.r-project.org/web/packages/tidybayes/vignettes/tidybayes-residuals.html
plot_data$res.med <- apply(e, 2, median )
res_plot<-plot_data %>%
  ggplot(aes(x = col.med, y = res.med)) +
  stat_pointinterval()
res_plot
res_plot+facet_wrap(~ gear2, ncol=2, scales = 'free')

plot_data  %>%
  ggplot(aes(sample = res.med)) +
  geom_qq() +
  geom_qq_line()

#### PREDICTION to historical ##### 
hist_dat<-read.csv("Data/historical_lmb2_with_secchi.csv") %>% 
  drop_na(secchi_combo) #only use samples with secchi from contemp data 

#fix julian day 
#hist_dat$julian<-ifelse(hist_dat$new_key == '65-63', 181, hist_dat$julian)

#set up data 
hist_dat$logeffort<-log(hist_dat$effort_sum)
## standardize and transform predictor values ##
hist_dat$z_lake_area<-as.numeric(scale(log(hist_dat$lake_area_m2))) #good 
hist_dat$z_max_depth<-as.numeric(scale(log(hist_dat$max_depth_m))) #good 
hist_dat$z_dd_year-as.numeric(scale(log(hist_dat$dd_year))) #good
hist_dat$z_surface_temp_year<-as.numeric(scale(log(hist_dat$surface_temp_year))) #good 
hist_dat$z_secchi<-as.numeric(scale(log(hist_dat$secchi_combo)))
hist_dat$z_bottom_do<-as.numeric(scale(log(hist_dat$bottom_do + 0.001))) ##has 0s so added 0.001 
hist_dat$z_julian<-as.numeric(scale(hist_dat$julian)) ##standardize julian date 
hist_dat$z_ws_forest<-as.numeric(scale(asin(sqrt(hist_dat$ws_forest_prop)))) #good
hist_dat$z_ws_wetland<-as.numeric(scale(asin(sqrt(hist_dat$ws_wetland_prop)))) #good

# need group number to match initial data 
#MCM =1 Central Lake Michigan
#MER = 2 Lake Erie
# MES = 3 Eastern Lake Superior
#MNH = 4 Northern Lake Huron
# MNM = 5 Northern Lake Michigan
#MSH = 6 Southern Lake Huron
#MSM = 7 Southern Lake Michigan
#MWS= 8 Western Lake Superior
hist_dat<-hist_dat %>% 
  mutate(
    group = case_when( MGMT_UNIT == "Central Lake Michigan" ~ 1,
                       MGMT_UNIT == "Lake Erie" ~ 2,
                       MGMT_UNIT == "Eastern Lake Superior" ~ 3, 
                       MGMT_UNIT == "Northern Lake Huron" ~ 4,
                       MGMT_UNIT == "Northern Lake Michigan" ~ 5,
                       MGMT_UNIT== "Southern Lake Huron" ~ 6,
                       MGMT_UNIT== "Southern Lake Michigan" ~ 7,
                       MGMT_UNIT== "Western Lake Superior" ~ 8)
  )

#need gear number to match original data (no shock) 
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
hist_dat<-hist_dat %>% 
  mutate(
    gear_index = case_when( gear == "FT_net" ~ 1,
                            gear == "gill" ~ 2,
                            gear == "seine" ~ 3)
  )

#create an indicator variable “IND”, with value 0 for every sample using the reference gear and value 1 otherwise
hist_dat$IND<-ifelse(hist_dat$gear == "FT_net", 0, 1) #FT_NET is the reference gear

#### *multiple test sets for each fold ####
#pull out about 50 lakes to match contemporary test dataset 
hist1<-partition(hist_dat$new_key, p=c(train=0.4,test=0.6), seed = 1, type =c("grouped"))
hist_test1<-hist_dat[hist1$test,]
hist2<-partition(hist_dat$new_key, p=c(train=0.4,test=0.6), seed = 2, type =c("grouped"))
hist_test2<-hist_dat[hist2$test,]
hist3<-partition(hist_dat$new_key, p=c(train=0.4,test=0.6), seed = 3, type =c("grouped"))
hist_test3<-hist_dat[hist3$test,]
hist4<-partition(hist_dat$new_key, p=c(train=0.4,test=0.6), seed = 4, type =c("grouped"))
hist_test4<-hist_dat[hist4$test,]
hist5<-partition(hist_dat$new_key, p=c(train=0.4,test=0.6), seed = 5, type =c("grouped"))
hist_test5<-hist_dat[hist5$test,]

hist_test<-hist_test1# just change this input each time 
n_distinct(hist_test$new_key)

# Container for predicted values
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
predict_hist<- array(NA, c(nsim,length(hist_test[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predict_hist)

#use regional mu.alphas instead of site specific, since dif sites (this takes like 5 mins)
for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(hist_test[[1]])){ #loop over covariate data (obs)
    predict_hist[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',hist_test[['group']][t],']')] + coefs[ID[i],paste0('b1[',hist_test[['gear_index']][t],']')]*hist_test[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',hist_test[['gear_index']][t],']')]*hist_test[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',hist_test[['gear_index']][t],']')]*hist_test[['z_surface_temp_year']][t]  + 
                               coefs[ID[i],paste0('b4[',hist_test[['gear_index']][t],']')]*hist_test[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',hist_test[['gear_index']][t],']')]*hist_test[['z_ws_forest']][t] + coefs[ID[i],paste0('b6[',hist_test[['gear_index']][t],']')]*hist_test[['z_ws_wetland']][t] +
                               coefs[ID[i],paste0('b7[',hist_test[['gear_index']][t],']')]*hist_test[['z_julian']][t] + 
                               coefs[ID[i], paste0('logq[',hist_test[['gear_index']][t],']')]*hist_test[['IND']][t] + hist_test[['logeffort']][t]) # exp the predictions because we used a log-link in the Poisson
    
  }
}


#container for p values 
#parametrization is by the dispersion parameter, where prob = size/(size+mu).
#p is the success parameter and r is the dispersion parameter.
p_hist <- array(NA, c(nsim,length(hist_test[[1]])) )
exp_catch_hist <- array(NA, c(nsim,length(hist_test[[1]])) )

for(i in 1:nrow(predict_hist)){ # loop over rows (sims)
  for(t in 1:ncol(predict_hist) ){ #loop over columns(obs)
    p_hist[i,t] <- coefs[ID[i],'r']/(coefs[ID[i],'r'] + predict_hist[i,t])
    exp_catch_hist[i,t] <- rnbinom(1, prob=p_hist[i,t], size= coefs[ID[i],'r'])
  }
}

#From this distribution, you can extract the median and the 2.5% and 97.5% percentiles (95% credible interval),
hist.med <- apply(exp_catch_hist, 2, median )
hist.means <- apply(exp_catch_hist, 2, mean )
# 95% CIs for fitted values
upperPI <- apply(exp_catch_hist, 2, quantile, probs=c(0.975) )
lowerPI <- apply(exp_catch_hist, 2, quantile, probs=c(0.025) )

#prep data for plotting
med_data=data.frame(hist.med, hist.means, upperPI, lowerPI) 
med_data$row <- as.numeric(row.names(med_data))

hist_test$site <- as.numeric(as.factor(hist_test$new_key))
hist_catch<-dplyr::select(hist_test, new_key, site, gear, largemouthbass_sum) %>% 
  mutate(row = row_number())
plot_data<-left_join(med_data, hist_catch) 

plot( med_data$hist.med , hist_catch$largemouthbass_sum)
plot( med_data$hist.means, hist_catch$largemouthbass_sum)

ypred<-exp_catch_hist
y<-plot_data$largemouthbass_sum
e <- -1 * sweep(ypred, 2, y) #residual distribution 
var_ypred <- apply(ypred, 1, var) # variance of modelled predictive means
var_e <- apply(e, 1, var) #variance of modeled residuals 
bayesR2res<- var_ypred / (var_ypred + var_e) #R^2
round(median(bayesR2res), 2) #0.36
sqrt(mean((plot_data$hist.means-plot_data$largemouthbass_sum)^2))


#remove fyke before plotting 
plot_data_nofyke=filter(plot_data, gear != 'FT_net')
#which can be compared to the observed catch for that combination of lake and gear. 
#One suggestion for visualizing the model fit is to generate a separate graph for each gear, 
#plot the median predicted catches (+- the 95% ci) of different lakes on the x-axis versus the observed catches of lakes on the y-axis.
pred_plot<-ggplot()+
  geom_pointrange(data=plot_data_nofyke, aes(x=hist.means, y=largemouthbass_sum, xmax=upperPI, xmin=lowerPI))+ 
  xlab('predicted catch')+
  ylab('historical catch')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)

pred_plot+facet_wrap(~ gear, ncol=2, scales = 'free') #allow scales to vary 

#plot log scale 
pred_plot_log<-ggplot(data=plot_data_nofyke, aes(x=log(hist.means+1), y=log(largemouthbass_sum + 1), xmax=log(upperPI +1), xmin=log(lowerPI +1) ))+
  geom_pointrange(color = "grey")+ 
  geom_point(color="black")+ 
  xlab('predicted catch (log)')+
  ylab('historical catch (log)')+
  theme_bw()+theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text = element_text(size=14)) + 
  geom_abline(intercept = 0, slope = 1)

model2_hist_pred<-pred_plot_log+facet_wrap(~ gear, ncol=2, scales = 'free') #allow scales to vary 
model2_hist_pred

ggsave(plot=model2_hist_pred, 
       device = "png", 
       filename = "figures/model2_hist_pred.png", 
       dpi = 600, height = 8, width = 12, units = "in",
       bg="#ffffff") #sets background to white


#### RESIDUALS #### 
#*Bayes-R2 using residuals
#https://avehtari.github.io/bayes_R2/bayes_R2.html#2_Functions_for_Bayesian_R-squared_for_stan_glm_models
ypred<-exp_catch_hist
y<-plot_data$largemouthbass_sum
e <- -1 * sweep(ypred, 2, y) #residual distribution 
var_ypred <- apply(ypred, 1, var) # variance of modelled predictive means
var_e <- apply(e, 1, var) #variance of modeled residuals 
bayesR2res<- var_ypred / (var_ypred + var_e) #R^2
round(median(bayesR2res), 2) #0.36

#* RMSE historical #### 
sqrt(mean((plot_data$hist.means-plot_data$largemouthbass_sum)^2))

#plot the R2 posterior and median
mcmc_hist(data.frame(bayesR2res), binwidth=0.02)+xlim(c(0,1))+
  xlab(('R^2'))+
  geom_vline(xintercept=median(bayesR2res))+
  ggtitle('Bayesian R squared posterior and median')

#plot fits 
color_scheme_set("brightblue")
ppc_dens_overlay(y, ypred[1:50, ])
ppc_dens_overlay(y, ypred[1:50, ]) + xlim(0, 50)

#group 
ppc_dens_overlay_grouped(y, ypred[1:50, ], group = hist_dat$gear) + xlim(0, 20)

#residual plot
#https://cran.r-project.org/web/packages/tidybayes/vignettes/tidybayes-residuals.html
plot_data$res.med <- apply(e, 2, median )
res_plot<-plot_data %>%
  ggplot(aes(x = hist.med, y = res.med)) +
  stat_pointinterval()

res_plot+facet_wrap(~ gear, ncol=2, scales = 'free')

plot_data  %>%
  ggplot(aes(sample = res.med)) +
  geom_qq() +
  geom_qq_line()

#### relative abundance hist ####
#index of density that is the sum of expected catches of all gears, using a common effort (E=1)
#extract just lake data (need a single row for each lake) 
hist_just_lake<-hist_dat %>% 
 distinct(new_key, .keep_all = TRUE) %>% #83 lakes 
  select(-c(gear, effort_sum, logeffort, gear_index, site)) %>% 
  mutate(effort = 1) 

dup<-hist_just_lake[rep(c(1:83),2),] # duplicate all rows 

hist_just_lake<-dup %>% 
  group_by(new_key) %>% 
  mutate(gear = if_else(row_number() %% 2 ==1, "seine", "gill")) %>% #add just seine and gill 
  ungroup %>% 
  mutate(gear_index = ifelse(gear == "gill", 2, 3))

# Container for predicted values
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
predict_hist<- array(NA, c(nsim,length(hist_just_lake[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predict_hist)

#use regional mu.alphas instead of site specific, since dif sites, use estimated q and effort = 1
for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(hist_just_lake[[1]])){ #loop over covariate data (obs)
    predict_hist[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',hist_just_lake[['group']][t],']')] + coefs[ID[i],paste0('b1[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_surface_temp_year']][t]  + 
                               coefs[ID[i],paste0('b4[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_ws_forest']][t] + coefs[ID[i],paste0('b6[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_ws_wetland']][t] +
                               coefs[ID[i],paste0('b7[',hist_just_lake[['gear_index']][t],']')]*hist_just_lake[['z_julian']][t] + 
                               coefs[ID[i], paste0('logq[',hist_just_lake[['gear_index']][t],']')] + hist_just_lake[['effort']][t]) # exp the predictions because we used a log-link in the Poisson
    
  }
}

#container for p values 
#parametrization is by the dispersion parameter, where prob = size/(size+mu).
#p is the success parameter and r is the dispersion parameter.
p_hist <- array(NA, c(nsim,length(hist_just_lake[[1]])) )
exp_abund_hist <- array(NA, c(nsim,length(hist_just_lake[[1]])) )

for(i in 1:nrow(predict_hist)){ # loop over rows (sims)
  for(t in 1:ncol(predict_hist) ){ #loop over columns(obs)
    p_hist[i,t] <- coefs[ID[i],'r']/(coefs[ID[i],'r'] + predict_hist[i,t])
    exp_abund_hist[i,t] <- rnbinom(1, prob=p_hist[i,t], size= coefs[ID[i],'r'])
  }
}

#From this distribution, you can extract the median and the 2.5% and 97.5% percentiles (95% credible interval),
hist.med <- apply(exp_abund_hist, 2, median )
hist.means <- apply(exp_abund_hist, 2, mean )

#prep data for plotting
abund_data=data.frame(hist.med, hist.means) 
abund_data$row <- as.numeric(row.names(abund_data))
abund_data$log_mean<-log(abund_data$hist.means)

hist_just_lake$site <- as.numeric(as.factor(hist_just_lake$new_key))
hist_ids<-dplyr::select(hist_just_lake, site, new_key) %>% 
  mutate(row = row_number())
map_data<-left_join(abund_data, hist_ids)

#add up the densities from both gears 
hist_sum_dens<-aggregate(map_data$hist.means, by=list(new_key=map_data$new_key), FUN=sum) %>% 
  rename(sum_dens = x)

#### relative abundance contemporary ####
#extract just lake data (need a single row for each lake) 
cont_just_lake<-test4 %>% 
  distinct(new_key, .keep_all = TRUE) %>% #83 lakes 
  select(-c(gear2, gear, logeffort, IND, site)) %>% 
  mutate(effort = 1) 

dup<-cont_just_lake[rep(c(1:142),2),] # duplicate all rows 

cont_just_lake<-dup %>% 
  group_by(new_key) %>% 
  mutate(gear = if_else(row_number() %% 2 ==1, "seine", "gill")) %>% #add just seine and gill 
  ungroup %>% 
  mutate(gear_index = ifelse(gear == "gill", 2, 3))


# Container for predicted values
#logq1 is FT which is neutral, 2 is gill, 3 is seine, 4 is shock 
predict_cont<- array(NA, c(nsim,length(cont_just_lake[[1]]))) # 2 dimensions - length of data as rows and sims as columns 
dim(predict_cont)

#use lake/site specific alphas and environment, no gear data 
for(i in 1:nsim ){  #loop over sims 
  for(t in 1:length(cont_just_lake[[1]])){ #loop over covariate data (obs)
    predict_cont[i,t] <- exp(coefs[ID[i],paste0('mu.alpha[',cont_just_lake[['group']][t],']')] + coefs[ID[i],paste0('b1[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_secchi']][t]  + coefs[ID[i],paste0('b2[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_lake_area']][t]  + coefs[ID[i],paste0('b3[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_surface_temp_year']][t]  + 
                               coefs[ID[i],paste0('b4[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_max_depth']][t]  + coefs[ID[i],paste0('b5[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_ws_forest']][t] + coefs[ID[i],paste0('b6[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_ws_wetland']][t] +
                               coefs[ID[i],paste0('b7[',cont_just_lake[['gear_index']][t],']')]*cont_just_lake[['z_doy']][t] + 
                               coefs[ID[i], paste0('logq[',cont_just_lake[['gear_index']][t],']')] + cont_just_lake[['effort']][t]) # exp the predictions because we used a log-link in the Poisson
    
  }
}

#container for p values 
#parametrization is by the dispersion parameter, where prob = size/(size+mu).
#p is the success parameter and r is the dispersion parameter.
p_cont <- array(NA, c(nsim,length(cont_just_lake[[1]])) )
exp_abund_cont <- array(NA, c(nsim,length(cont_just_lake[[1]])) )

for(i in 1:nrow(predict_cont)){ # loop over rows (sims)
  for(t in 1:ncol(predict_cont) ){ #loop over columns(obs)
    p_cont[i,t] <- coefs[ID[i],'r']/(coefs[ID[i],'r'] + predict_cont[i,t])
    exp_abund_cont[i,t] <- rnbinom(1, prob=p_cont[i,t], size= coefs[ID[i],'r'])
  }
}

#From this distribution, you can extract the median and the 2.5% and 97.5% percentiles (95% credible interval),
cont.med.abund <- apply(exp_abund_cont, 2, median )
cont.means.abund <- apply(exp_abund_cont, 2, mean )

#prep data for plotting
cont_abund_data=data.frame(cont.med.abund, cont.means.abund) 
cont_abund_data$row <- as.numeric(row.names(cont_abund_data))
cont_abund_data$log_mean<-log(cont_abund_data$cont.means.abund)

cont_just_lake$site <- as.numeric(as.factor(cont_just_lake$new_key))
cont_ids<-dplyr::select(cont_just_lake, site, new_key) %>% 
  mutate(row = row_number())
cont_map_data<-left_join(cont_abund_data, cont_ids) 

#add up the densities from both gears 
cont_sum_dens<-aggregate(cont_map_data$cont.means.abund, by=list(new_key=cont_map_data$new_key), FUN=sum) %>% 
  rename(sum_dens = x)


#*map predicted abundance ####
library(ggmap)
library(stringr)
library(devtools)
library(mapdata)
library(dplyr)

#get the lat/long points - need lat/lon from Humphries table 
points<-read.csv("Data/MI_data/Humphries_table.csv") %>%
  dplyr::select( New_Key, LAT_DD, LONG_DD) %>%
  rename(new_key = New_Key) 

hist_abund_map_dat<-left_join(hist_sum_dens, points)
cont_abund_map_dat<-left_join(cont_sum_dens, points)

summary(log(hist_abund_map_dat$sum_dens))
summary(log(cont_abund_map_dat$sum_dens ))

#get MI basemap 
MI_basemap<-map_data("state") %>%
  subset(region %in% c("michigan")) # select michigan 
p<-ggplot(data = MI_basemap) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + #this fill MI white and outline is black
  coord_fixed(1.3) 
#map of lmb abund historical 
p+ geom_point(data=hist_abund_map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(sum_dens)), alpha = 0.7) + 
  scale_colour_gradientn(colours = c("tan", "green","skyblue", "midnightblue"),
                         limits = c(0,16))  + 
  labs(color="lmb abund", tag = "c") + #changes the labels on the legend 
  theme(legend.position = "bottom", 
        plot.tag = element_text(), 
        plot.tag.position = c(0.1,0.95)) + 
  ylab(NULL) + xlab(NULL)
#log
hist_map<-p+ geom_point(data=hist_abund_map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(log(sum_dens + 1))), size = 3) + 
  scale_colour_gradient(low = "gray", high = "darkblue", limits = c(0,5))  + 
  labs(color="log relative density")  +#changes the labels on the legend 
  theme(legend.position = "bottom") + ylab(NULL) + xlab(NULL)
hist_map
#map of lmb abund contemporary 
p+ geom_point(data=cont_abund_map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(sum_dens)), alpha = 0.7) + 
  scale_colour_gradientn(colours = c("tan", "green","skyblue", "midnightblue"),
                         limits = c(0,16))  + 
  labs(color="lmb abund", tag = "d")  +  #changes the labels on the legend
  theme(legend.position = "bottom", 
        plot.tag = element_text(), 
        plot.tag.position = c(0.1,0.95)) + 
  ylab(NULL) + xlab(NULL)


#log
cont_map<-p+ geom_point(data=cont_abund_map_dat, aes(x = LONG_DD, y = LAT_DD, colour = c(log(sum_dens + 1))), size = 3) + 
  scale_colour_gradient(low = "gray", high = "darkblue")  + 
  labs(color="log relative density") +  #changes the labels on the legend
  theme(legend.position = "bottom") + ylab(NULL) + xlab(NULL)

#* save as 2-panel plot 
map_model2<-cowplot::plot_grid(hist_map, cont_map, labels = c("a", "b"))

ggsave(plot=map_model2, 
       device = "png", 
       filename = "figures/model2_map_abund.png", 
       dpi = 600, height = 8, width = 12, units = "in",
       bg="#ffffff") #sets background to white 

##### other stats ####
mean(output$BUGSoutput$sims.list$b1[,1] > 0) # Posterior probability that Secchi effect is positive in  fyke 1 (look at all_ests for mean and CRI)
mean(output$BUGSoutput$sims.list$b1[,2] > 0) # Posterior probability that Secchi effect is positive in  gill 2(look at all_ests for mean and CRI)
mean(output$BUGSoutput$sims.list$b1[,4] > 0) # Posterior probability that Secchi effect is positive in  shock 4(look at all_ests for mean and CRI)

#save predicted abundance in each lake 
write.csv(hist_abund_map_dat, "Data/output/model2_hist_abund_pred.csv", row.names = FALSE)
write.csv(cont_abund_map_dat, "Data/output/model2_cont_abund_pred.csv", row.names = FALSE)

### lat plot#### 
cont_dat<-cont_just_lake %>% 
  distinct(new_key, .keep_all = TRUE) %>% 
  left_join(cont_abund_map_dat) %>% 
  select(new_key, z_surface_temp_year, sum_dens, LAT_DD, LONG_DD) %>% 
  mutate(type='contemporary')
hist_dat<-hist_just_lake %>% 
  distinct(new_key, .keep_all = TRUE)%>% 
  left_join(hist_abund_map_dat) %>% 
  select(new_key, z_surface_temp_year, sum_dens, LAT_DD, LONG_DD)%>% 
  mutate(type='historical')
plot_dat<-rbind(cont_dat, hist_dat)

ggplot(data=plot_dat, aes(x=z_surface_temp_year, y=log(sum_dens), color=type)) + 
  geom_point()

lat_regression<-ggplot(data=plot_dat, aes(x=LAT_DD, y=log(sum_dens), color=type)) + 
  geom_point() + 
  geom_smooth(method='lm') + 
  scale_color_manual(values=c("darkblue", "lightsalmon")) + 
  xlab("latitude") + ylab("density (log)")

ggsave("figures/regression_plot.png", lat_regression)

