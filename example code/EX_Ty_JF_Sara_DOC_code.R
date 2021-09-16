#eample code from project with Sarah and JF on trend in DOC – give the R2 for level 1 and 2 predictors 
#In model change minus part
#Need to track this variable 
#It calculates the posterior residual for every observation and for regional scale, so takes a while 
#Aj is intercept 
#Neg value is possible 

# rm(list=ls())
library(data.table)
library(lme4)
library(jagsUI)
library(MCMCpack)
library(arm)

library(car)
library(ggmcmc)
library(gridExtra)
library(ggthemes)
library(coda)

dat <- readRDS('../Datasets//doc_july18.rds')
dat$ly.med=log(dat$ly.med)
head(dat)
dim(dat)
length(unique(dat$lagoslakeid))

#### Read in Covariate data
covs <- readRDS('../Datasets/DOCCovars_Mar19.rds')
head(covs)
dim(covs)

covs <- data.table(covs)

covs[, lagoslakeid:=as.factor(lagoslakeid)]

# Select covariates for lakes in dat
cov_subset <- covs[covs$lagoslakeid %in% dat$lagoslakeid,]
dim(cov_subset)
length(unique(cov_subset$lagoslakeid))
head(cov_subset)
summary(cov_subset)

# Remove lakes with missing covariate data
cov_subset <- cov_subset[!is.na(urban)]
summary(cov_subset)
cov_subset <- cov_subset[!is.na(maxdepth)]
cov_subset <- cov_subset[!is.na(so4changepct)]
cov_subset <- cov_subset[!is.na(lakemedp)]
summary(cov_subset)
dim(cov_subset)

# Remove lakes from dat that have missing covariate data
dat <- dat[dat$lagoslakeid %in% cov_subset$lagoslakeid,]
length(unique(dat$lagoslakeid))
# Sort both data sets by lakd id
dat <- dat[order(dat$lagoslakeid),]
cov_subset <- cov_subset[order(cov_subset$lagoslakeid),]

# Transform few covariates
cov_subset[, maxdepth:=log(maxdepth)]
cov_subset[, lake_area_ha:=log(lake_area_ha)]
cov_subset[, wetlands:=logit(wetlands)]
cov_subset[, agri:=logit(agri)]
cov_subset[, urban:=logit(urban)]
cov_subset[, forest:=logit(forest)]
cov_subset[, lakemedp:=log(lakemedp)]

# Standardize covariates
changeCols <- colnames(cov_subset)[c(2,4:8,11:16,18,20,22)]
cov_subset[,(changeCols):= lapply(.SD, scale), .SDcols = changeCols]

# Covariate matrix for modeling slopes
x_mat <- as.matrix(cov_subset[,c(2,18,22,4:8,11:16,20)])
colnames(x_mat)
cor(x_mat)
dim(x_mat)


################################################################
########## BUGS CODE ############################################
#################################################################
sink("Model.txt")
cat("
    model {
    for (i in 1:n){
    y[i] ~ dnorm (y.hat[i], tau.y)
   
    y.hat[i] <- alpha[group[i]] + beta[group[i]] * x[i] 
    e.y[i] <- y[i] - y.hat[i]
   
    }
   
    tau.y <- pow(sigma.y, -2)
    sigma.y ~ dunif (0, 10)
   
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] <- BB[j,1]
    beta[j] <- BB[j,2]
   
    BB[j,1:K] ~ dmnorm (BB.hat[j,], Tau.B[,])
    BB.hat[j,1] <- mu.a
    e.a[j] <- alpha[j] - BB.hat[j,1]
    BB.hat[j,2] <- mu.b + b[1] * z[j, 1] + b[2] * z[j, 2] + b[3] * z[j, 3] +  b[4] * z[j, 4] +  b[5] * z[j, 5] +  b[6] * z[j, 6] +
    b[7] * z[j, 7] + b[8] * z[j, 8] + b[9] * z[j, 9] +  b[10] * z[j, 10] +  b[11] * z[j, 11] +  b[12] * z[j, 12]  +
    b[13] * z[j, 13] +  b[14] * z[j, 14]  + b[15] * z[j, 15]
    e.b[j] <- beta[j] - BB.hat[j,2]
   
    
    }
   
    
    
    # Bayesian LASSO -  a Laplace (double exponential) prior
    for(k in 1:ncovs){
    b[k] ~ ddexp(0, lambda1)
    }
   
    # Hyper-prior on lambda
    lambda1 ~ dexp(10)
    mu.a ~ dnorm(0,0.0001)
    mu.b ~ dnorm(0,0.0001)
   
    
    # Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
   
    }
    ",fill=TRUE)
sink()

# Number of sites
J <-length(unique(dat$lagoslakeid))


# Create identity matrix for Wishart dist'n
#!!!!!!!Number of parameters to estimate (K)
K <- 2

W <- diag(K)

# Make sure ID is numeric and goes from 1-number of fisj
dat$ID2 <- as.numeric(as.factor(as.numeric(dat$lagoslakeid)))

# load data
data <- list(y = dat$ly.med , x = dat$year1, group = dat$ID2, n = dim(dat)[1],
             J = J, W = W, K = K, ncovs = dim(x_mat)[2], z = x_mat)


# Initial values
inits <- function (){
  list (BB=array(c(rep(rnorm(1,0,1),J),rep(rnorm(1,0,1),J)), c(J,K)),
        mu.a=rnorm(1,0,1),mu.b=rnorm(1,0,1),
        sigma.y=runif(1,0,10),
        Tau.B=rwish(K+1,diag(K))    )
}

# Parameters monitored
params1 <- c("BB","mu.a","mu.b", "sigma.y","sigma.B","rho.B","b","alpha","beta","e.a","e.b","e.y")


# MCMC settings
ni <- 10000
nt <- 1
nb <- 5000
nc <- 3


out1 <- jags(data, inits, params1, "Model.txt", n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb, parallel = T)

saveRDS(out1, file='DOC_jags.rds')

# Check model convergence
# out.mcmc <- as.mcmc(out1)
# S <- ggs(out.mcmc$samples)
# ggs_traceplot(S, family="mu.a")
# ggs_traceplot(S, family="mu.b")
# ggs_traceplot(S, family="^b")
# ggs_traceplot(S, family="sigma.B")
# ggs_traceplot(S, family="sigma.y")
# ggs_traceplot(S,family="BB\\[1,.\\]")

# Summarize posteriors
print(out1, dig = 3)

y = dat$ly.med
# data level summaries
rsquared.y <- 1 - mean (apply (out1$sims.list$e.y, 1, var))/ var (y)
# summaries for the intercept model
rsquared.a <- 1 - mean (apply (out1$sims.list$e.a, 1, var))/ mean (apply (out1$sims.list$alpha, 1, var))
# summaries for the slope model
rsquared.b <- 1 - mean (apply (out1$sims.list$e.b, 1, var))/ mean (apply (out1$sims.list$beta, 1, var))
print (round (c (rsquared.y, rsquared.a, rsquared.b), 2))


# BugsOut <- out1$summary
# write.csv(BugsOut, "BUGSutputSummary.csv", row.names = T)
# Export MCMC samples
# mcmcOut <- out1$sims.list
# saveRDS(mcmcOut, file="jags_out.rds")
# read back in mcmcOut
# out <- readRDS("jags_out.rds")
# str(out)

###############################
## PLOT COVARIATE EFFECTS #####
# Select random slopes 
mean.beta <- out1$mean$BB[,2]

# Fake data to predict
fake1 <- matrix(NA, nrow=50, ncol=dim(x_mat)[2])
for(j in 1:dim(x_mat)[2]){
  fake1[,j] <- seq(min(x_mat[,j]), max(x_mat[,j]), length=50)
}


# Obtain  CIs and fitted lines
est.lineA <- array(NA, dim=c( out1$mcmc.info$n.samples, dim(fake1)[1], dim(fake1)[2])  ) #container for predicted values
dim(est.lineA)

for(i in 1:out1$mcmc.info$n.samples ){
  for(j in 1:dim(fake1)[1] ){
    for(k in 1:dim(fake1)[2]){
      est.lineA[i,j,k] <- out1$sims.list$mu.b[i] + out1$sims.list$b[i, k] * fake1[j, k]
    }
  }
}


# CIs for fitted values
fit1 <- apply(est.lineA, c(2,3), mean )
upper.CIA <- apply(est.lineA, c(2,3), quantile, 0.975)
lower.CIA <- apply(est.lineA, c(2,3), quantile, 0.025)


## Grab 90% CIs for beta's
CRIs <- matrix(NA, nrow=2, ncol=length(mean.beta))
for(i in 1:length(mean.beta)) {
  CRIs[,i] <- quantile(out1$sims.list$BB[,i,2],probs=c(0.025, 0.975) )
}

# Significant indicator for plotting
beta_sig <- apply(out1$sims.list$b,2,quantile, c(0.025,0.975))
beta_sig_indicator <- numeric()
for(j in 1:dim(x_mat)[2]){
  beta_sig_indicator[j] <- as.numeric(beta_sig[1,j] * beta_sig[2,j] > 0)
}
# color <- rep("black", dim(x_mat)[2])
# color[beta_sig_indicator==1] <- "blue"

# Calculate covariate-specific probability that it is in
# direction of the posterior mean
# lake specific slopes
cov.slopes <- out1$sims.list$b
str(cov.slopes)
# sign of slope
cov.slopeSign <- out1$mean$b > 0    
# calc
cov.slopeProbs <- numeric()
for(i in 1:15){
  if(cov.slopeSign[i] > 0){
    cov.slopeProbs[i] <- mean(cov.slopes[,i] > 0)
  } else {
    cov.slopeProbs[i] <- mean(cov.slopes[,i] < 0)
  }
}
# Trend probs for each lake
cov.slopeProbs

color <- rep("black", dim(x_mat)[2])
color[cov.slopeProbs >= 0.9] <- "blue"

#################
# Create figure
# Panel names
pan.names <- c(expression(paste(SO[4],"change")),"Precipitation change", "Temperature change",expression(paste(NO[3],"Deposition")),
               "Runoff", "MAP", "MAT", "Baseflow", "Wetlands", "Agriculture", "Urban", "Forest", "Lake depth", "Lake area", "Median TP")
# X-axis xoordinate for placing text
x.text.loc <- c(2.1, 1.0, 0.5, 1.7, 1.3, 1.7, 2.3, 1.85, 2.25, 3.9, 6.4, 1.5, 1.3, 2.1, 2.8)


pdf("DOC_cov_plot_mar19.pdf", height = 6, width = 6)
def.par <- par(no.readonly = TRUE)

nf <- layout(matrix(c(1:15),nrow=5,ncol=3,byrow=TRUE),  TRUE)
# layout.show(nf)
par(oma=c(2,3.0,0,1),mai=c(0.10,0.1,0.05,0) )


size.labels = 1
size.text = 1
axissize <- 1
# ylabnums <- seq(min(lower.CIA),max(upper.CIA),length=5)
ylabnums <- seq(min(CRIs[1,]),0.18,length=5)

x.label = 'Standardized covariate'
y.label = expression(paste('Lake-specific DOC trend (',beta[j], ')'))

for(j in 1:dim(x_mat)[2]){
  plot(x_mat[,j], mean.beta ,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n',
       ylim=c(min(CRIs[1,]), 0.18) )
  
  axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  
  if( j ==1 | j==4 | j==7 | j==10 | j==13){
    axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1,at=format(pretty(ylabnums), digits=2),
         labels=format(pretty(ylabnums), digits=2))
  } else {
    axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
  }   
  
  polygon(x=c(fake1[,j],rev(fake1[,j])),y=c(lower.CIA[,j],rev(upper.CIA[,j])),
          col=adjustcolor('gray',0.8),border=NA)
  
  points(x_mat[,j], mean.beta,pch=16,cex=0.5, col='black')
  
  segments(x0=x_mat[,j], x1=x_mat[,j],
           y0=CRIs[1,], y1=CRIs[2,], lwd=0.5, col='black')
  
  
  lines(fake1[,j],fit1[,j], lwd = 2, col="black", lty = 1)
  
  mtext(x.label, line = 0.4, side = 1, cex = size.text, outer=T)
  mtext(y.label, line = 1.4, side = 2, cex = size.text, outer=T)
  
  text(x.text.loc[j],0.15, pan.names[j], col=color[j])
  text(x.text.loc[j],0.11, round(cov.slopeProbs[j],2), col=color[j])
  
  box()
} # end for loop
dev.off()

# END PLOT



# Calculate lake-specific probability of trend being in the
# direction of the posterior mean
# lake specific slopes
slopes <- out1$sims.list$BB[, ,2]
str(slopes)
# sign of slope
slopeSign <- out1$mean$BB[,2] > 0    
# calc
slopeProbs <- numeric()
for(i in 1:J){
  if(slopeSign[i] > 0){
    slopeProbs[i] <- mean(slopes[,i] > 0)
  } else {
    slopeProbs[i] <- mean(slopes[,i] < 0)
  }
}
# Trend probs for each lake
