##Credit to m-clark.github Michael Clark : consultung for statistics, computing, and analytics research at UofM
#Bayesian approach - all parameters are random draws from a distribution.
install.packages('brms')
install.packages('rstan')
library(brms)
library(rstan)

corr_res = lme(gpa ~ occasion,
  data = gpa,
  random = ~ 1 | student,
  correlation = corAR1(form = ~ occasion) #specify spatial autocorrelation 
)


# bays package brms 
#bayes_R2 function to get R2

# Note that scaling numeric predictors benefits here just like it does with lme4, and makes specifying the prior easier as well.
brms::brm(gpa ~ occasion, data = gpa)
brms::brm(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy) #left of the bar is the model formula, right of the bar is varying slope

brm(BOOMSHK ~ lake_area_ac + z + (1|g), 
    data=dat, 
    family = lognormal(link = "identity", link_sigma = "log"), #family agrgument to specify distribution of the response 
    prior = prior(normal(0, 1), class = b,
                  cauchy(0, 5), class = sd),
    chains = 2,
    cores = 2, 
    iter = 2000)

#test run 
brm(fish_count ~ (1|gear), 
    data=count_data, 
    family = poisson(link = "log"), #family agrgument to specify distribution of the response 
    chains = 2,
    cores = 2, 
    iter = 2000)

BRM1 <- brm(weight ~ Diet, data = ChickWeight) 
summary(BRM1)

#this command blocks the popup 
options("buildtools.check" = FALSE)