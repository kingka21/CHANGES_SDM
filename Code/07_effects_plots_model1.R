#### effects of the predictors plots model 1####

#*load library
library(dotwhisker)

####read in model output ####
output<-readRDS("Data/output/output_model1_secchi_fold3_lakes.rds")
jagsfit.mcmc <- as.mcmc(output)
model1tranformed <- ggs(jagsfit.mcmc)

#betas
bEst <- matrix(NA, nrow=7,ncol=3)
for(i in 1:7){ #parameters
  bEst[i,1] <- mean(output$BUGSoutput$sims.list$b[,i])
  bEst[i,2:3] <- quantile(output$BUGSoutput$sims.list$b[,i],c(0.05,0.95))
}
bEst<-as.data.frame(bEst)
bEst$variable<-c("Secchi depth",  "lake area", "surface temp", "max depth",  
                 "ws forest", "ws wetland", "day-of-year")

colnames(bEst)<- c("estimate", "conf.low", "conf.high", "term")
bEst$Parameter <-c("b[1]", "b[2]", "b[3]", "b[4]","b[5]","b[6]", "b[7]")


#add colors for sig different than 0
bEst$color <- as.numeric(bEst[,2] * bEst[,3] > 0 )

# plot using the dotwhisker package 
lake_plot<-dwplot(bEst, style = "dotwhisker", 
                  dot_args = list(aes(colour = factor(color))),
                  whisker_args = list(aes(colour = factor(color))), 
                  vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)) + # plot line at zero _behind_ coefs
  scale_colour_manual(values= c("black", "blue")) +    
  theme_bw() + 
  xlab("Estimated effect") + ylab("Variable") +
  ggtitle("lakes") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position = "none")
lake_plot

#PLOT WITH DISTRIBUTIONS 
model_1_covariates<-model1tranformed %>%
  filter(grepl("b",Parameter))%>%
  left_join(bEst) %>% #get parameter names and colors 
  ggplot(aes(x = value, y = term, fill= factor(color))) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 2) + 
  scale_fill_manual(values= c("gray80", "skyblue")) +
  ggdist::stat_halfeye(.width = c(.05, .95)) + 
  theme_bw() + 
  theme(legend.position = "none") +
  ylab("predictor")

ggsave(plot=model_1_covariates, 
       device = "png", 
       filename = "figures/model_1_covariates.png", 
       dpi = 600, height = 8, width = 5, units = "in")

#### other stats ####
mean(output$BUGSoutput$sims.list$b[,7] > 0) # Posterior probability that DOY effect is positive 
quantile(output$BUGSoutput$sims.list$b[,7],c(0.05,0.95))
mean(output$BUGSoutput$sims.list$b[,7])
