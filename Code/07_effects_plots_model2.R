#### model effects of the predictor plots model2####

#*load library
library(dotwhisker)

####read in model output ####
output<-readRDS("Data/output/output_model2_secchi_fold5_lakes.rds")

jagsfit.mcmc <- as.mcmc(output)
model1tranformed <- ggs(jagsfit.mcmc) 

#### dot and whisket plots ####

#*pull out betas ####
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

 

#* plot using the dotwhisker package ####
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
