#### Deviance Residuals  #### 


#### compare model 2 contemporary and historical #### 
mod2_f1<-read.csv("Data/output/model2_cont_res_dev_fold1.csv")
mod2_f2<-read.csv("Data/output/model2_cont_res_dev_fold2.csv")
mod2_f3<-read.csv("Data/output/model2_cont_res_dev_fold3.csv")
mod2_f4<-read.csv("Data/output/model2_cont_res_dev_fold4.csv")
mod2_f5<-read.csv("Data/output/model2_cont_res_dev_fold5.csv")

#put all test outputs into one column 
mod2_hist_f1<-read.csv("Data/output/model2_hist_res_dev_fold1.csv") %>% 
  pivot_longer(cols=c(test1, test2, test3, test4, test5), 
               names_to = "test", 
               values_to = "deviance_res")
hist(mod2_hist_f1$deviance_res)
mod2_hist_f2<-read.csv("Data/output/model2_hist_res_dev_fold1.csv")%>% 
  pivot_longer(cols=c(test1, test2, test3, test4, test5), 
               names_to = "test", 
               values_to = "deviance_res")

mod2_hist_f3<-read.csv("Data/output/model2_hist_res_dev_fold1.csv")%>% 
  pivot_longer(cols=c(test1, test2, test3, test4, test5), 
               names_to = "test", 
               values_to = "deviance_res")
mod2_hist_f4<-read.csv("Data/output/model2_hist_res_dev_fold1.csv")%>% 
  pivot_longer(cols=c(test1, test2, test3, test4, test5), 
               names_to = "test", 
               values_to = "deviance_res")
mod2_hist_f5<-read.csv("Data/output/model2_hist_res_dev_fold1.csv")%>% 
  pivot_longer(cols=c(test1, test2, test3, test4, test5), 
               names_to = "test", 
               values_to = "deviance_res")

#### HISTOGRAMS #### 
cont_histogram<-mod2_f1 %>% 
  select(deviance_res) %>% 
  mutate(time = "contemporary")

hist(cont_histogram$deviance_res)

hist_histogram<-mod2_hist_f1 %>% 
  select(deviance_res) %>% 
  mutate(time = "historical")

hist(hist_histogram$deviance_res)

dev_res_histogram<-rbind(cont_histogram, hist_histogram)

dev_res<-ggplot(dev_res_histogram, aes(x=deviance_res, fill=time)) +
  geom_histogram(alpha=0.5, position="identity", aes(y = ..density..), color="black", bins = 10) +
  labs(x="deviance residuals", y = "frequency") + 
  scale_fill_manual(values=c("lightblue", "lightsalmon")) +
  guides(fill=guide_legend(title='')) + 
  theme(legend.position = c(0.9,0.8), legend.text=element_text(size=14), axis.title = element_text(size=16), axis.text = element_text(size=14)) +
  theme_bw()
dev_res

ggsave(plot=dev_res, 
       device = "png", 
       filename = "figures/fig6_dev_res_histogram.png", 
       dpi = 600, height = 8, width = 12, units = "in",
       bg="#ffffff") #sets background to white

#### Kolmogorov-Smirnov Test #### 
#test the two distributions from each fold
#do the two datasets come from the same distribution 
ks.test(mod2_f1, mod2_hist_f1$deviance_res)
ks.test(mod2_f2, mod2_hist_f2$deviance_res)
ks.test(mod2_f3, mod2_hist_f3$deviance_res)
ks.test(mod2_f4, mod2_hist_f4$deviance_res)
ks.test(mod2_f5, mod2_hist_f5$deviance_res)
(0.1315+ 0.05412+ 0.06793+ 0.05384 + 0.05221) /5
