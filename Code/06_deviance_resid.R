#### Deviance and Residual Deviance #### 


#### compare model 1 and model 2 #### 
mod1_f1<-read.csv("Data/output/model1_deviance_fold1.csv")
mod1_f2<-read.csv("Data/output/model1_deviance_fold2.csv")
mod1_f3<-read.csv("Data/output/model1_deviance_fold3.csv")
mod1_f4<-read.csv("Data/output/model1_deviance_fold4.csv")
mod1_f5<-read.csv("Data/output/model1_deviance_fold5.csv")

mod2_f1<-read.csv("Data/output/model2_deviance_fold1.csv")
mod2_f2<-read.csv("Data/output/model2_deviance_fold2.csv")
mod2_f3<-read.csv("Data/output/model2_deviance_fold3.csv")
mod2_f4<-read.csv("Data/output/model2_deviance_fold4.csv")
mod2_f5<-read.csv("Data/output/model2_deviance_fold5.csv")


(sum(mod1_f1$deviance) + sum(mod1_f2$deviance) + sum(mod1_f3$deviance) + sum(mod1_f4$deviance) +
      sum(mod1_f5$deviance)) / 5

(sum(mod2_f1$deviance) + sum(mod2_f2$deviance) + sum(mod2_f3$deviance) + sum(mod2_f4$deviance) +
    sum(mod2_f5$deviance)) / 5


#### compare model 1 contemporary and historical #### 
mod1_f1<-read.csv("Data/output/model1_cont_res_dev_fold1.csv")
mod1_f2<-read.csv("Data/output/model1_cont_res_dev_fold2.csv")
mod1_f3<-read.csv("Data/output/model1_cont_res_dev_fold3.csv")
mod1_f4<-read.csv("Data/output/model1_cont_res_dev_fold4.csv")
mod1_f5<-read.csv("Data/output/model1_cont_res_dev_fold5.csv")

#put all test outputs into one column 
mod1_hist_f1<-read.csv("Data/output/model1_hist_res_dev_fold1.csv") %>% 
  pivot_longer(cols=c(test1, test2, test3, test4, test5), 
               names_to = "test", 
               values_to = "deviance_res")
hist(mod1_hist_f1$res_dev)
mod1_hist_f2<-read.csv("Data/output/model1_hist_res_dev_fold1.csv")%>% 
  pivot_longer(cols=c(test1, test2, test3, test4, test5), 
               names_to = "test", 
               values_to = "deviance_res")

mod1_hist_f3<-read.csv("Data/output/model1_hist_res_dev_fold1.csv")%>% 
  pivot_longer(cols=c(test1, test2, test3, test4, test5), 
               names_to = "test", 
               values_to = "deviance_res")
mod1_hist_f4<-read.csv("Data/output/model1_hist_res_dev_fold1.csv")%>% 
  pivot_longer(cols=c(test1, test2, test3, test4, test5), 
               names_to = "test", 
               values_to = "deviance_res")
mod1_hist_f5<-read.csv("Data/output/model1_hist_res_dev_fold1.csv")%>% 
  pivot_longer(cols=c(test1, test2, test3, test4, test5), 
               names_to = "test", 
               values_to = "deviance_res")

#### HISTOGRAMS #### 
cont_histogram<-mod1_f1 %>% 
  select(deviance_res) %>% 
  mutate(time = "contemporary")

hist(cont_histogram$deviance_res)

hist_histogram<-mod1_hist_f1 %>% 
  select(deviance_res) %>% 
  mutate(time = "historical")

hist(hist_histogram$deviance_res)

dev_res_histogram<-rbind(cont_histogram, hist_histogram)

dev_res<-ggplot(dev_res_histogram, aes(x=deviance_res, fill=time)) +
  geom_histogram(alpha=0.25, position="identity", aes(y = ..density..), color="black", bins = 10) +
  labs(x="deviance residuals", y = "frequency") + 
  scale_fill_manual(values=c("lightblue", "lightsalmon")) +
  guides(fill=guide_legend(title='')) + 
  theme(legend.position = c(0.9,0.8), legend.text=element_text(size=14))

ggsave(plot=dev_res, 
       device = "png", 
       filename = "figures/dev_res_histogram.png", 
       dpi = 600, height = 8, width = 12, units = "in",
       bg="#ffffff") #sets background to white
