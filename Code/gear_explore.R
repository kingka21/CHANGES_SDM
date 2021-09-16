### comparing gears for relative abundance measures #### 
#Paired study - relative (e.g. 10x as much ) catchability estimates 
#four gear types - gillnet, fyke/trap, electroshock, seining 

#summary function 
summaryfun <- function(x)list(N=length(x),Mean=mean(x),Median=median(x),SD=sd(x),Min=min(x),Max=max(x))


#### LMB data #### 
#lmb_format from script has cpue values only 
#count_lmb_format has count and effort 
lmb_dat<-left_join(count_lmb_format, lmb_format)

#select out pairs 
#gill net : electro fishing if > 1 means gill net caught more  
#large range, median under, mean over 
gill_elec<-filter(lmb_dat, effort_GILL > 0 & effort_SHOCK > 0 )
gill_elec$ratio<-gill_elec$GILL/gill_elec$SHOCK
table<-as.data.frame(summaryfun(gill_elec$ratio))

#gill net to fyke/trap #still a large range, but mean and median both <1 meaning trap more than FT
gill_trap<-filter(lmb_dat, effort_GILL > 0 & effort_FT_NET > 0 )
gill_trap$ratio<-gill_trap$GILL/gill_trap$FT_NET
table[2,]<-as.data.frame(summaryfun(gill_trap$ratio))

#gillnet to seine: med 1, mean 2.5 - gill more than seine
gill_seine<-filter(lmb_dat, effort_GILL > 0 & effort_SEINE > 0 )
gill_seine$ratio<-gill_seine$GILL/gill_seine$SEINE
table[3,]<-as.data.frame(summaryfun(gill_seine$ratio))

#electro to fyke : mean and median one above and one below 1 and big range 
elec_fyke<-filter(lmb_dat, effort_SHOCK > 0 & effort_FT_NET > 0 )
elec_fyke$ratio<-elec_fyke$SHOCK/elec_fyke$FT_NET
table[4,]<-as.data.frame(summaryfun(elec_fyke$ratio))

#electro to seine: both mean and med >1 so elect more than seine 
elec_seine<-filter(lmb_dat, effort_SHOCK > 0 & effort_SEINE > 0 )
elec_seine$ratio<-elec_seine$SHOCK/elec_seine$SEINE
table[5,]<-as.data.frame(summaryfun(elec_seine$ratio))

#fyke to seine: mean and med both >1 so fyke more than seine 
ft_seine<-filter(lmb_dat, effort_FT_NET > 0 & effort_SEINE > 0 )
ft_seine$ratio<-ft_seine$FT_NET/ft_seine$SEINE
table[6,]<-as.data.frame(summaryfun(ft_seine$ratio))


#### WAE data #### 
#join cpue values  and count and effort 
wae_dat<-left_join(count_wae_format, wae_format)
# note: no seine data 

#select out pairs 
#gill net : electro fishing if > 1 means gill net caught more  
gill_elec<-filter(wae_dat, effort_GILL > 0 & effort_SHOCK > 0 )
gill_elec$ratio<-gill_elec$GILL/gill_elec$SHOCK
table1<-as.data.frame(summaryfun(gill_elec$ratio))

#gill net to fyke/trap 
gill_trap<-filter(wae_dat, effort_GILL > 0 & effort_FT_NET > 0 )
gill_trap$ratio<-gill_trap$GILL/gill_trap$FT_NET
table1[2,]<-as.data.frame(summaryfun(gill_trap$ratio))

#electro to fyke 
elec_fyke<-filter(wae_dat, effort_SHOCK > 0 & effort_FT_NET > 0 )
elec_fyke$ratio<-elec_fyke$SHOCK/elec_fyke$FT_NET
table1[3,]<-as.data.frame(summaryfun(elec_fyke$ratio))

#### CIS data #### 
#cpue values and count and effort 
cis_dat<-left_join(count_cis_format, cis_format)
#note no seine data 

#gill net : electro fishing if > 1 means gill net caught more  
gill_elec<-filter(cis_dat, effort_GILL > 0 & effort_SHOCK > 0 )
gill_elec$ratio<-gill_elec$GILL/gill_elec$SHOCK
table2<-as.data.frame(summaryfun(gill_elec$ratio))

#gill net to fyke/trap 
gill_trap<-filter(cis_dat, effort_GILL > 0 & effort_FT_NET > 0 )
gill_trap$ratio<-gill_trap$GILL/gill_trap$FT_NET
table2[2,]<-as.data.frame(summaryfun(gill_trap$ratio))

#electro to fyke 
elec_fyke<-filter(cis_dat, effort_SHOCK > 0 & effort_FT_NET > 0 )
elec_fyke$ratio<-elec_fyke$SHOCK/elec_fyke$FT_NET
table2[3,]<-as.data.frame(summaryfun(elec_fyke$ratio))
