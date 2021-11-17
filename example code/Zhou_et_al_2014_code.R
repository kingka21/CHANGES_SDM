###########################
#9/2/2014, flathead, for paper
# 3 gears
# data: structure columns: catch in year.month.grid by swept.area and gear type (specified as 1,2,3)
#    catch[]	ymg[]	swpA[]	gear[]

model { 
  # for Negative binomial distribution, activate p,r,mu
  #p ~ dbeta(1,1)
  #r ~ dlnorm(0,0.1)
  #mu<- r*(1-p)/p
  mu ~ dlnorm(1, 0.1)	
  
  for (k in 1:3) { 
    q[k] ~ dbeta(1,1)   
    } 
  
  for (i in 1: 218) {
    #	Den[i] ~ dnegbin( p, r)
    Den[i] ~ dpois(mu)   
  } 
  
  for (i in 1:699) {
    #dummy[i]<- catch_wt[i]
    den2[i] ~ dpois(Den[ymg[i]])  #ymg: year.month.grid index
    N[i] <- den2[i] * swpA[i]     #swpA: swept.area = effort 
    catch[i] ~ dbin(q[gear[i]], N[i] )
    #	pred.c[i] ~dbin(q[gear[i]], N[i])
    est.c[i]<- q[gear[i]] * N[i]
    
    #	SPE[i] <- pow( pred.c[i] - catch[i], 2)
    #	X2[i] <-       pow( catch[i] -  est.c[i], 2) / pow(sd(catch[]), 2)
    #	X2pred[i] <- pow( pred.c[i] - est.c[i], 2) / pow(sd(pred.c[]), 2)
    #  ppp[i] <- step(X2pred[i] - X2[i])
    
    #	SPE.lg[i]<- pow( log( pred.c[i]+0.01) - log(catch[i]+0.01) , 2)
    #	X2.lg[i]   <- pow( log( catch[i] +0.01) - log( est.c[i] +0.01), 2) #/ pow( sd(est.c[]), 2)
    #	X2.lg.pred[i] <- pow( log(pred.c[i] +0.01) - log( est.c[i] +0.01), 2) #/ pow(sd(est.c[]) ,2)
    #ppp.lg[i] <- step(X2.lg.pred[i] - X2.lg[i])   
  }
  
  Den.mean<- mean(Den[])
  
  #MSPE[1] <- mean(SPE[1:176])
  #MSPE[2]<- mean(SPE[177:627])
  #MSPE[3]<- mean(SPE[628:699])
  
  #X2.sum[1]<- sum(X2[1:176])
  #X2.sum[2]<- sum(X2[177:627])
  #X2.sum[3]<- sum(X2[628:699])
  
  #X2pred.sum[1]<- sum( X2pred[1:176]  )
  #X2pred.sum[2]<- sum( X2pred[177:627]  )
  #X2pred.sum[3]<- sum( X2pred[628:699]  )
  
  #pppval[1] <- mean(ppp[1:176])
  #pppval[2] <- mean(ppp[177:627])
  #pppval[3] <- mean(ppp[628:699])
  
  #sdPredC[1]<-sd(pred.c[1:176])
  #sdPredC[2]<-sd(pred.c[177:627])
  #sdPredC[3]<-sd(pred.c[628:699])
  
  #MSPE.lg [1]<-  mean(SPE.lg[1:176]) 
  #MSPE.lg [2]<-  mean(SPE.lg[177:627]) 
  #MSPE.lg [3]<-  mean(SPE.lg[628:699]) 
  
  #minMaxCpred[1]<- minmax(pred.c[1:176]) 
  #minMaxCpred[2]<- minmax(pred.c[177:627]) 
  #minMaxCpred[3]<- minmax(pred.c[628:699]) 
  
  #X2.lg.sum<- sum(X2.lg[])
  #X2.lg.pred.sum <- sum(X2.lg.pred[])   
  #pppval.lg <- mean(ppp.lg[]) 
}


        


##### general model #### 
#Appendix 5: WinBUGS code for NB-PS model

#negative binomial and poisson distribution 
#Since the negative binomial distribution has one more parameter than the Poisson, 
#the second parameter can be used to adjust the variance independently of the mean

#NB.PS_Model
{ 
  # priors and constants
  nc <- 100	#number of cells
  ns <- 4		# number of samples per gear
  ng <- 2		# number of gears
  a <- 0.1		# fraction of cell fished
  p ~ dbeta(1, 1) 
  r ~ dlnorm(0, 0.01)  
  mu <- r * (1 - p) / p
  
  for (k in 1 : ng) { 
    Q[k] ~ dbeta(1,1)   
  }
  
  for (i in 1 : nc) { 				# for cells
    N1[i] ~ dnegbin(p, r)   		# between cell abundance/density
    for (j in 1 : ns) { 			# for samples
      N2[i, j] ~ dpois(N1[i])     	# within cell abundance
      n2[i, j] <- N2[i, j] * a	# local available abundance for catch
      for (k in 1 : ng) {
        C[i, j, k] ~ dbin(Q[k] , n2[i, j] )  
      }
    }	
  }
}

