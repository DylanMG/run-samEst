#Exercise 2: fit time varying models to data using samEst

library(samEst)
library(gsl)
library(ggplot)



#Select your data
stock_dat<-read.csv("data/psalmon_sr.csv")

unique(paste(stock_dat$stock, stock_dat$species))

#choose a stock of your preference from the list   
unique(paste(stock_dat$stock, stock_dat$species))

#eg
sr<-subset(stock_dat,stock=="Cowichan"&species=="Chinook")

#or pick a stock at random
sr<-subset(stock_dat,stock.id==round(runif(1,1,242)))

#modify
srdat<-data.frame(by=sr$broodyear,
  S=sr$spawners,
  R=sr$recruits,
  logRS=log(sr$recruits/sr$spawners))


Spred<-seq(0,max(srdat$S)*2,length.out=nrow(srdat))


#Simple model
p <- ricker_TMB(data=srdat)

Rpred_p<-Spred*p$alpha*exp(-p$beta*Spred)

bac <- ricker_stan(data=srdat,iter = 2000)
Rpred_b<-Spred*exp(b$alpha-b$beta*Spred)


# AR1 model
#TMB -MLE
pac<-ricker_TMB(data=srdat, AC=TRUE)
Rpred_pac<-Spred*exp(pac$alpha-pac$beta*Spred)
plot(Spred,Rpred_pac)
#stan
bac <- ricker_stan(data=srdat,iter = 2000,AC=TRUE)
Rpred_bac<-Spred*exp(bac$alpha-bac$beta*Spred)
plot(Spred,Rpred_bac)

#=====================================================================
# rw in alpha
ptva<- ricker_rw_TMB(data=srdat,tv.par="a")
Rpred_ptva<- matrix(0, nrow=length(Spred),ncol=length(ptva$alpha))
for(i in seq_along(ptva$alpha)){
  Rpred_ptva[,i]<-Spred*exp(ptva$alpha[i]-ptva$beta*Spred)
}

matplot(Spred,Rpred_ptva)

btva <- ricker_rw_stan(data=srdat, par="a",iter = 800)
Rpred_btva<- matrix(0, nrow=length(Spred),ncol=length(btva$alpha))
for(i in seq_along(btva$alpha)){
  Rpred_btva[,i]<-Spred*exp(btva$alpha[i]-btva$beta*Spred)
}
matplot(Spred,Rpred_btva)
#stan
#=====================================================================
#rw in log_beta
ptvb <- ricker_rw_TMB(data=srdat,tv.par="b")
Rpred_ptvb<- matrix(0, nrow=length(Spred),ncol=length(ptvb$beta))
for(i in seq_along(ptvb$beta)){
  Rpred_ptvb[,i]<-Spred*exp(ptvb$alpha-ptvb$beta[i]*Spred)
}

matplot(Spred,Rpred_ptvb)
btvb <- ricker_rw_stan(data=df, par="b",iter = 800)
Rpred_btvb<- matrix(0, nrow=length(Spred),ncol=length(btvb$beta))
for(i in seq_along(btvb$beta)){
  Rpred_btvb[,i]<-Spred*exp(btvb$alpha-btvb$beta[i]*Spred)
}

matplot(Spred,Rpred_btvb)
#=====================================================================
#rw in alpha and log beta
ptvab <- ricker_rw_TMB(data=srdat,tv.par="both")
Rpred_ptvab<- matrix(0, nrow=length(Spred),ncol=length(ptvab$beta))
for(i in seq_along(ptvab$beta)){
  Rpred_ptvab[,i]<-Spred*exp(ptvab$alpha[i]-ptvab$beta[i]*Spred)
}

btvab <- ricker_rw_stan(data=df, par="both",iter = 800, mod=rwab_mod) 
Rpred_btvab<- matrix(0, nrow=length(Spred),ncol=length(btvab$beta))
for(i in seq_along(btvab$beta)){
  Rpred_btvab[,i]<-Spred*exp(btvab$alpha[i]-btvab$beta[i]*Spred)
}

#=====================================================================
#regime shift in alpha
phmma <- ricker_hmm_TMB(data=srdat, tv.par='a')
Rpred_phmma<- matrix(0, nrow=length(Spred),ncol=length(phmma$alpha))
for(i in seq_along(phmma$alpha)){
  Rpred_phmma[,i]<-Spred*exp(phmma$alpha[i]-phmma$beta*Spred)
}


bhmma <- ricker_hmm_stan(data=df, par="a",iter = 800)

Rpred_bhmma<- matrix(0, nrow=length(Spred),ncol=length(bhmma$alpha))
for(i in seq_along(bhmma$alpha)){
  Rpred_bhmma[,i]<-Spred*exp(bhmma$alpha[i]-bhmma$beta*Spred)
}

#=====================================================================
#regime shift in beta
phmmb <- ricker_hmm_TMB(data=srdat, tv.par='b')
Rpred_phmmb<- matrix(0, nrow=length(Spred),ncol=length(phmmb$beta))
for(i in seq_along(phmmb$beta)){
  Rpred_phmmb[,i]<-Spred*exp(phmmb$alpha-phmmb$beta[i]*Spred)
}

bhmmb <- ricker_hmm_stan(data=srdat, par="b",iter = 800)
Rpred_bhmmb<- matrix(0, nrow=length(Spred),ncol=length(bhmmb$beta))
for(i in seq_along(bhmmb$beta)){
  Rpred_bhmmb[,i]<-Spred*exp(bhmmb$alpha-bhmmb$beta[i]*Spred)
}
#=====================================================================
#regime shift in alpha and beta
phmmab <- ricker_hmm_TMB(data=srdat, tv.par='both')
Rpred_phmmab<- matrix(0, nrow=length(Spred),ncol=length(phmmab$alpha))
for(i in seq_along(phmmab$alpha)){
  Rpred_phmmab[,i]<-Spred*exp(phmmab$alpha[i]-phmmab$beta[i]*Spred)
}

bhmmab <- ricker_hmm_stan(data=srdat, par="both",iter = 800) 
Rpred_bhmmab<- matrix(0, nrow=length(Spred),ncol=length(bhmmab$alpha))
for(i in seq_along(bhmmab$alpha)){
  Rpred_bhmmab[,i]<-Spred*exp(bhmmab$alpha[i]-bhmmab$beta[i]*Spred)
}
Spred<-seq(0,max(srdat$S)*2,length.out=nrow(srdat))


#model selection

#MLE
AICdf<-data.frame(AIC=c(p$AICc,
                        pac$AICc,
                        ptva$AICc,
                        ptvb$AICc, 
                        ptvab$AICc,
                        phmma$AICc,
                        phmmb$AICc,
                        phmm$AICc),
                  model=C("simple",
                           "autocorrelation",
                           "random walk alpha",
                           "random walk log beta",
                           "random walk alpha and log beta",
                           "regime shift alpha",
                           "regime shift log beta",
                           "regime shift alpha and log beta")                
                  )                          

BICdf<-data.frame(AIC=c(p$BIC,
                        pac$BIC,
                        ptva$BIC,
                        ptvb$BIC, 
                        ptvab$BIC,
                        phmma$BIC,
                        phmmb$BIC,
                        phmm$BIC),
                  model=C("simple",
                           "autocorrelation",
                           "random walk alpha",
                           "random walk log beta",
                           "random walk alpha and log beta",
                           "regime shift alpha",
                           "regime shift log beta",
                           "regime shift alpha and log beta")                
                  )                          


#plot results
data.frame(model=,
  method=,
  by=,
  a=,
  Smax=
  )
