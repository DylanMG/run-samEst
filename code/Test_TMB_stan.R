#----------------------------------------------------------------
#Code to test the installation of TMB and Stan
#January 2023
#----------------------------------------------------------------


#install TMB
install.packages("TMB")
library(TMB)

#test TMB
compile("code/TMB/Ricker_simple.cpp")
dyn.load(dynlib(here("code/TMB/Ricker_simple")))

#run the TMB code
s<-read.csv("data/example.csv")
srm <- lm(s$logRS~ s$S)
SRdata<-list(obs_logRS=s$logR_S,obs_S=s$S)
parameters<- list(
    alpha=srm$coefficients[1],
    logbeta = log(-srm$coefficients[2]),
    logsigobs=log(.4)
    )
obj_simple <- MakeADFun(SRdata,parameters_simple,DLL="Ricker_simple")


#

