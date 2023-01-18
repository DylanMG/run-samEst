#----------------------------------------------------------------
#Code to test the installation of TMB and Stan  before installing samEst
#January 2023
#----------------------------------------------------------------

#dependencies
install.packages('gsl')



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


#install latest version of stan
remove.packages(c("StanHeaders", "rstan"))
install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

library(rstan)
standata<-list(R_S = s$logRS,
            N=nrow(s),
            S=s$S)

fit1 <- stan(
  file = "code/stan/ricker_linear.stan",  # Stan program
  data = standata,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1,              # number of cores (could use one per chain)
  refresh = 0             # no progress shown
  )

#install samEst
remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')
