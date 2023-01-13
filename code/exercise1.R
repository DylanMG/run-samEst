#Exercise 1: Classic Ricker, sampling error, reference points, and autocorrelation

#Dependencies: please install these
#install.packages('gsl')
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')
library(samEst)
library(gsl)

#Function: this function gives a stochastic draw from a Ricker function with defined parameters
ricker_RS=function(alpha,beta,S,sigma){
  RS=rnorm(1,alpha - beta*S,sigma) #this draws the realized productivity (log RS) with the specified level of error
  return(RS)
}

#E1.1: Classic Ricker####

#Let's simulate a super simplified S-R time-series with these parameters:

#Key parameters: alpha, beta, and sigma
alpha<- 2
beta<- 1/(1000) #often easier to work in Smax (1/beta) as its far more interpretable
sigma<- 0.6 #ranges from ~0.2 to 1.5 for most Pacific salmon stocks
Seq<-alpha/beta

#Simulation parameters:
L=34#length of time-series - here 30 years
RS=numeric(L) #productivity in each year
S=numeric(L);S[1:4]=Seq #spawners in each year, we start with an initial escapement of 600 individuals
R=numeric(L) #recruits in each year
U=runif(L,0.4,0.6) #Annual harvest rate with low annual variance 
#(note U should be within the range of 0 to 1)

#Modify the above parameters to see how this changes the realized S-R data and our inferences
#eg. time-series length, alpha, beta, sigma, U bounds

for(t in 5:L){
  RS[t]=ricker_RS(alpha,beta,sigma,S=S[t-4]) #draw productivity in each year
  R[t]=exp(RS[t])*S[t-4] #transform into recruits by converting logRS to RS (recruits per spawner) times spawners
  S[t]=R[t]*U[t] #escapement left after harvest
}
S=S[-((L-3):L)] #Chop off the time seres so the obs match
R<-R[5:L]
RS<-RS[5:L]
par(mfrow=c(2,1))
plot(R~S,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5)) #Spawner Recruit curve
plot(RS~S,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5)) #Spawner Recruit curve

#Okay, we have our dataset - let's try fitting a model to recover the Ricker parameters

#Simplest route: good old linear regression
m=lm(RS~S)
summary(m)

#More complex routes - samEst in TMB or Stan


#E1.2 - Sampling Error####
#What if we don't have a perfect census?

#Let's define the level of error (coefficient of variation here - eg. sd/mean)
CV=0.1
#Adjust the CV to different levels (in reality we would... hope.. it usually falls in the 0.1 to 0.3 realm)

#Let's add in sampling error to just our estimate of spawners returning to natal streams
S_obs=rnorm(length(S),S,CV*S)

#Similarly this will change our estimates of productivity from the sampled data:
RS_obs=log(R/S_obs)

#Plot out S-R fits
par(mfrow=c(2,1))
plot(R~S_obs,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5)) #Spawner Recruit curve
plot(RS_obs~S_obs,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5)) #Spawner Recruit curve

m=lm(RS~S_obs)
summary(m)

#Error in our census of recruits (real world example: missing components of catch)
R_obs=rnorm(length(R),R,CV*R)
#And what our estimates of productivity would then be from sampled data:
RS_obs=log(R_obs/S_obs)

#Plot out S-R fits
par(mfrow=c(2,1))
plot(R_obs~S_obs,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5)) #Spawner Recruit curve
plot(RS_obs~S_obs,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5)) #Spawner Recruit curve

#Does this change the parameters?
m2=lm(RS_obs~S_obs)
summary(m2)


#What if the sampling error changes over time?
CV1=0.4 #coefficient of variation - period 1
CV2=0.1 #coefficient of variation - period 2
P1=round((L-4)*0.3) #Sampling period 1

#generate observations for each period:
S_obs[1:P1]=rnorm(P1,S[1:P1],CV1*S[1:P1]) #Estimates in the first sampling period
S_obs[c(P1+1):(L-4)]=rnorm((L-4)-P1,S[c(P1+1):(L-4)],CV2*S[c(P1+1):(L-4)]) #Estimates in the second sampling perido
R_obs[1:P1]=rnorm(P1,R[1:P1],CV1*R[1:P1]) #Estimates in the first sampling period
R_obs[c(P1+1):(L-4)]=rnorm((L-4)-P1,R[c(P1+1):(L-4)],CV2*R[c(P1+1):(L-4)]) #Estimates in the second sampling perido

RS_obs=log(R_obs/S_obs)

par(mfrow=c(2,1))
plot(R_obs~S_obs,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5)) #Spawner Recruit curve
plot(RS_obs~S_obs,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5)) #Spawner Recruit curve

m3=lm(RS_obs~S_obs)
summary(m3)

#Stochastic sampling errors do not necessarily create biases, just change precision...

#E1.3 - Reference points####

#Ultimately these parameters from the Ricker curve are important for estimating stock reference points

#Smax - Spawner abundance that maximizes surplus production (AKA production capacity)
#Simply the inverse of the per capita density-dependence parameter (AKA capacity rate parameter)
Smax = 1/beta

#Smsy - Spawner abundance that maximizes sustainable yield
Smsy = (1 - gsl::lambert_W0(exp(1 - alpha)))/beta

#lambert function allows for an explicit solution for calculating Smsy (Scheuerell 2016 PeerJ)

#Umsy - the corresponding harvest rate that would achieve Smsy
Umsy = 1-gsl::lambert_W0(exp(1-alpha))

#Check how estimated parameters relate to the real ones

Smax_obs=1/-m3$coefficients[2]

Smax
Smax_obs

Smsy_obs=(1 - gsl::lambert_W0(exp(1-m3$coefficients[1])))/-m3$coefficients[2]

Smsy
Smsy_obs

Umsy_obs=1-gsl::lambert_W0(exp(1-m3$coefficients[1]))

Umsy
Umsy_obs


#E1.4 - Autocorrelated residuals####
