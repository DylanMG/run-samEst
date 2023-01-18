#Function: this function gives a stochastic draw from a Ricker function with defined parameters
ricker_RS=function(alpha,beta,S,sigma){
  RS=rnorm(1,alpha - beta*S,sigma) #this draws the realized productivity (log RS) with the specified level of error
  return(RS)
}

#E3.1: Regime shifts####

#Let's simulate a super simplified S-R time-series with these parameters:

#Simulation parameters:
L=30 #length of time-series - here 30 years
RS=numeric(L) #productivity in each year
S=numeric(L);S[1]=600 #spawners in each year, we start with an initial escapement of 600 individuals
R=numeric(L) #recruits in each year
U=runif(L,0.4,0.6) #Annual harvest rate with low annual variance 
#(note U should be within the range of 0 to 1)

#Key parameters: alpha, beta, and sigma
#Now let's make distinct sets of parameters - eg. alpha
alpha<- c(0.8,2) #two regimes: very low productivity and high productivity
alpha_t<- numeric(L) #we'll set the durations for each regime through the time-series now
alpha_t[1:round(L*0.33)]=alpha[2] #high productivity for first 1/3rd
alpha_t[c(round(L*0.33)+1):round(L*0.66)]=alpha[1] #low productivity for middle 1/3rd
alpha_t[c(round(L*0.66)+1):L]=alpha[2] #and a return to high productivity
alpha_t
beta<- 1/(1000) #often easier to work in Smax (1/beta) as its far more interpretable
sigma<- 0.6 #ranges from ~0.2 to 1.5 for most Pacific salmon stocks


#Modify the above parameters to see how this changes the realized S-R data and our inferences
#eg. time-series length, alpha, beta, sigma, U bounds

for(t in 1:L){
  RS[t]=ricker_RS(alpha,beta,sigma,S=S[t]) #draw productivity in each year
  R[t]=exp(RS[t])*S[t] #transform into recruits by converting logRS to RS (recruits per spawner) times spawners
  S[t+1]=R[t]*U[t] #escapement left after harvest
}
S=S[-L+1] #Chop off the final spawner estimate to retain just the 40 years of observations

par(mfrow=c(2,1))
plot(R~S,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5)) #Spawner Recruit curve
plot(RS~S,bty='l',pch=21,bg=adjustcolor('black',alpha.f=0.5)) #Spawner Recruit curve

#Let's see if we can capture this using a hidden markov model approach implemented in samEst
df=data.frame(S=S,R=R,logRS=RS,by=seq(1:L))

#This may take a bit initially as it has to compile the stan code... run this line and take a break
m_hmm=samEst::sr_mod(type='hmm',par='a')

fit_hmm=rstan::sampling(m_hmm, data=df,
                 warmup = 200, 
                chains = 6, iter = 700)

hmm_fit=samEst::ricker_hmm_stan(data=df,par='a',k_regime=2)
