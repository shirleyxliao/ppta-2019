##Paper 3 data generation
##Fixed and time varying covariates, variable treatment

require(MASS)

#necessary side functions
expit = function(x){
  exp(x)/(1+exp(x))
}

logit = function(x){
  log(x/(1-x))
}


M = 90 #number of datasets

D = 5 #number of time points for data collection
n = 5000 #sample size at beginning of study
fixedp = 3 #number of fixed covariates
timep_incl = 3 #number of time varying covariates included in the model
timep = 3 #number of real time varying covariates

fixedalpha = matrix(rep(0.3,fixedp),ncol=1) #PS model coefficients for baseline covariates
alphaintercept = rep(0,n) #controls prevalence

fixedbeta = matrix(rep(0.3,fixedp),ncol=1) #outcome coefficients for baseline covariates
timebeta = matrix(c(rep(seq(0.1,0,length.out=D),each=timep)),ncol=1) #outcome coefficients for time-varying covariates
betaintercept = rep(-1,n) #outcome model intercept

#coefficients for time-varying covariates and past treatments in simulating new time-varying covariates
gamma = c(0.3,0.5) 
tau_t = c(0.1,0.2)
tau_x1 = 0.5
tau_x2 = 0.5

delta = 0.5 

savedata = list()

for(m in 1:M){
  savedata[[m]] = list()
  
  #simulate fixed covariates
  fixedX = matrix(mvrnorm(n,mu=rep(0,fixedp),Sigma=diag(fixedp)), ncol = fixedp)
  
  #save fixed covariates
  savedata[[m]]$fixed.X = fixedX
  savedata[[m]]$time.X = list()
  
  #data structures
  treatment = matrix(NA,nrow=n,ncol=D)
  truePS = matrix(NA,nrow=n,ncol=D)
  timeX = matrix(0,nrow=n,ncol=D*timep)
  
  for(d in 1:D){
    
    ##re-build matrix of past treatments and values for time-varying covariates
    if(d==1){
      gammax = matrix(0)
      tau_tx = matrix(0) 
      past.treatment = matrix(rep(0,n),nrow=n)
      past.x1 = matrix(0,nrow=n,ncol=timep)
      past.x2 = matrix(0,nrow=n,ncol=timep)
    } else if (d==2){
      gammax = matrix(gamma[1])
      tau_tx = matrix(tau_t[1])
      past.treatment = matrix(treatment[,1],nrow=n) #n by 1
      past.x1 = matrix(timeX[,1:timep],ncol=timep)
      past.x2 = matrix(0,nrow=n,ncol=timep)
    } else {
      gammax = matrix(gamma,ncol=2)
      tau_tx = matrix(tau_t,ncol=2)
      past.treatment = matrix(treatment[,(d-2):(d-1)],nrow=n) #n by 2
      past.x1 = matrix(timeX[,((d-2)*(timep)+1):((d-1)*timep)],ncol=timep) #n by timep
      past.x2 = matrix(timeX[,((d-3)*(timep)+1):((d-2)*timep)],ncol=timep) #n by timep
    }
    
    #simulate time-varying covariates depending on past values of those covariates and past treatment values
    #auto-regressive correlation structure, decreasing correlation over time
    thistimeX = matrix(NA,nrow=n,ncol = timep)
    for(p in 1:timep){
      thistimeX[,p] = past.treatment%*%t(tau_tx) + tau_x1*past.x1[,p] + tau_x2*past.x2[,p] + rnorm(n)
    }

    #insert time-varying covariate values into larger matrix of all values + save in list
    timeX[,((d-1)*(timep)+1):(d*timep)] = thistimeX
    
    if(timep<timep_incl){
      savedata[[m]]$time.X[[d]] = cbind(thistimeX,matrix(rnorm(n*(timep_incl-timep)),nrow=n))
    } else{
      savedata[[m]]$time.X[[d]] = thistimeX[,1:timep_incl] }
    
    #coefficients of time-varying covariates in PS model (dependence up to two time points previous)
    timealpha = matrix(0,nrow=timep*D,ncol=1) 
    timealpha[((d-1)*timep+1):(d*timep),1] = rep(1,timep)
    if(d>1){
      timealpha[((d-2)*timep+1):((d-1)*timep),1] = rep(0.5,timep)}
    
    #calculate PS from fixed and time varying covariates
    thistruePS = expit(alphaintercept + t(fixedalpha)%*%t(fixedX) + t(timealpha)%*%t(timeX)+t(past.treatment%*%t(gammax))) 
    truePS[,d] = thistruePS
    
    #simulate treatment
    treatment[,d] = rbinom(n,1,thistruePS)
    
  }
  
  #calculate fixed outcome
  error = rnorm(n)
  
  outcome = betaintercept + delta*apply(treatment,1,sum) + t(fixedbeta)%*%t(fixedX) + t(timebeta)%*%t(timeX) + error
  
  
  #Save data
  savedata[[m]]$outcome = t(outcome)
  savedata[[m]]$treatment = treatment
  savedata[[m]]$truePS = truePS
  
}

#save(savedata,file="/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper2/Code/Data generation/data_5.R")

#########################
