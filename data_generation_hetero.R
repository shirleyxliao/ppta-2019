require(MASS)

#necessary side functions
expit = function(x){
  exp(x)/(1+exp(x))
}

logit = function(x){
  log(x/(1-x))
}

#returns indicator function of whether individual is in area of high/low overlap
check.overlap = function(PS,treatment,cutoff=0.1,slice=20){ 
  slice.boundaries = seq(0,1,length.out=slice) #slice PS distribution into 20 quartiles
  indicator = rep(NA,length(PS))
  
  for(i in 2:slice){ #iterate through slices and see whether observations are in an area of high or low overlap
    
    if(i==1){
      lower = 0
      upper = slice.boundaries[1]
    } else{
      lower = slice.boundaries[i-1]
      upper = slice.boundaries[i]
    }
    
    indices = which(PS>lower & PS<upper)
    
    perc.treated = mean(treatment[indices])
    perc.control = 1-perc.treated
    
    if(perc.treated<cutoff | perc.control<cutoff){
      indicator[indices] = 0
    } else{
      indicator[indices] = 1
    }
  }
  
  indicator
}

#psi = -0.2
M = 90 #number of datasets

D = 5 #number of time points for data collection
n = 5000 #sample size at beginning of study
fixedp = 3 #number of fixed covariates
timep_incl = 3 #number of time varying covariates included in the model
timep = 3 #number of real time varying covariates

#coefficients for baseline covariates in PS model
fixedalpha = matrix(rep(0.3,fixedp),ncol=1)
alphaintercept = rep(0,n) #controls prevalence

#coefficients for baseline covariates in outcome model
fixedbeta = matrix(rep(0.3,fixedp),ncol=1)
#coefficents for time-varying covariates in outcome model
timebeta = matrix(c(rep(seq(0.1,0,length.out=D),each=timep)),ncol=1)
#outcome model intercept
betaintercept = rep(0,n)

#ATO
delta = 0.5

#coefficients for simulating time-varying covariate values conditional on past values + past treatments
gamma = c(0.3,0.5)
tau_t = c(0.1,0.2)
tau_x1 = 0.5
tau_x2 = 0.5

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
  overlap = matrix(NA,nrow=n,ncol=D)
  
  for(d in 1:D){
    
    ##re-create matrices for past treatment and past values of covariates
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
    
    #simulate time-varying covariates
    #dependent on past values + past treatments
    thistimeX = matrix(NA,nrow=n,ncol = timep)
    
    for(p in 1:timep){
      
      thistimeX[,p] = past.treatment%*%t(tau_tx) + tau_x1*past.x1[,p] + tau_x2*past.x2[,p] + rnorm(n)
    }
    
  #load current values of covariates into larger matrix + save it to list
    timeX[,((d-1)*(timep)+1):(d*timep)] = thistimeX
    
    if(timep<timep_incl){
      savedata[[m]]$time.X[[d]] = cbind(thistimeX,matrix(rnorm(n*(timep_incl-timep)),nrow=n))
    } else{
      savedata[[m]]$time.X[[d]] = thistimeX[,1:timep_incl] }
    
    #create coefficient matrix for time-varying covariates in PS model
    timealpha = matrix(0,nrow=timep*D,ncol=1) 
    timealpha[((d-1)*timep+1):(d*timep),1] = rep(1,timep)
    if(d>1){
      timealpha[((d-2)*timep+1):((d-1)*timep),1] = rep(0.5,timep)}
    
    #calculate PS from fixed and time varying covariates
    thistruePS = expit(alphaintercept + t(fixedalpha)%*%t(fixedX) + t(timealpha)%*%t(timeX)+t(past.treatment%*%t(gammax))) 
    truePS[,d] = thistruePS
    
    #simulate treatment
    treatment[,d] = rbinom(n,1,thistruePS)
    
    #calculate overlap
    overlap[,d] = check.overlap(thistruePS,treatment[,d])
  }
  
  #calculate fixed outcome - only members of the COP experience a treatment effect
  
  #identify members of the COP
  times.in.overlap = apply(overlap,1,sum)
  cop = rep(0,n)
  cop[times.in.overlap==D] = 1

  error = rnorm(n,0,1)
  
  outcome = betaintercept + delta*cum.treat*cop + t(fixedbeta)%*%t(fixedX) + t(timebeta)%*%t(timeX) + error
  
  #Save data
  savedata[[m]]$outcome = t(outcome)
  savedata[[m]]$treatment = treatment
  savedata[[m]]$truePS = truePS
  savedata[[m]]$overlap = overlap
}


  #save(savedata,file="/Users/shirleyliao/Dropbox/Shirley's Dissertation/Paper2/Code/Data generation/data_overlap_3.R")
##################################

#calculating the ATO for bias in Section 4.3.2

if(FALSE){
cum.treat = apply(treatment,1,sum)
#data.s = data.frame(outcome,cum.treat)
outcome = c(outcome)

ipw = apply((treatment/truePS + (1-treatment)/(1-truePS)),1,prod)
ow = apply((treatment*(1-truePS) + (1-treatment)*truePS),1,prod)

summary(lm(outcome~cum.treat,weights=ipw))

summary(lm(outcome~cum.treat,weights=ow))

##################################
props = table(times.in.overlap)/n

}