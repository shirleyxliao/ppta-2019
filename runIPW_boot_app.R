#loading data
load("~/simulated_app_data.R")
currentdata = new.data

#loading required packages
require(survey)
require(MCMCpack)
expit = function(x){exp(x)/(1+exp(x))}

#parameter specifications

bootstrap=TRUE
K=1000
D = dim(currentdata$treatment)[2]
n = dim(currentdata$fixed.X)[1]
W = 100

# initializing data structures 
result = list()

#ipw variables
save.weights = matrix(NA,nrow=n,ncol=D)
save.boot.weights = matrix(1,nrow=n,ncol=W)

#overlap variables 
save.overlap = matrix(NA,nrow=n,ncol=D)
save.boot.overlap = matrix(1,nrow=n,ncol=W)

#ppta variables
save.membership = array(NA,dim = c(n,K,D))
ppta.boot.save.membership = array(1,dim=c(n,K,W))
cross.probs = array(NA,dim = c(n,K,D))
ppta.size = matrix(NA,nrow=D,ncol=K)
ppta.treated.size = matrix(NA,nrow=D,ncol=K)
ppta.control.size = matrix(NA,nrow=D,ncol=K)

#stabilized variables
save.stable = matrix(NA,nrow=n,ncol=D)
save.boot.stable = matrix(1,nrow=n,ncol=W)

#specifying input data
treatment = data.frame(currentdata$treatment)
fixedcovar = data.frame(currentdata$fixed.X)
timecovar = currentdata$time.X

#naming input data
names(fixedcovar) = paste("fixed.",names(fixedcovar),sep="")
names(timecovar) = paste("time.",names(timecovar),sep="")
treatment.considered = data.frame(treatment[,1])
names(treatment.considered) = "treat.1"
time.covars.considered = data.frame(timecovar[[1]])

#sampling bootstrap draws
if(bootstrap){
  boot.indicies = matrix(sample(1:n,n*W,replace=TRUE),nrow=n,ncol=W)

  } else{
    W = 1
}

###################################################################
# iterate through time points

for(d in 1:D){
    
    data.matrix = cbind(treatment.considered,fixedcovar,time.covars.considered)
    
    #create datasets with bootstrapped indicies
    if(bootstrap){
      ppta.boot.data = apply(boot.indicies,2,function(x){data.matrix[x,]})
    }
    
    #specify formulas for propensity score model and model for stabilized weights
    if(d==1){
    no.treat = names(data.matrix)[names(data.matrix)!="treat.1"]
    form = paste("treat.1 ~",paste(no.treat,collapse="+"))
    form2 = "treat.1~1"
    } else {
      current.treat = paste("treat.",d,sep="")
      no.treat = names(data.matrix)[names(data.matrix)!=current.treat]
      form = paste(current.treat,"~", paste(no.treat,collapse="+"))
      form2 =  paste(current.treat,"~", paste(no.treat[grep("treat*",no.treat)],collapse="+"))
    }
    
####################################################
    
    # ALL WEIGHTING SCHEMES
    PS.fit =  glm(form,data=data.matrix,family="binomial")
    sum.ps.fit = summary(PS.fit)
    
    ipw.PS.est = predict(PS.fit,type="response") 
    stable.fit = glm(form2,data=data.matrix,family="binomial")
    stable.est = predict(stable.fit,type="response") 
    
    ###calculate ipw/overlap/stable weights
    save.weights[,d] = treatment[,d]/ipw.PS.est + (1-treatment[,d])/(1-ipw.PS.est)
    save.overlap[,d] = treatment[,d]*(1- ipw.PS.est) + (1-treatment[,d])*ipw.PS.est
    save.stable[,d] = treatment[,d]*(stable.est/ipw.PS.est) + (1-treatment[,d])*((1-stable.est)/(1-ipw.PS.est))
    
    ### bootstrapped MLE
    if(bootstrap){
      for (w in 1:W){
        boot.treatment = treatment[boot.indicies[,w],]
        PS.fit =  glm(form,data=ppta.boot.data[[w]],family="binomial")
        sum.ps.fit = summary(PS.fit)
        
        PS.est = predict(PS.fit,type="response") 
        stable.est = predict(glm(form2,data=ppta.boot.data[[w]],family="binomial"),type="response")
        
        save.boot.weights[,w] = save.boot.weights[,w]*(boot.treatment[,d]/PS.est + (1-boot.treatment[,d])/(1-PS.est))
        save.boot.stable[,w] = save.boot.stable[,w]*(boot.treatment[,d]*(stable.est/PS.est) + (1-boot.treatment[,d])*((1-stable.est)/(1-PS.est)))
        save.boot.overlap[,w] = save.boot.overlap[,w]*(boot.treatment[,d]*(1- PS.est) + (1-boot.treatment[,d])*PS.est)
       
      }
    
    }

###############################################
# PPTA
    
    #draws from PS distribution
    posteriors = as.matrix(MCMClogit(form,data=data.matrix,burnin=K/2,mcmc=K*20,thin=20))
    save(posteriors,file=paste("~/final_app_",d,"_posteriors.R")) 
    covars = as.matrix(cbind(rep(1,n),data.matrix[,-1]))
    ppta.PS.est = expit(covars%*%t(posteriors)) #PS estimated from draws of alpha
    
    # boostrapping for PPTA
    if(bootstrap){
        ppta.boot.posteriors = lapply(ppta.boot.data,function(x){as.matrix(MCMClogit(form,data=x,burnin=K/10,mcmc=K*10,thin=10))})
        ppta.boot.covars = lapply(ppta.boot.data,function(x){as.matrix(cbind(rep(1,n),x[,-1]))})
        ppta.boot.PS = mapply(function(x,y){as.matrix(expit(x%*%t(y)))},x = ppta.boot.covars,y=ppta.boot.posteriors,SIMPLIFY = FALSE)
        ppta.boot.PS = simplify2array(ppta.boot.PS)
        }

    #calculate in/out group
      
    for(k in 1:K){
    cross.prob = treatment[,d]*(1-ppta.PS.est[,k]) + (1-treatment[,d])*ppta.PS.est[,k]
    cross.probs[,k,d] = cross.prob    
    
    if(bootstrap){
      ppta.boot.treatment = apply(boot.indicies,2,function(x){treatment[x,d]})
      ppta.boot.cross.prob = apply(matrix(1:W,ncol=W),2,function(x){ppta.boot.treatment[,x]*(1-ppta.boot.PS[,k,x])+(1-ppta.boot.treatment[,x])*ppta.boot.PS[,k,x]})
      ppta.boot.save.membership[,k,] =ppta.boot.save.membership[,k,]*apply(ppta.boot.cross.prob,2,function(x){rbinom(n,1,x)})
    }
    
    }
    
    #updating covariates if not last time point
    if(d<D){
      treatment.considered = data.frame(treatment[,(d+1)],treatment[,d])
      time.covars.considered = data.frame(timecovar[[d+1]],timecovar[[d]])
      names(treatment.considered) = paste("treat.",seq(from=(d+1),to=(d)),sep="")
      names(time.covars.considered) = paste(rep(names(data.frame(timecovar[[d+1]])),2),rep(seq(from=(d+1),to=(d)),each=2),sep=".")
      }
    
    }


save(save.weights,file = "~/save.weights.R")
save(boot.indicies,file="~/boot.indicies.R")
save(save.boot.weights,file="~/save.boot.weights.R")

###############################################
# effect estimation

#creating necessary data structures 
outcome = currentdata$outcome
offset = currentdata$offset
cum.treat = apply(treatment,1,sum,na.rm=TRUE)
treat.out = data.frame(cum.treat,outcome,offset)
treat.out$outcome.rate = outcome/offset
treat.out.list = apply(boot.indicies,2,function(x,data){data[x,]},data=treat.out)

# function for calculating bootstrapped SE

coef.fxn = function(x,weights,data){
  weights = weights[,x]
  dat = data[[x]]
  fit = glm(outcome~cum.treat+offset(log(offset)),data=dat,weights=weights,family="poisson")

  xx = summary(fit)$coef
  
  if(dim(xx)[1]>1){
    return(xx[2,1])
  }else{
    return(NA)
  }
}


#####IPW
print("IPW")

outcome.weights = apply(save.weights,1,prod,na.rm=TRUE)
ipw.weights = outcome.weights  
ipw.fit = glm(outcome~cum.treat+offset(log(offset)),data=treat.out,family="poisson",weights=outcome.weights)

ipw.boots = apply(matrix(1:W,nrow=1),2,coef.fxn,weights=save.boot.weights,data=treat.out.list)
  
ipw.preddelta = summary(ipw.fit)$coef[2,1]
ipw.predse = summary(ipw.fit)$coef[2,2]
ipw.bootse = sd(ipw.boots,na.rm=TRUE)
ipw.lower = quantile(ipw.boots,prob=0.025)
ipw.upper = quantile(ipw.boots,prob=0.975)

##################
#overlap weighted regression
print("overlap")

outcome.weights = apply(save.overlap,1,prod,na.rm=TRUE)
overlap.weights = outcome.weights

overlap.fit = glm(outcome~cum.treat+offset(log(offset)),data=treat.out,family="poisson",weights=overlap.weights)

overlap.boots = apply(matrix(1:W,nrow=1),2,coef.fxn,weights=save.boot.overlap,data=treat.out.list)
  
overlap.preddelta = summary(overlap.fit)$coef[2,1]
overlap.predse = summary(overlap.fit)$coef[2,2]
overlap.bootse = sd(overlap.boots,na.rm=TRUE)
overlap.lower = quantile(overlap.boots,prob=0.025)
overlap.upper = quantile(overlap.boots,prob=0.975)

###############################
#stablized weighting
print("stabilized")

outcome.weights = apply(save.stable,1,prod,na.rm=TRUE)
stabilized.weights = outcome.weights
stable.fit2 = glm(outcome~cum.treat+offset(log(offset)),data=treat.out,family="poisson",weights=outcome.weights)

stable.boots = apply(matrix(1:W,nrow=1),2,coef.fxn,weights=save.boot.stable,data=treat.out.list)
  
stabilized.preddelta = summary(stable.fit2)$coef[2,1]
stabilized.predse = summary(stable.fit2)$coef[2,2]
stabilized.bootse = sd(stable.boots,na.rm=TRUE)
stabilized.lower = quantile(stable.boots,prob=0.025)
stabilized.upper = quantile(stable.boots,prob=0.975)

####################
#ppta
print("PPTA")

mean.boots = NA
ppta.fail = rep(0,K)

ppta.fit = rep(NA,K)
ppta.boot.fit = matrix(NA,nrow = W,ncol = K)
final.cross.probs = apply(cross.probs,c(1,2),prod)
final.membership = apply(final.cross.probs,2,function(x){rbinom(n,1,prob=x)})
ppta.final.size = apply(final.membership,2,sum)

for(k in 1:K){
  is.in = final.membership[,k]
  if (sum(is.in)>5){
    treat.out.in = treat.out[is.in==1,]
   ppta.fit[k] = coef(glm(outcome~cum.treat+offset(log(offset)),data = treat.out.in,family="poisson"))[2]    

    } else {
    ppta.fit[k] = NA
    ppta.fail[k] = 1} 
}

ppta.fxn = function(x,boot.data){
  if(sum(x)>10){
    boot.data.in = boot.data[x==1,]
    return(coef(glm(outcome ~ cum.treat + offset(log(offset)),data = boot.data.in,family="poisson"))[2])

  } else{return(NA)}
}

for(w in 1:W){
  #estimate ATE
  boot.data = treat.out[boot.indicies[,w],]

  boot.is.in = matrix(ppta.boot.save.membership[,,w],nrow=n)
  ppta.boot.fit[w,] = apply(boot.is.in,2,ppta.fxn,boot.data = boot.data)
  
}


ppta.fail.perc = mean(ppta.fail)
ppta.marginal = apply(final.membership,1,mean,na.rm=TRUE)
ppta.once = mean(apply(final.membership,1,max,na.rm=TRUE),na.rm=TRUE)
ppta.preddelta = NA
ppta.predse = NA
ppta.boot.predse = NA
mean.boots = NA

if(any(!is.na(ppta.fit))){
  
  ppta.preddelta = mean(ppta.fit,na.rm=TRUE)
  ppta.predse = sd(apply(ppta.boot.fit,1,mean,na.rm=TRUE),na.rm=TRUE) 
}

if(any(!is.na(ppta.boot.fit))){
  ppta.boot.predse = sd(ppta.boot.fit,na.rm=TRUE)
  mean.boots = apply(ppta.boot.fit,1,mean,na.rm=TRUE) 
ppta.lower = quantile(mean.boots,0.025)
ppta.upper = quantile(mean.boots,0.975)} 

###########################
#save results
result = 
  list(list(ipw.preddelta,ipw.predse,ipw.bootse,ipw.weights),list(ppta.preddelta,ppta.predse,
       ppta.marginal,ppta.once,ppta.size,
       ppta.treated.size,ppta.control.size,ppta.fail.perc),list(overlap.preddelta,
  overlap.predse,overlap.bootse,overlap.weights),list(stabilized.preddelta,stabilized.predse,stabilized.bootse,stabilized.weights))

save(result,file="~/final_app_result.R")


