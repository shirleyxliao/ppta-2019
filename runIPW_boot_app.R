
#install.packages("survey",lib="/n/home10/silverfoxflower/apps",repos="http://cran.mtu.edu/")
#require(survey,lib.loc="/n/home10/silverfoxflower/apps")
require(survey)
require(MCMCpack)

#library(parallel)

# Calculate the number of cores
#no_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
#no_cores = 4

# Initiate cluster
#cl <- makeCluster(no_cores,type="FORK")

##Read in data
#load("~/Dropbox/Shirley's Dissertation/Paper2/Code/Data generation/data_highbeta.R")


#runtimevarying = function(savedata){

  bootstrap=TRUE
  #bootstrap = FALSE
  expit = function(x){exp(x)/(1+exp(x))}
  
  K=1000
  M = length(savedata)
#  print(M)
  D = dim(savedata[[1]]$treatment)[2]
  n = dim(savedata[[1]]$fixed.X)[1]
  W = 100
  
result = list()
#result.frame = data.frame(method = NA,delta=NA,se = NA)
  
#ipw variables
#ipw.preddelta = rep(NA,M)
#ipw.predse = rep(NA,M)
#ipw.weights = matrix(NA,nrow = n,ncol = M)

#overlap.preddelta = rep(NA,M)
#overlap.predse = rep(NA,M)
#overlap.weights = matrix(NA,nrow = n,ncol = M)

#cefalu.preddelta = rep(NA,M)
#cefalu.predse = rep(NA,M)
cefalu.size = matrix(NA,nrow=D,ncol=K)
cefalu.treated.size = matrix(NA,nrow=D,ncol=K)
cefalu.control.size = matrix(NA,nrow=D,ncol=K)
#cefalu.fail.perc = matrix(NA,nrow = D,ncol = M)
#cefalu.boot.fail.perc = matrix(NA,nrow = D,ncol = M)

#stabilized.preddelta = rep(NA,M)
#stabilized.predse = rep(NA,M)
#stabilized.weights = matrix(NA,nrow = n,ncol = M)

#parSapply(cl,as.character(1:M),function(m,savedata){

#for(m in 1:M){  
#  m = as.numeric(m)
  
m = 1  
#  print(m)
  currentdata = savedata[[m]]
  
  treatment = data.frame(currentdata$treatment)
  fixedcovar = data.frame(currentdata$fixed.X)
  timecovar = currentdata$time.X
  names(fixedcovar) = paste("fixed.",names(fixedcovar),sep="")
  names(timecovar) = paste("time.",names(timecovar),sep="")
  
  if(bootstrap){
    boot.indicies = matrix(sample(1:n,n*W,replace=TRUE),nrow=n,ncol=W)
    #boot.indicies[,1] = 1:n
  } else{
    W = 1
  }
  
    #cefalu.fail = matrix(0,nrow = D,ncol = K)
  #cefalu.boot.fail = array(NA,dim = c(D,K,W))
  
  ##Calculate D sets of propensity scores
  
  treatment.considered = data.frame(treatment[,1])
  names(treatment.considered) = "treat.1"
  
  time.covars.considered = data.frame(timecovar[[1]])
  
  #ipw variables
  save.weights = matrix(NA,nrow=n,ncol=D)
  save.boot.weights = matrix(1,nrow=n,ncol=W)
  
  #overlap variables 
  save.overlap = matrix(NA,nrow=n,ncol=D)
  #save.bayes.overlap = array(NA,dim = c(n,K,D))
  #save.bayes.overlap = matrix(1,nrow=n,ncol=K)
  save.boot.overlap = matrix(1,nrow=n,ncol=W)
  
  #cefalu variables
  save.membership = array(NA,dim = c(n,K,D))
  #save.membership = matrix(1,nrow=n,ncol=K)
  #cefalu.boot.save.membership = array(NA,dim = c(n,K,D,W))
  cefalu.boot.save.membership = array(1,dim=c(n,K,W))
  
  #stabilized variables
  save.stable = matrix(NA,nrow=n,ncol=D)
  save.boot.stable = matrix(1,nrow=n,ncol=W)
  
  #calculate weights
  
  for(d in 1:D){
    
    data.matrix = cbind(treatment.considered,fixedcovar,time.covars.considered)
    
    if(bootstrap){
      cefalu.boot.data = apply(boot.indicies,2,function(x){data.matrix[x,]})
    }
    
    if(d==1){
    no.treat = names(data.matrix)[names(data.matrix)!="treat.1"]
    form = paste("treat.1 ~",paste(no.treat,collapse="+"))
    form2 = "treat.1~1"
    } else {
      current.treat = paste("treat.",d,sep="")
      no.treat = names(data.matrix)[names(data.matrix)!=current.treat]
      form = paste(current.treat,"~", paste(no.treat,collapse="+"))
      form2 =  paste(current.treat,"~", paste(no.treat[grep("treat*",no.treat)],collapse="+"))
      #form2 = paste(current.treat,"~ 1",sep="")
    }
    

    print(form)
    print(form2)

    PS.fit =  glm(form,data=data.matrix,family="binomial")
    sum.ps.fit = summary(PS.fit)
    
    ipw.PS.est = predict(PS.fit,type="response") 
    stable.fit = glm(form2,data=data.matrix,family="binomial")
    stable.est = predict(stable.fit,type="response") 
    
    
    if(bootstrap){
      for (w in 1:W){
        PS.fit =  glm(form,data=cefalu.boot.data[[w]],family="binomial")
        sum.ps.fit = summary(PS.fit)
        
        PS.est = predict(PS.fit,type="response") 
        stable.est = predict(glm(form2,data=cefalu.boot.data[[w]],family="binomial"),type="response")
        
        save.boot.weights[,w] = save.boot.weights[,w]*(treatment[,d]/PS.est + (1-treatment[,d])/(1-PS.est))
        save.boot.stable[,w] = save.boot.stable[,w]*(treatment[,d]*(stable.est/PS.est) + (1-treatment[,d])*((1-stable.est)/(1-PS.est)))
        save.boot.overlap[,w] = save.boot.overlap[,w]*(treatment[,d]*(1- PS.est) + (1-treatment[,d])*PS.est)
       
      }
    
    }
    
    posteriors = as.matrix(MCMClogit(form,data=data.matrix,burnin=K,mcmc=K*50,thin=50))
    save(posteriors,file=paste("~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_",d,"_posteriors.R")) 
      
    if(bootstrap){
        cefalu.boot.posteriors = lapply(cefalu.boot.data,function(x){as.matrix(MCMClogit(form,data=x,burnin=K/10,mcmc=K*10,thin=10))})
        cefalu.boot.covars = lapply(cefalu.boot.data,function(x){as.matrix(cbind(rep(1,n),x[,-1]))})
        cefalu.boot.PS = mapply(function(x,y){as.matrix(expit(x%*%t(y)))},x = cefalu.boot.covars,y=cefalu.boot.posteriors,SIMPLIFY = FALSE)
        cefalu.boot.PS = simplify2array(cefalu.boot.PS)
        }
      
      covars = as.matrix(cbind(rep(1,n),data.matrix[,-1]))
      cefalu.PS.est = expit(covars%*%t(posteriors)) #PS estimated from draws of alpha
    
    #calculate ipw/overlap/stable weights
    save.weights[,d] = treatment[,d]/ipw.PS.est + (1-treatment[,d])/(1-ipw.PS.est)
    save.overlap[,d] = treatment[,d]*(1- ipw.PS.est) + (1-treatment[,d])*ipw.PS.est
    save.stable[,d] = treatment[,d]*(stable.est/ipw.PS.est) + (1-treatment[,d])*((1-stable.est)/(1-ipw.PS.est))
    
    for(k in 1:K){
    #calculate cefalu membership
    cross.prob = treatment[,d]*(1-cefalu.PS.est[,k]) + (1-treatment[,d])*cefalu.PS.est[,k]
    save.membership[,k,d] = rbinom(n,1,cross.prob) 
    #save.bayes.overlap[,k] = save.bayes.overlap[,k]*cross.prob
    
    cefalu.size[d,k] = sum(save.membership[,k,d],na.rm=TRUE)
    cefalu.treated.size[d,k] = sum(save.membership[,k,d]*treatment[,d],na.rm=TRUE)
    cefalu.control.size[d,k] =sum(save.membership[,k,d]*(1-treatment[,d]),na.rm=TRUE)
    
    if(cefalu.treated.size[d,k]==0 | cefalu.control.size[d,k] == 0){
    cefalu.fail[d,k] = 1 }
    
    
    if(bootstrap){
      cefalu.boot.treatment = apply(boot.indicies,2,function(x){treatment[x,d]})
      cefalu.boot.cross.prob = apply(matrix(1:W,ncol=W),2,function(x){cefalu.boot.treatment[,x]*(1-cefalu.boot.PS[,k,x])+(1-cefalu.boot.treatment[,x])*cefalu.boot.PS[,k,x]})
      #cefalu.boot.cross.prob = mapply(function(x,y){x*(1-y)+(1-x)*y},x=cefalu.boot.treatment,y=cefalu.boot.PS[,k,],SIMPLIFY = FALSE)
      cefalu.boot.save.membership[,k,] =cefalu.boot.save.membership[,k,]*apply(cefalu.boot.cross.prob,2,function(x){rbinom(n,1,x)})
    }
    
    }
    
    #updating covariates if not last time point
    if(d<D){
      treatment.considered = data.frame(cbind(treatment[,(d+1)],treatment.considered))
      
      time.covars.considered = data.frame(cbind(data.frame(timecovar[[d+1]]),time.covars.considered))
      
      p_time = dim(data.frame(timecovar[[d+1]]))[2]
      if(d>=2){ #remove old treatments
        treatment.considered = data.frame(treatment.considered[,1:2])
        names(treatment.considered) = paste("treat.",seq(from=(d+1),to=(d)),sep="")
        #names(treatment.considered) = paste("treat.",d+1,sep="")
        time.covars.considered = time.covars.considered[,-((p_time*(4-1)+1):(p_time*4))]
      } else{
        names(treatment.considered) = paste("treat.",seq(from=(d+1),to=1),sep="")
      }
      }
    
    }



  
  outcome = currentdata$outcome
  cum.treat = apply(treatment,1,sum,na.rm=TRUE)

    #####ip weighted regression
  outcome.weights = apply(save.weights,1,prod,na.rm=TRUE)
  ipw.weights = outcome.weights
save(ipw.weights,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_ipw_weights.R")
save(save.boot.weights,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_ipw_boot_weights.R")

  dataset = data.frame(outcome,cum.treat,outcome.weights)
  des = svydesign(~0,data=dataset,weights=outcome.weights)
  #ipw.fit = svyglm(outcome~as.factor(cum.treat),design=des)
  
  ipw.fit = svyglm(outcome~cum.treat,design=des)
  #ipw.preddelta[m,] = summary(ipw.fit)$coef[2:(D+1),1]
  #ipw.predse[m,] = summary(ipw.fit)$coef[2:(D+1),2]

  coef.fxn = function(x){
  des = svydesign(~0,data=data.frame(outcome,cum.treat,x),weights=x)
  xx = summary(svyglm(outcome~cum.treat,design=des))$coef
  if(dim(xx)[1]>1){
    return(xx[2,1])
  }else{
    return(NA)
  }
  }
  
  ipw.boots = apply(save.boot.weights,2,coef.fxn)
  
  save(ipw.boots,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_ipw_boots.R")

  
  ipw.preddelta = summary(ipw.fit)$coef[2,1]
  ipw.predse = summary(ipw.fit)$coef[2,2]
  ipw.bootse = sd(ipw.boots,na.rm=TRUE)
  
  ##################
  #overlap weighted regression
save(save.overlap,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_overlap_weights.R")
save(save.boot.overlap,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_overlap_boot_weights.R")

  outcome.weights = apply(save.overlap,1,prod,na.rm=TRUE)
  overlap.weights = outcome.weights
#save(overlap.weights,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_overlap_weights.R")

  dataset = data.frame(outcome,cum.treat,outcome.weights)
  des = svydesign(~0,data=dataset,weights=outcome.weights)
  #overlap.fit = svyglm(outcome~as.factor(cum.treat),design=des)
  
  #overlap.preddelta[m,] = summary(overlap.fit)$coef[2:(D+1),1]
  #overlap.predse[m,] = summary(overlap.fit)$coef[2:(D+1),2]
  
  overlap.fit = svyglm(outcome~cum.treat,design=des)

overlap.boots = apply(save.boot.overlap,2,coef.fxn)

save(overlap.boots,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_overlap_boots.R")

  
  overlap.preddelta = summary(overlap.fit)$coef[2,1]
  overlap.predse = summary(overlap.fit)$coef[2,2]
  overlap.bootse = sd(overlap.boots,na.rm=TRUE)
  
  ###############################
  #stablized weighting
  
  outcome.weights = apply(save.stable,1,prod,na.rm=TRUE)
  stabilized.weights = outcome.weights
save(stabilized.weights,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_stable_weights.R")
save(save.boot.stable,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_stable_boot_weights.R")

  dataset = data.frame(outcome,cum.treat,outcome.weights)
  des = svydesign(~0,data=dataset,weights=outcome.weights)
  #ipw.fit = svyglm(outcome~as.factor(cum.treat),design=des)
  
  #ipw.preddelta[m,] = summary(ipw.fit)$coef[2:(D+1),1]
  #ipw.predse[m,] = summary(ipw.fit)$coef[2:(D+1),2]
  
  stable.fit2 = svyglm(outcome~cum.treat,design=des)

stable.boots = apply(save.boot.stable,2,coef.fxn)

save(stable.boots,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_stable_boots.R")

  
  stabilized.preddelta = summary(stable.fit2)$coef[2,1]
  stabilized.predse = summary(stable.fit2)$coef[2,2]
stabilized.bootse = sd(stable.boots,na.rm=TRUE)
  
  
  ####################
  #cefalu
#save(cefalu.treated.size,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_cefalu_treated_size.R")
save(cefalu.size,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_cefalu_size.R")
#  cefalu.deltas = rep(NA,K)
#  cefalu.ses = rep(NA,K)
  #bayes.outcome.weights = apply(save.bayes.overlap,c(1,2),prod) #n by k matrix of overlap probs

  once.matrix = matrix(NA,nrow=n,ncol=K)
 # cefalu.strat.var = rep(NA,K)
  cefalu.fail = rep(0,K)
  cefalu.final.size = rep(NA,K)
  cefalu.fit = rep(NA,K)
  cefalu.boot.fit = matrix(NA,nrow = W,ncol = K)
  save(save.membership,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_membership.R")
  
  for(k in 1:K){
  #print(k)
  #estimate ATE
  k.save.membership = save.membership[,k,]
  
  is.in = apply(k.save.membership,1,prod,na.rm=TRUE)
  once.matrix[,k] = is.in
  cefalu.final.size[k] = sum(is.in)
  
  if (sum(is.in)>10){
    
  outcome.in = outcome[is.in==1]
  cum.treat.in = cum.treat[is.in==1]
  
  cefalu.fit[k] = coef(lm(outcome.in ~ cum.treat.in))[2] } else {
    cefalu.fit[k] = NA
    cefalu.fail[k] = 1}
  
  #estimate SE
  if(bootstrap){
  #k.boot.save.membership = cefalu.boot.save.membership[,k,,] 
  #boot.is.in = apply(k.boot.save.membership,c(1,3),prod,na.rm=TRUE) #n by W matrix of who is in
  boot.is.in = cefalu.boot.save.membership[,k,]
  cefalu.boot.fit[,k] = apply(boot.is.in,2,function(x){coef(lm(outcome[x==1] ~ cum.treat[x==1]))[2]})
  }
  }
  
  save(cefalu.final.size,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_final.size.R")
  
  cefalu.marginal = apply(once.matrix,1,mean,na.rm=TRUE)
  cefalu.once = mean(apply(once.matrix,1,max,na.rm=TRUE),na.rm=TRUE)
  cefalu.preddelta = mean(cefalu.fit,na.rm=TRUE)
  cefalu.predse = sd(apply(cefalu.boot.fit,1,mean,na.rm=TRUE),na.rm=TRUE)
  cefalu.fail.perc = mean(cefalu.fail)

#result[[m]] =
 #return(
  
 result = 
  list(list(ipw.preddelta,ipw.predse,ipw.bootse,ipw.weights),list(cefalu.preddelta,cefalu.predse,
       cefalu.marginal,cefalu.once,cefalu.size,
       cefalu.treated.size,cefalu.control.size,cefalu.fail.perc),list(overlap.preddelta,
  overlap.predse,overlap.bootse,overlap.weights),list(stabilized.preddelta,stabilized.predse,stabilized.bootse,stabilized.weights))
#}
#},savedata=savedata)
save(result,file="~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_result")
#return(result)
#}

