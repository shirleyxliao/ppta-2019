
#install.packages("survey",lib="/n/home10/silverfoxflower/apps",repos="http://cran.mtu.edu/")
#require(survey,lib.loc="/n/home10/silverfoxflower/apps")
require(survey)
require(MCMCpack)

library(parallel)

# Calculate the number of cores
no_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
#no_cores = 4

# Initiate cluster
cl <- makeCluster(no_cores,type="FORK",outfile="")

##Read in data
#load("~/Dropbox/Shirley's Dissertation/Paper2/Code/Data generation/data_highbeta.R")

print("SOURCING ...")

runtimevarying = function(savedata,data.type){
  print("RUNNING FUNCTION ...")
  start.time = Sys.time()

  print("outer shell")
  print(data.type)
  
  bootstrap=TRUE
  expit = function(x){exp(x)/(1+exp(x))}
  

  #M = 200
  K=1500
  M = length(savedata)

  D = dim(savedata[[1]]$treatment)[2]
  n = dim(savedata[[1]]$fixed.X)[1]
  W = 200
  
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

#cefalu.fail.perc = matrix(NA,nrow = D,ncol = M)
#cefalu.boot.fail.perc = matrix(NA,nrow = D,ncol = M)

#stabilized.preddelta = rep(NA,M)
#stabilized.predse = rep(NA,M)
#stabilized.weights = matrix(NA,nrow = n,ncol = M)

parSapply(cl,as.character(1:M),function(m,savedata,data.type){
  print("inner shell")
  print(data.type)
#for(m in 1:M){  
  m = as.numeric(m)
  print("Iteration")
  print(m)
  print(Sys.time() - start.time)
  m.time = Sys.time()
  currentdata = savedata[[m]]
  
  treatment = data.frame(currentdata$treatment)
  fixedcovar = data.frame(currentdata$fixed.X)
  timecovar = currentdata$time.X
  #print("fixedcovar")
  names(fixedcovar) = paste("fixed.",names(fixedcovar),sep="")
  #print("timecovar")
  #names(timecovar) = paste("time.",names(timecovar),sep="")
  
  if(bootstrap){
    boot.indicies = matrix(sample(1:n,n*W,replace=TRUE),nrow=n,ncol=W)
    #boot.indicies[,1] = 1:n
  } else{
    W = 1
  }
  
    #cefalu.fail = matrix(0,nrow = D,ncol = K)
  #cefalu.boot.fail = array(NA,dim = c(D,K,W))
  
  ##Calculate D sets of propensity scores
  
  #print("treatment.considered")
  treatment.considered = data.frame(treatment[,1])
  names(treatment.considered) = "treat.1"
  
  time.covars.considered = data.frame(timecovar[[1]])
  
  #ipw variables
  save.weights = matrix(NA,nrow=n,ncol=D)
  save.boot.weights = matrix(1,nrow=n,ncol=W)
  
  #overlap variables 
  save.overlap = matrix(NA,nrow=n,ncol=D)
  #save.bayes.overlap = array(NA,dim = c(n,K,D))
  save.boot.overlap = matrix(1,nrow=n,ncol=W)
  
  #cefalu variables
  save.membership = array(NA,dim = c(n,K,D))
  #cefalu.boot.save.membership = array(0,dim=c(n,K,W))
  cefalu.boot.save.membership = list()
  
  for(w in 1:W){
    cefalu.boot.save.membership[[w]] = list()
    for(k in 1:K){
      cefalu.boot.save.membership[[w]][[k]] = 1:n
    }
  }
  #stabilized variables
  save.stable = matrix(NA,nrow=n,ncol=D)
  save.boot.stable = matrix(1,nrow=n,ncol=W)

  
  
  #cefalu variables
  #save.membership = array(NA,dim = c(n,K,D))
  #cefalu.boot.save.membership = array(NA,dim = c(n,K,D,W))
  cefalu.size = matrix(NA,nrow=D,ncol=K)
  cefalu.treated.size = matrix(NA,nrow=D,ncol=K)
  cefalu.control.size = matrix(NA,nrow=D,ncol=K)
  
  #stabilized variables
  #save.stable = matrix(NA,nrow=n,ncol=D)
  
  #calculate weights
  
for(d in 1:D){
    print("Time from start to time point ")
    print(d)
    print(Sys.time() - m.time)
    
    data.matrix = cbind(treatment.considered,fixedcovar,time.covars.considered)
    cefalu.boot.data = apply(boot.indicies,2,function(x){data.matrix[x,]})
    
    
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
    
    ################## IPW, SW, OW
    
    PS.fit =  glm(form,data=data.matrix,family="binomial")
    ipw.PS.est = predict(PS.fit,type="response") 
    stable.fit = glm(form2,data=data.matrix,family="binomial")
    stable.est = predict(stable.fit,type="response") 

    #calculate ipw/overlap/stable weights
    save.weights[,d] = treatment[,d]/ipw.PS.est + (1-treatment[,d])/(1-ipw.PS.est)
    save.overlap[,d] = treatment[,d]*(1- ipw.PS.est) + (1-treatment[,d])*ipw.PS.est
    save.stable[,d] = treatment[,d]*(stable.est/ipw.PS.est) + (1-treatment[,d])*((1-stable.est)/(1-ipw.PS.est))
    
    # BOOTSTRAPPING
    
    print("Make bootstrap weights for ipw/stable/ow")
    print(Sys.time() - m.time)
    
      for (w in 1:W){
        boot.treatment = treatment[boot.indicies[,w],d]
        PS.fit =  glm(form,data=cefalu.boot.data[[w]],family="binomial")
        sum.ps.fit = summary(PS.fit)
        
        PS.est = predict(PS.fit,type="response") 
        stable.est = predict(glm(form2,data=cefalu.boot.data[[w]],family="binomial"),type="response")
        
        save.boot.weights[,w] = save.boot.weights[,w]*(boot.treatment/PS.est + (1-boot.treatment)/(1-PS.est))
        save.boot.stable[,w] = save.boot.stable[,w]*(boot.treatment*(stable.est/PS.est) + (1-boot.treatment)*((1-stable.est)/(1-PS.est)))
        save.boot.overlap[,w] = save.boot.overlap[,w]*(boot.treatment*(1- PS.est) + (1-boot.treatment)*PS.est)
        
      }
      
#################### PPTA
    
    print("start ppta")
    print(Sys.time() - m.time)
    
    posteriors = as.matrix(MCMClogit(form,data=data.matrix,burnin=K/10,mcmc=K*10,thin=10))
    covars = as.matrix(cbind(rep(1,n),data.matrix[,-1]))
    cefalu.PS.est = expit(covars%*%t(posteriors)) #PS estimated from draws of alpha
    
    print("starting k iterations for ppta")
    print(Sys.time() - m.time)
    for(k in 1:K){
    #calculate cefalu membership
    cross.prob = treatment[,d]*(1-cefalu.PS.est[,k]) + (1-treatment[,d])*cefalu.PS.est[,k]
    save.membership[,k,d] = rbinom(n,1,cross.prob) 

    cefalu.size[d,k] = sum(save.membership[,k,d])
    cefalu.treated.size[d,k] = sum(save.membership[,k,d]*treatment[,d])
    cefalu.control.size[d,k] =sum(save.membership[,k,d]*(1-treatment[,d]))
    
    if(cefalu.treated.size[d,k]==0 | cefalu.control.size[d,k] == 0){
    cefalu.fail[d,k] = 1 }
    }
   
    ## BOOTSTRAP
    
    print("PPTA BOOTSTRAP")
    print(Sys.time() - m.time)
    
    cefalu.boot.posteriors = lapply(cefalu.boot.data,function(x){as.matrix(MCMClogit(form,data=x,burnin=K/10,mcmc=K*10,thin=10))})
    cefalu.boot.covars = lapply(cefalu.boot.data,function(x){as.matrix(cbind(rep(1,n),x[,-1]))})
            
    for(w in 1:W){
      boot.treatment = treatment[boot.indicies[,w],d]
      cefalu.boot.PS = as.matrix(expit(cefalu.boot.covars[[w]]%*%t(cefalu.boot.posteriors[[w]])))
        
      cefalu.boot.cross.prob = apply(cefalu.boot.PS,2,function(x){boot.treatment*(1-x)+(1-boot.treatment)*x})
      xx = apply(cefalu.boot.cross.prob,2,function(x){rbinom(n,1,x)})
      
      for(k in 1:K){
        
        cefalu.boot.save.membership[[w]][[k]]= intersect(which(xx[,k]==1),cefalu.boot.save.membership[[w]][[k]])
      }
    }
    
    print("updating covariates")
    print(Sys.time() - m.time)
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
  treat.out = data.frame(cum.treat,outcome = currentdata$outcome)
  treat.out.list = apply(boot.indicies,2,function(x,data){data[x,]},data=treat.out)
  
    #####ip weighted regression
  print("IPW")
  
  outcome.weights = apply(save.weights,1,prod,na.rm=TRUE)
  ipw.weights = outcome.weights
  
  if(is.na(any(outcome.weights))){
    ipw.boots = NA
    
    ipw.preddelta =NA
    ipw.predse = NA
    ipw.bootse = NA
  } else{
  #dataset = data.frame(outcome,cum.treat,outcome.weights)
  des = svydesign(~0,data=treat.out,weights=outcome.weights)
  #ipw.fit = svyglm(outcome~as.factor(cum.treat),design=des)
  
  ipw.fit = svyglm(outcome~cum.treat,design=des)
  #ipw.preddelta[m,] = summary(ipw.fit)$coef[2:(D+1),1]
  #ipw.predse[m,] = summary(ipw.fit)$coef[2:(D+1),2]
  
  coef.fxn = function(x,weights,data){
    weights = weights[,x]
    dat = data[[x]]
    des = svydesign(~0,data=dat,weights=weights)
    
    xx = summary(svyglm(outcome~cum.treat,design=des))$coef
    if(dim(xx)[1]>1){
      return(xx[2,1])
    }else{
      return(NA)
    }
  }
  
  ipw.boots = apply(matrix(1:W,nrow=1),2,coef.fxn,weights=save.boot.weights,data=treat.out.list)

  
  ipw.preddelta = summary(ipw.fit)$coef[2,1]
  ipw.predse = summary(ipw.fit)$coef[2,2]
  ipw.bootse = sd(ipw.boots,na.rm=TRUE)
  }
  
  
  ##################
  #overlap weighted regression
  print("overlap")
  
  outcome.weights = apply(save.overlap,1,prod,na.rm=TRUE)
  overlap.weights = outcome.weights
  
  if(any(is.na(outcome.weights))){
    overlap.boots = NA
    overlap.preddelta = NA
    overlap.predse =NA
    overlap.bootse = NA
   } else{
  #dataset = data.frame(outcome,cum.treat,outcome.weights)
  des = svydesign(~0,data=treat.out,weights=outcome.weights)
  
  overlap.fit = svyglm(outcome~cum.treat,design=des)
  
  overlap.boots = apply(matrix(1:W,nrow=1),2,coef.fxn,weights=save.boot.overlap,data=treat.out.list)

  overlap.preddelta = summary(overlap.fit)$coef[2,1]
  overlap.predse = summary(overlap.fit)$coef[2,2]
  overlap.bootse = sd(overlap.boots,na.rm=TRUE)
  }
  ###############################
  #stablized weighting
  print("stabilized")
  
  outcome.weights = apply(save.stable,1,prod,na.rm=TRUE)
  stabilized.weights = outcome.weights
  
  if(any(is.na(outcome.weights))){
    stable.boots = NA
    stabilized.preddelta = NA
    stabilized.predse = NA
    stabilized.bootse = NA
    }else{
  #dataset = data.frame(outcome,cum.treat,outcome.weights)
  des = svydesign(~0,data=treat.out,weights=outcome.weights)

  stable.fit2 = svyglm(outcome~cum.treat,design=des)
  
  stable.boots = apply(matrix(1:W,nrow=1),2,coef.fxn,weights=save.boot.stable,data=treat.out.list)

  stabilized.preddelta = summary(stable.fit2)$coef[2,1]
  stabilized.predse = summary(stable.fit2)$coef[2,2]
  stabilized.bootse = sd(stable.boots,na.rm=TRUE)}
  ####################
  #cefalu
  print("PPTA")
  
  mean.boots = NA
  once.matrix = matrix(NA,nrow=n,ncol=K)
  cefalu.fail = rep(0,K)
  cefalu.final.size = rep(NA,K)
  cefalu.fit = rep(NA,K)
  cefalu.boot.fit = matrix(NA,nrow = W,ncol = K)
  
  listtovector = function(sublist,n){
    result.vector = rep(0,n)
    result.vector[sublist] = 1
    result.vector
  }
  
  for(k in 1:K){
   
    #estimate ATE
    k.save.membership = save.membership[,k,]
    
    is.in = apply(k.save.membership,1,prod,na.rm=TRUE)
    cefalu.final.size[k] = sum(is.in)
    
    if (sum(is.in)>10){
      
      outcome.in = outcome[is.in==1]
      cum.treat.in = cum.treat[is.in==1]
      once.matrix[,k] = is.in
      cefalu.fit[k] = coef(lm(outcome.in ~ cum.treat.in))[2]
      
      #boot.is.in = matrix(unlist(lapply(cefalu.boot.save.membership[[k]],FUN=listtovector,n=n)), ncol = W, byrow = TRUE)
      #cefalu.boot.fit[,k] = apply(boot.is.in,2,function(x){if(sum(x)>10){coef(lm(boot.outcome[x==1] ~ boot.cum.treat[x==1]))[2]} else{NA}})
      
      } else {
        cefalu.fit[k] = NA
        cefalu.fail[k] = 1} 
    

    
  }
  
  for(w in 1:W){
    #estimate ATE
    boot.outcome = outcome[boot.indicies[,w]]
    boot.cum.treat = cum.treat[boot.indicies[,w]]
    
    boot.is.in = matrix(unlist(lapply(cefalu.boot.save.membership[[w]],FUN=listtovector,n=n)), ncol = K, byrow = FALSE)
    cefalu.boot.fit[w,] = apply(boot.is.in,2,function(x){if(sum(x)>10){coef(lm(boot.outcome[x==1] ~ boot.cum.treat[x==1]))[2]} else{NA}})
    
  }
  
  
  cefalu.fail.perc = mean(cefalu.fail)
  cefalu.marginal = apply(once.matrix,1,mean,na.rm=TRUE)
  cefalu.once = mean(apply(once.matrix,1,max,na.rm=TRUE),na.rm=TRUE)
  cefalu.preddelta = NA
  cefalu.predse = NA
  cefalu.boot.predse = NA
  mean.boots = NA
  
if(any(!is.na(cefalu.fit))){
    
    cefalu.preddelta = mean(cefalu.fit,na.rm=TRUE)
    cefalu.predse = sd(apply(cefalu.boot.fit,1,mean,na.rm=TRUE),na.rm=TRUE) }
    
if(any(!is.na(cefalu.boot.fit))){
    cefalu.boot.predse = sd(cefalu.boot.fit,na.rm=TRUE)
    mean.boots = apply(cefalu.boot.fit,1,mean,na.rm=TRUE) } 
  
  
  
#result[[m]] =
    return.thing = list(list(ipw.preddelta,ipw.predse,ipw.bootse,ipw.weights),list(cefalu.preddelta,cefalu.predse,cefalu.boot.predse,
        cefalu.marginal,cefalu.once,cefalu.size,
        cefalu.treated.size,cefalu.control.size,cefalu.fail.perc,mean.boots,cefalu.boot.fit,cefalu.final.size),list(overlap.preddelta,
      overlap.predse,overlap.bootse,overlap.weights,overlap.boots),list(stabilized.preddelta,stabilized.predse,stabilized.bootse,stabilized.weights))
    
    save(return.thing,file=paste("/n/home10/silverfoxflower/April19/result_",data.type,"_",m,".R",sep=""))
    return(return.thing)
},savedata=savedata,data.type=data.type)

#return(result)
}

