### Paper 2 graph code

########################################

# Table 4: Average baseline covariates at various levels of exposure, and post-analysis
#contribution status to the COS

load('~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/cut72_data.R')
fixed.X = cut72.data$fixed.X
cum.treat = apply(cut72.data$treatment,1,sum)
fixed.X$cum.treat = cum.treat

for(i in 1:length(names(fixed.X))){
  form = as.formula(paste(names(fixed.X)[i],"~cum.treat",sep=""))
  print(aggregate(form,fixed.X,FUN=mean))
}
#################################################################

# Table 5: Estimates of the causal rate ratio of an additional season of high exposure
#to air pollution on IHD admissions rates (per 10,000 person-years) at the zip-code
#level

modify.result = function(result){
  ipw_result = result[[1]]
  ppta_result = result[[2]]
  overlap_result = result[[3]]
  stabilized_result = result[[4]]
  
  result_frame = data.frame(method = c("IPW","PPTA","OW","SW"),
                            TE = c(ipw_result[[1]],ppta_result[[1]],overlap_result[[1]],stabilized_result[[1]]),
                            SE = c(ipw_result[[3]],ppta_result[[2]],overlap_result[[3]],stabilized_result[[3]])
                            #SE.type = c(rep("Asymp",4),rep("Boot",4))
  )
  
  
  result_frame$lower = exp(result_frame$TE - 1.96*result_frame$SE)
  result_frame$upper = exp(result_frame$TE+1.96*result_frame$SE) 
  result_frame$TE = exp(result_frame$TE)
  
  return(result_frame)}

load("~/shared_space/ci3_analysis/zigler_lab/projects/Time_varying_PPTA/final_app_result")
result_frame = modify.result(result)

library(xtable)
xtable(result_frame,digits=5)

######################################################################

# Table 7: exposure pattern prevalence after dichotomization

patterns = c()
count_patterns = data.frame(pattern=NA,count=NA)
count_patterns$count = as.numeric(count_patterns$count)
for (t in 1:dim(treatment)[1]){
  t_chara = paste(treatment[t,],collapse="") 
  
  if(t_chara%in% patterns){
    count_patterns[which(count_patterns$pattern==t_chara),2] = count_patterns[which(count_patterns$pattern==t_chara),2] + 1
  } else{
    patterns = c(patterns,t_chara)
    count_patterns = rbind(count_patterns,c(t_chara,1))
    count_patterns$count = as.numeric(count_patterns$count)
  }
}

count_patterns = count_patterns[-1,]
count_patterns

library(xtable)
xtable(count_patterns)
