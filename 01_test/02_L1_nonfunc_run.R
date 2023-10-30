setwd("~/Projects/Sanger_OT_CARMA_refactor/01_test/")
load("test.RData")
library(CARMA)
library(Matrix)

source("CARMA_funcs_L1.R")

#ASSUMPTIONS 
w.list=NULL
effect.size.prior='Spike-slab'
EM.dist='Logistic'
prior.prob.computation='Logistic'
label.list=NULL

#Init of variables
output.labels='.'

rho.index=0.99
BF.index=10
Max.Model.Dim=2e+5
all.iter=3
all.inner.iter=10
input.alpha=0
epsilon.threshold=1e-5
printing.log=F
num.causal=10
y.var=1
tau=0.04
outlier.switch=T
outlier.BF.index=1/3.2

model.prior='Poisson'

###input
z=z.list[[2]]
p_snp=length(z)
ld=ld.list[[2]]

epsilon.list<-epsilon.threshold*p_snp
all.epsilon.threshold<-epsilon.threshold*p_snp

######## Burning step###########
previous.result<-list()
########Run fine-mapping step (module function) for each locus included in the analysis
t0=Sys.time()
all.C.list<-Module.Cauchy.Shotgun(z,ld,epsilon=epsilon.list,
                           Max.Model.Dim=Max.Model.Dim,lambda = 1,
                           outlier.switch=outlier.switch,tau=tau,
                           num.causal = num.causal,y.var=y.var,
                           label = "l",output.labels = output.labels,
                           effect.size.prior=effect.size.prior,
                           model.prior=model.prior,inner.all.iter = all.inner.iter)
t1=Sys.time()-t0
print(paste0('This is locus burning time'))
print((t1))
  
########Running CARMA######## 
for(g in 1:all.iter){ 
 
  delete.list<-integer(0)
  if(outlier.switch & (nrow(all.C.list[[4]])!=0)){
      delete.list<-all.C.list[[4]]$Index
  }
  
  ac1=all.C.list[[1]][[1]]
  previous.result<-mean(ac1[1:round(length(ac1)/4)])
  
  prior.prob.list<-list(NULL)
  
  #######Fine-mapping step for each locus, i.e., the E-step in the EM algorithm
  t0=Sys.time()
  all.C.list<-Module.Cauchy.Shotgun(z=z,ld,input.conditional.S.list = all.C.list[[4]],
                                         Max.Model.Dim=Max.Model.Dim,y.var=y.var,num.causal = num.causal,epsilon=epsilon.list,
                                         C.list = all.C.list[[2]],prior.prob = prior.prob.list,
                                         outlier.switch=outlier.switch,tau=tau,
                                         lambda = 1, label = "l",output.labels = output.labels,
                                         effect.size.prior=effect.size.prior,model.prior=model.prior,inner.all.iter = all.inner.iter)
  t1=Sys.time()-t0
  print(paste0('This is locus computing time'))
  print((t1))

  ac1=all.C.list[[1]][[1]]
  difference<-abs(previous.result-mean(ac1[1:round(length(ac1)/4)]))
  print(paste0('This is difference; ',difference))
  if(difference<all.epsilon.threshold){
    break
  }
}

results.list<-list()
pip=all.C.list[[3]]
credible.set<-credible.set.fun.improved(pip,ld,rho=rho.index)
credible.model<-credible.model.fun(all.C.list[[1]][[1]],all.C.list[[1]][[2]],bayes.threshold = BF.index)
results.list[[1]]<-pip
results.list[[2]]<-credible.set
results.list[[3]]<-credible.model
results.list[[4]]<-all.C.list[[4]]

results.list[[5]]<-all.C.list

names(results.list)<-c('PIPs','Credible_set','Credible_model','Outliers','all.C.list')
  
#PIPs
#[1] 0.89254484 0.33198906 0.07088181 0.48519930 0.38178145 0.00000000 0.04397452
#[8] 0.04137350 0.03143955 0.03647640 0.02621831 0.03552681 0.03274398 0.02607696
#[15] 0.03059451 0.03359447 0.01415915 0.02176282 0.01689036 0.01430290 0.01447391




