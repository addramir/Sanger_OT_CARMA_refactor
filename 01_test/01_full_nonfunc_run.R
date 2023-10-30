setwd("~/Projects/Sanger_OT_CARMA_refactor/01_test/")
load("test.RData")
library(CARMA)
library(Matrix)

source("CARMA_funcs.R")

#Init of variables

w.list=NULL
output.labels='.'
label.list=NULL
effect.size.prior='Spike-slab'
rho.index=0.99
BF.index=10
EM.dist='Logistic'
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
prior.prob.computation='Logistic'

##########First part
log.2pi<-log(2*pi)
L<-length(z.list)
p.list<-list()
for(i in 1:L){
  z.list[[i]]<-as.matrix(z.list[[i]])
  p.list[[i]]<-nrow(z.list[[i]])
}
B<-Max.Model.Dim
all.B.list<-list()
for(i in 1:L){
  all.B.list[[i]]<-list()
  all.B.list[[i]][[1]]<-integer(0)
  all.B.list[[i]][[2]]<-Matrix(nrow = 0,ncol=p.list[[i]],data=0,sparse = T)
}
q.list<-list()
if(!is.null(w.list)){
  for(i in 1:L){
    
    q.list[[i]]<-ncol(w.list[[i]])
    invariant.var.index<-which((apply(w.list[[i]][,-1],2,sd))==0)
    if(length(invariant.var.index)!=0){
      invariant.var<-w.list[[i]][,invariant.var.index+1]
      w.list[[i]]<-as.matrix(cbind(1,scale(w.list[[i]][,-1])))
      w.list[[i]][,invariant.var.index+1]<-invariant.var
    }else{
      w.list[[i]]<-as.matrix(cbind(1,scale(w.list[[i]][,-1])))
    }
    
  }
}
if(is.null(label.list)){
  for(i in 1:L){
    label.list[[i]]=paste0('locus_',i)
  }
}
Sigma.list<-list()
for(i in 1:L){
  Sigma.list[[i]]<-as.matrix(ld.list[[i]])
}
S.list<-list()
for(i in 1:L){
  S.list[[i]]<-integer(0)
}

all.C.list<-list()
for(i in 1:L){
  all.C.list[[i]]<-list()
  all.C.list[[i]][[1]]<-integer(0)
  all.C.list[[i]][[2]]<-Matrix(nrow = 0,ncol=p.list[[i]],data=0,sparse = T)
}

all.epsilon.threshold<-0
epsilon.list<-list()
for(i in 1:L){
  epsilon.list[[i]]<-epsilon.threshold*p.list[[i]]
  all.epsilon.threshold<-all.epsilon.threshold+  epsilon.threshold*p.list[[i]]
}
model.prior='Poisson'
standardize.model.space=T

######## Burning step###########
previous.result<-list()
########Run fine-mapping step (module function) for each locus included in the analysis
for(i in 1:L){
  t0=Sys.time()
  all.C.list[[i]]<-Module.Cauchy.Shotgun(z.list[[i]],ld.list[[i]],epsilon=epsilon.list[[i]],
                                         Max.Model.Dim=Max.Model.Dim,lambda = lambda.list[[i]],
                                         outlier.switch=outlier.switch,tau=tau,
                                         num.causal = num.causal,y.var=y.var,
                                         label = label.list[[i]],output.labels = output.labels,
                                         effect.size.prior=effect.size.prior,model.prior=model.prior,inner.all.iter = all.inner.iter)
  t1=Sys.time()-t0
  print(paste0('This is locus ',i,' burning time'))
  print((t1))
  
}
########Running CARMA######## 
for(g in 1:all.iter){ 
  if(outlier.switch){
    delete.list<-list()
    for(i in 1:L){
      delete.list[[i]]<-integer(0)
      if(nrow(all.C.list[[i]][[4]])!=0){
        temp.delete.list<-c(all.C.list[[i]][[4]]$Index)
        delete.list[[i]]<-temp.delete.list
      }
    }
  }else{
    delete.list<-list()
    for(i in 1:L){
      delete.list[[i]]<-integer(0)
    }
  }
  
  for(i in 1:L){
    previous.result[[i]]<-mean(all.C.list[[i]][[1]][[1]][1:round(quantile(1:length(all.C.list[[i]][[1]][[1]]),probs = 0.25))])
  }
  prior.prob.list<-list()
  for(i in 1:L){
    prior.prob.list[[i]]<-list(NULL)
  }
  #######Fine-mapping step for each locus, i.e., the E-step in the EM algorithm
  for(i in 1:L){
    t0=Sys.time()
    all.C.list[[i]]<-Module.Cauchy.Shotgun(z=z.list[[i]],ld.list[[i]],input.conditional.S.list = all.C.list[[i]][[4]],
                                           Max.Model.Dim=Max.Model.Dim,y.var=y.var,num.causal = num.causal,epsilon=epsilon.list[[i]],
                                           C.list = all.C.list[[i]][[2]],prior.prob = prior.prob.list[[i]],
                                           outlier.switch=outlier.switch,tau=tau,
                                           lambda = lambda.list[[i]], label = label.list[[i]],output.labels = output.labels,
                                           effect.size.prior=effect.size.prior,model.prior=model.prior,inner.all.iter = all.inner.iter)
    t1=Sys.time()-t0
    print(paste0('This is locus ',i,' computing time'))
    print((t1))
  }
  
  difference<-0
  for(i in 1:L){
    difference<-difference+abs(previous.result[[i]]-mean(all.C.list[[i]][[1]][[1]][1:round(quantile(1:length(all.C.list[[i]][[1]][[1]]),probs = 0.25))]))
  }
  print(paste0('This is difference; ',difference))
  if(difference<all.epsilon.threshold){
    break
  }
}
results.list<-list()
for(i in 1:L){
  results.list[[i]]<-list()
  pip=all.C.list[[i]][[3]]
  credible.set<-credible.set.fun.improved(pip,ld.list[[i]],rho=rho.index)
  credible.model<-credible.model.fun(all.C.list[[i]][[1]][[1]],all.C.list[[i]][[1]][[2]],bayes.threshold = BF.index)
  results.list[[i]][[1]]<-pip
  results.list[[i]][[2]]<-credible.set
  results.list[[i]][[3]]<-credible.model
  results.list[[i]][[4]]<-all.C.list[[i]][[4]]
  
  results.list[[i]][[5]]<-all.C.list[[i]]
  
  names(results.list[[i]])<-c('PIPs','Credible set','Credible model','Outliers','all.C.list')
  
}
print("MODIFIED!")





