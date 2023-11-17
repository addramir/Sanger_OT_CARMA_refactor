setwd("~/Projects/Sanger_OT_CARMA_refactor/01_test/")
load("Zero_MCH_run.RData")
library(CARMA)
library(Matrix)
#source("CARMA_MCShotgun.R")
source("new_MCS_fun.R")
####
z=z
ld.matrix=ld
epsilon=epsilon.list
Max.Model.Dim=Max.Model.Dim
lambda = 1
outlier.switch=outlier.switch
tau=tau
num.causal = num.causal
y.var=y.var
label = "l"
output.labels = output.labels
effect.size.prior=effect.size.prior
inner.all.iter = all.inner.iter

#input.conditional.S.list = all.C.list[[4]]
C.list = NULL
input.conditional.S.list = NULL

#### START

ridge.fun<-function(x){
  temp.ld.S<-x*modi.ld.S+(1-x)*diag(nrow(modi.ld.S))
  temp.Sigma[test.S,test.S]<-temp.ld.S
  return( outlier_likelihood(test.S,temp.Sigma,z,outlier.tau,length(test.S),1) )
} 

#######The prior distributions on the model space#########
p<-length(z)

#Poisson.prior.dist
prior.dist<-function(t){
  l=unlist(strsplit(t,split=","))
  dim.model<-length(l)
  result<-dim.model*log(lambda)+lfactorial(p-dim.model)-lfactorial(p)
  return(result)
}
marginal_likelihood=ind_Normal_fixed_sigma_marginal
tau.sample<-tau
if(outlier.switch){
  outlier_likelihood=outlier_ind_Normal_marginal
  outlier.tau=tau.sample
}


########Feature learning for the fine-mapping step, such as learning the visited model space from the previous iterations#######

B<-Max.Model.Dim
stored.result.prob<-rep(0,p)
stored.bf<-0
Sigma<-as.matrix(ld.matrix)


S<-NULL

null.model<-""
null.margin<-prior.dist(null.model)
if(is.null(C.list)){
  C.list<-vector("list",2)
}
B.list<-list()
B.list[[1]]<-prior.dist(null.model)
B.list[[2]]<-""
if(length(input.conditional.S.list)==0){
  conditional.S.list<-list()
  conditional.S=NULL
}else{
  conditional.S.list<-input.conditional.S.list
  conditional.S<-input.conditional.S.list$Index
  conditional.S<-unique(conditional.S)
  S<-conditional.S
}

for(l in 1:inner.all.iter){
  for(h in 1:10){
    ##############Shotgun COMPUTATION ############
    
    set.gamma<-set.gamma.func(S,conditional.S,p)  
    
    if(is.null(conditional.S)){
      working.S=S
    }else{
      working.S=S[-match(conditional.S,S)]
    }
    
    set.gamma.margin<-list()
    set.gamma.prior<-list()
    matrix.gamma<-list()
    
    if(length(working.S)!=0){
      S.model<-paste0(sort(working.S),collapse = ",")
      p_S=length(working.S)
      current.log.margin<-marginal_likelihood(working.S,Sigma,z,tau=tau.sample,p_S=p_S,y.var)+prior.dist(S.model)
    }else{
      current.log.margin<-prior.dist(null.model)
    }
    
    for(i in  1:length(set.gamma)){
      if (length(set.gamma[[i]])>0){
        matrix.gamma[[i]]<-index.fun(set.gamma[[i]],p=p)
        p_S=dim(set.gamma[[i]])[2]
        set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
        set.gamma.prior[[i]]<-sapply(matrix.gamma[[i]],prior.dist)
        set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
      } else {
        set.gamma.margin[[i]]=null.margin
        set.gamma.prior[[i]]=0
        matrix.gamma[[i]]=null.model
      }
    }
    add.B<-vector("list",2)
    for(i in 1:length(set.gamma)){
      if (!is.null(set.gamma.margin[[i]])){
        add.B[[1]]<-c(add.B[[1]],set.gamma.margin[[i]])
        add.B[[2]]<-c(add.B[[2]],matrix.gamma[[i]])
      }
    }
    names(add.B[[1]])=NULL
    names(add.B[[2]])=NULL
    ########## add visited models into the storage space of models###############
    
    
    B.list[[1]]<-c((B.list[[1]]),(add.B[[1]]))
    B.list[[2]]<-c(B.list[[2]],add.B[[2]])
    
    ind1=B.list[[2]]
    ind1=which(duplicated(ind1))
    if (length(ind1)>0){
      B.list[[1]]=B.list[[1]][-ind1]
      B.list[[2]]=B.list[[2]][-ind1]
    }
    
    B.list[[2]]<-B.list[[2]][order(B.list[[1]],decreasing = T)]
    B.list[[1]]<-B.list[[1]][order(B.list[[1]],decreasing = T)]
    
    ###################Select next visiting model###############
    
    if(length(working.S)!=0){
      set.star<-data.frame(set.index=1:3,gamma.set.index=rep(NA,3),margin=rep(NA,3))
      for(i in 1:3){
        aa<-set.gamma.margin[[i]]-current.log.margin
        aa<-aa-aa[which.max(aa)]
        if(length(which(is.nan(aa)))!=0){
          aa[which(is.nan(aa))]<-min(aa)
        }
        set.star$gamma.set.index[i] <-sample(1:length(set.gamma.margin[[i]]),1,prob=exp(aa))
        set.star$margin[i]<-set.gamma.margin[[i]][set.star$gamma.set.index[i]]
      }
      
      #######The Bayesian hypothesis testing for Z-scores/LD discrepancies########
      if(outlier.switch){
        for(i in 2:length(set.gamma)){
          repeat{
            
            aa<-set.gamma.margin[[i]]-current.log.margin
            aa<-aa-aa[which.max(aa)]
            if(length(which(is.nan(aa)))!=0){
              aa[which(is.nan(aa))]<-min(aa[!is.na(aa)])
            }
            
            set.star$gamma.set.index[i]<-c(sample((1:length(set.gamma.margin[[i]])),
                                                  1,prob=exp(aa)))
            set.star$margin[i]<-set.gamma.margin[[i]][  set.star$gamma.set.index[i]]
            
            test.S<-set.gamma[[i]][set.star$gamma.set.index[i],]
            
            modi.Sigma=temp.Sigma=Sigma
            if(length(test.S)>1){
              
              modi.ld.S<- modi.Sigma[test.S,test.S]
              
              opizer<-optimize(ridge.fun,interval=c(0,1),maximum = T)
              modi.ld.S<-opizer$maximum*modi.ld.S+(1-opizer$maximum)*diag(nrow(modi.ld.S)) 
              
              
              modi.Sigma[test.S,test.S]<-modi.ld.S
              
              test.log.BF<-outlier_likelihood(test.S,Sigma,z,outlier.tau,length(test.S),1)-outlier_likelihood(test.S,modi.Sigma,z,outlier.tau,length(test.S),1)
              test.log.BF<--abs(test.log.BF)

            }
            
            if(exp(test.log.BF)<outlier.BF.index){
              set.gamma[[i]]<-set.gamma[[i]][-set.star$gamma.set.index[i],]
              set.gamma.margin[[i]]<-set.gamma.margin[[i]][-set.star$gamma.set.index[i]]
              conditional.S<-c(conditional.S,setdiff(test.S,working.S))
              conditional.S<-unique(conditional.S)
            }else{
              break
            }
          }
        }
      }

      if(length(working.S)==num.causal){
        set.star<-set.star[-2,]
        aa<-set.star$margin-current.log.margin-max(set.star$margin-current.log.margin)
        sec.sample<-sample(c(1,3),1,prob=exp(aa))
        S<-set.gamma[[sec.sample]][set.star$gamma.set.index[[which(sec.sample==set.star$set.index)]] ,]
      }else{
        aa<-set.star$margin-current.log.margin-max(set.star$margin-current.log.margin)
        sec.sample<-sample(1:3,1,prob=exp(aa) )
        S<-set.gamma[[sec.sample]][set.star$gamma.set.index[[sec.sample]] ,]
      }
      
    }else{
      set.star<-data.frame(set.index=rep(1,3),gamma.set.index=rep(NA,3),margin=rep(NA,3))
      aa<-set.gamma.margin[[2]]-current.log.margin
      aa<-aa-aa[which.max(aa)]
      if(length(which(is.nan(aa)))!=0){
        aa[which(is.nan(aa))]<-min(aa)
      }
      
      set.star$gamma.set.index[2] <-c(sample((1:length(set.gamma.margin[[2]]))[order(exp(aa),decreasing = T)[1:(min(length(aa),floor(p/2)))]],
                                             1,prob=exp(aa)[order(exp(aa),decreasing = T)[1:(min(length(aa),floor(p/2)))]]))
      set.star$margin[2]<-set.gamma.margin[[2]][  set.star$gamma.set.index[2]]
      
      S<-set.gamma[[2]][set.star$gamma.set.index[2],]
      
    }
    
    S<-unique(c(S,conditional.S))
    
    
    
  }
  ######Output of the results of the module function######
  
  
  
  if(!is.null(conditional.S)){
    all.c.index<-c()
    for(tt in conditional.S){
      l=strsplit(B.list[[2]],split=",")
      ind=which(sapply(l,FUN=function(x){tt%in%x}))
      all.c.index<-c(all.c.index,ind)
    }
  
    all.c.index<-unique(all.c.index)
    
    temp.B.list<-B.list
    temp.B.list[[1]]<-B.list[[1]][-all.c.index]
    temp.B.list[[2]]<-B.list[[2]][-all.c.index]
  }else{
    temp.B.list<-B.list
  }
  result.B.list<-list()
  result.B.list[[1]]<-temp.B.list[[1]][(1:min(B,length(temp.B.list[[2]])))]
  result.B.list[[2]]<-temp.B.list[[2]][(1:min(B,length(temp.B.list[[2]])))]
  
  result.prob=NULL
  
  conditional.S.list<-data.frame(Index=conditional.S,Z=z[conditional.S])
 
  rb1=result.B.list[[1]]
  
  difference<-abs(mean(rb1[1:round(length(rb1)/4)])-stored.bf)
  
  if(difference<epsilon){
    break
  }else{
    stored.bf<-mean(rb1[1:round(length(rb1)/4)])
  }
}
return(list(result.B.list,C.list,result.prob,conditional.S.list))



