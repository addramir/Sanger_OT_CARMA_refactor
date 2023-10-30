setwd("~/Projects/Sanger_OT_CARMA_refactor/01_test/")
load("Zero_MCH_run.RData")
library(CARMA)

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

input.conditional.S.list = NULL
C.list = NULL


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
  dim.model<-sum(t)
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

null.model<-Matrix(nrow = 1,ncol=p,data=0,sparse = T)
null.margin<-prior.dist(null.model)
if(is.null(C.list)){
  C.list<-list()
  C.list[[2]]<-list()
  C.list[[1]]<-list()
}
B.list<-list()
B.list[[1]]<-prior.dist(null.model)
B.list[[2]]<-Matrix(nrow = 1,ncol=p,data=0,sparse = T)
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
      base.model<-null.model
      base.model.margin<-null.margin
    }else{
      working.S=S[-match(conditional.S,S)]
      if(length(working.S)!=0){
        base.model<-as(Matrix(nrow=1,ncol=p,sparse = T,data=0),'CsparseMatrix')
        base.model[,working.S]<-1
        p_S=length(working.S);
        base.model.margin<-marginal_likelihood(working.S,Sigma,z,tau=tau.sample,p_S=p_S,y.var)+prior.dist(base.model)
      }else{
        base.model<-null.model
        base.model.margin<-null.margin
      }
    }
    set.gamma.margin<-list()
    set.gamma.prior<-list()
    matrix.gamma<-list()
    if(length(working.S)!=0){
      S.model<-as(Matrix(nrow=1,ncol=p,sparse = T,data=0),'CsparseMatrix')
      S.model[,working.S]<-1
      p_S=length(working.S);
      current.log.margin<-marginal_likelihood(working.S,Sigma,z,tau=tau.sample,p_S=p_S,y.var)+prior.dist(S.model)
    }else{
      current.log.margin<-prior.dist(null.model)
    }
    
    if(length(working.S)>1){
      for(i in  1:length(set.gamma)){
        t0=Sys.time()
        matrix.gamma[[i]]<-index.fun(set.gamma[[i]],p=p)
        
        if(length(C.list[[2]])<ncol(set.gamma[[i]])){
          C.list[[2]][[ncol(set.gamma[[i]])]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
          C.list[[1]][[ncol(set.gamma[[i]])]]<-integer(0)
          computed.index<-integer(0)
        }else{
          computed.index<-match.dgCMatrix(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
        }
        
        p_S=dim(set.gamma[[i]])[2]
        if(length(na.omit(computed.index))==0){
          set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
          C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]])
          C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
          set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
          set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
        }else{
          set.gamma.margin[[i]]<-rep(NA,nrow(matrix.gamma[[i]]))
          set.gamma.margin[[i]][!is.na(computed.index)]<-C.list[[1]][[ncol(set.gamma[[i]])]][na.omit(computed.index)]
          if(sum(is.na(computed.index))!=0){
            set.gamma.margin[[i]][is.na(computed.index)]<-apply(set.gamma[[i]][is.na(computed.index),,drop=F],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
          }
          C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]][is.na(computed.index)])
          C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]][is.na(computed.index),,drop=F])
          set.gamma.margin[[i]]<- set.gamma.margin[[i]]+apply(matrix.gamma[[i]],1,prior.dist)
        }
        t1=Sys.time()-t0
      }
      
      add.B<-list()
      add.B[[1]]<-c(set.gamma.margin[[1]],
                    set.gamma.margin[[2]],
                    set.gamma.margin[[3]])
      add.B[[2]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
      for(i in 1:3){
        add.B[[2]]<-rbind(add.B[[2]],matrix.gamma[[i]])
      }
      
    }
    
    if(length(working.S)==1){
      set.gamma.margin[[1]]<-null.margin
      matrix.gamma[[1]]<-null.model
      for(i in 2:3){
        
        matrix.gamma[[i]]<-index.fun(set.gamma[[i]],p=p)
        
        if(length(C.list[[2]])<ncol(set.gamma[[i]])){
          C.list[[2]][[ncol(set.gamma[[i]])]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
          C.list[[1]][[ncol(set.gamma[[i]])]]<-integer(0)
          computed.index<-integer(0)
        }else{
          computed.index<-match.dgCMatrix(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
        }
        
        
        
        p_S=dim(set.gamma[[i]])[2]
        if(length(na.omit(computed.index))==0){
          set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
          C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]])
          C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
          set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
          set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
        }else{
          set.gamma.margin[[i]]<-rep(NA,nrow(matrix.gamma[[i]]))
          set.gamma.margin[[i]][!is.na(computed.index)]<-C.list[[1]][[ncol(set.gamma[[i]])]][na.omit(computed.index)]
          if(sum(is.na(computed.index))!=0){
            set.gamma.margin[[i]][is.na(computed.index)]<-apply(set.gamma[[i]][is.na(computed.index),,drop=F],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
          }
          C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]][is.na(computed.index)])
          C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]][is.na(computed.index),,drop=F])
          set.gamma.margin[[i]]<- set.gamma.margin[[i]]+apply(matrix.gamma[[i]],1,prior.dist)
        }
      }
      
      
      add.B<-list()
      add.B[[1]]<-c(set.gamma.margin[[1]],
                    set.gamma.margin[[2]],
                    set.gamma.margin[[3]])
      add.B[[2]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
      for(i in 1:3){
        add.B[[2]]<-rbind(add.B[[2]],matrix.gamma[[i]])
      }
    }
    if(length(working.S)==0){
      
      i=2
      matrix.gamma[[i]]<-index.fun(set.gamma[[i]],p=p)
      
      if(length(C.list[[2]])<ncol(set.gamma[[i]])){
        C.list[[2]][[ncol(set.gamma[[i]])]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
        C.list[[1]][[ncol(set.gamma[[i]])]]<-integer(0)
        computed.index<-integer(0)
      }else{
        computed.index<-match.dgCMatrix(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
      }
      
      
      p_S=dim(set.gamma[[i]])[2]
      if(length(na.omit(computed.index))==0){
        set.gamma.margin[[i]]<-apply(set.gamma[[i]],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
        C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]])
        C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]])
        set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
        set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
      }else{
        set.gamma.margin[[i]]<-rep(NA,nrow(matrix.gamma[[i]]))
        set.gamma.margin[[i]][!is.na(computed.index)]<-C.list[[1]][[ncol(set.gamma[[i]])]][na.omit(computed.index)]
        if(sum(is.na(computed.index))!=0){
          set.gamma.margin[[i]][is.na(computed.index)]<-apply(set.gamma[[i]][is.na(computed.index),,drop=F],1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
        }
        C.list[[1]][[ncol(set.gamma[[i]])]]<-c(C.list[[1]][[ncol(set.gamma[[i]])]],set.gamma.margin[[i]][is.na(computed.index)])
        C.list[[2]][[ncol(set.gamma[[i]])]]<-rbind(C.list[[2]][[ncol(set.gamma[[i]])]],matrix.gamma[[i]][is.na(computed.index),,drop=F])
        set.gamma.margin[[i]]<- set.gamma.margin[[i]]+ apply(matrix.gamma[[i]],1,prior.dist)
      }
      
      
      add.B<-list()
      add.B[[1]]<-c(set.gamma.margin[[2]])
      add.B[[2]]<-matrix.gamma[[2]]
    }
    ########## add visited models into the storage space of models###############
    
    
    add.index<-match.dgCMatrix(B.list[[2]],add.B[[2]])
    if(length(which(!is.na(add.index)))>10){
      check.index<-sample(which(!is.na(add.index)),10)
      
    }
    if(length(na.omit(add.index))!=0){
      B.list[[1]]<-c((B.list[[1]]),(add.B[[1]][is.na(add.index)]))
      B.list[[2]]<-rbind(B.list[[2]],add.B[[2]][is.na(add.index),,drop=F])
    }else{
      B.list[[1]]<-c((B.list[[1]]),(add.B[[1]]))
      B.list[[2]]<-rbind(B.list[[2]],add.B[[2]])
    }
    B.list[[2]]<-B.list[[2]][order(B.list[[1]],decreasing = T),]
    B.list[[1]]<-B.list[[1]][order(B.list[[1]],decreasing = T)]
    ###################Select next visiting model###############
    
    if(length(working.S)!=0){
      set.star<-data.frame(set.index=1:3,gamma.set.index=rep(NA,3),margin=rep(NA,3))
      for(i in 1){
        aa<-set.gamma.margin[[i]]-current.log.margin
        aa<-aa-aa[which.max(aa)]
        if(length(which(is.nan(aa)))!=0){
          aa[which(is.nan(aa))]<-min(aa)
        }
        set.star$gamma.set.index[i] <-c(sample(1:length(set.gamma.margin[[i]]),1,prob=exp(aa)))
        set.star$margin[i]<-set.gamma.margin[[i]][  set.star$gamma.set.index[i]]
        rm(aa)
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
            
            modi.Sigma<-Sigma
            temp.Sigma<-Sigma
            if(length(test.S)>1){
              
              modi.ld.S<- modi.Sigma[test.S,test.S]
              
              opizer<-optimize(ridge.fun,interval=c(0,1),maximum = T)
              modi.ld.S<-opizer$maximum*modi.ld.S+(1-opizer$maximum)*diag(nrow(modi.ld.S)) 
              
              
              modi.Sigma[test.S,test.S]<-modi.ld.S
              
              test.log.BF<-outlier_likelihood(test.S,Sigma,z,outlier.tau,length(test.S),1)-outlier_likelihood(test.S,modi.Sigma,z,outlier.tau,length(test.S),1)
              test.log.BF<--abs(test.log.BF)
              print(paste0('Outlier BF: ', test.log.BF))
              print(test.S)
              print(paste0('This is xi hat: ', opizer))
            }
            
            if(exp(test.log.BF)<outlier.BF.index){
              set.gamma[[i]]<-set.gamma[[i]][-set.star$gamma.set.index[i],]
              set.gamma.margin[[i]]<-set.gamma.margin[[i]][-set.star$gamma.set.index[i]]
              conditional.S<-c(conditional.S,test.S[is.na(match(test.S,working.S))])
              conditional.S<-unique(conditional.S)
            }else{
              break
            }
          }
          rm(aa)
        }
      }else{
        for(i in 2:length(set.gamma)){
          
          aa<-set.gamma.margin[[i]]-current.log.margin
          aa<-aa-aa[which.max(aa)]
          if(length(which(is.nan(aa)))!=0){
            aa[which(is.nan(aa))]<-min(aa[!is.na(aa)])
          }
          
          set.star$gamma.set.index[i]<-c(sample((1:length(set.gamma.margin[[i]])),
                                                1,prob=exp(aa)))
          set.star$margin[i]<-set.gamma.margin[[i]][  set.star$gamma.set.index[i]]
          rm(aa)
        }
      }
      print(set.star)
      
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
      print(set.star)
    }
    print(paste0('this is running S: ',paste0(S,collapse = ',')))
    S<-unique(c(S,conditional.S))
    
    
    
  }
  ######Output of the results of the module function######
  
  
  result.B.list<-list()
  if(!is.null(conditional.S)){
    all.c.index<-c()
    
    
    for(tt in conditional.S){
      c.index<-(B.list[[2]]@i[min(length(B.list[[2]]@i),(B.list[[2]]@p[tt]+1)):B.list[[2]]@p[tt+1]])+1
      all.c.index<-c(all.c.index,c.index)
    }
    
    all.c.index<-unique(all.c.index)
    temp.B.list<-list()
    temp.B.list[[1]]<-B.list[[1]][-all.c.index]
    temp.B.list[[2]]<-B.list[[2]][-all.c.index,]
  }else{
    temp.B.list<-list()
    temp.B.list[[1]]<-B.list[[1]]
    temp.B.list[[2]]<-B.list[[2]]
  }
  result.B.list<-list()
  result.B.list[[1]]<-temp.B.list[[1]][(1:min(B,nrow(temp.B.list[[2]])))]
  result.B.list[[2]]<-temp.B.list[[2]][(1:min(B,nrow(temp.B.list[[2]]))),]
  
  if(num.causal==1){
    single.set<-matrix(1:p,p,1)
    single.marginal<-apply(single.set,1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=p_S,y_sigma=y.var)
    aa<-single.marginal-max(single.marginal,na.rm=T)
    
    prob.sum<-sum(exp(aa))
    result.prob<-(exp(aa))/prob.sum
  }else{
    result.prob=PIP.func(result.B.list[[1]],result.B.list[[2]],p)
  }
  conditional.S.list<-data.frame(Index=conditional.S,Z=z[conditional.S])
  if(!is.null(output.labels)){
    if(dir.exists(output.labels)==F ){
      dir.create(output.labels,recursive = T)
    }
    write.table(result.B.list[[1]],file=paste0(output.labels,'/post_',label,'_poi_likeli','.txt'),row.names = F,col.names = F)
    writeMM(result.B.list[[2]],file=paste0(output.labels,'/post_',label,'_poi_gamma','.mtx'))
    write.table((result.prob),file=paste0(output.labels,'/post_', label,'.txt'),row.names = F,append = F,col.names = F)
    if(outlier.switch){
      saveRDS(conditional.S.list,file=paste0(output.labels,'/post_', label,'_','outliers.RData'))
    }
  }
  
  
  difference<-abs(mean(result.B.list[[1]][1:round(quantile(1:length(result.B.list[[1]]),probs = 0.25))])-stored.bf)
  print(difference)
  if(difference<epsilon){
    break
  }else{
    stored.bf<-mean(result.B.list[[1]][1:round(quantile(1:length(result.B.list[[1]]),probs = 0.25))])
  }
}




