###The function that defines neighborhood model space
set.gamma.func<-function(input.S,condition.index=NULL,p){
  
  #S - vector of causal variants in index configuration
  #condition.index - outliers
  
  add.function<-function(S_sub,y){
    #this fucntion combines S_sub and y and sorts each row
    results<-apply(as.matrix(S_sub),1,function(x){return(sort(c(x,y)))})
    results=t(results)
    return(results)
  }
  
  set.gamma.func.base<-function(S,p){
    
    set.gamma<-list()
    
    #set of gamma-
    if(length(S)==0){
      set.gamma[[1]]<-NULL
      set.gamma[[2]]<-as.matrix(1:p)
      set.gamma[[3]]<-NULL
    }
    
    if(length(S)==1){
      S_sub<-(1:p)[-S]
      set.gamma[[1]]<-NULL
      set.gamma[[2]]<-add.function(S_sub,S)
      set.gamma[[3]]<-as.matrix(S_sub)
    }
    
    if(length(S)>1){
      S_sub<-(1:p)[-S]
      
      S=sort(S)
      set.gamma[[1]]=t(combn(S,length(S)-1))
      
      #set of gamma+
      set.gamma[[2]]<-add.function(S_sub,S)
      
      #set of gamma=
      xs<-NULL
      for(i in 1:nrow(set.gamma[[1]])){
        xs<-rbind(xs,add.function(S_sub,set.gamma[[1]][i,]))  
      }
      set.gamma[[3]]=xs
      
    }
    
    return(set.gamma)
  }
  
  set.gamma.func.conditional<-function(input.S,condition.index,p){
    
    set.gamma<-list()
    S=input.S[-match(condition.index,input.S)]
    
    #set of gamma-
    if(length(S)==0){
      S_sub<-(1:p)[-condition.index]
      set.gamma[[1]]<-NULL
      set.gamma[[2]]<-as.matrix(S_sub)
      set.gamma[[3]]<-NULL
    }
    
    if(length(S)==1){
      S_sub<-(1:p)[-input.S]
      set.gamma[[1]]<-NULL
      set.gamma[[2]]<-add.function(S_sub,S)
      set.gamma[[3]]<-as.matrix(S_sub)
    }
    
    if(length(S)>1){
      S_sub<-(1:p)[-input.S]
      
      S=sort(S)
      set.gamma[[1]]=t(combn(S,length(S)-1))
      
      #set of gamma+
      set.gamma[[2]]<-add.function(S_sub,S)
      
      #set of gamma=
      xs<-NULL
      for(i in 1:nrow(set.gamma[[1]])){
        xs<-rbind(xs,add.function(S_sub,set.gamma[[1]][i,]))  
      }
      set.gamma[[3]]=xs
      
      
    }
    
    return(set.gamma)
  }
  
  if(is.null(condition.index)){
    results<-set.gamma.func.base(input.S,p)
  }else{
    results<-set.gamma.func.conditional(input.S,condition.index,p)
  }
  return(results)
}

match.dgCMatrix <- function (dgCMat1,dgCMat2) {
  #  dgCMat1=B.list[[2]]
  # dgCMat2=add.B[[2]]
  n1 <- nrow(dgCMat1)
  p1 <- ncol(dgCMat1)
  J1 <- rep(1:p1, diff(dgCMat1@p))
  I1 <- dgCMat1@i + 1
  x1 <- dgCMat1@x
  n2 <- nrow(dgCMat2)
  p2 <- ncol(dgCMat2)
  J2 <- rep(1:p2, diff(dgCMat2@p))
  I2 <- dgCMat2@i + 1
  x2 <- dgCMat2@x
  ## check duplicated rows
  names(x1) <- J1
  RowLst1 <- split(J1, I1)
  is_empty1<- setdiff(1:n1, I1)
  ## check duplicated rows
  names(x2) <- J2
  RowLst2 <- split(J2, I2)
  is_empty2 <- setdiff(1:n2, I2)
  result<-(match(RowLst2,RowLst1))
  if(any(which(result>is_empty1))){
    result[which(result>=is_empty1)]=result[which(result>=is_empty1)]+1
  }
  if(any(is_empty1)){
    if(any(is_empty2)){
      result<-c(is_empty1,result)  
    }
  }else{
    if(any(is_empty2)){
      result<-c(NA,result)  
    }
  }
  return(result)
  
}
####Function that computes posterior inclusion probability based on the marginal likelihood and model space
PIP.func<-function(likeli,model.space,p){
  infi.index<-which(is.infinite(likeli))
  if(length(infi.index)!=0){
    likeli<-likeli[-infi.index]
    model.space<-model.space[-infi.index,]
  }
  na.index<-which(is.na(likeli))
  if(length(na.index)!=0){
    likeli<-likeli[-na.index]
    model.space<-model.space[-na.index,]
  }
  aa<-likeli-max(likeli,na.rm=T)
  prob.sum<-sum(exp(aa))
  result.prob<-rep(NA,p)
  for(i in 1:p){
    result.prob[i]<-sum(exp(aa[ which(model.space[,i]==1)]))/prob.sum
  }
  return(result.prob)
}

index.fun<-function(x,p){
  m=as(matrix(0,nrow=nrow(x),ncol=p), "TsparseMatrix")
  m@i<-as.integer(rep(1:nrow(x)-1,each=ncol(x)))
  m@j<-as.integer(c(t(x))-1)
  m@x=rep(1,nrow(x)*ncol(x))
  m<-as(m,"CsparseMatrix")
  return(m)
}




#########The module function of the CARMA fine-mapping step for each locus included in the analysis##########

MCS_modified<-function(z,ld.matrix,Max.Model.Dim=1e+4,lambda,label,
                                num.causal=10,output.labels,y.var=1,
                                effect.size.prior=effect.size.prior,
                                outlier.switch,input.conditional.S.list=NULL,tau=1/0.05^2,
                                epsilon=1e-3,inner.all.iter=10){
  
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
      }else{
        working.S=S[-match(conditional.S,S)]
      }
      
      set.gamma.margin<-list()
      set.gamma.prior<-list()
      matrix.gamma<-list()
      
      if(length(working.S)!=0){
        S.model<-as(Matrix(nrow=1,ncol=p,sparse = T,data=0),'CsparseMatrix')
        S.model[,working.S]<-1
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
          set.gamma.prior[[i]]<-apply(matrix.gamma[[i]],1,prior.dist)
          set.gamma.margin[[i]]=set.gamma.prior[[i]]+set.gamma.margin[[i]]
        } else {
          set.gamma.margin[[i]]=null.margin
          set.gamma.prior[[i]]=NULL
          matrix.gamma[[i]]=null.model
        }
      }
      add.B<-list()
      add.B[[1]]<-NULL
      add.B[[2]]<-as(Matrix(nrow=0,ncol=p,sparse = T,data=0),'CsparseMatrix')
      for(i in 1:length(set.gamma)){
        if (!is.null(set.gamma.margin[[i]])){
          add.B[[1]]<-c(add.B[[1]],set.gamma.margin[[i]])
          add.B[[2]]<-rbind(add.B[[2]],matrix.gamma[[i]])
        }
      }
      ########## add visited models into the storage space of models###############
      
      
      add.index<-match.dgCMatrix(B.list[[2]],add.B[[2]])
      
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
      single.marginal<-apply(single.set,1,marginal_likelihood,Sigma=Sigma,z=z,tau=tau.sample,p_S=1,y_sigma=y.var)
      aa<-single.marginal-max(single.marginal,na.rm=T)
      
      prob.sum<-sum(exp(aa))
      result.prob<-(exp(aa))/prob.sum
    }else{
      result.prob=PIP.func(result.B.list[[1]],result.B.list[[2]],p)
    }
    conditional.S.list<-data.frame(Index=conditional.S,Z=z[conditional.S])
    
    rb1=result.B.list[[1]]
    
    difference<-abs(mean(rb1[1:round(length(rb1)/4)])-stored.bf)
    
    if(difference<epsilon){
      break
    }else{
      stored.bf<-mean(rb1[1:round(length(rb1)/4)])
    }
  }
  
  out=list(result.B.list,C.list=NULL,result.prob,conditional.S.list)
  return(out)
}
