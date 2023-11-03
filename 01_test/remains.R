if(is.null(conditional.S)){
  working.S=S
  base.model<-null.model
  base.model.margin<-null.margin
}else{
  working.S=S[-match(conditional.S,S)]
  if(length(working.S)!=0){
    base.model<-as(Matrix(nrow=1,ncol=p,sparse = T,data=0),'CsparseMatrix')
    base.model[,working.S]<-1
    p_S=length(working.S)
    base.model.margin<-marginal_likelihood(working.S,Sigma,z,tau=tau.sample,p_S=p_S,y.var)+prior.dist(base.model)
  }else{
    base.model<-null.model
    base.model.margin<-null.margin
  }
}


if(length(working.S)>1){
  for(i in  1:length(set.gamma)){
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




duplicated.dgCMatrix <- function (dgCMat, MARGIN) {
  MARGIN <- as.integer(MARGIN)
  n <- nrow(dgCMat)
  p <- ncol(dgCMat)
  J <- rep(1:p, diff(dgCMat@p))
  I <- dgCMat@i + 1
  x <- dgCMat@x
  if (MARGIN == 1L) {
    ## check duplicated rows
    names(x) <- J
    RowLst <- split(x, I)
    is_empty <- setdiff(1:n, I)
    result <- duplicated.default(RowLst)
  } else if (MARGIN == 2L) {
    ## check duplicated columns
    names(x) <- I
    ColLst <- split(x, J)
    is_empty <- setdiff(1:p, J)
    result <- duplicated.default(ColLst)
  } else {
    warning("invalid MARGIN; return NULL")
    result <- NULL
  }
  
  if(any(is_empty)){
    out <- logical(if(MARGIN == 1L) n else p)
    out[-is_empty] <- result
    if(length(is_empty) > 1)
      out[is_empty[-1]] <- TRUE
    result <- out
  }
  
  result
}






index.fun.inner<-function(x,p){
  #m<-as(as(as(matrix(0,nrow=nrow(x),ncol=p), "TsparseMatrix"), "generalMatrix"), "TsparseMatrix")
  m=as(matrix(0,nrow=nrow(x),ncol=p), "TsparseMatrix")
  #m=matrix(0,nrow=nrow(x),ncol=p)
  m@i<-as.integer(rep(1:nrow(x)-1,each=ncol(x)))
  m@j<-as.integer(c(t(x))-1)
  m@x=rep(1,nrow(x)*ncol(x))
  m<-as(m,"CsparseMatrix")
  return(m)
}
index.fun<-function(outer.x,Max.Model.Dimins=10,p){
  if(nrow(outer.x)>1000){
    index.bins<-which((1:nrow(outer.x))%%floor(nrow(outer.x)/Max.Model.Dimins)==0)
    result.m<-index.fun.inner(outer.x[1:index.bins[1],,drop=F],p)
    for(b in 1:(length(index.bins)-1)){
      result.m<-rbind(result.m,index.fun.inner(outer.x[(index.bins[b]+1):index.bins[b+1],,drop=F]),p)
    }
    if(index.bins[length(index.bins)]!=nrow(outer.x)){
      result.m<-rbind(result.m,index.fun.inner(outer.x[(index.bins[length(index.bins)]+1):nrow(outer.x),,drop=F]),p)
    }
  }else{
    result.m<-index.fun.inner(outer.x,p)
  }
  return(result.m)
}