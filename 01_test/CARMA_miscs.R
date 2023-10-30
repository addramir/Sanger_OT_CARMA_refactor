###The computation of the credible set based on the final results of the fine-mapping step.
credible.set.fun.improved<-function(pip,ld,true.beta=NULL,rho=0.99){
  candidate.r<-c((seq(from=.5,0.95,by=0.05)),seq(0.96,0.99,0.01))
  
  snp.list<-list()
  colnames(ld)<-rownames(ld)<-1:nrow(ld)
  for(r in 1:length(candidate.r)){
    working.ld<-ld
    cor.threshold<-candidate.r[r]
    pip.order<-order(pip,decreasing = T)
    snp.list[[r]]<-list()
    group.index<-1
    s<-1
    while(sum(pip[pip.order[s:length(pip.order)]])>rho){
      
      cor.group<-as.numeric(names(which(abs(working.ld[which(pip.order[s]==as.numeric(colnames(working.ld))),])>cor.threshold)))
      if(sum(pip[cor.group])>rho){
        group.pip<- pip[cor.group]
        snp.index<- cor.group[order(group.pip,decreasing = T)][1: min(which(cumsum( sort(group.pip,decreasing = T))>rho))]
        snp.list[[r]][[group.index]]<-snp.index
        group.index<-group.index+1
        pip.order<-pip.order[-match(snp.index,pip.order)]
        working.ld<-working.ld[-match(snp.index,as.numeric(colnames(working.ld))),-match(snp.index,as.numeric(colnames(working.ld)))]
      }else{
        s=s+1
      }
    }
    
  }
  if(sum(sapply(snp.list,length))!=0){
    group.index<-max(which(sapply(snp.list,length)==max(sapply(snp.list,length))))
    credible.set.list<-snp.list[[group.index]]
    if(!is.null(true.beta)){
      purity<-c()
      for(s in 1:length(credible.set.list)){
        purity<-c(purity, (mean(ld[ credible.set.list[[s]], credible.set.list[[s]]]^2)))
      }
      causal.snp<-ifelse(length(na.omit(match(unlist(credible.set.list),true.beta)))!=0,
                         length(na.omit(match(unlist(credible.set.list),true.beta))),0)
      length.credible<-length(credible.set.list)
      return(list(c(causal.snp,
                    ifelse(length.credible==0,NA,length.credible),
                    mean(sapply(credible.set.list,length)),
                    mean(purity)),credible.set.list))
      
    }else{
      purity<-c()
      for(s in 1:length(credible.set.list)){
        purity<-c(purity, (mean(ld[ credible.set.list[[s]], credible.set.list[[s]]]^2)))
      }
      length.credible<-length(credible.set.list)
      return(list(c(ifelse(length.credible==0,NA,length.credible),
                    mean(sapply(credible.set.list,length)),
                    mean(purity)),
                  credible.set.list))
    }
  }else{
    return(list(rep(0,4),list()))
  }
}
###The computation of the credible models based on the final results of the fine-mapping step.
credible.model.fun<-function(likelihood,model.space,bayes.threshold=10){
  na.index<-which(is.na(likelihood))
  if(length(na.index)!=0){
    likelihood<-likelihood[-na.index]
    model.space<-model.space[-na.index,]
  }
  post.like.temp<-likelihood-likelihood[1]
  post.prob<-exp(post.like.temp)/(sum(exp(post.like.temp)))
  bayes.f<-post.prob[1]/post.prob
  candidate.model<-1
  
  credible.model.list<-list()
  credible.model.list[[1]]<-list()
  input.rs<-c()
  for(ss in 1:length(which(bayes.f<bayes.threshold))){
    credible.model.list[[1]][[ss]]<-which(model.space[ss,]==1)
    input.rs<-c(input.rs,which(model.space[ss,]==1))
  }
  credible.model.list[[2]]<-data.frame(Posterior.Prob=post.prob[which(bayes.f<bayes.threshold)])
  credible.model.list[[3]]<-unique(input.rs)
  return(credible.model.list )
  
  
}


