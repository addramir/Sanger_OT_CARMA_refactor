setwd("~/Projects/Sanger_OT_CARMA_refactor/03_unit_R/")
source("../02_R_src/00_CARMA_function.R")
library(data.table)
source("../01_test/CARMA_fixed_sigma.R")

library(CARMA)


#test1

z_full=fread("01_test_z.csv",header=T,sep=",",data.table=F)
ld=fread("01_test_ld.csv",header=F,sep=",",data.table=F)
z=z_full[,"z"]

R1=CARMA_spike_slab_noEM_rcpp(
        z,   
        ld,
        lambda=1,
        path_to_src="../02_R_src/",
        sparse_fun=FALSE)

R2=CARMA_spike_slab_noEM_rcpp(
  z,   
  ld,
  lambda=1,
  path_to_src="../02_R_src/",
  sparse_fun=TRUE)

#R3=CARMA_fixed_sigma(z.list = list(z),ld.list = list(ld),lambda.list = list(1),output.labels = NULL)
R3=CARMA(z.list = list(z),ld.list = list(ld),lambda.list = list(1),output.labels = NULL)
R3=R3[[1]]

R1$Outliers
R2$Outliers
R3$Outliers
cor(cbind(R1$PIPs,R2$PIPs,R3$PIPs))

#test2

z_full=fread("02_APOE_z.csv",header=T,sep=",",data.table=F)
ld=fread("02_APOE_ld.csv",header=F,sep=",",data.table=F)
z=z_full[,"z"]

R1=CARMA_spike_slab_noEM_rcpp(
  z,   
  ld,
  lambda=1,
  path_to_src="../02_R_src/",
  sparse_fun=FALSE)

R2=CARMA_spike_slab_noEM_rcpp(
  z,   
  ld,
  lambda=1,
  path_to_src="../02_R_src/",
  sparse_fun=TRUE)

#R3=CARMA_fixed_sigma(z.list = list(z),ld.list = list(ld),lambda.list = list(1),output.labels = NULL)
R3=CARMA(z.list = list(z),ld.list = list(ld),lambda.list = list(1),output.labels = NULL)
R3=R3[[1]]

R1$Outliers
R2$Outliers
R3$Outliers
cor(cbind(R1$PIPs,R2$PIPs,R3$PIPs))

out=array(NA,c(p_snp,5))
i=1
for (i in 1:5){
  R3=CARMA(z.list = list(z),ld.list = list(ld),lambda.list = list(2),output.labels = NULL)
  R3=R3[[1]]
  out[,i]=R3$PIPs
}

#test3

z_full=fread("03_test_z.csv",header=T,sep=",",data.table=F)
ld=fread("03_test_ld.csv",header=F,sep=",",data.table=F)
z=z_full[,"z"]

R1=CARMA_spike_slab_noEM_rcpp(
  z,   
  ld,
  lambda=1,
  path_to_src="../02_R_src/",
  sparse_fun=FALSE)

R2=CARMA_spike_slab_noEM_rcpp(
  z,   
  ld,
  lambda=1,
  path_to_src="../02_R_src/",
  sparse_fun=TRUE)

#R3=CARMA_fixed_sigma(z.list = list(z),ld.list = list(ld),lambda.list = list(1),output.labels = NULL)
R3=CARMA(z.list = list(z),ld.list = list(ld),lambda.list = list(1),output.labels = NULL)
R3=R3[[1]]

R1$Outliers
R2$Outliers
R3$Outliers
cor(cbind(R1$PIPs,R2$PIPs,R3$PIPs))

out=array(NA,c(length(z),5))
i=1
for (i in 1:5){
  R3=CARMA_fixed_sigma(z.list = list(z),ld.list = list(ld),lambda.list = list(2),output.labels = NULL)
  R3=R3[[1]]
  out[,i]=R3$PIPs
}
cor(out)


out=array(NA,c(p_snp,5))
i=1
for (i in 1:5){
  R1=CARMA_spike_slab_noEM_rcpp(
    z,   
    ld,
    lambda=1,
    path_to_src="../02_R_src/",
    sparse_fun=FALSE)
  out[,i]=R1$PIPs
}

x=R1$all.C.list[[1]][[2]]
l=apply(x,1,FUN=paste0,collapse="")
table(duplicated(l))

x=R3$all.C.list[[1]][[2]]
l=apply(x,1,FUN=paste0,collapse="")
table(duplicated(l))

