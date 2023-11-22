library(data.table)

load("test.RData")

z=z.list[[2]]
p_snp=length(z)
ld=ld.list[[2]]

snp_names=paste0("snp",1:p_snp)
z=cbind(snp_names,z=z)

fwrite(x=z,file="../03_unit_R/01_test_z.csv",sep=",",col.names = T,row.names = F,quote=F)
fwrite(x=ld,file="../03_unit_R/01_test_ld.csv",sep=",",col.names = F,row.names = F,quote=F)



ld=fread("../02_R_src/APOE_locus_ld.txt.gz",data.table=F)
z=fread("../02_R_src/APOE_locus_sumstats.txt.gz",data.table=F)
z=z$Z
p_snp=length(z)
snp_names=paste0("snp",1:p_snp)
z=cbind(snp_names,z=z)

fwrite(x=z,file="../03_unit_R/02_APOE_z.csv",sep=",",col.names = T,row.names = F,quote=F)
fwrite(x=ld,file="../03_unit_R/02_APOE_ld.csv",sep=",",col.names = F,row.names = F,quote=F)





ld=fread("../02_R_src/matrix.csv",data.table=F)
z=fread("../02_R_src/zScores.txt",data.table=F)
z=z[,1]
p_snp=length(z)
snp_names=paste0("snp",1:p_snp)
z=cbind(snp_names,z=z)

fwrite(x=z,file="../03_unit_R/03_test_z.csv",sep=",",col.names = T,row.names = F,quote=F)
fwrite(x=ld,file="../03_unit_R/03_test_ld.csv",sep=",",col.names = F,row.names = F,quote=F)


library(RcppCNPy)
ld=npyLoad("../01_test/dc16_etl_finemapping_output_-1335409592772313897_ld.npy")
z_full=fread("../01_test/02_test_df.csv",data.table=F)
z=z_full[,"beta"]/z_full[,"sebeta"]

z[1064]=-z[1064]

p_snp=length(z)
snp_names=paste0("snp",1:p_snp)
z=cbind(snp_names,z=z)

fwrite(x=z,file="../03_unit_R/02_test_z.csv",sep=",",col.names = T,row.names = F,quote=F)
npySave(filename = "../03_unit_R/02_test_ld.npy",object = ld)

#fwrite(x=ld,file="../03_unit_R/02_test_ld.csv",sep=",",col.names = F,row.names = F,quote=F)





