#########################################################
##################### mse ###############################
#########################################################
path="/home/bc758/research/decals/"
load(file=paste0(path,"data/rosmap/sig.RData"))
load(file=paste0(path,"data/rosmap/bulk.RData"))
#load(file=paste0(path,"data/rosmap/frac0.RData"))
load(file=paste0(path,"data/sc/bisque/sc.ref.RData"))
load(file=paste0(path,"data/sc/bisque/res.RData"))
load(file=paste0(path,"data/sc/bisque/estFrac_RNA_5CT.Robj"))
frac0=estFrac_RNA_5CT

match_id=match(colnames(res$bulk.props),rownames(frac0))
match_id1=match(colnames(res$sc.props),rownames(frac0))
match_id1_nna=which(!is.na(match_id1))
frac_bisque=matrix(0,541,5)
frac_bisque[match_id,]=t(res$bulk.props)
frac_bisque[match_id1[match_id1_nna],]=t(res$sc.props[,match_id1_nna])
rownames(frac_bisque)=rownames(frac0)
colnames(frac_bisque)=rownames(res$bulk.props)
frac_bisque=frac_bisque[,c(4,5,1,3,2)]


bulk_est_decals=t(as.matrix(frac0)%*%t(sig))
bulk_est_bisque=bulk_est_bisque0=t(frac_bisque%*%t(sc.ref))
for(j in 1:159){
  lm0=lm(bulk[j,match_id]~res$transformed.bulk[j,])
  bulk_est_bisque[j,]=lm0$coefficients[1]+bulk_est_bisque0[j,]*lm0$coefficients[2]
}
mse_exp=matrix(0,2,159)
for(j in 1:159){
  mse_exp[1,j]=mean((bulk_est_decals[j,]-bulk[j,])^2)
  mse_exp[2,j]=mean((bulk_est_bisque[j,]-bulk[j,])^2)
}
rownames(mse_exp)=c("DECALS","Bisque")
log10(rowMeans(mse_exp))
# DECALS   Bisque 
# 4.577820 8.820123 

#########################################################
##################### IHC ###############################
#########################################################
cell_names=c("neu","oli","ast","mic","end")
IHC_data=vector("list")
for(k in 1:5){
  IHC_k=read.delim(file=paste0(path,"data/sc/bisque/IHC_",cell_names[k],".txt"), header=FALSE)
  IHC_data[[k]]=IHC_k
}

sample_id=as.character(IHC_data[[1]][1,])
for(k in 2:5){
  sample_id=intersect(sample_id,as.character(IHC_data[[k]][1,]))
}

sample_frac=matrix(0,49,5)
for(k in 1:5){
  id0=match(sample_id,IHC_data[[k]][1,])
  frac_use=as.numeric(IHC_data[[k]][2,id0])
  sample_frac[,k]=frac_use
}
sample_frac0=sample_frac/rowSums(sample_frac)
rownames(sample_frac0)=sample_id



ROSMAP_clinical <- read.csv(file=paste0(path,"data/sc/bisque/ROSMAP_clinical.csv"))
id0=match(sample_id,ROSMAP_clinical$projid)
rownames(sample_frac0)=ROSMAP_clinical$individualID[id0]


#load("~/research/decals/data/ROSMAP/estFrac_RNA_5CT.Robj")
frac0=estFrac_RNA_5CT

id_frac=match(rownames(sample_frac0),rownames(frac0))
decals_frac=frac0[id_frac,]
bisque_frac=frac_bisque[id_frac,]








mse=matrix(0,2,5)
for(k in 1:5){
  mse[1,k]=mean((sample_frac0[,k]-decals_frac[,k])^2)
  mse[2,k]=mean((sample_frac0[,k]-bisque_frac[,k])^2)
}

rownames(mse)=c("DECALS","Bisque")
colnames(mse)=cell_names

colnames(mse)=c("Neu","Oli","Ast","Mic","End")
barplot(mse[1:2,],beside = TRUE,col=c("lightcoral","lightgreen"),ylab="MSE")
legend(12,0.045, legend=rownames(scor[1:2,]),fill=c("lightcoral","lightgreen"))

