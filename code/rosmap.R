path="/home/bc758/research/decals/"
source(paste0(path,"code/DECALS.R"))

########################################################
############### ROSMAP Data Analysis ###################
########################################################
load(file=paste0(path,"data/rosmap/bulk.RData")) # bulk expression data
load(file=paste0(path,"data/rosmap/sig.RData")) # signature matrix

########### CTS proportion estimation ##############
CTS_proportion=constraint_ols(sig=sig,bulk=bulk)
dim(CTS_proportion) # samples * CTS_proportions
cell_names=c("Neu","Oli","Ast","Mic","End")
colnames(CTS_proportion)=cell_names
boxplot(CTS_proportion,outline=FALSE,ylab="cell type proportions",ylim=c(0,0.8))

########### CTS proportion variance estimation ##############
lambda_ols <- tune_select(y_all=bulk,W=t(sig),propor_est=t(CTS_proportion),lambda_set=seq(0,0.5,by=0.025),tune_method="sequential")
# the tuning set can be changed, e.g. lambda_set=seq(0,0.6,0.02), the values should be in [0,1] 
sigma_ols0<-sigma_upt0(propor_all=t(CTS_proportion),y_all=bulk,W=t(sig),lambda_all=lambda_ols[[1]])


CTS_proportion_var=propor_cov0(W=t(sig),sigma_all=sigma_ols,propor_all=t(CTS_proportion))
dim(CTS_proportion_var) # CTS_proportions_var * samples

########## Simulate 100 CTS proportions ##################
# this returns a list of covariance matrix
# each element is one K*K  covariance matrix for cts proportions
CTS_proportion_var=propor_cov(W=t(sig),sigma_all=sigma_ols,propor_all=t(CTS_proportion))

seedN=1
library(mvtnorm)
set.seed(seedN)
CTS_proportion_seedN <- sapply(1:n,function(i) { # generate cts proportion by estimated mean and covariance
  frac1=rmvnorm(1,mean=CTS_proportion[i,],sigma=CTS_proportion_var[[i]]) 
  frac1[frac1<0]=0;frac1[frac1>1]=1 # control them between 0 and 1
  frac1/sum(frac1)
} )
frac_one=t(frac_one)

####################################################
################ bMIND Analysis ####################
####################################################
load(file=paste0(path,"data/rosmap/frac0.RData")) # cts proportion
load(file=paste0(path,"data/rosmap/y.RData")) # group label
library(MIND)
colnames(bulk) = rownames(frac0) = 1:nrow(frac0)
deconv = bmind_de(bulk1[1:5,], frac=frac0, y = y, np = T,ncore=1)
#deconv$pval
#we did this analysis for all genes and save the related p values.
#The following code shows how we obtain the DE gene sets.

gene_over <- function(cut_p=0.05,cut_prop=0.1){ # select signature genes from 100 replication
  name_over=vector("list")
  freq_all <- list()
  p_ave=vector("list")
  for(k in 1:K){
    pval_allk=NULL
    p_sel_all=NULL
    for(i in 1:100){
      load(paste(path,"data/rosmap/pval_all",i,".RData",sep=""))
      # load(paste(path,"bmind/result/protein/pval_all_",i,".RData",sep=""))
      name0=rownames(pval_all)
      name_sel = name0[which(pval_all[,k]<cut_p)]
      p_sel_all=c(p_sel_all,pval_all[,k][which(pval_all[,k]<cut_p)])
      pval_allk=c(pval_allk,name_sel) # get the genes with p value smaller than cut-off value
    }
    name_over0=table(pval_allk)
    name_over[[k]]=names(name_over0)[which(name_over0>=cut_prop*100)] # select the gene over 10 times
    p_ave0=NULL
    for(i in 1:length(name_over[[k]])){
      id_gene=which(pval_allk==name_over[[k]][i])
      p_ave0=c(p_ave0,mean(p_sel_all[id_gene]))
    }
    p_ave[[k]]=p_ave0 # average p values for selected genes
    freq_all[[k]] <- as.vector(name_over0)[which(name_over0>=cut_prop*100)]
    # frequency of selected genes with p value < 0.05
  }
  return(list(name_over, freq_all,p_ave))
}



K=5

########## select de genes by bMIND without considering decals ##############
cut_p=0.05;cut_prop=0.05
load(paste(path,"data/rosmap/pval_all0.RData",sep=""))
# p value for testing the difference in two groups by bMIND
# number of genes * number of cell types

gene_names=rownames(pval_all)
cts=colnames(pval_all)
gene_sel1=vector("list")
gene_sel1_pval <- list()
for(k in 1:ncol(pval_all)){ # select DE genes for one specific cell type
  pval_allk=pval_all[,k]
  gene_sel1[[k]]=gene_names[which(pval_allk<cut_p)] # set the alpha value here
  gene_sel1_pval[[k]] <- pval_allk[which(pval_allk<cut_p)]
}     
summary(gene_sel1_pval[[1]])
names(gene_sel1) <- names(gene_sel1_pval) <- cts

n_sel=sapply(1:K,function(i)length(gene_sel1[[i]]))
n_sel # number of signature genes in each cell type

gene_sel2_list=gene_over(cut_p=cut_p,cut_prop=cut_prop) # select DE gene with DECALS
gene_sel2 <- gene_sel2_list[[1]]
gene_sel2_freq <- gene_sel2_list[[2]]
gene_sel2_pval=gene_sel2_list[[3]]
names(gene_sel2_pval)= names(gene_sel2_freq) <- names(gene_sel2) <- cts
n_sel2=sapply(1:K,function(i)length(gene_sel2[[i]]))

ordered_gene_sel1 <- list() # order the DE genes with frequency number
for(ct in cts){
  ordered_gene_sel1[[ct]] <- gene_sel1[[ct]][order(gene_sel1_pval[[ct]])]
}
ordered_gene_sel2 <- list()
freq_val=NULL
for(k in 1:5){ # order the frequency in each cell type
  id0=order(gene_sel2_freq[[k]], decreasing = T)
  gene_freq0=gene_sel2_freq[[k]][id0]
  freq_val=c(freq_val,gene_freq0[n_sel[k]])
  id1=max(which(gene_freq0==freq_val[k]))
  ordered_gene_sel2[[k]] <- gene_sel2[[k]][id0[1:(id1)]]#gene_sel2[[k]][id0]#
}
sapply(1:5,function(i)length(ordered_gene_sel2[[i]]))

gene_before=ordered_gene_sel1 # selected DE genes without DECALS
gene_after=ordered_gene_sel2 # selected DE genes with DECALS

# save the DE gene set by cell types
Neu=list(before=gene_before[[1]],after=gene_after[[1]])
Oli=list(before=gene_before[[2]],after=gene_after[[2]])
Ast=list(before=gene_before[[3]],after=gene_after[[3]])
Mic=list(before=gene_before[[4]],after=gene_after[[4]])
End=list(before=gene_before[[5]],after=gene_after[[5]])

AD_DEG1.1 = list(Neu=Neu,Oli=Oli,Ast=Ast,Mic=Mic,End=End)


#### We use these DE genes to do the IPA enrichment analysis to obtain the table.



############### CTS proportions #####################
load(file=paste0(path,"data/rosmap/y.RData")) # group label
load(file=paste0(path,"data/rosmap/frac0.RData")) # cts proportions
library(vioplot)
par(mfrow=c(1,5))
for(k in 1:5){ #vioplot for cts proportion in two groups
  vioplot(frac0[,k]~y,col=c("lightcoral","lightgreen"),xlab=NULL,ylab=NULL,
          main=cell_names[k],ylim=c(0,0.8),outline=FALSE)
}
par(mfrow=c(1,1))

bulk_est=t(frac0%*%t(sig))
sample_cor=sapply(1:541,function(i)cor(bulk[,i],bulk_est[,i]))
boxplot(sample_cor,main="Correlations in ROSMAP Data")

bulk_log=log(bulk+1)
bulk_est_log=log(bulk_est+1)
par(mfrow=c(3,3))
for(i in 1:9){ # compare the observed and estimated bulk expression 
  plot(bulk_log[,i],bulk_est_log[,i],xlab="Observed",ylab="Estiamted",main=colnames(bulk_log)[i])
  lim_val=min(c(max(bulk_log[,i])),max(bulk_est_log[,i]))
  lines(x=c(0,lim_val),y=c(0,lim_val),col=2,lty=2)
  legend("bottomright", bty="n", text.col=2, legend=paste("cor=",round(cor(bulk_log[,i],bulk_est_log[,i]),2)))
}
par(mfrow=c(1,1))
