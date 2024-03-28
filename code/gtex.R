path="/home/bc758/research/decals/"
source(paste0(path,"code/DECALS.R"))

########################################################
################# GTEx Data Analysis ###################
########################################################
load(file=paste0(path,"data/gtex/bulk.RData"))
load(file=paste0(path,"data/gtex/sig.RData"))

########### CTS proportion estimation ##############
CTS_proportion=constraint_ols(sig=sig,bulk=bulk)
dim(CTS_proportion) # samples * CTS_proportions
cell_names=c("Ast","End","Mic","Ext","Inh","Oli")
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
CTS_proportion_var=propor_cov(W=t(sig),sigma_all=sigma_ols,propor_all=t(CTS_proportion))

seedN=1
library(mvtnorm)
set.seed(seedN)
CTS_proportion_seedN <- sapply(1:n,function(i) {
  frac1=rmvnorm(1,mean=CTS_proportion[i,],sigma=CTS_proportion_var[[i]])
  frac1[frac1<0]=0;frac1[frac1>1]=1
  frac1/sum(frac1)
} )
frac_one=t(frac_one)

####################################################
################ bMIND Analysis ####################
####################################################
load(file=paste0(path,"data/gtex/frac0.RData"))
load(file=paste0(path,"data/gtex/y.RData"))
library(MIND)
colnames(bulk) = rownames(frac0) = 1:nrow(frac0)
deconv = bmind_de(bulk1[1:5,], frac=frac0, y = y, np = T,ncore=1)
#deconv$pval
#we did this analysis for all genes and save the related p values.
#The following code shows how we obtain the DE gene sets and map them on chrosomes X and Y.

load(file=paste0(path,"data/gtex/datProbes.RData"))

chr_names=c(1:22,"X","Y")

chr_genes_table <- function(gene_sets){
  gene_chr_num=gene_chr_ratio=NULL
  for(chr in chr_names){
    gene_sets_chr=probes$external_gene_name[probes$chromosome_name==chr]
    gene_chr=match(gene_sets,gene_sets_chr)
    gene_chr_num=c(gene_chr_num,sum(!is.na(gene_chr)) )
    gene_chr_ratio=c(gene_chr_ratio,sum(!is.na(gene_chr))/length(gene_sets_chr))
  }
  rbind(gene_chr_num,gene_chr_ratio)
}

##################### map DE genes on chromosomes ######################
p_thr=0.05

load(file=paste0(path,"data/gtex/pval_all0.RData"))

CTS_chr=lapply(1:6,function(k)array(0,dim=c(3,24,100)))
gene_sets_sel=vector("list")
for(k in 1:6){
  pval_adjust=p.adjust(pval_all[,k],method="BH")
  gene_sets=rownames(pval_all)[pval_adjust<p_thr]
  gene_sets_sel[[k]]=gene_sets
  result=chr_genes_table(gene_sets)
  ratio2=result[1,]/length(gene_sets)
  result=rbind(result,ratio2)
  colnames(result)=chr_names
  CTS_chr[[k]]=result
}

gene_num0=sapply(1:6,function(k)length(gene_sets_sel[[k]]))


CTS_chr1=lapply(1:6,function(k)array(0,dim=c(3,24,100)))
gene_sets_all=vector("list")
load(file=paste0(path,"data/gtex/pval_all",1,".RData"))
for(k in 1:6){
  pval_adjust=p.adjust(pval_all[,k],method="BH")
  gene_sets=rownames(pval_all)[pval_adjust<p_thr]
  gene_sets_all[[k]]=gene_sets
}

for(i in 2:100){
  load(file=paste0(path,"data/gtex/pval_all",i,".RData"))
  for(k in 1:6){
    pval_adjust=p.adjust(pval_all[,k],method="BH")
    gene_sets=rownames(pval_all)[pval_adjust<p_thr]
    gene_sets_all[[k]]=c(gene_sets_all[[k]],gene_sets)
  }
}

gene_sets_sel1=gene_sets_all
for(k in 1:6){
  gene_hist=table(gene_sets_all[[k]])
  gene_sel1_id=order(gene_hist,decreasing=TRUE)
  thr=gene_hist[gene_sel1_id[gene_num0[k]]]
  gene_sets=names(gene_hist)[which(gene_hist>=thr)]
  gene_sets_sel1[[k]]=gene_sets
  result=chr_genes_table(gene_sets)
  ratio2=result[1,]/length(gene_sets)
  if(length(gene_sets)==0) ratio2=rep(0,24)
  result=rbind(result,ratio2)
  CTS_chr1[[k]]=result
}


x_data=matrix(0,2,6)
for(k in 1:6){
  x_data[,k]=c(CTS_chr[[k]][3,23],CTS_chr1[[k]][3,23])
}
rownames(x_data) = c("bMIND","bMIND+DECALS")
colnames(x_data)=cell_names
y_data=matrix(0,2,6)
for(k in 1:6){
  y_data[,k]=c(CTS_chr[[k]][3,24],CTS_chr1[[k]][3,24])
}
rownames(y_data) = c("bMIND","bMIND+DECALS")
colnames(y_data)=cell_names
x_genes=x_data
y_genes=y_data


load(file=paste0(path,"data/gtex/genes_x.RData"))
load(file=paste0(path,"data/gtex/genes_y.RData"))
load(file=paste0(path,"data/gtex/pval_all0.RData"))
xy_over <- function(genes_sets,genes_x,genes_y){
  xy_info=NULL
  for(k in 1:6){
    #pval_adjust=p.adjust(pval_all[,k],method="BH")
    gene_sets=genes_sets[[k]]
    gene_sets_xmale=sum(duplicated(c(names(genes_x$male_gene),gene_sets)))
    gene_sets_xfemale=sum(duplicated(c(names(genes_x$female_gene),gene_sets)))
    gene_sets_ymale=sum(duplicated(c(names(genes_y$male_gene),gene_sets)))
    gene_sets_yfemale=sum(duplicated(c(names(genes_y$female_gene),gene_sets)))
    xy_info=rbind(xy_info,c(gene_sets_xmale,gene_sets_xfemale,gene_sets_ymale,gene_sets_yfemale))
  }
  xy_info
}

tab0=xy_over(genes_sets=gene_sets_sel,genes_x,genes_y)
tab0_x=tab0[,2]/rowSums(tab0[,1:2])
tab0_y=tab0[,3]/rowSums(tab0[,3:4])

tab1=xy_over(genes_sets=gene_sets_sel1,genes_x,genes_y)
tab1_x=tab1[,2]/rowSums(tab1[,1:2])
tab1_y=tab1[,3]/rowSums(tab1[,3:4])

x_data=rbind(tab0_x,tab1_x)
rownames(x_data) = c("bMIND","bMIND+DECALS")
colnames(x_data)=cell_names

y_data=rbind(tab0_y,tab1_y)
rownames(y_data) = c("bMIND","bMIND+DECALS")
colnames(y_data)=cell_names
x_sex=x_data
y_sex=y_data

################# add 2 sd interval ################
CTS_chr2=lapply(1:6,function(k)array(0,dim=c(3,24,100)))
gene_sets_all=vector("list")
load(file=paste0(path,"data/gtex/pval_all",1,".RData"))
for(k in 1:6){
  pval_adjust=p.adjust(pval_all[,k],method="BH")
  gene_sets=rownames(pval_all)[pval_adjust<p_thr]
  gene_sets_sel[[k]]=gene_sets
  result=chr_genes_table(gene_sets)
  ratio2=result[1,]/length(gene_sets)
  result=rbind(result,ratio2)
  colnames(result)=chr_names
  CTS_chr2[[k]][,,1]=result
}

for(i in 2:100){
  load(file=paste0(path,"data/gtex/pval_all",i,".RData"))
  for(k in 1:6){
    pval_adjust=p.adjust(pval_all[,k],method="BH")
    gene_sets=rownames(pval_all)[pval_adjust<p_thr]
    gene_sets_sel[[k]]=gene_sets
    result=chr_genes_table(gene_sets)
    ratio2=result[1,]/length(gene_sets)
    result=rbind(result,ratio2)
    colnames(result)=chr_names
    CTS_chr2[[k]][,,i]=result
  }
}

CTS_chr_sd=array(0,dim=c(6,2))
for(k in 1:6){
  CTS_chr_sd[k,1]=sd(CTS_chr2[[k]][3,23,!is.nan(CTS_chr2[[k]][3,23,])])
  CTS_chr_sd[k,2]=sd(CTS_chr2[[k]][3,24,!is.nan(CTS_chr2[[k]][3,24,])])
}
CTS_chr_sd=CTS_chr_sd/10




tab1_x=tab1_y=matrix(0,100,6)
for(i in 1:100){
  load(file=paste0(path,"data/gtex/pval_all",i,".RData"))
  gene_sets_sel1=vector("list")
  for(k in 1:6){
    pval_adjust=p.adjust(pval_all[,k],method="BH")
    gene_sets=rownames(pval_all)[pval_adjust<p_thr]
    gene_sets_sel1[[k]]=gene_sets
  }
  
  tab1=xy_over(genes_sets=gene_sets_sel1,genes_x,genes_y)
  tab1_x[i,]=tab1[,2]/rowSums(tab1[,1:2])
  tab1_y[i,]=tab1[,3]/rowSums(tab1[,3:4])
}

sex_de_sd=matrix(0,6,2)
for(k in 1:6){
  sex_de_sd[k,1]=sd(tab1_x[!is.nan(tab1_x[,k]),k])
  sex_de_sd[k,2]=sd(tab1_y[!is.nan(tab1_y[,k]),k])
}
sex_de_sd=sex_de_sd/10





############### plot ######################
par(mfrow=c(2,2))
barplot(x_genes, beside=TRUE,
        main = "Proportions of DE transcripts on Chromosome X",
        #xlab = "Class",
        col = c("lightgreen","lightblue"),ylim=c(0,0.18)
)
for(k in 1:6){
  lines(x=c((3*k-1)+0.3,(3*k-1)+0.7),y=rep(x_genes[2,k]-2*CTS_chr_sd[k,1],2),lty=1,col="red")
  lines(x=c((3*k-1)+0.3,(3*k-1)+0.7),y=rep(x_genes[2,k]+2*CTS_chr_sd[k,1],2),lty=1,col="red")
  lines(x=c((3*k-1)+0.5,(3*k-1)+0.5),y=c(x_genes[2,k]-2*CTS_chr_sd[k,1],x_genes[2,k]+2*CTS_chr_sd[k,1]),lty=2,col="red")
}

barplot(y_genes, beside=TRUE,
        main = "Ratio of DE genes in Chromosome Y",
        #xlab = "Class",
        col = c("lightgreen","lightblue"),ylim=c(0,0.55)
)
for(k in 1:6){
  lines(x=c((3*k-1)+0.3,(3*k-1)+0.7),y=rep(y_genes[2,k]-2*CTS_chr_sd[k,2],2),lty=1,col="red")
  lines(x=c((3*k-1)+0.3,(3*k-1)+0.7),y=rep(y_genes[2,k]+2*CTS_chr_sd[k,2],2),lty=1,col="red")
  lines(x=c((3*k-1)+0.5,(3*k-1)+0.5),y=c(y_genes[2,k]-2*CTS_chr_sd[k,2],y_genes[2,k]+2*CTS_chr_sd[k,2]),lty=2,col="red")
}


barplot(x_sex, beside=TRUE,
        #main = "Chromosome X: the ratio of overexpressed DE genes in female",
        #xlab = "Class",
        col = c("lightgreen","lightblue"),ylim=c(0,1.02)
)
for(k in 1:6){
  lines(x=c((3*k-1)+0.3,(3*k-1)+0.7),y=rep(x_sex[2,k]-2*sex_de_sd[k,1],2),lty=1,col="red")
  lines(x=c((3*k-1)+0.3,(3*k-1)+0.7),y=rep(x_sex[2,k]+2*sex_de_sd[k,1],2),lty=1,col="red")
  lines(x=c((3*k-1)+0.5,(3*k-1)+0.5),y=c(x_sex[2,k]-2*sex_de_sd[k,1],x_sex[2,k]+2*sex_de_sd[k,1]),lty=2,col="red")
}
barplot(y_sex, beside=TRUE,
        #main = "Chromosome Y: the ratio of overexpressed DE genes in male",
        col = c("lightgreen","lightblue")
)
for(k in 1:6){
  lines(x=c((3*k-1)+0.3,(3*k-1)+0.7),y=rep(y_sex[2,k]-2*sex_de_sd[k,2],2),lty=1,col="red")
  lines(x=c((3*k-1)+0.3,(3*k-1)+0.7),y=rep(y_sex[2,k]+2*sex_de_sd[k,2],2),lty=1,col="red")
  lines(x=c((3*k-1)+0.5,(3*k-1)+0.5),y=c(y_sex[2,k]-2*sex_de_sd[k,2],y_sex[2,k]+2*sex_de_sd[k,2]),lty=2,col="red")
}
par(mfrow=c(1,1))




############### CTS proportions #####################
load(file=paste0(path,"data/gtex/y.RData"))
load(file=paste0(path,"data/gtex/frac0.RData"))
y[y==1]="Female"
y[y==0]="Male"
library(vioplot)
par(mfrow=c(1,6))
for(k in 1:6){
  vioplot(frac0[,k]~y,col=c("lightcoral","lightgreen"),xlab=NULL,ylab=NULL,
          main=cell_names[k],ylim=c(0,0.8),outline=FALSE)
}
par(mfrow=c(1,1))

bulk_est=t(frac0%*%t(sig))
sample_cor=sapply(1:1671,function(i)cor(bulk[,i],bulk_est[,i]))
boxplot(sample_cor,outline=FALSE,main="Correlations in GTEx Data")

bulk_log=log(bulk+1)
bulk_est_log=log(bulk_est+1)
par(mfrow=c(3,3))
for(i in 1:9){
  plot(bulk_log[,i],bulk_est_log[,i],xlab="Observed",ylab="Estiamted",main=colnames(bulk_log)[i])
  lim_val=min(c(max(bulk_log[,i])),max(bulk_est_log[,i]))
  lines(x=c(0,lim_val),y=c(0,lim_val),col=2,lty=2)
  legend("topleft", bty="n", text.col=2, legend=paste("cor=",round(cor(bulk[,i],bulk_est[,i]),2)))
}
par(mfrow=c(1,1))
