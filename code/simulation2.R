path="/home/bc758/research/decals/"
source(paste0(path,"code/DECALS.R"))

########################################################
################### Simulation 2 #######################
########################################################
load(file=paste0(path,"data/sim2/sig.RData"))
load(file=paste0(path,"data/sim2/frac0.RData"))
load(file=paste0(path,"data/sim2/sigma_ols0.RData"))

sim_data2 = data_gen2(frac0,sig,sigma_ols0,seedN=1)

bulk=sim_data2$bulk
sig=sim_data2$sig
CTS_proportion=constraint_ols(sig=sig,bulk=bulk)
lambda_ols <- tune_select(y_all=bulk,W=t(sig),propor_est=t(CTS_proportion),lambda_set=seq(0,0.5,by=0.025),tune_method="sequential")
sigma_ols0<-sigma_upt0(propor_all=t(CTS_proportion),y_all=bulk,W=t(sig),lambda_all=lambda_ols[[1]])

CTS_proportion_var=propor_cov0(W=t(sig),sigma_all=sigma_ols,propor_all=t(CTS_proportion))
dim(CTS_proportion_var) # CTS_proportions_var * samples
decals_res=list(CTS_proportion=CTS_proportion,CTS_proportion_var=CTS_proportion_var)

ols_res=ols_result(bulk=bulk,sig=sig)

sigma1=cor_est(sigma_ols0[[1]] )
sigma2=cor_est(sigma_ols0[[2]] )
sigma3=cor_est(sigma_ols0[[3]] )
sigma4=cor_est(sigma_ols0[[4]] )
sigma5=cor_est(sigma_ols0[[5]] )

library(MEAD)
p=159
R01=matrix(0,p,p)
net0=sigma1+sigma2+sigma3+sigma4+sigma5
R01[which(net0!=0)]=1
mead_res=MEAD_est(y=bulk,X=sig,Vg=matrix(0,159,5),R01=R01)


# RNA-Sieve is achieved by python, the code and result is in the folder "sim2-rna-sieve". 

##################################################################

n=541
k=5
sd_ols_all=sd_decals_all=sd_mead_all=array(0,dim=c(k,n,100))
Pi_ols_all=Pi_decals_all=Pi_mead_all=array(0,dim=c(k,n,100))
for(i in 1:100){
  load(paste(path,"data/sim2/decals_res_",i,".RData",sep=""))
  Pi_decals_all[,,i]=t(decals_res$CTS_proportion)
  sd_decals_all[,,i]=sqrt(decals_res$CTS_proportion_var)
  
  load(paste(path,"data/sim2/ols_res_",i,".RData",sep=""))
  Pi_ols_all[,,i]=t(ols_res$CTS_proportion)
  sd_ols_all[,,i]=sqrt(ols_res$CTS_proportion_var)
  
  load(paste(path,"data/sim2/mead_res_",i,".RData",sep=""))
  Pi_mead_all[,,i]=mead_res$p_hat
  sd_mead_all[,,i]=mead_res$p_hat_se
}


Pi=t(frac0)


alpha=0.975
cp_decals=NULL
for(j in 1:541){
  int_lower=Pi_decals_all[,j,]-qnorm(alpha)*sd_decals_all[,j,]
  int_upper=Pi_decals_all[,j,]+qnorm(alpha)*sd_decals_all[,j,]
  
  tab1=sapply(1:5,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  })
  
  
  cp_decals=rbind(cp_decals,tab1)
}

cp_ols=NULL
for(j in 1:541){
  int_lower=Pi_ols_all[,j,]-qnorm(alpha)*sd_ols_all[,j,]
  int_upper=Pi_ols_all[,j,]+qnorm(alpha)*sd_ols_all[,j,]
  
  tab1=sapply(1:5,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  })
  
  
  cp_ols=rbind(cp_ols,tab1)
}


cp_mead=NULL
for(j in 1:541){
  int_lower=Pi_mead_all[,j,]-qnorm(alpha)*sd_mead_all[,j,]
  int_upper=Pi_mead_all[,j,]+qnorm(alpha)*sd_mead_all[,j,]
  
  tab1=sapply(1:5,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  })
  
  
  cp_mead=rbind(cp_mead,tab1)
}



ci_lower <- function(x){
  x=as.character(x)
  n_x=nchar(x)
  x1=substr(x,2,n_x-1)
  x2=strsplit(x1, ", ")[[1]][[1]]
  as.numeric(x2)
}

ci_upper <- function(x){
  x=as.character(x)
  n_x=nchar(x)
  x1=substr(x,2,n_x-1)
  x2=strsplit(x1, ", ")[[1]][[2]]
  as.numeric(x2)
}
id0=NULL
for(i in 1:100){
  id0=c(id0,file.exists(paste0(path,"data/sim2-rna-sieve/results/set2/ci95_",i,".csv")))
}

id1=which(id0)
ci_lower_all=ci_upper_all=array(0,dim=c(5,541,length(id1)))
for(i in 1:length(id1)){
  ci0=read.csv(paste0(path,"data/sim2-rna-sieve/results/set2/ci95_",id1[i],".csv"))
  ci=ci0[,c(5,6,2,4,3)]
  ci_1_lower=apply(ci,c(1,2),ci_lower)
  ci_lower_all[,,i]=t(ci_1_lower)
  ci_1_upper=apply(ci,c(1,2),ci_upper)
  ci_upper_all[,,i]=t(ci_1_upper)
}


cp_rna_sieve=NULL
for(j in 1:541){
  int_lower=ci_lower_all[,j,]
  int_upper=ci_upper_all[,j,]
  
  tab1=sapply(1:5,function(k){
    sum(sapply(1:length(id1),function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/length(id1)
  })
  cp_rna_sieve=rbind(cp_rna_sieve,tab1)
}
cp_rna_sieve[is.na(cp_rna_sieve)]=0



library(vioplot)
par(mfrow=c(2,3))
for(k in 1:5){
  result=cbind(cp_ols[,k],cp_rna_sieve[,k],cp_mead[,k],cp_decals[,k])
  colnames(result)=c("OLS","RNA-Sieve","MEAD","DECALS")
  vioplot(result,col=c("grey","lightgreen","lightcoral","skyblue"),#axes=FALSE,
          main=paste("Cell",k),ylim=c(0,1),outline=FALSE,
          ylab="coverage probability")
  lines(c(0.5,3.5),c(0.95,0.95),lty=2)
}
par(mfrow=c(1,1))


###############################################################
##################### Gaussian Setting ########################
###############################################################

load(file=paste0(path,"data/sim2/gaussian/sig.RData"))
load(file=paste0(path,"data/sim2/gaussian/frac0.RData"))
load(file=paste0(path,"data/sim2/gaussian/sigma_ols0.RData"))

sim_data3 = data_gen3(frac0,sig,sigma_ols0,seedN=1)

bulk=sim_data3$bulk
sig=sim_data3$sig
CTS_proportion=constraint_ols(sig=sig,bulk=bulk)
lambda_ols <- tune_select(y_all=bulk,W=t(sig),propor_est=t(CTS_proportion),lambda_set=seq(0,0.5,by=0.025),tune_method="sequential")
sigma_ols0<-sigma_upt0(propor_all=t(CTS_proportion),y_all=bulk,W=t(sig),lambda_all=lambda_ols[[1]])

CTS_proportion_var=propor_cov0(W=t(sig),sigma_all=sigma_ols,propor_all=t(CTS_proportion))
dim(CTS_proportion_var) # CTS_proportions_var * samples
decals_res=list(CTS_proportion=CTS_proportion,CTS_proportion_var=CTS_proportion_var)

ols_res=ols_result(bulk=bulk,sig=sig)

sigma1=cor_est(sigma_ols0[[1]] )
sigma2=cor_est(sigma_ols0[[2]] )
sigma3=cor_est(sigma_ols0[[3]] )
sigma4=cor_est(sigma_ols0[[4]] )
sigma5=cor_est(sigma_ols0[[5]] )

p=159
R01=matrix(0,p,p)
net0=sigma1+sigma2+sigma3+sigma4+sigma5
R01[which(net0!=0)]=1
mead_res=MEAD_est(y=bulk,X=sig,Vg=matrix(0,159,5),R01=R01)


# RNA-Sieve is achieved by python, the code and result is in the folder "sim2-rna-sieve". 


##################################################################

n=541
k=5
sd_ols_all=sd_decals_all=sd_mead_all=array(0,dim=c(k,n,100))
Pi_ols_all=Pi_decals_all=Pi_mead_all=array(0,dim=c(k,n,100))
for(i in 1:100){
  load(paste(path,"data/sim2/gaussian/decals_res_",i,".RData",sep=""))
  Pi_decals_all[,,i]=t(decals_res$CTS_proportion)
  sd_decals_all[,,i]=sqrt(decals_res$CTS_proportion_var)
  
  load(paste(path,"data/sim2/gaussian/ols_res_",i,".RData",sep=""))
  Pi_ols_all[,,i]=t(ols_res$CTS_proportion)
  sd_ols_all[,,i]=sqrt(ols_res$CTS_proportion_var)
  
  load(paste(path,"data/sim2/gaussian/mead_res_",i,".RData",sep=""))
  Pi_mead_all[,,i]=mead_res$p_hat
  sd_mead_all[,,i]=mead_res$p_hat_se
}


Pi=t(frac0)


alpha=0.975
cp_decals=NULL
for(j in 1:541){
  int_lower=Pi_decals_all[,j,]-qnorm(alpha)*sd_decals_all[,j,]
  int_upper=Pi_decals_all[,j,]+qnorm(alpha)*sd_decals_all[,j,]
  
  tab1=sapply(1:5,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  })
  
  
  cp_decals=rbind(cp_decals,tab1)
}

cp_ols=NULL
for(j in 1:541){
  int_lower=Pi_ols_all[,j,]-qnorm(alpha)*sd_ols_all[,j,]
  int_upper=Pi_ols_all[,j,]+qnorm(alpha)*sd_ols_all[,j,]
  
  tab1=sapply(1:5,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  })
  
  
  cp_ols=rbind(cp_ols,tab1)
}


cp_mead=NULL
for(j in 1:541){
  int_lower=Pi_mead_all[,j,]-qnorm(alpha)*sd_mead_all[,j,]
  int_upper=Pi_mead_all[,j,]+qnorm(alpha)*sd_mead_all[,j,]
  
  tab1=sapply(1:5,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  })
  
  
  cp_mead=rbind(cp_mead,tab1)
}



ci_lower <- function(x){
  x=as.character(x)
  n_x=nchar(x)
  x1=substr(x,2,n_x-1)
  x2=strsplit(x1, ", ")[[1]][[1]]
  as.numeric(x2)
}

ci_upper <- function(x){
  x=as.character(x)
  n_x=nchar(x)
  x1=substr(x,2,n_x-1)
  x2=strsplit(x1, ", ")[[1]][[2]]
  as.numeric(x2)
}
id0=NULL
for(i in 1:100){
  id0=c(id0,file.exists(paste0(path,"data/sim2-rna-sieve/results/set1/ci95_",i,".csv")))
}

id1=which(id0)
ci_lower_all=ci_upper_all=array(0,dim=c(5,541,length(id1)))
for(i in 1:length(id1)){
  ci0=read.csv(paste0(path,"data/sim2-rna-sieve/results/set1/ci95_",id1[i],".csv"))
  ci=ci0[,c(5,6,2,4,3)]
  ci_1_lower=apply(ci,c(1,2),ci_lower)
  ci_lower_all[,,i]=t(ci_1_lower)
  ci_1_upper=apply(ci,c(1,2),ci_upper)
  ci_upper_all[,,i]=t(ci_1_upper)
}


cp_rna_sieve=NULL
for(j in 1:541){
  int_lower=ci_lower_all[,j,]
  int_upper=ci_upper_all[,j,]
  
  tab1=sapply(1:5,function(k){
    sum(sapply(1:length(id1),function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/length(id1)
  })
  cp_rna_sieve=rbind(cp_rna_sieve,tab1)
}
cp_rna_sieve[is.na(cp_rna_sieve)]=0



library(vioplot)
par(mfrow=c(2,3))
for(k in 1:5){
  result=cbind(cp_ols[,k],cp_rna_sieve[,k],cp_mead[,k],cp_decals[,k])
  colnames(result)=c("OLS","RNA-Sieve","MEAD","DECALS")
  vioplot(result,col=c("grey","lightgreen","lightcoral","skyblue"),#axes=FALSE,
          main=paste("Cell",k),ylim=c(0,1),outline=FALSE,
          ylab="coverage probability")
  lines(c(0.5,3.5),c(0.95,0.95),lty=2)
}
par(mfrow=c(1,1))
