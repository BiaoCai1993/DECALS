path="/home/bc758/research/decals/"
source(paste0(path,"code/DECALS.R"))

########################################################
################### Simulation 1 #######################
########################################################
seedN=1
sim_data1=data_gen1(n=500,p=300,k=3,seedN=seedN)
bulk=sim_data1$bulk
sig=sim_data1$sig
CTS_proportion=constraint_ols(sig=sig,bulk=bulk)
lambda_ols <- tune_select(y_all=bulk,W=t(sig),propor_est=t(CTS_proportion),lambda_set=seq(0,0.5,by=0.025),tune_method="sequential")
sigma_ols0<-sigma_upt0(propor_all=t(CTS_proportion),y_all=bulk,W=t(sig),lambda_all=lambda_ols[[1]])

CTS_proportion_var=propor_cov0(W=t(sig),sigma_all=sigma_ols,propor_all=t(CTS_proportion))
dim(CTS_proportion_var) # CTS_proportions_var * samples
decals_res=list(CTS_proportion=CTS_proportion,CTS_proportion_var=CTS_proportion_var)

ols_res=ols_result(bulk=bulk,sig=sig)


p=300
id1=1:(p/3);id2=(p/3+1):(2*p/3);id3=(2*p/3+1):p
sigma0_1=sigma_gen(p=100,rho=0.7,type="AR")
sigma0_2=matrix(0.3,p/3,p/3);diag(sigma0_2)=1


sigma1=sigma2=sigma3=matrix(0,p,p)
sigma1[id1,id1]=sigma0_2;sigma1[id2,id2]=sigma1[id3,id3]=sigma0_1
sigma2[id2,id2]=sigma0_2;sigma2[id1,id1]=sigma2[id3,id3]=sigma0_1
sigma3[id3,id3]=sigma0_2;sigma3[id1,id1]=sigma3[id2,id2]=sigma0_1


library(MEAD)
R01=matrix(0,p,p)
net0=sigma1+sigma2+sigma3
R01[which(net0!=0)]=1
mead_res=MEAD_est(y=bulk,X=sig,Vg=matrix(0,300,3),R01=R01)


#repeat this 100 times with seedN from 1 to 100

n=500
k=3
sd_ols_all=sd_decals_all=sd_mead_all=array(0,dim=c(k,n,100))
Pi_ols_all=Pi_decals_all=Pi_mead_all=array(0,dim=c(k,n,100))
for(i in 1:100){
  load(paste(path,"data/sim1/decals_res_",i,".RData",sep=""))
  Pi_decals_all[,,i]=t(decals_res$CTS_proportion)
  sd_decals_all[,,i]=sqrt(decals_res$CTS_proportion_var)
  
  load(paste(path,"data/sim1/ols_res_",i,".RData",sep=""))
  Pi_ols_all[,,i]=t(ols_res$CTS_proportion)
  sd_ols_all[,,i]=sqrt(ols_res$CTS_proportion_var)
  
  load(paste(path,"data/sim1/mead_res_",i,".RData",sep=""))
  Pi_mead_all[,,i]=mead_res$p_hat
  sd_mead_all[,,i]=mead_res$p_hat_se
}


library(mvtnorm)
p=300
id1=1:(p/3);id2=(p/3+1):(2*p/3);id3=(2*p/3+1):p

set.seed(1)
W<-matrix(rnorm(k*p,0,1),k,p)
Pi<-t(rdirichlet(n,c(3,2,1)))


alpha=0.975
cp_decals=NULL
for(j in 1:500){
  int_lower=Pi_decals_all[,j,]-qnorm(alpha)*sd_decals_all[,j,]
  int_upper=Pi_decals_all[,j,]+qnorm(alpha)*sd_decals_all[,j,]
  
  tab1=sapply(1:3,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  })
  
  
  cp_decals=rbind(cp_decals,tab1)
}

cp_ols=NULL
for(j in 1:500){
  int_lower=Pi_ols_all[,j,]-qnorm(alpha)*sd_ols_all[,j,]
  int_upper=Pi_ols_all[,j,]+qnorm(alpha)*sd_ols_all[,j,]
  
  tab1=sapply(1:3,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  })
  
  
  cp_ols=rbind(cp_ols,tab1)
}


cp_mead=NULL
for(j in 1:500){
  int_lower=Pi_mead_all[,j,]-qnorm(alpha)*sd_mead_all[,j,]
  int_upper=Pi_mead_all[,j,]+qnorm(alpha)*sd_mead_all[,j,]
  
  tab1=sapply(1:3,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  })
  
  
  cp_mead=rbind(cp_mead,tab1)
}

library(vioplot)
par(mfrow=c(1,3))
for(k in 1:3){
  result=cbind(cp_ols[,k],cp_mead[,k],cp_decals[,k])
  colnames(result)=c("OLS","MEAD","DECALS")
  vioplot(result,col=c("grey","lightcoral","skyblue"),#axes=FALSE,
          main=paste("Cell",k),ylim=c(0.6,1),outline=FALSE,
          ylab="coverage probability")
  lines(c(0.5,3.5),c(0.95,0.95),lty=2)
}
par(mfrow=c(1,1))

##################### GLS #############################
source(paste0(path,"code/jointalg.R"))
y_all=bulk
W=t(sig)
joint_est=joint_model(W,y_all,lambda_all)
CTS_proportion=joint_est$pi_est
CTS_proportion_var=propor_cov0(W=t(sig),sigma_all=joint_est$sigma_est,propor_all=t(CTS_proportion))

load(file=paste0(path,"data/sim1/gls_cp.RData"))
par(mfrow=c(1,3))
for(k in 1:3){
  result=gls_cp[[k]]
  colnames(result)=c("OLS true","GLS true","OLS est","GLS est")
  vioplot(result,col=c("mistyrose","skyblue"),axes=FALSE,
          main=paste("Cell",k),ylim=c(0.3,1),#outline=FALSE,
          ylab="coverage probability")
  lines(c(0.5,4.5),c(0.95,0.95),lty=2)
}
par(mfrow=c(1,1))

##################### Noise #############################
seedN=1
sim_data1=data_gen1_noise(n=500,p=300,k=3,seedN=seedN,noise=0.1)
#a = noise, it can be specified as 0.1,0.2,...,1
bulk=sim_data1$bulk
sig=sim_data1$sig
CTS_proportion=constraint_ols(sig=sig,bulk=bulk)
lambda_ols <- tune_select(y_all=bulk,W=t(sig),propor_est=t(CTS_proportion),lambda_set=seq(0,0.5,by=0.025),tune_method="sequential")
sigma_ols0<-sigma_upt0(propor_all=t(CTS_proportion),y_all=bulk,W=t(sig),lambda_all=lambda_ols[[1]])

CTS_proportion_var=propor_cov0(W=t(sig),sigma_all=sigma_ols,propor_all=t(CTS_proportion))
dim(CTS_proportion_var) # CTS_proportions_var * samples
decals_res=list(CTS_proportion=CTS_proportion,CTS_proportion_var=CTS_proportion_var)



n=500
k=3
sd_decals_all=array(0,dim=c(k,n,100))
Pi_decals_all=array(0,dim=c(k,n,100))
for(i in 1:100){
  load(paste(path,"data/sim1/noise1/decals_res_",i,".RData",sep=""))
  Pi_decals_all[,,i]=t(decals_res$CTS_proportion)
  sd_decals_all[,,i]=sqrt(decals_res$CTS_proportion_var)
}
load(paste0(path,"data/sim1/noise1/Pi.RData"))

alpha=0.975
cp_decals=NULL
for(j in 1:500){
  int_lower=Pi_decals_all[,j,]-qnorm(alpha)*sd_decals_all[,j,]
  int_upper=Pi_decals_all[,j,]+qnorm(alpha)*sd_decals_all[,j,]
  
  tab1=sapply(1:3,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  })
  
  
  cp_decals=rbind(cp_decals,tab1)
}


cell1=cell2=cell3=NULL
for(j in 1:10){
  load(file=paste0(path,"data/sim1/noise",j,"/cp_decals.RData"))
  cell1=cbind(cell1,cp_decals[,1])
  cell2=cbind(cell2,cp_decals[,2])
  cell3=cbind(cell3,cp_decals[,3])
}
colnames(cell1)=colnames(cell2)=colnames(cell3)=seq(0.1,1,by=0.1)
noise_cp=list(cell1=cell1,cell2=cell2,cell3=cell3)
par(mfrow=c(3,1))
for(k in 1:3){
  result=noise_cp[[k]]
  vioplot(result,col=c("skyblue"),#axes=FALSE,
          main=paste("Cell",k),ylim=c(0,1),outline=FALSE,
          ylab="coverage probability")
  lines(c(0.5,10.5),c(0.95,0.95),lty=2)
}
par(mfrow=c(1,1))


########################### heatplot ##########################
require(latex2exp)
require(pheatmap)
require(RColorBrewer)



# color for heatmap
palPos <- colorRampPalette(c("white", "red"), space = "rgb")
palNeg <- colorRampPalette(c("blue", "white"), space = "rgb")
coul <- c(palNeg(20), palPos(20))


pheatmap(sigma1, color = coul,#rev(coul(20))
         cluster_rows = F, cluster_cols = F,
         breaks = seq(-1, 1, by = 0.05),
         na_col = 'grey',
         show_rownames = F, show_colnames = F)

pheatmap(sigma2, color = coul,#rev(coul(20))
         cluster_rows = F, cluster_cols = F,
         breaks = seq(-1, 1, by = 0.05),
         na_col = 'grey',
         show_rownames = F, show_colnames = F)

pheatmap(sigma3, color = coul,#rev(coul(20))
         cluster_rows = F, cluster_cols = F,
         breaks = seq(-1, 1, by = 0.05),
         na_col = 'grey',
         show_rownames = F, show_colnames = F)