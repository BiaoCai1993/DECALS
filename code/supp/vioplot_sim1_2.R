########################################################
################### Simulation 1 #######################
########################################################
path="/home/bc758/research/decals/result/sim/"
n=500
k=3
sd_decals_all=array(0,dim=c(k,n,100))
Pi_decals_all=array(0,dim=c(k,n,100))
for(i in 1:100){
  load(paste(path,"sim1/decals_res_",i,".RData",sep=""))
  Pi_decals_all[,,i]=t(decals_res$CTS_proportion)
  sd_decals_all[,,i]=sqrt(decals_res$CTS_proportion_var)
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
colnames(cp_decals)=c("Cell 1", "Cell 2", "Cell 3")

library(vioplot)
vioplot(cp_decals, col="skyblue",ylim=c(0.6,1),ylab="coverage probability",main="DECALS")
lines(c(0.5,3.5),c(0.95,0.95),lty=2)


########################################################
################### Simulation 2 #######################
########################################################
load(file="/home/bc758/research/decals/code/jasa/data/sim2/frac0.RData")

n=541
k=5
sd_decals_all=array(0,dim=c(k,n,100))
Pi_decals_all=array(0,dim=c(k,n,100))
for(i in 1:100){
  load(paste(path,"data/sim2/decals_res_",i,".RData",sep=""))
  Pi_decals_all[,,i]=t(decals_res$CTS_proportion)
  sd_decals_all[,,i]=sqrt(decals_res$CTS_proportion_var)
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
colnames(cp_decals)=c("Cell 1", "Cell 2", "Cell 3", "Cell 4", "Cell 5")


library(vioplot)
vioplot(cp_decals, col="skyblue",ylim=c(0,1),ylab="coverage probability",main="DECALS")
lines(c(0.5,5.5),c(0.95,0.95),lty=2)
