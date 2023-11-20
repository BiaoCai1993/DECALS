##############################################################
################################ WLS #########################
##############################################################
args = (commandArgs(TRUE))
eval(parse(text = args[[1]]))
cat("Iteration Number = ", seedN)

path="/home/bc758/research/decals/"
source(paste(path,"code/wls.R",sep=""))
load(file=paste0(path,"data/sim2/frac0.RData"))
library(mvtnorm)
library(gtools)
library(quadprog)

constraint_ols<-function(sig, bulk,wt){ 
  p=dim(sig)[1];K=dim(sig)[2];n=dim(bulk)[2]
  A=cbind(rep(1,K),diag(K))
  b0=c(1,rep(0,K))
  
  
  frac0=NULL
  for(i in 1:n){
    D=t(sig)%*%diag(wt[i,])%*%sig/n
    d=bulk[,i]%*%diag(wt[i,])%*%sig/n
    max_lambda=max(eigen(D)$values)
    frac0=rbind(frac0,solve.QP(Dmat=D/(max_lambda),dvec=d/(max_lambda),Amat=A,bvec=b0,meq=1)$solution)
  }
  frac0[frac0<0]=0
  return(frac0)
}



p=159
n<-541
k<-5

load(file=paste0(path,"data/sim2/res_all_",seedN,".RData"))

y_all=res_all$y_all
W=res_all$W
Pi=res_all$Pi
wt=matrix(1,n,p)
var_cell=W^2*100
wt=(frac0^2%*%var_cell)
wt=1/log(wt+1)
up_val=quantile(wt,0.99)
low_val=quantile(wt,0.01)
wt[wt>up_val]=up_val
wt[wt<low_val]=low_val
#wt[wt<0.1]=0.1
pi_wls_bf=matrix(1/5,n,k)

pi_wls=constraint_ols(sig=t(W),bulk=y_all,wt=wt)

lambda_ols <- tune_select(y_all=y_all,W=W,wt=wt,propor_est=t(pi_wls),lambda_set=seq(0,0.5,by=0.025),tune_method="sequential")
sigma_ols<-sigma_upt0(propor_all=t(pi_wls),y_all=y_all,W=W,wt=wt,lambda_all=lambda_ols[[1]])

pi_wls_var=propor_cov0(W=W,sigma_all=sigma_ols,propor_all=t(pi_wls),wt=wt)
res_wls=list(pi_wls=pi_wls,pi_wls_var=pi_wls_var)
#save(res_wls,file=paste0("/home/bc758/research/decals/result/sim/wls/8/res_wls_",seedN,".RData"))


#################### Data Analysis ###################
n=541
k=5
sd_wls_all=array(0,dim=c(k,n,100))
Pi_wls_all=array(0,dim=c(k,n,100))
for(i in 1:100){
  load(paste(path,"data/sim2/res_wls_",i,".RData",sep=""))
  Pi_wls_all[,,i]=t(res_wls$pi_wls)
  var0=res_wls$pi_wls_var
  var0[var0<0]=0
  sd_wls_all[,,i]=sqrt(var0)
}

load(file=paste0(path,"data/sim2/res_all_",seedN,".RData"))
Pi=res_all$Pi


alpha=0.975
cp_wls=NULL
for(j in 1:541){
  int_lower=Pi_wls_all[,j,]-qnorm(alpha)*sd_wls_all[,j,]
  int_upper=Pi_wls_all[,j,]+qnorm(alpha)*sd_wls_all[,j,]
  
  tab1=sapply(1:5,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  })
  
  
  cp_wls=rbind(cp_wls,tab1)
}


