#######################################################################
########################## Simulation #################################
#######################################################################
path="/home/bc758/research/decals/"
source(paste0(path,"code/DECALS.R"))


propor_adjust <- function(x,method="all"){
  K=nrow(x);n=ncol(x)
  if(method=="all"){
    propor_adj=matrix(0,K*(K+1)/2,n)
    i1=1
    for(k1 in 1:K){
      for(k2 in k1:K){
        propor_adj[i1,]=propor_all[k1,]*propor_all[k2,]
        i1=i1+1
      }
    }
  }
  if(method=="ind"){
    propor_adj=propor_all^2
  }
  propor_adj
}

scad_filter <- function(x,lambda,a=3.7){
  signx=sign(x)
  if(abs(x)<= 2*lambda){
    x0=signx*max(abs(x)-lambda,0)
  }else if(abs(x)<a*lambda){
    x0=signx*(lambda+(a-1)/(a-2)*(abs(x)-2*lambda) )
  }else{
    x0=signx*abs(x)
  }
  return(x0)
}


cor_filter <- function(mat0,lambda0,a=3.7){
  D0=diag(mat0)
  id0=which(D0>0)
  if(length(id0)>1){
    mat1=diag(D0[id0]^{-1/2})%*%mat0[id0,id0]%*%diag(D0[id0]^{-1/2})
    mat1[mat1>1]=0;mat1[mat1< -1]=0;diag(mat1)=1
    mat2=apply(mat1,c(1,2),scad_filter,lambda=lambda0,a=a)
    mat3=diag(D0[id0]^{1/2})%*%mat2%*%diag(D0[id0]^{1/2})
    mat=matrix(0,length(D0),length(D0))
    mat[id0,id0]=mat3
  }else{
    mat=mat0
  }
  mat
}

sigma_upt2 <- function(propor_all,y_all,W=W,method="all"){
  K<-nrow(propor_all); p<-nrow(y_all); n=ncol(propor_all)
  # centralization
  y_center <- y_all-t(W)%*%propor_all
  #var_y <- apply(y_center,2,function(x) x%*%t(x) )
  dis=n%/%20
  id=vector("list")
  for(j in 1:19){
    id[[j]]=(j*dis-dis+1):(j*dis)
  }
  id[[20]]=(19*dis+1):n
  
  b=apply(y_center,2,function(x) x%*%t(x) )
  
  propor_adj=propor_adjust(propor_all,method)
  var_y0=matrix(0,nrow(propor_adj),p^2)
  for(j in 1:20){
    var_y0=var_y0+propor_adj[,id[[j]]]%*%t(b[,id[[j]]])
  }
  
  G=propor_adj%*%t(propor_adj)
  
  var_yk<-solve(G)%*%var_y0
  
  sigma_all=lapply(1:nrow(propor_adj),function(k){
    mat0=matrix(var_yk[k,],p,p)
  })
  
  return(sigma_all)
}



load(file=paste0(path,"data/sc/linear/ct_props_all1.RData"))
load(file=paste0(path,"data/rosmap/sig.RData"))
load(file=paste0(path,"data/sc/linear/cell_cov.RData"))
library(MASS)

pseudo_bulk0=sig%*%t(ct_props_all)
sigmai_all=array(0,dim=c(159,159,436))
for(i in 1:436){
  sigmai0=matrix(0,159,159)
  for(k1 in 1:5){
    for(k2 in 1:5){
      sigmai0=sigmai0+ct_props_all[i,k1]*ct_props_all[i,k2]*cell_cov[,,k1,k2]
    }
  }
  sigmai_all[,,i]=sigmai0
}
#save(sigmai_all,file="/home/bc758/research/decals/data/ROSMAP/sc/new/sigmai_all.RData")

pseudo_bulk=pseudo_bulk0
for(i in 1:436){
  pseudo_bulk[,i]=pseudo_bulk0[,i]+mvrnorm(n=1,mu=rep(0,159),Sigma=sigmai_all[,,i]/10^5)
}
cell_names=c("Neu","Oli","Ast","Mic","End")
W=sig
frac0=constraint_ols(sig=W,bulk=pseudo_bulk)
colnames(frac0)=cell_names
boxplot(frac0)


par(mfrow=c(2,3))
for(k in 1:5){
  plot(ct_props_all[,k],frac0[,k],main=cell_names[k],xlab="TRUE",ylab="EST")
}
par(mfrow=c(1,1))



sigma_all1=sigma_upt2(propor_all=t(frac0),y_all=pseudo_bulk,W=t(W),method="all")
sigma_all2=sigma_upt2(propor_all=t(frac0),y_all=pseudo_bulk,W=t(W),method="ind")


prop_cov_all=matrix(0,436,5)
prop_cov_diag=matrix(0,436,5)
for(sam_id in 1:436){
  cell_cov_all=cell_cov_diag=matrix(0,159,159)
  i1=1
  for(k1 in 1:5){
    for(k2 in k1:5){
      cell_cov_all=cell_cov_all+frac0[sam_id,k1]*frac0[sam_id,k2]*sigma_all1[[i1]]
      #if(k2==k1) cell_cov_diag=cell_cov_diag+frac0[sam_id,k1]*frac0[sam_id,k2]*sigma_all1[[i1]]
      i1=i1+1
    }
  }
  for(k in 1:5){
    cell_cov_diag=cell_cov_diag+frac0[sam_id,k]^2*sigma_all2[[k]]
  }
  prop_cov_all[sam_id,]=diag(var_b(x=W,sigma=cell_cov_all))
  prop_cov_diag[sam_id,]=diag(var_b(x=W,sigma=cell_cov_diag))
}

par(mfrow=c(2,3))
for(k in 1:5){
  plot(log(abs(prop_cov_all[,k])),log(abs(prop_cov_diag[,k])),main=cell_names[k],xlab="ALL",ylab="Within")
  cor_val=round(cor(log(abs(prop_cov_all[,k])),log(abs(prop_cov_diag[,k]))),3)
  legend("topleft",legend=cor_val)
}
par(mfrow=c(1,1))


n=436
k=5
var_all1=var_all2=array(0,dim=c(5,n,100))
Pi_ols_all=array(0,dim=c(5,n,100))
for(seedN in 1:100){
  set.seed(seedN)
  pseudo_bulk=pseudo_bulk0
  for(i in 1:436){
    pseudo_bulk[,i]=pseudo_bulk0[,i]+mvrnorm(n=1,mu=rep(0,159),Sigma=sigmai_all[,,i]/10^6)
  }
  
  W=sig
  frac0=constraint_ols(sig=W,bulk=pseudo_bulk)
  colnames(frac0)=cell_names
  Pi_ols_all[,,seedN]=t(frac0)
  frac0=ct_props_all
  sigma_all1=sigma_upt2(propor_all=t(frac0),y_all=pseudo_bulk,W=t(W),method="all")
  sigma_all2=sigma_upt2(propor_all=t(frac0),y_all=pseudo_bulk,W=t(W),method="ind")
  
  prop_cov_all=matrix(0,436,5)
  prop_cov_diag=matrix(0,436,5)
  for(sam_id in 1:436){
    cell_cov_all=cell_cov_diag=matrix(0,159,159)
    i1=1
    for(k1 in 1:5){
      for(k2 in k1:5){
        cell_cov_all=cell_cov_all+frac0[sam_id,k1]*frac0[sam_id,k2]*sigma_all1[[i1]]
        #if(k2==k1) cell_cov_diag=cell_cov_diag+frac0[sam_id,k1]*frac0[sam_id,k2]*sigma_all1[[i1]]
        i1=i1+1
      }
    }
    for(k in 1:5){
      cell_cov_diag=cell_cov_diag+frac0[sam_id,k]^2*sigma_all2[[k]]
    }
    prop_cov_all[sam_id,]=diag(var_b(x=W,sigma=cell_cov_all))
    prop_cov_diag[sam_id,]=diag(var_b(x=W,sigma=cell_cov_diag))
  }
  var_all1[,,seedN]=t(prop_cov_all)
  var_all2[,,seedN]=t(prop_cov_diag)
}


var_all1s=var_all1
var_all2s=var_all2
Pi=t(ct_props_all)

#var_all1[var_all1<0]=0
#var_all2[var_all2<0]=0
var_all1=abs(var_all1s)
var_all2=abs(var_all2s)

alpha=0.975
K=5
tab1_all=array(0,dim=c(K,2,n))
for(j in 1:n){
  int_lower=Pi_ols_all[,j,]-qnorm(alpha)*sqrt(var_all1[,j,])
  int_upper=Pi_ols_all[,j,]+qnorm(alpha)*sqrt(var_all1[,j,])
  
  tab1=sapply(1:K,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  })
  
  int_lower=Pi_ols_all[,j,]-qnorm(alpha)*sqrt(var_all2[,j,])
  int_upper=Pi_ols_all[,j,]+qnorm(alpha)*sqrt(var_all2[,j,])
  
  tab1=cbind(tab1,sapply(1:K,function(k){
    sum(sapply(1:100,function(i) sum(Pi[k,j]<int_upper[k,i] & Pi[k,j]>int_lower[k,i]) ))/100
  }))

  tab1_all[,,j]=tab1
}

load(file=paste0(path,"data/sc/linear/tab1_all.RData"))

cell_names=c("Neu","Oli","Ast","Mic","End")
library(vioplot)
par(mfrow=c(1,5))
for(k in 1:5){
  result=t(tab1_all[k,c(1:2),])#c(5,4,2,3)
  colnames(result)=c("Within","ALL")
  vioplot(result,col=c("lightcoral","skyblue"),
          axes=FALSE,main=cell_names[k],ylim=c(0.4,1),
          ylab="coverage probability")
  lines(c(0.5,2.5),c(0.95,0.95),lty=2)
}
par(mfrow=c(1,1))

par(mfrow=c(2,3))
for(k in 1:5){
  plot(tab1_all[k,2,],tab1_all[k,1,],main=cell_names[k],xlim=c(0.8,1),ylim=c(0.8,1),xlab="ALL",ylab="Within")
  cor_val=cor(tab1_all[k,2,],tab1_all[k,1,])
  text(0.85,0.98,labels=paste("cor =",round(cor_val,3)))
}
par(mfrow=c(1,1))
