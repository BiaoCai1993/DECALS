############################################################
############ Constraint OLS ################################
############################################################
library(quadprog)
library(gtools)
# sig: signature matrix, genes * cell types
# bulk: bulk, genes * samples
constraint_ols <- function(sig, bulk){ 
  p=dim(sig)[1];K=dim(sig)[2];n=dim(bulk)[2]
  A=cbind(rep(1,K),diag(K))
  b0=c(1,rep(0,K))
  D=t(sig)%*%sig/n
  
  frac0=NULL
  for(i in 1:n){
    d=bulk[,i]%*%sig/n
    frac0=rbind(frac0,solve.QP(Dmat=D/abs(max(d)),dvec=d/abs(max(d)),Amat=A,bvec=b0,meq=1)$solution)
  }
  #frac0[frac0<0]=0
  return(frac0)
}

###################################################################
################### Proportion Variance ###########################
###################################################################
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

near_posmat <- function(x,e0=0.001){
  x0=(x+t(x))/2
  dimx=nrow(x0)
  va=eigen(x0)$values
  va[va<0]=0
  ve=eigen(x0)$vectors
  z=ve%*%diag(va)%*%t(ve)
  z+e0*diag(rep(1,dimx))
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

sigma_error0 <- function(sigma_all1,sigma_all2){
  K=length(sigma_all1)
  err0=NULL
  for(k in 1:K){
    err0=c(err0,sqrt(sum((sigma_all1[[k]]-sigma_all2[[k]])^2)) )
  }
  err0
}

cor_est <-function(mat0){
  D0=diag(mat0)
  id0=which(D0>0)
  mat1=diag(D0[id0]^{-1/2})%*%mat0[id0,id0]%*%diag(D0[id0]^{-1/2})
  mat=matrix(0,length(D0),length(D0))
  mat[id0,id0]=mat1
  mat
}

var_b <- function(x,sigma){
  p=nrow(x); K=ncol(x)
  A=matrix(1,1,K)
  G_inv=solve(t(x)%*%x/p)
  B=diag(K)-G_inv%*%t(A)%*%solve(A%*%G_inv%*%t(A))%*%A
  var0=B%*%G_inv%*%t(x)%*%sigma%*%x%*%G_inv%*%t(B)/p
  var0/p
}

var_b1 <-function(x,sigma0){
  var0=sigma0*solve(t(x)%*%x/p)
  var0/p
}

propor_cov0 <- function(W,sigma_all=NULL,propor_all,sigma_val=NULL){ # the estimation of var(pi)
  K=nrow(propor_all)
  n=ncol(propor_all)
  p=ncol(W)
  if(is.null(sigma_val)){
    all_cov=sapply(1:n,function(i){
      sigma0=sigma_all[[1]]*propor_all[1,i]^2
      for(k in 2:K) sigma0=sigma0+sigma_all[[k]]*propor_all[k,i]^2
      #sigma1=near_posmat(sigma0,e0=0)
      cov0=var_b(x=t(W),sigma=sigma0)
      diag(cov0)
    })
  }else{
    all_cov=sapply(1:n,function(i){
      sigma0=sigma_val[i]*diag(p)
      cov0=var_b(x=t(W),sigma=sigma0)
      diag(cov0)
    })
  }
  all_cov
}


propor_cov_1 <- function(W,propor_all,sigma_val=NULL){ # the estimation of var(pi)
  K=nrow(propor_all)
  n=ncol(propor_all)
  p=ncol(W)
  all_cov=sapply(1:n,function(i){
    sigma0=sigma_val[i]
    cov0=var_b1(x=t(W),sigma0=sigma0)
    diag(cov0)
  })
  all_cov
}


propor_cov <- function(W,sigma_all=NULL,propor_all,sigma_val=NULL){ # the estimation of var(pi)
  K=nrow(propor_all)
  n=ncol(propor_all)
  p=ncol(W)
  if(is.null(sigma_val)){
    all_cov=lapply(1:n,function(i){
      sigma0=sigma_all[[1]]*propor_all[1,i]^2
      for(k in 2:K) sigma0=sigma0+sigma_all[[k]]*propor_all[k,i]^2
      
      cov0=var_b(x=t(W),sigma=sigma0)
      cov0
    })
  }else{
    all_cov=lapply(1:n,function(i){
      sigma0=sigma_val[i]*diag(p)
      
      cov0=var_b(x=t(W),sigma=sigma0)
      cov0
    })
  }
  all_cov
}

var_cal <-function(x){
  n=length(x)
  sum(x^2)/(n-1)
}



sigma_gen <- function(p,rho,type){
  if(type=="MA"){
    sigma0=matrix(rho,p,p)
    sigma0[abs(row(sigma0)-col(sigma0))>1]<-0
    diag(sigma0)=1
  }else if(type=="AR"){
    sigma0=diag(p)
    for(j in 1:(p-1)){
      sigma0[abs(row(sigma0)-col(sigma0))==j]<-rho*0.9^(j-1)
    }
  }else if(type=="full"){
    sigma0=matrix(rho,p,p)
    diag(sigma0)=1
  }
  return(sigma0)
}

cor_est <-function(mat0){
  D0=diag(mat0)
  id0=which(D0>0)
  mat1=diag(D0[id0]^{-1/2})%*%mat0[id0,id0]%*%diag(D0[id0]^{-1/2})
  mat=matrix(0,length(D0),length(D0))
  mat[id0,id0]=mat1
  mat
}

pi_random_cor <- function(pi_est,pi_cov_est){
  K=length(pi_est)
  term1=matrix(0,K,K)
  for(i in 1:K){
    for(j in 1:K){
      term1[i,j]=pi_cov_est[i,i]*pi_cov_est[j,j]+2*pi_cov_est[i,j]^2
    }
  }
  #diag(term1)=3*diag(pi_cov_est)
  
  term2=4*pi_est%*%t(pi_est)*pi_cov_est
  
  term3=pi_est^2%*%t(diag(pi_cov_est))+t(pi_est^2%*%t(diag(pi_cov_est)))
  term=pi_est^2%*%t(pi_est^2)-(term1+term2+term3)
  term
}

sigma_upt0 <- function(propor_all,y_all,W,lambda_all=lambda_all){
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
  
  ########### initialization ##############
  var0=apply(y_center,2,var_cal)
  pi_cov_all=propor_cov(W,propor_all=propor_all,sigma_val=var0)
  pi_cov_all0=sapply(1:n,function(i)diag(pi_cov_all[[i]]))
  pi_cov_all_before=matrix(1,K,n)
  err=sum((pi_cov_all0-pi_cov_all_before)^2)
  t=0
  while(err>10^{-5} & t<20){
    
    ########## correction of B1 ########
    term_all=lapply(1:n,function(i)pi_random_cor(pi_est=propor_all[,i],pi_cov_est=pi_cov_all[[i]]))
    G1=term_all[[1]]
    for(j in 2:n) G1=G1+term_all[[j]]
    
    ######### correction of B2 ########
    propor_all_adj = propor_all^2-pi_cov_all0
    
    b=apply(y_center,2,function(x) x%*%t(x) )
    
    var_y0=matrix(0,K,p^2)
    for(j in 1:20){
      var_y0=var_y0+propor_all_adj[,id[[j]]]%*%t(b[,id[[j]]])
    }
    
    var_yk<-solve(G1)%*%var_y0
    
    sigma_all0=lapply(1:K,function(k){
      mat0=matrix(var_yk[k,],p,p)
    })
    
    sigma_all=lapply(1:K,function(k){
      cor_filter(sigma_all0[[k]],lambda0=lambda_all[k])
    })
    
    pi_cov_all=propor_cov(W=W,sigma_all=sigma_all,propor_all=propor_all)
    pi_cov_all_before=pi_cov_all0
    pi_cov_all0=sapply(1:n,function(i)diag(pi_cov_all[[i]]))
    err=sum((pi_cov_all0-pi_cov_all_before)^2)
    t=t+1
  }
  
  return(sigma_all)
}


tune_select <- function(y_all,W,propor_est,lambda_set,tune_method="sequential"){
  p=nrow(y_all)
  n=ncol(y_all)
  K=nrow(W)
  n1=n%/%2
  id1=sample(1:n,n1)
  id2=(1:n)[-id1]
  y_all1=y_all[,id1]
  y_all2=y_all[,id2]
  if(tune_method=="sequential"){
    n_lam=length(lambda_set)
    err=matrix(0,K,n_lam)
    lambda_opt=rep(lambda_set[1],K)
    for(k in 1:K){
      for(j in 1:n_lam){
        lambda_opt[k]=lambda_set[j]
        res1=sigma_upt0(propor_all=propor_est[,id1],y_all=y_all1,W,lambda_all=rep(0,K))
        res2=sigma_upt0(propor_all=propor_est[,id2],y_all=y_all2,W,lambda_all=lambda_opt)
        err1=sigma_error0(sigma_all1=res1,sigma_all2=res2)[k]
        res1=sigma_upt0(propor_all=propor_est[,id1],y_all1,W,lambda_all=lambda_opt)
        res2=sigma_upt0(propor_all=propor_est[,id2],y_all2,W,lambda_all=rep(0,K))
        err2=sigma_error0(sigma_all1=res1,sigma_all2=res2)[k]
        err[k,j]=err1+err2
      }
      j_id=which.min(err[k,])[1]
      lambda_opt[k]=lambda_set[j_id]
    }
    return(list(lambda_opt,err))
  }
}
############################################################
############ Data Setting  #################################
############################################################
data_gen1 <- function(n,p,k,seedN){
  library(mvtnorm)
  id1=1:(p/3);id2=(p/3+1):(2*p/3);id3=(2*p/3+1):p
  
  set.seed(1)
  W<-matrix(rnorm(k*p,0,1),k,p)
  Pi<-t(rdirichlet(n,c(3,2,1)))
  D=W%*%t(W)
  A=cbind(rep(1,3), diag(3))
  b0=c(1,rep(0,3))
  
  sigma0_1=sigma_gen(p=100,rho=0.7,type="AR")
  sigma0_2=matrix(0.3,p/3,p/3);diag(sigma0_2)=1
  
  
  sigma1=sigma2=sigma3=matrix(0,p,p)
  sigma1[id1,id1]=sigma0_2;sigma1[id2,id2]=sigma1[id3,id3]=sigma0_1
  sigma2[id2,id2]=sigma0_2;sigma2[id1,id1]=sigma2[id3,id3]=sigma0_1
  sigma3[id3,id3]=sigma0_2;sigma3[id1,id1]=sigma3[id2,id2]=sigma0_1
  
  sigma1=10*sigma1
  sigma2=10*sigma2
  sigma3=10*sigma3
  sigma_all=list(sigma1,sigma2,sigma3)
  set.seed(seedN)
  y_all=matrix(0,p,n)
  Pi_wt<-matrix(0,k,n)
  for(i in 1:n){
    sigma<-Pi[1,i]^2*sigma1+Pi[2,i]^2*sigma2+Pi[3,i]^2*sigma3
    epsilon<-matrix(rmvnorm(1,mean=rep(0,p),sigma=sigma),p,1)
    y_all[,i]<-t(W)%*%Pi[,i]+epsilon
  }
  return(list(bulk=y_all,sig=t(W)))
}


data_gen1_noise<- function(n,p,k,seedN,noise){
  library(mvtnorm)
  id1=1:(p/3);id2=(p/3+1):(2*p/3);id3=(2*p/3+1):p
  
  set.seed(1)
  W0<-matrix(rnorm(k*p,0,1),k,p)
  W=W0+matrix(rnorm(k*p,0,noise),k,p)
  Pi<-t(rdirichlet(n,c(3,2,1)))
  D=W%*%t(W)
  A=cbind(rep(1,3), diag(3))
  b0=c(1,rep(0,3))
  
  sigma0_1=sigma_gen(p=100,rho=0.7,type="AR")
  sigma0_2=matrix(0.3,p/3,p/3);diag(sigma0_2)=1
  
  
  sigma1=sigma2=sigma3=matrix(0,p,p)
  sigma1[id1,id1]=sigma0_2;sigma1[id2,id2]=sigma1[id3,id3]=sigma0_1
  sigma2[id2,id2]=sigma0_2;sigma2[id1,id1]=sigma2[id3,id3]=sigma0_1
  sigma3[id3,id3]=sigma0_2;sigma3[id1,id1]=sigma3[id2,id2]=sigma0_1
  
  sigma1=10*sigma1
  sigma2=10*sigma2
  sigma3=10*sigma3
  sigma_all=list(sigma1,sigma2,sigma3)
  set.seed(seedN)
  y_all=matrix(0,p,n)
  Pi_wt<-matrix(0,k,n)
  for(i in 1:n){
    sigma<-Pi[1,i]^2*sigma1+Pi[2,i]^2*sigma2+Pi[3,i]^2*sigma3
    epsilon<-matrix(rmvnorm(1,mean=rep(0,p),sigma=sigma),p,1)
    y_all[,i]<-t(W)%*%Pi[,i]+epsilon
  }
  return(list(bulk=y_all,sig=t(W)))
}


rmultigamma <- function(n,p,cor_mat,shape,scale){
  cor_mat_up <- chol(cor_mat)
  copula <- matrix(rnorm(n*p),nrow = n)
  copula <- pnorm(t(cor_mat_up)%*%t(copula))
  cor_exp_matrix <- matrix(NA, nrow=n, ncol=p)
  for (i in 1:n){
    cor_exp_matrix[i,] <- qgamma(copula[,i], shape=shape,scale=scale)
    id0=which(scale==0)
    cor_exp_matrix[i,id0]=0
  }
  cor_exp_matrix
}


data_gen2 <- function(frac0,sig,sigma_ols0,seedN){
  p=nrow(sig)
  W=t(sig)
  Pi=t(frac0)
  n=nrow(frac0)
  
  set.seed(1)
  D=W%*%t(W)
  A=cbind(rep(1,5), diag(5))
  b0=c(1,rep(0,5))
  
  
  sigma1=cor_est(sigma_ols0[[1]]+10^(-6)*diag(rep(1,159)))
  sigma2=cor_est(sigma_ols0[[2]]+10^(-6)*diag(rep(1,159)))
  sigma3=cor_est(sigma_ols0[[3]]+10^(-6)*diag(rep(1,159)))
  sigma4=cor_est(sigma_ols0[[4]]+10^(-6)*diag(rep(1,159)))
  sigma5=cor_est(sigma_ols0[[5]]+10^(-6)*diag(rep(1,159)))
  
  set.seed(seedN)
  y_all1=rmultigamma(n,p,cor_mat=sigma1,shape=1/100,scale=W[1,]*100)
  y_all2=rmultigamma(n,p,cor_mat=sigma2,shape=1/100,scale=W[2,]*100)
  y_all3=rmultigamma(n,p,cor_mat=sigma3,shape=1/100,scale=W[3,]*100)
  y_all4=rmultigamma(n,p,cor_mat=sigma4,shape=1/100,scale=W[4,]*100)
  y_all5=rmultigamma(n,p,cor_mat=sigma5,shape=1/100,scale=W[5,]*100)
  y_all<-matrix(0,p,n)
  for(i in 1:n){
    y_all[,i]<-Pi[1,i]*y_all1[i,]+Pi[2,i]*y_all2[i,]+Pi[3,i]*y_all3[i,]+Pi[4,i]*y_all4[i,]+Pi[5,i]*y_all5[i,]
  }
  return(list(bulk=y_all,sig=t(W)))
}


data_gen3 <- function(frac0,sig,sigma_ols0,seedN){
  k=ncol(sig)
  W=t(sig)
  p=nrow(sig)
  n=nrow(frac0)
  
  Pi=t(frac0)
  sigma_all=sigma_ols0
  set.seed(seedN)
  
  y_all<-matrix(0,p,n)
  Pi_ols<-matrix(0,k,n)
  for(i in 1:n){
    sigma<-Pi[1,i]^2*sigma_all[[1]]+Pi[2,i]^2*sigma_all[[2]]+Pi[3,i]^2*sigma_all[[3]]+Pi[4,i]^2*sigma_all[[4]]+Pi[5,i]^2*sigma_all[[5]]
    epsilon<-matrix(rmvnorm(1,mean=rep(0,p),sigma=sigma),p,1)
    y_all[,i]<-t(W)%*%Pi[,i]+epsilon
  }
  return(list(bulk=y_all,sig=t(W)))
}



############################################################
##################### OLS  #################################
############################################################
sigma_upt <- function(propor_all,y_all,W,lambda_all,a=3.7,e0=0){
  K<-nrow(propor_all); p<-nrow(y_all); n=ncol(propor_all)
  # centralization
  y_center <- y_all-t(W)%*%propor_all
  #var_y <- apply(y_center,2,function(x) x%*%t(x) )
  dis=n%/%5
  id=vector("list")
  for(j in 1:4){
    id[[j]]=(j*dis-dis+1):(j*dis)
  }
  id[[5]]=(4*dis+1):n
  var_y0=matrix(0,K,p^2)
  for(j in 1:5){
    var_y0=var_y0+(propor_all^2)[,id[[j]]]%*%t(apply(y_center[,id[[j]]],2,function(x) x%*%t(x)))
  }
  
  var_yk<-solve(propor_all^2%*%t(propor_all^2))%*%var_y0
  
  sigma_all=lapply(1:K,function(k){
    mat0=matrix(var_yk[k,],p,p)
    cor_filter(mat0,lambda0=lambda_all[k])
  })
  return(sigma_all)
}

tune_select0 <- function(y_all,W,propor_est,lambda_set,tune_method="sequential"){
  p=nrow(y_all)
  n=ncol(y_all)
  K=nrow(W)
  n1=n%/%2
  id1=sample(1:n,n1)
  id2=(1:n)[-id1]
  y_all1=y_all[,id1]
  y_all2=y_all[,id2]
  if(tune_method=="sequential"){
    n_lam=length(lambda_set)
    err=matrix(0,K,n_lam)
    lambda_opt=rep(lambda_set[1],K)
    for(k in 1:K){
      for(j in 1:n_lam){
        lambda_opt[k]=lambda_set[j]
        res1=sigma_upt(propor_all=propor_est[,id1],y_all1,W,lambda_all=rep(0,K),a,e0)
        res2=sigma_upt(propor_all=propor_est[,id2],y_all2,W,lambda_all=lambda_opt,a,e0)
        err1=sigma_error0(sigma_all1=res1,sigma_all2=res2)[k]
        res1=sigma_upt(propor_all=propor_est[,id1],y_all1,W,lambda_all=lambda_opt,a,e0)
        res2=sigma_upt(propor_all=propor_est[,id2],y_all2,W,lambda_all=rep(0,K),a,e0)
        err2=sigma_error0(sigma_all1=res1,sigma_all2=res2)[k]
        err[k,j]=err1+err2
      }
      j_id=which.min(err[k,])[1]
      lambda_opt[k]=lambda_set[j_id]
    }
    return(list(lambda_opt,err))
  }
}

ols_result<-function(bulk,sig){
  p=dim(sig)[1];K=dim(sig)[2];n=dim(bulk)[2]
  A=cbind(rep(1,K),diag(K))
  b0=c(1,rep(0,K))
  D=t(sig)%*%sig/n
  W=t(sig)
  y_all=bulk
  Pi_ols=NULL
  for(i in 1:n){
    d=bulk[,i]%*%sig/n
    Pi_ols=rbind(Pi_ols,solve.QP(Dmat=D/abs(max(d)),dvec=d/abs(max(d)),Amat=A,bvec=b0,meq=1)$solution)
  }
  
  lambda_ols<- tune_select0(y_all=bulk,W=t(sig),propor_est=t(Pi_ols),lambda_set=seq(0,0.5,by=0.025),tune_method="sequential")

  sigma_ols<-sigma_upt(propor_all=t(Pi_ols),y_all=y_all,W=t(sig),lambda_all=lambda_ols[[1]]) 
  
  y_center <- y_all-sig%*%t(Pi_ols)
  var0=apply(y_center,2,var_cal)
  CTS_proportion_var=propor_cov_1(W,propor_all=t(Pi_ols),sigma_val=var0)
  return(CTS_proportion=Pi_ols,CTS_proportion_var=CTS_proportion_var)
}








