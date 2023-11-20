############################################################
############ Constraint OLS ################################
############################################################
library(quadprog)
# sig: signature matrix, genes * cell types
# bulk: bulk, genes * samples
# wt: withts, samples * genes
constraint_ols <- function(sig, bulk,wt){ 
  p=dim(sig)[1];K=dim(sig)[2];n=dim(bulk)[2]
  A=cbind(rep(1,K),diag(K))
  b0=c(1,rep(0,K))
  
  
  frac0=NULL
  for(i in 1:n){
    D=t(sig)%*%diag(wt[i,])%*%sig/n
    d=bulk[,i]%*%diag(wt[i,])%*%sig/n
    frac0=rbind(frac0,solve.QP(Dmat=D/100,dvec=d/100,Amat=A,bvec=b0,meq=1)$solution)
  }
  frac0[frac0<0]=0
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

var_b <- function(x,sigma,wt){
  p=nrow(x); K=ncol(x)
  A=matrix(1,1,K)
  G_inv=solve(t(x)%*%diag(wt)%*%x/p)
  B=diag(K)-G_inv%*%t(A)%*%solve(A%*%G_inv%*%t(A))%*%A
  var0=B%*%G_inv%*%t(x)%*%diag(wt)%*%sigma%*%diag(wt)%*%x%*%G_inv%*%t(B)/p
  var0/p
}

var_b1 <-function(x,sigma0,wt){
  var0=sigma0*solve(t(x)%*%diag(wt)%*%x/p)
  var0/p
}

propor_cov0 <- function(W,sigma_all=NULL,propor_all,sigma_val=NULL,wt=wt){ # the estimation of var(pi)
  K=nrow(propor_all)
  n=ncol(propor_all)
  p=ncol(W)
  if(is.null(sigma_val)){
    all_cov=sapply(1:n,function(i){
      sigma0=sigma_all[[1]]*propor_all[1,i]^2
      for(k in 2:K) sigma0=sigma0+sigma_all[[k]]*propor_all[k,i]^2
      #sigma1=near_posmat(sigma0,e0=0)
      cov0=var_b(x=t(W),sigma=sigma0,wt=wt[i,])
      diag(cov0)
    })
  }else{
    all_cov=sapply(1:n,function(i){
      sigma0=sigma_val[i]*diag(p)
      cov0=var_b(x=t(W),sigma=sigma0,wt=wt[i,])
      diag(cov0)
    })
  }
  all_cov
}


propor_cov_1 <- function(W,propor_all,sigma_val=NULL,wt=wt){ # the estimation of var(pi)
  K=nrow(propor_all)
  n=ncol(propor_all)
  p=ncol(W)
  all_cov=sapply(1:n,function(i){
    sigma0=sigma_val[i]
    cov0=var_b1(x=t(W),sigma0=sigma0,wt=wt[i,])
    diag(cov0)
  })
  all_cov
}


propor_cov <- function(W,sigma_all=NULL,propor_all,sigma_val=NULL,wt=wt){ # the estimation of var(pi)
  K=nrow(propor_all)
  n=ncol(propor_all)
  p=ncol(W)
  if(is.null(sigma_val)){
    all_cov=lapply(1:n,function(i){
      sigma0=sigma_all[[1]]*propor_all[1,i]^2
      for(k in 2:K) sigma0=sigma0+sigma_all[[k]]*propor_all[k,i]^2
      
      cov0=var_b(x=t(W),sigma=sigma0,wt=wt[i,])
      cov0
    })
  }else{
    all_cov=lapply(1:n,function(i){
      sigma0=sigma_val[i]*diag(p)
      
      cov0=var_b(x=t(W),sigma=sigma0,wt=wt[i,])
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

sigma_upt0 <- function(propor_all,y_all,W,lambda_all=lambda_all,wt=wt){
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
  pi_cov_all=propor_cov(W,propor_all=propor_all,sigma_val=var0,wt=wt)
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
    
    pi_cov_all=propor_cov(W=W,sigma_all=sigma_all,propor_all=propor_all,wt=wt)
    pi_cov_all_before=pi_cov_all0
    pi_cov_all0=sapply(1:n,function(i)diag(pi_cov_all[[i]]))
    err=sum((pi_cov_all0-pi_cov_all_before)^2)
    t=t+1
  }
  
  return(sigma_all)
}


tune_select <- function(y_all,W,wt=wt,propor_est,lambda_set,tune_method="sequential"){
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
        res1=sigma_upt0(propor_all=propor_est[,id1],y_all=y_all1,W,lambda_all=rep(0,K),wt=wt)
        res2=sigma_upt0(propor_all=propor_est[,id2],y_all=y_all2,W,lambda_all=lambda_opt,wt=wt)
        err1=sigma_error0(sigma_all1=res1,sigma_all2=res2)[k]
        res1=sigma_upt0(propor_all=propor_est[,id1],y_all1,W,lambda_all=lambda_opt,wt=wt)
        res2=sigma_upt0(propor_all=propor_est[,id2],y_all2,W,lambda_all=rep(0,K),wt=wt)
        err2=sigma_error0(sigma_all1=res1,sigma_all2=res2)[k]
        err[k,j]=err1+err2
      }
      j_id=which.min(err[k,])[1]
      lambda_opt[k]=lambda_set[j_id]
    }
    return(list(lambda_opt,err))
  }
}

