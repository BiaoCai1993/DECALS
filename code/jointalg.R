###############################################################################
###### The algorithm to jointly estimate proportion and tsparse variance ######
###############################################################################

propor_est0 <- function(W,y,sigma){
  solve(W%*%solve(sigma)%*%t(W))%*%W%*%solve(sigma)%*%y
}

propor_est <- function(W,y,sigma){
  aa=solve(W%*%solve(sigma)%*%t(W))
  propor0=aa%*%W%*%solve(sigma)%*%y
  A=matrix(1,1,nrow(W))
  propor0-aa%*%t(A)%*%solve(A%*%aa%*%t(A))%*%(A%*%propor0-1)
}

propor_est <- function(W,y,sigma){
  bb=solve(sigma)
  aa=solve(W%*%bb%*%t(W))
  propor0=aa%*%W%*%bb%*%y
  A=matrix(1,1,nrow(W))
  propor0-aa%*%t(A)%*%solve(A%*%aa%*%t(A))%*%(A%*%propor0-1)
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
  mat1=diag(D0[id0]^{-1/2})%*%mat0[id0,id0]%*%diag(D0[id0]^{-1/2})
  mat2=apply(mat1,c(1,2),scad_filter,lambda=lambda0,a=a)
  mat3=diag(D0[id0]^{1/2})%*%mat2%*%diag(D0[id0]^{1/2})
  mat=matrix(0,length(D0),length(D0))
  mat[id0,id0]=mat3
  mat
}

sigma_upt <- function(propor_all,y_all,W,lambda_all,a,e0){
  # centralization
  y_center <- y_all-t(W)%*%propor_all
  var_y <- apply(y_center,2,function(x) x%*%t(x) )
  #var_yk <- sapply(1:nrow(var_y),function(l){
  #  lm0=lm(var_y[l,]~t(propor_all^2)-1)
  #  lm0$coefficients
  #})
  var_yk=solve(propor_all^2%*%t(propor_all^2))%*%propor_all^2%*%t(var_y)
  K=nrow(propor_all)
  p=nrow(y_all)
  sigma_all=lapply(1:K,function(k){
    mat0=matrix(var_yk[k,],p,p)
    cor_filter(mat0,lambda0=lambda_all[k])
  })
  return(sigma_all)
}



joint_model <- function(W,y_all,lambda_all,a=3.7,e0=0.001){ # proposed algorithm
  p=dim(y_all)[1]
  n=dim(y_all)[2]
  K=dim(W)[1]
  sigma_est=lapply(1:K,function(i)diag(1,p))
  pi_est=sapply(1:n,function(i)rep(1/K,K))
  #pi_est=propor_all
  #sigma_est=sigma#cov_est(W,y_all,propor_all=pi_est)
  A=cbind(rep(1,3), diag(3))
  b0=c(1,rep(0,3))
  
  delta=10
  t=1
  while(delta>0.001 & t<=10){
    pi_est1=sapply(1:n,function(i){
      sigma0=sigma_est[[1]]*pi_est[1,i]^2
      for(k in 2:K) sigma0=sigma0+sigma_est[[k]]*pi_est[k,i]^2
      sigma0=near_posmat(sigma0,e0)
      D=W%*%solve(sigma0)%*%t(W)
      d=y_all[,i]%*%solve(sigma0)%*%t(W)
      solve.QP(Dmat=D/p,dvec=d/p,Amat=A,bvec=b0,meq=1)$solution
    })
    
    sigma_est1=sigma_upt(propor_all=pi_est1,y_all,W,lambda_all=lambda_all,a,e0)
    delta0=sum((pi_est1-pi_est)^2)
    delta1=sigma_error0(sigma_all1=sigma_est1,sigma_all2=sigma_est)#sapply(1:K,function(i)sum((sigma_est1[[i]]-sigma_est[[i]])^2)/sum(sigma_est[[i]]^2))
    delta=delta0+sum(delta1)
    pi_est=pi_est1
    sigma_est=sigma_est1
    #if(eigen(sigma_est)$values[p]<10^(-12)) return(list(pi_est=pi_est1,sigma_est=sigma_est1,iteration=t))
    t=t+1
  }
  return(list(pi_est=pi_est1,sigma_est=sigma_est1,iteration=t))
}


sigma_error <- function(sigma_all1,sigma_all2){
  K=length(sigma_all1)
  err0=NULL
  for(k in 1:K){
    err0=c(err0,sqrt(sum((sigma_all1[[k]]-sigma_all2[[k]])^2)) )
  }
  err0
}

tune_select <- function(y_all,W,lambda_set,tune_method="sequential"){
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
        res1=joint_model(W,y_all1,lambda_all=rep(0,K))
        res2=joint_model(W,y_all2,lambda_all=lambda_opt)
        err1=sigma_error(sigma_all1=res1$sigma_est,sigma_all2=res2$sigma_est)[k]
        res1=joint_model(W,y_all1,lambda_all=lambda_opt)
        res2=joint_model(W,y_all2,lambda_all=rep(0,K))
        err2=sigma_error(sigma_all1=res1$sigma_est,sigma_all2=res2$sigma_est)[k]
        err[k,j]=err1+err2
      }
      j_id=which.min(err[k,])[1]
      lambda_opt[k]=lambda_set[j_id]
    }
    return(list(lambda_opt,err))
  }
}
