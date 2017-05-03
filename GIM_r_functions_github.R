### Function to get nlminb

gim<-function(model,iv=get(paste0("iv.",model)),negll=get(paste0("negll.",model)),
              lb=get(paste0("lb.",model)),...){
  fit<-nlminb(start=iv,objective=negll,lower=lb)
  fit$par<-listfunc(model,fit)
  fit
}

### Function to simulate data from the GIM model:

sim.GIM<-function(a,theta,b,c1,c2,tau1,tau0,M1,M2,M1c,M2c,N){ #n is the nr of simulations, d the initial state \in {1,2,3,4}
  
  for (s in 1:3){
    n<-N[s]
    coal.times<-vector(mode="numeric",length=n)
    for (i in 1:n){ 
      d<-s
      repeat{
        if(d==1){
          coal.times[i]<-coal.times[i]+rexp(1,1/c1+M1c)
          if(coal.times[i]>tau1){
            coal.times[i]<-tau1
            break
          }
          if(runif(1)<=(1/c1)/(1/c1+M1c)){
            d<-4
            break
          }else{
            d<-3
          }
        }
        if(d==2){
          coal.times[i]<-coal.times[i]+rexp(1,1/c2+M2c)
          if(coal.times[i]>tau1){
            coal.times[i]<-tau1
            break
          }
          if(runif(1)<=(1/c2)/(1/c2+M2c)){
            d<-4
            break
          }else{
            d<-3
          }  
        }
        if(d==3){
          if(M1c>0 | M2c>0){
            coal.times[i]<-coal.times[i]+rexp(1,(M1c+M2c)/2)
            if(coal.times[i]>tau1){
              coal.times[i]<-tau1
              break
            }
            if(runif(1)<=M2c/(M1c+M2c)){
              d<-1
            }else{
              d<-2        
            }
          }else if(M1c==0 & M2c==0){
            coal.times[i]<-tau1
            break
          }
        }
        
      }
      
      if(d==4){next}
      repeat{
        if(d==1){
          coal.times[i]<-coal.times[i]+rexp(1,1+M1)
          if(coal.times[i]>tau0){
            coal.times[i]<-tau0
            break
          }
          if(runif(1)<=(1)/(1+M1)){
            d<-4
            break
          }else{
            d<-3
          }
        }
        if(d==2){
          coal.times[i]<-coal.times[i]+rexp(1,1/b+M2)
          if(coal.times[i]>tau0){
            coal.times[i]<-tau0
            break
          }
          if(runif(1)<=(1/b)/(1/b+M2)){
            d<-4
            break
          }else{
            d<-3
          }  
        }
        if(d==3){
          if(M1>0 | M2>0){
            coal.times[i]<-coal.times[i]+rexp(1,(M1+M2)/2)  
            
            if(coal.times[i]>tau0){
              coal.times[i]<-tau0
              break
            }
            if(runif(1)<=M2/(M1+M2)){
              d<-1
            }else{
              d<-2        
            }
          }else if(M1==0 & M2==0){
            coal.times[i]<-tau0
            break
          }
        }
        
        
      }
      
      if(d==4){next}else{
        coal.times[i]<-coal.times[i]+rexp(1,1/a)
      }
      
      
    }
    
    mutations<-rpois(length(coal.times),theta*coal.times)
    if(s==1){
      x1<-mutations
      t1<-coal.times
      N1<-N[s]
    }
    if(s==2){
      x2<-mutations
      t2<-coal.times
      N2<-N[s]
    }
    if(s==3){
      x3<-mutations
      t3<-coal.times
      N3<-N[s]
    }
  }
  
  list(x1=x1,x2=x2,x3=x3,t1=t1,t2=t2,t3=t3,N1=N1,N2=N2,N3=N3)
} 


### Function to list and label the maximum-likelihood estimates:

listfunc<-function(model,fit){
  if(model=="iso.1"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                   V=fit$par[4])
  }
  
  if(model=="IM.1"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                   V=fit$par[4],M1=fit$par[5],M2=fit$par[6])
  }
  
  if(model=="IIM.1"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                   T1=fit$par[4],V=fit$par[5],M1=fit$par[6],M2=fit$par[7])
  }
  
  if(model=="IIM.2"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                   theta_c1=fit$par[4],theta_c2=fit$par[5],T1=fit$par[6],V=fit$par[7],
                   M1=fit$par[8],M2=fit$par[9])
  }
  
  if(model=="IIM.3"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                   theta_c1=fit$par[4],theta_c2=fit$par[5],T1=fit$par[6],V=fit$par[7],
                   M2=fit$par[8])
  }
  
  if(model=="IM.2"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                   T1=fit$par[4],M2=fit$par[5])
  }
  
  if(model=="IM.3"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                   T1=fit$par[4],M1=fit$par[5])
  }
  
  if(model=="IM.4"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                   T1=fit$par[4],M=fit$par[5])
  }
  
  if(model=="IIM.4"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                   theta_c1=fit$par[4],theta_c2=fit$par[5],T1=fit$par[6],V=fit$par[7],
                   M1=fit$par[8])
  }
    if(model=="IIM.5"){
      par.list<-list(theta_a=fit$par[1],theta=fit$par[2],
                     theta_c1=fit$par[3],theta_c2=fit$par[4],T1=fit$par[5],V=fit$par[6],
                     M=fit$par[7])
    }
    
  if(model=="IIM.6"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                   theta_c1=fit$par[4],theta_c2=fit$par[5],T1=fit$par[6],V=fit$par[7],
                   M=fit$par[8])
  }
  
    if(model=="iso.2"){
      par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                     theta_c1=fit$par[4],theta_c2=fit$par[5],T1=fit$par[6],V=fit$par[7])
    }
    
    if(model=="GIM.1"){
      par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                     theta_c1=fit$par[4],theta_c2=fit$par[5],T1=fit$par[6],V=fit$par[7],
                     M1=fit$par[8],M2=fit$par[9],M1_prime=fit$par[10],M2_prime=fit$par[11])
    }
  
    if(model=="GIM.2"){
      par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                     theta_c1=fit$par[4],theta_c2=fit$par[5],T1=fit$par[6],V=fit$par[7],
                     M1_prime=fit$par[8],M2_prime=fit$par[9])
    }
  
  if(model=="GIM.3"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                   theta_c1=fit$par[4],theta_c2=fit$par[5],T1=fit$par[6],V=fit$par[7],
                   M=fit$par[8],M_prime=fit$par[9])
  }
  
  if(model=="GIM.4"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],theta_b=fit$par[3],
                   theta_c1=fit$par[4],theta_c2=fit$par[5],T1=fit$par[6],V=fit$par[7],
                   M_prime=fit$par[8])
  }
  
  
  if(model=="IIM.5"){
    par.list<-list(theta_a=fit$par[1],theta=fit$par[2],
                   theta_c1=fit$par[3],theta_c2=fit$par[4],T1=fit$par[5],V=fit$par[6],
                   M=fit$par[7])
  }
  
  return(par.list)
  
}

### Vectors of lower bounds for each model:

lb.iso.1<-0.000001
lb.IM.1<-c(0.000001,0.000001,0.000001,0.000001,0,0)
lb.IIM.1<-c(0.000001,0.000001,0.000001,0.000001,0.000001,0,0)
lb.IIM.2<-c(0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0,0)
lb.IIM.3<-c(0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0)
lb.IM.2<-c(0.000001,0.000001,0.000001,0.000001,0)
lb.IM.3<-c(0.000001,0.000001,0.000001,0.000001,0)
lb.IM.4<-c(0.000001,0.000001,0.000001,0.000001,0)
lb.IIM.4<-c(0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0)
lb.IIM.5<-c(0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0)
lb.IIM.6<-c(0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0)
lb.iso.2<-0.000001
lb.GIM.1<-c(0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0,0,0,0)
lb.GIM.2<-c(0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0,0)
lb.GIM.3<-c(0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0,0)
lb.GIM.4<-c(0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0)

### Vectors of initial values
iv.iso.1<-c(3.2,4)
iv.IM.1<-rep(3.2,6)
iv.IIM.1<-rep(3.2,7)
iv.IIM.2<-rep(3.2,9)
iv.IIM.3<-rep(3.2,8)
iv.IM.2<-rep(3.2,5)
iv.IM.3<-rep(3.2,5)
iv.IM.4<-rep(3.2,5)
iv.IIM.4<-rep(3.2,8)
iv.IIM.5<-rep(3.2,7)
iv.IIM.6<-rep(3.2,8)
iv.iso.2<-rep(3.2,7)
iv.GIM.1<-rep(3.2,11)
iv.GIM.2<-rep(3.2,9)
iv.GIM.3<-rep(3.2,9)
iv.GIM.4<-rep(3.2,8)



### GIM MODELS


# full GIM model

negll.GIM.1<-function(params,model){
  
  theta<-params[2]
  a<-params[1]/theta
  b<-params[3]/theta
  c<-c(params[4]/theta,params[5]/theta)
  tau1<-params[6]/theta
  tau0<-tau1+params[7]/theta
  M<-vector(length=2,mode="numeric")
  M[1]<-params[8]
  M[2]<-params[9]
  Mc<-vector(length=2,mode="numeric")
  Mc[1]<-params[10]
  Mc[2]<-params[11]
  
  R1<-matrix(rep(r1,4),nrow=4,byrow=TRUE)
  R2<-matrix(rep(r2,4),nrow=4,byrow=TRUE)
  R3<-matrix(rep(r3,4),nrow=4,byrow=TRUE)
  
  X1<-matrix(rep(x1,4),nrow=4,byrow=TRUE)
  X2<-matrix(rep(x2,4),nrow=4,byrow=TRUE)
  X3<-matrix(rep(x3,4),nrow=4,byrow=TRUE)
  
  Qt<-matrix(ncol=4,nrow=4)
  Qt[,1]<-c(-(1/c[1]+Mc[1]),0,Mc[1],1/c[1])
  Qt[,2]<-c(0,-(1/c[2]+Mc[2]),Mc[2],1/c[2])
  Qt[,3]<-c(Mc[2]/2,Mc[1]/2,-(Mc[1]+Mc[2])/2,0)
  Qt[,4]<-c(0,0,0,0)
  G<-t(eigen(Qt)$vectors)
  Ginv<-solve(G)
  nu<-eigen(Qt)$values
  
  Qt2<-matrix(ncol=4,nrow=4)
  Qt2[,1]<-c(-(1+M[1]),0,M[1],1)
  Qt2[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
  Qt2[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
  Qt2[,4]<-c(0,0,0,0)
  C<-t(eigen(Qt2)$vectors)
  Cinv<-solve(C)
  nu2<-eigen(Qt2)$values
  
  Qt3<-matrix(ncol=4,nrow=4)
  Qt3[,1]<-c(-1/a,0,0,1/a)
  Qt3[,2]<-c(0,-1/a,0,1/a)
  Qt3[,3]<-c(0,0,-1/a,1/a)
  Qt3[,4]<-c(0,0,0,0)
  D<-t(eigen(Qt3)$vectors)
  Dinv<-solve(D)
  nu3<-eigen(Qt3)$values
  
  EgX.1<-(theta*R1/(-nu+theta*R1))^X1*(-nu/(-nu+theta*R1))
  EgtauX.1<-exp(-theta*R1*tau1)*(theta*R1/(-nu+theta*R1))^X1*(-nu/(-nu+theta*R1))*exp(tau1*(-nu+theta*R1))*ppois(X1,tau1*(-nu+theta*R1))
  Egtau1Y.1<-exp(-theta*R1*tau1)*(theta*R1/(-nu2+theta*R1))^X1*(-nu2/(-nu2+theta*R1))*exp(tau1*(-nu2+theta*R1))*ppois(X1,tau1*(-nu2+theta*R1))
  Egtau0Y.1<-exp(-theta*R1*tau0)*(theta*R1/(-nu2+theta*R1))^X1*(-nu2/(-nu2+theta*R1))*exp(tau0*(-nu2+theta*R1))*ppois(X1,tau0*(-nu2+theta*R1))
  Egtau0Z.1<-exp(-theta*R1*tau0)*(theta*R1/(-nu3+theta*R1))^X1*(-nu3/(-nu3+theta*R1))*exp(tau0*(-nu3+theta*R1))*ppois(X1,tau0*(-nu3+theta*R1))
  
  EgX.2<-(theta*R2/(-nu+theta*R2))^X2*(-nu/(-nu+theta*R2))
  EgtauX.2<-exp(-theta*R2*tau1)*(theta*R2/(-nu+theta*R2))^X2*(-nu/(-nu+theta*R2))*exp(tau1*(-nu+theta*R2))*ppois(X2,tau1*(-nu+theta*R2))
  Egtau1Y.2<-exp(-theta*R2*tau1)*(theta*R2/(-nu2+theta*R2))^X2*(-nu2/(-nu2+theta*R2))*exp(tau1*(-nu2+theta*R2))*ppois(X2,tau1*(-nu2+theta*R2))
  Egtau0Y.2<-exp(-theta*R2*tau0)*(theta*R2/(-nu2+theta*R2))^X2*(-nu2/(-nu2+theta*R2))*exp(tau0*(-nu2+theta*R2))*ppois(X2,tau0*(-nu2+theta*R2))
  Egtau0Z.2<-exp(-theta*R2*tau0)*(theta*R2/(-nu3+theta*R2))^X2*(-nu3/(-nu3+theta*R2))*exp(tau0*(-nu3+theta*R2))*ppois(X2,tau0*(-nu3+theta*R2))
  
  EgX.3<-(theta*R3/(-nu+theta*R3))^X3*(-nu/(-nu+theta*R3))
  EgtauX.3<-exp(-theta*R3*tau1)*(theta*R3/(-nu+theta*R3))^X3*(-nu/(-nu+theta*R3))*exp(tau1*(-nu+theta*R3))*ppois(X3,tau1*(-nu+theta*R3))
  Egtau1Y.3<-exp(-theta*R3*tau1)*(theta*R3/(-nu2+theta*R3))^X3*(-nu2/(-nu2+theta*R3))*exp(tau1*(-nu2+theta*R3))*ppois(X3,tau1*(-nu2+theta*R3))
  Egtau0Y.3<-exp(-theta*R3*tau0)*(theta*R3/(-nu2+theta*R3))^X3*(-nu2/(-nu2+theta*R3))*exp(tau0*(-nu2+theta*R3))*ppois(X3,tau0*(-nu2+theta*R3))
  Egtau0Z.3<-exp(-theta*R3*tau0)*(theta*R3/(-nu3+theta*R3))^X3*(-nu3/(-nu3+theta*R3))*exp(tau0*(-nu3+theta*R3))*ppois(X3,tau0*(-nu3+theta*R3))
  
  pij1<-Ginv%*%diag(exp(nu*tau1))%*%G
  pij2<-Cinv%*%diag(exp(nu2*(tau0-tau1)))%*%C
  
  
  GG<--Ginv%*%diag(G[,4])
  CC<--Cinv%*%diag(C[,4])
  DD<--Dinv%*%diag(D[,4])
  
  loglike1<-sum(log(colSums(GG[1,]*(EgX.1-EgtauX.1*exp(nu*tau1)))+
                      colSums(pij1[1,]*(CC%*%(Egtau1Y.1-Egtau0Y.1*exp(nu2*(tau0-tau1)))))+
                      colSums(pij1[1,]*(pij2%*%(DD%*%Egtau0Z.1)))))
  
  loglike2<-sum(log(colSums(GG[2,]*(EgX.2-EgtauX.2*exp(nu*tau1)))+
                      colSums(pij1[2,]*(CC%*%(Egtau1Y.2-Egtau0Y.2*exp(nu2*(tau0-tau1)))))+
                      colSums(pij1[2,]*(pij2%*%(DD%*%Egtau0Z.2)))))
  
  loglike3<-sum(log(colSums(GG[3,]*(EgX.3-EgtauX.3*exp(nu*tau1)))+
                      colSums(pij1[3,]*(CC%*%(Egtau1Y.3-Egtau0Y.3*exp(nu2*(tau0-tau1)))))+
                      colSums(pij1[3,]*(pij2%*%(DD%*%Egtau0Z.3)))))
  
  -(loglike1+loglike2+loglike3)
  
}


# secondary contact model (ancestral migration rates fixed at zero)

negll.GIM.2<-function(params,model){
  
  theta<-params[2]
  a<-params[1]/theta
  b<-params[3]/theta
  c<-c(params[4]/theta,params[5]/theta)
  tau1<-params[6]/theta
  tau0<-tau1+params[7]/theta
  M<-vector(length=2,mode="numeric")
  M[1]<-0
  M[2]<-0
  Mc<-vector(length=2,mode="numeric")
  Mc[1]<-params[8]
  Mc[2]<-params[9]
  
  R1<-matrix(rep(r1,4),nrow=4,byrow=TRUE)
  R2<-matrix(rep(r2,4),nrow=4,byrow=TRUE)
  R3<-matrix(rep(r3,4),nrow=4,byrow=TRUE)
  
  X1<-matrix(rep(x1,4),nrow=4,byrow=TRUE)
  X2<-matrix(rep(x2,4),nrow=4,byrow=TRUE)
  X3<-matrix(rep(x3,4),nrow=4,byrow=TRUE)
  
  Qt<-matrix(ncol=4,nrow=4)
  Qt[,1]<-c(-(1/c[1]+Mc[1]),0,Mc[1],1/c[1])
  Qt[,2]<-c(0,-(1/c[2]+Mc[2]),Mc[2],1/c[2])
  Qt[,3]<-c(Mc[2]/2,Mc[1]/2,-(Mc[1]+Mc[2])/2,0)
  Qt[,4]<-c(0,0,0,0)
  G<-t(eigen(Qt)$vectors)
  Ginv<-solve(G)
  nu<-eigen(Qt)$values
  
  Qt2<-matrix(ncol=4,nrow=4)
  Qt2[,1]<-c(-(1+M[1]),0,M[1],1)
  Qt2[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
  Qt2[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
  Qt2[,4]<-c(0,0,0,0)
  C<-t(eigen(Qt2)$vectors)
  Cinv<-solve(C)
  nu2<-eigen(Qt2)$values
  
  Qt3<-matrix(ncol=4,nrow=4)
  Qt3[,1]<-c(-1/a,0,0,1/a)
  Qt3[,2]<-c(0,-1/a,0,1/a)
  Qt3[,3]<-c(0,0,-1/a,1/a)
  Qt3[,4]<-c(0,0,0,0)
  D<-t(eigen(Qt3)$vectors)
  Dinv<-solve(D)
  nu3<-eigen(Qt3)$values
  
  EgX.1<-(theta*R1/(-nu+theta*R1))^X1*(-nu/(-nu+theta*R1))
  EgtauX.1<-exp(-theta*R1*tau1)*(theta*R1/(-nu+theta*R1))^X1*(-nu/(-nu+theta*R1))*exp(tau1*(-nu+theta*R1))*ppois(X1,tau1*(-nu+theta*R1))
  Egtau1Y.1<-exp(-theta*R1*tau1)*(theta*R1/(-nu2+theta*R1))^X1*(-nu2/(-nu2+theta*R1))*exp(tau1*(-nu2+theta*R1))*ppois(X1,tau1*(-nu2+theta*R1))
  Egtau0Y.1<-exp(-theta*R1*tau0)*(theta*R1/(-nu2+theta*R1))^X1*(-nu2/(-nu2+theta*R1))*exp(tau0*(-nu2+theta*R1))*ppois(X1,tau0*(-nu2+theta*R1))
  Egtau0Z.1<-exp(-theta*R1*tau0)*(theta*R1/(-nu3+theta*R1))^X1*(-nu3/(-nu3+theta*R1))*exp(tau0*(-nu3+theta*R1))*ppois(X1,tau0*(-nu3+theta*R1))
  
  EgX.2<-(theta*R2/(-nu+theta*R2))^X2*(-nu/(-nu+theta*R2))
  EgtauX.2<-exp(-theta*R2*tau1)*(theta*R2/(-nu+theta*R2))^X2*(-nu/(-nu+theta*R2))*exp(tau1*(-nu+theta*R2))*ppois(X2,tau1*(-nu+theta*R2))
  Egtau1Y.2<-exp(-theta*R2*tau1)*(theta*R2/(-nu2+theta*R2))^X2*(-nu2/(-nu2+theta*R2))*exp(tau1*(-nu2+theta*R2))*ppois(X2,tau1*(-nu2+theta*R2))
  Egtau0Y.2<-exp(-theta*R2*tau0)*(theta*R2/(-nu2+theta*R2))^X2*(-nu2/(-nu2+theta*R2))*exp(tau0*(-nu2+theta*R2))*ppois(X2,tau0*(-nu2+theta*R2))
  Egtau0Z.2<-exp(-theta*R2*tau0)*(theta*R2/(-nu3+theta*R2))^X2*(-nu3/(-nu3+theta*R2))*exp(tau0*(-nu3+theta*R2))*ppois(X2,tau0*(-nu3+theta*R2))
  
  EgX.3<-(theta*R3/(-nu+theta*R3))^X3*(-nu/(-nu+theta*R3))
  EgtauX.3<-exp(-theta*R3*tau1)*(theta*R3/(-nu+theta*R3))^X3*(-nu/(-nu+theta*R3))*exp(tau1*(-nu+theta*R3))*ppois(X3,tau1*(-nu+theta*R3))
  Egtau1Y.3<-exp(-theta*R3*tau1)*(theta*R3/(-nu2+theta*R3))^X3*(-nu2/(-nu2+theta*R3))*exp(tau1*(-nu2+theta*R3))*ppois(X3,tau1*(-nu2+theta*R3))
  Egtau0Y.3<-exp(-theta*R3*tau0)*(theta*R3/(-nu2+theta*R3))^X3*(-nu2/(-nu2+theta*R3))*exp(tau0*(-nu2+theta*R3))*ppois(X3,tau0*(-nu2+theta*R3))
  Egtau0Z.3<-exp(-theta*R3*tau0)*(theta*R3/(-nu3+theta*R3))^X3*(-nu3/(-nu3+theta*R3))*exp(tau0*(-nu3+theta*R3))*ppois(X3,tau0*(-nu3+theta*R3))
  
  pij1<-Ginv%*%diag(exp(nu*tau1))%*%G
  pij2<-Cinv%*%diag(exp(nu2*(tau0-tau1)))%*%C
  
  
  GG<--Ginv%*%diag(G[,4])
  CC<--Cinv%*%diag(C[,4])
  DD<--Dinv%*%diag(D[,4])
  
  loglike1<-sum(log(colSums(GG[1,]*(EgX.1-EgtauX.1*exp(nu*tau1)))+
                      colSums(pij1[1,]*(CC%*%(Egtau1Y.1-Egtau0Y.1*exp(nu2*(tau0-tau1)))))+
                      colSums(pij1[1,]*(pij2%*%(DD%*%Egtau0Z.1)))))
  
  loglike2<-sum(log(colSums(GG[2,]*(EgX.2-EgtauX.2*exp(nu*tau1)))+
                      colSums(pij1[2,]*(CC%*%(Egtau1Y.2-Egtau0Y.2*exp(nu2*(tau0-tau1)))))+
                      colSums(pij1[2,]*(pij2%*%(DD%*%Egtau0Z.2)))))
  
  loglike3<-sum(log(colSums(GG[3,]*(EgX.3-EgtauX.3*exp(nu*tau1)))+
                      colSums(pij1[3,]*(CC%*%(Egtau1Y.3-Egtau0Y.3*exp(nu2*(tau0-tau1)))))+
                      colSums(pij1[3,]*(pij2%*%(DD%*%Egtau0Z.3)))))
  
  -(loglike1+loglike2+loglike3)
  
}


# GIM model with symmetric migration rates in both periods

negll.GIM.3<-function(params,model){
  
  theta<-params[2]
  a<-params[1]/theta
  b<-params[3]/theta
  c<-c(params[4]/theta,params[5]/theta)
  tau1<-params[6]/theta
  tau0<-tau1+params[7]/theta
  M<-vector(length=2,mode="numeric")
  M[1]<-params[8]
  M[2]<-params[8]
  Mc<-vector(length=2,mode="numeric")
  Mc[1]<-params[9]
  Mc[2]<-params[9]
  
  R1<-matrix(rep(r1,4),nrow=4,byrow=TRUE)
  R2<-matrix(rep(r2,4),nrow=4,byrow=TRUE)
  R3<-matrix(rep(r3,4),nrow=4,byrow=TRUE)
  
  X1<-matrix(rep(x1,4),nrow=4,byrow=TRUE)
  X2<-matrix(rep(x2,4),nrow=4,byrow=TRUE)
  X3<-matrix(rep(x3,4),nrow=4,byrow=TRUE)
  
  Qt<-matrix(ncol=4,nrow=4)
  Qt[,1]<-c(-(1/c[1]+Mc[1]),0,Mc[1],1/c[1])
  Qt[,2]<-c(0,-(1/c[2]+Mc[2]),Mc[2],1/c[2])
  Qt[,3]<-c(Mc[2]/2,Mc[1]/2,-(Mc[1]+Mc[2])/2,0)
  Qt[,4]<-c(0,0,0,0)
  G<-t(eigen(Qt)$vectors)
  Ginv<-solve(G)
  nu<-eigen(Qt)$values
  
  Qt2<-matrix(ncol=4,nrow=4)
  Qt2[,1]<-c(-(1+M[1]),0,M[1],1)
  Qt2[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
  Qt2[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
  Qt2[,4]<-c(0,0,0,0)
  C<-t(eigen(Qt2)$vectors)
  Cinv<-solve(C)
  nu2<-eigen(Qt2)$values
  
  Qt3<-matrix(ncol=4,nrow=4)
  Qt3[,1]<-c(-1/a,0,0,1/a)
  Qt3[,2]<-c(0,-1/a,0,1/a)
  Qt3[,3]<-c(0,0,-1/a,1/a)
  Qt3[,4]<-c(0,0,0,0)
  D<-t(eigen(Qt3)$vectors)
  Dinv<-solve(D)
  nu3<-eigen(Qt3)$values
  
  EgX.1<-(theta*R1/(-nu+theta*R1))^X1*(-nu/(-nu+theta*R1))
  EgtauX.1<-exp(-theta*R1*tau1)*(theta*R1/(-nu+theta*R1))^X1*(-nu/(-nu+theta*R1))*exp(tau1*(-nu+theta*R1))*ppois(X1,tau1*(-nu+theta*R1))
  Egtau1Y.1<-exp(-theta*R1*tau1)*(theta*R1/(-nu2+theta*R1))^X1*(-nu2/(-nu2+theta*R1))*exp(tau1*(-nu2+theta*R1))*ppois(X1,tau1*(-nu2+theta*R1))
  Egtau0Y.1<-exp(-theta*R1*tau0)*(theta*R1/(-nu2+theta*R1))^X1*(-nu2/(-nu2+theta*R1))*exp(tau0*(-nu2+theta*R1))*ppois(X1,tau0*(-nu2+theta*R1))
  Egtau0Z.1<-exp(-theta*R1*tau0)*(theta*R1/(-nu3+theta*R1))^X1*(-nu3/(-nu3+theta*R1))*exp(tau0*(-nu3+theta*R1))*ppois(X1,tau0*(-nu3+theta*R1))
  
  EgX.2<-(theta*R2/(-nu+theta*R2))^X2*(-nu/(-nu+theta*R2))
  EgtauX.2<-exp(-theta*R2*tau1)*(theta*R2/(-nu+theta*R2))^X2*(-nu/(-nu+theta*R2))*exp(tau1*(-nu+theta*R2))*ppois(X2,tau1*(-nu+theta*R2))
  Egtau1Y.2<-exp(-theta*R2*tau1)*(theta*R2/(-nu2+theta*R2))^X2*(-nu2/(-nu2+theta*R2))*exp(tau1*(-nu2+theta*R2))*ppois(X2,tau1*(-nu2+theta*R2))
  Egtau0Y.2<-exp(-theta*R2*tau0)*(theta*R2/(-nu2+theta*R2))^X2*(-nu2/(-nu2+theta*R2))*exp(tau0*(-nu2+theta*R2))*ppois(X2,tau0*(-nu2+theta*R2))
  Egtau0Z.2<-exp(-theta*R2*tau0)*(theta*R2/(-nu3+theta*R2))^X2*(-nu3/(-nu3+theta*R2))*exp(tau0*(-nu3+theta*R2))*ppois(X2,tau0*(-nu3+theta*R2))
  
  EgX.3<-(theta*R3/(-nu+theta*R3))^X3*(-nu/(-nu+theta*R3))
  EgtauX.3<-exp(-theta*R3*tau1)*(theta*R3/(-nu+theta*R3))^X3*(-nu/(-nu+theta*R3))*exp(tau1*(-nu+theta*R3))*ppois(X3,tau1*(-nu+theta*R3))
  Egtau1Y.3<-exp(-theta*R3*tau1)*(theta*R3/(-nu2+theta*R3))^X3*(-nu2/(-nu2+theta*R3))*exp(tau1*(-nu2+theta*R3))*ppois(X3,tau1*(-nu2+theta*R3))
  Egtau0Y.3<-exp(-theta*R3*tau0)*(theta*R3/(-nu2+theta*R3))^X3*(-nu2/(-nu2+theta*R3))*exp(tau0*(-nu2+theta*R3))*ppois(X3,tau0*(-nu2+theta*R3))
  Egtau0Z.3<-exp(-theta*R3*tau0)*(theta*R3/(-nu3+theta*R3))^X3*(-nu3/(-nu3+theta*R3))*exp(tau0*(-nu3+theta*R3))*ppois(X3,tau0*(-nu3+theta*R3))
  
  pij1<-Ginv%*%diag(exp(nu*tau1))%*%G
  pij2<-Cinv%*%diag(exp(nu2*(tau0-tau1)))%*%C
  
  
  GG<--Ginv%*%diag(G[,4])
  CC<--Cinv%*%diag(C[,4])
  DD<--Dinv%*%diag(D[,4])
  
  loglike1<-sum(log(colSums(GG[1,]*(EgX.1-EgtauX.1*exp(nu*tau1)))+
                      colSums(pij1[1,]*(CC%*%(Egtau1Y.1-Egtau0Y.1*exp(nu2*(tau0-tau1)))))+
                      colSums(pij1[1,]*(pij2%*%(DD%*%Egtau0Z.1)))))
  
  loglike2<-sum(log(colSums(GG[2,]*(EgX.2-EgtauX.2*exp(nu*tau1)))+
                      colSums(pij1[2,]*(CC%*%(Egtau1Y.2-Egtau0Y.2*exp(nu2*(tau0-tau1)))))+
                      colSums(pij1[2,]*(pij2%*%(DD%*%Egtau0Z.2)))))
  
  loglike3<-sum(log(colSums(GG[3,]*(EgX.3-EgtauX.3*exp(nu*tau1)))+
                      colSums(pij1[3,]*(CC%*%(Egtau1Y.3-Egtau0Y.3*exp(nu2*(tau0-tau1)))))+
                      colSums(pij1[3,]*(pij2%*%(DD%*%Egtau0Z.3)))))
  
  -(loglike1+loglike2+loglike3)
  
}


# secondary contact model (ancestral migration rates fixed at zero)
# with symmetric migration

negll.GIM.4<-function(params,model){
  
  theta<-params[2]
  a<-params[1]/theta
  b<-params[3]/theta
  c<-c(params[4]/theta,params[5]/theta)
  tau1<-params[6]/theta
  tau0<-tau1+params[7]/theta
  M<-vector(length=2,mode="numeric")
  M[1]<-0
  M[2]<-0
  Mc<-vector(length=2,mode="numeric")
  Mc[1]<-params[8]
  Mc[2]<-params[8]
  
  R1<-matrix(rep(r1,4),nrow=4,byrow=TRUE)
  R2<-matrix(rep(r2,4),nrow=4,byrow=TRUE)
  R3<-matrix(rep(r3,4),nrow=4,byrow=TRUE)
  
  X1<-matrix(rep(x1,4),nrow=4,byrow=TRUE)
  X2<-matrix(rep(x2,4),nrow=4,byrow=TRUE)
  X3<-matrix(rep(x3,4),nrow=4,byrow=TRUE)
  
  Qt<-matrix(ncol=4,nrow=4)
  Qt[,1]<-c(-(1/c[1]+Mc[1]),0,Mc[1],1/c[1])
  Qt[,2]<-c(0,-(1/c[2]+Mc[2]),Mc[2],1/c[2])
  Qt[,3]<-c(Mc[2]/2,Mc[1]/2,-(Mc[1]+Mc[2])/2,0)
  Qt[,4]<-c(0,0,0,0)
  G<-t(eigen(Qt)$vectors)
  Ginv<-solve(G)
  nu<-eigen(Qt)$values
  
  Qt2<-matrix(ncol=4,nrow=4)
  Qt2[,1]<-c(-(1+M[1]),0,M[1],1)
  Qt2[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
  Qt2[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
  Qt2[,4]<-c(0,0,0,0)
  C<-t(eigen(Qt2)$vectors)
  Cinv<-solve(C)
  nu2<-eigen(Qt2)$values
  
  Qt3<-matrix(ncol=4,nrow=4)
  Qt3[,1]<-c(-1/a,0,0,1/a)
  Qt3[,2]<-c(0,-1/a,0,1/a)
  Qt3[,3]<-c(0,0,-1/a,1/a)
  Qt3[,4]<-c(0,0,0,0)
  D<-t(eigen(Qt3)$vectors)
  Dinv<-solve(D)
  nu3<-eigen(Qt3)$values
  
  EgX.1<-(theta*R1/(-nu+theta*R1))^X1*(-nu/(-nu+theta*R1))
  EgtauX.1<-exp(-theta*R1*tau1)*(theta*R1/(-nu+theta*R1))^X1*(-nu/(-nu+theta*R1))*exp(tau1*(-nu+theta*R1))*ppois(X1,tau1*(-nu+theta*R1))
  Egtau1Y.1<-exp(-theta*R1*tau1)*(theta*R1/(-nu2+theta*R1))^X1*(-nu2/(-nu2+theta*R1))*exp(tau1*(-nu2+theta*R1))*ppois(X1,tau1*(-nu2+theta*R1))
  Egtau0Y.1<-exp(-theta*R1*tau0)*(theta*R1/(-nu2+theta*R1))^X1*(-nu2/(-nu2+theta*R1))*exp(tau0*(-nu2+theta*R1))*ppois(X1,tau0*(-nu2+theta*R1))
  Egtau0Z.1<-exp(-theta*R1*tau0)*(theta*R1/(-nu3+theta*R1))^X1*(-nu3/(-nu3+theta*R1))*exp(tau0*(-nu3+theta*R1))*ppois(X1,tau0*(-nu3+theta*R1))
  
  EgX.2<-(theta*R2/(-nu+theta*R2))^X2*(-nu/(-nu+theta*R2))
  EgtauX.2<-exp(-theta*R2*tau1)*(theta*R2/(-nu+theta*R2))^X2*(-nu/(-nu+theta*R2))*exp(tau1*(-nu+theta*R2))*ppois(X2,tau1*(-nu+theta*R2))
  Egtau1Y.2<-exp(-theta*R2*tau1)*(theta*R2/(-nu2+theta*R2))^X2*(-nu2/(-nu2+theta*R2))*exp(tau1*(-nu2+theta*R2))*ppois(X2,tau1*(-nu2+theta*R2))
  Egtau0Y.2<-exp(-theta*R2*tau0)*(theta*R2/(-nu2+theta*R2))^X2*(-nu2/(-nu2+theta*R2))*exp(tau0*(-nu2+theta*R2))*ppois(X2,tau0*(-nu2+theta*R2))
  Egtau0Z.2<-exp(-theta*R2*tau0)*(theta*R2/(-nu3+theta*R2))^X2*(-nu3/(-nu3+theta*R2))*exp(tau0*(-nu3+theta*R2))*ppois(X2,tau0*(-nu3+theta*R2))
  
  EgX.3<-(theta*R3/(-nu+theta*R3))^X3*(-nu/(-nu+theta*R3))
  EgtauX.3<-exp(-theta*R3*tau1)*(theta*R3/(-nu+theta*R3))^X3*(-nu/(-nu+theta*R3))*exp(tau1*(-nu+theta*R3))*ppois(X3,tau1*(-nu+theta*R3))
  Egtau1Y.3<-exp(-theta*R3*tau1)*(theta*R3/(-nu2+theta*R3))^X3*(-nu2/(-nu2+theta*R3))*exp(tau1*(-nu2+theta*R3))*ppois(X3,tau1*(-nu2+theta*R3))
  Egtau0Y.3<-exp(-theta*R3*tau0)*(theta*R3/(-nu2+theta*R3))^X3*(-nu2/(-nu2+theta*R3))*exp(tau0*(-nu2+theta*R3))*ppois(X3,tau0*(-nu2+theta*R3))
  Egtau0Z.3<-exp(-theta*R3*tau0)*(theta*R3/(-nu3+theta*R3))^X3*(-nu3/(-nu3+theta*R3))*exp(tau0*(-nu3+theta*R3))*ppois(X3,tau0*(-nu3+theta*R3))
  
  pij1<-Ginv%*%diag(exp(nu*tau1))%*%G
  pij2<-Cinv%*%diag(exp(nu2*(tau0-tau1)))%*%C
  
  
  GG<--Ginv%*%diag(G[,4])
  CC<--Cinv%*%diag(C[,4])
  DD<--Dinv%*%diag(D[,4])
  
  loglike1<-sum(log(colSums(GG[1,]*(EgX.1-EgtauX.1*exp(nu*tau1)))+
                      colSums(pij1[1,]*(CC%*%(Egtau1Y.1-Egtau0Y.1*exp(nu2*(tau0-tau1)))))+
                      colSums(pij1[1,]*(pij2%*%(DD%*%Egtau0Z.1)))))
  
  loglike2<-sum(log(colSums(GG[2,]*(EgX.2-EgtauX.2*exp(nu*tau1)))+
                      colSums(pij1[2,]*(CC%*%(Egtau1Y.2-Egtau0Y.2*exp(nu2*(tau0-tau1)))))+
                      colSums(pij1[2,]*(pij2%*%(DD%*%Egtau0Z.2)))))
  
  loglike3<-sum(log(colSums(GG[3,]*(EgX.3-EgtauX.3*exp(nu*tau1)))+
                      colSums(pij1[3,]*(CC%*%(Egtau1Y.3-Egtau0Y.3*exp(nu2*(tau0-tau1)))))+
                      colSums(pij1[3,]*(pij2%*%(DD%*%Egtau0Z.3)))))
  
  -(loglike1+loglike2+loglike3)
  
}



### MODELS IN FIGURE 4, IIM paper


# "ISO" 

negll.iso.1<-function(params){
  
  a<-params[1]/params[2]
  theta<-params[2]
  b<-params[3]/theta
  tau<-params[4]/theta
  
  
  R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
  R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
  R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
  
  X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
  X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
  X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
  
  loglike1<-log(((theta*R1[1,])^X1[1,])/    		
                  ((1+theta*R1[1,])^(X1[1,]+1))*							
                  (1-ppois(X1[1,],tau*(1+theta*R1[1,])))+  	
                  exp(-tau*(1+theta*R1[1,]))*
                  (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp(tau*(1/a+theta*R1[1,]))*ppois(X1[1,],tau*(1/a+theta*R1[1,])))
  
  loglike2<-log(((b*theta*R2[1,])^X2[1,])/				
                  ((1+b*theta*R2[1,])^(X2[1,]+1))*							
                  (1-ppois(X2[1,],tau*(1/b+theta*R2[1,])))+  	
                  exp(-tau*(1/b+theta*R2[1,]))*
                  (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp(tau*(1/a+theta*R2[1,]))*ppois(X2[1,],tau*(1/a+theta*R2[1,])))
  
  loglike3<-log(exp(-tau*theta*R3[1,])*			
                  (a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau*(1/a+theta*R3[1,]))*ppois(X3[1,],tau*(1/a+theta*R3[1,])))
  
  -sum(c(loglike1,loglike2,loglike3))
}


# "IM1" 

negll.IM.1<-function(params){
  
  a<-params[1]/params[2]
  theta<-params[2]
  b<-params[3]/theta
  tau<-params[4]/theta
  M<-vector(length=2,mode="numeric")
  M[1]<-params[5]
  M[2]<-params[6]
  
  
  
  
  if (M[1]>0 & M[2]>0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    Qt<-matrix(ncol=4,nrow=4)
    Qt[,1]<-c(-(1+M[1]),0,M[1],1)
    Qt[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
    Qt[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
    Qt[,4]<-c(0,0,0,0)
    X<-eigen(Qt)$vectors
    nu<-eigen(Qt)$values[-4]
    P.0<-diag(4) 	#boundary conditions for all starting states form an identity matrix
    C<-solve(X) #(%*%P.0)
    A<--{X[4,-4]*C[-4,]} #minus sign precedes because of differentiation wrt t;
    
    
    loglike1<-log(colSums(A[,1]*((-nu*(theta*R1)^X1)/
                                   ((theta*R1-nu)^(X1+1))*
                                   (1-ppois(X1,tau*(-nu+theta*R1)))+
                                   exp(-tau*(theta*R1-nu))*
                                   (a*theta*R1)^X1/
                                   (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A[,2]*((-nu*(theta*R2)^X2)/
                                   ((theta*R2-nu)^(X2+1))*
                                   (1-ppois(X2,tau*(-nu+theta*R2)))+
                                   exp(-tau*(theta*R2-nu))*
                                   (a*theta*R2)^X2/
                                   (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A[,3]*((-nu*(theta*R3)^X3)/
                                   ((theta*R3-nu)^(X3+1))*
                                   (1-ppois(X3,tau*(-nu+theta*R3)))+
                                   exp(-tau*(theta*R3-nu))*
                                   (a*theta*R3)^X3/
                                   (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
  }
  
  if (M[1]==0 & M[2]>0){
    
    R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    
    loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
                                (-nu1+theta*R1)^(X1+1)*
                                (1-ppois(X1,tau*(-nu1+theta*R1)))+
                                exp(-tau*(theta*R1-nu1))*
                                (a*theta*R1)^X1/
                                (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
                                (-nu2+theta*R2)^(X2+1)*
                                (1-ppois(X2,tau*(-nu2+theta*R2)))+
                                exp(-tau*(theta*R2-nu2))*
                                (a*theta*R2)^X2/
                                (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
                                (-nu3+theta*R3)^(X3+1)*
                                (1-ppois(X3,tau*(-nu3+theta*R3)))+
                                exp(-tau*(theta*R3-nu3))*
                                (a*theta*R3)^X3/
                                (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
  }
  
  if (M[1]>0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    A11<-(b^2*M[1]^2)/((b*M[1]-2)*(b-1+b*M[1]))
    nu11<--1/b
    A12<-4*M[1]/((2+M[1])*(2-b*M[1]))
    nu12<--M[1]/2
    A13<-1/(1+M[1])+M[1]^2/((2+M[1])*(b-1+b*M[1])*(1+M[1]))
    nu13<--(1+M[1])
    
    A1<-c(A11,A12,A13)
    nu1<-c(nu11,nu12,nu13)
    
    A2<-1
    nu2<--1/b
    
    A31<-b*M[1]/(b*M[1]-2)
    nu31<--1/b
    A32<-2/(2-b*M[1])
    nu32<--M[1]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
                                (-nu1+theta*R1)^(X1+1)*
                                (1-ppois(X1,tau*(-nu1+theta*R1)))+
                                exp(-tau*(theta*R1-nu1))*
                                (a*theta*R1)^X1/
                                (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
                                (-nu2+theta*R2)^(X2+1)*
                                (1-ppois(X2,tau*(-nu2+theta*R2)))+
                                exp(-tau*(theta*R2-nu2))*
                                (a*theta*R2)^X2/
                                (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
                                (-nu3+theta*R3)^(X3+1)*
                                (1-ppois(X3,tau*(-nu3+theta*R3)))+
                                exp(-tau*(theta*R3-nu3))*
                                (a*theta*R3)^X3/
                                (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
    
  }
  
  if (M[1]==0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    loglike1<-log(((theta*R1[1,])^X1[1,])/				
                    ((1+theta*R1[1,])^(X1[1,]+1))*							
                    (1-ppois(X1[1,],tau*(1+theta*R1[1,])))+  	
                    exp(-tau*(1+theta*R1[1,]))*
                    (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp(tau*(1/a+theta*R1[1,]))*ppois(X1[1,],tau*(1/a+theta*R1[1,])))
    
    loglike2<-log(((b*theta*R2[1,])^X2[1,])/				
                    ((1+b*theta*R2[1,])^(X2[1,]+1))*							
                    (1-ppois(X2[1,],tau*(1/b+theta*R2[1,])))+  	
                    exp(-tau*(1/b+theta*R2[1,]))*
                    (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp(tau*(1/a+theta*R2[1,]))*ppois(X2[1,],tau*(1/a+theta*R2[1,])))
    
    loglike3<-log(exp(-tau*theta*R3[1,])*			
                    (a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau*(1/a+theta*R3[1,]))*ppois(X3[1,],tau*(1/a+theta*R3[1,])))
    
  }
  
  -sum(c(loglike1,loglike2,loglike3),na.rm=TRUE)
}


# "IIM1" 

negll.IIM.1<-function(params){
  
  a<-params[1]/params[2]
  theta<-params[2]
  b<-params[3]/theta
  c<-c(1,params[3]/theta)
  tau1<-params[4]/theta
  tau0<-params[5]/theta+tau1
  M<-vector(length=2,mode="numeric")
  M[1]<-params[6]
  M[2]<-params[7]
  
  
  if (M[1]>0 & M[2]>0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    
    Qt<-matrix(ncol=4,nrow=4)
    Qt[,1]<-c(-(1+M[1]),0,M[1],1)
    Qt[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
    Qt[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
    Qt[,4]<-c(0,0,0,0)
    X<-eigen(Qt)$vectors
    nu<-eigen(Qt)$values[-4]
    P.0<-diag(4) 	#boundary conditions for all starting states form an identity matrix
    C<-solve(X) #(%*%P.0)
    A<--{X[4,-4]*C[-4,]} #minus sign precedes because of differentiation wrt t;
    
    
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A[,1]*(-nu*(R1*theta)^X1/(-nu+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu+theta*R1)*tau1)*ppois(X1,(-nu+theta*R1)*tau1)-
                                                                                                  exp(nu*(tau0-tau1)-theta*R1*tau0)*exp((-nu+theta*R1)*tau0)*ppois(X1,(-nu+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A[,1]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A[,2]*(-nu*(R2*theta)^X2/(-nu+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu+theta*R2)*tau1)*ppois(X2,(-nu+theta*R2)*tau1)-
                                                                                                  exp(nu*(tau0-tau1)-theta*R2*tau0)*exp((-nu+theta*R2)*tau0)*ppois(X2,(-nu+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A[,2]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A[,3]*(-nu*(R3*theta)^X3/(-nu+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu+theta*R3)*tau1)*ppois(X3,(-nu+theta*R3)*tau1)-
                                                                                exp(nu*(tau0-tau1)-theta*R3*tau0)*exp((-nu+theta*R3)*tau0)*ppois(X3,(-nu+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A[,3]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
  }
  
  if (M[1]==0 & M[2]>0){
    
    R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                                 exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*A1*exp(nu1*(tau0-tau1))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                                 exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                               exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
    
  }
  
  if (M[1]>0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    A11<-(b^2*M[1]^2)/((b*M[1]-2)*(b-1+b*M[1]))
    nu11<--1/b
    A12<-4*M[1]/((2+M[1])*(2-b*M[1]))
    nu12<--M[1]/2
    A13<-1/(1+M[1])+M[1]^2/((2+M[1])*(b-1+b*M[1])*(1+M[1]))
    nu13<--(1+M[1])
    
    A1<-c(A11,A12,A13)
    nu1<-c(nu11,nu12,nu13)
    
    A2<-1
    nu2<--1/b
    
    A31<-b*M[1]/(b*M[1]-2)
    nu31<--1/b
    A32<-2/(2-b*M[1])
    nu32<--M[1]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                                 exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A1*exp(nu1*(tau0-tau1)))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                                 exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                               exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
    
  }
  
  if (M[1]==0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    loglike1<-sum(log(
      (c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)+
        exp(-tau1*(1/c[1]+theta*R1[1,]))*(   (theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
                                               exp(tau1*(1+theta*R1[1,]))*ppois(X1[1,],tau1*(1+theta*R1[1,]))-(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*
                                               exp(tau1*(1/c[1]+theta*R1[1,]))*ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,]))  )
      +
        exp(-(tau1/c[1]+tau0-tau1+theta*R1[1,]*tau0))*(   (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*
                                                            exp(tau0*(1/a+theta*R1[1,]))*ppois(X1[1,],tau0*(1/a+theta*R1[1,]))-(theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
                                                            exp(tau0*(1+theta*R1[1,]))*ppois(X1[1,],tau0*(1+theta*R1[1,]))  )
    ))
    
    
    loglike2<-sum(log(
      (c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)+
        exp(-tau1*(1/c[2]+theta*R2[1,]))*(   (b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
                                               exp(tau1*(1/b+theta*R2[1,]))*ppois(X2[1,],tau1*(1/b+theta*R2[1,]))-(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*
                                               exp(tau1*(1/c[2]+theta*R2[1,]))*ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,]))  )
      +
        exp(-(tau1/c[2]+(tau0-tau1)/b+theta*R2[1,]*tau0))*(   (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*
                                                                exp(tau0*(1/a+theta*R2[1,]))*ppois(X2[1,],tau0*(1/a+theta*R2[1,]))-(b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
                                                                exp(tau0*(1/b+theta*R2[1,]))*ppois(X2[1,],tau0*(1/b+theta*R2[1,]))  )
    ))
    
    
    loglike3<-sum(log(
      exp(-tau0*theta*R3[1,])*(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau0*(1/a+theta*R3[1,]))*ppois(X3[1,],tau0*(1/a+theta*R3[1,]))
    ))
    
  }
  
  
  -(loglike1+loglike2+loglike3)
}


# "IIM2" 

negll.IIM.2<-function(params){
  
  a<-params[1]/params[2]
  theta<-params[2]
  b<-params[3]/theta
  c<-c(params[4]/theta,params[5]/theta)
  tau1<-params[6]/theta
  tau0<-params[7]/theta+tau1
  M<-vector(length=2,mode="numeric")
  M[1]<-params[8]
  M[2]<-params[9]
  
  if (M[1]>0 & M[2]>0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    
    Qt<-matrix(ncol=4,nrow=4)
    Qt[,1]<-c(-(1+M[1]),0,M[1],1)
    Qt[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
    Qt[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
    Qt[,4]<-c(0,0,0,0)
    X<-eigen(Qt)$vectors
    nu<-eigen(Qt)$values[-4]
    P.0<-diag(4) 	#boundary conditions for all starting states form an identity matrix
    C<-solve(X) #(%*%P.0)
    A<--{X[4,-4]*C[-4,]} #minus sign precedes because of differentiation wrt t;
    
    
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A[,1]*(-nu*(R1*theta)^X1/(-nu+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu+theta*R1)*tau1)*ppois(X1,(-nu+theta*R1)*tau1)-
                                                                                                  exp(nu*(tau0-tau1)-theta*R1*tau0)*exp((-nu+theta*R1)*tau0)*ppois(X1,(-nu+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A[,1]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A[,2]*(-nu*(R2*theta)^X2/(-nu+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu+theta*R2)*tau1)*ppois(X2,(-nu+theta*R2)*tau1)-
                                                                                                  exp(nu*(tau0-tau1)-theta*R2*tau0)*exp((-nu+theta*R2)*tau0)*ppois(X2,(-nu+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A[,2]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A[,3]*(-nu*(R3*theta)^X3/(-nu+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu+theta*R3)*tau1)*ppois(X3,(-nu+theta*R3)*tau1)-
                                                                                exp(nu*(tau0-tau1)-theta*R3*tau0)*exp((-nu+theta*R3)*tau0)*ppois(X3,(-nu+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A[,3]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
  }
  
  if (M[1]==0 & M[2]>0){
    
    R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                                 exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*A1*exp(nu1*(tau0-tau1))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                                 exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                               exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
    
  }
  
  if (M[1]>0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    A11<-(b^2*M[1]^2)/((b*M[1]-2)*(b-1+b*M[1]))
    nu11<--1/b
    A12<-4*M[1]/((2+M[1])*(2-b*M[1]))
    nu12<--M[1]/2
    A13<-1/(1+M[1])+M[1]^2/((2+M[1])*(b-1+b*M[1])*(1+M[1]))
    nu13<--(1+M[1])
    
    A1<-c(A11,A12,A13)
    nu1<-c(nu11,nu12,nu13)
    
    A2<-1
    nu2<--1/b
    
    A31<-b*M[1]/(b*M[1]-2)
    nu31<--1/b
    A32<-2/(2-b*M[1])
    nu32<--M[1]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                                 exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A1*exp(nu1*(tau0-tau1)))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                                 exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                               exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
    
  }
  
  if (M[1]==0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    loglike1<-sum(log(
      (c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)+
        exp(-tau1*(1/c[1]+theta*R1[1,]))*(   (theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
                                               exp(tau1*(1+theta*R1[1,]))*ppois(X1[1,],tau1*(1+theta*R1[1,]))-(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*
                                               exp(tau1*(1/c[1]+theta*R1[1,]))*ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,]))  )
      +
        exp(-(tau1/c[1]+tau0-tau1+theta*R1[1,]*tau0))*(   (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*
                                                            exp(tau0*(1/a+theta*R1[1,]))*ppois(X1[1,],tau0*(1/a+theta*R1[1,]))-(theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
                                                            exp(tau0*(1+theta*R1[1,]))*ppois(X1[1,],tau0*(1+theta*R1[1,]))  )
    ))
    
    
    loglike2<-sum(log(
      (c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)+
        exp(-tau1*(1/c[2]+theta*R2[1,]))*(   (b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
                                               exp(tau1*(1/b+theta*R2[1,]))*ppois(X2[1,],tau1*(1/b+theta*R2[1,]))-(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*
                                               exp(tau1*(1/c[2]+theta*R2[1,]))*ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,]))  )
      +
        exp(-(tau1/c[2]+(tau0-tau1)/b+theta*R2[1,]*tau0))*(   (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*
                                                                exp(tau0*(1/a+theta*R2[1,]))*ppois(X2[1,],tau0*(1/a+theta*R2[1,]))-(b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
                                                                exp(tau0*(1/b+theta*R2[1,]))*ppois(X2[1,],tau0*(1/b+theta*R2[1,]))  )
    ))
    
    
    loglike3<-sum(log(
      exp(-tau0*theta*R3[1,])*(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau0*(1/a+theta*R3[1,]))*ppois(X3[1,],tau0*(1/a+theta*R3[1,]))
    ))
    
  }
  
  
  -(loglike1+loglike2+loglike3)
}


# "IIM3"

negll.IIM.3<-function(params){
  
  a<-params[1]/params[2]
  theta<-params[2]
  b<-params[3]/theta
  c<-c(params[4]/theta,params[5]/theta)
  tau1<-params[6]/theta
  tau0<-params[7]/theta+tau1
  M2<-params[8]
  
  
  X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
  X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
  X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
  
  R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
  R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
  R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
  
  A1<-1
  nu1<--1
  
  A21<-(b*M2^2)/((-2 + M2)*(1-b+b*M2))
  nu21<--1
  A22<-4*b*M2/((2-M2)*(2+b*M2))
  nu22<--M2/2
  A23<-(1/b)/(1/b+M2)+b^2*M2^2/((2+b*M2)*(1-b+b*M2)*(1/b+M2))
  nu23<--(1/b+M2)
  
  A2<-c(A21,A22,A23)
  nu2<-c(nu21,nu22,nu23)
  
  A31<-M2/(M2-2)
  nu31<--1
  A32<-2/(2-M2)
  nu32<--M2/2
  
  A3<-c(A31,A32)
  nu3<-c(nu31,nu32)
  
  loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                      exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                               exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                      exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*A1*exp(nu1*(tau0-tau1))),na.rm=TRUE)
  
  loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                      exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                               exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                      exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
  
  loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                             exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                      exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
  
  
  
  -(loglike1+loglike2+loglike3)
}





### OTHER MODELS

# "IM1" model with M1=0 ("IM2")

negll.IM.2<-function(params){
  
  a<-params[1]/params[2]
  theta<-params[2]
  b<-params[3]/theta
  tau<-params[4]/theta
  M2<-params[5]
  
  
  X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
  X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
  X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
  
  R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
  R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
  R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
  
  A1<-1
  nu1<--1
  
  A21<-(b*M2^2)/((-2 + M2)*(1-b+b*M2))
  nu21<--1
  A22<-4*b*M2/((2-M2)*(2+b*M2))
  nu22<--M2/2
  A23<-(1/b)/(1/b+M2)+b^2*M2^2/((2+b*M2)*(1-b+b*M2)*(1/b+M2))
  nu23<--(1/b+M2)
  
  A2<-c(A21,A22,A23)
  nu2<-c(nu21,nu22,nu23)
  
  A31<-M2/(M2-2)
  nu31<--1
  A32<-2/(2-M2)
  nu32<--M2/2
  
  A3<-c(A31,A32)
  nu3<-c(nu31,nu32)
  
  
  loglike1<-log(t(A1)%*%((-nu1*(theta*R1)^X1)/
                           ((theta*R1-nu1)^(X1+1))*
                           (1-ppois(X1,tau*(-nu1+theta*R1)))+
                           exp(-tau*(theta*R1-nu1))*
                           (a*theta*R1)^X1/
                           (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1))))
  
  loglike2<-log(t(A2)%*%((-nu2*(theta*R2)^X2)/
                           ((theta*R2-nu2)^(X2+1))*
                           (1-ppois(X2,tau*(-nu2+theta*R2)))+
                           exp(-tau*(theta*R2-nu2))*
                           (a*theta*R2)^X2/
                           (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2))))
  
  loglike3<-log(t(A3)%*%((-nu3*(theta*R3)^X3)/
                           ((theta*R3-nu3)^(X3+1))*
                           (1-ppois(X3,tau*(-nu3+theta*R3)))+
                           exp(-tau*(theta*R3-nu3))*
                           (a*theta*R3)^X3/
                           (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3))))
  
  -sum(c(loglike1,loglike2,loglike3),na.rm=TRUE)
}


# "IM1" model with M2=0 ("IM3")

negll.IM.3<-function(params){
  
  a<-params[1]/params[2]
  theta<-params[2]
  b<-params[3]/theta
  tau<-params[4]/theta
  M1<-params[5]
  
  
  X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
  X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
  X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
  
  R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
  R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
  R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
  
  
  
  A11<-(b^2*M1^2)/((b*M1-2)*(b-1+b*M1))
  nu11<--1/b
  A12<-4*M1/((2+M1)*(2-b*M1))
  nu12<--M1/2
  A13<-1/(1+M1)+M1^2/((2+M1)*(b-1+b*M1)*(1+M1))
  nu13<--(1+M1)
  
  A1<-c(A11,A12,A13)
  nu1<-c(nu11,nu12,nu13)
  
  A2<-1
  nu2<--1/b
  
  A31<-b*M1/(b*M1-2)
  nu31<--1/b
  A32<-2/(2-b*M1)
  nu32<--M1/2
  
  A3<-c(A31,A32)
  nu3<-c(nu31,nu32)
  
  
  
  loglike1<-log(t(A1)%*%((-nu1*(theta*R1)^X1)/
                           ((theta*R1-nu1)^(X1+1))*
                           (1-ppois(X1,tau*(-nu1+theta*R1)))+
                           exp(-tau*(theta*R1-nu1))*
                           (a*theta*R1)^X1/
                           (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1))))
  
  loglike2<-log(t(A2)%*%((-nu2*(theta*R2)^X2)/
                           ((theta*R2-nu2)^(X2+1))*
                           (1-ppois(X2,tau*(-nu2+theta*R2)))+
                           exp(-tau*(theta*R2-nu2))*
                           (a*theta*R2)^X2/
                           (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2))))
  
  loglike3<-log(t(A3)%*%((-nu3*(theta*R3)^X3)/
                           ((theta*R3-nu3)^(X3+1))*
                           (1-ppois(X3,tau*(-nu3+theta*R3)))+
                           exp(-tau*(theta*R3-nu3))*
                           (a*theta*R3)^X3/
                           (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3))))
  
  -sum(c(loglike1,loglike2,loglike3),na.rm=TRUE)
}


# "IM1" model with symmetric migration rates ("IM4")

negll.IM.4<-function(params){
  
  a<-params[1]/params[2]
  theta<-params[2]
  b<-params[3]/theta
  tau<-params[4]/theta
  M<-vector(length=2,mode="numeric")
  M[1]<-params[5]
  M[2]<-params[5]
  
  
  
  
  if (M[1]>0 & M[2]>0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    Qt<-matrix(ncol=4,nrow=4)
    Qt[,1]<-c(-(1+M[1]),0,M[1],1)
    Qt[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
    Qt[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
    Qt[,4]<-c(0,0,0,0)
    X<-eigen(Qt)$vectors
    nu<-eigen(Qt)$values[-4]
    P.0<-diag(4) 	#boundary conditions for all starting states form an identity matrix
    C<-solve(X) #(%*%P.0)
    A<--{X[4,-4]*C[-4,]} #minus sign precedes because of differentiation wrt t;
    
    
    loglike1<-log(colSums(A[,1]*((-nu*(theta*R1)^X1)/
                                   ((theta*R1-nu)^(X1+1))*
                                   (1-ppois(X1,tau*(-nu+theta*R1)))+
                                   exp(-tau*(theta*R1-nu))*
                                   (a*theta*R1)^X1/
                                   (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A[,2]*((-nu*(theta*R2)^X2)/
                                   ((theta*R2-nu)^(X2+1))*
                                   (1-ppois(X2,tau*(-nu+theta*R2)))+
                                   exp(-tau*(theta*R2-nu))*
                                   (a*theta*R2)^X2/
                                   (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A[,3]*((-nu*(theta*R3)^X3)/
                                   ((theta*R3-nu)^(X3+1))*
                                   (1-ppois(X3,tau*(-nu+theta*R3)))+
                                   exp(-tau*(theta*R3-nu))*
                                   (a*theta*R3)^X3/
                                   (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
  }
  
  if (M[1]==0 & M[2]>0){
    
    R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    
    loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
                                (-nu1+theta*R1)^(X1+1)*
                                (1-ppois(X1,tau*(-nu1+theta*R1)))+
                                exp(-tau*(theta*R1-nu1))*
                                (a*theta*R1)^X1/
                                (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
                                (-nu2+theta*R2)^(X2+1)*
                                (1-ppois(X2,tau*(-nu2+theta*R2)))+
                                exp(-tau*(theta*R2-nu2))*
                                (a*theta*R2)^X2/
                                (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
                                (-nu3+theta*R3)^(X3+1)*
                                (1-ppois(X3,tau*(-nu3+theta*R3)))+
                                exp(-tau*(theta*R3-nu3))*
                                (a*theta*R3)^X3/
                                (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
  }
  
  if (M[1]>0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    A11<-(b^2*M[1]^2)/((b*M[1]-2)*(b-1+b*M[1]))
    nu11<--1/b
    A12<-4*M[1]/((2+M[1])*(2-b*M[1]))
    nu12<--M[1]/2
    A13<-1/(1+M[1])+M[1]^2/((2+M[1])*(b-1+b*M[1])*(1+M[1]))
    nu13<--(1+M[1])
    
    A1<-c(A11,A12,A13)
    nu1<-c(nu11,nu12,nu13)
    
    A2<-1
    nu2<--1/b
    
    A31<-b*M[1]/(b*M[1]-2)
    nu31<--1/b
    A32<-2/(2-b*M[1])
    nu32<--M[1]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
                                (-nu1+theta*R1)^(X1+1)*
                                (1-ppois(X1,tau*(-nu1+theta*R1)))+
                                exp(-tau*(theta*R1-nu1))*
                                (a*theta*R1)^X1/
                                (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
                                (-nu2+theta*R2)^(X2+1)*
                                (1-ppois(X2,tau*(-nu2+theta*R2)))+
                                exp(-tau*(theta*R2-nu2))*
                                (a*theta*R2)^X2/
                                (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
                                (-nu3+theta*R3)^(X3+1)*
                                (1-ppois(X3,tau*(-nu3+theta*R3)))+
                                exp(-tau*(theta*R3-nu3))*
                                (a*theta*R3)^X3/
                                (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
    
  }
  
  if (M[1]==0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    loglike1<-log(((theta*R1[1,])^X1[1,])/				
                    ((1+theta*R1[1,])^(X1[1,]+1))*							
                    (1-ppois(X1[1,],tau*(1+theta*R1[1,])))+  	
                    exp(-tau*(1+theta*R1[1,]))*
                    (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp(tau*(1/a+theta*R1[1,]))*ppois(X1[1,],tau*(1/a+theta*R1[1,])))
    
    loglike2<-log(((b*theta*R2[1,])^X2[1,])/				
                    ((1+b*theta*R2[1,])^(X2[1,]+1))*							
                    (1-ppois(X2[1,],tau*(1/b+theta*R2[1,])))+  	
                    exp(-tau*(1/b+theta*R2[1,]))*
                    (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp(tau*(1/a+theta*R2[1,]))*ppois(X2[1,],tau*(1/a+theta*R2[1,])))
    
    loglike3<-log(exp(-tau*theta*R3[1,])*			
                    (a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau*(1/a+theta*R3[1,]))*ppois(X3[1,],tau*(1/a+theta*R3[1,])))
    
  }
  
  -sum(c(loglike1,loglike2,loglike3),na.rm=TRUE)
}


# "IIM2" model with M2=0  ("IIM4")

negll.IIM.4<-function(params){
  
  a<-params[1]/params[2]
  theta<-params[2]
  b<-params[3]/theta
  c<-c(params[4]/theta,params[5]/theta)
  tau1<-params[6]/theta
  tau0<-params[7]/theta+tau1
  M1<-params[8]
  
  
  X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
  X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
  X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
  
  R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
  R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
  R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
  
  
  A11<-(b^2*M1^2)/((b*M1-2)*(b-1+b*M1))
  nu11<--1/b
  A12<-4*M1/((2+M1)*(2-b*M1))
  nu12<--M1/2
  A13<-1/(1+M1)+M1^2/((2+M1)*(b-1+b*M1)*(1+M1))
  nu13<--(1+M1)
  
  A1<-c(A11,A12,A13)
  nu1<-c(nu11,nu12,nu13)
  
  A2<-1
  nu2<--1/b
  
  A31<-b*M1/(b*M1-2)
  nu31<--1/b
  A32<-2/(2-b*M1)
  nu32<--M1/2
  
  A3<-c(A31,A32)
  nu3<-c(nu31,nu32)
  
  loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                      exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                               exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                      exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A1*exp(nu1*(tau0-tau1)))),na.rm=TRUE)
  
  loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                      exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                               exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                      exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
  
  loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                             exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                      exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
  
  
  -(loglike1+loglike2+loglike3)
} 


# "IIM2" model with M1=M2 and equal pop sizes during migration  ("IIM5")

negll.IIM.5<-function(params){
  
  a<-params[1]/params[2]
  theta<-params[2]
  b<-1
  c<-c(params[3]/theta,params[4]/theta)
  tau1<-params[5]/theta
  tau0<-params[6]/theta+tau1
  M<-vector(length=2,mode="numeric")
  M[1]<-params[7]
  M[2]<-params[7]
  
  if (M[1]>0 & M[2]>0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    
    Qt<-matrix(ncol=4,nrow=4)
    Qt[,1]<-c(-(1+M[1]),0,M[1],1)
    Qt[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
    Qt[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
    Qt[,4]<-c(0,0,0,0)
    X<-eigen(Qt)$vectors
    nu<-eigen(Qt)$values[-4]
    P.0<-diag(4)   #boundary conditions for all starting states form an identity matrix
    C<-solve(X) #(%*%P.0)
    A<--{X[4,-4]*C[-4,]} #minus sign precedes because of differentiation wrt t;
    
    
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A[,1]*(-nu*(R1*theta)^X1/(-nu+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu+theta*R1)*tau1)*ppois(X1,(-nu+theta*R1)*tau1)-
                                                                                                  exp(nu*(tau0-tau1)-theta*R1*tau0)*exp((-nu+theta*R1)*tau0)*ppois(X1,(-nu+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A[,1]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A[,2]*(-nu*(R2*theta)^X2/(-nu+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu+theta*R2)*tau1)*ppois(X2,(-nu+theta*R2)*tau1)-
                                                                                                  exp(nu*(tau0-tau1)-theta*R2*tau0)*exp((-nu+theta*R2)*tau0)*ppois(X2,(-nu+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A[,2]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A[,3]*(-nu*(R3*theta)^X3/(-nu+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu+theta*R3)*tau1)*ppois(X3,(-nu+theta*R3)*tau1)-
                                                                                exp(nu*(tau0-tau1)-theta*R3*tau0)*exp((-nu+theta*R3)*tau0)*ppois(X3,(-nu+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A[,3]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
  }
  
  
  if (M[1]==0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    loglike1<-sum(log(
      (c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)+
        exp(-tau1*(1/c[1]+theta*R1[1,]))*(   (theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
                                               exp(tau1*(1+theta*R1[1,]))*ppois(X1[1,],tau1*(1+theta*R1[1,]))-(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*
                                               exp(tau1*(1/c[1]+theta*R1[1,]))*ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,]))  )
      +
        exp(-(tau1/c[1]+tau0-tau1+theta*R1[1,]*tau0))*(   (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*
                                                            exp(tau0*(1/a+theta*R1[1,]))*ppois(X1[1,],tau0*(1/a+theta*R1[1,]))-(theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
                                                            exp(tau0*(1+theta*R1[1,]))*ppois(X1[1,],tau0*(1+theta*R1[1,]))  )
    ))
    
    
    loglike2<-sum(log(
      (c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)+
        exp(-tau1*(1/c[2]+theta*R2[1,]))*(   (b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
                                               exp(tau1*(1/b+theta*R2[1,]))*ppois(X2[1,],tau1*(1/b+theta*R2[1,]))-(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*
                                               exp(tau1*(1/c[2]+theta*R2[1,]))*ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,]))  )
      +
        exp(-(tau1/c[2]+(tau0-tau1)/b+theta*R2[1,]*tau0))*(   (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*
                                                                exp(tau0*(1/a+theta*R2[1,]))*ppois(X2[1,],tau0*(1/a+theta*R2[1,]))-(b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
                                                                exp(tau0*(1/b+theta*R2[1,]))*ppois(X2[1,],tau0*(1/b+theta*R2[1,]))  )
    ))
    
    
    loglike3<-sum(log(
      exp(-tau0*theta*R3[1,])*(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau0*(1/a+theta*R3[1,]))*ppois(X3[1,],tau0*(1/a+theta*R3[1,]))
    ))
    
  }
  
  
  -(loglike1+loglike2+loglike3)
}


# "IIM2" model with M1=M2 and allowing for unequal population sizes ("IIM6")

negll.IIM.6<-function(params){
  
  a<-params[1]/params[2]
  theta<-params[2]
  b<-params[3]/theta
  c<-c(params[4]/theta,params[5]/theta)
  tau1<-params[6]/theta
  tau0<-params[7]/theta+tau1
  M<-vector(length=2,mode="numeric")
  M[1]<-params[8]
  M[2]<-params[8]
  
  if (M[1]>0 & M[2]>0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    
    Qt<-matrix(ncol=4,nrow=4)
    Qt[,1]<-c(-(1+M[1]),0,M[1],1)
    Qt[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
    Qt[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
    Qt[,4]<-c(0,0,0,0)
    X<-eigen(Qt)$vectors
    nu<-eigen(Qt)$values[-4]
    P.0<-diag(4)   #boundary conditions for all starting states form an identity matrix
    C<-solve(X) #(%*%P.0)
    A<--{X[4,-4]*C[-4,]} #minus sign precedes because of differentiation wrt t;
    
    
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A[,1]*(-nu*(R1*theta)^X1/(-nu+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu+theta*R1)*tau1)*ppois(X1,(-nu+theta*R1)*tau1)-
                                                                                                  exp(nu*(tau0-tau1)-theta*R1*tau0)*exp((-nu+theta*R1)*tau0)*ppois(X1,(-nu+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A[,1]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A[,2]*(-nu*(R2*theta)^X2/(-nu+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu+theta*R2)*tau1)*ppois(X2,(-nu+theta*R2)*tau1)-
                                                                                                  exp(nu*(tau0-tau1)-theta*R2*tau0)*exp((-nu+theta*R2)*tau0)*ppois(X2,(-nu+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A[,2]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A[,3]*(-nu*(R3*theta)^X3/(-nu+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu+theta*R3)*tau1)*ppois(X3,(-nu+theta*R3)*tau1)-
                                                                                exp(nu*(tau0-tau1)-theta*R3*tau0)*exp((-nu+theta*R3)*tau0)*ppois(X3,(-nu+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A[,3]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
  }
  
  if (M[1]==0 & M[2]>0){
    
    R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                                 exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*A1*exp(nu1*(tau0-tau1))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                                 exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                               exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
    
  }
  
  if (M[1]>0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    A11<-(b^2*M[1]^2)/((b*M[1]-2)*(b-1+b*M[1]))
    nu11<--1/b
    A12<-4*M[1]/((2+M[1])*(2-b*M[1]))
    nu12<--M[1]/2
    A13<-1/(1+M[1])+M[1]^2/((2+M[1])*(b-1+b*M[1])*(1+M[1]))
    nu13<--(1+M[1])
    
    A1<-c(A11,A12,A13)
    nu1<-c(nu11,nu12,nu13)
    
    A2<-1
    nu2<--1/b
    
    A31<-b*M[1]/(b*M[1]-2)
    nu31<--1/b
    A32<-2/(2-b*M[1])
    nu32<--M[1]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                                 exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A1*exp(nu1*(tau0-tau1)))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                                 exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                               exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
    
  }
  
  if (M[1]==0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    loglike1<-sum(log(
      (c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)+
        exp(-tau1*(1/c[1]+theta*R1[1,]))*(   (theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
                                               exp(tau1*(1+theta*R1[1,]))*ppois(X1[1,],tau1*(1+theta*R1[1,]))-(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*
                                               exp(tau1*(1/c[1]+theta*R1[1,]))*ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,]))  )
      +
        exp(-(tau1/c[1]+tau0-tau1+theta*R1[1,]*tau0))*(   (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*
                                                            exp(tau0*(1/a+theta*R1[1,]))*ppois(X1[1,],tau0*(1/a+theta*R1[1,]))-(theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
                                                            exp(tau0*(1+theta*R1[1,]))*ppois(X1[1,],tau0*(1+theta*R1[1,]))  )
    ))
    
    
    loglike2<-sum(log(
      (c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)+
        exp(-tau1*(1/c[2]+theta*R2[1,]))*(   (b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
                                               exp(tau1*(1/b+theta*R2[1,]))*ppois(X2[1,],tau1*(1/b+theta*R2[1,]))-(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*
                                               exp(tau1*(1/c[2]+theta*R2[1,]))*ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,]))  )
      +
        exp(-(tau1/c[2]+(tau0-tau1)/b+theta*R2[1,]*tau0))*(   (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*
                                                                exp(tau0*(1/a+theta*R2[1,]))*ppois(X2[1,],tau0*(1/a+theta*R2[1,]))-(b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
                                                                exp(tau0*(1/b+theta*R2[1,]))*ppois(X2[1,],tau0*(1/b+theta*R2[1,]))  )
    ))
    
    
    loglike3<-sum(log(
      exp(-tau0*theta*R3[1,])*(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau0*(1/a+theta*R3[1,]))*ppois(X3[1,],tau0*(1/a+theta*R3[1,]))
    ))
    
  }
  
  
  -(loglike1+loglike2+loglike3)
}


# "IIM2" model with M1=M2=0 ("ISO2")

negll.iso.2<-function(params){
  
  a<-params[1]/params[2]
  theta<-params[2]
  b<-params[3]/theta
  c<-c(params[4]/theta,params[5]/theta)
  tau1<-params[6]/theta
  tau0<-params[7]/theta+tau1
  
  X1<-matrix(x1,nrow=1,byrow=TRUE)
  X2<-matrix(x2,nrow=1,byrow=TRUE)
  X3<-matrix(x3,nrow=1,byrow=TRUE)
  
  R1<-matrix(r1,nrow=1,byrow=TRUE)
  R2<-matrix(r2,nrow=1,byrow=TRUE)
  R3<-matrix(r3,nrow=1,byrow=TRUE)
  
  A1<-1
  nu1<--1
  
  A2<-1
  nu2<--1/b
  
  A3<-1
  nu3<-0
  
  loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                      exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                               exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                      exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A1*exp(nu1*(tau0-tau1)))),na.rm=TRUE)
  
  loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                      exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                               exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                      exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
  
  loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                             exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                      exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
  
  
  -(loglike1+loglike2+loglike3)
}




### Functions to calculate profile lik confidence intervals



CI.pl<-function(gim.list,parameter,model,quantile,ub=20){
  ci<-list(NA,NA)
  names(ci)<-c("lower", "upper")
  i<-which(names(gim.list$par)==parameter)
  params1<-unlist(gim.list$par)
  params<-params1[-i]
  lik<-gim.list$objective
  negll<-get(paste0("negll.",model,".pl"))
  
  ##lower bound search
  
  if(params1[i] <0.001){
    return("parameter too close to zero")
  }else if(params1[i]>=0.1){
    l1<-1
    l2<-1
    l3<-1
  }else if(params1[i]>=0.01){
    l1<-0
    l2<-1
    l3<-1
  }else if(params1[i]>=0.001){
    l1<-0
    l2<-0
    l3<-1
  }
  
  
  if(l1==1){
    pointvec<-seq(params1[i],0.0001,-0.1)
    w<-vector(length=0)
    for (j in 1:length(pointvec)){
      x<-nlminb(params,
                target=pointvec[j],i=i,negll,lower=0.00001)
      if (x$objective==0){
        cat("Lower level 1","j=",j,"no convergence","\n")
        break
      }
      w<-c(w,2*(x$objective-lik))
      cat("Lower level 1","j=",j,"lik.ratio=", w[j],"\n")
      if (w[j]>quantile){
        break
      }
      params<-x$par
    }
    
  }
  
  if(l2==1){
    if(l1==1){
      pointvec<-seq(pointvec[j-1],0.0001,-0.01)
    }
    if(l1==0){
      pointvec<-seq(params1[i],0.0001,-0.01)
    }
    w<-vector(length=0)
    for (j in 1:length(pointvec)){
      x<-nlminb(params,
                target=pointvec[j],i=i,negll,lower=0.000001)
      if (x$objective==0){
        cat("Lower level 2","j=",j,"no convergence","\n")
        break
      }
      w<-c(w,2*(x$objective-lik))
      cat("Lower level 2","j=",j,"lik.ratio=", w[j],"\n")
      if (w[j]>quantile){break}
      params<-x$par
    }
  }
  
  if(l3==1){
    if(l2==1){
      pointvec<-seq(pointvec[j-1],0.0001,-0.001)
    }
    if(l2==0){
      pointvec<-seq(params1[i],0.0001,-0.001)
    }
    
    w<-vector(length=0)
    for (j in 1:length(pointvec)){
      x<-nlminb(params,
                target=pointvec[j],i=i,negll,lower=0.00001)
      if (x$objective==0){
        cat("Lower level 3","j=",j,"no convergence","\n")
        break
      }
      w<-c(w,2*(x$objective-lik))
      cat("Lower level 3","j=",j,"lik.ratio=", w[j],"\n")
      if (w[j]>quantile){break}
      params<-x$par
    }
    
  }
  
  if (w[j]>quantile){
    # print(w[j])
    # print(x$message)
    ci$lower<-pointvec[j-1]
  }else{
    ci$lower<-paste0("< ",pointvec(j))
  }
  
  names(ci$lower)<-NULL
  
  ##upper bound search:    
  
  pointvec<-seq(params1[i],ub,0.1)
  params<-params1[-i]
  
  w<-vector(length=0)
  for (j in 1:length(pointvec)){
    x<-nlminb(params,
              target=pointvec[j],i=i,negll,lower=0.00001)
    if (x$objective==0){
      cat("Upper level 1","j=",j,"no convergence","\n")
      break
    }
    w<-c(w,2*(x$objective-lik))
    cat("Upper level 1","j=",j,"lik.ratio=", w[j],"\n")
    if (w[j]>quantile){break}
    params<-x$par
  }
  
  if(w[j]>quantile){
    pointvec<-seq(pointvec[j-1],ub,0.01)
    
  }else{
    ci$upper<-paste0("> ",pointvec[j])
    names(ci$upper)<-NULL
    return(list(lower=ci$lower,upper=ci$upper))
  }
  
  w<-vector(length=0)
  for (j in 1:length(pointvec)){
    x<-nlminb(params,
              target=pointvec[j],i=i,negll,lower=0.00001)
    if (x$objective==0){
      cat("Upper level 2","j=",j,"no convergence","\n")
      break
    }
    w<-c(w,2*(x$objective-lik))
    cat("Upper level 2","j=",j,"lik.ratio=", w[j],"\n")
    if (w[j]>quantile){break}
    params<-x$par
  }
  
  
  if(w[j]>quantile){
    pointvec<-seq(pointvec[j-1],ub,0.001)
    
  }else{
    ci$upper<-paste0("> ",pointvec[j])
    names(ci$upper)<-NULL
    return(list(lower=ci$lower,upper=ci$upper))
  }
  
  w<-vector(length=0)
  for (j in 1:length(pointvec)){
    x<-nlminb(params,
              target=pointvec[j],i=i,negll,lower=0.00001)
    if (x$objective==0){
      cat("Upper level 3","j=",j,"no convergence","\n")
      break
    }
    w<-c(w,2*(x$objective-lik))
    cat("Upper level 3","j=",j,"lik.ratio=", w[j],"\n")
    if (w[j]>quantile){
      break}
    params<-x$par
  }
  if(w[j]>quantile){
    ci$upper<-pointvec[j-1]
    names(ci$upper)<-NULL
    return(list(lower=ci$lower,upper=ci$upper))
    
  }else{
    ci$upper<-paste0("> ",pointvec[j])
    names(ci$upper)<-NULL
    return(list(lower=ci$lower,upper=ci$upper))
  }
} 

negll.IIM.3.pl<-function(params,target,i){
  
  M<-vector(length=2,mode="numeric")
  M[1]<-0
  
  if(i==1){
    a<-target/params[1]
    theta<-params[1]
    b<-params[2]/theta
    c<-c(params[3]/theta,params[4]/theta)
    tau1<-params[5]/theta
    tau0<-params[6]/theta+tau1
    M[2]<-params[7]
  }
  
  if(i==2){
    a<-params[1]/target
    theta<-target
    b<-params[2]/theta
    c<-c(params[3]/theta,params[4]/theta)
    tau1<-params[5]/theta
    tau0<-params[6]/theta+tau1
    M[2]<-params[7]
  }
  
  if(i==3){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-target/theta
    c<-c(params[3]/theta,params[4]/theta)
    tau1<-params[5]/theta
    tau0<-params[6]/theta+tau1
    M[2]<-params[7]
  }
  
  if(i==4){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-params[3]/theta
    c<-c(target/theta,params[4]/theta)
    tau1<-params[5]/theta
    tau0<-params[6]/theta+tau1
    M[2]<-params[7]
  }
  
  
  if(i==5){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-params[3]/theta
    c<-c(params[4]/theta,target/theta)
    tau1<-params[5]/theta
    tau0<-params[6]/theta+tau1
    M[2]<-params[7]
  }
  
  if(i==6){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-params[3]/theta
    c<-c(params[4]/theta,params[5]/theta)
    tau1<-target/theta
    tau0<-params[6]/theta+tau1
    M[2]<-params[7]
  }
  
  if(i==7){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-params[3]/theta
    c<-c(params[4]/theta,params[5]/theta)
    tau1<-params[6]/theta
    tau0<-target/theta+tau1
    M[2]<-params[7]
  }
  
  if(i==8){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-params[3]/theta
    c<-c(params[4]/theta,params[5]/theta)
    tau1<-params[6]/theta
    tau0<-params[7]/theta+tau1
    M[2]<-target
  }
  
  
  
  if (M[1]>0 & M[2]>0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    
    Qt<-matrix(ncol=4,nrow=4)
    Qt[,1]<-c(-(1+M[1]),0,M[1],1)
    Qt[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
    Qt[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
    Qt[,4]<-c(0,0,0,0)
    X<-eigen(Qt)$vectors
    nu<-eigen(Qt)$values[-4]
    P.0<-diag(4)   #boundary conditions for all starting states form an identity matrix
    C<-solve(X) #(%*%P.0)
    A<--{X[4,-4]*C[-4,]} #minus sign precedes because of differentiation wrt t;
    
    
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A[,1]*(-nu*(R1*theta)^X1/(-nu+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu+theta*R1)*tau1)*ppois(X1,(-nu+theta*R1)*tau1)-
                                                                                                  exp(nu*(tau0-tau1)-theta*R1*tau0)*exp((-nu+theta*R1)*tau0)*ppois(X1,(-nu+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A[,1]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A[,2]*(-nu*(R2*theta)^X2/(-nu+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu+theta*R2)*tau1)*ppois(X2,(-nu+theta*R2)*tau1)-
                                                                                                  exp(nu*(tau0-tau1)-theta*R2*tau0)*exp((-nu+theta*R2)*tau0)*ppois(X2,(-nu+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A[,2]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A[,3]*(-nu*(R3*theta)^X3/(-nu+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu+theta*R3)*tau1)*ppois(X3,(-nu+theta*R3)*tau1)-
                                                                                exp(nu*(tau0-tau1)-theta*R3*tau0)*exp((-nu+theta*R3)*tau0)*ppois(X3,(-nu+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A[,3]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
  }
  
  if (M[1]==0 & M[2]>0){
    
    R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                                 exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*A1*exp(nu1*(tau0-tau1))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                                 exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                               exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
    
  }
  
  if (M[1]>0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    A11<-(b^2*M[1]^2)/((b*M[1]-2)*(b-1+b*M[1]))
    nu11<--1/b
    A12<-4*M[1]/((2+M[1])*(2-b*M[1]))
    nu12<--M[1]/2
    A13<-1/(1+M[1])+M[1]^2/((2+M[1])*(b-1+b*M[1])*(1+M[1]))
    nu13<--(1+M[1])
    
    A1<-c(A11,A12,A13)
    nu1<-c(nu11,nu12,nu13)
    
    A2<-1
    nu2<--1/b
    
    A31<-b*M[1]/(b*M[1]-2)
    nu31<--1/b
    A32<-2/(2-b*M[1])
    nu32<--M[1]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                                 exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A1*exp(nu1*(tau0-tau1)))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                                 exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                               exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
    
  }
  
  if (M[1]==0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    loglike1<-sum(log(
      (c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)+
        exp(-tau1*(1/c[1]+theta*R1[1,]))*(   (theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
                                               exp(tau1*(1+theta*R1[1,]))*ppois(X1[1,],tau1*(1+theta*R1[1,]))-(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*
                                               exp(tau1*(1/c[1]+theta*R1[1,]))*ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,]))  )
      +
        exp(-(tau1/c[1]+tau0-tau1+theta*R1[1,]*tau0))*(   (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*
                                                            exp(tau0*(1/a+theta*R1[1,]))*ppois(X1[1,],tau0*(1/a+theta*R1[1,]))-(theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
                                                            exp(tau0*(1+theta*R1[1,]))*ppois(X1[1,],tau0*(1+theta*R1[1,]))  )
    ))
    
    
    loglike2<-sum(log(
      (c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)+
        exp(-tau1*(1/c[2]+theta*R2[1,]))*(   (b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
                                               exp(tau1*(1/b+theta*R2[1,]))*ppois(X2[1,],tau1*(1/b+theta*R2[1,]))-(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*
                                               exp(tau1*(1/c[2]+theta*R2[1,]))*ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,]))  )
      +
        exp(-(tau1/c[2]+(tau0-tau1)/b+theta*R2[1,]*tau0))*(   (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*
                                                                exp(tau0*(1/a+theta*R2[1,]))*ppois(X2[1,],tau0*(1/a+theta*R2[1,]))-(b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
                                                                exp(tau0*(1/b+theta*R2[1,]))*ppois(X2[1,],tau0*(1/b+theta*R2[1,]))  )
    ))
    
    
    loglike3<-sum(log(
      exp(-tau0*theta*R3[1,])*(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau0*(1/a+theta*R3[1,]))*ppois(X3[1,],tau0*(1/a+theta*R3[1,]))
    ))
    
  }
  
  
  -(loglike1+loglike2+loglike3)
}

negll.IIM.5.pl<-function(params,target,i){
  
  M<-vector(length=2,mode="numeric")
  
  if(i==1){
    a<-target/params[1]
    theta<-params[1]
    b<-1
    c<-c(params[2]/theta,params[3]/theta)
    tau1<-params[4]/theta
    tau0<-params[5]/theta+tau1
    M[1]<-params[6]
    M[2]<-params[6]
  }
  
  if(i==2){
    a<-params[1]/target
    theta<-target
    b<-1
    c<-c(params[2]/theta,params[3]/theta)
    tau1<-params[4]/theta
    tau0<-params[5]/theta+tau1
    M[1]<-params[6]
    M[2]<-params[6]
  }
  
  if(i==3){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-1
    c<-c(target/theta,params[3]/theta)
    tau1<-params[4]/theta
    tau0<-params[5]/theta+tau1
    M[1]<-params[6]
    M[2]<-params[6]
  }
  
  if(i==4){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-1
    c<-c(params[3]/theta,target/theta)
    tau1<-params[4]/theta
    tau0<-params[5]/theta+tau1
    M[1]<-params[6]
    M[2]<-params[6]
  }
  
  
  if(i==5){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-1
    c<-c(params[3]/theta,params[4]/theta)
    tau1<-target/theta
    tau0<-params[5]/theta+tau1
    M[1]<-params[6]
    M[2]<-params[6]
  }
  
  if(i==6){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-1
    c<-c(params[3]/theta,params[4]/theta)
    tau1<-params[5]/theta
    tau0<-target/theta+tau1
    #tau0<-target/theta
    M[1]<-params[6]
    M[2]<-params[6]
  }
  
  if(i==7){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-1
    c<-c(params[3]/theta,params[4]/theta)
    tau1<-params[5]/theta
    tau0<-params[6]/theta+tau1
    M[1]<-target
    M[2]<-target
  }
  
  
  
  
  if (M[1]>0 & M[2]>0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    
    Qt<-matrix(ncol=4,nrow=4)
    Qt[,1]<-c(-(1+M[1]),0,M[1],1)
    Qt[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
    Qt[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
    Qt[,4]<-c(0,0,0,0)
    X<-eigen(Qt)$vectors
    nu<-eigen(Qt)$values[-4]
    P.0<-diag(4)   #boundary conditions for all starting states form an identity matrix
    C<-solve(X) #(%*%P.0)
    A<--{X[4,-4]*C[-4,]} #minus sign precedes because of differentiation wrt t;
    
    
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A[,1]*(-nu*(R1*theta)^X1/(-nu+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu+theta*R1)*tau1)*ppois(X1,(-nu+theta*R1)*tau1)-
                                                                                                  exp(nu*(tau0-tau1)-theta*R1*tau0)*exp((-nu+theta*R1)*tau0)*ppois(X1,(-nu+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A[,1]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A[,2]*(-nu*(R2*theta)^X2/(-nu+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu+theta*R2)*tau1)*ppois(X2,(-nu+theta*R2)*tau1)-
                                                                                                  exp(nu*(tau0-tau1)-theta*R2*tau0)*exp((-nu+theta*R2)*tau0)*ppois(X2,(-nu+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A[,2]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A[,3]*(-nu*(R3*theta)^X3/(-nu+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu+theta*R3)*tau1)*ppois(X3,(-nu+theta*R3)*tau1)-
                                                                                exp(nu*(tau0-tau1)-theta*R3*tau0)*exp((-nu+theta*R3)*tau0)*ppois(X3,(-nu+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A[,3]*exp(nu*(tau0-tau1)))),na.rm=TRUE)
  }
  
  if (M[1]==0 & M[2]>0){
    
    R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                                 exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*A1*exp(nu1*(tau0-tau1))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                                 exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                               exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
    
  }
  
  if (M[1]>0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    A11<-(b^2*M[1]^2)/((b*M[1]-2)*(b-1+b*M[1]))
    nu11<--1/b
    A12<-4*M[1]/((2+M[1])*(2-b*M[1]))
    nu12<--M[1]/2
    A13<-1/(1+M[1])+M[1]^2/((2+M[1])*(b-1+b*M[1])*(1+M[1]))
    nu13<--(1+M[1])
    
    A1<-c(A11,A12,A13)
    nu1<-c(nu11,nu12,nu13)
    
    A2<-1
    nu2<--1/b
    
    A31<-b*M[1]/(b*M[1]-2)
    nu31<--1/b
    A32<-2/(2-b*M[1])
    nu32<--M[1]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-sum(log((c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*(1-ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,])))+
                        exp(-tau1/c[1])*colSums(A1*(-nu1*(R1*theta)^X1/(-nu1+theta*R1)^(X1+1)*(exp(-theta*R1*tau1)*exp((-nu1+theta*R1)*tau1)*ppois(X1,(-nu1+theta*R1)*tau1)-
                                                                                                 exp(nu1*(tau0-tau1)-theta*R1*tau0)*exp((-nu1+theta*R1)*tau0)*ppois(X1,(-nu1+theta*R1)*tau0))))+
                        exp(-tau1/c[1]-theta*R1[1,]*tau0)*(a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp((1/a+theta*R1[1,])*tau0)*ppois(X1[1,],(1/a+theta*R1[1,])*tau0)*sum(A1*exp(nu1*(tau0-tau1)))),na.rm=TRUE)
    
    loglike2<-sum(log((c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*(1-ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,])))+
                        exp(-tau1/c[2])*colSums(A2*(-nu2*(R2*theta)^X2/(-nu2+theta*R2)^(X2+1)*(exp(-theta*R2*tau1)*exp((-nu2+theta*R2)*tau1)*ppois(X2,(-nu2+theta*R2)*tau1)-
                                                                                                 exp(nu2*(tau0-tau1)-theta*R2*tau0)*exp((-nu2+theta*R2)*tau0)*ppois(X2,(-nu2+theta*R2)*tau0))))+
                        exp(-tau1/c[2]-theta*R2[1,]*tau0)*(a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp((1/a+theta*R2[1,])*tau0)*ppois(X2[1,],(1/a+theta*R2[1,])*tau0)*sum(A2*exp(nu2*(tau0-tau1)))),na.rm=TRUE)
    
    loglike3<-sum(log(colSums(A3*(-nu3*(R3*theta)^X3/(-nu3+theta*R3)^(X3+1)*(exp(-theta*R3*tau1)*exp((-nu3+theta*R3)*tau1)*ppois(X3,(-nu3+theta*R3)*tau1)-
                                                                               exp(nu3*(tau0-tau1)-theta*R3*tau0)*exp((-nu3+theta*R3)*tau0)*ppois(X3,(-nu3+theta*R3)*tau0))))+
                        exp(-R3[1,]*theta*tau0)*(a*R3[1,]*theta)^X3[1,]/(1+a*R3[1,]*theta)^(X3[1,]+1)*exp((1/a+theta*R3[1,])*tau0)*ppois(X3[1,],(1/a+theta*R3[1,])*tau0)*sum(A3*exp(nu3*(tau0-tau1)))),na.rm=TRUE)
    
  }
  
  if (M[1]==0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    loglike1<-sum(log(
      (c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)+
        exp(-tau1*(1/c[1]+theta*R1[1,]))*(   (theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
                                               exp(tau1*(1+theta*R1[1,]))*ppois(X1[1,],tau1*(1+theta*R1[1,]))-(c[1]*theta*R1[1,])^X1[1,]/(1+c[1]*theta*R1[1,])^(X1[1,]+1)*
                                               exp(tau1*(1/c[1]+theta*R1[1,]))*ppois(X1[1,],tau1*(1/c[1]+theta*R1[1,]))  )
      +
        exp(-(tau1/c[1]+tau0-tau1+theta*R1[1,]*tau0))*(   (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*
                                                            exp(tau0*(1/a+theta*R1[1,]))*ppois(X1[1,],tau0*(1/a+theta*R1[1,]))-(theta*R1[1,])^X1[1,]/(1+theta*R1[1,])^(X1[1,]+1)*
                                                            exp(tau0*(1+theta*R1[1,]))*ppois(X1[1,],tau0*(1+theta*R1[1,]))  )
    ))
    
    
    loglike2<-sum(log(
      (c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)+
        exp(-tau1*(1/c[2]+theta*R2[1,]))*(   (b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
                                               exp(tau1*(1/b+theta*R2[1,]))*ppois(X2[1,],tau1*(1/b+theta*R2[1,]))-(c[2]*theta*R2[1,])^X2[1,]/(1+c[2]*theta*R2[1,])^(X2[1,]+1)*
                                               exp(tau1*(1/c[2]+theta*R2[1,]))*ppois(X2[1,],tau1*(1/c[2]+theta*R2[1,]))  )
      +
        exp(-(tau1/c[2]+(tau0-tau1)/b+theta*R2[1,]*tau0))*(   (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*
                                                                exp(tau0*(1/a+theta*R2[1,]))*ppois(X2[1,],tau0*(1/a+theta*R2[1,]))-(b*theta*R2[1,])^X2[1,]/(1+b*theta*R2[1,])^(X2[1,]+1)*
                                                                exp(tau0*(1/b+theta*R2[1,]))*ppois(X2[1,],tau0*(1/b+theta*R2[1,]))  )
    ))
    
    
    loglike3<-sum(log(
      exp(-tau0*theta*R3[1,])*(a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau0*(1/a+theta*R3[1,]))*ppois(X3[1,],tau0*(1/a+theta*R3[1,]))
    ))
    
  }
  
  
  -(loglike1+loglike2+loglike3)
}

negll.iso.1.pl<-function(params,target,i){
  
  if(i==1){
    a<-target/params[1]
    theta<-params[1]
    b<-params[2]/theta
    tau<-params[3]/theta
  }
  
  if(i==2){
    a<-params[1]/target
    theta<-target
    b<-params[2]/theta
    tau<-params[3]/theta
  }
  
  if(i==3){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-target/theta
    tau<-params[3]/theta
  }
  
  if(i==4){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-params[3]/theta
    tau<-target/theta
  }
  
  R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
  R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
  R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
  
  X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
  X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
  X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
  
  loglike1<-log(((theta*R1[1,])^X1[1,])/      	
                  ((1+theta*R1[1,])^(X1[1,]+1))*							
                  (1-ppois(X1[1,],tau*(1+theta*R1[1,])))+  	
                  exp(-tau*(1+theta*R1[1,]))*
                  (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp(tau*(1/a+theta*R1[1,]))*ppois(X1[1,],tau*(1/a+theta*R1[1,])))
  
  loglike2<-log(((b*theta*R2[1,])^X2[1,])/				
                  ((1+b*theta*R2[1,])^(X2[1,]+1))*							
                  (1-ppois(X2[1,],tau*(1/b+theta*R2[1,])))+  	
                  exp(-tau*(1/b+theta*R2[1,]))*
                  (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp(tau*(1/a+theta*R2[1,]))*ppois(X2[1,],tau*(1/a+theta*R2[1,])))
  
  loglike3<-log(exp(-tau*theta*R3[1,])*			
                  (a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau*(1/a+theta*R3[1,]))*ppois(X3[1,],tau*(1/a+theta*R3[1,])))
  
  -sum(c(loglike1,loglike2,loglike3))
}

negll.IM.4.pl<-function(params,target,i){
  
 if(i==1){
  a<-target/params[1]
  theta<-params[1]
  b<-params[2]/theta
  tau<-params[3]/theta
  M<-vector(length=2,mode="numeric")
  M[1]<-params[4]
  M[2]<-params[4]
 }
  
  if(i==2){
    a<-params[1]/target
    theta<-target
    b<-params[2]/theta
    tau<-params[3]/theta
    M<-vector(length=2,mode="numeric")
    M[1]<-params[4]
    M[2]<-params[4]
  }
  
  if(i==3){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-target/theta
    tau<-params[3]/theta
    M<-vector(length=2,mode="numeric")
    M[1]<-params[4]
    M[2]<-params[4]
  }
  
  if(i==4){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-params[3]/theta
    tau<-target/theta
    M<-vector(length=2,mode="numeric")
    M[1]<-params[4]
    M[2]<-params[4]
  }
  
  if(i==5){
    a<-params[1]/params[2]
    theta<-params[2]
    b<-params[3]/theta
    tau<-params[4]/theta
    M<-vector(length=2,mode="numeric")
    M[1]<-target
    M[2]<-target
  }
  
  
  if (M[1]>0 & M[2]>0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    Qt<-matrix(ncol=4,nrow=4)
    Qt[,1]<-c(-(1+M[1]),0,M[1],1)
    Qt[,2]<-c(0,-(1/b+M[2]),M[2],1/b)
    Qt[,3]<-c(M[2]/2,M[1]/2,-(M[1]+M[2])/2,0)
    Qt[,4]<-c(0,0,0,0)
    X<-eigen(Qt)$vectors
    nu<-eigen(Qt)$values[-4]
    P.0<-diag(4) 	#boundary conditions for all starting states form an identity matrix
    C<-solve(X) #(%*%P.0)
    A<--{X[4,-4]*C[-4,]} #minus sign precedes because of differentiation wrt t;
    
    
    loglike1<-log(colSums(A[,1]*((-nu*(theta*R1)^X1)/
                                   ((theta*R1-nu)^(X1+1))*
                                   (1-ppois(X1,tau*(-nu+theta*R1)))+
                                   exp(-tau*(theta*R1-nu))*
                                   (a*theta*R1)^X1/
                                   (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A[,2]*((-nu*(theta*R2)^X2)/
                                   ((theta*R2-nu)^(X2+1))*
                                   (1-ppois(X2,tau*(-nu+theta*R2)))+
                                   exp(-tau*(theta*R2-nu))*
                                   (a*theta*R2)^X2/
                                   (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A[,3]*((-nu*(theta*R3)^X3)/
                                   ((theta*R3-nu)^(X3+1))*
                                   (1-ppois(X3,tau*(-nu+theta*R3)))+
                                   exp(-tau*(theta*R3-nu))*
                                   (a*theta*R3)^X3/
                                   (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
  }
  
  if (M[1]==0 & M[2]>0){
    
    R1<-matrix(rep(r1,1),nrow=1,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,1),nrow=1,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    
    A1<-1
    nu1<--1
    
    A21<-(b*M[2]^2)/((-2 + M[2])*(1-b+b*M[2]))
    nu21<--1
    A22<-4*b*M[2]/((2-M[2])*(2+b*M[2]))
    nu22<--M[2]/2
    A23<-(1/b)/(1/b+M[2])+b^2*M[2]^2/((2+b*M[2])*(1-b+b*M[2])*(1/b+M[2]))
    nu23<--(1/b+M[2])
    
    A2<-c(A21,A22,A23)
    nu2<-c(nu21,nu22,nu23)
    
    A31<-M[2]/(M[2]-2)
    nu31<--1
    A32<-2/(2-M[2])
    nu32<--M[2]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    
    loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
                                (-nu1+theta*R1)^(X1+1)*
                                (1-ppois(X1,tau*(-nu1+theta*R1)))+
                                exp(-tau*(theta*R1-nu1))*
                                (a*theta*R1)^X1/
                                (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
                                (-nu2+theta*R2)^(X2+1)*
                                (1-ppois(X2,tau*(-nu2+theta*R2)))+
                                exp(-tau*(theta*R2-nu2))*
                                (a*theta*R2)^X2/
                                (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
                                (-nu3+theta*R3)^(X3+1)*
                                (1-ppois(X3,tau*(-nu3+theta*R3)))+
                                exp(-tau*(theta*R3-nu3))*
                                (a*theta*R3)^X3/
                                (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
  }
  
  if (M[1]>0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,1),nrow=1,byrow=TRUE)
    R3<-matrix(rep(r3,2),nrow=2,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,1),nrow=1,byrow=TRUE)
    X3<-matrix(rep(x3,2),nrow=2,byrow=TRUE)
    
    A11<-(b^2*M[1]^2)/((b*M[1]-2)*(b-1+b*M[1]))
    nu11<--1/b
    A12<-4*M[1]/((2+M[1])*(2-b*M[1]))
    nu12<--M[1]/2
    A13<-1/(1+M[1])+M[1]^2/((2+M[1])*(b-1+b*M[1])*(1+M[1]))
    nu13<--(1+M[1])
    
    A1<-c(A11,A12,A13)
    nu1<-c(nu11,nu12,nu13)
    
    A2<-1
    nu2<--1/b
    
    A31<-b*M[1]/(b*M[1]-2)
    nu31<--1/b
    A32<-2/(2-b*M[1])
    nu32<--M[1]/2
    
    A3<-c(A31,A32)
    nu3<-c(nu31,nu32)
    
    loglike1<-log(colSums(A1*((-nu1*(R1*theta)^X1)/
                                (-nu1+theta*R1)^(X1+1)*
                                (1-ppois(X1,tau*(-nu1+theta*R1)))+
                                exp(-tau*(theta*R1-nu1))*
                                (a*theta*R1)^X1/
                                (1+a*theta*R1)^(X1+1)*exp(tau*(1/a+theta*R1))*ppois(X1,tau*(1/a+theta*R1)))))
    
    loglike2<-log(colSums(A2*((-nu2*(R2*theta)^X2)/
                                (-nu2+theta*R2)^(X2+1)*
                                (1-ppois(X2,tau*(-nu2+theta*R2)))+
                                exp(-tau*(theta*R2-nu2))*
                                (a*theta*R2)^X2/
                                (1+a*theta*R2)^(X2+1)*exp(tau*(1/a+theta*R2))*ppois(X2,tau*(1/a+theta*R2)))))
    
    loglike3<-log(colSums(A3*((-nu3*(theta*R3)^X3)/
                                (-nu3+theta*R3)^(X3+1)*
                                (1-ppois(X3,tau*(-nu3+theta*R3)))+
                                exp(-tau*(theta*R3-nu3))*
                                (a*theta*R3)^X3/
                                (1+a*theta*R3)^(X3+1)*exp(tau*(1/a+theta*R3))*ppois(X3,tau*(1/a+theta*R3)))))
    
    
  }
  
  if (M[1]==0 & M[2]==0){
    
    R1<-matrix(rep(r1,3),nrow=3,byrow=TRUE)
    R2<-matrix(rep(r2,3),nrow=3,byrow=TRUE)
    R3<-matrix(rep(r3,3),nrow=3,byrow=TRUE)
    
    X1<-matrix(rep(x1,3),nrow=3,byrow=TRUE)
    X2<-matrix(rep(x2,3),nrow=3,byrow=TRUE)
    X3<-matrix(rep(x3,3),nrow=3,byrow=TRUE)
    
    loglike1<-log(((theta*R1[1,])^X1[1,])/				
                    ((1+theta*R1[1,])^(X1[1,]+1))*							
                    (1-ppois(X1[1,],tau*(1+theta*R1[1,])))+  	
                    exp(-tau*(1+theta*R1[1,]))*
                    (a*theta*R1[1,])^X1[1,]/(1+a*theta*R1[1,])^(X1[1,]+1)*exp(tau*(1/a+theta*R1[1,]))*ppois(X1[1,],tau*(1/a+theta*R1[1,])))
    
    loglike2<-log(((b*theta*R2[1,])^X2[1,])/				
                    ((1+b*theta*R2[1,])^(X2[1,]+1))*							
                    (1-ppois(X2[1,],tau*(1/b+theta*R2[1,])))+  	
                    exp(-tau*(1/b+theta*R2[1,]))*
                    (a*theta*R2[1,])^X2[1,]/(1+a*theta*R2[1,])^(X2[1,]+1)*exp(tau*(1/a+theta*R2[1,]))*ppois(X2[1,],tau*(1/a+theta*R2[1,])))
    
    loglike3<-log(exp(-tau*theta*R3[1,])*			
                    (a*theta*R3[1,])^X3[1,]/(1+a*theta*R3[1,])^(X3[1,]+1)*exp(tau*(1/a+theta*R3[1,]))*ppois(X3[1,],tau*(1/a+theta*R3[1,])))
    
  }
  
  -sum(c(loglike1,loglike2,loglike3),na.rm=TRUE)
}
