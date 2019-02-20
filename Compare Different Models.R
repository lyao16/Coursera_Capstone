
################################################################
############   All the packages we need in the code ############
################################################################

library(doParallel)
library(flare)
library(mvtnorm)
library(glmnet)
library(foreach)
library(lars)
library(quantreg)
library(pheatmap)
library(plotrix)
cl <- makeCluster(2)
registerDoParallel(cl)


#####################################################################################
########## a) Elastic net with alpha and labda decided via cross-validation #########
#####################################################################################


enet<-function(Data)
{
  require(glmnet)
  k=10
  fold<- sample(1:k, size=length(Data$y), replace=TRUE)
  alpha = as.vector(seq(0,1, length=100))
  best.lambda = vector(length=length(alpha))
  min.mse = vector(length=length(alpha))
  lowest.coef = list()
  for(i in 1:length(alpha))
  {
    cvout= cv.glmnet(Data$X, Data$y, foldid=fold,alpha=alpha[i],parallel=TRUE)
    best.lambda[i] = cvout$lambda.min
    lowest.coef[[i]]=as.numeric(coef(cvout,s=cvout$lambda.min))
    min.mse[i] = cvout$cvm[which(cvout$lambda==cvout$lambda.min)]
  }
  index = which.min(min.mse)
  return(lowest.coef[[index]])
}


############################################################
########b) Lad lasso based on flare package ################
############################################################


ladlasso <-function(Data){
  x1=Data$X
  y1=Data$y
  out1=slim(x1, y1,method="lq", q=1,nlambda =5,verbose=FALSE)
  return(out1$beta0[[5]])
}

#################################################################
################ c) sqrt lasso based on flare package ###########
#################################################################

sqrtlasso <-function(Data)
{
  x1=Data$X
  y1=Data$y
  out2=slim(x1, y1, method="lq", q=2, nlambda=5, verbose=FALSE, lambda.min.value=sqrt(log(ncol(x1))/nrow(x1)))
  return(out2$beta0[[5]])
}

#### coef(out2,lambda.idx = (1:5), beta.idx = (1:10))

#################################################################
############## d) l_q lasso based on flare package  #############
#################################################################


lqlasso <-function(Data)
{
  x1=Data$X
  y1=Data$y
  out3=slim(x1, y1, method="lq", q=1.5, nlambda=5,verbose=FALSE)
  return(out3$beta0[[5]])
}

##################################################################
############ e)  Dantzig selector based on flare package #########
##################################################################


dantziglasso <-function(Data)
{
  require(flare)
  x1=Data$X
  y1=Data$y
  out4=slim(x1, y1, method='dantzig', nlambda=5, verbose=FALSE,lambda.min.ratio=0.35)
  return(out4$beta0[[5]])
}

#########################################################################
########### f) adaptive lasso based on the LS approximation #############
#########################################################################

require(lars)
lsa <- function(obj) 
{
  intercept <- attr(obj$terms,'intercept')
  if(class(obj)[1]=='coxph') intercept <- 0
  n <- length(obj$residuals)  
  Sigma <- vcov(obj) 
  SI <- solve(Sigma)
  beta.ols <- coef(obj)
  l.fit <- lars.lsa(SI, beta.ols, intercept, n)  
  t1 <- sort(l.fit$BIC, ind=T)
  t2 <- sort(l.fit$AIC, ind=T)
  beta <- l.fit$beta
  if(intercept) {
    beta0 <- l.fit$beta0+beta.ols[1]
    beta.bic <- c(beta0[t1$ix[1]],beta[t1$ix[1],])
    beta.aic <- c(beta0[t2$ix[1]],beta[t2$ix[1],]) 
  }
  else { 
    beta0 <- l.fit$beta0
    beta.bic <- beta[t1$ix[1],]
    beta.aic <- beta[t2$ix[1],]
  }
  obj <- list(beta.ols=beta.ols, beta.bic=beta.bic, beta.aic = beta.aic)  
  obj 
}

## lars variant for LSA

lars.lsa <- function (Sigma0, b0, intercept,  n, type = c("lasso", "lar"),    
                     eps = .Machine$double.eps,max.steps) 
  
{
  type <- match.arg(type)
  TYPE <- switch(type, lasso = "LASSO", lar = "LAR")
  n1 <- dim(Sigma0)[1]
  ## handle intercept
  
  if (intercept) {
    a11 <- Sigma0[1,1]
    a12 <- Sigma0[2:n1,1]
    a22 <- Sigma0[2:n1,2:n1]
    Sigma <- a22-outer(a12,a12)/a11
    b <- b0[2:n1]
    beta0 <- crossprod(a12,b)/a11
  }
  else {
    Sigma <- Sigma0
    b <- b0 
  }
  Sigma <- diag(abs(b))%*%Sigma%*%diag(abs(b))
  b <- sign(b)
  nm <- dim(Sigma)
  m <- nm[2]
  im <- inactive <- seq(m)
  Cvec <- drop(t(b)%*%Sigma)
  ssy <- sum(Cvec*b)
  if (missing(max.steps)) 
    max.steps <- 8 * m
  beta <- matrix(0, max.steps + 1, m)
  Gamrat <- NULL
  arc.length <- NULL
  R2 <- 1
  RSS <- ssy
  first.in <- integer(m)
  active <- NULL
  actions <- as.list(seq(max.steps))
  drops <- FALSE
  Sign <- NULL
  R <- NULL
  k <- 0
  ignores <- NULL
  while ((k < max.steps) & (length(active) < m)) {
    action <- NULL
    k <- k + 1
    C <- Cvec[inactive]
    Cmax <- max(abs(C))
    if (!any(drops)) {
      new <- abs(C) >= Cmax - eps
      C <- C[!new]
      new <- inactive[new]
      for (inew in new) {
        R <- updateR(Sigma[inew, inew], R, drop(Sigma[inew, active]),           
                     Gram = TRUE,eps=eps)
        if(attr(R, "rank") == length(active)) { 
          ##singularity; back out
          nR <- seq(length(active))
          R <- R[nR, nR, drop = FALSE]
          attr(R, "rank") <- length(active)
          ignores <- c(ignores, inew)
          action <- c(action,  - inew) 
        }
        else { 
          if(first.in[inew] == 0)
            first.in[inew] <- k
          active <- c(active, inew)
          Sign <- c(Sign, sign(Cvec[inew])) 
          action <- c(action, inew)   
        }   
      } 
    }
    else action <- -dropid
    Gi1 <- backsolve(R, backsolvet(R, Sign))
    dropouts <- NULL
    A <- 1/sqrt(sum(Gi1 * Sign))
    w <- A * Gi1    
    if (length(active) >= m) { 
      gamhat <- Cmax/A        
    } 
    else {           
      a <- drop(w %*% Sigma[active, -c(active,ignores), drop = FALSE])
      gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
      gamhat <- min(gam[gam > eps], Cmax/A)
    }
    if (type == "lasso") {
      dropid <- NULL
      b1 <- beta[k, active]
      z1 <- -b1/w
      zmin <- min(z1[z1 > eps], gamhat)
      # cat('zmin ',zmin, ' gamhat ',gamhat,'\n') 
      if (zmin < gamhat) {
        gamhat <- zmin
        drops <- z1 == zmin  
      }
      else drops <- FALSE  
    }
    beta[k + 1, ] <- beta[k, ]
    beta[k + 1, active] <- beta[k + 1, active] + gamhat * w
    Cvec <- Cvec - gamhat * Sigma[, active, drop = FALSE] %*% w    
    Gamrat <- c(Gamrat, gamhat/(Cmax/A))
    arc.length <- c(arc.length, gamhat)
    if (type == "lasso" && any(drops)) {
      dropid <- seq(drops)[drops]
      for (id in rev(dropid)) {
        R <- downdateR(R,id) 
      }
      dropid <- active[drops] 
      beta[k + 1, dropid] <- 0
      active <- active[!drops]
      Sign <- Sign[!drops] 
    }
    actions[[k]] <- action
    inactive <- im[-c(active)] 
  }
  beta <- beta[seq(k + 1), ]
  dff <- b-t(beta)
  RSS <- diag(t(dff)%*%Sigma%*%dff)  
  if(intercept)  
    beta <- t(abs(b0[2:n1])*t(beta))
  else 
    beta <- t(abs(b0)*t(beta))
  if (intercept) { 
    beta0 <- beta0-drop(t(a12)%*%t(beta))/a11 
  }
  else { 
    beta0 <- rep(0,k+1)  
  }
  dof <- apply(abs(beta)>eps,1,sum)
  BIC <- RSS+log(n)*dof
  AIC <- RSS+2*dof
  object <- list(AIC = AIC, BIC = BIC, 
  beta = beta, beta0 = beta0)  
  object 
}



lsa.linear<-function(x,y){
  # adaptive lasso for linear reg, tuning parameter by bic 
  # calls software from Wang and Leng (2007, JASA).
  m<-ncol(x)
  n<-nrow(x)
  x<-as.matrix(x)
  out<-lm(y~x)
  out.lsa<-lsa(out)
  coeff<-out.lsa$beta.bic
  coeff2<-coeff[2:(m+1)]               # get rid of intercept
  pred<-x%*%coeff2+coeff[1]
  st<-sum(coeff2 !=0)                                          # number nonzero
  mse<-sum((y-pred)^2)/(n-st-1)
  if(st>0) x.ind<-as.vector(which(coeff2 !=0)) else x.ind<-0
  return(list(mse=mse,coeff=coeff2,intercept=coeff[1]))
}

adaptivelasso <- function(Data)
{
  
  x=Data$X
  y1=Data$y
  out5=lsa.linear(x=x,y=y1)
  return(out5$coeff)
}



###############################################################
############ g) 4 variants of adaptive lad lasso ##############
###############################################################


adapladlasso1<- function(Data)
{
  require(quantreg)
  
  QR1<- coef(rq(y ~ X, tau=0.5, data=Data, method="lasso", lambda=1 ))
  lam= sqrt(1.5*nrow(Data$X)*log(ncol(Data$X)))/abs(QR1)
  lam[1]=0
  QR.coef1<- coef(rq(y~ X, tau=.5, data=Data, method="lasso", lambda=lam))
  return(QR.coef1)
}

adapladlasso2<- function(Data)
{
  require(quantreg)
  QR1<- coef(rq(y ~ X, tau=0.5, data=Data, method="lasso", lambda=1 ))
  lam= sqrt(2*nrow(Data$X)*log(ncol(Data$X)))/abs(QR1)
  lam[1]=0
  QR.coef2<- coef(rq(y~ X, tau=.5, data=Data, method="lasso", lambda=lam))
  return(QR.coef2)
  
}
adapladlasso3<- function(Data)
{
  require(quantreg)
  QR1<- coef(rq(y ~ X, tau=0.5, data=Data, method="lasso", lambda=1 ))
  lam= sqrt(4*nrow(Data$X)*log(ncol(Data$X)))/abs(QR1)
  lam[1]=0
  QR.coef3<- coef(rq(y~ X, tau=.5, data=Data, method="lasso", lambda=lam))
  return(QR.coef3)
  
}

adapladlasso4<- function(Data)
{
  require(quantreg)
  QR1<- coef(rq(y ~ X, tau=0.5, data=Data, method="lasso", lambda=1 ))
  lam= sqrt(10*nrow(Data$X)*log(ncol(Data$X)))/abs(QR1)
  lam[1]=0
  QR.coef4<- coef(rq(y~ X, tau=.5, data=Data, method="lasso", lambda=lam))
  return(QR.coef4)
}


lasso <- function(Data)
{
  require(glmnet)
  cvfit <- cv.glmnet(Data$X,Data$y, parallel=TRUE)
  return(as.numeric(coef(cvfit,s=cvfit$lambda.min)))
}



###############################################
##############  Generating the data ###########
###############################################


genData<-function(n, p, beta){
  
  x<- matrix(rmvt(n, delta=rep(0,length(beta)), sigma=diag(p), df=3 ),ncol=p,nrow=n)
  y<- rnorm(n, x%*%beta, 2)
  return(list(X=x,y=y))
}
genData2 <- function( n, p, beta, rho )
{
  s= c(1, rho^seq(1:(length(beta)-1)))
  s1= toeplitz(s)
  
  x<- matrix(rmvt(n, delta=rep(0,length(beta)), sigma=s1, df= 3), nrow=n, ncol=p)
  y<- rnorm(n, x%*%beta, 2)
  return(list(X=x,y=y))
}

###  x=rt(100000, df=3, ncp=0)

###  n=100000

###  hist(x, breaks=500,prob=TRUE, xlim=c(-5, 5),ylim=c(0,0.45),main="t vs. normal")

###  hist(rmvt(10000,sigma=diag(1),df=10),breaks=100,prob=TRUE)
###  curve(dnorm(x,0,1), add=TRUE,col="red", lwd=2)


##############################################################
######################  Scenario 1 ###########################
##############################################################

n=100
N=150
beta=c(2, -2, 1, -1,rep(0,16))
p=length(beta)
results = array(NA, dim=c(N, p, 12),
                dimnames =list(1:N, 1:p,c("Elastic Net", "Lad Lasso", "Sqrt Lasso", "lq Lasso",
                                          "Dantzig Selector", "Adaptive Lasso", "Ad.Lad.Lasso 1",
                                          "Ad.Lad.Lasso 2", "Ad.Lad.Lasso 3",
                                          "Ad.Lad.Lasso 4", "OLS","Lasso")))

dimnames =c("Elastic Net", "Lad Lasso", "Sqrt Lasso", "lq Lasso",
            "Dantzig Selector", "Adaptive Lasso", "Ad.Lad.Lasso 1",
            "Ad.Lad.Lasso 2", "Ad.Lad.Lasso 3",
            "Ad.Lad.Lasso 4", "OLS","Lasso")

for(i in 1:N){
  Data = genData(n, p, beta)
  results[i,,1]<- enet(Data)[-1]
  results[i,,2]<- ladlasso(Data) 
  results[i,,3]<- sqrtlasso(Data)
  print(i)
  results[i,,4]<- lqlasso(Data)
  results[i,,5]<- dantziglasso(Data)
  results[i,,6]<- adaptivelasso(Data)
  results[i,,7]<- adapladlasso1(Data)[-1]
  results[i,,8]<- adapladlasso2(Data)[-1]
  results[i,,9]<- adapladlasso3(Data)[-1]
  results[i,,10]<- adapladlasso4(Data)[-1]
  results[i,,11]<- coef(lm(y~X,Data))[-1]
  results[i,,12]<- lasso(Data)[-1]
  
}

B1<- apply(results, 2:3, mean)-beta
V1<- apply(results, 2:3, var)
MSE1<- B1^2+V1
s1 = round(apply(MSE1, 2, sum),5)

pheatmap(B1, cluster_rows=F, cluster_cols=F, cellwidth =30, cellheight = 15, labels_row = c(1:20),labels_col =dimnames, main="Bias Scenario 1")

pheatmap(V1, cluster_rows=F, cluster_cols=F, cellwidth =30, cellheight = 15, labels_row = c(1:20),labels_col =dimnames, main="Variance Scenario 1")

write.csv(s1,row.names=dimnames,file="C:/Users/Lanlan/Downloads/STP 598 Computational Statistics/hwk/HWK4/s1s.csv")

ylim=range(results)
par(mfrow=c(2:3))
for(i in 1:12)
{
  boxplot(results[,,i], col="gray", ylim=ylim, main=dimnames[i], xlab="Scenario 1")
  abline(h=c(2,-2,1,-1),lty=2, col="red")
}



####################################################################
########################  Scenario 2  ##############################
####################################################################

n=100
N=150
beta2=rep(0.2,20)
p=length(beta2)
results2 = array(NA, dim=c(N, p, 12),
                 dimnames =list(1:N, 1:p,c("Elastic Net", "Lad Lasso", "Sqrt Lasso", "lq Lasso",
                                           "Dantzig Selector", "Adaptive Lasso", "Ad.Lad.Lasso 1",
                                           "Ad.Lad.Lasso 2", "Ad.Lad.Lasso 3",
                                           "Ad.Lad.Lasso 4", "OLS","Lasso")))


for(i in 1:N)
{
  Data = genData(n, p, beta2)
  results2[i,,1]<- enet(Data)[-1]
  results2[i,,2]<- ladlasso(Data) 
  results2[i,,3]<- sqrtlasso(Data)
  print(i)
  results2[i,,4]<- lqlasso(Data)
  results2[i,,5]<- dantziglasso(Data)
  results2[i,,6]<- adaptivelasso(Data)
  results2[i,,7]<- adapladlasso1(Data)[-1]
  results2[i,,8]<- adapladlasso2(Data)[-1]
  results2[i,,9]<- adapladlasso3(Data)[-1]
  results2[i,,10]<- adapladlasso4(Data)[-1]
  results2[i,,11]<- coef(lm(y~X,Data))[-1]
  results2[i,,12]<- lasso(Data)[-1]
  
}

B2<- apply(results2, 2:3, mean)-beta2
V2<- apply(results2, 2:3, var)
MSE2<- B2^2+V2
s2 = round(apply(MSE2, 2, sum),5)

    

pheatmap(B2, cluster_rows=F, cluster_cols=F, cellwidth =30, cellheight = 15, labels_row = c(1:20),labels_col =dimnames, main="Bias Scenario 2")

pheatmap(V2, cluster_rows=F, cluster_cols=F, cellwidth =30, cellheight = 15, labels_row = c(1:20),labels_col =dimnames, main="Variance Scenario 2")

write.csv(s2,row.names=dimnames,file="C:/Users/Lanlan/Downloads/STP 598 Computational Statistics/hwk/HWK4/s2.csv")


ylim=range(results2)
par(mfrow=c(2,3))
for(k in 1:12)
{
  boxplot(results2[,,k], col="gray", ylim=ylim, main=dimnames[k], xlab="Scenario 2")
  abline(h=0.2,lty=2, col="red")
}



###########################################################################
############################  Scenario 3  #################################
###########################################################################

n=100
N=150
beta3=c(2, -2, 1, -1, 0.5, 0.2, -0.3, -0.15, rep(0,12))
p=length(beta3)
rho3=0.9
results3 = array(NA, dim=c(N, p, 12),
                 dimnames =list(1:N, 1:p,c("Elastic Net", "Lad Lasso", "Sqrt Lasso", "lq Lasso",
                                           "Dantzig Selector", "Adaptive Lasso", "Ad.Lad.Lasso 1",
                                           "Ad.Lad.Lasso 2", "Ad.Lad.Lasso 3",
                                           "Ad.Lad.Lasso 4", "OLS","Lasso")))


for(i in 1:N)  
{
  
  Data = genData2(n, p, beta3, rho3)
  results3[i,,1]<- enet(Data)[-1]
  results3[i,,2]<- ladlasso(Data) 
  results3[i,,3]<- sqrtlasso(Data)
  print(i)
  results3[i,,4]<- lqlasso(Data)
  results3[i,,5]<- dantziglasso(Data)
  results3[i,,6]<- adaptivelasso(Data)
  results3[i,,7]<- adapladlasso1(Data)[-1]
  results3[i,,8]<- adapladlasso2(Data)[-1]
  results3[i,,9]<- adapladlasso3(Data)[-1]
  results3[i,,10]<- adapladlasso4(Data)[-1]
  results3[i,,11]<- coef(lm(y~X,Data))[-1]
  results3[i,,12]<- lasso(Data)[-1]
  
}

B3<- apply(results3, 2:3, mean)-beta3
V3<- apply(results3, 2:3, var)
MSE3<- B3^2+V3
s3 = round(apply(MSE3, 2, sum),5)


pheatmap(B3, cluster_rows=F, cluster_cols=F, cellwidth =30, cellheight = 15, labels_row = c(1:20),labels_col =dimnames, main="Bias Scenario 3")

pheatmap(V3, cluster_rows=F, cluster_cols=F, cellwidth =30, cellheight = 15, labels_row = c(1:20),labels_col =dimnames, main="Variance Scenario 3")

write.csv(s3,row.names=dimnames,file="C:/Users/Lanlan/Downloads/STP 598 Computational Statistics/hwk/HWK4/s3.csv")


par(mfrow=c(2,3))
for(k in 1:12)
{

  boxplot(results3[,,k], col="gray", main=dimnames[k], xlab="Scenario 3")
  abline(h=c(2, -2, 1, -1, 0.5, 0.2, -0.3, -0.15),lty=2, col="red")
}



####################################################################
########################  Scenario 4  ##############################
####################################################################

n=100
N=150
beta4=c(2, 1, 0.5, rep(0, 217))
p=length(beta4)
results4 = array(NA, dim=c(N, p, 10),
                 dimnames = list(1:N,1:p, c("Elastic Net", "Lad Lasso", "Sqrt Lasso", "lq Lasso",
                                            "Dantzig Selector",  "Ad.Lad.Lasso 1", "Ad.Lad.Lasso 2", "Ad.Lad.Lasso 3",
                                            "Ad.Lad.Lasso 4","Lasso")))

dimnames2 =c("Elastic Net", "Lad Lasso", "Sqrt Lasso", "lq Lasso",
             "Dantzig Selector",  "Ad.Lad.Lasso 1", "Ad.Lad.Lasso 2", "Ad.Lad.Lasso 3",
             "Ad.Lad.Lasso 4","Lasso")


for(i in 1:N)
{
  Data = genData(n, p, beta4)
  results4[i,,1]<- enet(Data)[-1]
  results4[i,,2]<- ladlasso(Data) 
  results4[i,,3]<- sqrtlasso(Data)
  print(i)
  results4[i,,4]<- lqlasso(Data)
  results4[i,,5]<- dantziglasso(Data)
# results4[i,,6]<- adaptivelasso(Data)
  results4[i,,6]<- adapladlasso1(Data)[-1]
  results4[i,,7]<- adapladlasso2(Data)[-1]
  results4[i,,8]<- adapladlasso3(Data)[-1]
  results4[i,,9]<- adapladlasso4(Data)[-1]
# results4[i,,11]<- coef(lm(y~X,Data))[-1]
  results4[i,,10]<- lasso(Data)[-1]
  
}



B4<- apply(results4, 2:3, mean)-beta4
V4<- apply(results4, 2:3, var)
MSE4<- B4^2+V4
s4 = round(apply(MSE4, 2, sum),5)


pheatmap(B4[1:20,], cluster_rows=F, cluster_cols=F, cellwidth =30, cellheight = 15, labels_row = c(1:20),labels_col =dimnames2, main="Bias Scenario 4")

pheatmap(V4[1:20,], cluster_rows=F, cluster_cols=F, cellwidth =30, cellheight = 15, labels_row = c(1:20),labels_col =dimnames2, main="Variance Scenario 4")

write.csv(s4,row.names=dimnames2,file="C:/Users/Lanlan/Downloads/STP 598 Computational Statistics/hwk/HWK4/s4.csv")


ylim=range(results4[,1:20,])
par(mfrow=c(2,3))
for(k in 1:10)
{
  boxplot(results4[,1:20,k], col="gray", main=dimnames2[k], xlab="Scenario 4")
  abline(h=c(2,1,0.5),lty=2, col="red")
}

####################################################################
########################  Scenario 5  ##############################
####################################################################


n=100
N=150
beta5=c(2, -2, 1, -1, 0.5, 0.2, -0.3, -0.15, rep(0,212))
p=length(beta5)
rho5=0.9
results5 = array(NA, dim=c(N, p, 10),  
                 dimnames = list(1:N,1:p, c("Elastic Net", "Lad Lasso", "Sqrt Lasso", "lq Lasso",
                                            "Dantzig Selector",  "Ad.Lad.Lasso 1", "Ad.Lad.Lasso 2", "Ad.Lad.Lasso 3",
                                            "Ad.Lad.Lasso 4","Lasso")))


for(i in 1:N) 
{
  Data = genData2(n, p, beta5, rho5)
  results5[i,,1]<- enet(Data)[-1]
  results5[i,,2]<- ladlasso(Data) 
  results5[i,,3]<- sqrtlasso(Data)
  print(i)
  results5[i,,4]<- lqlasso(Data)
  results5[i,,5]<- dantziglasso(Data)
# results5[i,,6]<- adaptivelasso(Data)
  results5[i,,6]<- adapladlasso1(Data)[-1]
  results5[i,,7]<- adapladlasso2(Data)[-1]
  results5[i,,8]<- adapladlasso3(Data)[-1]
  results5[i,,9]<- adapladlasso4(Data)[-1]
# results5[i,,11]<- coef(lm(y~X,Data))[-1]
  results5[i,,10]<- lasso(Data)[-1]
  
}

B5<- apply(results5, 2:3, mean)-beta5
V5<- apply(results5, 2:3, var)
MSE5<- B5^2+V5
s5 = round(apply(MSE5, 2, sum))



pheatmap(B5[1:20,], cluster_rows=F, cluster_cols=F, cellwidth =30, cellheight = 15, labels_row = c(1:20),labels_col =dimnames2, main="Bias Scenario 5")

pheatmap(V5[1:20,], cluster_rows=F, cluster_cols=F, cellwidth =30, cellheight = 15, labels_row = c(1:20),labels_col =dimnames2, main="Variance Scenario 5")

write.csv(s5,row.names=dimnames2,file="C:/Users/Lanlan/Downloads/STP 598 Computational Statistics/hwk/HWK4/s5a.csv")


ylim=range(results5[,1:20,])
par(mfrow=c(2,3))
for(k in 1:10)
{
  boxplot(results5[,1:20,k], col="gray", ylim=ylim, main=dimnames2[k], xlab="Scenario 5")
  abline(h=c(2, -2, 1, -1, 0.5, 0.2, -0.3, -0.15),lty=2, col="red")
}

ComMSE= cbind(s1,s2, s3, c(s4[1:5],NA,s4[6:9],NA,s4[10]), c(s5[1:5],NA,s5[6:9],NA,s5[10]))
rownames(ComMSE)=c("Elastic net","Lad Lasso", "Sqrt Lasso", "lq Lasso",
           "Dantzig Selector", "Adaptive Lasso", "Ad.Lad.Lasso 1",
           "Ad.Lad.Lasso 2", "Ad.Lad.Lasso 3",
           "Ad.Lad.Lasso 4", "OLS","Lasso" )

write.csv(ComMSE,row.names=TRUE,file="C:/Users/Lanlan/Downloads/STP 598 Computational Statistics/hwk/HWK4/ComMSE.csv")

save.image(file="C:/Users/Lanlan/Downloads/STP 598 Computational Statistics/hwk/HWK4/hwk4a.RData")

stopCluster(cl) 









