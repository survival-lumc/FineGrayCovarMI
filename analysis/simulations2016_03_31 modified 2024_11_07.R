library(mitools)
library(survival)
library(mice)
library(smcfcs)

cumhazest <- function(data) {
  coxmod <- coxph(Surv(t,d==1)~x1+x2+x3, data)
  survobj <- survfit(coxmod,newdata=data.frame(x1=0,x2=0,x3=0))
  cumhaz <- survobj$cumhaz[survobj$time<H0evaltime]
  cumhazse <- survobj$std.err[survobj$time<H0evaltime]
  c(tail(cumhaz,1),tail(cumhazse,1))*100
}

ciok <- function(lower,upper,true) {
  100*((lower<true) & (upper>true))
}

#set this to 1,2,3 to choose scenario
scenario <- 1

nSim <- 100
numImps <- 5

n <- 1000
beta11 <- 1
beta12 <- 1
beta13 <- 1
beta21 <- 0.5
beta22 <- -1
beta23 <- 0.75
if (scenario==3) {
  beta23 <- 0
}
a1 <- 4
a2 <- 1
lambda10 <- 0.2
lambda20 <- 1
H0evaltime <- 0.5
H01_100 <- 100*0.2*H0evaltime^a1

numParms <- 7

full_parm <- array(0, dim=c(nSim,numParms))
full_ci <- array(0, dim=c(nSim,numParms))

cca_parm <- array(0, dim=c(nSim,numParms))
cca_ci <- array(0, dim=c(nSim,numParms))

smcfcs1_parm <- array(0, dim=c(nSim,numParms))
smcfcs1_ci <- array(0, dim=c(nSim,numParms))
smcfcs2_parm <- array(0, dim=c(nSim,numParms))
smcfcs2_ci <- array(0, dim=c(nSim,numParms))

mice1_parm <- array(0, dim=c(nSim,numParms))
mice1_ci <- array(0, dim=c(nSim,numParms))
mice2_parm <- array(0, dim=c(nSim,numParms))
mice2_ci <- array(0, dim=c(nSim,numParms))
mice3_parm <- array(0, dim=c(nSim,numParms))
mice3_ci <- array(0, dim=c(nSim,numParms))


set.seed(6981233)

generateData <- function() {
  x1 <- rbinom(n,1,0.5)
  x2 <- rbinom(n,1,0.25+0.50*x1)
  x3 <- -1+x1+x2+rnorm(n)
  xb1 <- beta11*x1+beta12*x2+beta13*x3
  xb2 <- beta21*x1+beta22*x2+beta23*x3
  
  lambda1 <- lambda10*exp(xb1)
  b1 <- lambda1^(-1/a1)
  t1 <- rweibull(n,shape=a1,scale=b1)
  
  lambda2 <- lambda20*exp(xb2)
  b2 <- lambda2^(-1/a2)
  t2 <- rweibull(n,shape=a2,scale=b2)
  
  t <- t1
  t[t2<t1] <- t2[t2<t1]
  d <- rep(1,n)
  d[t2<t1] <- 2
  rm(t1,t2)
  if (scenario==1) {
    c <- 0.5+1.5*runif(n)
    d[c<t] <- 0
    t[c<t] <- c[c<t]
  }
  
  dataframe <- data.frame(x1,x2,x3,t,d)
  dataframe[order(dataframe$t),]
}

if (scenario==2) {
  makeMissing <- function(dataframe) {
    dataframe$x3[runif(n)<(0.25+0.5*(dataframe$d-1))] <- NA
    dataframe
  }
} else if (scenario==3){
  makeMissing <- function(dataframe) {
    missxb <- dataframe$x3/1.32
    dataframe$x3[runif(n)<(exp(missxb)/(1+exp(missxb)))] <- NA
    dataframe
  } 
} else {
  makeMissing <- function(dataframe) {
    dataframe$x3[runif(n)<(0.25+0.5*dataframe$x1)] <- NA
    dataframe
  } 
}

analyseBoth <- function(imputations) {
  impobj <- imputationList(imputations)
  models <- with(impobj, coxph(Surv(t,d==1)~x1+x2+x3))
  temp <- summary(MIcombine(models))
  parmvec <- array(0, dim=numParms)
  civec <- array(0, dim=numParms)
  
  parmvec[1:3] <- temp[,1]
  civec[1] <- ciok(temp[,3][1],temp[,4][1],beta11)
  civec[2] <- ciok(temp[,3][2],temp[,4][2],beta12)
  civec[3] <- ciok(temp[,3][3],temp[,4][3],beta13)
  
  models <- with(impobj, coxph(Surv(t,d==2)~x1+x2+x3))
  temp <- summary(MIcombine(models))
  parmvec[4:6] <- temp[,1]
  civec[4] <- ciok(temp[,3][1],temp[,4][1],beta21)
  civec[5] <- ciok(temp[,3][2],temp[,4][2],beta22)
  civec[6] <- ciok(temp[,3][3],temp[,4][3],beta23)
  
  cumhazestimates <- with(impobj, fun=cumhazest)
  temp <- summary(MIcombine(as.list(sapply(cumhazestimates, function(x) x[1])), variances=as.list((sapply(cumhazestimates, function(x) x[2]))^2)))
  parmvec[numParms] <- temp[,1]
  civec[numParms] <- ciok(temp[,3],temp[,4],H01_100)
  list(parmvec, civec)
}

analyseFull <- function(data) {
  Mod1 <- coxph(Surv(t,d==1)~x1+x2+x3, data)
  Mod2 <- coxph(Surv(t,d==2)~x1+x2+x3, data)
  CumHaz <- cumhazest(data)
  
  parmvec <- c(Mod1$coef, Mod2$coef, CumHaz[1])
  civec <- rep(0, numParms)
  civec[1] <- ciok(Mod1$coef[1]-1.96*summary(Mod1)$coefficients[1,3],Mod1$coef[1]+1.96*summary(Mod1)$coefficients[1,3],beta11)
  civec[2] <- ciok(Mod1$coef[2]-1.96*summary(Mod1)$coefficients[2,3],Mod1$coef[2]+1.96*summary(Mod1)$coefficients[2,3],beta12)
  civec[3] <- ciok(Mod1$coef[3]-1.96*summary(Mod1)$coefficients[3,3],Mod1$coef[3]+1.96*summary(Mod1)$coefficients[3,3],beta13)
  civec[4] <- ciok(Mod2$coef[1]-1.96*summary(Mod2)$coefficients[1,3],Mod2$coef[1]+1.96*summary(Mod2)$coefficients[1,3],beta21)
  civec[5] <- ciok(Mod2$coef[2]-1.96*summary(Mod2)$coefficients[2,3],Mod2$coef[2]+1.96*summary(Mod2)$coefficients[2,3],beta22)
  civec[6] <- ciok(Mod2$coef[3]-1.96*summary(Mod2)$coefficients[3,3],Mod2$coef[3]+1.96*summary(Mod2)$coefficients[3,3],beta23)
  civec[7] <- ciok(CumHaz[1]-1.96*CumHaz[2],CumHaz[1]+1.96*CumHaz[2],H01_100)
  list(parmvec, civec)
}


for (sim in 1:nSim) {
  print(sim)

  mydata <- generateData()
  
  #full data analysis
  fullDataResults <- analyseFull(mydata)
  full_parm[sim,] <- fullDataResults[[1]]
  full_ci[sim,] <- fullDataResults[[2]]
    
  #make data missing  
  mydata <- makeMissing(mydata)
  
  #complete case analysis
  ccaDataResults <- analyseFull(mydata)
  cca_parm[sim,] <- ccaDataResults[[1]]
  cca_ci[sim,] <- ccaDataResults[[2]]
  
  #now use standard mice approach with norm
  # margmod1 <- coxph(Surv(t,d==1)~1, mydata)
  # NA1 <- basehaz(margmod1, centered=FALSE)[,1]
  # margmod2 <- coxph(Surv(t,d==2)~1, mydata)
  # NA2 <- basehaz(margmod2, centered=FALSE)[,1]
  # 
  # micedata <- cbind(mydata,NA1,NA2,d1=1*(mydata$d==1),d2=1*(mydata$d==2))
  # predictorMatrix <- array(0, dim=c(ncol(micedata),ncol(micedata)))
  # predictorMatrix[3,c(1,2,6,7,8,9)] <- 1
  # miceimps <- mice(micedata,m=numImps,method="norm",predictorMatrix=predictorMatrix)
  # 
  # imps <- vector("list", numImps)
  # for (impnum in 1:numImps) {
  #   imps[[impnum]] <- complete(miceimps,impnum)
  # }
  # result <- analyseBoth(imps)
  # mice1_parm[sim,] <- result[[1]]
  # mice1_ci[sim,] <- result[[2]]
  # 
  # #mice including interactions
  # micedata <- cbind(mydata,NA1,NA2,mydata$x1*NA1,mydata$x1*NA2,mydata$x2*NA1,mydata$x2*NA2,d1=1*(mydata$d==1),d2=1*(mydata$d==2))
  # predictorMatrix <- array(0, dim=c(ncol(micedata),ncol(micedata)))
  # predictorMatrix[3,c(1,2,5,6,7,8,9,10,11,12,13)] <- 1
  # miceimps <- mice(micedata,m=numImps,method="norm",predictorMatrix=predictorMatrix)
  # 
  # imps <- vector("list", numImps)
  # for (impnum in 1:numImps) {
  #   imps[[impnum]] <- complete(miceimps,impnum)
  # }
  # result <- analyseBoth(imps)
  # mice2_parm[sim,] <- result[[1]]
  # mice2_ci[sim,] <- result[[2]]
  # 
  # #mice ignoring second cause of failure
  # igndata <- mydata
  # igndata$d[igndata$d==2] <- 0
  # margmod1 <- coxph(Surv(t,d)~1, igndata)
  # NA1 <- basehaz(margmod1, centered=FALSE)[,1]
  # 
  # micedata <- cbind(igndata,NA1)
  # predictorMatrix <- array(0, dim=c(ncol(micedata),ncol(micedata)))
  # predictorMatrix[3,c(1,2,5,6)] <- 1
  # miceimps <- mice(micedata,m=numImps,method="norm",predictorMatrix=predictorMatrix)
  # 
  # imps <- vector("list", numImps)
  # for (impnum in 1:numImps) {
  #   imps[[impnum]] <- complete(miceimps,impnum)
  #   imps[[impnum]]$d <- mydata$d
  # }
  # result <- analyseBoth(imps)
  # mice3_parm[sim,] <- result[[1]]
  # mice3_ci[sim,] <- result[[2]]
  # 
  #smcfcs 1
  # imps <- smcfcs(mydata,smtype="compet",smformula=c("Surv(t,d==1)~x1+x2+x3", "Surv(t,d==2)~x1+x2+x3"),method=c("","","norm","",""),m=numImps)
  # 
  # result <- analyseBoth(imps$impDatasets)
  # smcfcs1_parm[sim,] <- result[[1]]
  # smcfcs1_ci[sim,] <- result[[2]]
  
  #smcfcs ignoring second competing risk
  igndata <- mydata
  igndata$d[igndata$d==2] <- 0
  
  imps <- smcfcs(igndata,smtype="coxph",smformula="Surv(t,d)~x1+x2+x3",method=c("","","norm","",""),m=numImps)
  
  #regenerate full failure cause indicator
  for (impnum in 1:numImps) {
    imps$impDatasets[[impnum]]$d <- mydata$d
  }
  result <- analyseBoth(imps$impDatasets)
  smcfcs2_parm[sim,] <- result[[1]]
  smcfcs2_ci[sim,] <- result[[2]]
  
}

library(xtable)

meanprint <- function(parm,ci) {
  result <- array(0, dim=c(numParms))
  for (i in 1:numParms) {
    result[i] <- format(round(mean(parm[,i]),2),nsmall=2)
  }
  result
}

meantable <- rbind(meanprint(full_parm, full_ci), meanprint(cca_parm, cca_ci), meanprint(smcfcs2_parm, smcfcs2_ci))
meantable <- cbind(c("Full data","Complete case","SMC-FCS survival"), meantable)
meantable

t.test(smcfcs2_parm[,1],mu=1)
t.test(smcfcs2_parm[,2],mu=1)