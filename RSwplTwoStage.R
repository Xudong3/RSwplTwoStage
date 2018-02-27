

#setting: notation
N1=20 ## number of sampling cluster in the first stage (population level) 
N2=50 ##number of elements in each sampling cluster (population level)
latitude<-1:N2
longitude<-1:N1
population<-expand.grid(lat=latitude,long=longitude)
population$PSU<-population$long
overlap=ceiling(N2*3/4)


model_cluster<-function(population, overlap){
   population$cluster<-numeric(nrow(population))
   
   id<-ifelse(population$lat<=overlap, 
              population$long, 
              ((population$long+population$lat-overlap) %% N1)+1
   )
   population$cluster<-id
   population	
}

population<-model_cluster(population, overlap)
T=length(unique(population$cluster))

##check
table(population$cluster==population$strata)
table(table(population$cluster))
table(table(population$strata))

#Model: parameter from random slope model  model
truebeta1=1
truebeta2=3
truesigma2=4
truetau2_11=2
truetau2_12=0.4
truetau2_22=3
PairCov<-matrix(c(truetau2_11, truetau2_12, truetau2_12, truetau2_22), nrow=2, byrow=T)
###check positive definite 
#install.packages("matrixcalc")
library("matrixcalc")
is.positive.definite(PairCov)

truevalue<-c(truebeta1,truebeta2, truesigma2, truetau2_11, truetau2_12, truetau2_22)
names(truevalue)<-c("beta1", "beta2", "sigma2", "tau2_11", "tau2_12", "tau2_22")

##Population data
#install.packages("MASS")
library("MASS")
#install.packages("rockchalk")
library(rockchalk)
re=mvrnorm(n=T, mu = c(0,0), Sigma = PairCov) #generate vector of random effect (a, b)
population$a<-re[,1][population$cluster]
population$b<-re[,2][population$cluster]

population$x<-rnorm(N1*N2)+rnorm(T)[population$cluster]
population$y<-with(population, truebeta1+a+truebeta2*x+b*x+rnorm(N1*N2,s=sqrt(truesigma2)))
population$r=with(population, x*(y-truebeta1-truebeta2*x))
population$ID_unit=with(population, 1:(N1*N2))

#uninformative two-stage sampling design (first-stage: SRSWOR, Second-stage:SRSWOR)
n1=ceiling(N1/10) ##number of sampling cluster in the first stage (sample level)
n2=ceiling(N2/10) ##umber of elements in each sampling cluster ( sample level)

# Using sampling package for two-stage sampling (First-stage: SRSWOR, Second-stage: SRSWOR ) 
#install.packages("sampling")
library("sampling")

##uninformative two-stage  sampling design (First-stage: SRSWOR, Second-stage: SRSWOR) and extracts the observed data
##first-stage
FirststageSRSWOR=srswor(n1, N1)
FirststageSRSWORSample=subset(population, population$PSU%in% which(FirststageSRSWOR==1))

#second-stage
SecondstageSRSWOR=unlist(lapply(rep(n2,n1), function(v) return(srswor(v, N2))))
TwostageSRSWORSample<-FirststageSRSWORSample[c(which(SecondstageSRSWOR==1)),] 

#informative two-stage sampling design (first-stage: SRSWOR, Second-stage:SRSWOR)
##number of elements in each sampling cluster
param=c(0.05, 3.5)
n2informative= function(r, sc, param, N2){
   a=rep(NA, length=length(unique(population$sc)))
   b=rep(NA, length=length(unique(population$sc)))
   for (i in unique(sc)){
      a[i]=mean(r[sc==i])
      b[i]=2*ceiling((param[1]*exp(-param[2]*a[i]))/(1 +param[1]*exp(-param[2]*a[i]))*N2/2)
   }
   b
}

##informative two-stage  sampling design (SRSWOR)[second-stage is informative] and extracts the observed data
###second-stage
n2pop=n2informative(population$r,population$PSU, param ,N2)
n2is=n2pop*FirststageSRSWOR
SecondstageSRSWORis=unlist(lapply(n2is[c(which(n2is!=0))], function(v) return(srswor(v, N2))))
TwostageSRSWORSampleis=FirststageSRSWORSample[c(which(SecondstageSRSWORis==1)), ]

# Estimation: full-likelihood
#install.packages("lme4")
library(lme4)

##Census estimator 
lmer(y~(1+x|cluster)+x,data=population)


##uninformative two-stage sampling design (SRSWOR)
lmer(y~(1+x|cluster)+x,data=TwostageSRSWORSample)

##informative two-stage sampling design (SRSWOR)
lmer(y~(1+x|cluster)+x,data=TwostageSRSWORSampleis)

# Estimation: pairwise likelihood (without weight)
l2<-function(y1,y2, g1,g2, x1,x2, alpha, beta, sigma2, tau2_11, tau2_12, tau2_22){
   pc11<-tau2_11+2*x1*tau2_12+(x1^2)*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+(x2^2)*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   #(-(r1*r1*pc-2*r1*r2*iftau+r2*r2*pc22)/2/det-log(det)/2)
   
   -log(det)/2-(1/2)*(1/det)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)
   
}	


dalpha<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2,tau2_11, tau2_12, tau2_22){
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dr1<- -1
   dr2<- -1
   
   
   #((r1*pc22-r2*iftau-r1*iftau+r2*pc11)/det)
   (-1/2)*(1/det)*(2*r1*dr1*pc22-2*dr1*r2*pc12-2*r1*dr2*pc12+2*r2*dr2*pc11 )
   
}

dbeta<-function(y1,y2, g1,g2, x1,x2,alpha,beta, sigma2,tau2_11, tau2_12, tau2_22){
   
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dr1<- -x1
   dr2<- -x2
   
   (-1/2)*(1/det)*(2*r1*dr1*pc22-2*dr1*r2*pc12-2*r1*dr2*pc12+2*r2*dr2*pc11)
}	

dsigma2<-function(y1,y2, g1,g2, x1,x2,alpha, beta, sigma2,tau2_11, tau2_12, tau2_22){
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dpc11<-1
   dpc22<-1
   dpc12<-0
   
   ddet<-dpc11*pc22+pc11*dpc22-2*pc12*dpc12
   
   (-1/2)*(ddet/det)-1/2*(-ddet)/(det)^2*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
}

dtau2_11<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11, tau2_12, tau2_22){
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dpc11<-1
   dpc22<-1
   dpc12<-ifelse(g1==g2, 1, 0)
   ddet<- dpc11*pc22+pc11*dpc22-2*pc12*dpc12
   
   
   (-1/2)*(ddet/det)-1/2*(-ddet)/(det^2)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
   #ddet<-2*(st-iftau)
   #(-( (r1*r1-2*r1*r2*ifistau+r2*r2)/2/det - (r1*r1*st-2*r1*r2*iftau+r2*r2*st)*(ddet)/2/det/det) -ddet/det/2)
   
   
}	


dtau2_12<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11, tau2_12, tau2_22){
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dpc11<-2*x1
   dpc22<-2*x2
   dpc12<-ifelse(g1==g2,x1+x2, 0)
   ddet<- dpc11*pc22+pc11*dpc22-2*pc12*dpc12
   
   -1/2*ddet/det-1/2*(-ddet)/(det^2)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
   #   (-( (r1*r1-2*r1*r2*ifistau+r2*r2)/2/det - (r1*r1*st-2*r1*r2*iftau+r2*r2*st)*(ddet)/2/det/det) -ddet/det/2)
}

dtau2_22<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2_11, tau2_12, tau2_22){
   pc11<-tau2_11+2*x1*tau2_12+x1^2*tau2_22+sigma2 #pairwise covariance for 11
   pc22<-tau2_11+2*x2*tau2_12+x2^2*tau2_22+sigma2 #pairwise covariance for 22
   pc12<-ifelse(g1==g2, tau2_11+x1*tau2_12+x2*tau2_12+x1*x2*tau2_22, 0) #pairwise covariance for 12
   
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   
   det<-pc11*pc22-pc12^2
   
   dpc11<-x1^2
   dpc22<-x2^2
   dpc12<-ifelse(g1==g2, x1*x2, 0 )
   ddet<- dpc11*pc22+pc11*dpc22-2*pc12*dpc12
   
   -1/2*ddet/det-1/2*(-ddet)/(det^2)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
   #  -1/2*ddet/det-1/2*(-ddet)/(det^2)*(r1^2*pc22-2*r1*r2*pc12+r2^2*pc11)-1/2*1/det*(r1^2*dpc22-2*r1*r2*dpc12+r2^2*dpc11)
}

#optimization problem for pariwise likelihood estimation (without weight)
fast_pl<-function(y,g,x, theta){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau2_12=exp(theta[5]), tau2_22=exp(theta[6]))
   sum(increment)/T
}


fast_fit<-function(y,g,x, pars){
   n<-length(y)
   ij=expand.grid(1:n,1:n)
   ij<-ij[ij[,1]<ij[,2],]
   ij<-ij[g[ij[,1]]==g[ij[,2]],]
   i<-ij[,1]
   j<-ij[,2]
   
   func1<-function(theta){
      increment=l2(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                   sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau2_12=exp(theta[5]), tau2_22=exp(theta[6]))
      sum(increment)/T
   }
   gr<-function(theta){
      incrementda=dalpha(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                         sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau2_12=exp(theta[5]), tau2_22=exp(theta[6]))
      incrementdb=dbeta(y[i],y[j],g[i],g[j],x[i],x[j],alpha=theta[1],beta=theta[2],
                        sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau2_12=exp(theta[5]), tau2_22=exp(theta[6]))
      incrementds=exp(theta[3])*dsigma2(y[i],y[j],g[i],g[j],x[i],x[j],
                                        alpha=theta[1],beta=theta[2],sigma2=exp(theta[3]),tau2_11=exp(theta[4]),tau2_12=exp(theta[5]),
                                        tau2_22=exp(theta[6]))
      incrementdt_11=exp(theta[4])*dtau2_11(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                            sigma2=exp(theta[3]), tau2_11=exp(theta[4]),tau2_12=exp(theta[5]), tau2_22=exp(theta[6]))
      incrementdt_12=exp(theta[5])*dtau2_12(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                            sigma2=exp(theta[3]), tau2_11=exp(theta[4]),tau2_12=exp(theta[5]), tau2_22=exp(theta[6]))
      incrementdt_22=exp(theta[6])*dtau2_22(y[i],y[j],g[i],g[j],x[i],x[j], alpha=theta[1],beta=theta[2],
                                            sigma2=exp(theta[3]), tau2_11=exp(theta[4]),tau2_12=exp(theta[5]), tau2_22=exp(theta[6]))
      c(sum(incrementda), sum(incrementdb), sum(incrementds), sum(incrementdt_11),  sum(incrementdt_12),  sum(incrementdt_22))/T
   }
   optim(pars,func1, gr,  method="BFGS",control=list(fnscale=-1,parscale=c(1/n,1/n,1/n,1/n, 1/n, 1/n)))
}





dtau2<-function(y1,y2, g1,g2, x1,x2,alpha, beta,sigma2,tau2){
   st<-sigma2+tau2
   r1<-y1-alpha-beta*x1
   r2<-y2-alpha-beta*x2
   iftau<-ifelse(g1==g2,tau2,0)
   ifistau<-1*(g1==g2)
   det<-st^2-iftau^2
   ddet<-2*(st-iftau)
   (-( (r1*r1-2*r1*r2*ifistau+r2*r2)/2/det - (r1*r1*st-2*r1*r2*iftau+r2*r2*st)*(ddet)/2/det/det) -ddet/det/2)
}	







