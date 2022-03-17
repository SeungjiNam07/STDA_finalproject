#### STDA Final Project
#### Seung ji Nam

# load library
library(sf)
library(classInt)
library(maps)
library(sp)
library(classInt)
library(gstat)
library(fields)
library(stringr)
library(mapdata)
library(mapproj)
library(rgdal)
library(nlme)
library(batchmeans)
library(nimble)
library(mvtnorm)
library(spacetime)


# set your directory
directory = '/Users/seungji/Library/Mobile Documents/com~apple~CloudDocs/STDA/Project/'
sensor = readOGR(paste0(directory, '/data/censor.shp'))
sensor_oct = readOGR(paste0(directory, '/data/censor_oct.shp'))
sensor_march = readOGR(paste0(directory, '/data/censor_march.shp'))
sensor_mt = readOGR(paste0(directory, 'data/sensor_mt.shp'))
sensor_wk = readOGR(paste0(directory, 'data/sensor_wk.shp'))
################################################################################

# Function Save 
load(paste0(directory, 'modeling/function.RData'))
################################################################################
# 1. scale data
sensor_df = data.frame(sensor_march)
data_scale = sensor_df
col = c('subway','bus','building_c','building_a','pop')
for (c in col){
  data_scale[,c] = sensor_df[,c]/(max(sensor_df[,c])-min(sensor_df[,c]))
}

################################################################################
# 2. Split into Model + Cross-Validation
data_all = data_scale

lat = data_all$lat
lon = data_all$lon

ncv = round(nrow(data_all)*0.2)
n = nrow(data_all)-ncv

set.seed(730)
cv_idx = sample(1:nrow(data_all), ncv,replace = F)
idx = as.numeric(row.names(data_all[-cv_idx,]))

data = data_all[idx,]
datacv = data_all[cv_idx,]
row.names(data) = as.numeric(1:n)
row.names(datacv) = as.numeric(1:ncv)

gridLocation<-cbind(lat[idx],lon[idx])
CVgridLocation<-cbind(lat[cv_idx], lon[cv_idx])
comboLocation<-rbind(gridLocation,CVgridLocation)
distMatFull<-as.matrix(rdist(comboLocation))

# Create Indices
distMatMod<-distMatFull[idx,idx]

# Covariates
XMat<-data[, c('bus','subway','building_c','building_a','pop')]

dim(XMat)
length(beta)
XMatCV<- datacv[, c('bus','subway','building_c','building_a','pop')]

dim(XMatCV)

beta = c(1,1,1,1,1)
phi = 0.1
sigma2 = 0.01

# beta = c(1,1)
XB<-as.matrix(XMat)%*%as.matrix(beta)
cvXB<-as.matrix(XMatCV)%*%as.matrix(beta)
XBFull<-rbind(XB,cvXB)

# Covariance Matrix
CovMat<-sigma2*expcov(distMatFull,phi)

# Latent Gaussian Random Field
gpWFull <- as.numeric(rmvnorm(n=1,mean=rep(0,nrow(CovMat)),sigma = CovMat,method = "chol"))

pWFullPois<-exp(gpWFull+XBFull)

# Observations
obsFullPois<-sapply(pWFullPois,rpois,n=1)

################################################################################
# 3. Model Sample
# Latent Process
pWModPois<-pWFullPois[idx]
obsModPois<-obsFullPois[idx]
pWCVPois<-pWFullPois[cv_idx]
obsCVPois<-obsFullPois[cv_idx]
truthCVMSPEPois<-mean((pWCVPois-obsCVPois)^2)

################################################################################
# 4. Nimble Model
model_string <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    lambda[i] <- exp(W[i]+XB[i])
    Z[i] ~ dpois(lambda[i])
  }
  
  # Constant and Cov Matrix
  XB[1:n]<-beta1*X[,1] + beta2*X[,2]+ beta3*X[,3]+ beta4*X[,4]+ beta5*X[,5]
  covMat[1:n,1:n]<- expcov(dists[1:n,1:n],phi)
  fullCovMat[1:n,1:n]<- sigma2*covMat[1:n,1:n]
  
  # Process Model
  W[1:n] ~ dmnorm(mean = mn[1:n], cov = fullCovMat[1:n,1:n])
  
  # Parameter Model
  sd ~ dinvgamma(0.2,0.2)
  sigma2   ~  dinvgamma(0.2, 0.2)
  phi   ~  dunif(0,1)
  beta1 ~  dnorm(0, sd=sqrt(10))
  beta2 ~  dnorm(0, sd=sqrt(10))
  beta3 ~  dnorm(0, sd=sqrt(10))
  beta4 ~  dnorm(0, sd=sqrt(10))
  beta5 ~  dnorm(0, sd=sqrt(10))
})

niter=200000
consts   <- list(n=n,X=XMat,dists=distMatMod,mn=rep(0,n))
data     <- list(Z=obsModPois)
inits    <- list(beta1=rnorm(1),beta2=rnorm(1),beta3=rnorm(1),beta4=rnorm(1),beta4~rnorm(1),phi=0.5,sigma2=2, 
                 W=rnorm(n))

################################################################################
# 5. Run MCMC
pt<-proc.time()
samples_march  <- nimbleMCMC(model_string, data = data, inits = inits,
                             constants=consts,
                             monitors = c("beta1", "beta2","beta3", "beta4",'beta5',"phi","sigma2","W"),
                             samplesAsCodaMCMC=TRUE,WAIC=FALSE,summary=FALSE,
                             niter = niter, nburnin = 5000, nchains = 5)
ptFinal<-proc.time()-pt
ptFinal

################################################################################
# 6. Summary
pois_march1 = samples_march$chain1
pois_march2 = samples_march$chain2
pois_march3 = samples_march$chain3
pois_march4 = samples_march$chain4
pois_march5 = samples_march$chain5


# Summary
# Table
want = c("beta1", "beta2","beta3", "beta4",'beta5',"phi","sigma2")
summaryMarch<-list()
summaryMarch[[1]]<-round(summaryFunction(data.frame(pois_march1)[want],
                                         time=ptFinal[3]),3)
# summaryMarch[[2]]<-round(summaryFunction(pois_march1[,1:n],
#                                         time=ptFinal[3]),3)
# summaryMarch<-list()
summaryMarch[[2]]<-round(summaryFunction(data.frame(pois_march2)[want],
                                         time=ptFinal[3]),3)

# summaryMat3<-list()
summaryMarch[[3]]<-round(summaryFunction(data.frame(pois_march3)[want],
                                         time=ptFinal[3]),3)


# summaryMat4<-list()
summaryMarch[[4]]<-round(summaryFunction(data.frame(pois_march4)[want],
                                         time=ptFinal[3]),3)
# summaryMat5<-list()
summaryMarch[[5]]<-round(summaryFunction(data.frame(pois_march5)[want],
                                         time=ptFinal[3]),3)

summaryMarch[[1]]
apply(summaryMarch[[2]],1,mean)
save(summaryMarch,ptFinal,file="Pois_march1_MCMC.RData")


# plot
pdf(file = paste0(directory,"modelingPoisson_march1.pdf"),width=10,height=12.5)
par(mfrow=c(7,3),mar=c(3,3,3,3))
# plot.ts(pois_march1[5000:nrow(pois_march1),75]); abline(h=c(beta,phi,sigma2)[i],col="red",lwd=2)
sampInd<-floor(seq(5000,nrow(pois_march1),length.out = 1000))

for(i in 75:81){
  plot.ts(pois_march1[sampInd,i], main= colnames(pois_march1)[i]); abline(h=c(beta,phi,sigma2)[i],col="red",lwd=2)
  dens <- density(pois_march1[sampInd,i])
  plot(dens, main = colnames(pois_march1)[i]); abline(v=c(beta,phi,sigma2)[i],col="red",lwd=2 )
  acf(pois_march1[sampInd,i])
}
for(i in 75:81){
  plot.ts(pois_march2[sampInd,i], main= colnames(pois_march2)[i]); abline(h=c(beta,phi,sigma2)[i],col="red",lwd=2)
  dens <- density(pois_march2[sampInd,i])
  plot(dens, main =colnames(pois_march1)[i] ); abline(v=c(beta,phi,sigma2)[i],col="red",lwd=2, main = NULL )
  acf(pois_march2[sampInd,i])
}
for(i in 75:81){
  plot.ts(pois_march3[sampInd,i], main= colnames(pois_march3)[i]); abline(h=c(beta,phi,sigma2)[i],col="red",lwd=2)
  dens <- density(pois_march3[sampInd,i])
  plot(dens, main =colnames(pois_march1)[i] ); abline(v=c(beta,phi,sigma2)[i],col="red",lwd=2, main = NULL )
  acf(pois_march3[sampInd,i])
}
for(i in 75:81){
  plot.ts(pois_march4[sampInd,i], main= colnames(pois_march4)[i]); abline(h=c(beta,phi,sigma2)[i],col="red",lwd=2)
  dens <- density(pois_march4[sampInd,i])
  plot(dens, main =colnames(pois_march1)[i] ); abline(v=c(beta,phi,sigma2)[i],col="red",lwd=2, main = NULL )
  acf(pois_march4[sampInd,i])
}
for(i in 75:81){
  plot.ts(pois_march5[sampInd,i], main= colnames(pois_march5)[i]); abline(h=c(beta,phi,sigma2)[i],col="red",lwd=2)
  dens <- density(pois_march5[sampInd,i])
  plot(dens, main =colnames(pois_march1)[i] ); abline(v=c(beta,phi,sigma2)[i],col="red",lwd=2, main = NULL )
  acf(pois_march5[sampInd,i])
}
dev.off()


#1, 3

summarymar= summaryMarch[[1]]
save(summarymar,ptFinal,file=paste0(directory,"Pois_mar_MCMC_final.RData"))





