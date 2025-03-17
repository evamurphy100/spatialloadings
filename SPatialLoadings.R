library(nimble)
library(coda)
library(dplyr)
############################################################
###########             Load data               ############
############################################################
library(readxl)
## Read in the NC data
NCdata <- read_excel("NCdata.xlsx")
### Use only 2017 or later data
outcomedataNC = NCdata[which(NCdata$year>=2017),]
#outcomedataNC = outcomedataNC %>%
  #filter(year == 2019)

############################################################
###########        Extract the outcomes         ############
############################################################
y.D = outcomedataNC$`Illicit opioid overdose deaths`
y.E = outcomedataNC$`Drug overdose ED visits`
Pop = outcomedataNC$pop
Year = outcomedataNC$year
y.C = outcomedataNC$AcuteHepCpdf + outcomedataNC$ChronicHepC
y.I = outcomedataNC$HIV
y.T = outcomedataNC$`People served by treatment programs`
y.B = outcomedataNC$`Patients receiving buprenorphine`

############################################################
###########           Address census            ############
############################################################
cens = 1*(is.na(outcomedataNC$`People served by treatment programs`))
cenlb = rep(NA,length(y.T))
cenub = rep(NA,length(y.T))
cenlb[which(cens==1)]= 1
cenub[which(cens==1)]= 5

cenlb[which(is.na(cenlb))]=y.T[which(is.na(cenlb))]
cenub[which(is.na(cenub))]=y.T[which(is.na(cenub))]

cen=cens
bd=cbind(cenlb,cenub)

### plug in values for censored right now
### plug in values for censored right now
### plug in values for censored right now
y.T0 = y.T
y.T0[which(cens == 1)] = 5
## Load the neighborhood matrix
load("Adj_Matrix_NC.RData")
adj.mat = adj.mat
n<-nrow(adj.mat)
TT = length(unique(outcomedataNC$year))
num<-colSums(adj.mat)


# Expected count standardized to the state level rate 2017
ET = rep(0,n*TT)
EB = rep(0,n*TT)
ED = rep(0,n*TT)
EE = rep(0,n*TT)
EC = rep(0,n*TT)
EI = rep(0,n*TT)
for(t in 1:TT){
  for(i in 1:n){
    ET[((t-1)*n+i)] =  sum(y.T0[which(outcomedataNC$year==2017)])/sum(Pop[which(outcomedataNC$year==2017)])*Pop[((t-1)*n+i)]
    ED[((t-1)*n+i)] =  sum(y.D[which(outcomedataNC$year==2017)])/sum(Pop[which(outcomedataNC$year==2017)])*Pop[((t-1)*n+i)]
    EB[((t-1)*n+i)] =  sum(y.B[which(outcomedataNC$year==2017)])/sum(Pop[which(outcomedataNC$year==2017)])*Pop[((t-1)*n+i)]
    EE[((t-1)*n+i)] =  sum(y.E[which(outcomedataNC$year==2017)])/sum(Pop[which(outcomedataNC$year==2017)])*Pop[((t-1)*n+i)]
    EC[((t-1)*n+i)] =  sum(y.C[which(outcomedataNC$year==2017)])/sum(Pop[which(outcomedataNC$year==2017)])*Pop[((t-1)*n+i)]
    EI[((t-1)*n+i)] =  sum(y.I[which(outcomedataNC$year==2017)])/sum(Pop[which(outcomedataNC$year==2017)])*Pop[((t-1)*n+i)]
  }
}




############################################################
###########         Set up for NIMBLE          ############
############################################################

adj<-NULL
for(j in 1:n){
  adj<-c(adj,which(adj.mat[j,]==1)) ##takes only the neighbors
}
adj<-as.vector(adj)
num<-as.vector(num)
weights<-1+0*adj


model_code=nimbleCode({
  for (t in T0:T0){
    for (i in 1:n){
    y.D[n*(t-T0)+i] ~ dpois(ED[n*(t-T0)+i]*lambdaD[n*(t-T0)+i])
    cen[n*(t-T0)+i] ~ dinterval(y.T[n*(t-T0)+i],bd[n*(t-T0)+i,1:2])
    y.T[n*(t-T0)+i] ~ dpois(ET[n*(t-T0)+i]*lambdaT[n*(t-T0)+i])
    y.E[n*(t-T0)+i] ~ dpois(EE[n*(t-T0)+i]*lambdaE[n*(t-T0)+i])
    y.B[n*(t-T0)+i] ~ dpois(EB[n*(t-T0)+i]*lambdaB[n*(t-T0)+i])
    y.C[n*(t-T0)+i] ~ dpois(EC[n*(t-T0)+i]*lambdaC[n*(t-T0)+i])
    y.I[n*(t-T0)+i] ~ dpois(EI[n*(t-T0)+i]*lambdaI[n*(t-T0)+i])
    log(lambdaD[n*(t-T0)+i]) <- load.shiftedD[i]*(U[n*(t-T0)+i]+mu[n*(t-T0)+i]) + VD[n*(t-T0)+i]
    log(lambdaT[n*(t-T0)+i]) <- load.shiftedT[i]*(U[n*(t-T0)+i]+mu[n*(t-T0)+i]) + VT[n*(t-T0)+i]
    log(lambdaE[n*(t-T0)+i]) <- load.shiftedE[i]*(U[n*(t-T0)+i]+mu[n*(t-T0)+i]) + VE[n*(t-T0)+i]
    log(lambdaB[n*(t-T0)+i]) <- load.shiftedB[i]*(U[n*(t-T0)+i]+mu[n*(t-T0)+i]) + VB[n*(t-T0)+i]
    log(lambdaC[n*(t-T0)+i]) <- load.shiftedC[i]*(U[n*(t-T0)+i]+mu[n*(t-T0)+i]) + VC[n*(t-T0)+i]
    log(lambdaI[n*(t-T0)+i]) <- load.shiftedI[i]*(U[n*(t-T0)+i]+mu[n*(t-T0)+i]) + VI[n*(t-T0)+i]
    VD[n*(t-T0)+i] ~ dnorm(0,tau.VD) 
    VT[n*(t-T0)+i] ~ dnorm(0,tau.VT)
    VE[n*(t-T0)+i] ~ dnorm(0,tau.VE)
    VB[n*(t-T0)+i] ~ dnorm(0,tau.VB)
    VC[n*(t-T0)+i] ~ dnorm(0,tau.VC) 
    VI[n*(t-T0)+i] ~ dnorm(0,tau.VI) 
    mu[n*(t-T0)+i] <- mu0[t-T0+1]
  }
  
  # ICAR prior for the spatial factors
  U[(n*(t-T0)+1):(n*(t-T0)+n)] ~ dcar_normal(adj[], weights[], num[], tau = tau.U[t-(T0-1)], zero_mean=1)
  tau.U[t-(T0-1)] ~ dgamma(0.5,0.5)
  
  }
  
  for (t in (T0+1):(T0+TT-1)){
    for (i in 1:n){
      y.D[n*(t-T0)+i] ~ dpois(ED[n*(t-T0)+i]*lambdaD[n*(t-T0)+i])
      cen[n*(t-T0)+i] ~ dinterval(y.T[n*(t-T0)+i],bd[n*(t-T0)+i,1:2])
      y.T[n*(t-T0)+i] ~ dpois(ET[n*(t-T0)+i]*lambdaT[n*(t-T0)+i])
      y.E[n*(t-T0)+i] ~ dpois(EE[n*(t-T0)+i]*lambdaE[n*(t-T0)+i])
      y.B[n*(t-T0)+i] ~ dpois(EB[n*(t-T0)+i]*lambdaB[n*(t-T0)+i])
      y.C[n*(t-T0)+i] ~ dpois(EC[n*(t-T0)+i]*lambdaC[n*(t-T0)+i])
      y.I[n*(t-T0)+i] ~ dpois(EI[n*(t-T0)+i]*lambdaI[n*(t-T0)+i])
      log(lambdaD[n*(t-T0)+i]) <- load.shiftedD[i]*(U[n*(t-T0)+i]+mu[n*(t-T0)+i]) + VD[n*(t-T0)+i]
      log(lambdaT[n*(t-T0)+i]) <- load.shiftedT[i]*(U[n*(t-T0)+i]+mu[n*(t-T0)+i]) + VT[n*(t-T0)+i]
      log(lambdaE[n*(t-T0)+i]) <- load.shiftedE[i]*(U[n*(t-T0)+i]+mu[n*(t-T0)+i]) + VE[n*(t-T0)+i]
      log(lambdaB[n*(t-T0)+i]) <- load.shiftedB[i]*(U[n*(t-T0)+i]+mu[n*(t-T0)+i]) + VB[n*(t-T0)+i]
      log(lambdaC[n*(t-T0)+i]) <- load.shiftedC[i]*(U[n*(t-T0)+i]+mu[n*(t-T0)+i]) + VC[n*(t-T0)+i]
      log(lambdaI[n*(t-T0)+i]) <- load.shiftedI[i]*(U[n*(t-T0)+i]+mu[n*(t-T0)+i]) + VI[n*(t-T0)+i]
      VD[n*(t-T0)+i] ~ dnorm(0,tau.VD) 
      VT[n*(t-T0)+i] ~ dnorm(0,tau.VT)
      VE[n*(t-T0)+i] ~ dnorm(0,tau.VE)
      VB[n*(t-T0)+i] ~ dnorm(0,tau.VB)
      VC[n*(t-T0)+i] ~ dnorm(0,tau.VC) 
      VI[n*(t-T0)+i] ~ dnorm(0,tau.VI) 
      mu[n*(t-T0)+i] <- mu0[t-T0+1]+eta*(U[n*(t-(T0+1))+i]-mu0[t-T0])
    }
    
    # ICAR prior for the spatial factors
    U[(n*(t-T0)+1):(n*(t-T0)+n)] ~ dcar_normal(adj[], weights[], num[], tau = tau.U[t-(T0-1)], zero_mean=1)
    tau.U[t-(T0-1)] ~ dgamma(0.5,0.5)
    
    
  }
  
  
  #ICAR priors for the spatial loadings
  load.matD[1:n] ~ dcar_normal(adj[], weights[], num[], tau = 1, zero_mean=0)
  load.matE[1:n] ~ dcar_normal(adj[], weights[], num[], tau = 1, zero_mean=0)
  load.matB[1:n] ~ dcar_normal(adj[], weights[], num[], tau = 1, zero_mean=0)
  load.matC[1:n] ~ dcar_normal(adj[], weights[], num[], tau = 1, zero_mean=0)
  load.matI[1:n] ~ dcar_normal(adj[], weights[], num[], tau = 1, zero_mean=0)
  
  load.shiftedD[1:n] <- rep(1,n) + load.matD[1:n]
  load.shiftedE[1:n] <- rep(1,n) + load.matE[1:n]
  load.shiftedB[1:n] <- rep(1,n) + load.matB[1:n]
  load.shiftedC[1:n] <- rep(1,n) + load.matC[1:n]
  load.shiftedI[1:n] <- rep(1,n) + load.matI[1:n]
  
  for(t in T0:(T0+TT-1)){
    mu0[t-T0+1] ~ dflat()
  }
  bD ~ dflat()
  bT ~ dflat()
  bE ~ dflat()
  bB ~ dflat()
  bC ~ dflat()
  bI ~ dflat()
  
  tau.VD ~ dgamma(0.5,0.5)
  tau.VT ~ dgamma(0.5,0.5)
  tau.VE ~ dgamma(0.5,0.5)
  tau.VC ~ dgamma(0.5,0.5)
  tau.VI ~ dgamma(0.5,0.5)
  tau.VB ~ dgamma(0.5,0.5)
  eta ~ dunif(0,1)
  
})

############################################################
#########              Call NIMBLE                 ###########
############################################################

mod_constants=list(bd=bd,num=num,adj=adj,weights=weights,n=n,TT=TT,T0=min(Year))

mod_data=list(cen=cen,y.D=y.D,y.T=y.T,y.E=y.E,y.B=y.B,y.C=y.C,y.I=y.I,
              ET=ET,EE=EE,ED=ED,EB=EB,EC=EC,EI=EI)

mod_inits=list(U=rep(0,n*TT),load.shiftedT = rep(1,n),load.matD=rep(0,n),load.matE=rep(0,n),load.matB=rep(0,n),load.matC=rep(0,n),load.matI=rep(0,n),
               load.shiftedD=rep(0,n),load.shiftedE=rep(0,n),load.shiftedB=rep(0,n),load.shiftedC=rep(0,n),
               load.shiftedI=rep(0,n),mu0=rep(0,TT),
               VD=rep(0,n*TT),VE=rep(0,n*TT),VT=rep(0,n*TT),VB=rep(0,n*TT),VC=rep(0,n*TT),VI=rep(0,n*TT),
               y.T=floor(.5*(bd[,1]+bd[,2])))

# Build the model.
nim_model <- nimbleModel(model_code,mod_constants,mod_data,mod_inits)
compiled_model <- compileNimble(nim_model,resetFunctions = TRUE)

# Set up samplers.
mcmc_conf <- configureMCMC(nim_model,monitors=c("mu","load.shiftedD","load.shiftedE",
                                                "load.shiftedB","load.shiftedT","load.shiftedC","load.shiftedI",
                                                "U","VD","VT","VE","VB","VC","VI", "lambdaT","lambdaD","lambdaE","lambdaB",
                                                "lambdaC","lambdaI","eta"),
                           useConjugacy = TRUE)

mod_mcmc<-buildMCMC(mcmc_conf)
compiled_mcmc<-compileNimble(mod_mcmc, project = nim_model,resetFunctions = TRUE)


# Run the model 
MCS=500000
st<-Sys.time()
samples=runMCMC(compiled_mcmc,inits=mod_inits,
                nchains = 1, nburnin=floor(MCS/2),niter = MCS,samplesAsCodaMCMC = TRUE,thin=50,
                summary = FALSE, WAIC = FALSE,progressBar=TRUE) 
Sys.time()-st

save(samples,file="ModOutput500000_Spatial_Loadings.Rda")
