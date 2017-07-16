# script to sweep over all parameters and evaluate likelihood

rm(list=ls())

# Load packages -----------------------------------------------------------
library(deSolve)
library(abind)
library(coda)
library(pomp)

# Set work directory ------------------------------------------------------
setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# Set directory to save figures -------------------------------------------
path <- '../Figures/Simulated/'
cd <- getwd()

# Load source files  ---------------------------------------
source('Parameters.LTFU.SR.Exit.R')
source('SR.Exit.MCMC.Functions.R')
source('SR.Exit.Models.R')

# set parameters  ---------------------------------------------------------
n = 150                   # number of runs per parameter sweep
Mosquitoes = 1000
num.exp.days = 5
kd.interval.factor = 2
num.sim = 100             # number of different simulations
Names=c('movement 2 away','movement 1 away','movement SR', 
        'prop away from SR', 
        'exit 2 away', 'exit 1 away', 'exit SR hut',
        'KD 2 away', 'KD 1 away', 'KD SR hut', 'loss-to-follow-up' )

# Sweep across data sets --------------------------------------------------
Defaults =  exp(sample.sobol(num.sim,LTFU=TRUE,movement.tied.to.exit))
someData <- rep(NaN, (n*dim(Defaults)[2]));  
LL.50 <- array(someData,c(n,dim(Defaults)[2]))
Starters =  exp(sample.sobol(n,LTFU=TRUE,movement.tied.to.exit))

for (nn in 1:num.sim){
  defaults = as.numeric(Defaults[nn,])
  source('Main.SimulateData.LTFU.SR.Exit.R')
  interval = diff(Sim.Data$Time)[1]
  print('Data simulated')
  Sim.Data.org <- Sim.Data
  Sim.Data.end.org <- Sim.Data.end
  LL.50.temp <- LL.10.temp <- LL.5.temp <- matrix(data=NA, nrow = n, ncol = length(defaults))
  for (ii in 1:length(defaults)) {
    test.param = ii
    print(paste('next parameter',ii))
    for (jj in 1:n){
      temp = as.numeric(Starters[jj,])
      temp[-test.param] = as.numeric(defaults[-test.param])
      retrieve.and.attach(Data=Sim.Data, End = FALSE)
      retrieve.and.attach(Data=Sim.Data.end, End = TRUE)
      LL.50.temp[jj,ii] = likelihood.hut.movement.ltfu(params = temp,
                                                       t = Sim.Data$Time,interval = interval,K_exit = K_exit,K_kd = K_kd,k_huts = k_huts,
                                                       K_exit.end = K_exit.end,K_kd.end = K_kd.end,k_u = k_u,
                                                       K_huts = K_huts,Rel_loc = Rel_loc,Rel_loc.end = Rel_loc.end,end.time = end.time,
                                                       delta = delta,kd.interval.factor = kd.interval.factor, first.measurement = 30,movement.tied.to.exit)

      Sim.Data <- Sim.Data.org
    }
  }
  LL.50 = abind(LL.50,LL.50.temp,along = 3); 
}
LL.50 <- LL.50[,,-1]

save(list = ls(), file = 'workspace.parametersweeps.1000.RData')
# Figure of likelihoods------------------------------------------------------------------
setwd(path)

pdf('LL.sweep.pdf',width=4.75,height=4.75)
par(mfrow = c(4,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  #mar=c(1,3,1,1.9)

for (ii in 1:12){
  ind = complete.cases(LL.50[,ii,20])
  plot(Starters[ind,ii], LL.50[ind,ii,20])
  abline(v=Defaults[20,ii]) 
  lines(lowess(Starters[ind,ii],LL.50[ind,ii,20]), col='red', lwd = 2)
  mtext(side = 1, Names[ii], line = 2.25, cex =.8)
}
mtext(side = 2, 'LL', line = 2.25, cex = .8, outer = TRUE)

dev.off()
setwd(cd)

# Figure of likelihoods------------------------------------------------------------------
setwd(path)

pdf('LL.sweep.RealvsEstimated.pdf',width=4.75,height=4.75)
par(mfrow = c(4,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  #mar=c(1,3,1,1.9)

for (ii in 1:12) {
  max = which(LL.50[ind,ii,1] == max(LL.50[ind,ii,1]))
  if (length(max) ==1){
    plot(Defaults[1,ii],Starters[max,ii],xlim=range(Defaults[,ii]), ylim = range(Defaults[,ii]))
    lines(Defaults[,ii], Defaults[,ii], col = 'red', lwd = 2)
    mtext(side = 3, Names[ii], line = 1, cex =.8)
  }
  for (jj in 2:num.sim){
    max = which(LL.50[ind,ii,jj] == max(LL.50[ind,ii,jj]))
    if (length(max) ==1){
      points(Defaults[jj,ii],Starters[max,ii])
    }
  }
}
mtext(side = 1, text = 'Real values', line = 2.25, cex = .8, outer = TRUE)  
mtext(side = 2, text = 'Estimated values', line = 2.25, cex = .8, outer = TRUE)      

dev.off()
setwd(cd)


