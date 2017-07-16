# script to sweep over all parameters and evaluate likelihood

rm(list=ls())

# Load packages -----------------------------------------------------------
library(deSolve)
library(abind)
# library(tictoc)
library(coda)
library(pomp)

# Set work directory ------------------------------------------------------
setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# Set directory to save figures -------------------------------------------
path <- '../Figures'
cd <- getwd()

# set parameters  ---------------------------------------------------------
n = 15

# Load source files + simulate data ---------------------------------------
source('Parameters.LTFU.SR.Exit.R')
Mosquitoes = 1000
num.exp.days = 5
kd.interval.factor = 2
source('Main.SimulateData.LTFU.SR.Exit.R')
source('SR.Exit.MCMC.Functions.R')
source('SR.Exit.Models.R')

# one at a time (ltfu) -----------------------------------------------------------

test.param = 9
Starters =  exp(sample.sobol(n,LTFU=TRUE,movement.tied.to.exit))
Starters

LL = numeric()
for (jj in 1:n){
  print(jj)
  temp = as.numeric(Starters[jj,])
  temp[-test.param] = as.numeric(defaults[-test.param])
  LL = c(LL,likelihood.hut.movement.ltfu(params = temp,
                                         t = Sim.Data$Time,interval = interval,K_exit = K_exit,K_kd = K_kd,k_huts = k_huts,
                                         K_exit.end = K_exit.end,K_kd.end = K_kd.end,k_u = k_u,
                                         K_huts = K_huts,Rel_loc = Rel_loc,Rel_loc.end = Rel_loc.end,end.time = end.time,
                                         delta = delta,kd.interval.factor = kd.interval.factor, first.measurement = 30,movement.tied.to.exit))
}

ind = complete.cases(LL)
plot((Starters[ind,test.param]),LL[ind],col='red')
abline(v=(as.numeric(defaults[test.param])),col='blue')



# idem without ltfu -------------------------------------------------------
source('Parameters.SR.Exit.R')
kd.interval.factor = 1
source('Main.SimulateData.SR.Exit.R')
K_exit = cbind(k.ex.A, k.ex.B, k.ex.C, k.ex.D, k.ex.E )
K_kd = cbind(k.kd.A, k.kd.B, k.kd.C, k.kd.D, k.kd.E )
k_huts = k.hut
interval = diff(Sim.Data$Time)[1]
Init=Init 
Rel_loc = factor(Sim.Data$Rel_loc)
Locations <- Sim.Data$Rel_loc


# one at a time -----------------------------------------------------------
test.param = 11
Starters =  exp(sample.sobol(n,LTFU=FALSE))

LL = numeric()
for (jj in 1:n){
  print(jj)
  temp = as.numeric(Starters[jj,])
  temp[-test.param] = as.numeric(defaults[-test.param])
  LL = c(LL,likelihood.hut.movement(params = temp,
                                    t = Sim.Data$Time,interval,K_exit,K_kd,k_huts,Rel_loc,delta, kd.interval.factor = kd.interval.factor))
}

ind = complete.cases(LL)
plot((Starters[ind,test.param]),LL[ind])
abline(v=(as.numeric(defaults[test.param])),col='red')

# plot likelihoods --------------------------------------------------------
setwd(path)

Names=rownames(Starters)

pdf('LL.sweep.pdf',width=4.75,height=4.75)
par(mfrow = c(4,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  #mar=c(1,3,1,1.9)

for (ii in 1:12){
  ind = complete.cases(LL)
  plot(log(Starters[ind,ii]),LL[ind])
  abline(v=log(as.numeric(defaults[ii])))
  lines(lowess(log(Starters[ind,ii]),LL[ind]), col='red', lwd = 2)
  mtext(side = 1, Names[ii], line = 1, cex = 0.5)
  mtext(side = 2, LL, line = 1, cex = 0.5)
}

dev.off()
setwd(cd)

