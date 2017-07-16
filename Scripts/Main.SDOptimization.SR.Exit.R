# script to explore optimal values for standard deviation

# Load packages -----------------------------------------------------------
library(deSolve)
library(abind)
library(tictoc)
library(coda)

# # Set work directory ------------------------------------------------------
# setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# Source files ------------------------------------------------------------
source('Parameters.SR.Exit.R')
source('SR.Exit.Models.R')
source('SR.Exit.MCMC.Functions.R')

# set parameters ----------------------------------------------------------
SD <- c(0.001, 0.01,0.1)
iterations = 1000
# explore sd's ----------------------------------------------------------------
Acceptance <- numeric()
for ( i in 1:length(SD) ) {
  tic()
  chain = run_metropolis_MCMC(startvalue,iter = iterations,t = Sim.Data$Time, filename = mcmc.output.file, writing.freq = wr.freq,fix=fix, Qfix = TRUE,SD[i])
  toc()
  burnIn = 0
  acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
  Acceptance = c(Acceptance, acceptance)
}

# plot --------------------------------------------------------------------

plot(SD,Acceptance)

