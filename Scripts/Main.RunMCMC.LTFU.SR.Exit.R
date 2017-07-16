

# Load packages -----------------------------------------------------------
library(deSolve)
library(abind)
# library(tictoc)
library(coda)
library(MASS)
library(mvtnorm)

# # Set work directory ------------------------------------------------------
# setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# Source files ------------------------------------------------------------
#source('Parameters.LTFU.SR.Exit.R')
source('SR.Exit.Models.R')
source('SR.Exit.MCMC.Functions.R')

# run mcmc ----------------------------------------------------------------
# for local use set e.g.:
# iter = 1500
# wr.freq = iter +1
# Starters = sample.from.uniform(1,LTFU=TRUE,movement.tied.to.exit = TRUE )
# Starters[1,4:8] <- defaults[4:8]
# startvalue <-  Starters[1,]
# # startvalue[6:8]<- c(3.477827e-02, 4.071661e-02, 4.075752e-05)  # set equal to medians from control experiment
# draw.probs = c(rep(1/8,4),rep(0,4),rep(1/8,4))
# # p.prior = 'Dose.dependent'
# # q.prior = 'Gamma'
# # relative.sds = TRUE

# # draw.probs <- c(rep(0,3),1,rep(0,8))
# draw.probs = c(rep(1/11,4),0,rep(1/11,7))
# relative.sds = TRUE
# sd = c(rep(0.1,3),rep(0.1,2),rep(0.1,3),rep(0.1,3),0.1)
# dose = 0.1
# p.prior = 'Flat'
# q.prior = 'Gamma'
# tic()
chain = run_metropolis_MCMC(startvalue = startvalue,iterations = iter,t = Sim.Data$Time, draw.probabilities = draw.probs, 
                            filename = mcmc.output.file, writing.freq = wr.freq, 
                            Qfix = FALSE,sd = sd,one.at.a.time = TRUE, relative.sds = relative.sds, 
                            adaptive = adaptive, adapt.par = adaptive.parameters,kd.interval.factor = kd.interval.factor,
                            LTFU = LTFU, debugging.plot = FALSE, p.prior = p.prior, q.prior = q.prior, dose = dose, 
                            movement.tied.to.exit = movement.tied.to.exit)

# # toc()
# chain = mcmc(chain)
# burnIn = 1
# chain = chain[-1,]; chain = chain[complete.cases(chain),]; chain = mcmc(chain); 
# acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
# acceptance
# 












