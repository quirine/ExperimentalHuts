
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


# load data ---------------------------------------------------------------
# load('ExitData.ltfu.RData')

# Prepare the data  ---------------------------------------------------------------
# dose = 0
Ind <- which(Data$Dose == dose)    # for now, just looking at one dosage. add '.dose' to nll-fuction to assess dose-effects
Data <- Data[Ind,]

attach(Data)
K_exit = cbind(k.ex.A, k.ex.B, k.ex.C, k.ex.D, k.ex.E )
K_kd = cbind(k.kd.A, k.kd.B, k.kd.C, k.kd.D, k.kd.E )
k_huts = k.bb
interval = diff(Data$Time)[1]
Init <- diag(5) # create 5 sets of initial conditions, one for each release location
Init <- cbind(Init,matrix(0,5,10))
Rel_loc = factor(Data$Rel_loc)
Locations <- Data$Rel_loc

attach(Data.end)
K_exit.end = cbind(k.ex.A.end, k.ex.B.end, k.ex.C.end, k.ex.D.end, k.ex.E.end )
K_kd.end = cbind(k.kd.A.end, k.kd.B.end, k.kd.C.end, k.kd.D.end, k.kd.E.end )
k_u = k.u
K_huts = cbind(k.hut.A, k.hut.B, k.hut.C, k.hut.D, k.hut.E )
Rel_loc.end = Rel_loc.end
end.time = Data.end$Time.end[1]
# kd.interval.factor = 2

# run mcmc ----------------------------------------------------------------
# for local use set e.g.:
# iter = 50
# wr.freq = iter +1
# tic()
chain = run_metropolis_MCMC(startvalue = startvalue,iterations = iter,t = Data$Time, draw.probabilities = draw.probs, 
                            filename = mcmc.output.file, writing.freq = wr.freq, 
                            Qfix = FALSE,sd = sd,one.at.a.time = TRUE, relative.sds = relative.sds, 
                            adaptive = adaptive, adapt.par = adaptive.parameters,kd.interval.factor = kd.interval.factor,
                            LTFU = LTFU, debugging.plot = FALSE, p.prior = p.prior, q.prior = q.prior, dose = dose,
                            first.measurement = first.measurement, movement.tied.to.exit = movement.tied.to.exit)
# toc()
# chain = mcmc(chain)
# burnIn = 1
# chain = chain[-1,]; chain = chain[complete.cases(chain),]; chain = mcmc(chain); 
# acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
# acceptance
# plot(chain)









