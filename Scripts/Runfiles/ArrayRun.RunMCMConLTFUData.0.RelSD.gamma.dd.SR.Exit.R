rm(list=ls())

# Set work directory ------------------------------------------------------
# setwd("~/Box Sync/Research Backup/Research/Iquitos/Analysis/Movement/MarkovChains/Scripts")

# Run data fitting ---
load('ExitData.ltfu.RData')
source('Parameters.LTFU.SR.Exit.R')
source('SR.Exit.MCMC.Functions.R')

args <- commandArgs(trailingOnly = TRUE)
Starters <- sample.from.uniform(1000, LTFU = TRUE)
startvalue <- Starters[as.numeric(args[1]),]
sd = c(0.5, 0.5, 0.8, 0.4, 0.4, 0.4, 0.5, 0.6, 1, 1, 1, 0.4)
relative.sds = TRUE
p.prior = 'Dose.dependent'
q.prior = 'Gamma'
dose = 0
iter = 1e5                          # number of iterations
wr.freq = round(iter/1000)              # how frequently would you like the parameters to be printed to file?
name =  'output.ar.ltfudata.relsd.gamma.dd_'
directory.name = paste(name,Sys.Date(),sep='')
dir.create(directory.name) 
filename = paste(name,'dose_',dose,'_',sep='')  
if (as.numeric(args[1]) <10) {
  filename = paste(filename,'0',sep='')
}
mcmc.output.file = paste(directory.name,'/',filename,args[1],'_',Sys.Date(),'.csv',sep='')   
save(list = c('sd','startvalue','draw.probs'),file = paste(directory.name,'/',filename,dose,'_',Sys.Date(),'.RData',sep='') )
source('Main.RunMCMConLTFUData.SR.Exit.R')
