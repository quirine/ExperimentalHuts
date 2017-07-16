# script to simulate different scenarios and run mcmc on those 
# takes 2 args: one for which data to simulate, one for starting values
source('Parameters.LTFU.SR.Exit.R')
source('SR.Exit.MCMC.Functions.R')
args <- commandArgs(trailingOnly = TRUE)
set.seed(85)
rm(defaults)
RealValues <- sample.from.uniform(1000,LTFU=TRUE,movement.tied.to.exit)
defaults <- RealValues[as.numeric(args[1]),]
Qs <- RealValues[as.numeric(args[1]),c(1,2,3)]            # q_AorE, q_BorD, q_C: movement rates between huts
Ps <- RealValues[as.numeric(args[1]),c(4,5)]                       # p_BorD, p_C (p_AorE is always 1 and thus not defined) : probability of moving away from SR-hut (C)
Rs <- RealValues[as.numeric(args[1]),c(6,7,8)] #c(0.001, 0.001, 0.001)            # r_AorE, r_BorD, r_C : exit rates
Ks <- RealValues[as.numeric(args[1]),c(9,10,11)]         # k_AorE, k_BorD, k_C : knockdown rates
u <- RealValues[as.numeric(args[1]),12]

# Mosquito.vector<-runif(1000,50,1000)
# Mosquitoes <- Mosquito.vector[as.numeric(args[1])]                      # Number of mosquitoes released per hut
Mosquitoes <- 25
set.seed(89) # just to make sure we start with different initial conditions
Starters <- sample.from.uniform(1000,LTFU=TRUE,movement.tied.to.exit)
startvalue <- Starters[as.numeric(args[2]),]
sd = c(0.4, 0.4, 0.4, 0.4, 0.4, 0.3, 0.3, 0.4, 1, 1, 1, 0.1)
relative.sds = TRUE
iter = 1e5                          # number of iterations
wr.freq = round(iter/1000)              # how frequently would you like the parameters to be printed to file?
p.prior = 'Beta'
q.prior = 'Gamma'
name =  'simulated.data.ltfu.RelSD.25_'
directory.name = paste(name,Sys.Date(),sep='')
dir.create(directory.name) 
filename = name
mcmc.output.file = paste(directory.name,'/',filename,sprintf("%02d",as.numeric(args[1])),'_',sprintf("%02d",as.numeric(args[2])),'_',Sys.Date(),'.csv',sep='') 
source('Main.SimulateData.LTFU.SR.Exit.R')
save(list = c('sd','defaults','Starters','draw.probs','relative.sds','adaptive','p.prior','q.prior','Mosquitoes','Sim.Data', 'Sim.Data.end', 'name','directory.name'),file = paste(directory.name,'/',filename,'_',sprintf("%02d",as.numeric(args[1])),'.RData',sep=''))
source('Main.RunMCMC.LTFU.SR.Exit.R')