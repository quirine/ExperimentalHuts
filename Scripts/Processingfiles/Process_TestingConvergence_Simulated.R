rm(list=ls())
# load libraries ----------------------------------------------------------
library(coda)
library(IDPmisc)

# set wd ------------------------------------------------------------------
setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# load source files -------------------------------------------------------
source('SR.Exit.FiguresScripts.R')

# set data directory ------------------------------------------------------
data.path <- "../Output/Simulated/output.sim.ltfu.array_25_RelSD_2016-11-11/"

# Set directory to save figures -------------------------------------------
figure.path <- '../Figures/Simulated'
cd <- getwd()

# specify file names ------------------------------------------------------
setwd(data.path)

list.of.files = dir(pattern = '*.csv')
list.of.data = dir(pattern = '*.RData')

setwd(cd)

# set parameters -----------------------------------------------------
Names=c('movement 2 away','movement 1 away','movement SR', 
        'prop away from SR', 
        'exit 2 away', 'exit 1 away', 'exit SR hut',
        'KD 2 away', 'KD 1 away', 'KD SR hut', 'loss-to-follow-up' )

# load and prepare data ---------------------------------------------------
setwd(data.path)

load(list.of.data[1])

# first, get the appropriate length of the chains (shortest chain determines)
SIZE <- numeric()
for (ii in 1:length(list.of.files)) {
  chain = read.csv(list.of.files[ii],header=TRUE,colClasses = 'numeric')
  if (dim(chain)[2] > 12) { temp = ((dim(chain)[2] - 1) / 2); chain = chain[-1,1:temp]  }
  else {  chain = chain[-1,] }
  print(dim(chain)[1])
  SIZE = c(SIZE,dim(chain)[1])
}
size = min(SIZE)

for (ii in 1:length(list.of.files)){
  chain = read.csv(list.of.files[ii],header=TRUE,colClasses = 'numeric')
  if (dim(chain)[2] > 12) { temp = ((dim(chain)[2] - 1) / 2); chain = chain[-1,1:temp]  }
  else {  chain = chain[-1,] }
  chain = mcmc(chain[1:size,-5])
  name <- paste("chain.raw.", ii, sep = "")
  assign(name, chain)
}
chains = apropos("^chain.raw.")
combinedchains.raw = mcmc.list(chain.raw.1, chain.raw.2, chain.raw.3, chain.raw.4,chain.raw.5)
GEL = gelman.diag(combinedchains.raw,multivariate = FALSE)

setwd(cd)
# create gelman.plots to assess burnin period (control)------------------------------------------------------------
setwd(figure.path)

pdf('GelmanPlots.pdf',width=4.75,height=4.75)
par(mfrow = c(4,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (ii in 1:11) { 
  gelman.plot(combinedchains.raw[,ii],auto.layout = FALSE,ask = FALSE,ylim=c(0,10))
  title = paste(Names[ii])
  mtext(side = 3, Names[ii], line = 1, cex = 0.5)
}
mtext(side = 1, 'Iterations', outer = TRUE)
mtext(side = 2, 'Shrinking factor', outer = TRUE)

dev.off()
setwd(cd)

# Trace plots ---------------------------------------
setwd(figure.path)

pdf('TracePlots.pdf',width=4.75,height=4.75)
par(mfrow = c(4,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (ii in 1:dim(chain)[2]){ 
  temp = c(seq(1,4),seq(6,12))
  traceplot(combinedchains.raw[,ii])
  title = paste(Names[ii])
  abline(h = defaults[temp[ii]],col='red')
  mtext(side = 3, Names[ii], line = 1, cex = 0.5)
}
mtext(side = 1, 'Iterations', outer = TRUE)
mtext(side = 2, 'Rates', outer = TRUE)

dev.off()
setwd(cd)

# Set burnin and derive summary statistics --------------------------------------------
burnIn = 50000

chain.1 = mcmc(chain.raw.1[burnIn:size,]); 
chain.2 = mcmc(chain.raw.2[burnIn:size,])
chain.3 = mcmc(chain.raw.3[burnIn:size,])
chain.4 = mcmc(chain.raw.4[burnIn:size,])
chain.5 = mcmc(chain.raw.5[burnIn:size,])
combinedchains = mcmc.list(chain.1, chain.2, chain.3, chain.4,chain.5)

chain = combinedchains  

output.summary = summary(chain)

# Get acceptance rates ----------------------------------------------------

Acceptance <- matrix(data=NA,nrow=output.summary$nchain,ncol=dim(chain.1)[2])
for (ii in 1:output.summary$nchain) { 
  for (jj in 1:dim(chain.raw.1)[2]){
    name <- paste("chain.", ii, sep = "")
    assign(name, chain)
    acceptance = (1-mean(duplicated(chain[-(1),jj]))) * length(Names)
    print(acceptance)
    Acceptance[ii,jj] = acceptance
  }
}

# Density plots with real values ---------------------------------------
setwd(figure.path)

# Names=rownames(output.summary$statistics)

pdf('DensityPlots_simulated.pdf',width=4.75,height=4.75)
par(mfrow = c(4,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (ii in 1:dim(chain.raw.1)[2]){ 
  densplot(chain[,ii],show.obs = TRUE)
  temp = c(seq(1,4),seq(6,12))
  abline(v=defaults[temp[ii]])
  abline(v=output.summary$quantiles[ii,3],col='red')
  print(output.summary$quantiles[ii,3])
  title = paste(Names[ii])
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  
}
mtext(side = 1, 'Rates', outer = TRUE)
mtext(side = 2, 'Density', outer = TRUE)

dev.off()
setwd(cd)


