rm(list=ls())
# load libraries ----------------------------------------------------------
library(coda)
library(IDPmisc)

# set wd ------------------------------------------------------------------
setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# load source files -------------------------------------------------------
source('SR.Exit.FiguresScripts.R')

# set data directory ------------------------------------------------------
data.path <- "../Output/Final/output.ar.ltfudata.relsd.gamma.beta.fixedr_2016-12-13/"

# Set directory to save figures -------------------------------------------
figure.path <- '../Figures/Final'
cd <- getwd()

# specify file names ------------------------------------------------------
setwd(data.path)

list.of.files = dir(pattern = '*.csv')
output.file.0 = list.of.files[1]
output.file.625 = list.of.files[2]
output.file.125 = list.of.files[3]

setwd(cd)

# set parameters -----------------------------------------------------
burnIn = 5000

Names=c('movement 2 away','movement 1 away','movement SR', 
        'prop away from SR', 'prop L vs R from SR hut',
        'exit 2 away', 'exit 1 away', 'exit SR hut',
        'KD 2 away', 'KD 1 away', 'KD SR hut', 'loss-to-follow-up' )


# load and prepare data ---------------------------------------------------
setwd(data.path)

chain = read.csv(output.file.0,header=TRUE,colClasses = 'numeric')
chain = chain[-1,]
chain.raw.0 = mcmc(chain)
chain = chain[burnIn:dim(chain)[1],]
chain = mcmc(chain,start = burnIn)
acceptance.0 = 1-mean(duplicated(chain.raw.0[-(1:burnIn),])) 
output.summary.0 = summary(chain)
chain.0 = chain

chain = read.csv(output.file.625,header=TRUE,colClasses = 'numeric')
chain = chain[-1,]
chain.raw.625 = mcmc(chain)
chain = chain[burnIn:dim(chain)[1],]
chain = mcmc(chain,start = burnIn)
acceptance.625 = 1-mean(duplicated(chain.raw.625[-(1:burnIn),])) 
output.summary.625 = summary(chain)
chain.625 = chain

chain = read.csv(output.file.125,header=TRUE,colClasses = 'numeric')
chain = chain[-1,]
chain.raw.125 = mcmc(chain)
chain = chain[burnIn:dim(chain)[1],]
chain = mcmc(chain,start = burnIn)
acceptance.125 = 1-mean(duplicated(chain.raw.125[-(1:burnIn),])) 
output.summary.125 = summary(chain)
chain.125 = chain

setwd(cd)

# Density plots with real values - movement rates---------------------------------------
setwd(figure.path)

# Names=rownames(output.summary.0$statistics)

pdf('DensityPlots_movement.pdf',width=4.75,height=4.75)
par(mfrow = c(3,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (ii in 1:3){ 
  densplot(chain.0[,ii],show.obs = TRUE, xlim=c(0,.2))
  abline(v=output.summary.0$quantiles[ii,3],col='red')
  print(output.summary.0$quantiles[ii,3])
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 1) { mtext( side = 2, 'control', line = 2, cex = 0.5) }
 
}
for (ii in 1:3){ 
  
  densplot(chain.625[,ii],show.obs = TRUE, xlim=c(0,.2))
  abline(v=output.summary.625$quantiles[ii,3],col='red')
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 1) { mtext( side = 2, '0.0625 FAR', line = 2, cex = 0.5) }
}
for (ii in 1:3){ 
  densplot(chain.125[,ii],show.obs = TRUE, xlim=c(0,.2))
  abline(v=output.summary.125$quantiles[ii,3],col='red')
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 1) { mtext( side = 2, '0.125 FAR', line = 2, cex = 0.5) }
}
mtext(side = 1, 'Rates', outer = TRUE)
mtext(side = 2, 'Density', outer = TRUE, line = 2)
# 
dev.off()
setwd(cd)

# Density plots with real values : probabilities away from SR ---------------------------------------
setwd(figure.path)

# Names=rownames(output.summary.0$statistics)
ii = 4

pdf('DensityPlots_proportion_away.pdf',width=4.75,height=4.75)
par(mfrow = c(3,1),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1)) 

  densplot(chain.0[,ii],show.obs = TRUE, xlim=c(0,1))
  abline(v=output.summary.0$quantiles[ii,3],col='red')
  print(output.summary.0$quantiles[ii,3])
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  mtext( side = 2, 'control', line = 2, cex = 0.8)
  
  densplot(chain.625[,ii],show.obs = TRUE, xlim=c(0,1))
  abline(v=output.summary.625$quantiles[ii,3],col='red')
  mtext( side = 2, '0.0625 FAR', line = 2, cex = 0.8)

  densplot(chain.125[,ii],show.obs = TRUE, xlim=c(0,1))
  abline(v=output.summary.125$quantiles[ii,3],col='red')
  mtext( side = 2, '0.125 FAR', line = 2, cex = 0.8)

mtext(side = 1, 'Rates', outer = TRUE)
mtext(side = 2, 'Density', outer = TRUE, line = 2)
# 
dev.off()
setwd(cd)

# Density plots with real values : Loss to follow up ---------------------------------------
setwd(figure.path)

# Names=rownames(output.summary.0$statistics)
ii = 12

pdf('DensityPlots_loss_to_follow_up.pdf',width=4.75,height=4.75)
par(mfrow = c(3,1),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

densplot(chain.0[,ii],show.obs = TRUE, xlim=c(0,1e-3))
abline(v=output.summary.0$quantiles[ii,3],col='red')
print(output.summary.0$quantiles[ii,3])
mtext(side = 3, Names[ii], line = 1, cex =0.5)
mtext( side = 2, 'control', line = 2, cex = 0.8)

densplot(chain.625[,ii],show.obs = TRUE, xlim=c(0,1e-3))
abline(v=output.summary.625$quantiles[ii,3],col='red')
mtext( side = 2, '0.0625 FAR', line = 2, cex = 0.8)

densplot(chain.125[,ii],show.obs = TRUE, xlim=c(0,1e-3))
abline(v=output.summary.125$quantiles[ii,3],col='red')
mtext( side = 2, '0.125 FAR', line = 2, cex = 0.8)

mtext(side = 1, 'Rates', outer = TRUE)
mtext(side = 2, 'Density', outer = TRUE, line = 2)
# 
dev.off()
setwd(cd)


# Density plots with real values - exit rates---------------------------------------
setwd(figure.path)

# Names=rownames(output.summary.0$statistics)

pdf('DensityPlots_exit.pdf',width=4.75,height=4.75)
par(mfrow = c(3,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (ii in 6:8){ 
  densplot(chain.0[,ii],show.obs = TRUE, xlim=c(0,0.01))
  abline(v=output.summary.0$quantiles[ii,3],col='red')
  print(output.summary.0$quantiles[ii,3])
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 6) { mtext( side = 2, 'control', line = 2, cex = 0.5) }
  
}
for (ii in 6:8){ 
  
  densplot(chain.625[,ii],show.obs = TRUE, xlim=c(0,0.01))
  abline(v=output.summary.625$quantiles[ii,3],col='red')
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 6) { mtext( side = 2, '0.0625 FAR', line = 2, cex = 0.5) }
}
for (ii in 6:8){ 
  densplot(chain.125[,ii],show.obs = TRUE, xlim=c(0,0.01))
  abline(v=output.summary.125$quantiles[ii,3],col='red')
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 6) { mtext( side = 2, '0.125 FAR', line = 2, cex = 0.5) }
}
mtext(side = 1, 'Rates', outer = TRUE)
mtext(side = 2, 'Density', outer = TRUE, line = 2)

dev.off()
setwd(cd)

# Density plots with real values - kd rates---------------------------------------
setwd(figure.path)

# Names=rownames(output.summary.0$statistics)

pdf('DensityPlots_kd.pdf',width=4.75,height=4.75)
par(mfrow = c(3,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1)) 

for (ii in 9:11){ 
  densplot(chain.0[,ii],show.obs = TRUE, xlim=c(0,0.001))
  abline(v=output.summary.0$quantiles[ii,3],col='red')
  print(output.summary.0$quantiles[ii,3])
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 9) { mtext( side = 2, 'control', line = 2, cex = 0.5) }
  
}
for (ii in 9:11){ 
  
  densplot(chain.625[,ii],show.obs = TRUE, xlim=c(0,0.001))
  abline(v=output.summary.625$quantiles[ii,3],col='red')
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 6) { mtext( side = 2, '0.0625 FAR', line = 2, cex = 0.5) }
}
for (ii in 9:11){ 
  densplot(chain.125[,ii],show.obs = TRUE, xlim=c(0,0.001))
  abline(v=output.summary.125$quantiles[ii,3],col='red')
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 9) { mtext( side = 2, '0.125 FAR', line = 2, cex = 0.5) }
}
mtext(side = 1, 'Rates', outer = TRUE)
mtext(side = 2, 'Density', outer = TRUE, line = 2)

dev.off()
setwd(cd)

# Trace plots ---------------------------------------
setwd(figure.path)

# Names=rownames(output.summary.0$statistics)

pdf('TracePlots_0.pdf',width=4.75,height=4.75)
par(mfrow = c(4,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1)) 

for (ii in 1:dim(chain)[2]){ 
  traceplot(chain.raw.0[,ii])
  title = paste(Names[ii])
  mtext(side = 3, Names[ii], line = 1, cex = 0.5)
  acceptance = (1-mean(duplicated(chain.raw.0[-(1:burnIn),ii]))) * length(Names)   
  print(acceptance)
  mtext(side = 3 ,acceptance, line = 0.2, cex = 0.5)
}
mtext(side = 1, 'Iterations', outer = TRUE)
mtext(side = 2, 'Rates', outer = TRUE)

dev.off()
setwd(cd)

# Trace plots ---------------------------------------
setwd(figure.path)

# Names=rownames(output.summary.0$statistics)

pdf('TracePlots_0625.pdf',width=4.75,height=4.75)
par(mfrow = c(4,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (ii in 1:dim(chain)[2]){ 
  traceplot(chain.raw.625[,ii])
  title = paste(Names[ii])
  mtext(side = 3, Names[ii], line = 1, cex = 0.5)
  acceptance = (1-mean(duplicated(chain.raw.625[-(1:burnIn),ii]))) * length(Names)   
  print(acceptance)
  mtext(side = 3 ,acceptance, line = 0.2, cex = 0.5)
}
mtext(side = 1, 'Iterations', outer = TRUE)
mtext(side = 2, 'Rates', outer = TRUE)

dev.off()
setwd(cd)

# Trace plots ---------------------------------------
setwd(figure.path)

# Names=rownames(output.summary.0$statistics)

pdf('TracePlots_125.pdf',width=4.75,height=4.75)
par(mfrow = c(4,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (ii in 1:dim(chain)[2]){ 
  traceplot(chain.raw.125[,ii])
  title = paste(Names[ii])
  mtext(side = 3, Names[ii], line = 1, cex = 0.5)
  acceptance = (1-mean(duplicated(chain.raw.125[-(1:burnIn),ii]))) * length(Names)   
  print(acceptance)
  mtext(side = 3 ,acceptance, line = 0.2, cex = 0.5)
}
mtext(side = 1, 'Iterations', outer = TRUE)
mtext(side = 2, 'Rates', outer = TRUE)

dev.off()
setwd(cd)

# Correlation plots  ---------------------------------------
setwd(figure.path)

pdf('CorrelationPlots_Data_Dose625.pdf',width=7.5,height=7.5)
par(mfrow = c(11,11),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

betterPairs(data.frame(chain.625[,]))


dev.off()
setwd(cd)
