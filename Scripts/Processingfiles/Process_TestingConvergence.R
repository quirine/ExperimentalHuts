rm(list=ls())
# load libraries ----------------------------------------------------------
library(coda)
library(IDPmisc)
library(abind)

# set wd ------------------------------------------------------------------
setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# load source files -------------------------------------------------------
source('SR.Exit.FiguresScripts.R')

# set data directory ------------------------------------------------------
data.path <- "../Output/Final/output.ar.ltfudata.relsd.gamma.beta.fixedr_2016-12-13/"
# data.path <- "../Output/Final/output.ar.ltfudata.relsd.gamma.dd_2016-12-03/"

# Set directory to save figures -------------------------------------------
figure.path <- '../Figures/Final/'
cd <- getwd()

# specify file names ------------------------------------------------------
setwd(data.path)

list.of.files = dir(pattern = '*.csv')

files.0 = list.of.files[1:5]
files.l = list.of.files[6:10]
files.h = list.of.files[11:15]


setwd(cd)

# set parameters -----------------------------------------------------
Names=c('movement 2 away','movement 1 away','movement SR', 
        'prop away from SR', 
        'exit 2 away', 'exit 1 away', 'exit SR hut',
        'KD 2 away', 'KD 1 away', 'KD SR hut', 'loss-to-follow-up' )

Gelman.Names = c(expression('Movement (q' [2]*')'),expression('Movement (q' [1]*')'),
                 expression('Movement (q' [T]*')'),expression('Prop away from SR (p' [1]*')'),
                 expression('Exit (x' [2]*')'),expression('Exit (x' [1]*')'),
                 expression('Exit (x' [T]*')'),
                 expression('KD (k' [2]*')'),expression('KD (k' [1]*')'),
                 expression('KD (k' [T]*')'),'ltfu (u)')

Cor.Names = c('q2','q1','qT','p1','r2','r1','rT','k2','k1','kT','u')

fixed.r = FALSE

# load and prepare data ---------------------------------------------------
setwd(data.path)

SIZE.0 <- SIZE.l <-  SIZE.h <- numeric()
for (ii in 1:length(files.0)) {
  chain = read.csv(files.0[ii],header=TRUE,colClasses = 'numeric')
  if (dim(chain)[2] > 12) { temp = ((dim(chain)[2] - 1) / 2); chain = chain[-1,1:temp]  }
  else { chain = chain[-1,] }
  SIZE.0 = c(SIZE.0,dim(chain)[1])
  chain = read.csv(files.l[ii],header=TRUE,colClasses = 'numeric')
  if (dim(chain)[2] > 12) { temp = ((dim(chain)[2] - 1) / 2); chain = chain[-1,1:temp]  }
  else {  chain = chain[-1,] }
  SIZE.l = c(SIZE.l,dim(chain)[1])
  chain = read.csv(files.h[ii],header=TRUE,colClasses = 'numeric')
  if (dim(chain)[2] > 12) { temp = ((dim(chain)[2] - 1) / 2); chain = chain[-1,1:temp]  }
  else {  chain = chain[-1,] }
  SIZE.h = c(SIZE.h,dim(chain)[1])
}
size.0 = min(SIZE.0)
size.l = min(SIZE.l)
size.h = min(SIZE.h)

size.0 <- size.l <- size.h <- 9e4 #25000

for (ii in 1:length(files.0)){
  chain = read.csv(files.0[ii],header=TRUE,colClasses = 'numeric')
  if (dim(chain)[2] > 12) { temp = ((dim(chain)[2] - 1) / 2); 
                            LL = chain[-1,temp + 1]; props = chain[-1,(temp + 2):dim(chain)[2]];
                            chain = chain[-1,1:temp]
  }
  else {  chain = chain[-1,] }
  chain = mcmc(chain[1:size.0,-5])
  name <- paste("chain.raw.0.", ii, sep = "")
  assign(name, chain)
  name <- paste("LL.0.", ii, sep = "")
  assign(name, LL)
  name <- paste("props.0.", ii, sep = "")
  assign(name, props)
}
chains = apropos("^chain.raw.0.")
combinedchains.raw.0 = mcmc.list(chain.raw.0.1, chain.raw.0.2, chain.raw.0.3, chain.raw.0.4,chain.raw.0.5)
GEL.0 = gelman.diag(combinedchains.raw.0,multivariate = FALSE)
LL.0 = rbind(LL.0.1[1:size.0], LL.0.2[1:size.0], LL.0.3[1:size.0], LL.0.4[1:size.0], LL.0.5[1:size.0])
Props.0 = abind(props.0.1[1:size.0,], props.0.2[1:size.0,], props.0.3[1:size.0,], props.0.4[1:size.0,], props.0.5[1:size.0,], along = 3 )

for (ii in 1:length(files.l)){
  chain = read.csv(files.l[ii],header=TRUE,colClasses = 'numeric')
  if (dim(chain)[2] > 12) { 
    temp = ((dim(chain)[2] - 1) / 2); 
    LL = chain[-1,temp + 1]; props = chain[-1,(temp + 2):dim(chain)[2]];
    chain = chain[-1,1:temp]
  }
  else {  chain = chain[-1,] }
  chain = mcmc(chain[1:size.l,-5])
  name <- paste("chain.raw.l.", ii, sep = "")
  assign(name, chain)
  name <- paste("LL.l.", ii, sep = "")
  assign(name, LL)
  name <- paste("props.l.", ii, sep = "")
  assign(name, props)
}
chains = apropos("^chain.raw.l.")

combinedchains.raw.l = mcmc.list(chain.raw.l.1, chain.raw.l.2, chain.raw.l.3, chain.raw.l.4,chain.raw.l.5)
GEL.l = gelman.diag(combinedchains.raw.l,multivariate = FALSE)
LL.l = c(LL.l.1[1:size.l], LL.l.2[1:size.l], LL.l.3[1:size.l], LL.l.4[1:size.l], LL.l.5[1:size.l])
Props.l = abind(props.l.1[1:size.l,], props.l.2[1:size.l,], props.l.3[1:size.l,], props.l.4[1:size.l,], props.l.5[1:size.l,], along = 3 )

for (ii in 1:length(files.h)){
  chain = read.csv(files.h[ii],header=TRUE,colClasses = 'numeric')
  if (dim(chain)[2] > 12) { 
    temp = ((dim(chain)[2] - 1) / 2); 
    LL = chain[-1,temp + 1]; props = chain[-1,(temp + 2):dim(chain)[2]];
    chain = chain[-1,1:temp]
  }
  else {  chain = chain[-1,] }
  chain = mcmc(chain[1:size.h,-5])
  name <- paste("chain.raw.h.", ii, sep = "")
  assign(name, chain)
  name <- paste("LL.h.", ii, sep = "")
  assign(name, LL)
  name <- paste("props.h.", ii, sep = "")
  assign(name, props)
}
chains = apropos("^chain.raw.h.")
combinedchains.raw.h = mcmc.list(chain.raw.h.1, chain.raw.h.2, chain.raw.h.3, chain.raw.h.4,chain.raw.h.5)
GEL.h = gelman.diag(combinedchains.raw.h,multivariate = FALSE)
LL.h = c(LL.h.1[1:size.h], LL.h.2[1:size.h], LL.h.3[1:size.h], LL.h.4[1:size.h], LL.h.5[1:size.h])
Props.h = abind(props.h.1[1:size.h,], props.h.2[1:size.h,], props.h.3[1:size.h,], props.h.4[1:size.h,], props.h.5[1:size.h,], along = 3 )

setwd(cd)

# create gelman.plots to assess burnin period (control)------------------------------------------------------------
setwd(figure.path)

Iis = c(seq(1,4),seq(8,11))

pdf('FigureS7.pdf',width=4.75,height=3.00)
par(mfrow = c(2,4),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  #mar=c(1,3,1,1.9)

for (ii in Iis) { 
  gelman.plot(combinedchains.raw.0[,ii],auto.layout = FALSE,ask = FALSE, cex = 0.8, ylim = c(0,2))
  mtext(side = 3, Gelman.Names[ii], line = 1, cex = 0.6)
  if (ii == 4) legend('topright',c('median','97.5% CI'), col = c('black','red'), lty=c(1,2), cex = 0.5, bty = 'n')
}
mtext(side = 1, 'Iterations', outer = TRUE, cex = 0.8)
mtext(side = 2, 'Shrink factor', outer = TRUE, cex = 0.8)

dev.off()
setwd(cd)

# create gelman.plots to assess burnin period (low dosage)------------------------------------------------------------
setwd(figure.path)

Iis = c(seq(1,4),seq(8,11))

pdf('FigureS8.pdf',width=4.75,height=3.00)
par(mfrow = c(2,4),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (ii in Iis) { 
  gelman.plot(combinedchains.raw.l[,ii],auto.layout = FALSE,ask = FALSE, cex = 0.8, ylim = c(0,2))
  mtext(side = 3, Gelman.Names[ii], line = 1, cex = 0.6)
  if (ii == 4) legend('topright',c('median','97.5% CI'), col = c('black','red'), lty=c(1,2), cex = 0.5, bty = 'n')
}
mtext(side = 1, 'Iterations', outer = TRUE, cex = 0.8)
mtext(side = 2, 'Shrink factor', outer = TRUE, cex = 0.8)

dev.off()
setwd(cd)

# create gelman.plots to assess burnin period (high dosage)------------------------------------------------------------
setwd(figure.path)

Iis = c(seq(1,4),seq(8,11))

pdf('FigureS9.pdf',width=4.75,height=3.00)
par(mfrow = c(2,4),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (ii in Iis) { 
  gelman.plot(combinedchains.raw.h[,ii],auto.layout = FALSE,ask = FALSE, cex = 0.8, ylim = c(0,2))
  mtext(side = 3, Gelman.Names[ii], line = 1, cex = 0.6)
  if (ii == 4) legend('topright',c('median','97.5% CI'), col = c('black','red'), lty=c(1,2), cex = 0.5, bty = 'n')
}
mtext(side = 1, 'Iterations', outer = TRUE, cex = 0.8)
mtext(side = 2, 'Shrink factor', outer = TRUE, cex = 0.8)

dev.off()
setwd(cd)

# Trace plots ---------------------------------------
setwd(figure.path)

pdf('FigureS10.pdf',width=4.75,height=3.00)
par(mfrow = c(2,4),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (iii in 1:length(Iis)){ 
  ii = Iis[iii]
  traceplot(combinedchains.raw.0[,ii],cex = 0.8)
  mtext(side = 3, Gelman.Names[ii], line = 1, cex = 0.6)
}
mtext(side = 1, 'Iterations', outer = TRUE, cex = 0.8)
mtext(side = 2, 'Rates', outer = TRUE, cex = 0.8)

dev.off()
setwd(cd)

# Trace plots ---------------------------------------
setwd(figure.path)

pdf('FigureS11.pdf',width=4.75,height=3.00)
par(mfrow = c(2,4),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (iii in 1:length(Iis)){ 
  ii = Iis[iii]
  traceplot(combinedchains.raw.l[,ii], cex = 0.8)
  mtext(side = 3, Gelman.Names[ii], line = 1, cex = 0.6)
}
mtext(side = 1, 'Iterations', outer = TRUE, cex = 0.8)
mtext(side = 2, 'Rates', outer = TRUE, cex = 0.8)

dev.off()
setwd(cd)

# Trace plots ---------------------------------------
setwd(figure.path)

pdf('FigureS12.pdf',width=4.75,height=3.00)
par(mfrow = c(2,4),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (iii in 1:length(Iis)){ 
  ii = Iis[iii]
  traceplot(combinedchains.raw.h[,ii], cex = 0.8)
  mtext(side = 3, Gelman.Names[ii], line = 1, cex = 0.6)
}
mtext(side = 1, 'Iterations', outer = TRUE, cex = 0.8)
mtext(side = 2, 'Rates', outer = TRUE, cex = 0.8)

dev.off()
setwd(cd)


# Set burnin and derive summary statistics --------------------------------------------
burnIn = 10000

chain.0.1 = mcmc(chain.raw.0.1[burnIn:size.0,]) 
chain.0.2 = mcmc(chain.raw.0.2[burnIn:size.0,])
chain.0.3 = mcmc(chain.raw.0.3[burnIn:size.0,])
chain.0.4 = mcmc(chain.raw.0.4[burnIn:size.0,])
chain.0.5 = mcmc(chain.raw.0.5[burnIn:size.0,])
combinedchains.0 = mcmc.list(chain.0.1, chain.0.2, chain.0.3, chain.0.4,chain.0.5)

chain.l.1 = mcmc(chain.raw.l.1[burnIn:size.l,])
chain.l.2 = mcmc(chain.raw.l.2[burnIn:size.l,])
chain.l.3 = mcmc(chain.raw.l.3[burnIn:size.l,])
chain.l.4 = mcmc(chain.raw.l.4[burnIn:size.l,])
chain.l.5 = mcmc(chain.raw.l.5[burnIn:size.l,])
combinedchains.l = mcmc.list(chain.l.1, chain.l.2, chain.l.3, chain.l.4,chain.l.5)

chain.h.1 = mcmc(chain.raw.h.1[burnIn:size.h,])
chain.h.2 = mcmc(chain.raw.h.2[burnIn:size.h,])
chain.h.3 = mcmc(chain.raw.h.3[burnIn:size.h,])
chain.h.4 = mcmc(chain.raw.h.4[burnIn:size.h,])
chain.h.5 = mcmc(chain.raw.h.5[burnIn:size.h,])
combinedchains.h = mcmc.list(chain.h.1, chain.h.2, chain.h.3, chain.h.4,chain.h.5)

chain.0 = combinedchains.0  
chain.625 = combinedchains.l 
chain.125 = combinedchains.h 

output.summary.0 = summary(chain.0)
output.summary.625 = summary(chain.625)
output.summary.125 = summary(chain.125)

GEL.0.bi = gelman.diag(chain.0,multivariate = FALSE)
GEL.l.bi = gelman.diag(chain.625,multivariate = FALSE)
GEL.h.bi = gelman.diag(chain.125,multivariate = FALSE)

# get r’s for initialization ----------------------------------------------
Rs = output.summary.0$quantiles[5:7,]
Res.times = 1/(output.summary.0$quantiles[1:3,])

# Get acceptance rates ----------------------------------------------------

Acceptance.0 <- Acceptance.l <- Acceptance.h <- matrix(data=NA,nrow=output.summary.0$nchain,ncol=dim(chain)[2])
for (ii in 1:output.summary.0$nchain) { 
  for (jj in 1:dim(chain)[2]){
    name <- paste("chain.0.", ii, sep = "")
    assign(name, chain)
    acceptance = (1-mean(duplicated(chain[-(1),jj]))) * length(Names)
    print(acceptance)
    Acceptance.0[ii,jj] = acceptance
  }
}

for (ii in 1:output.summary.0$nchain) { 
  for (jj in 1:dim(chain)[2]){
    name <- paste("chain.l.", ii, sep = "")
    assign(name, chain)
    acceptance = (1-mean(duplicated(chain[-(1),jj]))) * length(Names)
    print(acceptance)
    Acceptance.l[ii,jj] = acceptance
  }
}

for (ii in 1:output.summary.0$nchain) { 
  for (jj in 1:dim(chain)[2]){
    name <- paste("chain.h.", ii, sep = "")
    assign(name, chain)
    acceptance = (1-mean(duplicated(chain[-(1),jj]))) * length(Names)
    print(acceptance)
    Acceptance.h[ii,jj] = acceptance
  }
}


# Correlations ------------------------------------------------------------

crosscorr(chain.0)
crosscorr(chain.625)
crosscorr(chain.125)
min(crosscorr(chain.0),na.rm=TRUE)
min(crosscorr(chain.625),na.rm=TRUE)
min(crosscorr(chain.125),na.rm=TRUE)


# Density plots with real values - movement rates ---------------------------------------

setwd(figure.path)

# Names=rownames(output.summary.0$statistics)

pdf('DensityPlots_movement.pdf',width=4.75,height=4.75)
par(mfrow = c(3,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (ii in 1:3){ 
  densplot(chain.0[,ii],show.obs = TRUE, xlim=c(0,.1))
  abline(v=output.summary.0$quantiles[ii,3],col='red')
  print(output.summary.0$quantiles[ii,3])
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 1) { mtext( side = 2, 'control', line = 2, cex = 0.5) }
  
}
for (ii in 1:3){ 
  
  densplot(chain.625[,ii],show.obs = TRUE, xlim=c(0,.1))
  abline(v=output.summary.625$quantiles[ii,3],col='red')
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 1) { mtext( side = 2, '0.0625 FAR', line = 2, cex = 0.5) }
}
for (ii in 1:3){ 
  densplot(chain.125[,ii],show.obs = TRUE, xlim=c(0,.1))
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
ii = 11

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


# Density plots with real values - exit proportions---------------------------------------
setwd(figure.path)

# Names=rownames(output.summary.0$statistics)

pdf('DensityPlots_exit_prop.pdf',width=4.75,height=4.75)
par(mfrow = c(3,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

for (ii in 5:7){ 
  densplot(chain.0[,ii],show.obs = TRUE, xlim=c(0,0.15))
  abline(v=output.summary.0$quantiles[ii,3],col='red')
  print(output.summary.0$quantiles[ii,3])
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 5) { mtext( side = 2, 'control', line = 2, cex = 0.5) }
  
}
for (ii in 5:7){ 
  
  densplot(chain.625[,ii],show.obs = TRUE, xlim=c(0,0.15))
  abline(v=output.summary.625$quantiles[ii,3],col='red')
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 5) { mtext( side = 2, '0.0625 FAR', line = 2, cex = 0.5) }
}
for (ii in 5:7){ 
  densplot(chain.125[,ii],show.obs = TRUE, xlim=c(0,0.15))
  abline(v=output.summary.125$quantiles[ii,3],col='red')
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 5) { mtext( side = 2, '0.125 FAR', line = 2, cex = 0.5) }
}
mtext(side = 1, 'Proportion', outer = TRUE)
mtext(side = 2, 'Density', outer = TRUE, line = 2)

dev.off()
setwd(cd)

# Density plots with real values - exit rates IN PROGRESS!!!---------------------------------------
setwd(figure.path)

# Names=rownames(output.summary.0$statistics)

pdf('DensityPlots_exit_rate.pdf',width=4.75,height=4.75)
par(mfrow = c(3,3),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  #mar=c(1,3,1,1.9)

for (ii in 5:7){ 
  densplot(chain.0[,ii]*chain.0[,ii-4],show.obs = TRUE, xlim=c(0,0.15))
  print(output.summary.0$quantiles[ii,3])
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 5) { mtext( side = 2, 'control', line = 2, cex = 0.5) }
  
}
for (ii in 5:7){ 
  
  densplot(chain.625[,ii]*chain.625[,ii-4],show.obs = TRUE, xlim=c(0,0.15))
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 5) { mtext( side = 2, '0.0625 FAR', line = 2, cex = 0.5) }
}
for (ii in 5:7){ 
  densplot(chain.125[,ii]*chain.125[,ii-4],show.obs = TRUE, xlim=c(0,0.15))
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 5) { mtext( side = 2, '0.125 FAR', line = 2, cex = 0.5) }
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

for (ii in 8:10){ 
  densplot(chain.0[,ii],show.obs = TRUE, xlim=c(0,0.001))
  abline(v=output.summary.0$quantiles[ii,3],col='red')
  print(output.summary.0$quantiles[ii,3])
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 8) { mtext( side = 2, 'control', line = 2, cex = 0.5) }
  
}
for (ii in 8:10){ 
  
  densplot(chain.625[,ii],show.obs = TRUE, xlim=c(0,0.001))
  abline(v=output.summary.625$quantiles[ii,3],col='red')
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 8) { mtext( side = 2, '0.0625 FAR', line = 2, cex = 0.5) }
}
for (ii in 8:10){ 
  densplot(chain.125[,ii],show.obs = TRUE, xlim=c(0,0.001))
  abline(v=output.summary.125$quantiles[ii,3],col='red')
  mtext(side = 3, Names[ii], line = 1, cex =0.5)
  if (ii == 8) { mtext( side = 2, '0.125 FAR', line = 2, cex = 0.5) }
}
mtext(side = 1, 'Rates', outer = TRUE)
mtext(side = 2, 'Density', outer = TRUE, line = 2)

dev.off()
setwd(cd)


# Create output matrix for correlation plots (just for the all parameter case. Otherwise done in Process_Results)------------------------------
someData = rep(NaN,((size-burnIn+1)*5)*11)  # num.chains is hard coded here!!!!
Output.all = array(someData, c(((size.0-burnIn+1)*5),11))
dimnames(Output.all)[[2]] <- Cor.Names

for (ii in 1:11) {
  temp.mult = numeric()
  for (kk in 1:5){
    temp = eval(parse(text=paste('chain.0.',kk,sep='')))
    temp.mult = c( temp.mult,temp[,ii] ) 
  }
  Output.all[,ii] = temp.mult
}

# # Correlation plots  ---------------------------------------
setwd(figure.path)

pdf('FigureS1.pdf',width=10,height=10)
par(mfrow = c(11,11),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

temp = Output.all
betterPairs(data.frame(temp))


dev.off()
setwd(cd)

