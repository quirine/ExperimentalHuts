rm(list=ls())

# load libraries ----------------------------------------------------------
library(coda)
library(IDPmisc)
library(abind)
library(ggplot2)

# set wd ------------------------------------------------------------------
setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# load source files -------------------------------------------------------
source('SR.Exit.FiguresScripts.R')

# set data directory ------------------------------------------------------
data.path <- "../Output/Final/output.ar.ltfudata.relsd.gamma.beta.fixedr_2016-12-13/"
# data.path <- "../Output/Final/output.ar.ltfudata.relsd.gamma.beta.fixedr.low_2016-12-13/"
# data.path <- "../Output/Final/output.ar.ltfudata.relsd.gamma.beta.fixedr.high_2016-12-13/"

# Set directory to save figures -------------------------------------------
figure.path <- '~/Box Sync/Research Backup/Research/Iquitos/Analysis/Movement/MarkovChains/Figures/Final/'
cd <- getwd()

# specify file names ------------------------------------------------------
setwd(data.path)

list.of.files = dir(pattern = '*.csv')

files.0 = list.of.files[1:5]
files.l = list.of.files[6:10]
files.h = list.of.files[11:15]

setwd(cd)

# set parameters -----------------------------------------------------
Names=c('exit rate 2 away','exit rate 1 away','exit rate SR', 
        'prop away from SR', 
        'exit 2 away', 'exit 1 away', 'exit SR hut',
        'KD 2 away', 'KD 1 away', 'KD SR hut', 'loss-to-follow-up' )

Cor.Names = c('x2','x1','xT','p1','k2','k1','kT','u')
num.chains = 5
Dosages = seq(0,2,0.001)
colorlist= colorRampPalette(c("lightcyan3", "darkorange","deeppink"))(length(Dosages))
colorlist = c(rgb(0,0,0), colorlist)

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
# size = min (size.0, size.l, size.h)
size =9e4


for (ii in 1:length(files.0)){
  chain = read.csv(files.0[ii],header=TRUE,colClasses = 'numeric')
  if (dim(chain)[2] > 12) { temp = ((dim(chain)[2] - 1) / 2); 
                            LL = chain[-1,temp + 1]; props = chain[-1,(temp + 2):dim(chain)[2]];
                            chain = chain[-1,1:temp]
  }
  else {  chain = chain[-1,] }
  chain = mcmc(chain[1:size,-5])
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
LL.0 = rbind(LL.0.1[1:size], LL.0.2[1:size], LL.0.3[1:size], LL.0.4[1:size], LL.0.5[1:size])
Props.0 = abind(props.0.1[1:size,], props.0.2[1:size,], props.0.3[1:size,], props.0.4[1:size,], props.0.5[1:size,], along = 3 )

for (ii in 1:length(files.l)){
  chain = read.csv(files.l[ii],header=TRUE,colClasses = 'numeric')
  if (dim(chain)[2] > 12) { 
    temp = ((dim(chain)[2] - 1) / 2); 
    LL = chain[-1,temp + 1]; props = chain[-1,(temp + 2):dim(chain)[2]];
    chain = chain[-1,1:temp]
  }
  else {  chain = chain[-1,] }
  chain = mcmc(chain[1:size,-5])
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
LL.l = c(LL.l.1[1:size], LL.l.2[1:size], LL.l.3[1:size], LL.l.4[1:size], LL.l.5[1:size])
Props.l = abind(props.l.1[1:size,], props.l.2[1:size,], props.l.3[1:size,], props.l.4[1:size,], props.l.5[1:size,], along = 3 )

for (ii in 1:length(files.h)){
  chain = read.csv(files.h[ii],header=TRUE,colClasses = 'numeric')
  if (dim(chain)[2] > 12) { 
    temp = ((dim(chain)[2] - 1) / 2); 
    LL = chain[-1,temp + 1]; props = chain[-1,(temp + 2):dim(chain)[2]];
    chain = chain[-1,1:temp]
  }
  else {  chain = chain[-1,] }
  chain = mcmc(chain[1:size,-5])
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
LL.h = c(LL.h.1[1:size], LL.h.2[1:size], LL.h.3[1:size], LL.h.4[1:size], LL.h.5[1:size])
Props.h = abind(props.h.1[1:size,], props.h.2[1:size,], props.h.3[1:size,], props.h.4[1:size,], props.h.5[1:size,], along = 3 )

setwd(cd)

# Set burnin and derive summary statistics --------------------------------------------
burnIn = 10000

chain.0.1 = mcmc(chain.raw.0.1[burnIn:size,])
chain.0.2 = mcmc(chain.raw.0.2[burnIn:size,])
chain.0.3 = mcmc(chain.raw.0.3[burnIn:size,])
chain.0.4 = mcmc(chain.raw.0.4[burnIn:size,])
chain.0.5 = mcmc(chain.raw.0.5[burnIn:size,])
combinedchains.0 = mcmc.list(chain.0.1, chain.0.2, chain.0.3, chain.0.4,chain.0.5)

chain.l.1 = mcmc(chain.raw.l.1[burnIn:size,])
chain.l.2 = mcmc(chain.raw.l.2[burnIn:size,])
chain.l.3 = mcmc(chain.raw.l.3[burnIn:size,])
chain.l.4 = mcmc(chain.raw.l.4[burnIn:size,])
chain.l.5 = mcmc(chain.raw.l.5[burnIn:size,])
combinedchains.l = mcmc.list(chain.l.1, chain.l.2, chain.l.3, chain.l.4,chain.l.5)

chain.h.1 = mcmc(chain.raw.h.1[burnIn:size,])
chain.h.2 = mcmc(chain.raw.h.2[burnIn:size,])
chain.h.3 = mcmc(chain.raw.h.3[burnIn:size,])
chain.h.4 = mcmc(chain.raw.h.4[burnIn:size,])
chain.h.5 = mcmc(chain.raw.h.5[burnIn:size,])
combinedchains.h = mcmc.list(chain.h.1, chain.h.2, chain.h.3, chain.h.4,chain.h.5)

chain.0 = combinedchains.0  
chain.625 = combinedchains.l 
chain.125 = combinedchains.h 

output.summary.0 = summary(chain.0)
output.summary.625 = summary(chain.625)
output.summary.125 = summary(chain.125)

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

# Main Figure - hut paper: MCMC-results by dose ---------------------------------------
setwd(figure.path)

someData = rep(NaN,((size-burnIn+1)*num.chains)*3*8)
Output.all = array(someData, c(((size-burnIn+1)*num.chains),3,8))
dimnames(Output.all)[[3]] <- Cor.Names

transparancy = 0.5
XLIMS = rbind(c(5e-4,4e-3),c(5e-4,4e-3),c(5e-4,4e-3),
              c(0.4,0.75), c(NA,NA), c(NA,NA),c(NA,NA),
              c(0,2e-4), c(0,2e-4),c(0,7e-4),
              c(5e-4,1e-3) )

YLIMS = rbind(c(0,9e-3),c(0,9e-3),c(0,11e-3),
              c(0,7e-3), c(NA,NA), c(NA,NA),c(NA,NA),
              c(0,1.6e-2), c(0,1.6e-2),c(0,1.6e-2),
              c(0,7e-3) )
# Names=rownames(output.summary.0$statistics)

pdf('Figure5.pdf',width=7.5,height=4.75)
par(mfrow = c(2,4),mai=c(0.3,0.3,0.3,0.15) ,oma=c(3,4,1,1))  

for (ii in 1:3){ 
  temp.0.mult <- temp.l.mult <- temp.h.mult <- numeric()
  for (kk in 1:num.chains){
    temp.0.mult = c( temp.0.mult,chain.0[[kk]][,ii] * chain.0[[kk]][,ii + 4] )
    temp.625.mult = c( temp.l.mult,chain.625[[kk]][,ii] * chain.625[[kk]][,ii + 4] )
    temp.125.mult = c( temp.h.mult,chain.125[[kk]][,ii] * chain.125[[kk]][,ii + 4] )
  } 
  len = min (length(temp.0.mult),length(temp.625.mult),length(temp.125.mult))
  Temp = rbind(temp.0.mult, temp.625.mult, temp.125.mult)
  if (ii == 1) Output.all[,,1] <- Exit_2 <- Temp; if (ii == 2) Output.all[,,2] <- Exit_1 <- Temp; if (ii == 3) Output.all[,,3] <- Exit_0 <- Temp; 
  p1 <- density(temp.0.mult);     p1$y <- p1$y*diff(p1$x)[1]
  p2 <- density(temp.625.mult);   p2$y <- p2$y*diff(p2$x)[1]
  p3 <- density(temp.125.mult);   p3$y <- p3$y*diff(p3$x)[1]
  plot(p1,col='white',main='',xlim=XLIMS[ii,],ylim=YLIMS[ii,],xaxs='i',yaxs='i',las=1)
  polygon( p1, col=alpha('black',transparancy),main='',border = 'black')  
  polygon( p2, col=alpha('darkorange',transparancy),border = 'darkorange')  
  polygon( p3, col=alpha('deeppink',transparancy),border = 'deeppink')  
  if (ii == 1) { mtext( bquote(bold("A")), side = 3, line = .2 ,at =6e-4,cex=1); mtext(text=expression('Exit rate (' ~x[2]* ')'),side=1,line=2.25,cex=0.8) } 
  if (ii == 2) { mtext( bquote(bold("B")), side = 3, line = .2 ,at =6e-4,cex=1); mtext(text=expression('Exit rate (' ~x[1]* ')'),side=1,line=2.25,cex=0.8) }  
  if (ii == 3) { mtext( bquote(bold("C")), side = 3, line = .2 ,at =6e-4,cex=1); mtext(text=expression('Exit rate (' ~x[T]* ')'),side=1,line=2.25,cex=0.8) }  
  
}
for (ii in c(4,seq(8,11))){ 
  
  temp.0.mult <- temp.l.mult <- temp.h.mult <- numeric()
  for (kk in 1:num.chains){
    temp.0.mult = c( temp.0.mult,chain.0[[kk]][,ii])
    temp.625.mult = c( temp.l.mult,chain.625[[kk]][,ii]  )
    temp.125.mult = c( temp.h.mult,chain.125[[kk]][,ii]  )
  } 
  Temp = rbind(temp.0.mult, temp.625.mult, temp.125.mult)
  if (ii == 4) Output.all[,,4] <- PBorD <- Temp; if (ii == 8) Output.all[,,5] <- Kd_2 <- Temp; if (ii == 9) Output.all[,,6] <- Kd_1 <- Temp; if (ii == 10) Output.all[,,7] <- Kd_0 <- Temp; if (ii == 11) Output.all[,,8] <- Us <- Temp; 
  p1 <- density(temp.0.mult);     p1$y <- p1$y*diff(p1$x)[1]
  p2 <- density(temp.625.mult);   p2$y <- p2$y*diff(p2$x)[1]
  p3 <- density(temp.125.mult);   p3$y <- p3$y*diff(p3$x)[1]
  plot(p1,col='white',main='',xlim=XLIMS[ii,],ylim=YLIMS[ii,],xaxs='i',yaxs='i',las=1)
  polygon( p1, col=alpha('black',transparancy),main='',border = 'black')  
  polygon( p2, col=alpha('darkorange',transparancy),border = 'darkorange')  
  polygon( p3, col=alpha('deeppink',transparancy),border = 'deeppink')  
  if (ii == 4) { mtext( bquote(bold("D")), side = 3, line = .2 ,at =0.42,cex=1); mtext(text=expression('Prop away from SR (' ~p[1]* ')'),side=1,line=2.25,cex=0.8) }  
  if (ii == 8) { mtext( bquote(bold("E")), side = 3, line = .2 ,at =6e-6,cex=1); mtext(text=expression('Knock down rate (' ~k[2]* ')'),side=1,line=2.25,cex=0.8) }   
  if (ii == 9) { mtext( bquote(bold("F")), side = 3, line = .2 ,at =6e-6,cex=1); mtext(text=expression('Knock down rate (' ~k[1]* ')'),side=1,line=2.25,cex=0.8) }   
  if (ii == 10){ mtext( bquote(bold("G")), side = 3, line = .2 ,at =1.5e-5,cex=1); mtext(text=expression('Knock down rate (' ~k[T]* ')'),side=1,line=2.25,cex=0.8) }   
  if (ii == 11){ mtext( bquote(bold("H")), side = 3, line = .2 ,at =5.2e-4,cex=1); mtext(text=expression('Loss to follow up rate (u)'),side=1,line=2.25,cex=0.8) }   
  
  if ( ii == 4) abline(v=0.5,col = 'gray', lwd = 2,lty = 3 )  
}

mtext(side = 2, 'Density', outer = TRUE, line = 1.3)
legend('topright',c('control','low','high'), fill=c(alpha('black',transparancy),alpha('darkorange',transparancy),alpha('deeppink',transparancy)),box.lty=0)
# 
dev.off()
setwd(cd)

# # Correlation plots  ---------------------------------------
setwd(figure.path)

pdf('FigureS2.pdf',width=7.5,height=7.5)
par(mfrow = c(11,11),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

temp = Output.all[,1,]
betterPairs(data.frame(temp))


dev.off()
setwd(cd)

# # Correlation plots  ---------------------------------------
setwd(figure.path)

pdf('FigureS3.pdf',width=7.5,height=7.5)
par(mfrow = c(11,11),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

temp = Output.all[,2,]
betterPairs(data.frame(temp))


dev.off()
setwd(cd)

# # Correlation plots  ---------------------------------------
setwd(figure.path)

pdf('FigureS4.pdf',width=7.5,height=7.5)
par(mfrow = c(11,11),mai=c(0.3,0.3,0.3,0.1) ,oma=c(3,4,1,1))  

temp = Output.all[,3,]
betterPairs(data.frame(temp))


dev.off()
setwd(cd)

# Main Figure - FOI paper: exit rates and repellency ---------------------------------------
setwd(figure.path)

transparancy = 0.5
XLIMS = rbind(c(5e-4,4e-3),c(5e-4,4e-3),c(5e-4,4e-3),
              c(0.4,0.75), c(NA,NA), c(NA,NA),c(NA,NA),
              c(0,2e-4), c(0,2e-4),c(0,7e-4),
              c(5e-4,1e-3) )
# Names=rownames(output.summary.0$statistics)

pdf('MainFig_Posteriors_MovementEstimates.pdf',width=7.5,height=4.75)
par(mfrow = c(1,2),mai=c(0.3,0.6,0.6,0.1) ,oma=c(3,4,1,1))  

temp.0.mult <- temp.l.mult <- temp.h.mult <- numeric()
for (kk in 1:num.chains){
  temp.0.mult = c( temp.0.mult,chain.0[[kk]][,4]);           
  temp.625.mult = c( temp.l.mult,chain.625[[kk]][,4]  )
  temp.125.mult = c( temp.h.mult,chain.125[[kk]][,4]  )
} 
Repellency = rbind( ( (temp.0.mult[runif(n=1e5,min=1,max=size)] - (temp.0.mult[runif(n=1e5,min=1,max=size)]) ) / (1-(temp.0.mult[runif(n=1e5,min=1,max=size)]) ) ),
                    ( (temp.625.mult[runif(n=1e5,min=1,max=size)] - (temp.0.mult[runif(n=1e5,min=1,max=size)]) ) / (1-(temp.0.mult[runif(n=1e5,min=1,max=size)]) ) ) ,
                    ( (temp.125.mult[runif(n=1e5,min=1,max=size)] - (temp.0.mult[runif(n=1e5,min=1,max=size)]) ) / (1-(temp.0.mult[runif(n=1e5,min=1,max=size)]) ) ) ) 
p1 <- density( ( (temp.0.mult[runif(n=1e5,min=1,max=size)] - (temp.0.mult[runif(n=1e5,min=1,max=size)]) ) / (1-(temp.0.mult[runif(n=1e5,min=1,max=size)]) ) ), na.rm = TRUE);   p1$y <- p1$y*diff(p1$x)[1]
p2 <- density( ( (temp.625.mult[runif(n=1e5,min=1,max=size)] - (temp.0.mult[runif(n=1e5,min=1,max=size)]) ) / (1-(temp.0.mult[runif(n=1e5,min=1,max=size)]) ) ) ,na.rm=TRUE);  p2$y <- p2$y*diff(p2$x)[1]  
p3 <- density( ( (temp.125.mult[runif(n=1e5,min=1,max=size)] - (temp.0.mult[runif(n=1e5,min=1,max=size)]) ) / (1-(temp.0.mult[runif(n=1e5,min=1,max=size)]) ) ) ,na.rm=TRUE);  p3$y <- p3$y*diff(p3$x)[1]
plot(p1,col='white',main='',xlim=c(-.4,.4), ylim=c(0,7.5e-3),xaxs='i',yaxs='i',las=1)
polygon( p1, col=alpha('black',transparancy),main='',border = 'black')  
polygon( p2, col=alpha('darkorange',transparancy),border = 'darkorange')  
polygon( p3, col=alpha('deeppink',transparancy),border = 'deeppink')  
mtext( bquote(bold("A")), side = 3, line = 0.1 ,at =-0.37,cex=1) 
mtext(text=expression('Repellency (' ~rho* ' )'),side=1,line=2.25)

temp.0.mult <- temp.l.mult <- temp.h.mult <- numeric()
for (kk in 1:num.chains){
  temp.0.mult = c( temp.0.mult,chain.0[[kk]][,3] * chain.0[[kk]][,3 + 4] )
  temp.625.mult = c( temp.l.mult,chain.625[[kk]][,3] * chain.625[[kk]][,3 + 4] )
  temp.125.mult = c( temp.h.mult,chain.125[[kk]][,3] * chain.125[[kk]][,3 + 4] )
} 
Expellency = rbind((temp.0.mult[runif(n=1e5,min=1,max=length(temp.0.mult))]/temp.0.mult[runif(n=1e5,min=1,max=length(temp.0.mult))]),
                   (temp.625.mult[runif(n=1e5,min=1,max=length(temp.625.mult))]/temp.0.mult[runif(n=1e5,min=1,max=length(temp.0.mult))]),
                   (temp.125.mult[runif(n=1e5,min=1,max=length(temp.125.mult))]/temp.0.mult[runif(n=1e5,min=1,max=length(temp.0.mult))]))
p1 <- density(temp.0.mult[runif(n=1e5,min=1,max=length(temp.0.mult))]/temp.0.mult[runif(n=1e5,min=1,max=length(temp.0.mult))]);   p1$y <- p1$y*diff(p1$x)[1]
p2 <- density(temp.625.mult[runif(n=1e5,min=1,max=length(temp.625.mult))]/temp.0.mult[runif(n=1e5,min=1,max=length(temp.0.mult))]); p2$y <- p2$y*diff(p2$x)[1]  
p3 <- density(temp.125.mult[runif(n=1e5,min=1,max=length(temp.125.mult))]/temp.0.mult[runif(n=1e5,min=1,max=length(temp.0.mult))]); p3$y <- p3$y*diff(p3$x)[1]
plot(p1,col='white',main='',xlim=c(0,2),ylim=c(0,1.2e-2),xaxs='i',yaxs='i',las=1)
polygon( p1, col=alpha('black',transparancy),main='',border = 'black')  
polygon( p2, col=alpha('darkorange',transparancy),border = 'darkorange')  
polygon( p3, col=alpha('deeppink',transparancy),border = 'deeppink')  
mtext( bquote(bold("B")), side = 3, line = 0.1 ,at =0.08,cex=1) 
mtext(text=expression('Expellency (' ~phi* ' )'),side=1,line=2.25)
legend('topright',c('control','low','high'), fill=c(alpha('black',transparancy),alpha('darkorange',transparancy),alpha('deeppink',transparancy)),box.lty=0)


mtext(side = 2, 'Density', outer = TRUE, line = 1)
# 
dev.off()
setwd(cd)


# Relative outcomes -------------------------------------------------------
Exit_2_rel = rbind( (Exit_2[1,runif(n=1e5,min=1,max=size)]/Exit_2[1,runif(n=1e5,min=1,max=size)]),
                    (Exit_2[2,runif(n=1e5,min=1,max=size)]/Exit_2[1,runif(n=1e5,min=1,max=size)]),
                    (Exit_2[3,runif(n=1e5,min=1,max=size)]/Exit_2[1,runif(n=1e5,min=1,max=size)]) )
Exit_1_rel = rbind( (Exit_1[1,runif(n=1e5,min=1,max=size)]/Exit_1[1,runif(n=1e5,min=1,max=size)]),
                    (Exit_1[2,runif(n=1e5,min=1,max=size)]/Exit_1[1,runif(n=1e5,min=1,max=size)]),
                    (Exit_1[3,runif(n=1e5,min=1,max=size)]/Exit_1[1,runif(n=1e5,min=1,max=size)]) )
Exit_0_rel = rbind( (Exit_0[1,runif(n=1e5,min=1,max=size)]/Exit_0[1,runif(n=1e5,min=1,max=size)]),
                    (Exit_0[2,runif(n=1e5,min=1,max=size)]/Exit_0[1,runif(n=1e5,min=1,max=size)]),
                    (Exit_0[3,runif(n=1e5,min=1,max=size)]/Exit_0[1,runif(n=1e5,min=1,max=size)]) )
PBorD_rel =  rbind( (PBorD[1,runif(n=1e5,min=1,max=size)]/PBorD[1,runif(n=1e5,min=1,max=size)]),
                    (PBorD[2,runif(n=1e5,min=1,max=size)]/PBorD[1,runif(n=1e5,min=1,max=size)]),
                    (PBorD[3,runif(n=1e5,min=1,max=size)]/PBorD[1,runif(n=1e5,min=1,max=size)]) )
Kd_2_rel =   rbind( (Kd_2[1,runif(n=1e5,min=1,max=size)]/Kd_2[1,runif(n=1e5,min=1,max=size)]),
                    (Kd_2[2,runif(n=1e5,min=1,max=size)]/Kd_2[1,runif(n=1e5,min=1,max=size)]),
                    (Kd_2[3,runif(n=1e5,min=1,max=size)]/Kd_2[1,runif(n=1e5,min=1,max=size)]) )
Kd_1_rel =   rbind( (Kd_1[1,runif(n=1e5,min=1,max=size)]/Kd_1[1,runif(n=1e5,min=1,max=size)]),
                    (Kd_1[2,runif(n=1e5,min=1,max=size)]/Kd_1[1,runif(n=1e5,min=1,max=size)]),
                    (Kd_1[3,runif(n=1e5,min=1,max=size)]/Kd_1[1,runif(n=1e5,min=1,max=size)]) )
Kd_0_rel =   rbind( (Kd_0[1,runif(n=1e5,min=1,max=size)]/Kd_0[1,runif(n=1e5,min=1,max=size)]),
                    (Kd_0[2,runif(n=1e5,min=1,max=size)]/Kd_0[1,runif(n=1e5,min=1,max=size)]),
                    (Kd_0[3,runif(n=1e5,min=1,max=size)]/Kd_0[1,runif(n=1e5,min=1,max=size)]) )
Us_rel =     rbind( (Us[1,runif(n=1e5,min=1,max=size)]/Us[1,runif(n=1e5,min=1,max=size)]),
                    (Us[2,runif(n=1e5,min=1,max=size)]/Us[1,runif(n=1e5,min=1,max=size)]),
                    (Us[3,runif(n=1e5,min=1,max=size)]/Us[1,runif(n=1e5,min=1,max=size)]) )

# Repellency specific outcomes ---------------------------------------
# proportion > 0.5
length(which(PBorD[1,]>0.5))/length(PBorD[1,])
length(which(PBorD[2,]>0.5))/length(PBorD[2,])
length(which(PBorD[3,]>0.5))/length(PBorD[3,])

# Odds
median(PBorD[1,]/(1-PBorD[1,]))
HDIofMCMC(mcmc(PBorD[1,]/(1-PBorD[1,])),credMass = 0.95)
median(PBorD[2,]/(1-PBorD[2,]))
HDIofMCMC(mcmc(PBorD[2,]/(1-PBorD[2,])),credMass = 0.95)
median(PBorD[3,]/(1-PBorD[3,]))
HDIofMCMC(mcmc(PBorD[3,]/(1-PBorD[3,])),credMass = 0.95)

# Loss to follow up specific outcomes -------------------------------------
Us.temp = Us[1,] - Us[2,]
Us.temp = rbind(Us.temp, Us[2,] - Us[1,])
Us.temp = rbind(Us.temp, Us.temp[1,] * Us.temp[2,])
Ind = which(Us.temp[3,]<0)
length(Ind)/length(Us[1,])

# Correlations ------------------------------------------------------------

crosscorr(chain.0)
crosscorr(chain.625)
crosscorr(chain.125)
min(crosscorr(chain.0),na.rm=TRUE)
min(crosscorr(chain.625),na.rm=TRUE)
min(crosscorr(chain.125),na.rm=TRUE)

# Outcomes for Experimental hut paper ------------------------------------------------------
# Absolute
Out.Exit2 =  rbind( c(median(mcmc(Exit_2[1,])), HDIofMCMC(mcmc(Exit_2[1,]),credMass = 0.95)),
        c(median(mcmc(Exit_2[2,])), HDIofMCMC(mcmc(Exit_2[2,]),credMass = 0.95)),
        c(median(mcmc(Exit_2[3,])), HDIofMCMC(mcmc(Exit_2[3,]),credMass = 0.95)))
row.names(Out.Exit2) = c('0','l','h');  colnames(Out.Exit2) = c('median','low','high') 

Out.Exit1 =  rbind( c(median(mcmc(Exit_1[1,])), HDIofMCMC(mcmc(Exit_1[1,]),credMass = 0.95)),
                    c(median(mcmc(Exit_1[2,])), HDIofMCMC(mcmc(Exit_1[2,]),credMass = 0.95)),
                    c(median(mcmc(Exit_1[3,])), HDIofMCMC(mcmc(Exit_1[3,]),credMass = 0.95)))
row.names(Out.Exit1) = c('0','l','h');  colnames(Out.Exit1) = c('median','low','high') 

Out.Exit0 =  rbind( c(median(mcmc(Exit_0[1,])), HDIofMCMC(mcmc(Exit_0[1,]),credMass = 0.95)),
                    c(median(mcmc(Exit_0[2,])), HDIofMCMC(mcmc(Exit_0[2,]),credMass = 0.95)),
                    c(median(mcmc(Exit_0[3,])), HDIofMCMC(mcmc(Exit_0[3,]),credMass = 0.95)))
row.names(Out.Exit0) = c('0','l','h');  colnames(Out.Exit0) = c('median','low','high') 

Out.P =  rbind( c(median(mcmc(PBorD[1,])), HDIofMCMC(mcmc(PBorD[1,]),credMass = 0.95)),
                    c(median(mcmc(PBorD[2,])), HDIofMCMC(mcmc(PBorD[2,]),credMass = 0.95)),
                    c(median(mcmc(PBorD[3,])), HDIofMCMC(mcmc(PBorD[3,]),credMass = 0.95)))
row.names(Out.P) = c('0','l','h');  colnames(Out.P) = c('median','low','high') 

Out.Kd2 =  rbind( c(median(mcmc(Kd_2[1,])), HDIofMCMC(mcmc(Kd_2[1,]),credMass = 0.95)),
                    c(median(mcmc(Kd_2[2,])), HDIofMCMC(mcmc(Kd_2[2,]),credMass = 0.95)),
                    c(median(mcmc(Kd_2[3,])), HDIofMCMC(mcmc(Kd_2[3,]),credMass = 0.95)))
row.names(Out.Kd2) = c('0','l','h');  colnames(Out.Kd2) = c('median','low','high') 

Out.Kd1 =  rbind( c(median(mcmc(Kd_1[1,])), HDIofMCMC(mcmc(Kd_1[1,]),credMass = 0.95)),
                    c(median(mcmc(Kd_1[2,])), HDIofMCMC(mcmc(Kd_1[2,]),credMass = 0.95)),
                    c(median(mcmc(Kd_1[3,])), HDIofMCMC(mcmc(Kd_1[3,]),credMass = 0.95)))
row.names(Out.Kd1) = c('0','l','h');  colnames(Out.Kd1) = c('median','low','high') 

Out.Kd0 =  rbind( c(median(mcmc(Kd_0[1,])), HDIofMCMC(mcmc(Kd_0[1,]),credMass = 0.95)),
                    c(median(mcmc(Kd_0[2,])), HDIofMCMC(mcmc(Kd_0[2,]),credMass = 0.95)),
                    c(median(mcmc(Kd_0[3,])), HDIofMCMC(mcmc(Kd_0[3,]),credMass = 0.95)))
row.names(Out.Kd0) = c('0','l','h');  colnames(Out.Kd0) = c('median','low','high') 

Out.Us =  rbind( c(median(mcmc(Us[1,])), HDIofMCMC(mcmc(Us[1,]),credMass = 0.95)),
                  c(median(mcmc(Us[2,])), HDIofMCMC(mcmc(Us[2,]),credMass = 0.95)),
                  c(median(mcmc(Us[3,])), HDIofMCMC(mcmc(Us[3,]),credMass = 0.95)))
row.names(Out.Us) = c('0','l','h');  colnames(Out.Us) = c('median','low','high') 

# Relative
Out.Exit2_rel =  rbind( c(median(mcmc(Exit_2_rel[1,])), HDIofMCMC(mcmc(Exit_2_rel[1,]),credMass = 0.95)),
                        c(median(mcmc(Exit_2_rel[2,])), HDIofMCMC(mcmc(Exit_2_rel[2,]),credMass = 0.95)),
                        c(median(mcmc(Exit_2_rel[3,])), HDIofMCMC(mcmc(Exit_2_rel[3,]),credMass = 0.95)))
row.names(Out.Exit2_rel) = c('0','l','h');  colnames(Out.Exit2_rel) = c('median','low','high') 

Out.Exit1_rel =  rbind( c(median(mcmc(Exit_1_rel[1,])), HDIofMCMC(mcmc(Exit_1_rel[1,]),credMass = 0.95)),
                        c(median(mcmc(Exit_1_rel[2,])), HDIofMCMC(mcmc(Exit_1_rel[2,]),credMass = 0.95)),
                        c(median(mcmc(Exit_1_rel[3,])), HDIofMCMC(mcmc(Exit_1_rel[3,]),credMass = 0.95)))
row.names(Out.Exit1_rel) = c('0','l','h');  colnames(Out.Exit1_rel) = c('median','low','high') 

Out.Exit0_rel =  rbind( c(median(mcmc(Exit_0_rel[1,])), HDIofMCMC(mcmc(Exit_0_rel[1,]),credMass = 0.95)),
                    c(median(mcmc(Exit_0_rel[2,])), HDIofMCMC(mcmc(Exit_0_rel[2,]),credMass = 0.95)),
                    c(median(mcmc(Exit_0_rel[3,])), HDIofMCMC(mcmc(Exit_0_rel[3,]),credMass = 0.95)))
row.names(Out.Exit0_rel) = c('0','l','h');  colnames(Out.Exit0_rel) = c('median','low','high') 

Out.P_rel =  rbind( c(median(mcmc(PBorD_rel[1,])), HDIofMCMC(mcmc(PBorD_rel[1,]),credMass = 0.95)),
                    c(median(mcmc(PBorD_rel[2,])), HDIofMCMC(mcmc(PBorD_rel[2,]),credMass = 0.95)),
                    c(median(mcmc(PBorD_rel[3,])), HDIofMCMC(mcmc(PBorD_rel[3,]),credMass = 0.95)))
row.names(Out.P_rel) = c('0','l','h');  colnames(Out.P_rel) = c('median','low','high') 

Out.Kd2_rel =  rbind( c(median(mcmc(Kd_2_rel[1,])), HDIofMCMC(mcmc(Kd_2_rel[1,]),credMass = 0.95)),
                      c(median(mcmc(Kd_2_rel[2,])), HDIofMCMC(mcmc(Kd_2_rel[2,]),credMass = 0.95)),
                      c(median(mcmc(Kd_2_rel[3,])), HDIofMCMC(mcmc(Kd_2_rel[3,]),credMass = 0.95)))
row.names(Out.Kd2_rel) = c('0','l','h');  colnames(Out.Kd2_rel) = c('median','low','high') 

Out.Kd1_rel =  rbind( c(median(mcmc(Kd_1_rel[1,])), HDIofMCMC(mcmc(Kd_1_rel[1,]),credMass = 0.95)),
                      c(median(mcmc(Kd_1_rel[2,])), HDIofMCMC(mcmc(Kd_1_rel[2,]),credMass = 0.95)),
                      c(median(mcmc(Kd_1_rel[3,])), HDIofMCMC(mcmc(Kd_1_rel[3,]),credMass = 0.95)))
row.names(Out.Kd1_rel) = c('0','l','h');  colnames(Out.Kd1_rel) = c('median','low','high') 

Out.Kd0_rel =  rbind( c(median(mcmc(Kd_0_rel[1,])), HDIofMCMC(mcmc(Kd_0_rel[1,]),credMass = 0.95)),
                      c(median(mcmc(Kd_0_rel[2,])), HDIofMCMC(mcmc(Kd_0_rel[2,]),credMass = 0.95)),
                      c(median(mcmc(Kd_0_rel[3,])), HDIofMCMC(mcmc(Kd_0_rel[3,]),credMass = 0.95)))
row.names(Out.Kd0_rel) = c('0','l','h');  colnames(Out.Kd0_rel) = c('median','low','high') 

Out.Us_rel =  rbind( c(median(mcmc(Us_rel[1,])), HDIofMCMC(mcmc(Us_rel[1,]),credMass = 0.95)),
                 c(median(mcmc(Us_rel[2,])), HDIofMCMC(mcmc(Us_rel[2,]),credMass = 0.95)),
                 c(median(mcmc(Us_rel[3,])), HDIofMCMC(mcmc(Us_rel[3,]),credMass = 0.95)))
row.names(Out.Us_rel) = c('0','l','h');  colnames(Out.Us_rel) = c('median','low','high') 

