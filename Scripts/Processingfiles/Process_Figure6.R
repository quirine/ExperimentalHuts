# This code uses the posteriors from the data fits to run the ode-model and 
# estimate what proportion of time is spent in each hut. 
# This takes a little while. 
# Alternatively, you can load 'timespent.fixedr.RData', which is present in the
# repository and go straight to the section 'double plot of total and alive time (for manuscript)'
# to create the plot 

rm(list=ls())
# load libraries ----------------------------------------------------------
library(coda)
library(IDPmisc)
library(deSolve)
library(ggplot2)
library(pomp)

# set wd ------------------------------------------------------------------
setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# load source files -------------------------------------------------------
source('SR.Exit.Models.R')

# set data directory ------------------------------------------------------
data.path <- "../Output/Final/output.ar.ltfudata.relsd.gamma.beta.fixedr_2016-12-13/"

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

delta = 0.1
Init = c(rep(0.2,5),rep(0,15))
Time <- seq(0,750,by=delta)
movement.tied.to.exit = TRUE
load('timespent.fixedr.RData')

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

setwd(cd)

# Run ode-model for selection of 500 parameter sets -----------------------
BurnIn = 10000
n = 2000
Prop.a.0 <- Prop.a.l <- Prop.a.h <- numeric()
Prop.b.0 <- Prop.b.l <- Prop.b.h <- numeric()
Prop.c.0 <- Prop.c.l <- Prop.c.h <- numeric()
Prop.tot.a.0 <- Prop.tot.a.l <- Prop.tot.a.h <- numeric()
Prop.tot.b.0 <- Prop.tot.b.l <- Prop.tot.b.h <- numeric()
Prop.tot.c.0 <- Prop.tot.c.l <- Prop.tot.c.h <- numeric()

Temp = combinedchains.raw.0[[3]][BurnIn:summary(combinedchains.raw.0)$end,]
Ind = round(runif(n,min=1, max = dim(Temp)[1]))
for (ii in 1: n) {
  print(ii)
  parameters = c(Temp[Ind[ii],1:4],0.5,Temp[Ind[ii],5:11])
  if (movement.tied.to.exit == TRUE){
    temp.Q <- parameters[1:3]; temp.R <- parameters[6:8]
    parameters[1:3] <- temp.Q * (1 - temp.R)
    parameters[6:8] <- temp.R * temp.Q
  }
  out = ode(Init,Time,hut.movement.model.ltfu,parms=parameters)
  temp.a = sum(out[,2]) / sum(out[,2:6])
  temp.b = sum(out[,3]) / sum(out[,2:6])
  temp.c = sum(out[,4]) / sum(out[,2:6])
  Prop.a.0 = c(Prop.a.0, temp.a); Prop.b.0 = c(Prop.b.0, temp.b); Prop.c.0 = c(Prop.c.0, temp.c)  
  temp.a = sum(out[,2]) / sum(out[,-1])
  temp.b = sum(out[,3]) / sum(out[,-1])
  temp.c = sum(out[,4]) / sum(out[,-1])
  Prop.tot.a.0 = c(Prop.tot.a.0, temp.a); Prop.tot.b.0 = c(Prop.tot.b.0, temp.b); Prop.tot.c.0 = c(Prop.tot.c.0, temp.c) 
}

Temp = combinedchains.raw.l[[3]][BurnIn:summary(combinedchains.raw.l)$end,]
Ind = round(runif(n,min=1, max = dim(Temp)[1]))
for (ii in 1: n) {
  print(ii)
  parameters = c(Temp[Ind[ii],1:4],0.5,Temp[Ind[ii],5:11])
  if (movement.tied.to.exit == TRUE){
    temp.Q <- parameters[1:3]; temp.R <- parameters[6:8]
    parameters[1:3] <- temp.Q * (1 - temp.R)
    parameters[6:8] <- temp.R * temp.Q
  }
  out = ode(Init,Time,hut.movement.model.ltfu,parms=parameters)
  temp.a = sum(out[,2]) / sum(out[,2:6])
  temp.b = sum(out[,3]) / sum(out[,2:6])
  temp.c = sum(out[,4]) / sum(out[,2:6])
  Prop.a.l = c(Prop.a.l, temp.a); Prop.b.l = c(Prop.b.l, temp.b); Prop.c.l = c(Prop.c.l, temp.c)  
  temp.a = sum(out[,2]) / sum(out[,-1])
  temp.b = sum(out[,3]) / sum(out[,-1])
  temp.c = sum(out[,4]) / sum(out[,-1])
  Prop.tot.a.l = c(Prop.tot.a.l, temp.a); Prop.tot.b.l = c(Prop.tot.b.l, temp.b); Prop.tot.c.l = c(Prop.tot.c.l, temp.c) 
}

Temp = combinedchains.raw.h[[3]][BurnIn:summary(combinedchains.raw.h)$end,]
Ind = round(runif(n,min=1, max = dim(Temp)[1]))
for (ii in 1: n) {
  print(ii)
  parameters = c(Temp[Ind[ii],1:4],0.5,Temp[Ind[ii],5:11])
  if (movement.tied.to.exit == TRUE){
    temp.Q <- parameters[1:3]; temp.R <- parameters[6:8]
    parameters[1:3] <- temp.Q * (1 - temp.R)
    parameters[6:8] <- temp.R * temp.Q
  }
  out = ode(Init,Time,hut.movement.model.ltfu,parms=parameters)
  temp.a = sum(out[,2]) / sum(out[,2:6])
  temp.b = sum(out[,3]) / sum(out[,2:6])
  temp.c = sum(out[,4]) / sum(out[,2:6])
  Prop.a.h = c(Prop.a.h, temp.a); Prop.b.h = c(Prop.b.h, temp.b); Prop.c.h = c(Prop.c.h, temp.c)
  temp.a = sum(out[,2]) / sum(out[,-1])
  temp.b = sum(out[,3]) / sum(out[,-1])
  temp.c = sum(out[,4]) / sum(out[,-1])
  Prop.tot.a.h = c(Prop.tot.a.h, temp.a); Prop.tot.b.h = c(Prop.tot.b.h, temp.b); Prop.tot.c.h = c(Prop.tot.c.h, temp.c) 
}

# double plot of total and alive time (for manuscript) -------------------------------------

# load('timespent.fixedr.RData')
setwd(figure.path)


transparancy = 0.5
pdf('Figure6.pdf',width=4.75,height=5)

par(mfrow = c(3,2),mai=c(0.2,.4,0,0.2) ,oma=c(3,3,3,.5))  

p1 = density(Prop.a.0[runif(1e5,min=1,max=n)] / Prop.a.0[runif(1e5,min=1,max=n)]); p1$y <- p1$y*diff(p1$x)[1]
p2 = density(Prop.a.l[runif(1e5,min=1,max=n)] / Prop.a.0[runif(1e5,min=1,max=n)]); p2$y <- p2$y*diff(p1$x)[1]
p3 = density(Prop.a.h[runif(1e5,min=1,max=n)] / Prop.a.0[runif(1e5,min=1,max=n)]); p3$y <- p3$y*diff(p1$x)[1]
plot(p1,col='white',main='',xlim=c(0.4,1.8),ylim=c(0,0.9e-2),xaxs='i',yaxs='i',las=1, pty="s",xaxt='n')
polygon( p1, col=alpha('black',transparancy),main='',border = 'black')  
polygon( p2, col=alpha('darkorange',transparancy),border = 'darkorange')  
polygon( p3, col=alpha('deeppink',transparancy),border = 'deeppink')  
abline(v = 1, col = 'gray', lty = 2, lwd = 2)
mtext(side = 3, 'Time spent before exit, kd, or ltfu', outer = FALSE, line = 1,cex=0.7)
mtext( bquote(bold("A")), side = 3, line = -2 ,at = 0.5,cex=1) 

p1 = density(Prop.tot.a.0[runif(1e5,min=1,max=n)] / Prop.tot.a.0[runif(1e5,min=1,max=n)]); p1$y <- p1$y*diff(p1$x)[1]
p2 = density(Prop.tot.a.l[runif(1e5,min=1,max=n)] / Prop.tot.a.0[runif(1e5,min=1,max=n)]); p2$y <- p2$y*diff(p1$x)[1]
p3 = density(Prop.tot.a.h[runif(1e5,min=1,max=n)] / Prop.tot.a.0[runif(1e5,min=1,max=n)]); p3$y <- p3$y*diff(p1$x)[1]
plot(p1,col='white',main='',xlim=c(0.4,1.8),ylim=c(0,0.9e-2),xaxs='i',yaxs='i',las=1, pty="s",yaxt='n',xaxt='n')
polygon( p1, col=alpha('black',transparancy),main='',border = 'black')  
polygon( p2, col=alpha('darkorange',transparancy),border = 'darkorange')  
polygon( p3, col=alpha('deeppink',transparancy),border = 'deeppink')  
abline(v = 1, col = 'gray', lty = 2, lwd = 2)
mtext(side = 3, 'Time spent during experiment', outer = FALSE, line = 1,cex=0.7)
mtext(side = 4, expression('H' [2]*''), outer = FALSE, line = 0.3,cex=0.7)
mtext( bquote(bold("B")), side = 3, line = -2 ,at = 0.5,cex=1)
legend('topright',c('control','low','high'), fill=c(alpha('black',transparancy),alpha('darkorange',transparancy),alpha('deeppink',transparancy)),box.lty=0)

p1 = density(Prop.b.0[runif(1e5,min=1,max=n)] / Prop.b.0[runif(1e5,min=1,max=n)]); p1$y <- p1$y*diff(p1$x)[1]
p2 = density(Prop.b.l[runif(1e5,min=1,max=n)] / Prop.b.0[runif(1e5,min=1,max=n)]); p2$y <- p2$y*diff(p1$x)[1]
p3 = density(Prop.b.h[runif(1e5,min=1,max=n)] / Prop.b.0[runif(1e5,min=1,max=n)]); p3$y <- p3$y*diff(p1$x)[1]
plot(p1,col='white',main='',xlim=c(0.4,1.8),ylim=c(0,0.9e-2),xaxs='i',yaxs='i',las=1, pty="s",xaxt='n')
polygon( p1, col=alpha('black',transparancy),main='',border = 'black')  
polygon( p2, col=alpha('darkorange',transparancy),border = 'darkorange')  
polygon( p3, col=alpha('deeppink',transparancy),border = 'deeppink')  
abline(v = 1, col = 'gray', lty = 2, lwd = 2)
mtext( bquote(bold("C")), side = 3, line = -2 ,at = 0.5,cex=1)

p1 = density(Prop.tot.b.0[runif(1e5,min=1,max=n)] / Prop.tot.b.0[runif(1e5,min=1,max=n)]); p1$y <- p1$y*diff(p1$x)[1]
p2 = density(Prop.tot.b.l[runif(1e5,min=1,max=n)] / Prop.tot.b.0[runif(1e5,min=1,max=n)]); p2$y <- p2$y*diff(p1$x)[1]
p3 = density(Prop.tot.b.h[runif(1e5,min=1,max=n)] / Prop.tot.b.0[runif(1e5,min=1,max=n)]); p3$y <- p3$y*diff(p1$x)[1]
plot(p1,col='white',main='',xlim=c(0.4,1.8),ylim=c(0,0.9e-2),xaxs='i',yaxs='i',las=1, pty="s",yaxt='n',xaxt='n')
polygon( p1, col=alpha('black',transparancy),main='',border = 'black')  
polygon( p2, col=alpha('darkorange',transparancy),border = 'darkorange')  
polygon( p3, col=alpha('deeppink',transparancy),border = 'deeppink')  
abline(v = 1, col = 'gray', lty = 2, lwd = 2)
mtext(side = 4, expression('H' [1]*''), outer = FALSE, line = 0.3,cex=0.7)
mtext( bquote(bold("D")), side = 3, line = -2 ,at = 0.5,cex=1)


p1 = density(Prop.c.0[runif(1e5,min=1,max=n)] / Prop.c.0[runif(1e5,min=1,max=n)]); p1$y <- p1$y*diff(p1$x)[1]
p2 = density(Prop.c.l[runif(1e5,min=1,max=n)] / Prop.c.0[runif(1e5,min=1,max=n)]); p2$y <- p2$y*diff(p1$x)[1]
p3 = density(Prop.c.h[runif(1e5,min=1,max=n)] / Prop.c.0[runif(1e5,min=1,max=n)]); p3$y <- p3$y*diff(p1$x)[1]
plot(p1,col='white',main='',xlim=c(0.4,1.8),ylim=c(0,0.9e-2),xaxs='i',yaxs='i',las=1, pty="s")
polygon( p1, col=alpha('black',transparancy),main='',border = 'black')  
polygon( p2, col=alpha('darkorange',transparancy),border = 'darkorange')  
polygon( p3, col=alpha('deeppink',transparancy),border = 'deeppink')  

abline(v = 1, col = 'gray', lty = 2, lwd = 2)
mtext( bquote(bold("E")), side = 3, line = -2 ,at = 0.5,cex=1)

p1 = density(Prop.tot.c.0[runif(1e5,min=1,max=n)] / Prop.tot.c.0[runif(1e5,min=1,max=n)]); p1$y <- p1$y*diff(p1$x)[1]
p2 = density(Prop.tot.c.l[runif(1e5,min=1,max=n)] / Prop.tot.c.0[runif(1e5,min=1,max=n)]); p2$y <- p2$y*diff(p1$x)[1]
p3 = density(Prop.tot.c.h[runif(1e5,min=1,max=n)] / Prop.tot.c.0[runif(1e5,min=1,max=n)]); p3$y <- p3$y*diff(p1$x)[1]
plot(p1,col='white',main='',xlim=c(0.4,1.8),ylim=c(0,0.9e-2),xaxs='i',yaxs='i',las=1, pty="s",yaxt='n')
polygon( p1, col=alpha('black',transparancy),main='',border = 'black')  
polygon( p2, col=alpha('darkorange',transparancy),border = 'darkorange')  
polygon( p3, col=alpha('deeppink',transparancy),border = 'deeppink')  
abline(v = 1, col = 'gray', lty = 2, lwd = 2)
mtext(side = 4, expression('H' [T]*''), outer = FALSE, line = 0.3,cex=0.7)
mtext( bquote(bold("F")), side = 3, line = -2 ,at = 0.5,cex=1)

mtext(side = 1, 'Proportion relative to control', outer = TRUE, line = 1,cex = 0.7)
mtext(side = 2, 'Density', outer = TRUE, line = 1,cex = 0.7)

dev.off()
setwd(cd)

