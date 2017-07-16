rm(list=ls())

# Load packages -----------------------------------------------------------
library(deSolve)
library(emdbook)
library(bbmle)
library(abind)

# Set work directory ------------------------------------------------------
setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# Load source files -------------------------------------------------------
source('SR.Exit.Models.R')
source('SR.Exit.FigureFunctions.R')

# Set directory to save figures -------------------------------------------
figure.path <- '../Figures/Illustrations'
cd <- getwd()

# Define parameters -------------------------------------------------------
Init <- diag(5) # create 5 sets of initial conditions, one for each release location
Init <- cbind(Init,matrix(0,5,15))
Time <- seq(1,12*60,by=0.1)          # run for 1 day
D <-c(0,0.0625,0.125)

# Store parameter table ---------------------------------------------------
Qs <- c(0.05, 0.05, 0.05)     # q_AorE, q_BorD, q_C: movement rates between huts
Ps <- c(0.6, 0.5)                # p_BorD, p_C (p_AorE is always 1 and thus not defined) : probability of moving away from SR-hut (C)
Rs <- c(0.005, 0.005, 0.005)     #r_AorE, r_BorD, r_C : exit rates
Ks <- c(0.0001, 0.0001, 0.0001)  #k_AorE, k_BorD, k_C : kd rates
u <- 7e-4                        # loss to follow up rate

# Run the model -----------------------------------------------------------
Init.4 <- Init/4                 # this is assuming you release a quarter of the total in all but the treatment hut
Init.4[3,3] <- 0
someData <- rep(NaN, length(Time)*21);  
Out <- array(someData,c(length(Time),21));   

Out = abind(Out,ode(Init.4[1,],Time,hut.movement.model.ltfu,parms=c(Qs,Ps,Rs,Ks,u)),along = 3)
Out = abind(Out,ode(Init.4[2,],Time,hut.movement.model.ltfu,parms=c(Qs,Ps,Rs,Ks,u)),along = 3)
Out = abind(Out,ode(Init.4[3,],Time,hut.movement.model.ltfu,parms=c(Qs,Ps,Rs,Ks,u)),along = 3)
Out = abind(Out,ode(Init.4[4,],Time,hut.movement.model.ltfu,parms=c(Qs,Ps,Rs,Ks,u)),along = 3)
Out = abind(Out,ode(Init.4[5,],Time,hut.movement.model.ltfu,parms=c(Qs,Ps,Rs,Ks,u)),along = 3)
Out <- Out[,,-1]

# Polygon plots normalised to be 1 in each hut -----------------------------------------------------------
cd <- getwd()
setwd(figure.path)

pdf('SR.Exit_HutPolygons_norm.pdf',width=4.5,height=3)
old.par <- par(mfrow=c(1,5),oma = c(5,4,4,2) + 0.1,
               mar = c(0,0,1,0) + 0.1)

stackplot.norm(Time,Out = Out,ylim=c(0,1)); par(xaxt="n")

mtext(text = "Time (min)", side = 1, line = 2, outer = TRUE )
mtext(text = "Proportion", side = 2, line = 2, outer = TRUE )
xaxt="y"
par(old.par)

dev.off()
setwd(cd)

# Polygon plots presence hut-----------------------------------------------------------
cd <- getwd()
setwd(figure.path)

pdf('SR.Exit_HutPolygons_presence.pdf',width=4.5,height=4.5)
old.par <- par(mfrow=c(1,5),oma = c(5,4,0,0) + 0.1,
               mar = c(0,0,1,1) + 0.1)
stackplot.loc(Time,Out = Out,ylim=c(0,0.4)); par(xaxt="n")

mtext(text = "Time (min)", side = 1, line = 2, outer = TRUE )
mtext(text = "Proportion", side = 2, line = 2, outer = TRUE )
xaxt="y"
par(old.par)

dev.off()
setwd(cd)

# Polygon plots exit -----------------------------------------------------------
cd <- getwd()
setwd(figure.path)

pdf('SR.Exit_HutPolygons_exit.pdf',width=4.5,height=4.5)
old.par <- par(mfrow=c(1,5),oma = c(5,4,0,0) + 0.1,
               mar = c(0,0,1,1) + 0.1)
stackplot.exit(Time,Out = Out,ylim=c(0,0.4))

mtext(text = "Time (min)", side = 1, line = 2, outer = TRUE )
mtext(text = "Proportion", side = 2, line = 2, outer = TRUE )
xaxt="y"
par(old.par)

dev.off()
setwd(cd)

# Polygon plots kd -----------------------------------------------------------
cd <- getwd()
setwd(figure.path)

pdf('SR.Exit_HutPolygons_kd.pdf',width=4.5,height=4.5)
old.par <- par(mfrow=c(1,5),oma = c(5,4,0,0) + 0.1,
               mar = c(0,0,1,1) + 0.1)
stackplot.kd(Time,Out = Out,ylim=c(0,0.03))

mtext(text = "Time (min)", side = 1, line = 2, outer = TRUE )
mtext(text = "Proportion", side = 2, line = 2, outer = TRUE )
xaxt="y"

dev.off()
setwd(cd)
