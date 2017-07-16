# Load packages -----------------------------------------------------------
library(deSolve)
library(abind)
# library(lhs)
# library(tictoc)

# Set work directory ------------------------------------------------------
# setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# Source files ------------------------------------------------------------
#source('Parameters.LTFU.SR.Exit.R')
source('SR.Exit.Models.R')
source('SR.Exit.SimulateData.Functions.R')

# Prepare Inputs  ---------------------------------------------------------------
Timepoints <- seq(0,Duration,by=Interval)#seq(0,750,by=Interval)
# half.of.int <- 0.5 * Interval
t <- seq(0,24*60,by=delta)          # run for 1 day
Huts <- LETTERS[1:5]
Init <- diag(5) # create 5 sets of initial conditions, one for each release location
Init <- cbind(Init,matrix(0,5,15))

# Simulate Data -----------------------------------------------------------
ii = 1
Sim.Data <- simulate.movement.data(Mosquitoes, Timepoints,defaults,model = hut.movement.model.ltfu, 
                                   kd.interval.factor = kd.interval.factor, first.measurement = first.measurement,
                                   day = ii, movement.tied.to.exit = movement.tied.to.exit )
Ind.c <- which(Sim.Data$Rel_loc == "C")
Sim.Data <- Sim.Data[-Ind.c,]  # get rid of C-releases, only done in control experiment

Ind <- which(Sim.Data$Time == Duration)
Sim.Data.end <- Sim.Data[Ind,]

Ind.a <- which(Sim.Data$Rel_loc == "A");  Ind.a.end <- which(Sim.Data.end$Rel_loc == "A")
k.u = sum(Sim.Data$k.u[Ind.a]);           Sim.Data.end$k.u[Ind.a.end] <- k.u 
Ind.b <- which(Sim.Data$Rel_loc == "B");  Ind.b.end <- which(Sim.Data.end$Rel_loc == "B")
k.u = sum(Sim.Data$k.u[Ind.b]);           Sim.Data.end$k.u[Ind.b.end] <- k.u 
Ind.d <- which(Sim.Data$Rel_loc == "D");  Ind.d.end <- which(Sim.Data.end$Rel_loc == "D")
k.u = sum(Sim.Data$k.u[Ind.d]);           Sim.Data.end$k.u[Ind.d.end] <- k.u 
Ind.e <- which(Sim.Data$Rel_loc == "E");  Ind.e.end <- which(Sim.Data.end$Rel_loc == "E")
k.u = sum(Sim.Data$k.u[Ind.e]);           Sim.Data.end$k.u[Ind.e.end] <- k.u 

colnames(Sim.Data.end)[9:18] <- c('k.ex.A.end', 'k.ex.B.end', 'k.ex.C.end', 'k.ex.D.end', 'k.ex.E.end', 
                                  'k.kd.A.end', 'k.kd.B.end', 'k.kd.C.end', 'k.kd.D.end', 'k.kd.E.end')
colnames(Sim.Data.end)[1] <- 'Time.end'
colnames(Sim.Data.end)[19] <- 'Rel_loc.end'

Sim.Data <- Sim.Data[-Ind,]
Sim.Data$k.bb <- Sim.Data$k.hut + Sim.Data$k.u
if (num.exp.days > 1){
  # print('Check')
  for (ii in 2:num.exp.days){ 
    Sim.Data.temp <- simulate.movement.data(Mosquitoes, Timepoints,defaults,model = hut.movement.model.ltfu, 
                                            kd.interval.factor = kd.interval.factor, first.measurement = first.measurement,
                                            day = ii, movement.tied.to.exit = movement.tied.to.exit )
    Ind.c <- which(Sim.Data.temp$Rel_loc == "C")
    Sim.Data.temp <- Sim.Data.temp[-Ind.c,]  # get rid of C-releases, only done in control experiment
    
    Ind <- which(Sim.Data.temp$Time == Duration)
    Sim.Data.end.temp <- Sim.Data.temp[Ind,]
    
    Ind.a <- which(Sim.Data.temp$Rel_loc == "A");  Ind.a.end <- which(Sim.Data.end.temp$Rel_loc == "A")
    k.u = sum(Sim.Data.temp$k.u[Ind.a]);           Sim.Data.end.temp$k.u[Ind.a.end] <- k.u 
    Ind.b <- which(Sim.Data.temp$Rel_loc == "B");  Ind.b.end <- which(Sim.Data.end.temp$Rel_loc == "B")
    k.u = sum(Sim.Data.temp$k.u[Ind.b]);           Sim.Data.end.temp$k.u[Ind.b.end] <- k.u 
    Ind.d <- which(Sim.Data.temp$Rel_loc == "D");  Ind.d.end <- which(Sim.Data.end.temp$Rel_loc == "D")
    k.u = sum(Sim.Data.temp$k.u[Ind.d]);           Sim.Data.end.temp$k.u[Ind.d.end] <- k.u 
    Ind.e <- which(Sim.Data.temp$Rel_loc == "E");  Ind.e.end <- which(Sim.Data.end.temp$Rel_loc == "E")
    k.u = sum(Sim.Data.temp$k.u[Ind.e]);           Sim.Data.end.temp$k.u[Ind.e.end] <- k.u 
    
    colnames(Sim.Data.end.temp)[9:18] <- c('k.ex.A.end', 'k.ex.B.end', 'k.ex.C.end', 'k.ex.D.end', 'k.ex.E.end', 
                                           'k.kd.A.end', 'k.kd.B.end', 'k.kd.C.end', 'k.kd.D.end', 'k.kd.E.end')
    colnames(Sim.Data.end.temp)[1] <- 'Time.end'
    colnames(Sim.Data.end.temp)[19] <- 'Rel_loc.end'
    
    Sim.Data.temp <- Sim.Data.temp[-Ind,]
    Sim.Data.temp$k.bb <- Sim.Data.temp$k.hut + Sim.Data.temp$k.u
    
    Sim.Data = rbind(Sim.Data,Sim.Data.temp)
    Sim.Data.end = rbind(Sim.Data.end,Sim.Data.end.temp)
  }
}

# Retrieve relevant vectors and attach ------------------------------------
# attach(Sim.Data)
K_exit = cbind(Sim.Data$k.ex.A, Sim.Data$k.ex.B, Sim.Data$k.ex.C, Sim.Data$k.ex.D, Sim.Data$k.ex.E )
K_kd = cbind(Sim.Data$k.kd.A, Sim.Data$k.kd.B, Sim.Data$k.kd.C, Sim.Data$k.kd.D, Sim.Data$k.kd.E )
k_huts = Sim.Data$k.bb
interval = diff(Sim.Data$Time)[1]
Init <- diag(5) # create 5 sets of initial conditions, one for each release location
Init <- cbind(Init,matrix(0,5,10))
Rel_loc = factor(Sim.Data$Rel_loc)
Locations <- Sim.Data$Rel_loc

# attach(Sim.Data.end)
K_exit.end = cbind(Sim.Data.end$k.ex.A.end, Sim.Data.end$k.ex.B.end, Sim.Data.end$k.ex.C.end, Sim.Data.end$k.ex.D.end, Sim.Data.end$k.ex.E.end )
K_kd.end = cbind(Sim.Data.end$k.kd.A.end, Sim.Data.end$k.kd.B.end, Sim.Data.end$k.kd.C.end, Sim.Data.end$k.kd.D.end, Sim.Data.end$k.kd.E.end )
k_u = Sim.Data.end$k.u
K_huts = cbind(Sim.Data.end$k.hut.A, Sim.Data.end$k.hut.B, Sim.Data.end$k.hut.C, Sim.Data.end$k.hut.D, Sim.Data.end$k.hut.E )
Rel_loc.end = factor(Sim.Data.end$Rel_loc.end)
end.time = Sim.Data.end$Time.end[1]





