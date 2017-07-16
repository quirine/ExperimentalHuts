# functions SR movement analysis

simulate.movement.data <- function(Mosquitoes.released, Timepoints,params,model = hut.movement.model, kd.interval.factor = 1,
                                   first.measurement = 30, day = 1, movement.tied.to.exit = FALSE){
  # browser()
  if (movement.tied.to.exit == TRUE){
    input.Qs <- params[1:3]; input.Rs <- params[6:8]
    params[1:3] <- input.Qs * (1 - input.Rs)
    params[6:8] <- input.Qs * (input.Rs)
  }
  
  K <- numeric();
  Rel_loc <- character(); 
  for (ii in 1:length(Huts)){
    Out = data.frame(ode(Init[ii,],t,model,parms=params),method="ode45")
    mos <- Mosquitoes.released
    for (jj in 2:length(Timepoints)){
      p.hut <- sum(Out[Timepoints[jj]/delta+1,2:6]) 
      p.ex <- ( as.numeric(Out[Timepoints[jj]/delta+1, 7:11]) -  as.numeric(Out[(Timepoints[jj]/delta) - (Interval/delta)+1,7:11]) ) 
      p.kd <- ( as.numeric(Out[Timepoints[jj]/delta+1, 12:16]) - as.numeric(Out[(Timepoints[jj]/delta) - (Interval/delta)+1,12:16]) ) 
      if (dim(Out)[2] == 22) { 
        p.u <- sum(as.numeric(Out[Timepoints[jj]/delta+1, 17:21]) - as.numeric(Out[(Timepoints[jj]/delta) - (Interval/delta)+1,17:21]))  
        p.huts <- Out[Timepoints[jj]/delta+1,2:6]
      }
      if (dim(Out)[2] == 17) { K.temp <- rmultinom(1,mos,c(p.hut,p.ex,p.kd)) 
                               K <- rbind(K, aperm(K.temp,c(2,1)))
      }
      else 
      { 
        K.temp <- rmultinom(1,mos,c(p.huts,p.u,p.ex,p.kd)) 
        K.temp <- c( sum(K.temp[1:5]), K.temp) 
        K <- rbind(K, K.temp)
      }
      
      Rel_loc <- c(Rel_loc, Huts[ii])
      mos <- K.temp[1]
      }
  }
  if (dim(Out)[2] == 17) { 

    colnames(K) <- c('k.hut', 'k.ex.A', 'k.ex.B', 'k.ex.C', 'k.ex.D', 'k.ex.E', 
                     'k.kd.A', 'k.kd.B', 'k.kd.C', 'k.kd.D', 'k.kd.E')
  }
  else 
  { 
    colnames(K) <- c('k.hut','k.hut.A','k.hut.B','k.hut.C','k.hut.D','k.hut.E', 'k.u',
                     'k.ex.A', 'k.ex.B', 'k.ex.C', 'k.ex.D', 'k.ex.E', 
                     'k.kd.A', 'k.kd.B', 'k.kd.C', 'k.kd.D', 'k.kd.E') 
  }
  Sim.Data <- data.frame(rep(Timepoints[-1],length(Huts)),K,Rel_loc,rep(0,length(Rel_loc)),rep(day,length(Rel_loc)))
  if (dim(Out)[2] == 17) { colnames(Sim.Data)[c(1,14,15)] <- c('Time', 'Dose','Exp.day') }
  else { colnames(Sim.Data)[c(1,20,21)] <- c('Time', 'Dose','Exp.day') }
  if (kd.interval.factor > 1) { 
    Ind <- which((Sim.Data$Time-first.measurement) %% (Interval*kd.interval.factor) == 0 & Sim.Data$Time > first.measurement ) # WATCH OUT: HARD CODED. SHOULD ALSO WORK FOR FACTORS > 2
    Sim.Data$k.kd.A[Ind] = Sim.Data$k.kd.A[Ind] + Sim.Data$k.kd.A[Ind-1];
    Sim.Data$k.kd.B[Ind] = Sim.Data$k.kd.B[Ind] + Sim.Data$k.kd.B[Ind-1]; 
    Sim.Data$k.kd.C[Ind] = Sim.Data$k.kd.C[Ind] + Sim.Data$k.kd.C[Ind-1]; 
    Sim.Data$k.kd.D[Ind] = Sim.Data$k.kd.D[Ind] + Sim.Data$k.kd.D[Ind-1]; 
    Sim.Data$k.kd.E[Ind] = Sim.Data$k.kd.E[Ind] + Sim.Data$k.kd.E[Ind-1]; 
    Sim.Data$k.huts[Ind-1] = Sim.Data$k.huts[Ind-1] + sum(Sim.Data$k.kd.A[Ind-1],Sim.Data$k.kd.B[Ind-1],
                                                          Sim.Data$k.kd.C[Ind-1],Sim.Data$k.kd.D[Ind-1],Sim.Data$k.kd.E[Ind-1] ) 
    Ind <- which((Sim.Data$Time-first.measurement) %% (Interval*kd.interval.factor) == 0 | Sim.Data$Time == first.measurement) 
    Sim.Data$k.kd.A[-Ind] = 0; Sim.Data$k.kd.B[-Ind] = 0; Sim.Data$k.kd.C[-Ind] = 0; Sim.Data$k.kd.D[-Ind] = 0; Sim.Data$k.kd.E[-Ind] = 0; 

  }
return(Sim.Data)
}

retrieve.and.attach <- function(Data, End = FALSE) {
  if (End == FALSE){
    K_exit = cbind(Data$k.ex.A, Data$k.ex.B, Data$k.ex.C, Data$k.ex.D, Data$k.ex.E )
    K_kd = cbind(Data$k.kd.A, Data$k.kd.B, Data$k.kd.C, Data$k.kd.D, Data$k.kd.E )
    k_huts = Data$k.bb
    interval = diff(Data$Time)[1]
    Rel_loc = factor(Data$Rel_loc)
    Locations <- Data$Rel_loc
  }
  else {
    K_exit.end = cbind(Data$k.ex.A.end, Data$k.ex.B.end, Data$k.ex.C.end, Data$k.ex.D.end, Data$k.ex.E.end )
    K_kd.end = cbind(Data$k.kd.A.end, Data$k.kd.B.end, Data$k.kd.C.end, Data$k.kd.D.end, Data$k.kd.E.end )
    k_u = Data$k.u
    K_huts = cbind(Data$k.hut.A, Data$k.hut.B, Data$k.hut.C, Data$k.hut.D, Data$k.hut.E )
    Rel_loc.end = factor(Data$Rel_loc.end)
    end.time = Data$Time.end[1]
  }
}
