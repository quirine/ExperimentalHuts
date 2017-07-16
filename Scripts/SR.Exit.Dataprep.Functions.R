
Prepare.data <- function(Dosages, Huts, Exp.day, ExitSR){

  Data.exit <- array(data=NA, dim=c(1,5)); Data.kd <- array(data=NA, dim=c(1,5));  Rel_loc <- character(); K.hut <- numeric();
  Dose <- numeric(); Day <- numeric();  Time <- numeric(); 
  for (ii in Dosages) {
    for (jj in Huts)  {
      for (kk in Exp.day)  {
        Total <- num.mos;
        Ind <- which(ExitSR$Loc_release==jj & ExitSR$Dose==ii & ExitSR$Day == kk);
        Timepoints <- levels(as.factor(ExitSR$TTE[Ind]))
        for (ll in Timepoints){
          Ind.Exit <- which(ExitSR$TTE==as.numeric(ll) & ExitSR$Event=="exit" & ExitSR$Loc_release==jj & ExitSR$Dose==ii & ExitSR$Day == kk)
          if( length(Ind.Exit) >0 ) {
            Data.temp <- Raw.counts (Ind.Exit,ExitSR) ;  Data.exit <- rbind(Data.exit,Data.temp[-1])
            exit <- sum(Data.temp[-1])
          } else { 
            Data.exit <- rbind(Data.exit,rep(0,length(Huts)))
            exit <- 0
          }
          Ind.KD <- which(ExitSR$TTE==as.numeric(ll) & ExitSR$Event=="KD" & ExitSR$Loc_release==jj & ExitSR$Dose==ii & ExitSR$Day == kk)
          if( length(Ind.KD) >0 ) {
            Data.temp <- Raw.counts (Ind.KD,ExitSR) ;  Data.kd <- rbind(Data.kd,Data.temp[-1])  
            kd <- sum(Data.temp[-1])
          } else { 
            Data.kd <- rbind(Data.kd,rep(0,length(Huts))) 
            kd <- 0
          }
          loc.temp <- rep(jj,1) ;  Rel_loc <- c(Rel_loc,loc.temp);
          dose.temp <- rep(ii,1);  Dose <- c(Dose,dose.temp) ;
          day.temp <- rep(kk,1); Day <- c(Day,day.temp);
          time.temp <- rep(ll,1); Time <- c(Time,time.temp); 
          event.tot <- exit + kd; 
          Total <- Total - event.tot; K.hut <- c(K.hut,Total)
          
        }
      }
    }
  }
  
  
  Data <- data.frame(Time,K.hut,Data.exit[-1,],Data.kd[-1,], Rel_loc,Dose, Day)
  Data <- upData(Data, rename = list(K.hut = "k.hut", X1 = 'k.ex.A', X2 = 'k.ex.B', X3 = 'k.ex.C', X4 = 'k.ex.D', X5 = 'k.ex.E',
                                     X1.1 = 'k.kd.A', X2.1 = 'k.kd.B', X3.1 = 'k.kd.C', X4.1 = 'k.kd.D', X5.1 = 'k.kd.E'))
  rm(Day, Time, Dose, Rel_loc, K.hut)
  Data$Time <- as.numeric(as.character(Data$Time))
  Ind <- which(Data$k.hut<0)
  Data$k.hut[Ind] <- 0
  return(Data)
}


# Add Censored Mosquitoes -------------------------------------------------

add.censored.mosquitoes <- function(ExitSR,num.mos,Dosages,Exp.day,ReleaseTime){
  
  # get number of dosages and experiment days to loop through
  Time <- character(); TTE <- character(); Event <- numeric(); Fed <- character(); Loc_release <- character(); ID <- numeric(); Dose <- numeric();
  Day <- numeric(); Loc_event <-character() ; ExitKD <-logical(); Loc_release <- character(); Loc_event <- character(); Treat <- character();
  
  for (ii in 1:length(Dosages)){
    Ind <- which(ExitSR$Dose==Dosages[ii])
    Temp <- ExitSR[Ind,]
    for (jj in 1:length(Exp.day)){
      Ind <- which(Temp$Day==Exp.day[jj])
      Data <- Temp[Ind,]
      Num.caught <- summary(Data$Loc_release)
      Num.censored <- c(num.mos,num.mos,0,num.mos,num.mos) - Num.caught[1:5] # HARD CODED, assuming 5 huts
      if (Dosages[ii]==0) {
        Num.censored <- num.mos - Num.caught[1:5] # HARD CODED, assuming 5 huts
      }
      print(Num.censored)
      Names <- names(Num.censored)
      ind <- which(Num.censored > 0) 
      Num.censored.pos <- Num.censored[ind]
      CS.tot <- sum(Num.censored.pos)
      ID <- c(ID,seq(max(Data$ID)+1,max(Data$ID)+CS.tot))
      Time.temp <- rep("18:00:00.00",CS.tot); as.character(Time.temp); Time <- c(Time,Time.temp); Time.temp <- strptime(Time.temp, format = '%H:%M:%S');
      TTE <- c(TTE,as.numeric(difftime(Time.temp,ReleaseTime,units='mins')))
      Event <- c(Event,rep("CS",CS.tot))
      ExitKD <- c(ExitKD,rep("NA",CS.tot))
      Loc_event <- c(Loc_event,rep("NA",CS.tot))
      Treat <- c(Treat,rep("SR",CS.tot))
      Dose <- c(Dose,rep(Dosages[ii],CS.tot))
      Day <- c(Day,rep(Exp.day[jj],CS.tot))
      Fed <- c(Fed,rep("NA",CS.tot))
      for (kk in 1:length(Names)){
        if(Num.censored[kk] > 0){
          Lcr <- rep(Names[kk],Num.censored[kk])
          Loc_release <- c(Loc_release,Lcr)
        }
      }
    }
  }
  ExitSR.cens <- data.frame(ID,Time,TTE,Event,ExitKD,Loc_release,Loc_event,Treat,Dose,Day,Fed)  # took out Loc_release
  ExitSR <- rbind(ExitSR, ExitSR.cens, make.row.names = FALSE)
  return(ExitSR)    
}



# Normalize Data ----------------------------------------------------------

Normalized.counts <- function(Ind,ExitSR){
  Rel <- ExitSR[Ind,]
  t <- table(Rel$TTE)
  t <- as.numeric(names(t))
  Events <- matrix(data=NA,nrow=length(t),ncol = 5)
  for (ii in 1:length(t)){
    Ind <- which(Rel$TTE==t[ii])
    temp <- summary(Rel$Loc_event[Ind])
    temp <- temp[1:5]
    Events[ii,] <- temp
  }
  Events.norm <- Events/rowSums(Events)
  Data <- cbind(t,Events.norm)
  return(Data)
}


# Non-normalized counts ---------------------------------------------------

Raw.counts <- function(Ind,ExitSR){
  Rel <- ExitSR[Ind,]
  t <- table(Rel$TTE)
  t <- as.numeric(names(t))
  Events <- matrix(data=NA,nrow=length(t),ncol = 5)
  for (ii in 1:length(t)){
    Ind <- which(Rel$TTE==t[ii])
    temp <- summary(Rel$Loc_event[Ind])
    temp <- temp[1:5]
    Events[ii,] <- temp
  }
  Data <- cbind(t,Events)
  return(Data)
}


# Normalize Data per event locations ----------------------------------------------------------

Normalized.counts.event <- function(Ind,ExitSR){
  Rel <- ExitSR[Ind,]
  t <- table(Rel$TTE)
  t <- as.numeric(names(t))
  Events <- matrix(data=NA,nrow=length(t),ncol = 5)
  for (ii in 1:length(t)){
    Ind <- which(Rel$TTE==t[ii])
    temp <- summary(Rel$Loc_release[Ind])
    temp <- temp[1:5]
    Events[ii,] <- temp
  }
  Events.norm <- Events/rowSums(Events)
  Data <- cbind(t,Events.norm)
  return(Data)
}


# Count Data per event locations -------------------------------------------------------

Raw.counts.event <- function(Ind,ExitSR){
  Rel <- ExitSR[Ind,]
  t <- table(Rel$TTE)
  t <- as.numeric(names(t))
  Events <- matrix(data=NA,nrow=length(t),ncol = 5)
  for (ii in 1:length(t)){
    Ind <- which(Rel$TTE==t[ii])
    temp <- summary(Rel$Loc_release[Ind])
    temp <- temp[1:5]
    Events[ii,] <- temp
  }
  Data <- cbind(t,Events)
  return(Data)
}

