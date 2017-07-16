# Copyright for adaptive sampling functions (see underneath)
# Author:
#   Copyright 2003-05 Korbinian Strimmer
#   Rank, condition, and positive definiteness of a matrix
#   GNU General Public License, Version 2

# assuming a single chain
run_metropolis_MCMC <- function(startvalue, iterations, t, draw.probabilities, filename = 'temp.csv', writing.freq = 100, Qfix = FALSE, sd = 1, 
                                one.at.a.time = TRUE, relative.sds = FALSE, adaptive = FALSE, adapt.par = NULL, kd.interval.factor = 1, 
                                LTFU = FALSE, debugging.plot = FALSE, p.prior = 'Flat', q.prior = 'Flat',dose = 0, first.measurement = 30, movement.tied.to.exit = FALSE) {
  # browser()
  if (debugging.plot == TRUE) {  plot(seq(0,1,by=1)); plot.count = 1 }
  if(relative.sds == TRUE && adaptive == TRUE) { print('ERROR: adaptive and relative sds do not go together'); break }
  if(LTFU == TRUE && length(startvalue) == 11) { print('ERROR: your parameter vector is too short. Did you forget the u-parameter for the spider-compartment?'); break }
  if(adaptive == FALSE){ chain <- chain.prop <- array(dim = c(writing.freq+1,length(startvalue))) }
  else if (adaptive == TRUE) { chain <- chain.prop <- array(dim = c(iterations+1,length(startvalue))) }
  chain[1,] <- chain.prop[1,] <- startvalue
  LLs = numeric()
  column_names <- c(names(startvalue),'LLs',names(startvalue))
  write.table(array(NA,dim = c(1,1+length(startvalue)*2)), file = filename, row.names = FALSE, 
              append = FALSE, col.names = column_names, sep = ", ", quote = TRUE)
  #write.table(array(NA,dim = c(1,2*length(startvalue)+1)), file = filename, row.names = FALSE, 
   #           append = FALSE, col.names = column_names, sep = ", ", quote = TRUE)
  param.vector <- seq(1,length(startvalue))
  kk = 0
  if (LTFU == FALSE) { former.ll <- posterior.hut.movement(startvalue,t,interval,K_exit,K_kd,k_huts,Rel_loc,delta,kd.interval.factor,p.prior,q.prior,dose,first.measurement,movement.tied.to.exit) }
  if (LTFU == TRUE) {  former.ll <- posterior.hut.movement(startvalue,t,interval,K_exit,K_kd,k_huts, Rel_loc,delta,kd.interval.factor, 
                                                           K_exit.end,K_kd.end,k_u,K_huts,Rel_loc.end,end.time,p.prior,q.prior,dose,first.measurement,movement.tied.to.exit) }
  if(adaptive == TRUE) {
    sd = diag(sd,nrow = length(startvalue), ncol = length(startvalue))
  }
  for (ii in 1:iterations){
    if (debugging.plot == TRUE) { plot.count = plot.count + 1/iterations }
    kk = kk + 1
    if(adaptive == TRUE) { kk = ii }
    if (ii %% writing.freq == 0){ 
      if (adaptive == FALSE) {
        write.table(cbind(chain[(kk-writing.freq+1):kk,], LLs[(ii-writing.freq+1):ii], chain.prop[(kk-writing.freq+1):kk,]), file = filename, row.names = FALSE, 
                    append = TRUE, col.names = FALSE, sep = ", ")
        temp = array(dim = dim(chain)); temp.prop = array(dim = dim(chain.prop));   
        temp[1,] = chain[writing.freq,]; temp.prop[1,] = chain.prop[writing.freq,] 
        rm(chain); rm(chain.prop)
        chain = temp; chain.prop = temp.prop
        kk = 1
      } else if (adaptive == TRUE){
        write.table(cbind(chain[(ii-writing.freq+1):ii,],LLs[(ii-writing.freq+1):ii], chain.prop[(ii-writing.freq+1):ii,]), file = filename, row.names = FALSE,
                    append = TRUE, col.names = FALSE, sep = ", ")
      }
    }
    former = chain[kk,]
    if(adaptive == TRUE){
      if(ii > adapt.par[1] && ii %% adapt.par[2] == 0 && ii < (adapt.par[4]*iterations) )          # adapt the proposal covariance structure
      {   
        # browser()
        len<-floor(ii*adapt.par[3]):ii
        x<-chain[len,]
        N<-length(len)
        p_sigma <- (N-1) * var(x)/N
        p_sigma <-makePositiveDefinite(p_sigma)   # To deal with rounding problems that can de-symmetrize
        print(p_sigma)
        if(!(0 %in% p_sigma) ) {
          sd<-p_sigma 
          print('updating sd matrix')
        }
        print('not updating sd matrix. p_sigma returns zeros')
      }
    }   
    proposal = proposalf(chain[kk,],SD=sd,relative.sds,adaptive,movement.tied.to.exit)
    if (one.at.a.time == TRUE) { 
      temp = rmultinom(1,1,draw.probabilities)
      proposal[ - which(temp==1)] = former[ - which(temp==1)]
      }
    if (Qfix==TRUE){
      proposal[c(1,2,3)] = proposal[1]  # NEED TO DO SOMETHING ABOUT THIS ONE. IF I GO BY PARAMETER, IT SHOULDN'T ALWAYS BE PROPOSAL[1]
    }
    if (LTFU == FALSE) { proposal.ll = posterior.hut.movement(proposal,t,interval,K_exit,K_kd,k_huts,Rel_loc,delta,kd.interval.factor,p.prior,q.prior,dose,first.measurement,movement.tied.to.exit) }
    if (LTFU == TRUE) {  proposal.ll = posterior.hut.movement(proposal,t,interval,K_exit,K_kd,k_huts, Rel_loc,delta,kd.interval.factor, 
                                                             K_exit.end,K_kd.end,k_u,K_huts,Rel_loc.end,end.time,p.prior,q.prior,dose,first.measurement,movement.tied.to.exit) }
#     print(paste('former ll: ', former.ll))
#     print(paste('proposal ll: ', proposal.ll))
    if (is.numeric(proposal.ll) && is.nan(proposal.ll)==FALSE) {
      ratio.posteriors = exp( proposal.ll - former.ll)
      
        proposal.probabilities = proposal.prob(proposal = proposal, former = former, SD = sd, relative.sds, adaptive,movement.tied.to.exit)
        proposal.probabilities.reciprocal = proposal.prob(proposal = former, former = proposal, SD = sd, relative.sds, adaptive,movement.tied.to.exit)
      if(adaptive == FALSE){  
      acceptance = ratio.posteriors * ( prod(proposal.probabilities.reciprocal[param.vector]) / prod(proposal.probabilities[param.vector]) )
      }
      else {
        acceptance = ratio.posteriors * ( proposal.probabilities.reciprocal[param.vector] / proposal.probabilities[param.vector] )
      }
      if (runif(1) < acceptance && acceptance > 0 && !is.nan(acceptance) ){  # CHECK: why would this ever be smaller than 0?!!!
        chain[kk+1,] = proposal
        former.ll = proposal.ll
        if (debugging.plot == TRUE) { points(plot.count,proposal[4],col='red') }
      }
      else {
        chain[kk+1,] = former
        former.ll = former.ll
        if (debugging.plot == TRUE) { points(plot.count,proposal[4],col='black') }
      }
    }
    else {
      chain[kk+1,] = former
      former.ll = former.ll
      if (debugging.plot == TRUE) { points(plot.count,proposal[4],col='black') }
    }
    LLs = c(LLs, proposal.ll)
    chain.prop[kk+1,] = proposal
  }
  return(chain)
}

posterior.hut.movement <- function(params,t,interval,K_exit,K_kd,k_huts,
                                   Rel_loc,delta,kd.interval.factor = 1,
                                   K_exit.end = NULL,K_kd.end = NULL,k_u = NULL,K_huts = NULL,Rel_loc.end = NULL,end.time = NULL,
                                   p.prior = 'Flat', q.prior = 'Flat', dose = 0,first.measurement = 30, movement.tied.to.exit = FALSE) {  
#   browser()
  if (is.null(k_u)) {
  return(likelihood.hut.movement(params,t,interval,K_exit,K_kd,k_huts,Rel_loc,delta,kd.interval.factor,first.measurement,movement.tied.to.exit) 
         + prior.hut.movement(params,p.prior, q.prior,dose,movement.tied.to.exit))
  }
  else {
    return(likelihood.hut.movement.ltfu(params,t,interval,K_exit,K_kd,k_huts,K_exit.end,K_kd.end,k_u,K_huts,Rel_loc,Rel_loc.end,end.time,delta,
                                        kd.interval.factor,first.measurement,movement.tied.to.exit) 
           + prior.hut.movement(params, p.prior, q.prior, dose, movement.tied.to.exit))
  }
}

prior.hut.movement <- function(params, p.prior = 'Flat', q.prior = 'Flat', dose = 0, movement.tied.to.exit = FALSE){
  # browser()
  q_AorE <- params[1];   q_BorD <- params[2];   q_C <- params[3]
  p_BorD <- params[4];   p_C <- params[5]
  r_AorE <- params[6];   r_BorD <- params[7];   r_C <- params[8]
  k_AorE <- params[9];   k_BorD <- params[10];  k_C <- params[11]
  if(length(params) == 12) {
    u <- params[12];   
  }
  small.step = 1e-8
  
  if (q.prior == 'Flat'){
    q_AorE.prior <- max ( -6000,log( punif(q_AorE + small.step, (1/1200), 1,log=F) - punif(q_AorE - small.step, (1/1200), 1,log=F) ) );     
    q_BorD.prior <- max ( -6000,log( punif(q_BorD + small.step, (1/1200), 1,log=F) - punif(q_BorD - small.step, (1/1200), 1,log=F) ) );     
    q_C.prior    <- max ( -6000,log( punif(q_C + small.step, (1/1200), 1,log=F)    - punif(q_C - small.step, (1/1200), 1,log=F) ) );     
  }
  else if (q.prior == 'Gamma'){
      q_AorE.prior <- max ( -6000,log( pgamma(q_AorE + small.step, shape = 1.5,rate = (1.5/2e-2), log=F) - pgamma(q_AorE - small.step, shape = 1.5, rate = (1.5/2e-2),log=F) ) );   
      q_BorD.prior <- max ( -6000,log( pgamma(q_BorD + small.step, shape = 1.5,rate = (1.5/2e-2), log=F) - pgamma(q_BorD - small.step, shape = 1.5, rate = (1.5/2e-2),log=F) ) );     
      q_C.prior    <- max ( -6000,log( pgamma(q_C + small.step, shape = 1.5,rate = (1.5/2e-2), log=F)    - pgamma(q_C - small.step, shape = 1.5, rate = (1.5/2e-2), log=F) ) ); 
  }
  if (p.prior == 'Flat') {
    p_BorD.prior <- max ( -6000,log( punif(p_BorD + small.step, 0, 1,log=F)      - punif(p_BorD - small.step, 0, 1,log=F) ) );    
  }
  else if (p.prior == 'Beta' ) {
    p_BorD.prior <- max ( -6000,log( pbeta(p_BorD + small.step, 4, 4,log=F)      - pbeta(p_BorD - small.step, 4, 4,log=F) ) );    
  }
  else if (p.prior == 'Dose.dependent')   {
    if (dose == 0) { 
      print('using dose = 0 beta')
      p_BorD.prior <- max ( -6000,log( pbeta(p_BorD + small.step, 18, 18,log=F)    - pbeta(p_BorD - small.step, 18, 18,log=F) ) );    
    }
    else {
      print('using dose > 0 beta')
      p_BorD.prior <- max ( -6000,log( pbeta(p_BorD + small.step, 15, 10,log=F)    - pbeta(p_BorD - small.step, 15, 10,log=F) ) );    
    }
  }
  p_C.prior    <- max ( -6000,log( punif(p_C + small.step, 0, 1,log=F)           - punif(p_C - small.step, 0, 1,log=F) ) );
  if (movement.tied.to.exit == FALSE) {
    r_AorE.prior <- max ( -6000,log( punif(r_AorE + small.step, (1/1200), (1/5),log=F) - punif(r_AorE - small.step, (1/1200), (1/5),log=F) ) );     
    r_BorD.prior <- max ( -6000,log( punif(r_BorD + small.step, (1/1200), (1/5),log=F) - punif(r_BorD - small.step, (1/1200), (1/5),log=F) ) );     
    r_C.prior    <- max ( -6000,log( punif(r_C + small.step, (1/1200), (1/5),log=F)    - punif(r_C - small.step, (1/1200), (1/5),log=F) ) );    
  }
  else {
    print('movement is tied to exit rates')
    # r_AorE.prior <- max ( -6000,log( punif(r_AorE + small.step, 0, .25,log=F)           - punif(r_AorE - small.step, 0, .25,log=F) ) );
    # r_BorD.prior <- max ( -6000,log( punif(r_BorD + small.step, 0, .25,log=F)           - punif(r_BorD - small.step, 0, .25,log=F) ) );
    # r_C.prior    <- max ( -6000,log( punif(r_C + small.step, 0, .25,log=F)           - punif(r_C - small.step, 0, .25,log=F) ) );
    r_AorE.prior <- max ( -6000,log( pbeta(r_AorE + small.step, 1.25, 3.75,log=F)    - pbeta(r_AorE - small.step,  1.25, 3.75,log=F) ) );    
    r_BorD.prior <- max ( -6000,log( pbeta(r_BorD + small.step, 1.25, 3.75,log=F)    - pbeta(r_BorD - small.step,  1.25, 3.75,log=F) ) );    
    r_C.prior    <- max ( -6000,log( pbeta(r_C + small.step, 1.25, 3.75,log=F)    - pbeta(r_C - small.step,  1.25, 3.75,log=F) ) );    
  }
  k_AorE.prior <- max ( -6000,log( punif(k_AorE + small.step, (1/28800), (1/60),log=F) - punif(k_AorE - small.step, (1/28800), (1/60),log=F) ) );     
  k_BorD.prior <- max ( -6000,log( punif(k_BorD + small.step, (1/28800), (1/60),log=F) - punif(k_BorD - small.step, (1/28800), (1/60),log=F) ) );      
  k_C.prior    <- max ( -6000,log( punif(k_C + small.step, (1/28800), (1/60),log=F) - punif(k_C - small.step,(1/28800), (1/60),log=F) ) );    
  
  if (length(params) == 12) { u.prior <- max ( -6000,log( punif(u + small.step, (1/144000), (1/30),log=F) - punif(u - small.step, (1/144000), (1/30),log=F) ) ) }
  if(length(params) == 11) { return(q_AorE.prior + q_BorD.prior + q_C.prior + p_BorD.prior + p_C.prior + r_AorE.prior + r_BorD.prior + r_C.prior + k_AorE.prior + k_BorD.prior + k_C.prior) }
  else { 
    return(q_AorE.prior + q_BorD.prior + q_C.prior + p_BorD.prior + p_C.prior + r_AorE.prior + r_BorD.prior + r_C.prior + k_AorE.prior + k_BorD.prior + k_C.prior + u.prior) 
  }
}

proposalf <- function(params, SD = 1,relative.sds = FALSE, adaptive = FALSE, movement.tied.to.exit = FALSE) {
#   browser()
  q_AorE <- params[1];   q_BorD <- params[2];   q_C <- params[3]
  p_BorD <- params[4];   p_C <- params[5]
  r_AorE <- params[6];   r_BorD <- params[7];   r_C <- params[8]
  k_AorE <- params[9];   k_BorD <- params[10];  k_C <- params[11]
  if (length(params) == 12 ) { u <- params[12] }
  if (length(SD)==1) {SD = rep(SD,length(params))}
  if (relative.sds == TRUE){
    SD[1] = q_AorE * SD[1];   SD[2] = q_BorD * SD[2]; SD[3] = q_C * SD[3];
    SD[4] = p_BorD * SD[4];   SD[5] = p_C * SD[5];
    SD[6] = r_AorE * SD[6];   SD[7] = r_BorD * SD[7]; SD[8 ] = r_C * SD[8];
    SD[9] = k_AorE * SD[9];   SD[10] = k_BorD * SD[10]; SD[11] = k_C * SD[11];
    if (length( params) == 12) { SD[12] = u *SD[12] }
  }
  if (adaptive == FALSE) {
    repeat{
      q_AorE.proposal <- q_AorE + rnorm(1,sd=SD[1]);   q_BorD.proposal <- q_BorD + rnorm(1,sd=SD[2]);   q_C.proposal <- q_C + rnorm(1,sd=SD[3]);
      p_BorD.proposal <- p_BorD + rnorm(1,sd=SD[4]);      p_C.proposal <- p_C + rnorm(1,sd=SD[5]); 
      r_AorE.proposal <- r_AorE + rnorm(1,sd=SD[6]);   r_BorD.proposal <- r_BorD + rnorm(1,sd=SD[7]);   r_C.proposal <- r_C + rnorm(1,sd=SD[8]);
      k_AorE.proposal <- k_AorE + rnorm(1,sd=SD[9]);   k_BorD.proposal <- k_BorD + rnorm(1,sd=SD[10]);   k_C.proposal <- k_C + rnorm(1,sd=SD[11]);
      if (length( params) == 12) { u.proposal <- u + rnorm(1,sd=SD[12]) }
      if (length(params) == 11){
        if (movement.tied.to.exit == FALSE) {
          if (all(c(q_AorE.proposal, q_BorD.proposal, q_C.proposal,   
                    p_BorD.proposal, p_C.proposal, 
                    r_AorE.proposal, r_BorD.proposal, r_C.proposal, 
                    k_AorE.proposal, k_BorD.proposal, k_C.proposal) > 0) & all(c(p_BorD.proposal, p_C.proposal)<1)) {
            break
          }
        }
        else {
          if (all(c(q_AorE.proposal, q_BorD.proposal, q_C.proposal,   
                    p_BorD.proposal, p_C.proposal, 
                    r_AorE.proposal, r_BorD.proposal, r_C.proposal, 
                    k_AorE.proposal, k_BorD.proposal, k_C.proposal) > 0) & all(c(p_BorD.proposal, p_C.proposal,r_AorE.proposal, r_BorD.proposal, r_C.proposal)<1)) {
            break
          }
        }
      }
      if (length(params) == 12){
        if (movement.tied.to.exit == FALSE) {
          if (all(c(q_AorE.proposal, q_BorD.proposal, q_C.proposal,   
                    p_BorD.proposal, p_C.proposal, 
                    r_AorE.proposal, r_BorD.proposal, r_C.proposal, 
                    k_AorE.proposal, k_BorD.proposal, k_C.proposal, u.proposal) > 0) & all(c(p_BorD.proposal, p_C.proposal)<1)) {
            break
          }
        }
        else {
          if (all(c(q_AorE.proposal, q_BorD.proposal, q_C.proposal,   
                    p_BorD.proposal, p_C.proposal, 
                    r_AorE.proposal, r_BorD.proposal, r_C.proposal, 
                    k_AorE.proposal, k_BorD.proposal, k_C.proposal, u.proposal) > 0) & all(c(p_BorD.proposal, p_C.proposal, r_AorE.proposal, r_BorD.proposal, r_C.proposal)<1)) {
            break
          }
        }
      }
    }
  }
  if (adaptive == TRUE){
    repeat{
      proposals<-mvrnorm(1,mu=params,Sigma=SD)  # Draw proposal point
      q_AorE.proposal <- proposals[1];   q_BorD.proposal <- proposals[2];   q_C.proposal <- proposals[3]
      p_BorD.proposal <- proposals[4];   p_C.proposal <- proposals[5]
      r_AorE.proposal <- proposals[6];   r_BorD.proposal <- proposals[7];   r_C.proposal <- proposals[8]
      k_AorE.proposal <- proposals[9];   k_BorD.proposal <- proposals[10];  k_C.proposal <- proposals[11]
      if (length( params) == 12) { u.proposal <-proposals[12] }
      if (length( params) == 11) { 
        if (movement.tied.to.exit == FALSE) {
          if (all(c(q_AorE.proposal, q_BorD.proposal, q_C.proposal,   
                    p_BorD.proposal, p_C.proposal, 
                    r_AorE.proposal, r_BorD.proposal, r_C.proposal, 
                    k_AorE.proposal, k_BorD.proposal, k_C.proposal) > 0) & all(c(p_BorD.proposal, p_C.proposal)<1)) {
            break
          }
        }
        else {
          if (all(c(q_AorE.proposal, q_BorD.proposal, q_C.proposal,   
                    p_BorD.proposal, p_C.proposal, 
                    r_AorE.proposal, r_BorD.proposal, r_C.proposal, 
                    k_AorE.proposal, k_BorD.proposal, k_C.proposal) > 0) & all(c(p_BorD.proposal, p_C.proposal,r_AorE.proposal, r_BorD.proposal, r_C.proposal)<1)) {
            break
          }
        }
      }
      if (length( params) == 12) { 
        if (movement.tied.to.exit == FALSE) {
          if (all(c(q_AorE.proposal, q_BorD.proposal, q_C.proposal,   
                    p_BorD.proposal, p_C.proposal, 
                    r_AorE.proposal, r_BorD.proposal, r_C.proposal, 
                    k_AorE.proposal, k_BorD.proposal, k_C.proposal, u.proposal) > 0) & all(c(p_BorD.proposal, p_C.proposal)<1)) {
            break
          }
        }
        else {
          if (all(c(q_AorE.proposal, q_BorD.proposal, q_C.proposal,   
                    p_BorD.proposal, p_C.proposal, 
                    r_AorE.proposal, r_BorD.proposal, r_C.proposal, 
                    k_AorE.proposal, k_BorD.proposal, k_C.proposal, u.proposal) > 0) & all(c(p_BorD.proposal, p_C.proposal, r_AorE.proposal, r_BorD.proposal, r_C.proposal)<1)) {
            break
          }
        }
      }
    }
  }
  if (length( params) == 12) { return(c(q_AorE.proposal, q_BorD.proposal, q_C.proposal, p_BorD.proposal, p_C.proposal, r_AorE.proposal, r_BorD.proposal, r_C.proposal, k_AorE.proposal, k_BorD.proposal, k_C.proposal,
                                        u.proposal)) }
  else { 
    return(c(q_AorE.proposal, q_BorD.proposal, q_C.proposal, p_BorD.proposal, p_C.proposal, r_AorE.proposal, r_BorD.proposal, r_C.proposal, k_AorE.proposal, k_BorD.proposal, k_C.proposal)) 
  }
}

proposal.prob <- function(proposal, former, SD = sd,relative.sds = FALSE, adaptive = FALSE, movement.tied.to.exit = FALSE) {
#   browser()
  q_AorE <- proposal[1];   q_BorD <- proposal[2];   q_C <- proposal[3]
  p_BorD <- proposal[4];   p_C <- proposal[5]
  r_AorE <- proposal[6];   r_BorD <- proposal[7];   r_C <- proposal[8]
  k_AorE <- proposal[9];   k_BorD <- proposal[10];  k_C <- proposal[11]
  if (length( proposal) == 12) { u <- proposal[12] }
  q_AorE.former <- former[1];   q_BorD.former <- former[2];   q_C.former <- former[3]
  p_BorD.former <- former[4];   p_C.former <- former[5]
  r_AorE.former <- former[6];   r_BorD.former <- former[7];   r_C.former <- former[8]
  k_AorE.former <- former[9];   k_BorD.former <- former[10];  k_C.former <- former[11]
  if (length( proposal) == 12) { u.former <- former[12] }
  small.step = 1e-8
  if (adaptive == FALSE){
    if (length(SD)==1) { SD = rep(SD,length(proposal)) }
    if (relative.sds == TRUE){
      SD[1] = q_AorE.former * SD[1];   SD[2] = q_BorD.former * SD[2]; SD[3] = q_C.former * SD[3];
      SD[4] = p_BorD.former * SD[4];   SD[5] = p_C.former * SD[5];
      SD[6] = r_AorE.former * SD[6];   SD[7] = r_BorD.former * SD[7]; SD[8 ] = r_C.former * SD[8];
      SD[9] = k_AorE.former * SD[9];   SD[10] = k_BorD.former * SD[10]; SD[11] = k_C.former * SD[11];
      if (length( proposal) == 12) { SD[12] = u.former * SD[12] }
    }
    q_AorE.prob <- ( pnorm(q_AorE+small.step,q_AorE.former,sd=SD[1]) - pnorm(q_AorE-small.step,q_AorE.former,SD[1]) ) / ( 1 - pnorm(0,q_AorE.former,SD[1]) )  ;   
    q_BorD.prob <- ( pnorm(q_BorD+small.step,q_BorD.former,sd=SD[2]) - pnorm(q_BorD-small.step,q_BorD.former,sd=SD[2]) ) / ( 1 - pnorm(0,q_BorD.former,SD[2]) ) ;     
    q_C.prob <- ( pnorm(q_C+small.step,q_C.former,sd=SD[3]) - pnorm(q_C-small.step,q_C.former,sd=SD[3]) ) / ( 1 - pnorm(0,q_C.former,SD[3]) ) ;  
    p_BorD.prob <- ( pnorm(p_BorD+small.step,p_BorD.former,sd=SD[4]) - pnorm(p_BorD-small.step,p_BorD.former,sd=SD[4]) ) / ( pnorm(1,p_BorD.former,SD[4]) - pnorm(0,p_BorD.former,SD[4]) )  ;     
    p_C.prob <- ( pnorm(p_C+small.step,p_C.former,sd=SD[5]) - pnorm(p_C-small.step,p_C.former,sd=SD[5]) ) / ( pnorm(1,p_C.former,SD[5]) - pnorm(0,p_C.former,SD[5]) )  ;     
    if (movement.tied.to.exit == FALSE) {
      r_AorE.prob <- ( pnorm(r_AorE+small.step,r_AorE.former,sd=SD[6]) - pnorm(r_AorE-small.step,r_AorE.former,sd=SD[6]) ) / ( 1 - pnorm(0,r_AorE.former,SD[6]) )  ;   
      r_BorD.prob <- ( pnorm(r_BorD+small.step,r_BorD.former,sd=SD[7]) - pnorm(r_BorD-small.step,r_BorD.former,sd=SD[7]) ) / ( 1 - pnorm(0,r_BorD.former,SD[7]) )  ;   
      r_C.prob <- ( pnorm(r_C+small.step,r_C.former,sd=SD[8]) - pnorm(r_C-small.step,r_C.former,sd=SD[8]) ) / ( 1 - pnorm(0,r_C.former,SD[8]) )  ;  
    }
    else {
      r_AorE.prob <- ( pnorm(r_AorE+small.step,r_AorE.former,sd=SD[6]) - pnorm(r_AorE-small.step,r_AorE.former,sd=SD[6]) ) / ( pnorm(1,r_AorE.former,SD[6]) - pnorm(0,r_AorE.former,SD[6]) )  ;   
      r_BorD.prob <- ( pnorm(r_BorD+small.step,r_BorD.former,sd=SD[6]) - pnorm(r_BorD-small.step,r_BorD.former,sd=SD[6]) ) / ( pnorm(1,r_BorD.former,SD[6]) - pnorm(0,r_BorD.former,SD[6]) )  ;   
      r_C.prob <- ( pnorm(r_C+small.step,r_C.former,sd=SD[6]) - pnorm(r_C-small.step,r_C.former,sd=SD[6]) ) / ( pnorm(1,r_C.former,SD[6]) - pnorm(0,r_C.former,SD[6]) )  ;   
    }
    k_AorE.prob <- ( pnorm(k_AorE+small.step,k_AorE.former,sd=SD[9]) - pnorm(k_AorE-small.step,k_AorE.former,sd=SD[9]) ) / ( 1 - pnorm(0,k_AorE.former,SD[9]) )  ;      
    k_BorD.prob <- ( pnorm(k_BorD+small.step,k_BorD.former,sd=SD[10]) - pnorm(k_BorD-small.step,k_BorD.former,sd=SD[10]) ) / ( 1 - pnorm(0,k_BorD.former,SD[10]) )  ;      
    k_C.prob <- ( pnorm(k_C+small.step,k_C.former,sd=SD[11]) - pnorm(k_C-small.step,k_C.former,sd=SD[11]) ) / ( 1 - pnorm(0,k_C.former,SD[11]) )  ;  
    if (length( proposal) == 12) { 
      u.prob <- ( pnorm(u+small.step,u.former,sd=SD[12]) - pnorm(u-small.step,u.former,sd=SD[12]) ) / ( 1 - pnorm(0,u.former,SD[12]) ) 
      return(c(q_AorE.prob, q_BorD.prob, q_C.prob, p_BorD.prob, p_C.prob, r_AorE.prob, r_BorD.prob, r_C.prob, k_AorE.prob, k_BorD.prob, k_C.prob, u.prob))
    }
    else { return(c(q_AorE.prob, q_BorD.prob, q_C.prob, p_BorD.prob, p_C.prob, r_AorE.prob, r_BorD.prob, r_C.prob, k_AorE.prob, k_BorD.prob, k_C.prob)) }
  }  
  if(adaptive == TRUE){
    if (length(SD)==1) { SD = diag(SD,length(proposal),length(proposal)) }
    if (length( proposal) == 12) {
      if (movement.tied.to.exit == FALSE) {
        prob = pmvnorm(lower = proposal-small.step, upper = proposal+small.step, mean = former, sigma = SD) / 
          pmvnorm(lower = rep(0,length(proposal)), upper = c(rep(1/0,3),rep(1,2),rep(1/0,7)), 
                  mean = former, sigma = SD)
      }
      else {
        prob = pmvnorm(lower = proposal-small.step, upper = proposal+small.step, mean = former, sigma = SD) / 
          pmvnorm(lower = rep(0,length(proposal)), upper = c(rep(1/0,3),rep(1,5),rep(1/0,4)), 
                  mean = former, sigma = SD)
      }
    }
    else {
      if (movement.tied.to.exit == FALSE) {
        prob = pmvnorm(lower = proposal-small.step, upper = proposal+small.step, mean = former, sigma = SD) / 
          pmvnorm(lower = rep(0,length(proposal)), upper = c(rep(1/0,3),rep(1,2),rep(1/0,6)), 
                  mean = former, sigma = SD)
      }
      else {
        prob = pmvnorm(lower = proposal-small.step, upper = proposal+small.step, mean = former, sigma = SD) / 
          pmvnorm(lower = rep(0,length(proposal)), upper = c(rep(1/0,3),rep(1,5),rep(1/0,3)), 
                  mean = former, sigma = SD)
      }
    }
    return(prob)
  }
}

likelihood.hut.movement <- function(params,t,interval,K_exit,K_kd,k_huts,Rel_loc,delta,kd.interval.factor = 1,first.measurement = 30, movement.tied.to.exit = FALSE) {
  # browser()
  
  if (movement.tied.to.exit == FALSE){
    q_AorE <- params[1];   q_BorD <- params[2];   q_C <- params[3]
    r_AorE <- params[6];   r_BorD <- params[7];   r_C <- params[8]
  }
  else {
    q_AorE <- params[1] * (1 - params[6]);   q_BorD <- params[2] * (1 - params[7]);   q_C <- params[3] * (1 - params[8])
    r_AorE <- params[1] * (params[6]);       r_BorD <- params[2] * (params[7]);       r_C <- params[3] * (params[8])
  }
  p_BorD <- params[4];   p_C <- params[5]
  k_AorE <- params[9];   k_BorD <- params[10];  k_C <- params[11]
  Qs <- c(q_AorE, q_BorD, q_C);   Ps <- c(p_BorD, p_C);  Rs <- c(r_AorE, r_BorD, r_C); Ks <- c(k_AorE, k_BorD, k_C);
  Init <- diag(5) # create 5 sets of initial conditions, one for each release location
  Init <- cbind(Init,matrix(0,5,10))
  Time <- seq(0,max(t)+100,by=delta)
  someData <- rep(NaN, length(Time)*16);  
  Out <- array(someData,c(length(Time),16))
  
  temp.a = ode(Init[1,],Time,hut.movement.model,parms=c(Qs,Ps,Rs,Ks))
  Out = abind(Out,temp.a,along = 3)
  if (is.nan(range(Out[,-1,2])[1])) return( -6000);
  temp.b = ode(Init[2,],Time,hut.movement.model,parms=c(Qs,Ps,Rs,Ks))
  Out = abind(Out,temp.b,along = 3)
  Out = abind(Out,ode(Init[3,],Time,hut.movement.model,parms=c(Qs,Ps,Rs,Ks)),along = 3)
  temp.d = temp.b[,c(1,6,5,4,3,2,11,10,9,8,7,16,15,14,13,12)]
  Out = abind(Out,temp.d,along = 3)
  temp.e = temp.a[,c(1,6,5,4,3,2,11,10,9,8,7,16,15,14,13,12)]
  Out = abind(Out,temp.e,along = 3)
  Out <- Out[,,-1]
  
  if (is.nan(range(Out)[1])) return( -6000 );
  # browser()
  Loc <- as.numeric(Rel_loc)

  Prob_exit <- Prob_kd <- prob_huts <- prob_bb <- numeric()
  for(ii in 1:length(Loc))
  {
    Prob_exit <- rbind(Prob_exit,
                       (Out[t[ii]/delta+1,7:11,Loc[ii]]-Out[(t[ii]/delta)-(interval/delta)+1,7:11,Loc[ii]]) )
    if((t[ii]-first.measurement) %% (interval*kd.interval.factor) == 0 ) {
      if (t[ii]==first.measurement){
        Prob_kd <- rbind(Prob_kd,
                         (Out[(t[ii]/delta+1),12:16,Loc[ii]]-Out[(t[ii]/delta)-((interval)/delta)+1,12:16,Loc[ii]]) )
      }
      else {
        Prob_kd <- rbind(Prob_kd,
                         (Out[(t[ii]/delta+1),12:16,Loc[ii]]-Out[(t[ii]/delta)-((interval*kd.interval.factor)/delta)+1,12:16,Loc[ii]]) )
        if(range(Prob_kd,na.rm=TRUE)[1]<0) {browser()}
      }
      prob_huts <- rbind(prob_huts,
                         sum(Out[(t[ii]/delta+1),2:6,Loc[ii]]) )
    }
    else {
      Prob_kd <- rbind(Prob_kd,NaN) # NaNs are to ensure an error if we're using the wrong prob at a specific interval
      prob_huts <- rbind(prob_huts,NaN)
    }
  }
  Ind <- which((t-first.measurement) %% (interval*kd.interval.factor) == 0 )
  Ind2 = which(t==max(t))
  prob_huts <- as.vector(prob_huts); 
  
  LL <- 0;
  temp <- dmnom(c(K_exit[,1],K_exit[,2], K_exit[,3], K_exit[,4], K_exit[,5], 
                  K_kd[Ind,1],K_kd[Ind,2], K_kd[Ind,3], K_kd[Ind,4], K_kd[Ind,5], k_huts[Ind2]),
                prob=c(Prob_exit[,1],Prob_exit[,2], Prob_exit[,3], Prob_exit[,4], Prob_exit[,5], 
                       Prob_kd[Ind,1], Prob_kd[Ind,2], Prob_kd[Ind,3], Prob_kd[Ind,4], Prob_kd[Ind,5],prob_huts[Ind2]),
                log=TRUE)
  LL <- LL+temp 
  
  return(LL) 
}

likelihood.hut.movement.ltfu <- function(params,t,interval,K_exit,K_kd,k_huts,K_exit.end,K_kd.end,k_u,K_huts,Rel_loc,Rel_loc.end,end.time,delta,
                                         kd.interval.factor = 1, first.measurement = 30, movement.tied.to.exit){
#   browser()
  if (movement.tied.to.exit == FALSE){
    q_AorE <- params[1];   q_BorD <- params[2];   q_C <- params[3]
    r_AorE <- params[6];   r_BorD <- params[7];   r_C <- params[8]
  }
  else {
    q_AorE <- params[1] * (1 - params[6]);   q_BorD <- params[2] * (1 - params[7]);   q_C <- params[3] * (1 - params[8])
    r_AorE <- (params[6]) * params[1] ;       r_BorD <- (params[7]) * params[2];       r_C <- (params[8]) * params[3]
  }
  p_BorD <- params[4];   p_C <- params[5]
  k_AorE <- params[9];   k_BorD <- params[10];  k_C <- params[11]
  u <- params[12];  
  Qs <- c(q_AorE, q_BorD, q_C);   Ps <- c(p_BorD, p_C);  Rs <- c(r_AorE, r_BorD, r_C); 
  Ks <- c(k_AorE, k_BorD, k_C);   
  Init <- diag(5) # create 5 sets of initial conditions, one for each release location
  Init <- cbind(Init,matrix(0,5,15))
  Time <- seq(0,max(t)+100,by=delta)
  someData <- rep(NaN, length(Time)*21);  
  Out <- array(someData,c(length(Time),21))
  
  temp.a = ode(Init[1,],Time,hut.movement.model.ltfu,parms=c(Qs,Ps,Rs,Ks,u))
  Out = abind(Out,temp.a,along = 3)
  if (is.nan(range(Out[,-1,2])[1])) return( -6000);
  temp.b = ode(Init[2,],Time,hut.movement.model.ltfu,parms=c(Qs,Ps,Rs,Ks,u))
  Out = abind(Out,temp.b,along = 3)
  Out = abind(Out,ode(Init[3,],Time,hut.movement.model.ltfu,parms=c(Qs,Ps,Rs,Ks,u)),along = 3)
  temp.d = temp.b[,c(1,6,5,4,3,2,11,10,9,8,7,16,15,14,13,12,21,20,19,18,17)]
  Out = abind(Out,temp.d,along = 3)
  temp.e = temp.a[,c(1,6,5,4,3,2,11,10,9,8,7,16,15,14,13,12,21,20,19,18,17)]
  Out = abind(Out,temp.e,along = 3)
  Out <- Out[,,-1]
  
  if (is.nan(range(Out)[1])) return( -6000 );
  
  Loc <- match(Rel_loc,LETTERS)
  Prob_exit <- Prob_kd <- numeric()
  for(ii in 1:length(Loc))
  {
    Prob_exit <- rbind(Prob_exit,
                       (Out[t[ii]/delta+1,7:11,Loc[ii]]-Out[(t[ii]/delta)-(interval/delta)+1,7:11,Loc[ii]]) )
    if((t[ii]-first.measurement) %% (interval*kd.interval.factor) == 0 ) {
      if (t[ii]==first.measurement){
        Prob_kd <- rbind(Prob_kd,
                         (Out[t[ii]/delta+1,12:16,Loc[ii]]-Out[(t[ii]/delta)-(interval/delta)+1,12:16,Loc[ii]]) )
      }
      else {
        Prob_kd <- rbind(Prob_kd,
                         (Out[t[ii]/delta+1,12:16,Loc[ii]]-Out[(t[ii]/delta)-((interval*kd.interval.factor)/delta)+1,12:16,Loc[ii]]) )
        if(range(Prob_kd,na.rm=TRUE)[1]<0) {browser()}
      }
    }
    else {
      Prob_kd <- rbind(Prob_kd,NaN) # NaNs are to ensure an error if we're using the wrong prob at a specific interval
    }
  }
  
  Ind <- which((t-first.measurement) %% (interval*kd.interval.factor) == 0 )
  Loc <- match(Rel_loc.end,LETTERS) # now the likelihoods for the end of the experiment
  Prob_exit.end <- Prob_kd.end <- prob_u <- Prob_huts <- numeric();  
  for(ii in 1:length(Loc))
  {
    Prob_exit.end <- rbind(Prob_exit.end,
                           ( Out[end.time/delta+1,7:11,Loc[ii]] - Out[(end.time/delta)-(interval/delta)+1,7:11,Loc[ii]]) )
    Prob_kd.end <- rbind(Prob_kd.end,
                         ( Out[(end.time/delta+1),12:16,Loc[ii]] - Out[(end.time/delta)-((interval*kd.interval.factor)/delta)+1,12:16,Loc[ii]]) )
    Prob_huts <- rbind(Prob_huts,
                       ( Out[(end.time/delta+1),2:6,Loc[ii]] ) )  
    prob_u <- rbind(prob_u,
                    ( sum(Out[(end.time/delta+1),17:21,Loc[ii]]) ) )  
  }
  prob_u <- as.vector(prob_u)
  
  Ind <- which((t-first.measurement) %% (interval*kd.interval.factor) == 0 )
  LL <- 0;
  
  temp <- 0;
  temp <- dmnom(c(K_exit[,1],K_exit[,2], K_exit[,3], K_exit[,4], K_exit[,5], 
                  K_kd[Ind,1],K_kd[Ind,2], K_kd[Ind,3], K_kd[Ind,4], K_kd[Ind,5], 
                  K_exit.end[,1],K_exit.end[,2], K_exit.end[,3], K_exit.end[,4], K_exit.end[,5], 
                  K_kd.end[,1],K_kd.end[,2], K_kd.end[,3], K_kd.end[,4], K_kd.end[,5], 
                  k_u[], 
                  K_huts[,1],K_huts[,2], K_huts[,3], K_huts[,4], K_huts[,5]), 
                prob=c(Prob_exit[,1],Prob_exit[,2], Prob_exit[,3], Prob_exit[,4], Prob_exit[,5], 
                       Prob_kd[Ind,1], Prob_kd[Ind,2], Prob_kd[Ind,3], Prob_kd[Ind,4], Prob_kd[Ind,5],
                       Prob_exit.end[,1],Prob_exit.end[,2], Prob_exit.end[,3], Prob_exit.end[,4], Prob_exit.end[,5], 
                       Prob_kd.end[,1], Prob_kd.end[,2], Prob_kd.end[,3], Prob_kd.end[,4], Prob_kd.end[,5],
                       prob_u[],
                       Prob_huts[,1], Prob_huts[,2], Prob_huts[,3], Prob_huts[,4], Prob_huts[,5]),
                log=TRUE)
  LL <- LL+temp 
  return(LL) 
}
  
sample.from.priors <- function(n, LTFU = FALSE, movement.tied.to.exit = FALSE){
  q_AorE.start <- rgamma(n, 0.5, 4);   q_BorD.start <- rgamma(n, 0.5, 4);   q_C.start <- rgamma(n, 0.5, 4)
  p_BorD.start <- rbeta(n, 1, 1);      p_C.start <- rbeta(n, 4, 4)
  if (movement.tied.to.exit == FALSE){
    r_AorE.start <- rgamma(n, 0.5, 4); r_BorD.start <- rgamma(n, 0.5, 4);   r_C.start <- rgamma(n, 0.5, 4)
  }
  else {
    r_AorE.start <- rbeta(n, 1, 1);    r_BorD.start <- rbeta(n, 1, 1);      r_C.start <- rbeta(n, 1, 1);
  }
  k_AorE.start <- rgamma(n, 0.5, 4);   k_BorD.start <- rgamma(n, 0.5, 4);   k_C.start <- rgamma(n, 0.5, 4)
  if (LTFU == TRUE) { 
    u.start <- rgamma(n, 0.5, 4) 
    return(cbind(q_AorE.start, q_BorD.start, q_C.start, p_BorD.start, p_C.start, r_AorE.start, r_BorD.start, r_C.start, k_AorE.start, k_BorD.start, k_C.start, u.start))
  }
  else { return(cbind(q_AorE.start, q_BorD.start, q_C.start, p_BorD.start, p_C.start, r_AorE.start, r_BorD.start, r_C.start, k_AorE.start, k_BorD.start, k_C.start)) }
}

sample.from.uniform <- function(n, LTFU = FALSE, movement.tied.to.exit = FALSE){
  q_AorE.start <- runif(n, 1/360, 1/30);      q_BorD.start <- runif(n, 1/360, 1/30);     q_C.start <- runif(n, 1/360, 1/30)
  p_BorD.start <- runif(n, 0.5, 1);          p_C.start <- 0.5; #runif(n, 0, 1)
  if (movement.tied.to.exit == FALSE){
    r_AorE.start <- runif(n, 1/900, 1/300);  r_BorD.start <- runif(n,1/900, 1/300);    r_C.start <- runif(n, 1/900, 1/300)
  }
  else {
    r_AorE.start <- runif(n, 0, 0.5);          r_BorD.start <- runif(n, 0, 0.5);           r_C.start <- runif(n, 0, 0.5)
  }
  k_AorE.start <- runif(n, 1/1400, 1/720);   k_BorD.start <- runif(n,1/1400, 1/720);   k_C.start <- runif(n, 1/1400, 1/720)
  if (LTFU == TRUE) { 
    u.start <- runif(n, 1/2000, 1/1000)
    return(cbind(q_AorE.start, q_BorD.start, q_C.start, p_BorD.start, p_C.start, r_AorE.start, r_BorD.start, r_C.start, k_AorE.start, k_BorD.start, k_C.start, u.start))
  }
  else { return(cbind(q_AorE.start, q_BorD.start, q_C.start, p_BorD.start, p_C.start, r_AorE.start, r_BorD.start, r_C.start, k_AorE.start, k_BorD.start, k_C.start)) }
}

sample.sobol <- function(n, LTFU = FALSE, movement.tied.to.exit = FALSE){
  if (LTFU == FALSE) {
    if (movement.tied.to.exit == FALSE) {
    Out = sobolDesign(lower = c(q_AorE.start = log(1/400), q_BorD.start = log(1/400), q_C.start = log(1/400),
                                p_BorD.start = log(0.1), p_C.start = log(0.5), 
                                r_AorE.start = log(1/1200), r_BorD.start = log(1/1200), r_C.start = log(1/1200), 
                                k_AorE.start = log(1/14400), k_BorD.start = log(1/14400), k_C.start = log(1/14400)),
                      upper = c(q_AorE.start = log(1/30), q_BorD.start = log(1/30), q_C.start = log(1/30),
                                p_BorD.start = 0, p_C.start = log(0.5), 
                                r_AorE.start =log(1/100), r_BorD.start = log(1/100), r_C.start = log(1/100), 
                                k_AorE.start = log(1/720), k_BorD.start = log(1/720), k_C.start = log(1/720)), 
                      n); 
    }
    else {
      Out = sobolDesign(lower = c(q_AorE.start = log(1/1200), q_BorD.start = log(1/1200), q_C.start = log(1/1200),
                                  p_BorD.start = log(0.1), p_C.start = log(0.5), 
                                  r_AorE.start = log(0.01), r_BorD.start = log(0.01), r_C.start = log(0.01), 
                                  k_AorE.start = log(1/14400), k_BorD.start = log(1/14400), k_C.start = log(1/14400)),
                        upper = c(q_AorE.start = log(1/100), q_BorD.start = log(1/100), q_C.start = log(1/100),
                                  p_BorD.start = 0, p_C.start = log(0.5), 
                                  r_AorE.start = 0, r_BorD.start = 0, r_C.start = 0, 
                                  k_AorE.start = log(1/720), k_BorD.start = log(1/720), k_C.start = log(1/720)), 
                        n); 
    }
  }
  else {
    if (movement.tied.to.exit == FALSE) {
      Out = sobolDesign(lower = c(q_AorE.start = log(1/1200), q_BorD.start = log(1/1200), q_C.start = log(1/1200),
                                  p_BorD.start = log(0.1), p_C.start = log(0.5), 
                                  r_AorE.start = log(1/1200), r_BorD.start = log(1/1200), r_C.start = log(1/1200), 
                                  k_AorE.start = log(1/144000), k_BorD.start = log(1/144000), k_C.start = log(1/144000), u.start = log(1/144000)),
                        upper = c(q_AorE.start = log(1/150), q_BorD.start = log(1/150), q_C.start = log(1/150),
                                  p_BorD.start = 0, p_C.start = log(0.5), 
                                  r_AorE.start = log(1/500), r_BorD.start = log(1/500), r_C.start = log(1/500), 
                                  k_AorE.start = log(1/720), k_BorD.start = log(1/720), k_C.start = log(1/720), u.start = log(1/400)), 
                        n);
    }
    else {
      Out = sobolDesign(lower = c(q_AorE.start = log(1/1200), q_BorD.start = log(1/1200), q_C.start = log(1/1200),
                                  p_BorD.start = log(0.1), p_C.start = log(0.5), 
                                  r_AorE.start = log(0.01), r_BorD.start = log(0.01), r_C.start = log(0.01), 
                                  k_AorE.start = log(1/144000), k_BorD.start = log(1/144000), k_C.start = log(1/144000), u.start = log(1/144000)),
                        upper = c(q_AorE.start = log(1/150), q_BorD.start = log(1/150), q_C.start = log(1/150),
                                  p_BorD.start = 0, p_C.start = log(0.5), 
                                  r_AorE.start = log(0.5), r_BorD.start = log(0.5), r_C.start = log(0.5), 
                                  k_AorE.start = log(1/6000), k_BorD.start = log(1/6000), k_C.start = log(1/6000), u.start = log(1/1000)), 
                        n);
    }
  }
  return(Out)
}

dmnom <- function(x,prob,log=FALSE) {
  r <- lgamma(sum(x) + 1) + sum(x * log(prob) - lgamma(x + 1))
  if (log) r else exp(r)
}


################################################################################
# Functions from MHadaptive package used when adaptive = TRUE
# FUNCTION:                 DESCRIPTION:
#  isPositiveDefinite        M  Checks if the matrix X is positive definite
#  makePositiveDefinite      M  Forces the matrix x to be positive definite
################################################################################


isPositiveDefinite <-
  function(x)
  {
    # Transform:
    x = as.matrix(x)
    
    # Check if matrix is positive definite:
    ans = .is.positive.definite(m = x)
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.is.positive.definite <-
  function (m, tol, method = c("eigen", "chol"))
  {
    # Author:
    #   Copyright 2003-05 Korbinian Strimmer
    #   Rank, condition, and positive definiteness of a matrix
    #   GNU General Public License, Version 2
    
    method = match.arg(method)
    if (!is.matrix(m)) {
      m = as.matrix(m)
    }
    if (method == "eigen") {
      eval = eigen(m, only.values = TRUE)$values
      if( missing(tol) ) {
        tol = max(dim(m))*max(abs(eval))*.Machine$double.eps
      }
      if (sum(eval > tol) == length(eval)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else if (method == "chol") {
      val = try(chol(m), silent = TRUE)
      if (class(val) == "try-error") {
        return(FALSE)
      } else {
        return(TRUE)
      }
    }
  }


# ------------------------------------------------------------------------------


makePositiveDefinite <-
  function(x)
  {
    # Make Positive Definite:
    ans = .make.positive.definite(m = x)
    
    # Return Value:
    ans
  }


# ------------------------------------------------------------------------------


.make.positive.definite <-
  function(m, tol)
  {
    # Author:
    #   Copyright 2003-05 Korbinian Strimmer
    #   Rank, condition, and positive definiteness of a matrix
    #   GNU General Public License, Version 2
    
    # Method by Higham 1988
    
    if (!is.matrix(m)) {
      m = as.matrix(m)
    }
    
    d = dim(m)[1]
    if ( dim(m)[2] != d ) {
      stop("Input matrix is not square!")
    }
    
    es = eigen(m)
    esv = es$values
    
    if (missing(tol)) {
      tol = d*max(abs(esv))*.Machine$double.eps
    }
    delta =  2*tol
    # factor two is just to make sure the resulting
    # matrix passes all numerical tests of positive definiteness
    
    tau = pmax(0, delta - esv)
    dm = es$vectors %*% diag(tau, d) %*% t(es$vectors)
    
    return( m + dm )
  }



