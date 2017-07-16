

stackplot <- function(Time, Out, ylim=NA, main=NA, colors=NA, xlab=NA, ylab=NA) {
  # stacked line plot
  if (is.na(ylim)) {
    ylim=c(0, 1)
  }
  if (is.na(colors)) {
    colors = c("yellow","blue","purple","green","pink")
  }
  for (jj in 1:5){
    summary = rep(0, nrow(Out))
    recent = summary
    if (jj == 1 ) {
      plot(c(-100), c(-100), xlim=c(min(Time, na.rm=T), max(Time, na.rm=T)), ylim=ylim, tck=0 )
    }
    else plot(c(-100), c(-100), xlim=c(min(Time, na.rm=T), max(Time, na.rm=T)), ylim=ylim, tck=0, yaxt='n' )
    for( ii in 1:5) {
      current = Out[,1+jj,ii]
      summary = summary + current
      polygon( x=c(Time,rev(Time)), c(summary,rev(recent)), col=colors[[ii]])
      recent = summary
      current = Out[,6+jj,ii]
      summary = summary + current
      polygon(x=c(Time,rev(Time)),c(summary,rev(recent)),col=colors[[ii]],density=10)
      recent = summary
      current = Out[,11+jj,ii]
      summary = summary + current
      polygon( x=c(Time,rev(Time)), c(summary,rev(recent)), col=colors[[ii]],density=0)
      recent = summary
    }
  }
}

stackplot.release <- function(Time, Out, ylim=NA, main=NA, colors=NA, xlab=NA, ylab=NA) {
  # stacked line plot
  if (is.na(ylim)) {
    ylim=c(0, 1)
  }
  if (is.na(colors)) {
    colors = c("yellow","blue","firebrick1","green","pink")
  }
  for (ii in 1:5){
    summary = rep(0, nrow(Out))
    recent = summary
    if (ii == 1 ) {
      plot(c(-100), c(-100), xlim=c(min(Time, na.rm=T), max(Time, na.rm=T)), ylim=ylim, tck=0 )
    }
    else plot(c(-100), c(-100), xlim=c(min(Time, na.rm=T), max(Time, na.rm=T)), ylim=ylim, tck=0, yaxt='n' )
    for( jj in 1:5) {
      current = Out[,1+jj,ii]
      summary = summary + current
      polygon( x=c(Time,rev(Time)), c(summary,rev(recent)), col=colors[[jj]])
      recent = summary
      current = Out[,6+jj,ii]
      summary = summary + current
      polygon(x=c(Time,rev(Time)),c(summary,rev(recent)),col=colors[[jj]], density=20)
      recent = summary
      current = Out[,11+jj,ii]
      summary = summary + current
      polygon( x=c(Time,rev(Time)), c(summary,rev(recent)), col=colors[[jj]],density=0)
      recent = summary
    }
  }
}


stackplot.norm <- function(Time, Out, ylim=NA, main=NA, colors=NA, xlab=NA, ylab=NA) {
  # stacked line plot
  if (is.na(ylim)) {
    ylim=c(0, 1)
  }
  if (is.na(colors)) {
    colors = c("yellow","blue","purple","green","pink")
  }
  for (jj in 1:5){
    summary = rep(0, nrow(Out))
    recent = summary
    if (jj == 1 ) {
      plot(c(-100), c(-100), xlim=c(min(Time, na.rm=T), max(Time, na.rm=T)), ylim=ylim, tck=0 )
    }
    else plot(c(-100), c(-100), xlim=c(min(Time, na.rm=T), max(Time, na.rm=T)), ylim=ylim, tck=0, yaxt='n' )
      for( ii in 1:5) {
      Total = 0
      for (kk in 1:5) Total = Total + Out[,1+jj,kk] + Out[,6+jj,kk] + Out[,11+jj,kk]
      current = Out[,1+jj,ii]/Total
      summary = summary + current
      polygon( x=c(Time,rev(Time)), c(summary,rev(recent)), col=colors[[ii]])
      recent = summary
      current = Out[,6+jj,ii]/Total
      summary = summary + current
      polygon(x=c(Time,rev(Time)),c(summary,rev(recent)),col=colors[[ii]],density=10)
      recent = summary
      current = Out[,11+jj,ii]/Total
      summary = summary + current
      polygon( x=c(Time,rev(Time)), c(summary,rev(recent)), col=colors[[ii]],density=0)
      recent = summary
    }
  }
}


stackplot.loc <- function(Time, Out, ylim=NA, main=NA, colors=NA, xlab=FALSE, ylab=FALSE) {
  # stacked line plot
  if (is.na(ylim)) {
    ylim=c(0, 1)
  }
  if (is.na(colors)) {
    colors = c("yellow","blue","purple","green","pink")
  }
  for (jj in 1:5){
    summary = rep(0, nrow(Out))
    recent = summary
    if (jj == 1 ) {
    plot(c(-100), c(-100), xlim=c(min(Time, na.rm=T), max(Time, na.rm=T)), ylim=ylim, tck=0 )
    }
    else plot(c(-100), c(-100), xlim=c(min(Time, na.rm=T), max(Time, na.rm=T)), ylim=ylim, tck=0, yaxt='n' )
    if (xlab == FALSE) axis(side = 1,xaxt='n')
    if (jj == 1 & ylab == TRUE) mtext('Proportion',2,line=2.75)
    if (jj == 3 & xlab == TRUE) mtext('Time (min)',1,line=2.5)
    for( ii in 1:5) {
      current = Out[,1+jj,ii]
      summary = summary + current
      polygon( x=c(Time,rev(Time)), c(summary,rev(recent)), col=colors[[ii]])
      recent = summary
    }
  }
}

stackplot.exit <- function(Time, Out, ylim=NA, main=NA, colors=NA, xlab=FALSE, ylab=FALSE) {
  # stacked line plot
  if (is.na(ylim)) {
    ylim=c(0, 1)
  }
  if (is.na(colors)) {
    colors = c("yellow","blue","purple","green","pink")
  }
  for (jj in 1:5){
    summary = rep(0, nrow(Out))
    recent = summary
    if (jj == 1 ) {
      plot(c(-100), c(-100), xlim=c(min(Time, na.rm=T), max(Time, na.rm=T)), ylim=ylim, tck=0 )
    }
    else plot(c(-100), c(-100), xlim=c(min(Time, na.rm=T), max(Time, na.rm=T)), ylim=ylim, tck=0, yaxt='n' )
    if (xlab == FALSE) axis(side = 1,xaxt='n')
    if (jj == 1 & ylab == TRUE) mtext('Proportion',2,line=2.75)
    if (jj == 3 & xlab == TRUE) mtext('Time (min)',1,line=2.5)
    for( ii in 1:5) {
      current = Out[,6+jj,ii]
      summary = summary + current
      polygon(x=c(Time,rev(Time)),c(summary,rev(recent)),col=colors[[ii]],density=10)
      recent = summary
    }
  }
}

stackplot.kd <- function(Time, Out, ylim=NA, main=NA, colors=NA, xlab=FALSE, ylab=FALSE) {
  # stacked line plot
  if (is.na(ylim)) {
    ylim=c(0, 1)
  }
  if (is.na(colors)) {
    colors = c("yellow","blue","purple","green","pink")
  }
  for (jj in 1:5){
    summary = rep(0, nrow(Out))
    recent = summary
    if (jj == 1 ) {
      plot(c(-100), c(-100), xlim=c(min(Time, na.rm=T), max(Time, na.rm=T)), ylim=ylim, tck=0 )
    }
    else plot(c(-100), c(-100), xlim=c(min(Time, na.rm=T), max(Time, na.rm=T)), ylim=ylim, tck=0, yaxt='n' )
    if (xlab == FALSE) axis(side = 1,xaxt='n')
    if (jj == 1 & ylab == TRUE) mtext('Proportion',2,line=2.75)
    if (jj == 3 & xlab == TRUE) mtext('Time (min)',1,line=2.5)
    for( ii in 1:5) {
      current = Out[,11+jj,ii]
      summary = summary + current
      polygon( x=c(Time,rev(Time)), c(summary,rev(recent)), col=colors[[ii]],density=0)
      recent = summary
    }
  }
}

piechart <- function(slices, main=NA, cols=NA, xlab=NA, ylab=NA) {
  if (is.na(cols)) {
    cols = c("yellow","blue","purple")
  }
  pct <- round(slices/sum(slices)*100)
  lbls <- paste(pct,"%",sep="") # ad % to labels 
  pie(slices,labels = lbls, col=cols[1:length(lbls)])
  mtext(main,2,line=2.75)
}

get.prop.timespent <- function(Out){
  Sum.A <- Sum.B <- Sum.C <- Sum.D <- Sum.E <- 0
  for(ii in 1:5) Sum.A <- Sum.A + sum(Out[,2,ii])
  for(ii in 1:5) Sum.B <- Sum.B + sum(Out[,3,ii])
  for(ii in 1:5) Sum.C <- Sum.C + sum(Out[,4,ii])
  for(ii in 1:5) Sum.D <- Sum.D + sum(Out[,5,ii])
  for(ii in 1:5) Sum.E <- Sum.E + sum(Out[,6,ii])
  SUM <- Sum.A + Sum.B + Sum.C + Sum.D + Sum.E
  slices <- c((Sum.A+Sum.E)/SUM,(Sum.B+Sum.D)/SUM,Sum.C/SUM)
    return(slices)
}

get.timespent <- function(Out){
  Sum.A <- Sum.B <- Sum.C <- Sum.D <- Sum.E <- 0
  for(ii in 1:5) Sum.A <- Sum.A + sum(Out[,2,ii])
  for(ii in 1:5) Sum.B <- Sum.B + sum(Out[,3,ii])
  for(ii in 1:5) Sum.C <- Sum.C + sum(Out[,4,ii])
  for(ii in 1:5) Sum.D <- Sum.D + sum(Out[,5,ii])
  for(ii in 1:5) Sum.E <- Sum.E + sum(Out[,6,ii])
  SUM <- Sum.A + Sum.B + Sum.C + Sum.D + Sum.E
  slices <- c((Sum.A+Sum.E),(Sum.B+Sum.D),Sum.C)
  return(slices)
}
  