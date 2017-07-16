rm(list=ls())

# load libraries ----------------------------------------------------------
library(coda)
library(pomp)
library(car)
library(IDPmisc)
library(vioplot)

# set wd ------------------------------------------------------------------
setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# load source files -------------------------------------------------------
source('SR.Exit.MCMC.Functions.R')
source('SR.Exit.FiguresScripts.R')
source('SR.Exit.ProcessOutputs.Functions.R')

# set data directory ------------------------------------------------------
data.path <- "../Output/Simulated/simulated.data.ltfu.RelSD.25_2016-12-16/"
# data.path <- "../Output/Simulated/simulated.data.ltfu.RelSD.1000_2016-12-17/"

# Set directory to save figures -------------------------------------------
figure.path <- '../Figures/Final/'
cd <- getwd()


# set parameters -----------------------------------------------------
burnIn = 10000
num.chains = 5

# specify file names ------------------------------------------------------
setwd(data.path)

list.of.files = dir(pattern = '*.csv')
list.of.datafiles = dir(pattern = '*.RData')

temp = length(list.of.files)
N = temp / num.chains

for (ii in 1:N){
  files = (list.of.files[((ii*5)-4):(ii*5)])
  name <- paste("files.", ii, sep = "")
  assign(name, files)
}

setwd(cd)

# set parameters -----------------------------------------------------
Names=c('movement 2 away','movement 1 away','movement SR', 
        'prop away from SR', 
        'exit 2 away', 'exit 1 away', 'exit SR hut',
        'KD 2 away', 'KD 1 away', 'KD SR hut', 'loss to follow-up' )

Violin.Names = c(expression('Exit (x' [2]*')'),expression('Exit (x' [1]*')'),
                 expression('Exit (x' ['T']*')'),expression('Repellency (p' [1]*')'),
                 expression('KD (k' [2]*')'),expression('KD (k' [1]*')'),
                 expression('KD (k' ['T']*')'),'ltfu (u)')
# load and prepare data ---------------------------------------------------
setwd(data.path)
size = 5e4
for (ii in 1:N){
  for (jj in 1:num.chains){
    files = paste('files.',ii,sep='')
    filename = eval(parse(text=paste(files,'[',jj,']',sep='')))
    chain = read.csv(filename,header=TRUE,colClasses = 'numeric')
    if (dim(chain)[2] > 12) { temp = ((dim(chain)[2] - 1) / 2); chain = chain[-1,1:temp]  }
    else { chain = chain[-1,] }
    chain = mcmc(chain)
    name = paste('chain.',jj,sep='')
    assign(name, chain)
    eval(parse(text=paste(name,'=mcmc(',name,')',sep='')))
  }
  temp = mcmc.list(mcmc(chain.1[burnIn:size,-5]),
                   mcmc(chain.2[burnIn:size,-5]),
                   mcmc(chain.3[burnIn:size,-5]),
                   mcmc(chain.4[burnIn:size,-5]),
                   mcmc(chain.5[burnIn:size,-5])) 
  name = paste('combined.chain',ii,sep='')
  assign(name,temp)
}

setwd(cd)


# Gelman rubin statistics -------------------------------------------------
GEL = numeric()
for (ii in 1:N){
  temp = eval(parse(text=paste('combined.chain',ii,sep='')))
  GEL = c(GEL,gelman.diag(temp,multivariate = FALSE))
}
write.xlsx(GEL[[19]],'test.xls', sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

# Parse outputs for figures -----------------------------------------------
setwd(data.path)

someData = rep(NaN,((size-burnIn+1)*num.chains)*8*N)
Output.all = array(someData, c(((size-burnIn+1)*num.chains),8,N))
Defaults.temp = numeric()
Iis = c(seq(1,4),seq(8,11))

for (jj in 1:length(list.of.datafiles)) {
  load(list.of.datafiles[jj])
  Defaults.temp = rbind(Defaults.temp,defaults[-5])
  for (ii in 1:length(Iis)) {
    iii = Iis[ii]
    temp = eval(parse(text=paste('combined.chain',jj,sep='')))
    if (ii < 4){
      temp.mult = numeric()
      for (kk in 1:num.chains){
        temp.mult = c( temp.mult,temp[[kk]][,iii] * temp[[kk]][,iii + 4] ) 
      }
    }
    else {
      temp.mult = numeric()
      for (kk in 1:num.chains){
        temp.mult = c( temp.mult,temp[[kk]][,iii] ) 
      }
    }
    Output.all[,ii,jj] = temp.mult
  }
}

Defaults = Defaults.temp[,Iis]
for (ii in 1:3){
   Defaults[,ii] = Defaults.temp[,ii]*Defaults.temp[,ii+4]
}

setwd(cd)


# correlations estimated vs median ----------------------------------------
Medians <- Rel.SDs <- matrix(data = NA, dim(Defaults)[1], dim(Defaults)[2])

for (ii in 1:length(list.of.datafiles)){
  for (jj in 1:dim(Defaults)[2]){
    Medians[ii,jj] = median(Output.all[,jj,ii])
    Rel.SDs[ii,jj] = sd(Output.all[,jj,ii])/Defaults[ii,jj]
  }
}

h = boxplot(Rel.SDs)
h$stats[3,]

Correlations = numeric()
for (ii in 1:dim(Defaults)[2]){
  Correlations = c(Correlations, cor(Defaults[,ii], Medians[,ii]) )
}

Correlations
# Violin plots ------------------------------------------------------------
setwd(figure.path)

pdf('Figure4.pdf',width=4.75,height=7.5)

par(mfrow = c(8,1),mai=c(0.3,0.3,0.1,0.1) ,oma=c(3,1,0.1,1),mar=c(0.5,3,1,1.9))  

for (ii in 1:8){
  for (jj in 1:N){
    temp = Output.all[,ii,jj]
    name = paste('x',jj,sep='')
    assign(name,temp)
  }
  par (xaxt='n')
  if (ii == 8) par(xaxt='s')
  vioplot(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10, names=seq(1,10),col = "lightcyan3",rectCol = 'white', colMed = 'black',border = FALSE, pchMed = 16)
  
  points(seq(1,10),Defaults[,ii],pch=23, col = "blue")
  if ( ii == 4) abline(h=0.5,col = 'gray', lwd = 2,lty = 3 )
  mtext(Violin.Names[ii], side = 2, line = 2.25,cex = 0.7)
  if ( ii == 1) legend('topleft',c('simulated true value','estimated median'),pch=c(23,16),col=c('blue','black'),box.lty=0)
}

mtext(side = 1, text = 'Simulated parameter set', line = 1.4, cex = .8, outer = TRUE)  

dev.off()
setwd(cd)

