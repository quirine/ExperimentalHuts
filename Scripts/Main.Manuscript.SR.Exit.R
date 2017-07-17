# Code to reproduce Manuscript:
#Title: "Model-based analysis of experimental hut data elucidates a multitude of effects of volatile chemicals on Aedes aegypti mosquitoes" 
#Authors: Quirine A. ten Bosch, Fanny Castro-Llanos, Hortance Manda, Amy C. Morrison, John Grieco, Nicole L. Achee, T. Alex Perkins 
#Journal: submitted to 'Parasites and Vectors' 
#Year: submission 2017s

rm(list=ls())

# Set work directory ------------------------------------------------------
setwd(dirname(list.files(pattern='Main.Manuscript.SR.Exit.R', recursive=TRUE, full.names=TRUE)))

# Run test run approach with recaptured and lost mosquitoes (ltfu case) ---
source('Parameters.LTFU.SR.Exit.R')
source('Main.SimulateData.LTFU.SR.Exit.R')
source('Main.RunMCMC.LTFU.SR.Exit.R')

# Run on real data -------------------------------------------
rm(list=ls())
source('Parameters.LTFU.SR.Exit.R')
load('ExitData.ltfu.RData')
dose = 0
source('Main.RunMCMConLTFUData.SR.Exit.R')
dose = 0.0625
source('Main.RunMCMConLTFUData.SR.Exit.R')
dose = 0.125
source('Main.RunMCMConLTFUData.SR.Exit.R')






