ExperimentalHuts_2017
====================

This repository contains code used in the following paper.

Quirine A. ten Bosch,Fanny Castro-Llanos, Hortance Manda, Amy C. Morrison, John P. Grieco, Nicole L. Achee, T. Alex Perkins (2017) **Model-based analysis of experimental hut data elucidates multifaceted effects of a volatile chemical on Aedes aegypti mosquitoes**. *BioRXivs* 

All code contained within this repository is released under the [CRAPL v0.1 License](http://matt.might.net/articles/crapl/). Data provided in this repository are sufficient to rerun all analyses. Because of the large sizes of the output files, we have deposited those on Open Science Framework [https://osf.io/xtmy7/](https://osf.io/xtmy7/) as part of project [https://osf.io/5hcpf/] (https://osf.io/5hcpf/).


====================

### Set up of code base: 

* 'Scripts' contains code to simulate and fit the model  
* './Runfiles' contains scripts to reproduce all runs described in the manuscript. Calls scripts from the 'Scripts' folder
* './Processingfiles' contains scripts to process the results and create figures

### Scripts folder

The scripts in this folder are used to fit the modeling framework to simulated data as well as the experimental hut data. The outputs of these exercises are given in [https://osf.io/xtmy7/](https://osf.io/xtmy7/) and are processed using the scripts in the Processing folder.  

These scripts for simulations and model fitting were performed on the University of Notre Dame's Center for Research Computing cluster [http://crc.nd.edu](http://crc.nd.edu). Processing of outputs was done on desktop comuputer (Mac OSX) 

To run a test on simulated data and to run the algorithm once on real data for baseline, high, and low dosage, run:

* Main.Manuscript.SR.R, which calls the files below

to simulate data with known parameters and fit the modeling framework to these data: 
* 'Parameters.LTFU.SR.Exit.R'  
* 'Main.SimulateData.LTFU.SR.Exit.R'
* 'Main.RunMCMC.LTFU.SR.Exit.R'

to fit the modeling framework to the real data:
* 'Parameters.LTFU.SR.Exit.R'  
* 'ExitData.ltfu.RData'
* 'Main.RunMCMConLTFUData.SR.Exit.R'

Note that for a quick run, you can adjust the number of iterations in Paramers.LTFU.SR.Exit.R 

Main.....: are drivers of specific procedures 
SR.....: function files 
....RData: data files to reproduce results:
* 'ExitData.ltfu.RData' contains the experimental hut data
* 'timespent.fixedr.RData' contains the results from running the odes on the posterior outcomes to estimate the proportion of time spent in each hut. Can be reproduced with ./Processingfiles/Process_Figure6.R

### Scripts/Runfiles folder

Contains the R-files to run the model framework on the experimental hut data for the three different dosages. 
For each dosage, the R-file is provided that runs at the default r (proportion of movement directed outside), and the low and high values used to assess sensitivity to this value. 

### Scripts/Runfiles/SimulationSets folder 

Contains python script to automate the creation of Run-files to simulate data at different values and run the modeling framework on these data. 


### Scripts/Processingfiles folder 
Contains the scripts to create all the figures in the manuscript. It uses output data deposited on [https://osf.io/xtmy7/](https://osf.io/xtmy7/)

* 'Process_SimulatedSet.R' to recreate figures 3 and 4
* 'Process_Results_MCMConData.R' to recreate figures 5 and S2, S3, and S4
* 'Process_Figure6.R' to recreate figure  6
* 'Process_TestingConvergence.R' for supplementary figures on traceplots and Gelman-Rubin statistics
* 'Process_Results_MCMConData.R' derives outcomes and summary statistics reported in manuscript






* `0_numfunctions.R`
* `1_raster_aligning.R`
* `run2job.pbs` calls `2_fit_attackrates_seroprev.R`
* `run2stat.pbs` calls `2_fit_attackrates_seroprev_stat.R`
* `2_pop_projection.R`
* `3_fill_blank_econs.R`
* `run4job.pbs` calls `4_output_grids_crc.R`
* `run4stat.pbs` calls `4_output_grids_stat_crc.R`
* `run5job.pbs` calls `5_min_max_mean_1st_round.R`
* `5_run_AR_R0.R`
* `run6job.pbs` calls `6_min_max_mean_2nd_round.R`
* `6_fig_maps.R`
* `run7job.pbs` calls `7_min_max_mean_surface.R`
* `run9job.pbs` calls `9_median_1st_round.R`
* `9_median_surface.R`
* `10_country_summary.py`


### data folder

Data included here pertain to the 13 sites listed in Table S1 from which we obtained seroprevalence estimates, 100 replicates of *Aedes aegypti* occurrence probabilities from those sites from [Kraemer et al. (2015)](https://elifesciences.org/content/4/e08347), and a generalized additive model object that describes the relationship between temperature and adult female *Aedes aegypti* mortality from [Brady et al. (2013)](https://parasitesandvectors.biomedcentral.com/articles/10.1186/1756-3305-6-351).


### generated folder

Files included here contain 1,000 replicates of parameterizations of the mechanistic and statistical models, respectively.


### outputs folder

The file included here contains the country-level sums of total infections and infections among childbearing women used in Table 1 and Figures 3 & 4.


### maps folder

This folder contains files that can be used to access raster data shown in the maps in Figures S2-S10. The numbers in each 5x5 km grid cell are only pertinent to that grid cell, and totals across multiple grid cells cannot be interpreted as sums across larger areas. The reason is that, for example, the map of median attack rates contains each grid cell's median across 1,000 replicates rather than a reflection of what might somehow be considered a median spatial layer of attack rates across all grid cells. By contrast, distributions of country- and continent-level totals in Table 1 and Figures 3 & 4 reflect each such quantity calculated in each of the 1,000 replicates and then examined as a distribution.