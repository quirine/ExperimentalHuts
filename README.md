repository to reproduce:

Title: "Model-based analysis of experimental hut data elucidates a multitude of effects of volatile chemicals on Aedes aegypti mosquitoes" 
Authors: Quirine A. ten Bosch, Fanny Castro-Llanos, Hortance Manda, Amy C. Morrison, John Grieco, Nicole L. Achee, T. Alex Perkins 
Journal: Submitted to 'Parasites and Vectors' 
Year: submission 2017

To run a test on simulated data, followed by running the algorithm on real data for baseline, high, and low dosage, run:

Main.Manuscript.SR.R

Overall set up of code base: 

Main.....: are drivers of specific procedures 
SR.....: function files 
./Processingfiles: contains script to process the results from iterations 
with: - Process_SimulatedSet.R to recreate figures 3 and 4
	  - Process_Results_MCMConData.R to recreate figures 5 and S2, S3, and S4
	  - Process_Figure6.R to recreate figure  6
	  - Process_TestingConvergence.R for supplementary figures on traceplots and Gelman-Rubin statistics
	  - Process_Results_MCMConData.R derives outcomes used in manuscript
./Runfiles: contains scripts to run on cluster 
./Runfiles/SimulationSets: contains python script to automate the creation of Run-files for cluster

RData-files:
- ExitData.ltfu.RData: experimental hut data used for analysis:  3 different dosages and 5 experimental days per dosage. 
- workspace.parametersweep...: results to reproduce figures 3 and 4 on simulated data for a scenario with 1000 or 25 mosquitoes per hut per experimental day
- timespent.fixedr.RData: data to reproduce figure 6 on the time spent in the huts by dosage.  

