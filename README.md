This is code to run an application of the 'gcipdr' package (see homonymous repository) on the open multi-center International Stroke Trial (IST) data. This analysis appeared in [Stat Med](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.8470). 

Notes: as by [version 6.1](https://data2knowledge.atlassian.net/wiki/spaces/DSDEV/pages/1600061445/Version+6.1.0) of DataSHIELD there are available functions to compute skewness and kurtosis. 

# Instructions

Simply clone or download the repository. The main analysis is in file 'gcipdr_IST_analysis.R'. Just execute the whole script. 
Additional simulations are in 'sim_1.R' and 'sim_2.R' (simulations on blood-pressure and aspirin effect respectively). Run these scripts separately. 
Files depend on 'simulation_functions_v9.R'. Beware: changing number of cores can yield different results. Computation time in the simulation study can take up to several hours under current settings. 

