## sim_2.R contains R commands to execute 'gcipdr application to multi-center IST data' (from omonimous repository). 
## Copyright (C) 2018 Federico Bonofiglio

    ## This Program is free software: you can redistribute it and/or modify
    ## it under the terms of the GNU General Public License as published by
    ## the Free Software Foundation, either version 3 of the License, or
    ## (at your option) any later version.

    ## This Program is distributed in the hope that it will be useful,
    ## but WITHOUT ANY WARRANTY; without even the implied warranty of
    ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    ## GNU General Public License for more details.

    ## You should have received a copy of the GNU General Public License
    ## along with This Program.  If not, see <https://www.gnu.org/licenses/>.


#### simulation study
############ MULTI-CENTER STRUCTURE ....
## Aspirin effect



# R version 3.6.1 (2019-07-05) -- "Action of the Toes"  



library(gcipdr)

library(lme4)

library(meta)
library(metafor)

library(xtable)


source("simulation_functions_v9.R")


## simulations parameters

ncores <- 28 # BEWARE: changing Nr cores can affect random seeds and give different results

 params <- data.frame(n = rep(c(100, 500, 1000), rep(2,3)), K = rep(c(2, 10), 3) )

R <- 100

## rxasp

 system.time( ppp <- TabulateSimulations(R, params, ncores = ncores, latextable = FALSE, beta0_I = -1.1, beta0_II = c(-2, 0), betarsbp_I = -0.3, betarsbp_II = c(-0.1, -0.4), binary_covariate = TRUE ) )  

# in alternative analsis was also set betarsbp_II = c(-0.04, -0.4)

### saveRDS(ppp, "simul_2_final.rds")



### re-arrange table as in publication

rxasp_simul <- matrix( nrow = dim(ppp$central$rxasp)[1], ncol = dim(ppp$central$rxasp)[2] )

rxasp_simul[which(ppp$central$rxasp["params"] != "$\\sigma_{b}$"), ] <- apply(ppp$central$rxasp[ppp$central$rxasp["params"] != "$\\sigma_{b}$", ], 2, as.character)

rxasp_simul[which(ppp$central$rxasp["params"] == "$\\sigma_{b}$"), ] <- apply(ppp$central$`(Intercept)`[ppp$central$`(Intercept)`["params"] == "$\\sigma_{b}$", ], 2, as.character)

 rxaspsimul <- data.frame(rxasp_simul); colnames(rxaspsimul) <- colnames(ppp$central$rxasp)
 
rxaspsimulsm <- rxaspsimul[-c(7:12, 25:30), ]

 levels(rxaspsimulsm$copula.data)[2:3] <- c("\\parbox[t]{2mm}{\\multirow{12}{*}{\\rotatebox[origin=c]{90}{Mismodeled}}}",  
 "\\parbox[t]{2mm}{\\multirow{12}{*}{\\rotatebox[origin=c]{90}{Well modeled}}}")

## small error with multirow em

levels(rxaspsimulsm$K)[2:3] <- c("\\multirow{3}{1em}{10}", "\\multirow{3}{1em}{2}")




