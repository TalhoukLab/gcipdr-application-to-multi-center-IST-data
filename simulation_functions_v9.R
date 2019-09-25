## simulation_functions_v9.R contains R commands to execute 'gcipdr application to multi-center IST data' (from homonymous repository). 
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


### SIMULATION MODULES


###

logistic <- function(x) 1 / (1 + exp(-x))

#
 statf <- function( x, formula ){
     
 model <- glm( formula, family = "binomial", data = as.data.frame(x) )

     out <- summary(model)

     if (length(model$coef) > length(out$coef[, 1])){  # summary dropped NAs ?
         nout <- as.matrix(data.frame(model$coef, NA))
     nout[, 2][-which(is.na(nout[, 1]))] <- out$coef[, 2]    
    out <- list(coef = nout, aic = out$aic)
     }
         
     effnames <- rownames(out$coef)
     nas <- rep(NA, length(effnames))
                aics <- rep(out$aic, length(effnames))
                names(aics) <- names(nas) <- effnames
  return(
      list( coef= out$coef[, 1], sd = out$coef[, 2], ranefsd= nas , ranefsdSE = nas, ranefcor = nas, aic = aics )
      )
    }

                                        #

statma <- function(TE, seTE, hubsnames, effect = "OR", statist = c("random", "fixed")){

    statist <- match.arg(statist)

    drop <- which(seTE >= 1e+02)
    if (length(drop) > 0)
        out <- try( metagen(TE[-drop], seTE[-drop], hubsnames[-drop], sm = effect, method.tau = "SJ"), silent = TRUE)
    else
            out <- try( metagen(TE, seTE, hubsnames, sm = effect, method.tau = "SJ"), silent = TRUE)

    if (class(out)[[1]] == "try-error")
       return( data.frame( coef= NA, sd = NA, ranefsd= NA, ranefsdSE = NA, ranefcor = NA, aic = NA ))
    else
    return( data.frame( coef= out[[paste0("TE.",statist)]], sd = out[[paste0("seTE.",statist)]], ranefsd = ifelse(statist[1] == "random", out$tau, NA), ranefsdSE = ifelse(statist[1] == "random", out$se.tau2, NA), ranefcor = NA, aic = NA ) )
       
    }


   
#


statm <- function( x, formula){

  model <-     try(
 glmer(  formula , data = x, family = "binomial", nAGQ = 0,
                         control = glmerControl(optimizer = "Nelder_Mead", boundary.tol = 1e-3) )
                    , silent = TRUE      )
    
   if (class(model) == "try-error")
       return( list( coef= NA, sd = NA, ranefsd= NA, ranefsdSE = NA, ranefcor = NA, aic = NA ))
    
     out <- summary(model)
    effnames <- rownames(out$coef)

    ranefsd <- attr(out$varcor[[1]],"stddev")
    ranefsd  <- c(ranefsd, rep(NA, length(effnames) - length(ranefsd) ) )
    ranefcor <- attr(out$varcor[[1]],"correlation")
      if (sum(dim(ranefcor)) > 2) ### BUG: will not work for length(lower.tri) > 1
          ranefcor <- rep(as.numeric(ranefcor[lower.tri(ranefcor)]),length(effnames) )
    else
        ranefcor <- rep(as.numeric(ranefcor),length(effnames) ) 
    
     nas <- rep(NA, length(effnames))
                aics <- rep(out$AICtab[1], length(effnames))
        names(ranefsd) <-  names(ranefcor) <- names(aics) <- names(nas) <- effnames
        
  return(
      list( coef= out$coef[, 1], sd = out$coef[, 2], ranefsd= ranefsd, ranefsdSE = nas, ranefcor = ranefcor,  aic = aics )
      )
    }


#### SCENARIO (I)

generateDataScenarioOne <- function(n = 100, beta0 = -6, betarsbp = 0.05, betaage = 0.02, K = 2, rdn_sd = c(1.3, 0.0, 0.0), meanrsbp = 150, sdrsbp = 35, meanage= 60, sdage = 10, cor_rsbpage = 0.35, binary_covariate = FALSE ){

    
Sigma <- matrix(c(sdage^2, rep(cor_rsbpage*(sdage*sdrsbp[1]), 2), sdrsbp[1]^2), ncol = 2)
    
simd <- data.frame(MASS::mvrnorm(n*K, c(meanage, meanrsbp[1]), Sigma) )

    if (binary_covariate)
    simd[, 2] <- rbinom(n, 1, 0.5) # random 50/50 treatment assignment
    
colnames(simd) <- c("age", "rsbp")

rnd_int <- rnorm(K, 0, rdn_sd[1])
rnd_slp <- rnorm(K, 0, rdn_sd[2])
rnd_slp2 <- rnorm(K, 0, rdn_sd[3])

    site <- gl(K, n, labels = 1:K)
    
  mus <- beta0[1] + rnd_int[site] + (betarsbp[1] + rnd_slp[site])*simd$rsbp + (betaage[1] + rnd_slp2[site]*simd$age)
    
simd$death <- rbinom(K*n, 1, logistic(mus))

    simd$site <- site 
    
    return(simd)
   
  }


#### SCENARIO (II)

generateDataScenarioTwo <- function(n = 100, beta0 = c(-2, -11), betarsbp = c(0.01, 0.08), betaage = c(0.02, 0.02), K = 2, rdn_sd = c(1.3, 0.0, 0.0), meanrsbp = c(90, 150), sdrsbp = c(10, 35), meanage= 60, sdage = 10, cor_rsbpage = 0.35, agecutoff_rsbp = 40, agecutoff_y = 50, binary_covariate = FALSE){

age <- rgmom(n*K, meanage, sdage)

# blood pressure as function of age
    if (binary_covariate)
        rsbp <- rbinom(n, 1, 0.5)
            else
 rsbp <- sapply(age, function(x) ifelse(x <= agecutoff_rsbp, rnorm(1, meanrsbp[1], sdrsbp[1]), rnorm(1, meanrsbp[2], sdrsbp[2]) ) )

simd <- data.frame(age, rsbp)

    site <- gl(K, n, labels = 1:K)
    
rnd_int <- rnorm(K, 0, rdn_sd[1])
rnd_slp <- rnorm(K, 0, rdn_sd[2])
rnd_slp2 <- rnorm(K, 0, rdn_sd[3])
    
groupA <-  which(age <= agecutoff_y)
groupB <-  which(age > agecutoff_y)

mus <- as.vector(matrix(nrow = n*K))
    
 mus[groupA] <- beta0[1] + rnd_int[site][groupA] + (betarsbp[1] + rnd_slp[site][groupA])*simd$rsbp[groupA] + (betaage[1] + rnd_slp2[site][groupA]*simd$age[groupA])

    mus[groupB] <- beta0[2] + rnd_int[site][groupB] + (betarsbp[2] + rnd_slp[site][groupB])*simd$rsbp[groupB] + (betaage[2] + rnd_slp2[site][groupB]*simd$age[groupB])

simd$death <- rbinom(K*n, 1, logistic(mus))

    simd$site <- site 
    
      return(simd)    

  }

##
                                      
generateArtificialDataOneSite <- function(simd, SI_k = 20000, meth = 3, avoid.zero.sd = TRUE){
    
 datchek <- Return.key.IPD.summaries( Return.IPD.design.matrix(simd) )
 zerosd <-  which( datchek$first.four.moments[, 2] == 0)
    if (length(zerosd) > 0 )
    simd2 <-  simd[, -zerosd]  # drop zero SD variables
    else
        simd2 <- simd    

    out <- try( Simulate.many.datasets(list(simd2), H = NULL, meth, stochastic.integration = TRUE, SI_k = SI_k, checkdata = TRUE, tabulate.similar.data = TRUE ), silent = TRUE)

    if (class(out[[1]]) == "try-error")
        return( list( list("similar.data" = lapply(1:ifelse(dim(simd2)[1] < 500, 300, 100), function(x) replace(simd2[1,], values = NA) )   ) ) )
    
    ## replace dropped zero SD variables with its mean value ...
    if (length(zerosd) > 0 & length(zerosd) < 2)
        out[[1]]$similar.data <- lapply(out[[1]]$similar.data, function(x){
       x <- as.data.frame(x)
       x[[names(zerosd)]] <- datchek$first.four.moments[zerosd, 1]
       if (avoid.zero.sd)
        x[[names(zerosd)]][sample(1:dim(simd2)[2], 1)]  <- ifelse(datchek$first.four.moments[zerosd, 1] == 0, 1, 0)  
      return(as.matrix(x))
  }  )
  
                         
    return(out)
    
  }

                                        #
PoolArtificialDataBySite <- function(obj, hubsvariable, hubsnames){
    
    H <- length(obj[[1]][[1]][["similar.data"]])

      listofpooled <-  lapply(1:H, function(h){      # by simulation

       merged <- na.omit( do.call("rbind",
         lapply(1:length(hubsnames), function(s){  ##  bind by country

           ad <- as.data.frame(obj[[s]][[1]][["similar.data"]][[h]])  # matrix simul
             ad <- data.frame(ad, hubsnames[s])
             colnames(ad)[dim(ad)[2]] <- hubsvariable
             return(ad) 

             } )  ) )
       
        return(merged)                     
                 
      } )
    
  return(listofpooled)
    
        }

                                        #

generateArtificialDataMultiCenter <- function(simd, site = "site", SI_k = 20000, meth = 3, ncores = 3, seed = 19, avoid.zero.sd = TRUE){

    sitenames <- levels(as.factor(simd[[site]]))
    
    datlist <- mclapply(sitenames, function(s){

   set.seed( seed + which(sitenames == s), "L'Ecuyer") 
           print( system.time(    
              res <-  generateArtificialDataOneSite(simd[simd[[site]] == s, -which(names(simd) == site)], SI_k, meth, avoid.zero.sd)
           ))
        return(res)
    }, mc.cores = ncores  )
              
                                        # pool data  
 out <-  PoolArtificialDataBySite(datlist, site, sitenames)

    return( list(pooled = out, details = datlist ) )
 
 }

##

                                     

trimOutliers <- function(x, trim = TRUE){   # discard outliers

   if (!trim)
       return(x)
    
 m <- median(x,na.rm=T)
   s <- mad(x, na.rm=T)
        
    if (all(is.na(m)) | all(is.na(s))) {
        warning("discard: cannot compute median or mad (all x values missing)")
        return(x)
     }
        
    if(max(x, na.rm=T) > (m+s*3.5) | min(x,na.rm=T) < (m-s*3.5)){
    
    min <- m-(2.5*s)
    max <- m+(2.5*s)

    res <- x[x>=min & x<=max]
        }else{
 res <- x
        }
    res
   }


####

NozeroSD <- function(simd, site = "site", target = "death"){

    sitenames <- levels(as.factor(simd[[site]]))
    
    out <- lapply(sitenames, function(s){
      res <-  simd[simd[[site]] == s, which(names(simd) == target)]
        if ( sd(res, na.rm = TRUE) == 0 )
        res[sample(1:dim(simd)[2], 1)]  <- ifelse(unique(res) == 0, 1, 0)

      return(res)
        }  )
    
   simd[[target]] <- unlist(out)
    return(simd)
  }

                                        #


    CollectMAsummaries <- function(simd, formula = NULL, site = "site", ncores = 3, noMA = FALSE){

        out <- NULL
                if (!is.null(formula)){
                    
sitenames <- levels(as.factor(simd[[site]]))
                    
                  out <- lapply(sitenames, function(s){            
 statf( simd[simd[[site]] == s, -which(names(simd) == site)], formula)
    } ); names(out) <- sitenames

       MAdat <- lapply(names(out[[1]]$coef), function(b) do.call("rbind", lapply(sitenames, function(s) data.frame( beta = out[[s]]$coef[b], betaSD = out[[s]]$sd[b], study = s) )))
    names(MAdat) <- names(out[[1]]$coef)

                }
        
        if (is.null(out)){
            sitenames <- names(simd)
            out <- simd
        }
        
      if (noMA)
     MAdat <- do.call( "rbind", lapply(c("coef", "sd"), function(i)  do.call( "rbind", lapply(names(out[[1]]$coef), function(b) do.call( "rbind", lapply(sitenames, function(s) data.frame( est = switch(i, "coef" = "$\\hat\\beta$", "sd" = "$\\hat\\sigma_{\\hat\\beta}$"), variab = b, site = s, ipd.ref =  out[[s]][[i]][b] )  )   ) ) ) ) )
      
                                         
        return(MAdat)

     }
    
                                        #

MultiplMA <- function(MAdat, MAeffect, statist){

        MAinfer <- lapply(MAdat, function(x) statma(x$beta, x$betaSD, x$study, MAeffect, statist))

        variab <- names(MAdat)
    params <- names(MAinfer[[1]])
    
    res <- lapply(params, function(b){
        out <- unlist(lapply(MAinfer, function(x) x[[b]]) )
        names(out) <- variab
        return(out)    }); names(res) <- params

    return(res)
    
   }

#

  BagInfer <- function(inferencesSample, operator = c("mean", "sd", "95qq"), trim.ol = FALSE ){

      operator <- match.arg(operator)
      
      op <- function(type, x) switch(type,
                                            "mean" = mean(x, na.rm = TRUE),
                                     "sd" = sd(x, na.rm = TRUE),
                   "95qq" = quantile(x, c(0.025, 0.975), na.rm = TRUE) )
     
      betas <- names(inferencesSample[[1]])
      
      out <- lapply(betas, function(b) apply(do.call("rbind", lapply(inferencesSample, function(y) y[[b]]) ), 2, function(x) op(operator, trimOutliers(x, trim.ol)))); names(out) <- betas

      return(out)
      
       }

                                        #
 
IPDboot <- function(simd, R = 100, distributed = TRUE, site = "site", statist = c("random", "fixed"), reg_formula_fixed = death ~ rsbp + age, reg_formula_random = death ~ rsbp + age + (1 + rsbp | site), ncores = 3, avoid.zero.sd = FALSE, trim.ol = FALSE ){

    statist <- match.arg(statist)

    if (distributed){  # distributed boot data 
        sitenames <- levels(as.factor(simd[[site]]))
        
            out <- lapply(sitenames, function(s) lapply(1:R, function(h){
      res <-  simd[simd[[site]] == s, -which(names(simd) == site)]
        resboot <- res[sample(1:dim(res)[1], dim(res)[1], TRUE), ]         
            return(resboot)
        }
                                                           )  )
    bootdat <- lapply(1:R, function(h) do.call("rbind", # merge by site 
         lapply(1:length(sitenames), function(s){  
           bd <- as.data.frame(out[[s]][[h]])  # matrix simul
             bd <- data.frame(bd, sitenames[s])
             colnames(bd)[dim(bd)[2]] <- site             
           return(bd)          } )  ) )

    } else  # pooled bood data
        bootdat <- lapply(1:R, function(h) simd[sample(1:dim(simd)[1], dim(simd)[1], TRUE), ]   )     


  if(avoid.zero.sd & distributed) # pooled version do not need this correction
      bootdat <- lapply(bootdat, function(x) NozeroSD(x, site))
    
# bootstrap statistics    
     IPDboot <-  mclapply(bootdat, function(x) switch(statist,
    "fixed" = statf(x, reg_formula_fixed),
    "random" = statm(x, reg_formula_random)                                                                  ), mc.cores = ncores)

        IPDbag <- BagInfer(IPDboot, "mean", trim.ol)
    IPDbagSE <- BagInfer(IPDboot, "sd", trim.ol)
    IPDbag95QQ <- BagInfer(IPDboot, "95qq", trim.ol)

   rm(bootdat)  # memory clean-up
    
    return( list( IPDbag = IPDbag, IPDbagSE = IPDbagSE, IPDbag95QQ = IPDbag95QQ) )

}


                                        #

     
    OneSimulationRun <- function(n, beta0, betarsbp, betaage, K, rdn_sd, meanrsbp, sdrsbp, meanage, sdage, cor_rsbpage, reg_formula_fixed = death ~ rsbp + age, reg_formula_random = death ~ rsbp + age + (1 | site), statist = c("random", "fixed"), MAeffect = "OR", scenario = c("one", "two"), agecutoff_rsbp = 40, agecutoff_y = 50, site = "site", SI_k = 20000, ncores = 3, seedAD = 19, avoid.zero.sd = FALSE, trim.ol = FALSE, binary_covariate = FALSE){


  statist <- match.arg(statist)
    scenario <- match.arg(scenario)

    ## simulate IPD
        
     simd <- switch(scenario,
                    "one" = generateDataScenarioOne(n, beta0, betarsbp, betaage, K, rdn_sd, meanrsbp, sdrsbp, meanage, sdage, cor_rsbpage, binary_covariate),
                    "two" = generateDataScenarioTwo(n, beta0, betarsbp, betaage, K, rdn_sd, meanrsbp, sdrsbp, meanage, sdage, cor_rsbpage, agecutoff_rsbp, agecutoff_y, binary_covariate) )
  
                    
                                        # generate artificial data
  ad <- lapply(3:4, function(j) generateArtificialDataMultiCenter(simd, site, SI_k, j, ncores, seedAD, avoid.zero.sd) )


 inferencesSample <- lapply(ad, function(j) mclapply(j$pooled, function(x) switch(statist,
                        "fixed" = statf(x, reg_formula_fixed),
                        "random" = statm(x, reg_formula_random)                                                                  ), mc.cores = ncores)
                        )
        
        GCinfer <- lapply( inferencesSample, function(j) BagInfer(j, "mean", trim.ol))
        GCinferSE <- lapply( inferencesSample, function(j) BagInfer(j, "sd", trim.ol))
        GCinfer95QQ <- lapply( inferencesSample, function(j) BagInfer(j, "95qq", trim.ol))
        names(GCinfer) <- names(GCinferSE) <- names(GCinfer95QQ) <- c("MV", "MVSK")
        
# two type of IPD bootstraps (distributed vs not distributed)
        IPDboots <- IPDboot(simd, length(ad[[1]]$pooled), distributed = FALSE, site, statist, reg_formula_fixed, reg_formula_random, ncores, trim.ol)
        
        IPDbootsDistr <- IPDboot(simd, length(ad[[1]]$pooled), distributed = TRUE, site, statist, reg_formula_fixed, reg_formula_random, ncores, avoid.zero.sd)
        
        if (avoid.zero.sd)
     simd <- NozeroSD(simd, site)
                                               
   IPDinfer <- switch(statist,   # original IPD infer
                        "fixed" = statf(simd, reg_formula_fixed),
                     "random" = statm(simd, reg_formula_random)  )
        # analytic 95% CIs                                             
        IPDinfer2.5QQ <- lapply(IPDinfer, function(x) replace(x, values = NA)); IPDinfer2.5QQ$coef <- CI(IPDinfer$coef, IPDinfer$sd, "low")
              IPDinfer97.5QQ <- IPDinfer2.5QQ; IPDinfer97.5QQ$coef <- CI(IPDinfer$coef, IPDinfer$sd)
  IPDinferSE <- lapply(IPDinfer, function(x) replace(x, values = NA))

     MAdat <- CollectMAsummaries(simd, reg_formula_fixed, site, ncores)
     
    MAinfer <- MultiplMA(MAdat, MAeffect, statist)

        MAinfer2.5QQ <- IPDinfer2.5QQ; MAinfer2.5QQ$coef <- CI( MAinfer$coef,  MAinfer$sd, "low") 
            MAinfer97.5QQ <- IPDinfer2.5QQ; MAinfer97.5QQ$coef <- CI( MAinfer$coef,  MAinfer$sd)  
     MAinferSE <- IPDinferSE; MAinferSE$ranefsd <- MAinfer$ranefsdSE

 rm(ad)  # memory clean-up !!!
        
        return( list(GCgamma = GCinfer[[1]], GCgammaSE = GCinferSE[[1]], GCgamma2.5QQ = lapply(GCinfer95QQ[[1]], function(x) x["2.5%", ]), GCgamma97.5QQ = lapply(GCinfer95QQ[[1]], function(x) x["97.5%", ]), GCjohn = GCinfer[[2]], GCjohnSE = GCinferSE[[2]], GCjohn2.5QQ = lapply(GCinfer95QQ[[2]], function(x) x["2.5%", ]), GCjohn97.5QQ = lapply(GCinfer95QQ[[2]], function(x) x["97.5%", ]), IPD = IPDinfer, IPDSE = IPDinferSE, IPD2.5QQ = IPDinfer2.5QQ, IPD97.5QQ = IPDinfer97.5QQ, IPDboot = IPDboots$IPDbag, IPDbootSE = IPDboots$IPDbagSE, IPDboot2.5QQ = lapply(IPDboots$IPDbag95QQ, function(x) x["2.5%", ]), IPDboot97.5QQ = lapply(IPDboots$IPDbag95QQ, function(x) x["97.5%", ]), IPDbootDistr = IPDbootsDistr$IPDbag, IPDbootDistrSE = IPDbootsDistr$IPDbagSE, IPDbootDistr2.5QQ = lapply(IPDbootsDistr$IPDbag95QQ, function(x) x["2.5%", ]), IPDbootDistr97.5QQ = lapply(IPDbootsDistr$IPDbag95QQ, function(x) x["97.5%", ]), MA = MAinfer, MASE = MAinferSE, MA2.5QQ = MAinfer2.5QQ, MA97.5QQ = MAinfer97.5QQ,
###  GCdetails =  list(MV = ad[[1]]$pooled, MVSK = ad[[2]]$pooled),
                     GCdetails = list (MV = NULL, MVSK = NULL ),
                     IPDdetails = simd ))
        
}


  # 

SimulationRuns <- function(nruns, n, K, beta0, betarsbp, betaage, rdn_sd, meanrsbp, sdrsbp, meanage, sdage, cor_rsbpage, reg_formula_fixed = death ~ rsbp + age, reg_formula_random = death ~ rsbp + age + (1 | site), poolruns = TRUE, statist = c("random", "fixed"), MAeffect = "OR", scenario = c("one", "two"), agecutoff_rsbp = 40, agecutoff_y = 50, site = "site", SI_k = 20000, ncores = 3, seedHyper = 98, seedAD = 19, avoid.zero.sd = FALSE, trim.ol = FALSE, binary_covariate = FALSE){

      statist <- match.arg(statist)
    scenario <- match.arg(scenario)

    
    out <-  mclapply(1:nruns, function(j){  
           set.seed(seedHyper + j)       
        OneSimulationRun(n, beta0, betarsbp, betaage, K, rdn_sd, meanrsbp, sdrsbp, meanage, sdage, cor_rsbpage, reg_formula_fixed, reg_formula_random, statist, MAeffect, scenario, agecutoff_rsbp, agecutoff_y, site, SI_k, ncores, (seedAD + j), avoid.zero.sd, trim.ol, binary_covariate )
                         }, mc.cores = ncores )


    MCmean <- NULL
    MCsd <- NULL
    if (poolruns){
     names <- names(out[[1]])   
  MCmean <- lapply(names[-c(length(names)-1,length(names))], function(i) BagInfer(lapply(out, function(x) x[[i]])) );  names(MCmean) <- names[-c(length(names)-1,length(names))]
        MCsd <- lapply(names[-c(length(names)-1,length(names))], function(i) BagInfer(lapply(out, function(x) x[[i]]), "sd")  );  names(MCsd) <- names[-c(length(names)-1,length(names))]

    }

    return( list( MCmean = MCmean, MCsd = MCsd, details = out ))

   }


## CREATE TABLES


LatexNotation <- function(betas){

    out <- ifelse(betas == "coef", "$\\beta$",
           ifelse(betas == "sd", "$\\sigma_{\\beta}$",
           ifelse(betas == "ranefsd", "$\\sigma_{b}$",
           ifelse(betas == "ranefsdSE", "SD$(\\sigma_{b})$",
           ifelse(betas == "ranefcor", "$\\rho_{b}$",
           ifelse(betas == "aic", "{\\Small AIC}", NA)))
                                       )     ))
   return(out)           
    
}

                                        #

CI <- function(m, s, bound = c("up", "low"), alpha = 0.05){

    bound <- match.arg(bound)
    
zq <- qnorm(alpha/2, lower = FALSE)
    
    out <- switch(bound,

                  "up"  =  m + (zq*s),
                  "low" =  m - (zq*s)
                  )
     return(out)

}

## ## 


TabulateSimulationRes <- function(res, n, K, truevalues, scenario = c("one", "two"), ranefsd = TRUE, ranefcor = FALSE, aic = FALSE, rounddec = 2, agecutoff_y = 50, binary_covariate = FALSE){

        scenario <- match.arg(scenario)
    analyses <- names(res$MCmean)
    cmres <- analyses[-c(grep("SE", analyses), grep("QQ", analyses))]       #
  cmres.ord <- c( which(cmres == "IPD"), which(cmres == "MA"), grep("boot", cmres), grep("GC", cmres) )
  
   betas <- names(res$MCmean[[1]])
       variab <- names(res$MCmean[[1]][[1]]) 
    
    out <- lapply(variab, function(vv){        
   line <-  as.data.frame(do.call("cbind", lapply(cmres, function(cc) unlist(lapply(betas, function(b) ifelse(!is.na(res$MCmean[[cc]][[b]][[vv]]), paste0(round(res$MCmean[[cc]][[b]][[vv]], rounddec), " (", round(res$MCsd[[cc]][[b]][[vv]], rounddec) ,")" ), NA) )  ) ) ) )
   colnames(line) <- cmres; rownames(line) <- betas
                                        # TODO: trueval wrongly reported in table ...   
     MeanRoughLFP <- round(apply(do.call("cbind", lapply(res$details, function(x) as.numeric(as.character(RoughLFP(truevalues[[vv]], vv, agecutoff_y, x$IPDdetails, rounddec, binary_covariate))))), 1, mean, na.rm = T), rounddec)
        
   return( make.multirow(data.frame(copula.data = switch(scenario, "one" = "Well modeled", "two" = "Mismodeled"), n, K, params = LatexNotation(betas),
       trueval = factor(MeanRoughLFP, levels = MeanRoughLFP[!is.na(MeanRoughLFP)]),     #  trueval = truevalues[[vv]],
                                    LatexBoldNumbers(line, MeanRoughLFP)[, cmres.ord])[c(T, T, ranefsd, F, ranefcor, aic), ], 2:3, rotate = c(T,F), em = c(NA, 2))[, c(2:6, 8:9, 7, 1, 10:11)] ) # small error with em here it does no insert unit measure "em"
   
    }    ); names(out) <- variab

###
    
     out2 <- lapply(variab, function(vv){
   line <-  as.data.frame(do.call("cbind", lapply(cmres, function(cc) do.call("rbind", lapply(betas, function(b) do.call("rbind", lapply(c("SE", "2.5QQ", "97.5QQ"), function(qq) ifelse(!is.na(res$MCmean[[paste0(cc, qq)]][[b]][[vv]]), paste0(round(res$MCmean[[paste0(cc, qq)]][[b]][[vv]], rounddec), " (", round(res$MCsd[[paste0(cc, qq)]][[b]][[vv]], rounddec) ,")" ), NA) ) ) ) ) ) ) )[, cmres.ord]
   colnames(line) <- cmres[cmres.ord]

  cis <- line[-seq(1, dim(line)[1], 3), ]
  lowci <- cis[seq(1, dim(cis)[1], 2), ]
   upcis <- cis[-seq(1, dim(cis)[1], 2), ]

   
   cismap <- apply(t(matrix( apply( cbind( do.call("rbind", lapply(1:dim(upcis)[1], function(i) t(upcis[i,]))),  do.call("rbind", lapply(1:dim(lowci)[1], function(i) t(lowci[i,]))) ), 1, function(y) ifelse(any(is.na(sapply(y, function(x) ExtractNum(x) ))), NA, abs(diff(sapply(y, function(x) ExtractNum(x) ) ) ) ) ), ncol = 6)), 1, function(x){
       if(all(is.na(x)))
           return(NA)
       else
      return(which(x == min(x, na.rm = T)))
   }  )
   
   cis <- apply(cis, 2, as.character)
   lookup <- matrix( 1:12, ncol = 6)
   
   for(i in 1:length(cismap)){
       for(j in 1:length(cismap[[i]])){
           if (!is.na(cismap[[i]][[j]]))
               cis[lookup[, i], cismap[[i]][[j]]] <- sapply(cis[lookup[, i], cismap[[i]][[j]]], function(x) makeBold(x))
         
           }
   }

   
   line <- apply(line, 2, as.character)
   line[-seq(1, dim(line)[1], 3), ] <- cis
   line[seq(1, dim(line)[1], 3), ] <- apply(LatexBoldNumbers( line[seq(1, dim(line)[1], 3), ], rep(NA, dim(line[seq(1, dim(line)[1], 3), ])[1])), 2, as.character)
   
 return( make.multirow( data.frame(copula.data = switch(scenario, "one" = "Well modeled", "two" = "Mismodeled"), n, K, params = rep(LatexNotation(betas), rep(3,6)), distrsumm = rep(c("SD", "Q2.5th", "Q97.5th"), 6), as.data.frame(line) )[rep(c(T, T, ranefsd, F, ranefcor, aic), rep(3,6)), ], 2:4, rotate = c(T,F,F), em = c(NA, 2, 2) )[, c(2:6, 8:9, 7, 1, 10:11)] ) # small error with em here it does no insert unit measure "em"  

       }    ); names(out2) <- variab

   rm(res) # memory-clean-up
    
    return( list(central = out, dispersion = out2) )
    
   }

# 

RoughLFP <- function(truevalues, namesvariab, agecutoff_y, IPDdetails, rounddec, binary_covariate = FALSE ){
    
    age <- IPDdetails$age
    rsbp <- IPDdetails$rsbp

averbeta <- function(x, weight1, weight2){

        stop <-  which(strsplit(as.character(x), "")[[1]] == ",")
        if (length(stop) > 0)
          return( as.character( round( ( (as.numeric(substr(as.character(x), 1, (stop-1) ))*weight1) + (as.numeric(substr(as.character(x), (stop+2), nchar(as.character(x)) ))*weight2) )/(weight1+weight2), rounddec) ) )
        else
            return(x)               
    }
# 
            weights <- switch(namesvariab,
          "(Intercept)" = c(agecutoff_y-min(age), max(age)-(agecutoff_y+1)),
          "rsbp" = { if(binary_covariate)
                         c(agecutoff_y-min(age), max(age)-(agecutoff_y+1))
                         else
          c( max(rsbp[which(age <= agecutoff_y)]) - min(rsbp[which(age <= agecutoff_y)]),  max(rsbp[which(age > agecutoff_y)]) - min(rsbp[which(age > agecutoff_y)])) },  
    "age" = c(agecutoff_y-min(age), max(age)-(agecutoff_y+1))       
    )
    
            out <- as.factor(sapply(as.character(truevalues), function(x) averbeta(x, weights[1], weights[2]) ))
            
            return(out)
    
  }

# small function to disentangle numeric value from character string with format num1 (num2)


ExtractNum <- function(x, strspl = " ", charslide = -1, leftterm = TRUE){

    stop <-  which(strsplit(as.character(x), "")[[1]] == strspl)
    if (length(stop) > 0){
        out <- substr(as.character(x), ifelse(leftterm, 1, stop), ifelse(leftterm, (stop + charslide), nchar(as.character(x))) )
        return( ifelse(leftterm, as.numeric(out), out) )
    }else
        return(x)
}

                                        #
    minDist <- function(x, trueval = NA){
                                       
        if (is.na(x))
            return(NA)
                
   if (is.na(trueval))
  return( ExtractNum(x) )      
       else
           return( abs(ExtractNum(x) - as.numeric(as.character(trueval)) ) )          
    }
    #
    makeBold <- function(x){

        if (is.na(x))
            return(NA)
                
  new <- paste0("\\textbf{",ExtractNum(x), "}", ExtractNum(x, leftterm = FALSE))
        return(new)
        }

# small function to bolden values with minimum distance from true value, or to bolden minimum values if true value is NA.

LatexBoldNumbers <- function(out, truevalues, bound = c("min", "max")){

    bound <- match.arg(bound)

    target <- apply(cbind(truevalues, out), 1, function(y){
        
   outi <- sapply( y[-1], function(x) minDist(x, y[1]) )
        
        if (all(is.na(outi)))
            return(NA)
        else
       return( which(outi == switch(bound, "min" = min(outi, na.rm = T), "max" = max(outi, na.rm = T)) ) )
        }
   )
    
    out <- apply(out, 2, as.character)
    
    for(i in 1:length(target)){ # target can be list due to multiple minima
        for( j in 1:length(target[[i]]) )
            out[i, target[[i]][j]] <- makeBold(ifelse(is.na(target[[i]][j]), NA, out[i, target[[i]][j]]))        
    }
   
   return(as.data.frame(out))
    
  }

##

TabulateSimulations <- function(nruns, params, beta0_I = -6, beta0_II = c(-2, -11), betarsbp_I = 0.05, betarsbp_II = c(0.01, 0.08), betaage = c(0.02, 0.02), rdn_sd = c(1.3, 0.0, 0.0), meanrsbp = c(90, 150), sdrsbp = c(10, 35), meanage = 60, sdage = 10, cor_rsbpage = 0.35, reg_formula_fixed = death ~ rsbp + age, reg_formula_random = death ~ rsbp + age + (1 | site), poolruns = TRUE, statist = c("random", "fixed"), MAeffect = "OR", scenario = c("one", "two"), agecutoff_rsbp = 40, agecutoff_y = 50, site = "site", SI_k = 20000, ncores = 3, seedHyper = 98, seedAD = 19, avoid.zero.sd = FALSE, trim.ol = FALSE, table.ranefcor = FALSE, table.aic = FALSE, rounddec = 2, latextable = FALSE, binary_covariate = FALSE){
    
       statist <- match.arg(statist)
       
      truevalues <- function(scenario) switch(scenario,
    "one" = { out <- data.frame( c(beta0_I, NA, rdn_sd[1], NA, NA, NA),  c(betarsbp_I, NA, rdn_sd[2], NA, NA, NA),  as.character(c(betaage[1], NA, rdn_sd[3], NA, NA, NA))); colnames(out) <- c("(Intercept)", "rsbp", "age" ); out},
    "two" = {out <- data.frame( c(paste0(beta0_II[1],", ", beta0_II[2]), NA, rdn_sd[1], NA, NA, NA),  c(paste0(betarsbp_II[1],", ", betarsbp_II[2]), NA, rdn_sd[2], NA, NA, NA),  as.character(c(switch(as.character(length(unique(betaage))), "1" = betaage[1], "2" = paste0(betaage[1],", ", betaage[2])), NA, rdn_sd[3], NA, NA, NA))); colnames(out) <- c("(Intercept)", "rsbp", "age" ); out}
        )
    
    out <-  lapply(scenario, function(sc)  mclapply(1:dim(params)[1], function(i){   
   #     lapply(1:dim(params)[1], function(i){ # TWEAK HERE            
        res <- SimulationRuns(nruns, params[i, 1], params[i, 2], switch(sc, "one" = beta0_I, "two" = beta0_II), switch(sc, "one" = betarsbp_I, "two" = betarsbp_II), betaage, rdn_sd, switch(sc, "one" = mean(meanrsbp), "two" = meanrsbp), switch(sc, "one" = mean(sdrsbp), "two" = sdrsbp), meanage, sdage, cor_rsbpage, reg_formula_fixed, reg_formula_random, poolruns, statist, MAeffect, scenario = sc, agecutoff_rsbp, agecutoff_y, site, SI_k, ncores, (seedHyper + which(scenario == sc) + i), (seedAD + which(scenario == sc) + i), avoid.zero.sd, trim.ol, binary_covariate)
              
     TabulateSimulationRes(res, params[i, 1], params[i, 2], truevalues(sc), sc, switch(statist, "random" = TRUE, "fixed" = FALSE), table.ranefcor, table.aic, rounddec, agecutoff_y, binary_covariate) 

   }, mc.cores = ncores ) ); names(out) <- scenario

       
   OUT <- lapply(names(out[[1]][[1]]), function(nn){
       res <- lapply(names(out[[1]][[1]][[nn]]), function(b) do.call("rbind", lapply(1:length(scenario), function(i) make.multirow(do.call("rbind", lapply(1:dim(params)[1], function(j) out[[i]][[j]][[nn]][[b]]  )  ), which(names(out[[i]][[1]][[nn]][[b]]) == "copula.data"), rotate = which(names(out[[i]][[1]][[nn]][[b]]) == "copula.data") )   )  )  )
       names(res) <- names(out[[1]][[1]][[nn]])
       if (binary_covariate)
           names(res)[2] <- "rxasp"
       return(res)
   }); names(OUT) <- names(out[[1]][[1]])
       
       
       if (latextable){
           OUT <- lapply(names(OUT), function(y){
               res <- lapply(names(OUT[[y]]), function(x)
        LatexTable(OUT, x, MAeffect, switch(statist, "random" = TRUE, "fixed" = FALSE), table.ranefcor, table.aic, y ) )
               names(res) <- names(OUT[[y]])
                      if (binary_covariate)
           names(res)[2] <- "rxasp"
               return(res)
    } ); names(OUT) <- names(out[[1]][[1]]) }
       
  
     return(OUT)     ### NOTE: u need to cat result when calling externally               
  }  


##


LatexTable <- function(res, variab, MAeffect, table.ranef, table.ranefcor = FALSE, table.aic = FALSE, type = c("central", "dispersion")){

    type <- match.arg(type)

    tab <- switch(type,
                  "central" = res$central[[variab]],
                  "dispersion" = res$dispersion[[variab]] )
                      
    nrules <- switch(type,
                     "central" =  seq(length(unique(tab$params)), nrow(tab), length(unique(tab$params))),

                     "dispersion" =  seq((length(unique(tab$params))-1)*3, nrow(tab), (length(unique(tab$params))-1)*3 ))

    headers <- switch(type,
                      "central" = c(  " \\toprule \\multicolumn{4}{c}{ }& \\multicolumn{3}{c}{IPD available} & \\multicolumn{4}{c}{IPD not available}\\\\ [0.3 cm]
      $n$ & $K$ & \\rotatebox[origin=c]{80}{Param.} & True & IPD Ref. & Bootstrap & Dist. Boot. & Meta anal. & \\rotatebox[origin=c]{80}{Artif. data} & NORTA-$\\Gamma$ & NORTA-J  \\\\[0.2 cm] \\midrule ", rep("[0.4 cm]", length(nrules)-1), "\\bottomrule" ),
      
      "dispersion" = c(  " \\toprule \\multicolumn{4}{c}{ }& \\multicolumn{3}{c}{IPD available} & \\multicolumn{4}{c}{IPD not available}\\\\ [0.3 cm]
      $n$ & $K$ & \\rotatebox[origin=c]{80}{Param.} & Est. & IPD Ref.* & Bootstrap & Dist. Boot. & Meta anal.* & \\rotatebox[origin=c]{80}{Artif. data} & NORTA-$\\Gamma$ & NORTA-J  \\\\[0.2 cm] \\midrule ", rep("[0.4 cm]", length(nrules)-1), "\\bottomrule" )  )

      #TODO: * clarify what is analystic and empirical estimate, clarify what is reference true value in mismodelled scenario ..
    
    captions <- switch(type,
                      "central" = paste0("Simulation results (100 runs). Mean (Standard Error) estimate for the following parameters (Param.): log ",MAeffect," ($\\beta$) for variable '",variab,"', log ",MAeffect," standard deviation ($\\sigma_{\\beta}$)", ifelse(table.ranef, ", log random effect standard deviation ($\\sigma_{b}$)", "."), ifelse(table.ranefcor, ", random effects correlation ($\\rho_{b}$)", "."), ifelse(table.aic, ", AIC.", "."), " Parameter estimates closest to the true value (True), or the minimum estimate, are printed in bold. Different estimation approaches according to IPD availability: original IPD regression estimates (IPD ref.), non-parametric bootstrap (Bootstrap), non-parametric bootstrap distributed across sites (Dist. Boot.), inverse-weighting meta-analysis (Meta anal.), proposed NORTA approaches (NORTA-$\\Gamma$, NORTA-J). IPD = individual person data. The true IPD generating mechanism is chosen to either agree or disagree with the NORTA model. Consequently NORTA artificial data (Artif. data) is respectively well- or mis-modeled.  $n$ = sample size at each site, $K$ = number of sites." ),
                      
                     "dispersion" =  paste0("Simulation results (100 runs). Mean (Standard Error) estimate of the empirical standard deviation (SD), 2.5th and 97.5th quantile (Q2.5th and Q97.5th) of the following parameters (Param.): log ",MAeffect," ($\\beta$) for variable '",variab,"', log ",MAeffect," standard deviation ($\\sigma_{\\beta}$)", ifelse(table.ranef, ", log random effect standard deviation ($\\sigma_{b}$)", "."), ifelse(table.ranefcor, ", random effects correlation ($\\rho_{b}$)", "."), ifelse(table.aic, ", AIC.", "."), " The smallest SD or quantile interval is printed in bold. Different estimation approaches according to IPD availability: non-parametric bootstrap (Bootstrap), non-parametric bootstrap distributed across sites (Dist. Boot.), inverse-weighting meta-analysis (Meta anal.), proposed NORTA approaches (NORTA-$\\Gamma$, NORTA-J). IPD = individual person data. The true IPD generating mechanism is chosen to either agree or disagree with the NORTA model. Consequently NORTA artificial data (Artif. data) is respectively well- or mis-modeled.  $n$ = sample size at each site, $K$ = number of sites.")                    
                      )
        
out <- xtable(tab, digits = rep(0, dim(tab)[2]+1),
             align = c( rep("c", 4 ), "c|", "c", "c", "c|", rep("c", 4)),  caption = captions, label = paste0("tab:", type, ":", variab))
    
  print(out,
      include.rownames = FALSE, 
      include.colnames = FALSE, 
      hline.after = NULL ,
      table.placement = "h", 
      floating=FALSE,
            caption.placement="top",
      sanitize.text.function = force,
      tabular.environment = "longtable",
      add.to.row = list(pos = lapply( c(0, nrules), function(x) x ),
      command = headers  )  )              
                           
    

  }



                                        # memory check from: https://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session

.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           format(utils::object.size(x), units = "auto") })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Length/Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}





