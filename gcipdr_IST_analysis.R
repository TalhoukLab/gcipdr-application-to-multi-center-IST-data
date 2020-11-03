## gcipdr_IST_analysis.R contains R commands to execute 'gcipdr application to multi-center IST data' (from omonimous repository).
## Copyright (C) 2019 Federico Bonofiglio

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

if (!require("devtools")) {
    install.packages("devtools")
    library(devtools)
}


## Install 'JohnsonDistribution' dependency (only available on CRAN archives)

url <- "https://cran.r-project.org/src/contrib/Archive/JohnsonDistribution/JohnsonDistribution_0.24.tar.gz"
 pkgFile <- "JohnsonDistribution_0.24.tar.gz"
 download.file(url = url, destfile = pkgFile)

 install.packages(pkgs=pkgFile, type="source", repos=NULL)

 unlink(pkgFile)


### INSTALL package (install other dependencies manually if needed)

install_github("bonorico/gcipdr")


# load libraries

library(gcipdr)

if (!require("cowplot")) {
    install.packages("cowplot")
    library(cowplot)
}

if (!require("meta")) {
    install.packages("meta")
    library(meta)
}

if (!require("metafor")) {
    install.packages("metafor")
    library(metafor)
}

if (!require("lme4")) {
    install.packages("lme4")
    library(lme4)
}

library(dplyr)
library(forcats)

## do not run !
## options("mc.cores") <- 3 ## set up your global 'mclapply()' forking options other than default (2 cores). Beware: changing default options might not reproduce results as shown in paper.


## DOWNLOAD IST DATASET


url2 <- "https://datashare.is.ed.ac.uk/bitstream/handle/10283/128/IST_corrected.csv?sequence=5&isAllowed=y"

ist <- read.csv( url2 )


 ##++++++++++++++++++++++++++++++++++++++++++++++++++++

dim(ist)

sort(names(ist))

# str(ist)

tokeep <- c( "RCONSC", "SEX", "AGE", "RATRIAL", "RXASP", "RXHEP", "RSBP", "DRSISC", "DRSISCD", "FDEAD", "FDEADD", "COUNTRY")

ist <- ist[, tokeep]


# levels(ist$COUNTRY)


world.macroregions <- c( "AMER-SOUTH", "SOUTH-PAC", "EU-WEST", "EU-WEST", "AMER-SOUTH", "AMER", "AMER-SOUTH", "EU-EAST", "EU-NORTH", "EU-NORTH", "EU-NORTH", "EU-WEST", "EU-SOUTH", "ASIA-SOUTH", "EU-EAST", "INDIA", "MID-EST", "EU-SOUTH", "ASIA", "EU-WEST", "SOUTH-PAC", "EU-NORTH", "EU-EAST", "EU-SOUTH", "EU-EAST", "ASIA-SOUTH", "EU-EAST", "EU-EAST", "ASIA", "EU-SOUTH", "ASIA-SOUTH", "EU-NORTH", "EU-WEST", "EU-SOUTH", "EU-NORTH", "AMER" )


# levels(ist$COUNTRY) <- world.macroregions
ist <- ist %>%
  mutate(COUNTRY = factor(case_when(
    COUNTRY %in% c("CANA", "USA") ~ "AMER",
    COUNTRY %in% c("ARGE", "BRAS", "CHIL") ~ "AMER-SOUTH",
    COUNTRY %in% c("JAPA", "SOUT") ~ "ASIA",
    COUNTRY %in% c("HONG", "SING", "SRI") ~ "ASIA-SOUTH",
    COUNTRY %in% c("CZEC", "HUNG", "POLA", "ROMA", "SLOK", "SLOV") ~ "EU-EAST",
    COUNTRY %in% c("DENM", "EIRE", "FINL", "NORW", "SWED", "UK") ~ "EU-NORTH",
    COUNTRY %in% c("GREE", "ITAL", "PORT", "SPAI", "TURK") ~ "EU-SOUTH",
    COUNTRY %in% c("AUST", "BELG", "FRAN", "NETH", "SWIT") ~ "EU-WEST",
    COUNTRY %in% c("INDI") ~ "INDIA",
    COUNTRY %in% c("ISRA") ~ "MID-EST",
    COUNTRY %in% c("AUSL", "NEW") ~ "SOUTH-PAC",
    TRUE ~ NA_character_
  )))

table(ist$COUNTRY)


## rename levels

ist$RCONSC <- factor(ist$RCONSC, labels = c(1, 0, 2))
# levels(ist$RCONSC) <- c( 1, 0, 2 )

ist$SEX <- factor(ist$SEX, labels = c(0, 1))
# levels(ist$SEX) <- c( 0, 1)

ist$RATRIAL <- factor(ist$RATRIAL, levels = c("N", "Y"), labels = c(0, 1))
# levels(ist$RATRIAL) <- c(NA, 0, 1)

ist$RXASP <- factor(ist$RXASP, labels = c(0, 1))
# levels(ist$RXASP) <- c(0, 1)

ist$RXHEP <- fct_collapse(
  ist$RXHEP,
  `0` = "N",
  `1` = c("H", "L", "M")
) %>%
  relevel(ref = "0")
# levels(ist$RXHEP) <- c(1, 1, 1, 0)
# ist$RXHEP <- relevel(ist$RXHEP, ref = "0")

ist$FDEAD <- factor(ist$FDEAD, levels = c("N", "Y"), labels = c(0, 1))
# levels(ist$FDEAD) <- c(NA, 0, NA, 1)

ist$is.EU <- grepl("EU", ist$COUNTRY)

## select only EU countries

ist <- ist[ist$is.EU, ]
ist <- ist[ , -c( 8,9, 11, 13 )]
ist[ ,-c(1,9)] <- apply(ist[ ,-c(1,9)], 2, as.integer)  ## convert to numeric

ist$COUNTRY <- as.factor(as.character(ist$COUNTRY))
ist$COUNTRY <- factor(ist$COUNTRY, levels(ist$COUNTRY)[c(4,2,1,3)])

ist$RCONSC <- relevel(ist$RCONSC, ref = "0")  #### IMPORTANT HERE !!! relevel RCONSC ..


ist <- ist[ which(!is.na(ist$RATRIAL)),  ]

head(ist)

hubs <- levels(ist$COUNTRY)


######### LOAD COPULA DATA ##########


 seed <- 49632

#
istsimul.by.region <- lapply( 3:4, function(i)
                     lapply(hubs, function(j){

                         data <- ist[ist$COUNTRY == j, -9]

   jiseed <- as.integer(paste(which(hubs == j), seed, i, sep = ""))
 set.seed( jiseed, "L'Ecuyer") # delete this line to assess stability
           print( system.time(
 artificial_data_object <-  Simulate.many.datasets(list(data), H = NULL, i,
      checkdata = TRUE, tabulate.similar.data = TRUE, NI_maxEval = 0 )
 ))
                         return( artificial_data_object[[1]])
          }  )         ); names(istsimul.by.region) <- paste0("meth_", 3:4)


  for(i in 1:2) names(istsimul.by.region[[i]]) <- hubs

### NOTE: One first computes the input IPD summaries and then feed them to 'DataRebuild'. 'Simulate.many.datasets' is just a handy wrapper for 'DataRebuild' that executes these two steps automatically (read 'gcipdr' documentation)

### NOTE: method 3 and 4 is a code for NORTA-$\Gamma$ and NORTA-J respectively


## ## ###  artificial POOLED IST function

    PoolArtifData <-function(istsimul.by.region){

    out <- lapply(1:2, function(i) lapply(1:100, function(h){

         merged <- do.call("rbind",
         lapply(1:length(hubs), function(j){  ##  bind by country

  data <- as.data.frame(istsimul.by.region[[i]][[j]]$similar.data[[h]])  # matrix simul
             data$COUNTRY <- hubs[j]
             return(data)
             } )  )

        return(merged)
    })    ); names(out) <- paste0("meth_", 3:4)

        return(out)
    }



##############################################
 ########## IPD INFERENCE RECOVERY  ##########
##############################################

### GLM RANDOM EFFECT (IST)


## original IST analysis

istformula <- FDEAD ~ RXASP*RXHEP + RSBP + RATRIAL + RCONSC + SEX + AGE + (1|COUNTRY)

 origIPDreinfer <- glmer(  istformula , data = ist, family = "binomial", nAGQ = 0,
                         control = glmerControl(optimizer = "Nelder_Mead", boundary.tol = 1e-3) )

oomist <- summary(origIPDreinfer)

  origbeta <- oomist$coef[, 1]
 origsd <- oomist$coef[, 2]
origsdranef <- sqrt(as.numeric(oomist$varcor))
origaic <- oomist$AICtab[1]

#

istform <- FDEAD ~ RXASP*RXHEP + RSBP + RATRIAL + RCONSC1 + RCONSC2 + SEX + AGE + (1|COUNTRY)



### inference calculation function

stat <- function( x ){

  model <-     try(
 glmer(  istform , data = x, family = "binomial", nAGQ = 0,
                         control = glmerControl(optimizer = "Nelder_Mead", boundary.tol = 1e-3) )
                    , silent = TRUE      )

   if (class(model) == "try-error")
       return( list( coef= NA, sd = NA, ranefsd= NA, aic = NA ))

     out <- summary(model)

     ranefsd <- sqrt(as.numeric(out$varcor))

  return(
      list( coef= out$coef[, 1], sd = out$coef[, 2], ranefsd= ranefsd, aic = out$AICtab[1] )
      )
    }


# IST GLM random effect regressin estimates recovery


  istmeth3 <- mclapply( PoolArtifData(istsimul.by.region)[["meth_3"]], function(x) stat(x) )


  istmeth4 <- mclapply( PoolArtifData(istsimul.by.region)[["meth_4"]], function(x) stat(x) )

#


istmeth3 <- lapply(c("coef", "sd", "ranefsd", "aic"), function(j){
do.call("rbind",
       lapply(1:100, function(i){

           istmeth3[[i]][[j]]

       }) )
    }
       )



istmeth4 <- lapply(c("coef", "sd", "ranefsd", "aic"), function(j){
do.call("rbind",
       lapply(1:100, function(i){

           istmeth4[[i]][[j]]

       }) )
    }
       )



 ###

istinferm3 <- lapply(istmeth3, function(x) apply(x, 2, mean, na.rm = T))

istinferm4 <- lapply(istmeth4, function(x) apply(x, 2, mean, na.rm = T))
#

 variables <- names(origbeta)



### MC standard deviation
# log OR SD

MCsdmeth3 <- apply(istmeth3[[1]], 2, sd, na.rm = T)
MCsdmeth4 <- apply(istmeth4[[1]], 2, sd, na.rm = T)

table3b <- data.frame( vars= variables, data = "IST", meth3 = MCsdmeth3, meth4 = MCsdmeth4 )

    # SD of log RE SD

MCsdREmeth3 <- apply(istmeth3[[3]], 2, sd, na.rm = T)
MCsdREmeth4 <- apply(istmeth4[[3]], 2, sd, na.rm = T)

sdranefsd <- c(MCsdREmeth3, MCsdREmeth4)


# TABLE 2 OF MAIN MANUSCRIPT


table2 <- data.frame( names= c(rep(variables,2), "EU-REGION", "") , est= c(rep("$\\hat\\beta$", 10), rep("$\\hat\\sigma_{\\hat\\beta}$", 10), "$\\hat\\sigma_{\\hat b}$",  "{\\Small AIC}") , original = c( origbeta, origsd, origsdranef, origaic), meth3= unlist(istinferm3), meth4 = unlist(istinferm4))


 table2 <- make.multirow( table2, 2, em = "1em" )


### TABLE 1 OF SUPPLEMENT PART I


 table3 <- make.multirow( table3b, 2, rotate = T, em = "1em")


############# ADDITIONAL ANALYSES (REVISOIN I) ###########################
################################################################


source("simulation_functions_v9.R")


### STEP ONE: report location specific ORs --> IPD inverse weighting meta-analysis

 istformulafixed <- FDEAD ~ RXASP*RXHEP + RSBP + RATRIAL + RCONSC + SEX + AGE

orig_ORs <- CollectMAsummaries(ist, istformulafixed, "COUNTRY", ncores = 2, noMA = TRUE)

# 'mc.cores' > 1 is not supported on Windows
if (.Platform$OS.type == "windows") {
  outbootORs <- lapply(hubs, function(s){ set.seed(4 + which(hubs == s)); IPDboot(ist[which(ist$COUNTRY == s), ], distributed = FALSE, statist = "fixed", reg_formula_fixed = istformulafixed, ncores = 1)$IPDbag }); names(outbootORs) <- hubs
} else {
  outbootORs <- lapply(hubs, function(s){ set.seed(4 + which(hubs == s)); IPDboot(ist[which(ist$COUNTRY == s), ], distributed = FALSE, statist = "fixed", reg_formula_fixed = istformulafixed, ncores = 2)$IPDbag }); names(outbootORs) <- hubs
}

bootORs <- CollectMAsummaries(outbootORs, noMA = TRUE)


 istformfixed <- FDEAD ~ RXASP*RXHEP + RSBP + RATRIAL + RCONSC1 + RCONSC2 + SEX + AGE

out <- lapply(istsimul.by.region, function(s) lapply(s, function(d) BagInfer( lapply(d$similar.data, function(x) statf(x, istformfixed)  ) )  ) )

nortagam <- CollectMAsummaries(out$meth_3, noMA = TRUE)
nortajohn <- CollectMAsummaries(out$meth_4, noMA = TRUE)

##

alphalev <- 0.005

countryORs <- data.frame(orig_ORs, IPDboot = bootORs[ ,4],  meth_3 = nortagam[, 4], meth_4 = nortajohn[, 4])

# bold

signif <- countryORs[1:40, -c(1:3)]/ countryORs[41:80, -c(1:3)]
countryORs[, -c(1:3) ] <- round(countryORs[, -c(1:3) ], 3)

boldrepl <- as.vector(matrix( nrow = prod(dim(countryORs[1:40, -c(1:3)]))))

boldrepl[which(abs(signif) > abs(qnorm(alphalev)))] <- paste0("\\textbf{",as.numeric(as.matrix(countryORs[1:40, -c(1:3)]))[which(abs(signif) > abs(qnorm(alphalev)))],"}")

boldrepl[-which(abs(signif) > abs(qnorm(alphalev)))] <- as.character(as.numeric(as.matrix(countryORs[1:40, -c(1:3)]))[-which(abs(signif) > abs(qnorm(alphalev)))])

countryORs[, -c(1:3)] <- apply(countryORs[, -c(1:3)], 2, as.character)

replacem <- as.data.frame(matrix(boldrepl, ncol = dim(countryORs[1:40, -c(1:3)])[2])); names(replacem) <- names(countryORs)[-c(1:3)]

  replacem <- rbind(replacem, countryORs[41:80, -c(1:3)] )

countryORs[, -c(1:3)] <- replacem

countryORs <- countryORs[countryORs$variab == "RXASP", ]

countryORs <- make.multirow(countryORs, 1:2)


# meta-analysis

MAres <- MultiplMA( CollectMAsummaries(ist, istformulafixed, "COUNTRY", ncores = 2), "OR", "random")

table2$MA <- c(MAres$coef, MAres$sd, MAres$ranefsd[1], NA)


### STEP TWO: report IPD bootstrap
if (.Platform$OS.type == "windows") {
  bootres <- IPDboot(ist, distributed = FALSE, site = "COUNTRY", reg_formula_random = istformula, ncores = 1)
  bootresDistr <- IPDboot(ist, distributed = TRUE, site = "COUNTRY", reg_formula_random = istformula, ncores = 1)
} else {
  bootres <- IPDboot(ist, distributed = FALSE, site = "COUNTRY", reg_formula_random = istformula, ncores = 2)
  bootresDistr <- IPDboot(ist, distributed = TRUE, site = "COUNTRY", reg_formula_random = istformula, ncores = 2)
}


table2$IPDboot <- c(bootres$IPDbag$coef, bootres$IPDbag$sd, bootres$IPDbag$ranefsd[1], bootres$IPDbag$aic[1])

table2$IPDbootDistr <- c(bootresDistr$IPDbag$coef, bootresDistr$IPDbag$sd, bootresDistr$IPDbag$ranefsd[1], bootresDistr$IPDbag$aic[1])

table2 <- table2[, c(1:3, 7:8, 6, 4:5)]

table2$est[22] <- "\\multirow{1}{1em}{{\\small AIC}}"

 signif <- table2$original[2:10]/ table2$original[12:20]

  sigvars <- as.character(table2$names[2:10][which(abs(signif) > abs(qnorm(alphalev)))])

 for( i in sigvars)
       levels(table2$names)[which(levels(table2$names) == i)] <- paste0("\\textbf{",i,"}")


### integrate new info into toreport object

table3$IPDboot <- bootres$IPDbagSE$coef

table3$IPDbootDistr <- bootresDistr$IPDbagSE$coef

table3 <- table3[, c(1,5,6,3,4)]

sdranefsd <- c(bootres$IPDbagSE$ranefsd[1], bootresDistr$IPDbagSE$ranefsd[1], sdranefsd); names(sdranefsd) <- c("IPDboot", "IPDbootDistr", "meth_3", "meth_4")


## include 95% quantiles of sdranefsd

qqranefsd <- data.frame(IPDboot = bootres$IPDbag95QQ$ranefsd[,1], IPDbootDistr = bootresDistr$IPDbag95QQ$ranefsd[,1], meth_3 = quantile(istmeth3[[3]], c(0.025, 0.975), na.rm = TRUE), met_4 = quantile(istmeth4[[3]], c(0.025, 0.975), na.rm = TRUE))


### STEP THREE: report location specific RXASP OR for different AGE cutoffs (leq 60; > 60)

cutoff <- 60

scatter(na.omit(subset(ist, select = c("AGE", "FDEAD")))[, -9])

orig_ORs65age <- CollectMAsummaries(ist[ist$AGE <= cutoff , ], istformulafixed, "COUNTRY", ncores = 2, noMA = TRUE)

orig_ORs65age <- orig_ORs65age[orig_ORs65age$variab == "RXASP", ]

orig_ORs65age$variab <- paste0("AGE \u2266 ",cutoff)
# levels(orig_ORs65age$variab)[2] <- paste0("AGE \u2266 ",cutoff)

orig_ORs65ageplus <- CollectMAsummaries(ist[ist$AGE > cutoff, ], istformulafixed, "COUNTRY", ncores = 2, noMA = TRUE)

orig_ORs65ageplus <- orig_ORs65ageplus[orig_ORs65ageplus$variab == "RXASP", ]
orig_ORs65ageplus$variab <- paste0("AGE > ", cutoff)
# levels(orig_ORs65ageplus$variab)[2] <- paste0("AGE > ", cutoff)

ORs65 <- rbind(orig_ORs65age, orig_ORs65ageplus)


## beware: occasionally RCONSC2 levels disappear during data bootstrapping and coef vector is shorter --> rbind Error in BagInfer. I change func statf, use here a design.matrix like data.frame for ist with appropriate reg_formula (factoring all RCONCS levels explicetely)

if (.Platform$OS.type == "windows") {
  outbootORs65 <- lapply(hubs, function(s){ set.seed(86 + which(hubs == s)); IPDboot(subset(data.frame(FDEAD = na.omit(ist$FDEAD), model.matrix(istformulafixed, ist[, -9])[ ,-c(1,10)], COUNTRY = na.omit(ist)$COUNTRY), AGE <= 65 & COUNTRY == s), distributed = FALSE, statist = "fixed", reg_formula_fixed = istformfixed, ncores = 1)$IPDbag }); names(outbootORs65) <- hubs
} else {
  outbootORs65 <- lapply(hubs, function(s){ set.seed(86 + which(hubs == s)); IPDboot(subset(data.frame(FDEAD = na.omit(ist$FDEAD), model.matrix(istformulafixed, ist[, -9])[ ,-c(1,10)], COUNTRY = na.omit(ist)$COUNTRY), AGE <= 65 & COUNTRY == s), distributed = FALSE, statist = "fixed", reg_formula_fixed = istformfixed, ncores = 2)$IPDbag }); names(outbootORs65) <- hubs
}


bootORs65 <- CollectMAsummaries(outbootORs65, noMA = TRUE)
bootORs65 <-  bootORs65[bootORs65$variab == "RXASP", ]

if (.Platform$OS.type == "windows") {
  outbootORs65plus <- lapply(hubs, function(s){ set.seed(86 + which(hubs == s)); IPDboot(subset(data.frame(FDEAD = na.omit(ist$FDEAD), model.matrix(istformulafixed, ist[, -9])[ ,-c(1,10)], COUNTRY = na.omit(ist)$COUNTRY), AGE > 65 & COUNTRY == s), distributed = FALSE, statist = "fixed", reg_formula_fixed = istformfixed, ncores = 1)$IPDbag }); names(outbootORs65plus) <- hubs
} else {
  outbootORs65plus <- lapply(hubs, function(s){ set.seed(86 + which(hubs == s)); IPDboot(subset(data.frame(FDEAD = na.omit(ist$FDEAD), model.matrix(istformulafixed, ist[, -9])[ ,-c(1,10)], COUNTRY = na.omit(ist)$COUNTRY), AGE > 65 & COUNTRY == s), distributed = FALSE, statist = "fixed", reg_formula_fixed = istformfixed, ncores = 2)$IPDbag }); names(outbootORs65plus) <- hubs
}

bootORs65plus <- CollectMAsummaries(outbootORs65plus, noMA = TRUE)
bootORs65plus <-  bootORs65plus[bootORs65plus$variab == "RXASP", ]

boot65 <- rbind(bootORs65, bootORs65plus)

ORs65$IPDboot <- boot65[, 4]

###

out65 <- lapply(istsimul.by.region, function(s) lapply(s, function(d) BagInfer( lapply(d$similar.data, function(x) statf(subset(as.data.frame(x), AGE <= 65), istformfixed)  ) )  ) )

out65plus <- lapply(istsimul.by.region, function(s) lapply(s, function(d) BagInfer( lapply(d$similar.data, function(x) statf(subset(as.data.frame(x), AGE > 65), istformfixed)  ) )  ) )

nortagam65 <- CollectMAsummaries(out65$meth_3, noMA = TRUE)
nortagam65 <-  nortagam65[nortagam65$variab == "RXASP", ]
nortajohn65 <- CollectMAsummaries(out65$meth_4, noMA = TRUE)
nortajohn65 <-  nortajohn65[nortajohn65$variab == "RXASP", ]
nortagam65plus <- CollectMAsummaries(out65plus$meth_3, noMA = TRUE)
nortagam65plus <-  nortagam65plus[nortagam65plus$variab == "RXASP", ]
nortajohn65plus <- CollectMAsummaries(out65plus$meth_4, noMA = TRUE)
nortajohn65plus <-  nortajohn65plus[nortajohn65plus$variab == "RXASP", ]

nortagam65 <- rbind(nortagam65, nortagam65plus)
ORs65$meth_3 <- nortagam65[, 4]
nortajohn65 <- rbind(nortajohn65, nortajohn65plus)
ORs65$meth_4 <- nortajohn65[, 4]


# bold

ORs65 <- ORs65[order(ORs65$est), ]

signif <- ORs65[ORs65$est == "$\\hat\\beta$", -c(1:3) ]/ ORs65[ORs65$est == "$\\hat\\sigma_{\\hat\\beta}$", -c(1:3) ]

ORs65[, -c(1:3) ] <- round(ORs65[, -c(1:3) ], 3)

boldrepl <- as.vector(matrix( nrow = prod(dim(ORs65[ORs65$est == "$\\hat\\beta$", -c(1:3) ]))))

boldrepl[which(abs(signif) > abs(qnorm(alphalev)))] <- paste0("\\textbf{",as.numeric(as.matrix(ORs65[ORs65$est == "$\\hat\\beta$", -c(1:3) ]))[which(abs(signif) > abs(qnorm(alphalev)))],"}")

boldrepl[-which(abs(signif) > abs(qnorm(alphalev)))] <- as.character(as.numeric(as.matrix(ORs65[ORs65$est == "$\\hat\\beta$", -c(1:3) ]))[-which(abs(signif) > abs(qnorm(alphalev)))])

ORs65[, -c(1:3) ] <- apply(ORs65[, -c(1:3) ], 2, as.character)

replacem <- as.data.frame(matrix(boldrepl, ncol = dim(ORs65[ORs65$est == "$\\hat\\beta$", -c(1:3) ])[2])); names(replacem) <- names(ORs65)[-c(1:3)]

  replacem <- rbind(replacem, ORs65[ORs65$est == "$\\hat\\sigma_{\\hat\\beta}$", -c(1:3) ] )

ORs65[, -c(1:3)] <- replacem

ORs65 <- ORs65[order(ORs65$site), c(3,1,2,4:7)]

cols <- matrix(c(1,4,5,8,9,12,13,16), nrow = 2)
newcol <- c()
for( i in 1:4  )
  newcol <- c(newcol, make.multirow(ORs65[seq(cols[, i][1], cols[, i][2], 1), ], 2, em = "1em")[, 2] )

ORs65[, 2] <- as.factor(newcol)

ORs65 <- make.multirow(ORs65, 1)

### JUST KEEP eu-sud results here ....


################ CALIBRATION OF EU-SOUTH RXASP EFFECT ##############

inp <- Return.key.IPD.summaries(Return.IPD.design.matrix(subset(ist, COUNTRY == "EU-SOUTH", select = -9)))

# calibrate RXASP effect by changing input FDEAD-RXASP correlation
inp[[5]][9,6] <- inp[[5]][6,9] <- inp[[5]][9,6] - 0.055

marg <- "gamma"  # "johnson"

set.seed(0681, "L'Ecuyer")
pp <-  DataRebuild(100, inp[[1]], inp[[5]], inp[[4]], inp[[8]] , marg.mod = marg, variable.names = inp[[3]], checkdata = TRUE, tabulate.similar.data = TRUE, NI_maxEval = 0 )


eusud <- lapply(pp$Xspace, function(x) statf(x, istformfixed) )

BagInfer(eusud)

### original IPD analysis
statf(subset(ist, COUNTRY == "EU-SOUTH"), FDEAD ~ RXASP*RXHEP + RSBP + RATRIAL + RCONSC + SEX + AGE  )



##########################################

### STEP FOUR: report pooled fixed-effect logistic regression (REDUNDAT !!! DO NOT REPORT)

table2B <- table2[-21, 1:2]

pooledIPDfix <- statf(ist, istformulafixed)

table2B$original <- c(pooledIPDfix$coef, pooledIPDfix$sd, pooledIPDfix$aic[1] )

 if (.Platform$OS.type == "windows") {
   bootresf <- IPDboot(ist, distributed = FALSE, site = "COUNTRY", statist = "fixed", reg_formula_fixed = istformulafixed, ncores = 1)$IPDbag

   bootresfDistr <- IPDboot(ist, distributed = TRUE, site = "COUNTRY", statist = "fixed", reg_formula_fixed = istformulafixed, ncores = 1)$IPDbag
 } else {
   bootresf <- IPDboot(ist, distributed = FALSE, site = "COUNTRY", statist = "fixed", reg_formula_fixed = istformulafixed, ncores = 2)$IPDbag

   bootresfDistr <- IPDboot(ist, distributed = TRUE, site = "COUNTRY", statist = "fixed", reg_formula_fixed = istformulafixed, ncores = 2)$IPDbag
 }

table2B$IPDboot <- c(bootresf$coef, bootresf$sd, bootresf$aic[1])

table2B$IPDbootDistr <- c(bootresfDistr$coef, bootresfDistr$sd, bootresfDistr$aic[1])

MAresfix <- MultiplMA( CollectMAsummaries(ist, istformulafixed, "COUNTRY", ncores = 2), "OR", "fixed")

table2B$MA <- c(MAresfix$coef, MAresfix$sd, NA)

GCpooledfix <- lapply(PoolArtifData(istsimul.by.region), function(m) BagInfer(lapply(m, function(x) statf(x, istformfixed) )) )

table2B$meth_3 <- c(GCpooledfix$meth_3$coef, GCpooledfix$meth_3$sd, GCpooledfix$meth_3$aic[1])

table2B$meth_4 <- c(GCpooledfix$meth_4$coef, GCpooledfix$meth_4$sd, GCpooledfix$meth_4$aic[1])


table2B$est[21] <- "\\multirow{1}{1em}{{\\small AIC}}"

signif <- table2B$original[2:10]/ table2B$original[12:20]

sigvars <- as.character(table2B$names[2:10][which(abs(signif) > qnorm(alphalev))])

for( i in sigvars)
  levels(table2B$names)[which(levels(table2B$names) == i)] <- paste0("\\textbf{",i,"}")


############ ARTIFICIAL DATA COMPATIBILITY DIAGNOSTICS #########



##################### GRAPHICAL CHECK AGAINST ORIG. IPD ################


# summary(apply(ist[ist$COUNTRY == "EU", ], 1, function(x) any(is.na(x))))

 vertical.istsimul <- do.call( "rbind", lapply(c("meth_3","meth_4"),
       function(j) data.frame( do.call( "rbind", lapply(PoolArtifData(istsimul.by.region)[[j]], function(x) x) ),
                              DATA = j     ) ) )

 vertical.ist <- rbind( vertical.istsimul,
    na.omit( data.frame(Return.IPD.design.matrix(ist[, -9], fill = T)[, -1], COUNTRY = ist$COUNTRY, DATA = "Orig. IPD") ) )

 # levels( vertical.ist$DATA)[1:2] <- c("Fewer moments", "All moments")  # c("NORTA-\u0393", "NORTA-J")
 vertical.ist$DATA <- factor(
   vertical.ist$DATA,
   labels = c("Fewer moments", "All moments", "Orig. IPD")
 )

 vertical.ist$COUNTRY <- factor(vertical.ist$COUNTRY, hubs)


 mp <- aes( x= AGE, y=..density.., fill = DATA, color = DATA)
pl <- ggplot( vertical.ist[, c("AGE","COUNTRY","DATA")] %>%
                filter(!is.na(DATA)), mp )
# hist <- geom_histogram(alpha=0.2, position="identity")
hist <- geom_density(alpha=0.2, position="identity", show.legend = FALSE)

panel <- facet_wrap(~ as.factor(COUNTRY), ncol=2)
 xl <- xlab("AGE")

 hist_age <- pl + hist + panel + xl + theme_bw()

#

mp <- aes( x= RSBP, y=..density.., fill = DATA, color = DATA)
pl <- ggplot( vertical.ist[, c("RSBP","COUNTRY","DATA")] %>%
                filter(!is.na(DATA)), mp )
hist <- geom_density(alpha=0.2, position="identity")

panel <- facet_wrap(~ as.factor(COUNTRY), ncol=2)
 xl <- xlab("RSBP")

 hist_rsbp <- pl + hist + panel + xl + theme_bw() + theme(legend.position="bottom", legend.title = element_blank(), legend.text = element_text(size = 12) )



hist_ist <- plot_grid( hist_age, hist_rsbp, ncol = 1, rel_heights = c(2.5, 2.7) )



### FIGURE 1 OF MAIN MANUSCRIPT

pdf("hist_ist.eps", width = 5.7, height = 6.8)
## png("hist_ist.png", width = 550, height = 650)

hist_ist

dev.off()


rm(hist_ist , hist_age, hist_rsbp,  vertical.istsimul)

### 2D CONTOURS


if (!require("reshape2")) {
    install.packages("reshape2")
    library(reshape2)
}



 contour.data <- function(D, cname1, cname2, breaks = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), cumulative = TRUE){


   kk <- MASS::kde2d(D[, cname1], D[, cname2], n = 100)

    dx <- diff(kk$x[1:2])
    dy <- diff(kk$y[1:2])

    sz <- sort(kk$z)
    c1 <- cumsum(sz) * dx * dy

    dimnames(kk$z) <- list(kk$x,kk$y)
  dc <- melt(kk$z)
    colnames(dc)[1:2] <- c(cname1, cname2)

   dcs <- dc[order(dc$value), ]  ## order grid along cumulative sum

  Cdat <- data.frame(dcs, cumul.dens = 1-c1)

# (not run: check against 'contour' function in browser mode)

     ## Cdats <- Cdat[order(Cdat$AGE, Cdat$RSBP, decreasing = TRUE), ]

  ##    mykk <- list()

  ## mykk$x <- kk$x
  ##    mykk$y <- kk$y
  ##       mykk$z  <- matrix(Cdats$cumul.dens, ncol = 100)

  ##    plot(D$AGE, D$RSBP )
  ##      contour(mykk, levels = c(breaks, 1), add = T)

     contourbreaks <- NULL

     if (!cumulative)

 contourbreaks <- do.call("rbind", lapply(breaks, function(b)
   Cdat[ Cdat$cumul.dens == max(Cdat$cumul.dens[Cdat$cumul.dens <= b]),  ]                           ))

     rm(kk, dcs, dc, c1, sz, dy, dx)

 return(list(contours= Cdat, breaks = contourbreaks))

    }


#

 contourplot <- function(D, mp, breaks = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), cumulative = TRUE, legend.command ){

cname1 <- as.character(mp[1][[1]][2])
     cname2 <- as.character(mp[2][[1]][2])


 Cd <- contour.data( subset( D, subset = DATA == "Orig. IPD"), cname1, cname2, breaks, cumulative )

   Ds <- subset(D, subset = DATA != "Orig. IPD")

 pl <- ggplot( Ds , mp )
points <- geom_point( aes(colour = DATA), alpha = 0.2, show.legend = legend.command )

     if (cumulative)
      contour <- geom_contour(aes( z=cumul.dens ), data = Cd$contours, breaks = breaks, colour = "black" )
  else
    contour <- geom_contour(aes( z=value ), data = Cd$contours, breaks = Cd$breaks$value, colour = "black" )


 age_rsbp <- pl + points + theme_bw() + guides(colour = guide_legend(override.aes = list(alpha = 1))) + contour

     }


#

br <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)

 contourlist <- lapply(hubs, function(j){

     Ds <-  subset(vertical.ist, subset = COUNTRY == j)
 gr <- contourplot( Ds, aes(x = AGE, y = RSBP), cumulative = TRUE , legend.command = ifelse(levels(factor(Ds$COUNTRY)) == "EU-SOUTH", TRUE, FALSE),
                   breaks = br)

     if ( j == "EU-EAST")
    xl <- xlab("AGE")
 else
 xl <- xlab(" ")

if ( j == "EU-WEST")
  yl <- ylab("RSBP")
 else
yl <- ylab(" ")

 if ( j == "EU-SOUTH")
    legend <- theme(legend.position="bottom", legend.title = element_blank(), legend.text = element_text(size = 12))
 else
    legend <- NULL

  title <- annotate("text", label = j, x = min(Ds$AGE) + 10, y = max(Ds$RSBP))

    out <- gr + xl + yl + title + legend

 }  )


   names( contourlist) <- hubs



 allcontour <- plot_grid( contourlist[[1]], contourlist[[2]], contourlist[[3]], contourlist[[4]] , ncol = 2 )


#

 contour.condit.dead_1 <- lapply(hubs, function(j){

     Ds <-  subset(vertical.ist, subset = COUNTRY == j & FDEAD == 1)
 gr <- contourplot( Ds, aes(x = AGE, y = RSBP), cumulative = TRUE , legend.command = ifelse(levels(factor(Ds$COUNTRY)) == "EU-SOUTH", TRUE, FALSE),
                   breaks = br)

     if ( j == "EU-EAST")
    xl <- xlab("AGE")
 else
 xl <- xlab(" ")

if ( j == "EU-WEST")
  yl <- ylab("RSBP")
 else
yl <- ylab(" ")

      if ( j == "EU-SOUTH")
    legend <- theme(legend.position="bottom", legend.title = element_blank(), legend.text = element_text(size = 12))
 else
    legend <- NULL


  title <- annotate("text", label = j, x = min(Ds$AGE) + 10, y = max(Ds$RSBP))

    out <- gr + xl + yl + title + legend

 }  )


 allcontour.fdead1 <- plot_grid(  contour.condit.dead_1[[1]],  contour.condit.dead_1[[2]],  contour.condit.dead_1[[3]],  contour.condit.dead_1[[4]] , ncol = 2 )

#

 contour.condit.dead_0 <- lapply(hubs, function(j){

     Ds <-  subset(vertical.ist, subset = COUNTRY == j & FDEAD == 0)
 gr <- contourplot( Ds, aes(x = AGE, y = RSBP), cumulative = TRUE , legend.command = ifelse(levels(factor(Ds$COUNTRY)) == "EU-SOUTH", TRUE, FALSE),
                   breaks = br)

     if ( j == "EU-EAST")
    xl <- xlab("AGE")
 else
 xl <- xlab(" ")

if ( j == "EU-WEST")
  yl <- ylab("RSBP")
 else
yl <- ylab(" ")

 if ( j == "EU-SOUTH")
    legend <- theme(legend.position="bottom", legend.title = element_blank(), legend.text = element_text(size = 12) )
 else
    legend <- NULL

  title <- annotate("text", label = j, x = min(Ds$AGE) + 10, y = max(Ds$RSBP))

    out <- gr + xl + yl + title + legend

 }  )

#

allcontour.fdead0 <- plot_grid(  contour.condit.dead_0[[1]],  contour.condit.dead_0[[2]],  contour.condit.dead_0[[3]],  contour.condit.dead_0[[4]] , ncol = 2 )


### FIGURE 2 OF MAIN MANUSCRIPT

 png("contours.png", height = 600, width = 800)

 allcontour

dev.off()

rm(allcontour, contourlist)





### FIGURE 1 and 2 OF SUPPLEMENT PART I

png("contours_fdead1.png", height = 600, width = 800)

 allcontour.fdead1

dev.off()

rm( allcontour.fdead1, contour.condit.dead_1)

#
png("contours_fdead0.png", height = 600, width = 800)

 allcontour.fdead0

dev.off()


rm( allcontour.fdead0, contour.condit.dead_0)


##### (not shown in paper) check dependence conditional to sex


 contour.condit.sex_1 <- lapply(hubs, function(j){

     Ds <-  subset(vertical.ist, subset = COUNTRY == j & SEX == 1)
 gr <- contourplot( Ds, aes(x = AGE, y = RSBP), cumulative = TRUE , legend.command = ifelse(levels(factor(Ds$COUNTRY)) == "EU-SOUTH", TRUE, FALSE),
                   breaks = br)

     if ( j == "EU-EAST")
    xl <- xlab("AGE")
 else
 xl <- xlab(" ")

if ( j == "EU-WEST")
  yl <- ylab("RSBP")
 else
yl <- ylab(" ")

      if ( j == "EU-SOUTH")
    legend <- theme(legend.position="bottom", legend.title = element_blank(), legend.text = element_text(size = 12))
 else
    legend <- NULL


  title <- annotate("text", label = j, x = min(Ds$AGE) + 10, y = max(Ds$RSBP))

    out <- gr + xl + yl + title + legend

 }  )



#

 contour.condit.sex_0 <- lapply(hubs, function(j){

     Ds <-  subset(vertical.ist, subset = COUNTRY == j & SEX == 0)
 gr <- contourplot( Ds, aes(x = AGE, y = RSBP), cumulative = TRUE , legend.command = ifelse(levels(factor(Ds$COUNTRY)) == "EU-SOUTH", TRUE, FALSE),
                   breaks = br)

     if ( j == "EU-EAST")
    xl <- xlab("AGE")
 else
 xl <- xlab(" ")

if ( j == "EU-WEST")
  yl <- ylab("RSBP")
 else
yl <- ylab(" ")

 if ( j == "EU-SOUTH")
    legend <- theme(legend.position="bottom", legend.title = element_blank(), legend.text = element_text(size = 12))
 else
    legend <- NULL

  title <- annotate("text", label = j, x = min(Ds$AGE) + 10, y = max(Ds$RSBP))

    out <- gr + xl + yl + title + legend

 }  )

#

allcontour.sex1 <- plot_grid(  contour.condit.sex_1[[1]],  contour.condit.sex_1[[2]],  contour.condit.sex_1[[3]],  contour.condit.sex_1[[4]] , ncol = 2 )

allcontour.sex0 <- plot_grid(  contour.condit.sex_0[[1]],  contour.condit.sex_0[[2]],  contour.condit.sex_0[[3]],  contour.condit.sex_0[[4]] , ncol = 2 )


rm( contour.condit.sex_1,  contour.condit.sex_0, allcontour.sex1, allcontour.sex0, vertical.ist)



############# SUPPLEMENTARY MATERIAL PART I ###########################
################################################################



#### IST REGIONAL SUMMARY TABULAR DESCRIPTION (SUPPLEMENTARY MATERIAL PART I)

      ist.similarity <- lapply(istsimul.by.region, function(i){  # fish info
    lapply(i, function(j){
  return(j$is.data.similar)
    })   })


###### GIVE SHORT DESCRIPTION OF 3rd and 4th moments of IST marginals (REGION SPECIFIC)


  istsummary <- lapply(ist.similarity[[1]], function(x){

      out <- data.frame(third = x$third.moment[, "ipd"], fourth = x$fourth.moment[, "ipd"])
      rownames(out) <- rownames(x$third.moment)

      list( highmoms = out,
           corr = x$lower.triangular.Rx[, "ipd"] )

        }); names(istsummary) <- hubs



##

istsimnorm <- lapply( ist.similarity, function(i){

    lapply(i, function(j){

     return( list("1st mom"  = j$first.moment,
          "2nd mom" = j$second.moment,
          "3rd mom" =  j$third.moment ,
          "4th mom" =  j$fourth.moment ,
          "corr" = j$lower.triangular.Rx
         )   )
    })
     })



###   # pool info

  istsimlong <- do.call("rbind",
  lapply(1:2, function(k){        # meth

      do.call("rbind",
   lapply( 1:length(hubs), function(j){     # country

          stats <- names(istsimnorm[[k]][[j]])
        data <- istsimnorm[[k]][[j]][stats]

      do.call("rbind",
   lapply(stats, function(i){      # stats
 subdata <- data[[i]]

       data.frame( var = rownames(subdata),
 bias = subdata$mc.mean - subdata$ipd, stat = i, meth = k + 2, country = hubs[j] )
  }) )

   }    )   )

     })  )


### relabel stats

istsimlong$stat <- factor(istsimlong$stat, labels = c("1st moment", "2nd moment", "3rd moment", "4th moment", "correlation"))
# levels(istsimlong$stat) <- c( "1st moment", "2nd moment", "3rd moment", "4th moment", "correlation")

names(istsimlong)[3] <- "statistic"


## succinct graphical diagnostic

 mp <- aes( x= meth, y = bias, shape = statistic)
pl <- ggplot( istsimlong, mp );  bxpl <- geom_point( size = 4, show.legend = FALSE )
  panels <- facet_wrap(  ~ country )
 xl <- xlab("Generative Method")
  yl <- ylab("Bias")
zero <- geom_hline(yintercept=0, linetype=2, colour="green")
up <- geom_hline(yintercept= 0.5, linetype=3, colour="red")
down <- geom_hline(yintercept=-0.5, linetype=3, colour="red")

 plot_two <- pl + bxpl + panels + xl + yl + theme_bw() + zero + up + down + scale_x_continuous(labels = c("", "NORTA-\u0393", "NORTA-J", ""), limits = c(2,5))


## EU-EAST 4th mom outlier under NORTA  (IMPORTANT !!!!)

  ##   ist.similarity[[4]][[3]]$fourth.moment  ## is Heparin assignement in EU-EAST that is well below 1%
## mean(ist$RHEP24[ist$COUNTRY=="EU-EAST"])


########## GRAPHICAL CHECK AGAINST ORIG. IPD SUMMARIES ############

md <- Return.IPD.design.matrix(ist[, -9])

 IPDsummary <- Return.key.IPD.summaries( md, "moment.corr" )

IPDsummary$johnson.parameters
istmoms <- IPDsummary$first.four.moments
istcorr <- IPDsummary$correlation.matrix[lower.tri(IPDsummary$correlation.matrix)]
sort(abs((round(istcorr, 4))))

istsumm <- list(istmoms[, 1], istmoms[, 2], istmoms[, 3], istmoms[, 4], istcorr )  # original pooled summaries
#
 istsimulsumm <- lapply(PoolArtifData(istsimul.by.region), function(j){

     lapply(j, function(i){

  data <- i[, -10]
    corr <- momcor(data)

    list(mx= apply(data, 2, mean, na.rm = T), sdx = apply(data, 2, sd, na.rm = T),
         skx = apply(data, 2, skewness, na.rm = T),
         ktx = apply(data, 2, kurtosis, na.rm = T), corr = corr[lower.tri(corr)])

     } )    })

#

 summmerg <- lapply(istsimulsumm, function(j){

     lapply( c("mx", "sdx", "skx", "ktx", "corr"), function(s){

         do.call("cbind",
   lapply(j, function(i){

       i[[s]]

     })  )     }   )                   })
#

aversumm <- lapply( summmerg, function(j){

    lapply(j, function(i) apply(i, 1, mean, na.rm = T))

   })

#

mergedsumm <- lapply(aversumm, function(j){

    lapply(1:5, function(i){

        data.frame(ipd = istsumm[[i]], mc.mean = j[[i]] )
    })        })



# redundant ....

istsummnorm <- lapply( mergedsumm, function(j){

     return( list("1st mom"  =  j[[1]],
          "2nd mom" =  j[[2]],
          "3rd mom" =  j[[3]] ,
          "4th mom" =  j[[4]] ,
          "corr" = j[[5]]
         )   )
    })



#

  istsummlong <- do.call("rbind",
  lapply(1:2, function(j){

      stats <- names(istsummnorm[[j]])
        data <- istsummnorm[[j]][stats]

      do.call("rbind",
   lapply(stats, function(i){
 subdata <- data[[i]]

       data.frame( var = rownames(subdata),
 bias = subdata$mc.mean - subdata$ipd, stat = i, meth = j + 2  )
  }) )

     })  )


# dispaly merged ist

 mp <- aes( x= meth, y = bias, shape = stat)
pl <- ggplot( istsummlong, mp );  bxpl <- geom_point( size = 4 )

 xl <- xlab("Generative Method")
  yl <- ylab("Bias")
zero <- geom_hline(yintercept=0, linetype=2, colour="green")
up <- geom_hline(yintercept=0.5, linetype=3, colour="red")
down <- geom_hline(yintercept=-0.5, linetype=3, colour="red")

 plot_three <- pl + bxpl + xl + yl + theme_bw() + zero + up + down + scale_x_continuous(labels = c("", "NORTA-\u0393", "NORTA-J", ""), limits = c(2,5)) +  labs(shape="Key IPD summaries")


### print FIGURE 5 of SUPPLEMENT PART I

png( file.path(getwd(),"fig1AB.png"), width = 700, height = 900)

plot_grid( plot_two, plot_three, ncol = 1, labels = c("A", "B"), rel_heights = c(2.5, 1.5) )

dev.off()


rm(istsimulsumm, summmerg)


### FIGURE 4 OF SUPPLEMENT PART I
### theoretical fit on original IPD

  john <- pick.density("johnson")
gamm <- pick.density("gamma")



pdf( file.path( getwd(),"SM_fig2.pdf"), width= 7, height = 7)

par(mfrow=c(3,3))

for(h in hubs){

hist(ist$AGE[ist$COUNTRY == h], freq=F, xlab = "AGE",
     main=paste0(h,"_AGE (", istsimul.by.region[[2]][[h]]$ipd.johnson.parameters[, "AGE"][1], ")")); curve(john(x, istsimul.by.region[[2]][[h]]$ipd.johnson.parameters[, "AGE"][-6]), add = T, col = 2, lwd = 2)
curve(gamm(x, t(istsimul.by.region[[2]][[h]]$ipd.moments)[, "AGE"]), add = T, lty = 2, lwd = 2)


hist(ist$RSBP[ist$COUNTRY == h], freq=F, xlab = "RSBP",
     main=paste0(h,"_RSBP (", istsimul.by.region[[2]][[h]]$ipd.johnson.parameters[, "RSBP"][1], ")")); curve(john(x, istsimul.by.region[[2]][[h]]$ipd.johnson.parameters[, "RSBP"][-6]), add = T, col = 2, lwd = 2)
curve(gamm(x, t(istsimul.by.region[[2]][[h]]$ipd.moments)[, "RSBP"]), add = T, lty = 2, lwd = 2)

  }

legend("topright", title = "marginal", c("J", "G"), col =c(2,1), lwd = c(2,2), lty=c(1,2), bty = "n")


dev.off()


#     rm(istsimul.by.region)


########### SCATTERPLOT WITH LINEAR CORRELATION


panel.fit <- function(x,y){
    b <- lm(y~x)$coef
    pred <- b[1] + b[2]*x
   points(x,y, pch = ".")
lines(lowess(x,y), col= 'blue', lty = 4, lwd = 2)
    lines(x, pred , col='red', lwd = 2)


}

#

panel.cor <- function(x,y){

    r <- cor(x,y)

    legend("center", paste(round(r,2)), bty = "n")

  }



## FIGURE 3 OF SUPPLEMENT PART I

dat2 <- as.data.frame(Return.IPD.design.matrix(ist[ist$COUNTRY=="EU-EAST", -9])[, -1])

pdf( file.path(getwd(),"scatter2C.pdf"), width= 9, height = 9)

pairs(dat2, upper.panel = panel.cor, lower.panel = panel.fit, main = "EU-EAST" )

dev.off()




# not shown figures

dat <- as.data.frame(Return.IPD.design.matrix(ist[ist$COUNTRY=="EU-WEST", -9])[, -1])

pdf( file.path(getwd(),"scatter2A.pdf"), width= 9, height = 9)

pairs(dat, upper.panel = panel.cor, lower.panel = panel.fit, main = "EU-WEST" )

dev.off()


dat1 <- as.data.frame(Return.IPD.design.matrix(ist[ist$COUNTRY=="EU-NORTH", -9])[, -1])

pdf( file.path(getwd(),"scatter2B.pdf"), width= 9, height = 9)

pairs(dat1, upper.panel = panel.cor, lower.panel = panel.fit, main = "EU-NORTH" )

dev.off()


dat3 <- as.data.frame(Return.IPD.design.matrix(ist[ist$COUNTRY=="EU-SOUTH", -9])[, -1])

pdf( file.path(getwd(),"scatter2D.pdf"), width= 9, height = 9)

pairs(dat3, upper.panel = panel.cor, lower.panel = panel.fit, main = "EU-SOUTH" )

dev.off()





#### IST POOLED SUMMARY TABULAR DESCRIPTION (SUPPLEMENTARY MATERIAL PART II)

### tables 1-10 of supplement part II

  stnames <- rownames(ist.similarity[[1]][[1]][[1]])

# build regional moments table

isttabs <-  lapply(hubs, function(h)

    do.call("rbind", lapply(stnames, function(j){     # variables

     meth.array <- do.call("rbind",
    lapply(ist.similarity, function(m){  # methods

           do.call("cbind", lapply(1:4, function(i) m[[h]][[i]][j, 2]
     )  )  # moments

    }
         )  )

  ipd.ref <- do.call("cbind",
   lapply(1:4, function(i) ist.similarity[[1]][[h]][[i]][j, 1]                    )  )   # ipd ref value (out of meth loop)

     main <- as.data.frame(rbind(ipd.ref, meth.array))  # main block

  out <- data.frame("Variable" = j, "Origin" = c("IPD", 3:4), main)
       colnames(out)[3:6] <- c("1st", "2nd", "3rd", "4th")
         return(out)
             }

      )  )  )



 ## highlight pm 0.5 departures from ipd value

 ## must round digits before ....


  istlatextabs <- lapply(isttabs, function(h){

      tabhl <- do.call("rbind", lapply(stnames, function(j)
         apply( h[h$Variable == j, 3:6], 2, function(x)
   ifelse(abs(x[1] - x) > 0.5, paste("\\bfseries", sprintf("%.3f", round(x, 3)), sep = " "), sprintf("%.3f", round(x, 3)))
         )  )  )


 out <- data.frame( h[1:2], tabhl)


  make.multirow( out, 1, rotate = T)


  }   )

   names(istlatextabs) <- hubs


# build regional correlations table


 cpnamesfull <- unlist(lapply( stnames, function(i) lapply(stnames, function(j){

        paste( i, "-", j, sep = "" )

    }  )) )

#
 cpnames <- cpnamesfull[as.vector(lower.tri(matrix( nrow = length(stnames), ncol=length(stnames))))]


#
 istlatextabs2 <-  lapply(hubs, function(h){

  corr.array <- do.call("cbind",
     lapply(ist.similarity, function(m){  # methods

      m[[h]][[5]][, 2] # correlation pairs

        }  ))

  ipd.corr <- ist.similarity[[1]][[h]][[5]][, 1]

# boldfaced departures pm 0.1

       corhl <- apply( cbind(ipd.corr, corr.array), 1, function(x)
   ifelse(abs(x[1] - x) > 0.1, paste("\\bfseries", sprintf("%.3f", round(x, 3)), sep = " "), sprintf("%.3f", round(x, 3)))  )

#
   data.frame( title = cpnames, t(corhl) )
      } )

  names(istlatextabs2) <- hubs



 # build pooled moments table

  tabS3 <- do.call("rbind", lapply(stnames, function(j){     # variables

     meth.array <- do.call("rbind", lapply(istsummnorm, function(m){  # methods

           do.call("cbind", lapply(1:4, function(i) m[[i]][j, 2]
     )  )  # moments

    }
         )  )

  ipd.ref <- do.call("cbind", lapply(1:4, function(i) istsummnorm[[1]][[i]][j, 1]                    )  )   # ipd ref value (out of meth loop)

     main <- as.data.frame(rbind(ipd.ref, meth.array))  # main block

  out <- data.frame("Variable" = j, "Origin" = c("IPD", 3:4), main)
       colnames(out)[3:6] <- c("1st", "2nd", "3rd", "4th")
         return(out)
             }

      )  )



## highlight pm 0.5 departures from ipd value

 ## must round digits before ....

 tabhl <- do.call("rbind", lapply(stnames, function(j)
         apply( tabS3[tabS3$Variable == j, 3:6], 2, function(x)
   ifelse(abs(x[1] - x) > 0.5, paste("\\bfseries", sprintf("%.3f", round(x, 3)), sep = " "), sprintf("%.3f", round(x, 3)))
         )  )  )




 tabS3 <- data.frame( tabS3[1:2], tabhl)


tabS3 <- make.multirow( tabS3, 1, rotate = T)


# build correlations table


  corr.array <- do.call("cbind", lapply(istsummnorm, function(m){  # methods

      m[[5]][, 2] # correlation pairs

        }  ))

  ipd.corr <- istsummnorm[[1]][[5]][, 1]


## highlight pm 0.1 departures from ipd value

 ## must round digits before ....


       corhl <- apply( cbind(ipd.corr, corr.array), 1, function(x)
   ifelse(abs(x[1] - x) > 0.1, paste("\\bfseries", sprintf("%.3f", round(x, 3)), sep = " "), sprintf("%.3f", round(x, 3)))  )


#
 tabS4 <- data.frame( title = cpnames, t(corhl) )



### OUT REPORT


 toreport <- list(ist = ist, istsummary = istsummary, tabletwo = table2, tablethree= table3, sdre = sdranefsd, qqsdre = qqranefsd)


saveRDS( toreport, file.path(getwd(), "toreport.rds"))

rm(toreport)


# SAVE SUPPLEMENTARY MATERIAL

suppmaterial <- list( istregiomoms =  istlatextabs, istregiocorr = istlatextabs2, istmoms = tabS3, istcorr = tabS4, countryORs = countryORs, ORs65 = ORs65, tab2B = table2B)


 saveRDS( suppmaterial, "suppmaterial.rds")

rm(suppmaterial)



### NOTE: contour plotting might cause memory-allocation errors. To re-run the script you might need to quit/restart R.
