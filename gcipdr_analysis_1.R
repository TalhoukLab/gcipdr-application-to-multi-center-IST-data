## gcipdr_analysis_1.R contains R commands to reproduce a particular analysis.
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


library(devtools)

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
library(cowplot)

library(xtable)


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


world.macroregions <- c( "AMER-SUD", "SOUTH-PAC", "EU", "EU", "AMER-SUD", "AMER", "AMER-SUD", "EU-OST", "EU-NORD", "EU-NORD", "EU-NORD", "EU", "EU-SUD", "ASIA-SUD", "EU-OST", "INDIA", "MID-EST", "EU-SUD", "ASIA", "EU", "SOUTH-PAC", "EU-NORD", "EU-OST", "EU-SUD", "EU-OST", "ASIA-SUD", "EU-OST", "EU-OST", "ASIA", "EU-SUD", "ASIA-SUD", "EU-NORD", "EU", "EU-SUD", "EU-NORD", "AMER" ) 


levels(ist$COUNTRY) <- world.macroregions


table(ist$COUNTRY)


## rename levels

levels(ist$RCONSC) <- c( 1, 0, 2 )

levels(ist$SEX) <- c( 0, 1)

levels(ist$RATRIAL) <- c(NA, 0, 1) 

    levels(ist$RXASP) <- c(0, 1)

     levels(ist$RXHEP) <- c(1, 1, 1, 0) 
ist$RXHEP <- relevel(ist$RXHEP, ref = "0")

      levels(ist$FDEAD) <- c(NA, 0, NA, 1) 
  

 ist$is.EU <- grepl("EU", ist$COUNTRY)


## select only EU countries

 ist <- ist[ist$is.EU, ]; ist <- ist[ , -c( 8,9, 11, 13 )]

  ist[ ,-c(1,9)] <- apply(ist[ ,-c(1,9)], 2, as.integer)  ## convert to numeric

 ist$COUNTRY <- as.factor(as.character(ist$COUNTRY))

ist$RCONSC <- relevel(ist$RCONSC, ref = "0")  #### IMPORTANT HERE !!! relevel RCONSC ..


 ist <- ist[ which(!is.na(ist$RATRIAL)),  ]

head(ist)

hubs <- levels(ist$COUNTRY)

#############################################################
######### GAUSSIAN COPULA ARTIFICIAL DATA GENERATION ########
#############################################################


 seed <- 49632

#
 lapply(hubs, function(j){  
 
  data <- ist[ist$COUNTRY == j, -9]
     
     lapply( 3:4, # method coding
       function(i){
   jiseed <- as.integer(paste(which(hubs == j), seed, i, sep = ""))
 set.seed( jiseed, "L'Ecuyer") # delete this line to assess stability
           print( system.time(
  simulation.list <- Simulate.many.datasets(list(data), H = NULL, i,
      checkdata = TRUE, tabulate.similar.data = TRUE )
  ))
  Robject.path <- paste( paste( file.path( getwd(), paste(j, "ist.IPDstar", sep = "_") ), i, sep = "_"), ".rds", sep="" )
   saveRDS(simulation.list, Robject.path)
       rm(simulation.list)          }  )       
      }) 
   

### NOTE: One first computes the input IPD summaries and then feed them to 'DataRebuild'. 'Simulate.many.datasets' is just a handy wrapper for 'DataRebuild' that executes these two steps automatically (read 'gcipdr' documentation)

### NOTE: method 3 and 4 is a code for NORTA-$\Gamma$ and NORTA-J respectively



## ## ###  merge country and method specific simulations

## ## first read data: for each method put country-spec data simulations in one bag 

istsimullist <- lapply(3:4, function(i){  # 

         lapply(hubs, function(j){  ##            
      
       Robject.path <- paste( paste( file.path( getwd(), paste(j, "ist.IPDstar", sep = "_") ), i, sep = "_"), ".rds", sep="" )

       datamenge <- readRDS(Robject.path)

             datamenge[[1]]$similar.data

   } )
    } )


## # access method, row-bind by country each simulation, and save

  lapply(3:4, function(i){  # by meth

   istsimul <- lapply(1:100, function(h){      # by simulation

         merged <- do.call("rbind",
         lapply(1:length(hubs), function(j){  ##  bind by country

           data <- as.data.frame(istsimullist[[(i-2)]][[j]][[h]])  # matrix simul
             data$COUNTRY <- hubs[j]
             return(data) 

             } )  )
       
        return(merged)
       })
      
       final.path <- paste( paste( file.path( getwd(), "ist.ISTstar"), i, sep = "_"), ".rds", sep="" )

      saveRDS(istsimul, final.path)
      rm(istsimul)
    })



 rm(istsimullist)




############ ARTIFICIAL DATA COMPATIBILITY DIAGNOSTICS #########


### check features of pooled artificial data 
## NOTE: THERE ARE NAs in IST so overall sample size is smaller in simulations.

folder <- lapply(3:4, function(i)
        paste( paste( file.path( getwd(), "ist.ISTstar"), i, sep = "_"), ".rds", sep="" )        )


   istsimul <- lapply(folder, function(i){ ## pooled simulated data objects 

     data <- readRDS(i)
  return(data)
   }); names(istsimul) <- paste0("meth_", 3:4)


 

##################### GRAPHICAL CHECK AGAINST ORIG. IPD ################

                               
# summary(apply(ist[ist$COUNTRY == "EU", ], 1, function(x) any(is.na(x))))  

 vertical.istsimul <- do.call( "rbind", lapply(c("meth_3","meth_4"),
       function(j) data.frame( do.call( "rbind", lapply(istsimul[[j]], function(x) x) ),
                              DATA = j     ) ) )

 vertical.ist <- rbind( vertical.istsimul,
    na.omit( data.frame(Return.IPD.design.matrix(ist[, -9], fill = T)[, -1], COUNTRY = ist$COUNTRY, DATA = "Orig. IPD") ) ) 

 levels( vertical.ist$DATA)[1:2] <- c("NORTA-\u0393", "NORTA-J")




 mp <- aes( x= AGE, y=..density.., fill = DATA, color = DATA)
pl <- ggplot( vertical.ist[, c("AGE","COUNTRY","DATA")], mp )
# hist <- geom_histogram(alpha=0.2, position="identity")
hist <- geom_density(alpha=0.2, position="identity", show.legend = FALSE)

panel <- facet_wrap(~ as.factor(COUNTRY), ncol=2)
 xl <- xlab("AGE")

 hist_age <- pl + hist + panel + xl + theme_bw() 

#

mp <- aes( x= RSBP, y=..density.., fill = DATA, color = DATA)
pl <- ggplot( vertical.ist[, c("RSBP","COUNTRY","DATA")], mp )
hist <- geom_density(alpha=0.2, position="identity")

panel <- facet_wrap(~ as.factor(COUNTRY), ncol=2)
 xl <- xlab("RSBP")

 hist_rsbp <- pl + hist + panel + xl + theme_bw() + theme(legend.position="bottom" ) 



hist_ist <- plot_grid( hist_age, hist_rsbp, ncol = 1, rel_heights = c(2.5, 2.7) )



### FIGURE 1 OF MAIN MANUSCRIPT

png("hist_ist.png", width = 550, height = 650)

hist_ist 

dev.off()


rm(hist_ist , hist_age, hist_rsbp,  vertical.istsimul)

### 2D CONTOURS


library(reshape2)


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
 gr <- contourplot( Ds, aes(x = AGE, y = RSBP), cumulative = TRUE , legend.command = ifelse(levels(factor(Ds$COUNTRY)) == "EU-SUD", TRUE, FALSE),
                   breaks = br) 

     if ( j == "EU-OST")
    xl <- xlab("AGE")
 else
 xl <- xlab(" ")
     
if ( j == "EU")
  yl <- ylab("RSBP")
 else      
yl <- ylab(" ")

 if ( j == "EU-SUD")
    legend <- theme(legend.position="bottom")
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
 gr <- contourplot( Ds, aes(x = AGE, y = RSBP), cumulative = TRUE , legend.command = ifelse(levels(factor(Ds$COUNTRY)) == "EU-SUD", TRUE, FALSE),
                   breaks = br) 

     if ( j == "EU-OST")
    xl <- xlab("AGE")
 else
 xl <- xlab(" ")
     
if ( j == "EU")
  yl <- ylab("RSBP")
 else      
yl <- ylab(" ")

      if ( j == "EU-SUD")
    legend <- theme(legend.position="bottom")
 else
    legend <- NULL

     
  title <- annotate("text", label = j, x = min(Ds$AGE) + 10, y = max(Ds$RSBP))

    out <- gr + xl + yl + title + legend
     
 }  )


 allcontour.fdead1 <- plot_grid(  contour.condit.dead_1[[1]],  contour.condit.dead_1[[2]],  contour.condit.dead_1[[3]],  contour.condit.dead_1[[4]] , ncol = 2 )

#

 contour.condit.dead_0 <- lapply(hubs, function(j){

     Ds <-  subset(vertical.ist, subset = COUNTRY == j & FDEAD == 0)
 gr <- contourplot( Ds, aes(x = AGE, y = RSBP), cumulative = TRUE , legend.command = ifelse(levels(factor(Ds$COUNTRY)) == "EU-SUD", TRUE, FALSE),
                   breaks = br) 

     if ( j == "EU-OST")
    xl <- xlab("AGE")
 else
 xl <- xlab(" ")
     
if ( j == "EU")
  yl <- ylab("RSBP")
 else      
yl <- ylab(" ")

 if ( j == "EU-SUD")
    legend <- theme(legend.position="bottom")
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
 gr <- contourplot( Ds, aes(x = AGE, y = RSBP), cumulative = TRUE , legend.command = ifelse(levels(factor(Ds$COUNTRY)) == "EU-SUD", TRUE, FALSE),
                   breaks = br) 

     if ( j == "EU-OST")
    xl <- xlab("AGE")
 else
 xl <- xlab(" ")
     
if ( j == "EU")
  yl <- ylab("RSBP")
 else      
yl <- ylab(" ")

      if ( j == "EU-SUD")
    legend <- theme(legend.position="bottom")
 else
    legend <- NULL

     
  title <- annotate("text", label = j, x = min(Ds$AGE) + 10, y = max(Ds$RSBP))

    out <- gr + xl + yl + title + legend
     
 }  )


 
#

 contour.condit.sex_0 <- lapply(hubs, function(j){

     Ds <-  subset(vertical.ist, subset = COUNTRY == j & SEX == 0)
 gr <- contourplot( Ds, aes(x = AGE, y = RSBP), cumulative = TRUE , legend.command = ifelse(levels(factor(Ds$COUNTRY)) == "EU-SUD", TRUE, FALSE),
                   breaks = br) 

     if ( j == "EU-OST")
    xl <- xlab("AGE")
 else
 xl <- xlab(" ")
     
if ( j == "EU")
  yl <- ylab("RSBP")
 else      
yl <- ylab(" ")

 if ( j == "EU-SUD")
    legend <- theme(legend.position="bottom")
 else
    legend <- NULL
           
  title <- annotate("text", label = j, x = min(Ds$AGE) + 10, y = max(Ds$RSBP))

    out <- gr + xl + yl + title + legend
     
 }  )

#

allcontour.sex1 <- plot_grid(  contour.condit.sex_1[[1]],  contour.condit.sex_1[[2]],  contour.condit.sex_1[[3]],  contour.condit.sex_1[[4]] , ncol = 2 )

allcontour.sex0 <- plot_grid(  contour.condit.sex_0[[1]],  contour.condit.sex_0[[2]],  contour.condit.sex_0[[3]],  contour.condit.sex_0[[4]] , ncol = 2 )


rm( contour.condit.sex_1,  contour.condit.sex_0, allcontour.sex1, allcontour.sex0, vertical.ist)

##############################################
 ########## IPD INFERENCE RECOVERY  ##########
##############################################

### GLM RANDOM EFFECT (IST)

    
 library("lme4")


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

 
  istmeth3 <- mclapply( istsimul[["meth_3"]], function(x) stat(x) )


  istmeth4 <- mclapply( istsimul[["meth_4"]], function(x) stat(x) )

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



############# SUPPLEMENTARY MATERIAL PART I ###########################
################################################################



#### IST REGIONAL SUMMARY TABULAR DESCRIPTION (SUPPLEMENTARY MATERIAL PART I)

############## region-specific data similarity (IST)  

 folder <- lapply(3:4, function(i){  # fish info

   lapply(hubs, function(j){
       paste( paste( file.path( getwd(), paste(j, "ist.IPDstar", sep = "_") ), i, sep = "_"), ".rds", sep="" )
    
    })
   
 })   


#
 
  istsimul.by.region <- lapply(folder, function(i){  # fish info
    lapply(i, function(j){
      data <- readRDS(j)
  return(data[[1]])  
    })   }); names(istsimul.by.region) <- paste0("meth_", 3:4)


  for(i in 1:2) names(istsimul.by.region[[i]]) <- hubs



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

levels(istsimlong$stat) <- c( "1st moment", "2nd moment", "3rd moment", "4th moment", "correlation")

names(istsimlong)[3] <- "statistic"


## succinct graphical diagnostic

 mp <- aes( x= meth, y = bias, shape = statistic)
pl <- ggplot( istsimlong, mp );  bxpl <- geom_point( show.legend = FALSE )
  panels <- facet_wrap(  ~ country )
 xl <- xlab("Simulation Method")
  yl <- ylab("Bias")
zero <- geom_hline(yintercept=0, linetype=2, colour="green")
up <- geom_hline(yintercept= 0.5, linetype=3, colour="red")
down <- geom_hline(yintercept=-0.5, linetype=3, colour="red")

 plot_two <- pl + bxpl + panels + xl + yl + theme_bw() + zero + up + down + scale_x_continuous(labels = c("", "NORTA-\u0393", "NORTA-J", ""), limits = c(2,5))  


## EU-OST 4th mom outlier under NORTA  (IMPORTANT !!!!)

  ##   ist.similarity[[4]][[3]]$fourth.moment  ## is Heparin assignement in EU-OST that is well below 1%
## mean(ist$RHEP24[ist$COUNTRY=="EU-OST"])


########## GRAPHICAL CHECK AGAINST ORIG. IPD SUMMARIES ############

md <- Return.IPD.design.matrix(ist[, -9])

 IPDsummary <- Return.key.IPD.summaries( md, "moment.corr" )

IPDsummary$johnson.parameters
istmoms <- IPDsummary$first.four.moments
istcorr <- IPDsummary$correlation.matrix[lower.tri(IPDsummary$correlation.matrix)]
sort(abs((round(istcorr, 4))))

istsumm <- list(istmoms[, 1], istmoms[, 2], istmoms[, 3], istmoms[, 4], istcorr )  # original pooled summaries
#
 istsimulsumm <- lapply(istsimul, function(j){

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
pl <- ggplot( istsummlong, mp );  bxpl <- geom_point( )

 xl <- xlab("Simulation Method")
  yl <- ylab("Bias")
zero <- geom_hline(yintercept=0, linetype=2, colour="green")
up <- geom_hline(yintercept=0.5, linetype=3, colour="red")
down <- geom_hline(yintercept=-0.5, linetype=3, colour="red")

 plot_three <- pl + bxpl + xl + yl + theme_bw() + zero + up + down + scale_x_continuous(labels = c("", "NORTA-\u0393", "NORTA-J", ""), limits = c(2,5))  



### print FIGURE 5 of SUPPLEMENT PART I

pdf( file.path(getwd(),"fig1AB.pdf"), width= 7, height = 9)

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

hist(ist$AGE[ist$COUNTRY == h], freq=F,
     main=paste0(h,"_AGE (", istsimul.by.region[[2]][[h]]$ipd.johnson.parameters[, "AGE"][1], ")")); curve(john(x, istsimul.by.region[[2]][[h]]$ipd.johnson.parameters[, "AGE"][-6]), add = T, col = 2, lwd = 2)
curve(gamm(x, t(istsimul.by.region[[2]][[h]]$ipd.moments)[, "AGE"]), add = T, lty = 2, lwd = 2)


hist(ist$RSBP[ist$COUNTRY == h], freq=F,
     main=paste0(h,"_RSBP (", istsimul.by.region[[2]][[h]]$ipd.johnson.parameters[, "RSBP"][1], ")")); curve(john(x, istsimul.by.region[[2]][[h]]$ipd.johnson.parameters[, "RSBP"][-6]), add = T, col = 2, lwd = 2)
curve(gamm(x, t(istsimul.by.region[[2]][[h]]$ipd.moments)[, "RSBP"]), add = T, lty = 2, lwd = 2)

  }

legend("topright", title = "marginal", c("J", "G"), col =c(2,1), lwd = c(2,2), lty=c(1,2), bty = "n")


dev.off()


rm(istsimul.by.region)


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

dat2 <- as.data.frame(Return.IPD.design.matrix(ist[ist$COUNTRY=="EU-OST", -9])[, -1])

pdf( file.path(getwd(),"scatter2C.pdf"), width= 9, height = 9)

pairs(dat2, upper.panel = panel.cor, lower.panel = panel.fit, main = "EU-OST" )

dev.off()




# not shown figures
 
dat <- as.data.frame(Return.IPD.design.matrix(ist[ist$COUNTRY=="EU", -9])[, -1])

pdf( file.path(getwd(),"scatter2A.pdf"), width= 9, height = 9)

pairs(dat, upper.panel = panel.cor, lower.panel = panel.fit, main = "EU" )

dev.off()


dat1 <- as.data.frame(Return.IPD.design.matrix(ist[ist$COUNTRY=="EU-NORD", -9])[, -1])

pdf( file.path(getwd(),"scatter2B.pdf"), width= 9, height = 9)

pairs(dat1, upper.panel = panel.cor, lower.panel = panel.fit, main = "EU-NORD" )

dev.off()


dat3 <- as.data.frame(Return.IPD.design.matrix(ist[ist$COUNTRY=="EU-SUD", -9])[, -1])

pdf( file.path(getwd(),"scatter2D.pdf"), width= 9, height = 9)

pairs(dat3, upper.panel = panel.cor, lower.panel = panel.fit, main = "EU-SUD" )

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

 toreport <- list(ist = ist, istsummary = istsummary, tabletwo = table2, tablethree= table3, sdre = sdranefsd)



saveRDS( toreport, file.path(getwd(), "toreport.rds"))

rm(toreport)


# SAVE SUPPLEMENTARY MATERIAL

suppmaterial <- list( istregiomoms =  istlatextabs, istregiocorr = istlatextabs2, istmoms = tabS3, istcorr = tabS4)


 saveRDS( suppmaterial, "suppmaterial.rds")

rm(suppmaterial)
