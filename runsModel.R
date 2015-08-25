# Code to fit the Siler mortality model to re-sighting data 
# of male and female lions of two populations, one from Hwange National Park
# and one from Serengeti National Park.

# Date: 25 August 2015
# Authors: Fernando Colchero (colchero@imada.sdu.dk) and Julia Barthold
#          (bartholdja@gmail.com)

# Comments:

# Download the this file, the functions file fcts01.R, and the data datSH.txt and
# save them to the same folder.
# Set your working directory to that folder.
# The model output will be automatically stored as an Rdata file to that folder
# There will be three output objectsn named:
# outFrom"model start age"_"version".Rdata

# This r code and the used data will be deposited at www.dryad.org after 
# the paper is accepted for publication.

# set the working directory
rm(list = ls())
setwd("...")

# load libraries
library(msm)
library(RColorBrewer)

# version of the analysis
vers <- "01"

# this loop runs the model from the three different starting ages 0, 0.5 and 1 year of age.

for (i in 1:3) {
  source("fcts01.R")
  
  # starting age for model fitting
  ageFitVec <- c(0, 0.5, 1)
  ageFit <- ageFitVec[i]
  
  # data
  dat <- read.table("datSH.txt", header = TRUE)
  dat <- dat[dat$xli >= ageFit, ]  
  
  # number of observations
  ni <- nrow(dat)
  
  # ages
  minDx <- 1.5 - ageFit
  xl <- dat$xli - ageFit
  xt <- rep(NA, nrow(dat))
  xt[dat$xti != 0] <- dat$xti[dat$xti != 0] - ageFit
  xt[dat$xti == 0] <- dat$xti[dat$xti == 0]
  xt[which(xt < 0)] <- 0
  ageMat <- matrix(c(2.5 - ageFit, 3.5 - ageFit), 2, 1,
                   dimnames =  list(c("f", "m"), "age at maturity"))
  
  # indicators and indices
  # sexes (1: female, 0: male)
  probFemVec <- c(0.51, 0.53, 0.55)
  probFem <- probFemVec[i]
  idNoSexS <- which(dat$unsex == 1 & dat$loc == "s")
  idNoSexH <- which(dat$unsex == 1 & dat$loc == "h")
  sexSStart <- rbinom(length(idNoSexS), 1, probFem)
  sexHStart <- rbinom(length(idNoSexH), 1, probFem)
  
  # Sex and population covariate matrix (sf, sm, hf, hm)
  # sexes: m for males, f for females
  # locations: s for serengeti, h for Hwange, 
  ncovs <- 4
  covs <- matrix(0, ni, ncovs, dimnames = list(NULL, c("f:s", "m:s", "f:h", "m:h")))
  covs[(dat$sex == "f" | dat$sex == "u" | dat$sex == "x") & dat$loc == "s", 1] <- 1
  covs[dat$sex == "m" & dat$loc == "s", 2] <- 1
  
  covs[(dat$sex == "f" | dat$sex == "u") & dat$loc == "h", 3] <- 1
  covs[dat$sex == "m" & dat$loc == "h", 4] <- 1
  
  # set up the covariate matrix according to proposed sexes
  covsStart <- CalcCov(covs, sexSStart, sexHStart)
  
  # dispersal covariate matrix (s, h)
  # locations: s for serengeti, h for Hwange, 
  nDispCovs <- 2
  dCovs <- matrix(0, nrow(covs), 2)
  colnames(dCovs) <- c("s", "h")
  dCovs[ ,1] <- rowSums(covsStart[ ,1:2])
  dCovs[ ,2] <- rowSums(covsStart[ ,3:4])
  
  # immigrant indicator
  idIm <- which(dat$im == 1)
  # exclude all immigrants that have truncation ages smaller than 
  # minimum age at dispersal (truncation caused by other mechanism than natal
  # dispersals)
  idIm <- idIm[which((xt-minDx)[idIm] > 0 )]
  
  # potential natal dispersal emigrants indicator
  idPotEm <- which(dat$potEm == 1 & xl > minDx)
  # likely secondary dispersal emigrants in the Serengeti indictor
  idPotSecEm <- which(dat$im == 1 & xt > minDx & dat$cens == 0 & dat$dead == 0 &
                        dat$outMigr == 1)
  # potential disperser starting index
  dispStateStart <- rep(0, ni)
  dispStateStart[c(idPotEm, idPotSecEm)] <- rbinom(length(c(idPotEm, 
                                                            idPotSecEm)), 1, 0.5)
  # potential dispersers starting indicator
  idEmStart <- which(dispStateStart == 1)
  
  # censored individuals indicator
  idCens <- which(dat$cens == 1)
  # censored native-born males = future potential dispersers indicator
  idCensPotEm <- which(dat$sex == "m" & dat$cens == 1 & dat$im == 0)
  
  # Siler mortality parameters: starting values and priors
  thf <- matrix(c(-1.4, 0.65, 0.07, -3.8, 0.2), 1, 5, 
                dimnames = list("f", c("a0", "a1", "c", "b0", "b1")))
  thm <- matrix(c(-1.2, 0.7, 0.16, -3.5, 0.23), 1, 5, 
                dimnames = list("m", c("a0", "a1", "c", "b0", "b1")))
  thetaStart <- matrix(rep(c(thf, thm), 2), ncovs, 5, 
                       byrow = TRUE, 
                       dimnames = list(colnames(covs), 
                                       c("a0", "a1", "c", "b0", "b1")))
  nthe <- length(thetaStart)
  
  thetaPriorMean <- matrix(c(-3, 0.2, 0, -4, 0.01), ncovs, 5, byrow = TRUE, 
                           dimnames = dimnames(thetaStart))
  thetaPriorSd <- matrix(rep(c(0.5, 0.25, 0.25, 0.5, 0.25), each = ncovs), ncovs, 5,  
                         dimnames = dimnames(thetaStart))
  thetaLow <- thetaStart * 0
  thetaLow[, c(1, 4)] <- -Inf
  
  # dispersal pars: starting values and priors
  gammaPrior1 <- matrix(rep(c(8, 2), 2), 2, 2, byrow = T)
  gammaPrior2 <- matrix(rep(c(2, 1), 2), 2, 2, byrow = T)
  gammaStart <- matrix(rep(c(4, 0.1), 2), 2, 2, byrow = T, dimnames =
                         list(c("disp pars s", "disp pars h"), c("gam1.1", "gam1.2")))
  
  gamma2Prior1 <- matrix(rep(c(8, 2), 2), 2, 2, byrow = T)
  gamma2Prior2 <- matrix(rep(c(2, 1), 2), 2, 2, byrow = T)
  gamma2Start <- matrix(c(4, 0.1, 0, 0), 2, 2, byrow = T, dimnames =
                          list(c("sec disp pars s", "sec disp pars h"), c("gam2.1", "gam2.2")))
  
  
  # run model
  niter <- 15000
  burn <- 5000
  thin <- 20
  keep <- seq(burn, niter, thin)
  
  nsim <- 4
  ncpus <- 4
  require(snowfall)
  sfInit(parallel = TRUE, cpus = ncpus)
  sfExport(list = c(ls(), ".Random.seed"))
  sfLibrary(msm, warn.conflicts = FALSE)
  outParallel <- sfClusterApplyLB(1:nsim, RunMCMC)
  sfStop() 
  
  # Process output
  out <- ExtractParalOut(outParallel)
  quantList <- lapply(1:length(out), function(pop){
    CalcDemoQuant(out[[pop]])
  })
  names(quantList) <- c("serengeti", "hwange")
  
  spOut <- list(Serengeti = list(theta = out$serengeti$theta,
                                 coeffs = out$serengeti$coeffs,
                                 mort = quantList$serengeti$mort,
                                 surv = quantList$serengeti$surv,
                                 pdf = quantList$serengeti$pdf,
                                 EH = quantList$serengeti$EH,
                                 cuts99 = quantList$serengeti$cuts99,
                                 cuts95 = quantList$serengeti$cuts95,
                                 cutsDisp99 = quantList$serengeti$cutsDisp99,
                                 cutsDisp95 = quantList$serengeti$cutsDisp95,
                                 x = quantList$serengeti$x,
                                 disp = quantList$serengeti$disp,
                                 disp2 = quantList$serengeti$disp2
  ),
  Hwange = list(theta = out$hwange$theta,
                coeffs = out$hwange$coeffs,
                mort = quantList$hwange$mort,
                surv = quantList$hwange$surv,
                pdf = quantList$hwange$pdf,
                EH = quantList$hwange$EH,
                cuts99 = quantList$hwange$cuts99,
                cuts95 = quantList$hwange$cuts95,
                cutsDisp99 = quantList$hwange$cutsDisp99,
                cutsDisp95 = quantList$hwange$cutsDisp95,
                x = quantList$hwange$x,
                disp = quantList$hwange$disp
  ))
  
  # Calculating the Kuhlback-Leibler divergences:
  # Calculate the KLc for comparison between sexes within locations
  # Serengeti females/males
  sFM <- lapply(1:5, function(i) {
    CalcKLc(spOut$Serengeti$theta[ , i], spOut$Serengeti$theta[ , (i + 5)], 
            low = rep(thetaLow[1, i], 2))
  })
  sFMmqKl <- do.call('rbind', lapply(1:length(sFM), function(x){
    sFM[[x]]$mqKl
  }))
  dimnames(sFMmqKl) <- list(dimnames(thetaStart)[[2]], "sFM")
  # Hwange females/males
  hFM <- lapply(1:5, function(i) {
    CalcKLc(spOut$Hwange$theta[ , i], spOut$Hwange$theta[ , (i + 5)], 
            low = rep(thetaLow[1, i], 2))
  })
  hFMmqKl <- do.call('rbind', lapply(1:length(hFM), function(x){
    hFM[[x]]$mqKl
  }))
  colnames(hFMmqKl) <- "hFM"
  
  # same sex between locations
  # Females Serengeti/Hwange
  fSH <- lapply(1:5, function(i) {
    CalcKLc(spOut$Serengeti$theta[ , i], spOut$Hwange$theta[ , i], 
            low = rep(thetaLow[1, i], 2))
  })
  fSHmqKl <- do.call('rbind', lapply(1:length(hFM), function(x){
    fSH[[x]]$mqKl
  }))
  colnames(fSHmqKl) <- "fSH"
  
  mSH <- lapply(1:5, function(i) {
    CalcKLc(spOut$Serengeti$theta[ , i+5], spOut$Hwange$theta[ , i+5], 
            low = rep(thetaLow[1, i], 2))
  })
  mSHmqKl <- do.call('rbind', lapply(1:length(hFM), function(x){
    mSH[[x]]$mqKl
  }))
  colnames(mSHmqKl) <- "mSH"
  
  mqKl <- round(cbind(sFMmqKl, hFMmqKl, fSHmqKl, mSHmqKl), 3)
  rm(list = setdiff(ls(), c("spOut", "niter", "burn", "thin", "keep", "ageFit", "outParallel", "mqKl", "vers")))
  #Sys.time()
  # save output
  save.image(as.character(paste("outFrom", ageFit, "_", vers, ".Rdata",sep = "")))
  rm(list = setdiff(ls(), "vers"))
}  
system("say Ich bin fertig")






