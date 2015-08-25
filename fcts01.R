# Version: 04
# Author: Julia Barthold
# Date: 04 Jun 2015
# Comment: 
# x: last seen age, xt: truncation age, xd: age at dispersal, xar: age at arrival (immigration)# x: last seen age, xt: truncation age, xd: age at dispersal, xar: age at arrival (immigration)
# allows secondary dispersal
# - difference to 03: only one postNow that alters the values for second and natal
#   dispersers through indices
# - difference to 05: fits both populations
#
{

### Functions to compute the Siler mortality hazard and surival:
  
# Siler mortality rates at age x
CalcMort <- function(th, ...) UseMethod("CalcMort")

CalcMort.matrix <- function(th, x) {
  exp(th[, 1] - th[, 2] * x) + th[, 3] + exp(th[, 4] + th[, 5] * x)
}

CalcMort.numeric <- function(th, x) {
  exp(th[1] - th[2] * x) + th[3] + exp(th[4] + th[5] * x)
}

# Siler survival from birth to age x
CalcSurv <- function(th, ...) UseMethod("CalcSurv")

CalcSurv.matrix <- function(th, x) {
  exp(exp(th[, 1])/th[, 2] * (exp(-th[, 2] * x) - 1) - th[, 3] * x + 
        exp(th[, 4])/th[, 5] * (1 - exp(th[, 5] * x)))
}

CalcSurv.numeric <- function(th, x) {
  exp(exp(th[1])/th[2] * (exp(-th[2] * x) - 1) - th[3] * x + 
        exp(th[4])/th[5] * (1 - exp(th[5] * x)))
}

### Functions to set up covariate matrices:

# function to populate covariate matrices with parameters
CalcCovPars <- function(pars, covs) {
  covs %*% pars
}

# function to update covariate matrix
CalcCov <- function(covs, sexS, sexH) {
  # arguments: covariate matrix, proposed sex Serengeti, proposed Sex Hwange
  covs[idNoSexS, 1] <- sexS
  covs[idNoSexS, 2] <- abs(sexS - 1)
  covs[idNoSexH, 3] <- sexH
  covs[idNoSexH, 4] <- abs(sexH - 1)
  return(covs)
}

### Functions to sample the Siler mortality parameters (theta):

# function to compute the mortality PDF of age at death
CalcMortPdf <- function(thCov, x) {
  log(CalcMort(thCov, x) * CalcSurv(thCov, x))
}

# function to compute the posterior for the mortality parameters
CalcMortPost <- function(th, thCov, mortPdf, idEm) {
  mortLike <- mortPdf
  mortLike[idCens] <- log(CalcSurv(thCov[idCens, ], xl[idCens]))
  mortLike[idEm] <- log(CalcSurv(thCov[idEm, ], xl[idEm]))
  mortLike <- mortLike - log(CalcSurv(thCov, xt))
  mortPost <- sum(mortLike) + sum(dtnorm(c(th), c(thetaPriorMean), 
                                         c(thetaPriorSd),low = c(thetaLow), 
                                         log = TRUE))
  return(mortPost)
}

### Functions to sample the dispersal parameters for natal and 
### secondary dispersal (gamma1 and gamma2) and impute dispersal state:

# function to compute the natal dispersal likelihood
CalcDisp1Like <- function(gamCov, idEm, thCov) {
  dispLike <- xl * 0
  idNatEm <- idEm[which(idEm %in% idPotEm)]
  # Immigrants
  dispLike[idIm] <- dgamma(xt[idIm] - minDx, gamCov[idIm,1], gamCov[idIm,2], 
                           log = TRUE) - log(CalcSurv(thCov[idIm, ], xt[idIm]))
  # Non-dispersers
  idNotEm <- idPotEm[!(idPotEm %in% idEm)]
  dispLike[idNotEm] <- pgamma(xl[idNotEm] - minDx, gamCov[idNotEm,1], 
                              gamCov[idNotEm,2], lower.tail = FALSE, log = TRUE)
  # Dispersers
  dispLike[idNatEm] <- dgamma(xl[idNatEm] - minDx, gamCov[idNatEm, 1],
                              gamCov[idNatEm, 2], log = TRUE)
  # Censored native born males
  dispLike[idCensPotEm] <-  pgamma(xl[idCensPotEm] - minDx,
                                   gamCov[idCensPotEm,1], gamCov[idCensPotEm,2],
                                   lower.tail = FALSE, log = TRUE)
  return(dispLike)
}

# function to compute the secondary dispersal likelihood
CalcDisp2Like <- function(gam2Cov, idEm, thCov) {
  disp2Like <- xl * 0
  idSecEm <- idEm[which(idEm %in% idPotSecEm)]
  # Non out-migrators
  idNotEm <- idPotSecEm[!(idPotSecEm %in% idEm)]
  disp2Like[idNotEm] <- pgamma(xl[idNotEm] - minDx, gam2Cov[idNotEm,1],
                               gam2Cov[idNotEm,2], 
                               lower.tail = FALSE, log = TRUE) - 
                                log(CalcSurv(thCov[idNotEm, ], xt[idNotEm]))
  # Out-migrators
  disp2Like[idSecEm] <- dgamma(xl[idSecEm] - minDx, gam2Cov[idSecEm, 1], 
                               gam2Cov[idSecEm, 2], 
                               log = TRUE) - 
                                log(CalcSurv(thCov[idSecEm, ], xt[idSecEm]))
  return(disp2Like)
}

# funciton to compute the natal dispersal gamma parameter posterior
CalcGamPost <- function(gam, dispLike) {
  sum(dispLike) + 
    sum(dtnorm(gam, gammaPrior1, gammaPrior2, log = TRUE, 
               low = 0))  
}

# function to compute the secondary dispersal gamma parameter posterior
CalcGam2Post <- function(gam2, disp2Like) {
  sum(disp2Like) + 
    sum(dtnorm(gam2, gamma2Prior1, gamma2Prior2, log = TRUE, 
               low = 0))  
}

# function to compute the prior for the dispersal state of potential natal 
# dispersers
CalcDisp1Prior <- function(idEm) {
  dispPr <- xl * 0
  idNatEm <- idEm[which(idEm %in% idPotEm)]
  dispPr[idNatEm] <- dgamma(xl[idNatEm] - minDx, gammaPrior1[1, 1], gammaPrior1[1, 2], 
                            log = TRUE)
  idNotEm <- idPotEm[!(idPotEm %in% idEm)]
  dispPr[idNotEm] <- pgamma(xl[idNotEm] - minDx, gammaPrior1[1, 1], gammaPrior1[1, 2], 
                            log = TRUE, lower.tail = FALSE)
  return(dispPr)
}

# function to compute the prior for the dispersal state of potential 
# secondary dispersers
CalcDisp2Prior <- function(idEm) {
  disp2Pr <- xl * 0
  idSecEm <- idEm[which(idEm %in% idPotSecEm)]
  disp2Pr[idSecEm] <- dgamma(xl[idSecEm] - minDx, gamma2Prior1[1, 1], gamma2Prior1[1, 2], 
                             log = TRUE)
  idNotEm <- idPotSecEm[!(idPotSecEm %in% idEm)]
  disp2Pr[idNotEm] <- pgamma(xl[idNotEm] - minDx, gamma2Prior1[1, 1], gamma2Prior1[1, 2], 
                             log = TRUE, lower.tail = FALSE)
  return(disp2Pr)
}

# function to compute the posterior for the dispersal state
# of potential natal dispersers
CalcDisp1Post <- function(thCov, mortPdf, idEm, disp1Like) {
  disp1Prior <- CalcDisp1Prior(idEm)
  disp1Post <- rep(0, ni)
  agePost <- rep(0, ni)
  
  # potential natal dispersers currently imputed to be natal dispersers
  idNatEm <- idEm[which(idEm %in% idPotEm)]
  # potential natal dispersers currently imputed to not be natal dispersers
  idNotEm <- idPotEm[!(idPotEm %in% idEm)]
  
  # dispersers survived to age xl and dispersed at age xl
  agePost[idNatEm] <- log(CalcSurv(thCov[idNatEm, ], xl[idNatEm])) + 
    log(CalcSurv(thetaPriorMean[1, ], xl[idNatEm]))
  disp1Post[idNatEm] <- disp1Like[idNatEm] + disp1Prior[idNatEm]
    
  # non-dispersers died at xl and would have disperserd at ages >xl.
  agePost[idNotEm] <- mortPdf[idNotEm] + 
    CalcMortPdf(thetaPriorMean[1, ], xl[idNotEm])
  disp1Post[idNotEm] <- disp1Like[idNotEm] + disp1Prior[idNotEm]
  return(disp1Post + agePost)
}
  
# function to compute the posterior for the dispersal state
# of potential secondary dispersers
CalcDisp2Post <- function(thCov, mortPdf, idEm, disp2Like) {
  disp2Prior <- CalcDisp2Prior(idEm)
  disp2Post <- rep(0, ni)
  agePost <- rep(0, ni)
  
  # potential secondary dispersers currently imputed to be secondary dispersers
  idSecEm <- idEm[which(idEm %in% idPotSecEm)]  
  # pot. secondary dispersers currently imputed to not be secondary dispersers
  idNotEm <- idPotSecEm[!(idPotSecEm %in% idEm)]
  
  # dispersers survived to age xl and dispersed at age xl
  agePost[idSecEm] <- log(CalcSurv(thCov[idSecEm, ], xl[idSecEm])) + 
    log(CalcSurv(thetaPriorMean[2, ], xl[idSecEm]))
  disp2Post[idSecEm] <- disp2Like[idSecEm] + disp2Prior[idSecEm]
  
  # non-dispersers died at xl and would have disperserd at ages >xl.
  agePost[idNotEm] <- mortPdf[idNotEm] + 
    CalcMortPdf(thetaPriorMean[2, ], xl[idNotEm])
  disp2Post[idNotEm] <- disp2Like[idNotEm] + disp2Prior[idNotEm]
  return(disp2Post + agePost)
}

### Functions to propose sex: 

# function to compute the prior for sex
CalcSexPrior <- function(covs) {
  (covs[, 1]) * probFem + (1 - covs[, 1]) * (1 - probFem)
}

# function to compute the posterior for sex
CalcSexPost <- function(covs, mortPdf) {
  mortPdf + log(CalcSexPrior(covs))
}

### Function to update jumps:
UpdateJumps <- function(jumps, updMat, iter, iterUpd, updTarg) {
  updRate <- apply(updMat[iter - ((iterUpd - 1):0), ], 2, sum) / iterUpd  
  updRate[updRate == 0] <- 1e-2
  jumps <- jumps * 
    matrix(updRate, nrow(jumps), ncol(jumps)) / updTarg
  return(jumps)
}
}

# MCMC
RunMCMC <- function(sim) {
  dispStateNow <- dispStateStart
  idEmNow <- idEmStart
  covsNow <- covsStart
  thNow <- thetaStart
  sexSNow <- sexSStart
  sexHNow <- sexHStart
  gamNow <- gammaStart
  gam2Now <- gamma2Start
  
  # start chains at random values for 3 of the 4 chains
  if (sim > 1) {
    rm(".Random.seed", envir = .GlobalEnv); runif(1)
    thNow[1:nthe] <- rtnorm(nthe, thetaStart[1:nthe], 0.1, 
                             low = thetaLow[1:nthe]) # adapt the low
  }

  if (sim > 1) {
    rm(".Random.seed", envir = .GlobalEnv); runif(1)
    gamNow <- matrix(rep(rtnorm(2, gammaStart[1, ], 0.1, low = 0), 2),
                     2, 2, byrow = T)
    gam2Now <- matrix(c(rtnorm(2, gamma2Start[1, ], 0.1, low = 0), rep(0, 2)),
                      2, 2, byrow = T)
  }

  gamCovNow <- CalcCovPars(gamNow, dCovs)
  gam2CovNow <- CalcCovPars(gam2Now, dCovs)
  thCovNow <- CalcCovPars(thNow, covsNow)
  mortPdfNow <- CalcMortPdf(thCovNow, xl)
  mortPostNow <- CalcMortPost(thNow, thCovNow, mortPdfNow, idEmNow)
  sexPostNow <- CalcSexPost(covsNow, mortPdfNow)
  disp1LikeNow <- CalcDisp1Like(gamCovNow, idEmNow, thCovNow)
  disp2LikeNow <- CalcDisp2Like(gam2CovNow, idEmNow, thCovNow)
  gamPostNow <- CalcGamPost(gamNow, disp1LikeNow)
  gam2PostNow <- CalcGam2Post(gam2Now, disp2LikeNow)
  disp1PostNow <- CalcDisp1Post(thCovNow, mortPdfNow, idEmNow, disp1LikeNow)
  disp2PostNow <- CalcDisp2Post(thCovNow, mortPdfNow, idEmNow, disp2LikeNow)
  # ojects to save outputs:
  namesThe <- paste(rep(colnames(thNow), ncovs), rep(rownames(thNow), each = 5), 
                    sep = '.')
  thMat <- matrix(0, niter, nthe, dimnames = list(NULL, namesThe))
  thUpdMat <- thMat
  thMat[1, ] <- c(t(thNow))
  gamMat <- matrix(0, niter, 4, dimnames = list(NULL, c("gam1.1:s", "gam1.2:s", "gam1.1:h", "gam1.2:h")))
  gamUpdMat <- gamMat
  gamMat[1, ] <- c(t(gamNow))
  gam2Mat <- matrix(0, niter, 2, dimnames = list(NULL, c("gam2.1:s", "gam2.2:s")))
  gam2UpdMat <- gam2Mat
  gam2Mat[1, ] <- c(t(gam2Now))[1:2]
  dispMat <- dispStateNow[idEmNow]
  sexMat <- matrix(c(sexSStart, sexHStart), 1, length(c(idNoSexS, idNoSexH)))
  
  # objects pertaining to updating the jumps
  theJumps <- matrix(rep(0.1, nthe), ncovs, 5, dimnames = dimnames(thNow))
  gamJumps <- matrix(0.01, 2, 2)
  gam2Jumps <- matrix(0.01, 1, 2)
  iterUpd <- 100
  updTarg <- 0.25
  UpdJumps <- TRUE
  
  # MCMC:
  Start <- Sys.time()
  for (iter in 2:niter) {
    
    # 1. Update mortality parameters:
    for (pp in 1:nthe) {
      thNew <- thNow
      thNew[pp] <- rtnorm(1, thNow[pp], theJumps[pp], low = thetaLow[pp])
      thCovNew <- CalcCovPars(thNew, covsNow)
      mortPdfNew <- CalcMortPdf(thCovNew, xl)
      mortPostNew <- CalcMortPost(thNew, thCovNew, mortPdfNew, idEmNow)
      r <- exp(mortPostNew - mortPostNow)
      if (!is.na(r) & r > runif(1)) {
        thNow <- thNew
        thCovNow <- thCovNew
        mortPdfNow <- mortPdfNew
        mortPostNow <- mortPostNew
        if (UpdJumps) thUpdMat[iter, pp] <- 1
      }
    }
    if (any(thUpdMat[iter, ] == 1)) {
      disp1LikeNow <- CalcDisp1Like(gamCovNow, idEmNow, thCovNow)
      disp2LikeNow <- CalcDisp2Like(gam2CovNow, idEmNow, thCovNow)
      gamPostNow <- CalcGamPost(gamNow, disp1LikeNow)
      gam2PostNow <- CalcGam2Post(gam2Now, disp2LikeNow)
      disp1PostNow <- CalcDisp1Post(thCovNow, mortPdfNow, idEmNow, disp1LikeNow)
      disp2PostNow <- CalcDisp2Post(thCovNow, mortPdfNow, idEmNow, disp2LikeNow)
    }
    
    # 2. Propose sex for unsexed individuals
    sexSNew <- rbinom(length(idNoSexS), 1, probFem)
    sexHNew <- rbinom(length(idNoSexH), 1, probFem)
    covsNew <-  CalcCov(covsNow, sexSNew, sexHNew)
    
    thCovNew <- CalcCovPars(thNow, covsNew)
    mortPdfNew <- CalcMortPdf(thCovNew, xl)
    mortPostNew <- CalcMortPost(thNow, thCovNew, mortPdfNew, idEmNow)
    sexPostNew <- CalcSexPost(thCovNew, mortPdfNew)
    
    r <- exp(sexPostNew - sexPostNow)[c(idNoSexS, idNoSexH)]
    
    idUpd2 <- c(idNoSexS, idNoSexH)[r > runif(length(c(idNoSexS, idNoSexH)))]
    idUpd2 <-idUpd2[!(is.na(idUpd2))]
    if (length(idUpd2) > 0 ) {
      mortPdfNow[idUpd2] <- mortPdfNew[idUpd2]
      covsNow[idUpd2, ] <- covsNew[idUpd2, ]
      thCovNow[idUpd2, ] <- thCovNew[idUpd2, ]
      sexPostNow[idUpd2] <- sexPostNew[idUpd2]
      sexNow <- c(covsNow[idNoSexS, 1], covsNow[idNoSexH, 3])     
    }
    
    mortPostNow <- CalcMortPost(thNow, thCovNow, mortPdfNow, idEmNow)
    disp1LikeNow <- CalcDisp1Like(gamCovNow,idEmNow, thCovNow)
    disp2LikeNow <- CalcDisp2Like(gam2CovNow,idEmNow, thCovNow)
    gamPostNow <- CalcGamPost(gamNow, disp1LikeNow)
    gam2PostNow <- CalcGam2Post(gam2Now, disp2LikeNow)
    disp1PostNow <- CalcDisp1Post(thCovNow, mortPdfNow, idEmNow, disp1LikeNow)
    disp2PostNow <- CalcDisp2Post(thCovNow, mortPdfNow, idEmNow, disp2LikeNow)
    
    # 3. update dispersal parameters:
    # a) natal
    for (pp in 1:4) {
      gamNew <- gamNow
      gamNew[pp] <- rtnorm(1, gamNow[pp], gamJumps[pp], lower = 0)
      gamCovNew <- CalcCovPars(gamNew, dCovs)
      disp1LikeNew <- CalcDisp1Like(gamCovNew,idEmNow, thCovNow)
      gamPostNew <- CalcGamPost(gamNew, disp1LikeNew)
      r <- exp(gamPostNew - gamPostNow)
      if (!is.na(r) & r > runif(1)) {
        gamNow <- gamNew
        gamCovNow <- gamCovNew
        disp1LikeNow <- disp1LikeNew
        gamPostNow <- gamPostNew
        if (UpdJumps) gamUpdMat[iter, pp] <- 1
      }
    }
    if (any(gamUpdMat[iter, ] == 1)) {
      disp1PostNow <- CalcDisp1Post(thCovNow, mortPdfNow, idEmNow, disp1LikeNow)
    }
    # b) secondary
    for (pp in 1:2) {
      gam2New <- gam2Now
      gam2New[c(1,3)[pp]] <- rtnorm(1, gam2Now[c(1,3)[pp]], gam2Jumps[pp], lower = 0)
      gam2CovNew <- CalcCovPars(gam2New, dCovs)
      disp2LikeNew <- CalcDisp2Like(gam2CovNew,idEmNow, thCovNow)
      gam2PostNew <- CalcGam2Post(gam2New, disp2LikeNew)
      r <- exp(gam2PostNew - gam2PostNow)
      if (!is.na(r) & r > runif(1)) {
        gam2Now <- gam2New
        gam2CovNow <- gam2CovNew
        disp2LikeNow <- disp2LikeNew
        gam2PostNow <- gam2PostNew
        if (UpdJumps) gam2UpdMat[iter, pp] <- 1
      }
    }
    if (any(gam2UpdMat[iter, ] == 1)) {
      disp2PostNow <- CalcDisp2Post(thCovNow, mortPdfNow, idEmNow, disp2LikeNow)
    }
    
    # 4. Update dispersal state, ages at death and ages at dispersal:
    # a) natal dispersal
    dispStateNew <- rep(0, ni)
    dispStateNew[c(idPotEm, idPotSecEm)] <- rbinom(length(c(idPotEm,
                                                          idPotSecEm)), 1, 0.5)
    idEmNew <- which(dispStateNew == 1)
    disp1LikeNew <- CalcDisp1Like(gamCovNow,idEmNew, thCovNow)
    disp1PostNew <- CalcDisp1Post(thCovNow, mortPdfNow, idEmNew, disp1LikeNew)
    r <- exp(disp1PostNew - disp1PostNow)
    idUpd <- which(r > runif(ni))
    idUpd <-idUpd[!(is.na(idUpd))] 
    if (length(idUpd) > 0) { 
      dispStateNow[idUpd] <- dispStateNew[idUpd]
      disp1LikeNow[idUpd] <- disp1LikeNew[idUpd]
      disp1PostNew[idUpd] <- disp1PostNew[idUpd]
    }
    gamPostNow <- CalcGamPost(gamNow, disp1LikeNow)
    idEmNow <- which(dispStateNow == 1)
    mortPostNow <- CalcMortPost(thNow, thCovNow, mortPdfNow, idEmNow)
    
    # b) secondary dispersal
    dispStateNew <- dispStateNow
    dispStateNew[idPotSecEm] <- rbinom(length(idPotSecEm), 1, 0.5)
    idEmNew <- which(dispStateNew == 1)
    disp2LikeNew <- CalcDisp2Like(gamCovNow,idEmNew, thCovNow)
    disp2PostNew <- CalcDisp2Post(thCovNow, mortPdfNow, idEmNew, disp2LikeNew)    
    r <- exp(disp2PostNew - disp2PostNow)
    z <- runif(ni)
    idUpd <- which(r > z)
    idUpd <-idUpd[!(is.na(idUpd))] 
    if (length(idUpd) > 0) { 
      dispStateNow[idUpd] <- dispStateNew[idUpd]
      disp2LikeNow[idUpd] <- disp1LikeNew[idUpd]
      disp2PostNew[idUpd] <- disp2PostNew[idUpd]
    }
    gam2PostNow <- CalcGam2Post(gam2Now, disp2LikeNow)
    idEmNow <- which(dispStateNow == 1)
    mortPostNow <- CalcMortPost(thNow, thCovNow, mortPdfNow, idEmNow)
    
    # 5. Dynamic Metropolis to update jumps:
    if (UpdJumps) {
      if (is.element(iter/iterUpd,c(1:10))) {
        theJumps <- UpdateJumps(theJumps, thUpdMat, iter, iterUpd, updTarg)
        gamJumps <- UpdateJumps(gamJumps, gamUpdMat, iter, iterUpd, updTarg)
        gam2Jumps <- UpdateJumps(gam2Jumps, gam2UpdMat, iter, iterUpd, updTarg)
      }
    }
    
    # 6. Fill in the output matrices:
    thMat[iter, ] <- c(t(thNow))
    gamMat[iter, ] <- c(t(gamNow))
    gam2Mat[iter, ] <- c(t(gam2Now))[1:2]
    if (iter %in% keep) {
      dispMat <- rbind(dispMat, dispStateNow[idPotEm])
      sexMat <- rbind(sexMat, sexNow)
    }
  }
  End <- Sys.time()
  out <- list(theta = thMat, gamma = gamMat, gamma2 = gam2Mat, disp = dispMat, 
              keep = keep, names = namesThe, potEm = idPotEm)
  return(out)
}

ExtractParalOut <- function(out) {
  # joint output objects
  thetaOut <- do.call('rbind', lapply(1:length(out), function(x){
    out[[x]]$theta[keep, ]
  }))
  gamOut <- do.call('rbind', lapply(1:length(out), function(x){
    out[[x]]$gamma[keep, ]
  }))
  gam2Out <- do.call('rbind', lapply(1:length(out), function(x){
    out[[x]]$gamma2[keep, ]
  })) 
  coeffs <- rbind(cbind(apply(thetaOut, 2, mean), apply(thetaOut, 2, sd), 
                        t(apply(thetaOut, 2, quantile, c(0.025, 0.975)))),
                  cbind(apply(gamOut, 2, mean), apply(gamOut, 2, sd), 
                        t(apply(gamOut, 2, quantile, c(0.025, 0.975)))),
                  cbind(apply(gam2Out, 2, mean), apply(gam2Out, 2, sd), 
                        t(apply(gam2Out, 2, quantile, c(0.025, 0.975)))))
  colnames(coeffs) <- c("Mean", "SE", "2.5%", "97.5%")
  parList <- list(serengeti = 
                    list(theta = thetaOut[ , grep("s", colnames(thetaOut))], 
                         coeffs = coeffs[grep("s", rownames(coeffs)), ],
                         pop = "serengeti", 
                         gam = gamOut[ , grep("s", colnames(gamOut))],
                         gam2 = gam2Out[ , grep("s", colnames(gam2Out))]),
                  hwange = 
                    list(theta = thetaOut[ , grep("h", colnames(thetaOut))], 
                         coeffs = coeffs[grep("h", rownames(coeffs)), ],
                         pop = "hwange", 
                         gam = gamOut[ , grep("h", colnames(gamOut))]))
  dimnames(parList$serengeti$theta) <- 
    list(NULL, sub(":s", "", 
                   colnames(thetaOut)[grep("s", 
                                           colnames(thetaOut))]))
  dimnames(parList$hwange$theta) <- dimnames(parList$serengeti$theta)
  rownames(parList$serengeti$coeffs) <- 
    sub(":s", "", rownames(coeffs)[grep("s", 
                                        rownames(coeffs))])
  rownames(parList$hwange$coeffs) <- rownames(parList$serengeti$coeffs)[1:12]
  return(parList)
}

# Functions to calculate life expectancy and entropy with CIs:
CalcHx <- function(Sx, dx) {
  Sx1 <- Sx[Sx > 0]
  -sum(Sx1 * log(Sx1) * dx) / sum(Sx1 * dx)
}

CalcEx <- function(Sx, dx) sum(Sx * dx)


# Function to calculate survival and mortality quantiles:
CalcDemoQuant <- function(out, ...) {
  
  # Ages
  dx <- 0.1
  xv <- seq(0, 100, dx)    
  # Mortality
  mortList <- lapply(c("f", "m"), function(ss) {
    ids <- grep(ss, colnames(out$theta))
    mort <- apply(out$theta[, ids], 1, function(th) CalcMort(th, xv))
    mortave <- apply(mort, 1, mean)
    mortci <- apply(mort, 1, quantile, c(0.025, 0.975))
    mortfin <- rbind(mortave, mortci)
    return(mortfin)
  })
  # Survival
  survList <- lapply(c("f", "m"), function(ss) {
    ids <- ids <- grep(ss, colnames(out$theta))
    surv <- apply(out$theta[, ids], 1, function(th) CalcSurv(th, xv))
    survave <- apply(surv, 1, mean)
    survci <- apply(surv, 1, quantile, c(0.025, 0.975))
    survfin <- rbind(survave, survci)
    return(survfin)
  })
  # pdf of ages at death
  pdfList <- lapply(c("f", "m"), function(ss) {
    ids <- ids <- grep(ss, colnames(out$theta))
    pdfm <- apply(out$theta[, ids], 1, 
                  function(th) CalcSurv(th, xv) * CalcMort(th, xv))
    pdfave <- apply(pdfm, 1, mean)
    pdfci <- apply(pdfm, 1, quantile, c(0.025, 0.975))
    pdffin <- rbind(pdfave, pdfci)
    return(pdffin)
  })
  # pdf of ages at natal out-migration
  dispList <- list()
  pdf <- apply(out$gam, 1, function(gam) {
    dgamma(xv, gam[1], gam[2])})
  pdfave <- apply(pdf, 1, mean)
  pdfci <- apply(pdf, 1, quantile, c(0.025, 0.975))
  dispList[[1]] <- rbind(pdfave, pdfci)
  # 1 - CDF to calculate cut off points for plotting
  OneMinusCDF <- apply(out$gam, 1, function(gam) {
    pgamma(xv, gam[1], gam[2], lower.tail = F)})
  OneMinusCDFave <- apply(OneMinusCDF, 1, mean)
  
  if(length(out$gam2 > 0)) {
    disp2List <- list()
    pdf <- apply(out$gam2, 1, function(gam2) {
      dgamma(xv, gam2[1], gam2[2])})
    pdfave <- apply(pdf, 1, mean)
    pdfci <- apply(pdf, 1, quantile, c(0.025, 0.975))
    disp2List[[1]] <- rbind(pdfave, pdfci)
  }
  
  EHqMat <- lapply(c("f", "m"), function(ss) {
    ids <- grep(ss, colnames(out$theta))
    xm <- ageMat[ss, 1] - ageFit
    xvm <- xv[xv >= xm]
    surv0 <- apply(out$the[, ids], 1, function(th) CalcSurv(th, xv))
    survM <- apply(out$the[, ids], 1, function(th) CalcSurv(th, xvm) / 
                     CalcSurv(th, xvm[1]))
    E0 <- apply(surv0, 2, CalcEx, dx = dx)
    EM <- apply(survM, 2, CalcEx, dx = dx)
    H0 <- apply(surv0, 2, CalcHx, dx = dx)
    HM <- apply(survM, 2, CalcHx, dx = dx)
    EHq <- rbind(c(mean(E0), quantile(E0, c(0.025, 0.975))),
                 c(mean(EM), quantile(EM, c(0.025, 0.975))),
                 c(mean(H0), quantile(H0, c(0.025, 0.975))),
                 c(mean(HM), quantile(HM, c(0.025, 0.975))))
    dimnames(EHq) <- list(c("E0", "EM", "H0", "HM"), 
                          c("Mean", "2.5%", "97.5%"))
    return(EHq)
  })
  
  plotCut99 <- list()
  for (i in 1:length(survList)) {
    plotCut99[[i]] <- which(survList[[i]][1, ] < 0.01)[1]
    if (is.na(plotCut99[[i]])) plotCut99[[i]] <- length(xv)
  }
  plotCut95 <- list()
  for (i in 1:length(survList)) {
    plotCut95[[i]] <- which(survList[[i]][1, ] < 0.05)[1]
    if (is.na(plotCut95[[i]])) plotCut95[[i]] <- length(xv)
  }
  
  plotCutDisp99 <- which(OneMinusCDFave < 0.01)[1]
  plotCutDisp95 <- which(OneMinusCDFave < 0.05)[1]
  
  names(EHqMat) <- names(survList) <- names(mortList) <- names(pdfList) <- 
    names(plotCut95) <- names(plotCut99) <- c("f", "m")
  if(length(out$gam2) > 0) {
    return(list(mort = mortList, surv = survList, pdf = pdfList, x = xv, 
                EH = EHqMat, cuts99 = plotCut99, cuts95 = plotCut95, 
                cutsDisp99 = plotCutDisp99, cutsDisp95 = plotCutDisp95,
                disp = dispList, disp2 = disp2List))
  } else {
    return(list(mort = mortList, surv = survList, pdf = pdfList, x = xv, 
                EH = EHqMat, cuts99 = plotCut99, cuts95 = plotCut95, 
                cutsDisp99 = plotCutDisp99, cutsDisp95 = plotCutDisp95,
                disp = dispList))
  }
}

# Calculate KLcs:
CalcKLc <- function(p1, p2, low) {
  require(msm)
  pMat <- cbind(p1, p2)
  parRan <- range(sapply(1:2, function(pp) qtnorm(c(0.001, 0.999), 
                                                  mean(pMat[, pp]), 
                                                  sd(pMat[, pp]), 
                                                  lower = low[pp])))
  parVec <- seq(parRan[1], parRan[2], length = 100)
  dp <- parVec[2] - parVec[1]
  parDens <- sapply(1:2, function(pp) 
    dtnorm(parVec, mean(pMat[, pp]), sd(pMat[, pp]), lower = low[pp]))
  # if one of the densities has 0s, restrict densities to the non-0 range
  kld1 <- sum(parDens[, 1] * log(parDens[, 1] / parDens[, 2]) * dp)
  kld2 <- sum(parDens[, 2] * log(parDens[, 2] / parDens[, 1]) * dp)
  qKlc1 <- (1 + (1 - exp(-2 * kld1)^(1 / 2))) / 2
  qKlc2 <- (1 + (1 - exp(-2 * kld2)^(1 / 2))) / 2
  mqKl <- (qKlc1 + qKlc2) / 2
  outList <- list(kl1 = kld1, kl2 = kld2, qkl1 = qKlc1, 
                  qkl2 = qKlc2, mqKl = mqKl)
  return(outList)
}



RunThis = F
if(RunThis) {
# function to extract thinned sequences from multiple 
# runs and calculate coefficients:
ExtractParalOut <- function(out) {
  # joint output objects
  thetaOut <- do.call('rbind', lapply(1:length(out), function(x){
    out[[x]]$theta[keep, ]
  }))
  gamOut <- do.call('rbind', lapply(1:length(out), function(x){
    out[[x]]$gamma[keep, ]
  }))
  gam2Out <- do.call('rbind', lapply(1:length(out), function(x){
    out[[x]]$gamma2[keep, ]
  }))
  coeffs <- rbind(cbind(apply(thetaOut, 2, mean), apply(thetaOut, 2, sd), 
                        t(apply(thetaOut, 2, quantile, c(0.025, 0.975)))),
                  cbind(apply(gamOut, 2, mean), apply(gamOut, 2, sd), 
                        t(apply(gamOut, 2, quantile, c(0.025, 0.975)))),
                  cbind(apply(gam2Out, 2, mean), apply(gam2Out, 2, sd), 
                        t(apply(gam2Out, 2, quantile, c(0.025, 0.975)))))
  colnames(coeffs) <- c("Mean", "SE", "2.5%", "97.5%")
  parList <- list(theta = thetaOut, coeffs = coeffs, gam = gamOut, gam2 = gam2Out)
  dimnames(parList$theta) <- list(NULL, colnames(thetaOut))
  return(parList)
}



# Functions to calculate life expectancy and entropy with CIs:
CalcHx <- function(Sx, dx) {
  Sx1 <- Sx[Sx > 0]
  -sum(Sx1 * log(Sx1) * dx) / sum(Sx1 * dx)
}

CalcEx <- function(Sx, dx) sum(Sx * dx)


# Function to calculate survival and mortality quantiles:
CalcDemoQuant <- function(out, ...) {
  
  # Ages
  dx <- 0.1
  xv <- seq(0, 100, dx)    
  # Mortality
  mortList <- lapply(c("f", "m"), function(ss) {
    ids <- grep(ss, colnames(out$theta))
    mort <- apply(out$theta[, ids], 1, function(th) CalcMort(th, xv))
    mortave <- apply(mort, 1, mean)
    mortci <- apply(mort, 1, quantile, c(0.025, 0.975))
    mortfin <- rbind(mortave, mortci)
    return(mortfin)
  })
  # Survival
  survList <- lapply(c("f", "m"), function(ss) {
    ids <- ids <- grep(ss, colnames(out$theta))
    surv <- apply(out$theta[, ids], 1, function(th) CalcSurv(th, xv))
    survave <- apply(surv, 1, mean)
    survci <- apply(surv, 1, quantile, c(0.025, 0.975))
    survfin <- rbind(survave, survci)
    return(survfin)
  })
  # pdf of ages at death
  pdfList <- lapply(c("f", "m"), function(ss) {
    ids <- ids <- grep(ss, colnames(out$theta))
    pdfm <- apply(out$theta[, ids], 1, 
                  function(th) CalcSurv(th, xv) * CalcMort(th, xv))
    pdfave <- apply(pdfm, 1, mean)
    pdfci <- apply(pdfm, 1, quantile, c(0.025, 0.975))
    pdffin <- rbind(pdfave, pdfci)
    return(pdffin)
  })
  # pdf of ages at natal out-migration
  dispList <- list()
  pdf <- apply(out$gam, 1, function(gam) {
    dgamma(xv, gam[1], gam[2])})
  pdfave <- apply(pdf, 1, mean)
  pdfci <- apply(pdf, 1, quantile, c(0.025, 0.975))
  dispList[[1]] <- rbind(pdfave, pdfci)
  # 1 - CDF to calculate cut off points for plotting
  OneMinusCDF <- apply(out$gam, 1, function(gam) {
    pgamma(xv, gam[1], gam[2], lower.tail = F)})
  OneMinusCDFave <- apply(OneMinusCDF, 1, mean)
  
  if(length(out$gam2 > 0)) {
  disp2List <- list()
  pdf <- apply(out$gam2, 1, function(gam2) {
    dgamma(xv, gam2[1], gam2[2])})
  pdfave <- apply(pdf, 1, mean)
  pdfci <- apply(pdf, 1, quantile, c(0.025, 0.975))
  disp2List[[1]] <- rbind(pdfave, pdfci)
  }
  
  EHqMat <- lapply(c("f", "m"), function(ss) {
    ids <- grep(ss, colnames(out$theta))
    xm <- ageMat[ss, 1] - ageFit
    xvm <- xv[xv >= xm]
    surv0 <- apply(out$the[, ids], 1, function(th) CalcSurv(th, xv))
    survM <- apply(out$the[, ids], 1, function(th) CalcSurv(th, xvm) / 
                     CalcSurv(th, xvm[1]))
    E0 <- apply(surv0, 2, CalcEx, dx = dx)
    EM <- apply(survM, 2, CalcEx, dx = dx)
    H0 <- apply(surv0, 2, CalcHx, dx = dx)
    HM <- apply(survM, 2, CalcHx, dx = dx)
    EHq <- rbind(c(mean(E0), quantile(E0, c(0.025, 0.975))),
                 c(mean(EM), quantile(EM, c(0.025, 0.975))),
                 c(mean(H0), quantile(H0, c(0.025, 0.975))),
                 c(mean(HM), quantile(HM, c(0.025, 0.975))))
    dimnames(EHq) <- list(c("E0", "EM", "H0", "HM"), 
                          c("Mean", "2.5%", "97.5%"))
    return(EHq)
  })
  
  plotCut99 <- list()
  for (i in 1:length(survList)) {
    plotCut99[[i]] <- which(survList[[i]][1, ] < 0.01)[1]
    if (is.na(plotCut99[[i]])) plotCut99[[i]] <- length(xv)
  }
  plotCut95 <- list()
  for (i in 1:length(survList)) {
    plotCut95[[i]] <- which(survList[[i]][1, ] < 0.05)[1]
    if (is.na(plotCut95[[i]])) plotCut95[[i]] <- length(xv)
  }
  
  plotCutDisp99 <- which(OneMinusCDFave < 0.01)[1]
  plotCutDisp95 <- which(OneMinusCDFave < 0.05)[1]
  
  names(EHqMat) <- names(survList) <- names(mortList) <- names(pdfList) <- 
    names(plotCut95) <- names(plotCut99) <- c("f", "m")
  if(length(out$gam2) > 0) {
    return(list(mort = mortList, surv = survList, pdf = pdfList, x = xv, 
                EH = EHqMat, cuts99 = plotCut99, cuts95 = plotCut95, 
                cutsDisp99 = plotCutDisp99, cutsDisp95 = plotCutDisp95,
                disp = dispList, disp2 = disp2List))
  } else {
    return(list(mort = mortList, surv = survList, pdf = pdfList, x = xv, 
                EH = EHqMat, cuts99 = plotCut99, cuts95 = plotCut95, 
                cutsDisp99 = plotCutDisp99, cutsDisp95 = plotCutDisp95,
                disp = dispList))
  }
}

# Calculate KLcs:
CalcKLc <- function(p1, p2, low) {
  require(msm)
  pMat <- cbind(p1, p2)
  parRan <- range(sapply(1:2, function(pp) qtnorm(c(0.001, 0.999), 
                                                  mean(pMat[, pp]), 
                                                  sd(pMat[, pp]), 
                                                  lower = low[pp])))
  parVec <- seq(parRan[1], parRan[2], length = 100)
  dp <- parVec[2] - parVec[1]
  parDens <- sapply(1:2, function(pp) 
    dtnorm(parVec, mean(pMat[, pp]), sd(pMat[, pp]), lower = low[pp]))
  # if one of the densities has 0s, restrict densities to the non-0 range
  kld1 <- sum(parDens[, 1] * log(parDens[, 1] / parDens[, 2]) * dp)
  kld2 <- sum(parDens[, 2] * log(parDens[, 2] / parDens[, 1]) * dp)
  qKlc1 <- (1 + (1 - exp(-2 * kld1)^(1 / 2))) / 2
  qKlc2 <- (1 + (1 - exp(-2 * kld2)^(1 / 2))) / 2
  mqKl <- (qKlc1 + qKlc2) / 2
  outList <- list(kl1 = kld1, kl2 = kld2, qkl1 = qKlc1, 
                  qkl2 = qKlc2, mqKl = mqKl)
  return(outList)
}
}