# sourced by LionMortXXX.R

# Functions:
# Basic mortality:
CalcBaseMort <- function(th, ...) UseMethod("CalcBaseMort")
CalcBaseMort.ct <- function(th, x) {
  rep(th, length(x))
}
CalcBaseMort.lo <- function(th, x) {
  (exp(th[, 1] + th[, 2] * x)) / 
    (1 + th[, 3] * exp(th[, 1]) / th[, 2] * (exp(th[, 2] * x) - 1))
}
CalcBaseMort.we <- function(th, x) {
  th[, 1] * th[, 2]^th[, 1] * x^(th[, 1] - 1)
}
CalcBaseMort.go <- function(th, x) {
  exp(th[, 1] + th[, 2] * x)
}

# Mortality by shape (simple, makeham, bathtub):
CalcMort <- function(th, ...) UseMethod("CalcMort")
CalcMort.si <- function(th, x) {   # si stands for simple
  thb <- th
  class(thb) <- class(th)
  CalcBaseMort(thb, x)
}
CalcMort.ma <- function(th, x) {   
  thb <- th[, -1]
  class(thb) <- class(th)
  th[, 1] + CalcBaseMort(thb, x)
}
CalcMort.bt <- function(th, x) {
  thb <- matrix(th[, -c(1:3)], nrow(th)) 
  class(thb) <- class(th)
  exp(th[, 1] - th[, 2] * x) + th[, 3] + CalcBaseMort(thb, x)
}

# Basic survival:
CalcBaseSurv <- function(th, ...) UseMethod("CalcBaseSurv")
CalcBaseSurv.ct <- function(th, x) {
  exp(-th * x)
}
CalcBaseSurv.we <- function(th, x) {
  exp(-(th[, 2] * x)^th[, 1])
}
CalcBaseSurv.lo <- function(th, x) {
  (1 + th[, 3] * exp(th[, 1]) / th[, 2] * (exp(th[, 2] * x) - 1))^(-1 / th[, 3])
}
CalcBaseSurv.go <- function(th, x) {
  Sx <- exp(exp(th[, 1]) / th[, 2] * (1 - exp(th[, 2] * x)))
  idNoTh12 <- which(th[, 1] == 0)
  Sx[idNoTh12] <- 1
  return(Sx)
}

# Survival by shape
CalcSurv <- function(th, ...) UseMethod("CalcSurv")
CalcSurv.si <- function(th, x) {
  thb <- th
  class(thb) <- class(th)
  CalcBaseSurv(thb, x)
}
CalcSurv.ma <- function(th, x) {
  thb <- th[, -1]
  class(thb) <- class(th)
  exp(-th[, 1] * x) * CalcBaseSurv(thb, x)
}
CalcSurv.bt <- function(th, x) {
  thb <- matrix(th[, -c(1:3)], nrow(th)) # avoids function failure if nrow(th) = 1
  class(thb) <- class(th)
  Sx0 <- exp(exp(th[, 1]) / th[, 2] * (exp(-th[, 2] * x) - 1) - th[, 3] * x)
  idNoTh12 <- which((th[, 1] != 0 | th[, 1] == 0)) 
  return(Sx0 * CalcBaseSurv(thb, x))
}

# Ages at death density:
CalcPdf <- function(th, x) CalcMort(th, x) * CalcSurv(th, x)

# Ages at death cdf:
CalcCdf <- function(th, x) {
  cdfx <- 1-CalcSurv(th, x)
  return(cdfx)
}

# Life expectancy:
CalcEx <- function(th, minAge = 0, dx = 0.1, xMax = 1000) {
  sum(CalcSurv(th, seq(minAge, xMax, dx)) * dx)
}

# Prior for age distribution:
CalcPriorAgeDist <- function(x, th, exPrior) {   # Prior age distr.: S(x, prior theta) / ex(prior theta)
  CalcSurv(th, x) / exPrior
}

# Likelihoods
# Define the mortality likelihood:
CalcLikeMort <- function(x, th, logPdf) {
  logPdf - log(CalcSurv(th, ageTrunc))#  denominator becomes 1 for non-truncated individuals
}

# Define the emigration likelihood:
CalcLikeDispPars <- function(lambda, ageToFirst, ageToLast, idIM, likeDispAge) {
  like <- likeDispAge
  like[idIM] <- plnorm(ageToFirst[idIM] - minDispAge, lambda[1], lambda[2], 
                       log = TRUE)
  return(like)
}

CalcLikeDispAge <- function(lambda, ageToFirst, ageToLast, dispState) {
  like <- dispState * dlnorm(ageToLast - minDispAge, lambda[1], 
                             lambda[2], log = TRUE)
  return(like)
}

# Likelihood for sexes:


# Define the ages at death likelihood:
CalcLikeAges <- function(x, th, logPdf, dispState) {
  like <- logPdf - dispState * (detectPar * (x - ageToLast))
  return(like)
}

# Construct theta by covariate matrix:
CalcCovTheta <- function(th, covars = NA) {
  if (is.na(covars[1])) {   # this might not work if there are no covariates
    theta <- matrix(th, n, defPars$length, byrow = TRUE, 
                    dimnames = list(NULL, defPars$name)) 
  } else {
    theta <- covars %*% th
  }
  class(theta) <- class(th)
  return(theta)
}

# Prior for lambda pars:
CalcLambdaPrior <- function(lambda, priorLamMean, priorLamSd) {
  dnorm(lambda[1], priorLamMean[1], priorLamMean[2], log = TRUE) +
    1 / dgamma(lambda[2]^2, priorLamSd[1], priorLamSd[2], log = TRUE)
}

# Set starting values:
SetDefaultTheta  <- function() {
  if (model == "ct") {
    nTh <- 1
    startTh <- 0.2 
    jumpTh <- 0.1
    priorMean <- 0.06
    priorSd <- 1
    nameTh <- "b0"
    lowTh <- 0
    jitter <- 0.5
  } else if (model == "go") {
    nTh <- 2 
    startTh <- c(-2, 0.01) 
    jumpTh <- c(0.1, 0.0375)
    priorMean <- c(-3, 0.01)
    priorSd <- c(1, 1)
    nameTh <- c("b0", "b1")
    lowTh <- c(-Inf, -Inf)
    jitter <- c(0.5, 0.2) 
    if (shape == "bt") {
      lowTh <- c(-Inf, 0)
    }
  } else if (model == "we") {
    nTh <- 2
    startTh <- c(1.5, 0.2) 
    jumpTh <- c(.01, 0.1)
    priorMean <- c(1.5, .05)
    priorSd <- c(1, 1)
    nameTh <- c("b0", "b1")
    lowTh <- c(0, 0)
    jitter <- c(0.5, 0.2) 
  } else if (model == "lo") {
    nTh <- 3 
    startTh <- c(-2, 0.01, 1e-04) 
    jumpTh <- c(0.1, 0.1, 0.1) 
    priorMean <- c(-3, 0.01, 1e-10)
    priorSd <- c(1, 1, 1)
    nameTh <- c("b0", "b1", "b2")
    lowTh <- c(-Inf, 0, 0)
    jitter <- c(0.5, 0.2, 0.5) 
  }
  if (shape == "ma") {
    nTh <- nTh + 1 
    startTh <- c(0, startTh) 
    jumpTh <- c(0.1, jumpTh) 
    priorMean <- c(0, priorMean)
    priorSd <- c(1, priorSd)
    nameTh <- c("c", nameTh)
    lowTh <- c(0, lowTh)
    jitter <- c(0.25, jitter) 
  } else if (shape == "bt") {
    nTh <- nTh + 3 
    startTh <- c(-0.1, 0.6, 0, startTh)
    jumpTh <- c(0.1, 0.1, 0.1, jumpTh) 
    priorMean <- c(-2, 0.01, 0, priorMean)
    priorSd <- c(1, 1, 1, priorSd)
    nameTh <- c("a0", "a1", "c", nameTh)
    lowTh <- c(-Inf, 0, 0, lowTh)
    jitter <- c(0.5, 0.2, 0.2, jitter) 
  }
  defaultTheta  <- list(length = nTh, start = startTh, jump = jumpTh, 
                        priorMean = priorMean, priorSd = priorSd, name = nameTh, 
                        low = lowTh, jitter = jitter)
  return(defaultTheta)
}

# MCMC:
RunMCMC <- function(sim) {
  xNow <- xStart
  sexFemNow <- sexFemStart
  idEMnow <- idEM
  idSTnow <- idST
  covarsNow <- covarsStart
  if (sim == 1) {
    thetaNow <- thetaStart
  } else {
    rm(".Random.seed", envir = .GlobalEnv); runif(1)
    thetaNow <- matrix(rtnorm(npars, t(thetaStart),
                              rep(defPars$jitter, each = ncovs), 
                              low = rep(defPars$low, ncovs)), 
                       2, 5, byrow = TRUE, 
                       dimnames = dimnames(thetaStart)) 
    class(thetaNow) <- class(thetaStart)
  }
  thetaMatNow <- CalcCovTheta(thetaNow, covarsNow)
  logPdfNow <- CalcPdf(thetaMatNow, xNow)
  likeMortNow <- CalcLikeMort(xNow, thetaMatNow, logPdfNow)
  dispStateNow <- rep(0, n)
  dispStateNow[idEMnow] <- 1
  likeAgesNow <- CalcLikeAges(xNow, thetaMatNow, logPdfNow, dispStateNow)
  parPostNow <- sum(likeMortNow) + 
                sum(dtnorm(c(thetaNow), 
                rep(defPars$priorMean, each = ncovs),
                rep(defPars$priorSd, each = ncovs), 
                low = rep(defPars$low, each = ncovs), log = TRUE))
  agePostNow <- likeAgesNow + CalcPriorAgeDist(xNow, thetaMatNow, exPrior)
  
  sexPostNow <- logPdfNow + (sexFemNow) * log(probFem) +
    (1 - sexFemNow) * log(1 - probFem)
  likeDispAgeNow <- CalcLikeDispAge(lambdaNow, ageToFirst, ageToLast, 
                                    dispStateNow)
  likeDispParNow <- CalcLikeDispPars(lambdaNow, ageToFirst, ageToLast, idIM, 
                                     likeDispAgeNow)
  priorLamNow <- CalcLambdaPrior(lambdaNow, priorLamMean, priorLamSd)
  dispPostParNow <- sum(likeDispParNow) + priorLamNow
  dispPostAgeNow <- likeDispAgeNow + dispStateNow * log(dispStatePrior) +
    (1 - dispStateNow) * log(1 - dispStatePrior)
  
  # Output matrices and vectors:
  parMat <- matrix(0, niter, npars, dimnames = list(NULL, thetaNames))
  namesMat <- matrix(thetaNames, ncovs, defPars$len, byrow = TRUE)
  parPostVec <- rep(0, niter)
  agePostMat <- matrix(0, niter, n)
  sexMat <- matrix(NA, niter, length(idNoSex))
  
  # Objects for Dynamic Metropolis:
  jumpMat <- jumpMatStart
  if (UpdJumps) {
    updMat <- parMat
    iterUpd <- 100
    updTarg <- 0.25
    jumpLargeMat <- parMat[0, ]
  }
  
  # Individual runs:
  for (iter in 1:niter) {
    
    # 1. Propose parameters:
    for (pp in 1:npars) {
      thetaNew <- thetaNow
      thetaNew[pp] <- rtnorm(1, thetaNow[pp], jumpMat[pp], 
                             low = rep(defPars$low, each = ncovs)[pp])
      thetaMatNew <- CalcCovTheta(thetaNew, covarsNow)
      logPdfNew <- CalcPdf(thetaMatNew, xNow)
      likeMortNew <- CalcLikeMort(xNow, thetaMatNew, logPdfNew)
      parPostNew <- sum(likeMortNew) + 
                    sum(dtnorm(c(thetaNew), 
                    rep(defPars$priorMean, each = ncovs),
                    rep(defPars$priorSd, each = ncovs), 
                    low = rep(defPars$low, each = ncovs), log = TRUE))
       
      r <- exp(parPostNew - parPostNow)
      if (!is.na(r)) {
        z <- runif(1)
        if (r > z) {
          thetaNow <- thetaNew
          thetaMatNow <- thetaMatNew
          logPdfNow <- logPdfNew
          likeMortNow <- likeMortNew
          likeAgesNow <- CalcLikeAges(xNow, thetaMatNow, logPdfNow, dispStateNow)
          agePostNow <- likeAgesNow + CalcPriorAgeDist(xNow, thetaMatNow, exPrior)
          sexPostNow <- logPdfNow + (sexFemNow) * log(probFem) +
            (1 - sexFemNow) * log(1 - probFem)
          if (UpdJumps) updMat[iter, namesMat[pp]] <- 1
        }
      }
    }
    
    # 2. Dispersal parameters:
    for (pp in 1:2) {
      lambdaNew <- lambdaNow
      lambdaNew[pp] <- rtnorm(1, lambdaNow[pp], lambdaJump[pp], 
                              lower = c(-Inf, 0)[pp])
      likeDispAgeNew <- CalcLikeDispAge(lambdaNew, ageToFirst, ageToLast, 
                                        dispStateNow)
      likeDispParNew <- CalcLikeDispPars(lambdaNew, ageToFirst, ageToLast, idIM, 
                                         likeDispAgeNew)
      priorLamNew <- CalcLambdaPrior(lambdaNew, priorLamMean, priorLamSd)
      dispPostParNew <- sum(likeDispParNew) + priorLamNew
      r <- exp(dispPostParNew - dispPostParNow)
      if (!is.na(r)) {
        if (r > runif(1)) {
          lambdaNow <- lambdaNew
          likeDispParNow <- likeDispParNew
          priorLambNow <- priorLamNew
          dispPostParNow <- dispPostParNew
          likeDispAgeNow <- CalcLikeDispAge(lambdaNow, ageToFirst, ageToLast, 
                                            dispStateNow)
          dispPostAgeNow <- likeDispAgeNow + dispStateNow * 
            log(dispStatePrior) + (1 - dispStateNow) * log(1 - dispStatePrior)
          
        }
      }
    }
    
    # 3. Update ages at death:
    xNew <- xNow
    xNew[idNoDeath] <- rtnorm(length(idNoDeath), xNow[idNoDeath], 0.1, 
                              low = ageToLast[idNoDeath])
    logPdfNew <- CalcPdf(thetaMatNow, xNew)
    likeAgesNew <- CalcLikeAges(xNew, thetaMatNow, logPdfNew, dispStateNow)
    
    agePostNew <- likeAgesNew + CalcPriorAgeDist(xNew, thetaMatNow, exPrior) 
    
    likeMortNew <- CalcLikeMort(xNew, thetaMatNow, logPdfNew)
    sexPostNew <- logPdfNew + (sexFemNow) * log(probFem) +
      (1 - sexFemNow) * log(1 - probFem)
    
    r <- exp(agePostNew - agePostNow)[idNoDeath] # infinities in agePostNew & Now
    z <- runif(length(idNoDeath))
    idUpd <- idNoDeath[r > z]
    idUpd <-idUpd[!(is.na(idUpd))] 
    if (length(idUpd) > 0) {
      logPdfNow[idUpd] <- logPdfNew[idUpd]
      likeMortNow[idUpd] <- likeMortNew[idUpd]
      likeAgesNow[idUpd] <- likeAgesNew[idUpd]
      agePostNow[idUpd] <- agePostNew[idUpd]
      xNow[idUpd] <- xNew[idUpd]
      sexPostNow[idUpd] <- sexPostNew[idUpd]
    }
    parPostNow <- sum(likeMortNow) + 
      sum(dtnorm(c(thetaNow),
                 rep(defPars$priorMean, each = ncovs),
                 rep(defPars$priorSd, each = ncovs), 
                 low = rep(defPars$low, each = ncovs), log = TRUE))
    
    # 4. Update unknown sexes:
    sexFemNew <- sexFemNow
    sexFemNew[idNoSex] <- rbinom(length(idNoSex), 1, 0.5)
    covarsNew <- cbind(sexFemNew, 1 - sexFemNew)
    colnames(covarsNew) <- names
    dispStateNew <- dispStateNow
    dispStateNew[sexFemNew == 1] <- 0
    thetaMatNew <- CalcCovTheta(thetaNow, covarsNew)
    logPdfNew <- CalcPdf(thetaMatNew, xNow)   
    likeAgesNew <- CalcLikeAges(xNew, thetaMatNew, logPdfNew, dispStateNew)    
    agePostNew <- likeAgesNew + CalcPriorAgeDist(xNow, thetaMatNew, exPrior) 
    likeMortNew <- CalcLikeMort(xNow, thetaMatNew, logPdfNew)
    sexPostNew <- logPdfNew + (sexFemNew) * log(probFem) +
      (1 - sexFemNew) * log(1 - probFem)

    r <- exp(sexPostNew - sexPostNow)[idNoSex]
    z <- runif(length(idNoSex))
    
    idUpd2 <- idNoSex[r > z]
    idUpd2 <-idUpd2[!(is.na(idUpd2))]
    if (length(idUpd2) > 0) {
      logPdfNow[idUpd2] <- logPdfNew[idUpd2]
      likeMortNow[idUpd2] <- likeMortNew[idUpd2]
      likeAgesNow[idUpd2] <- likeAgesNew[idUpd2]
      agePostNow[idUpd2] <- agePostNew[idUpd2]
      covarsNow[idUpd2, ] <- covarsNew[idUpd2, ]  # was covarsNow, Jul, 02 Nov
      thetaMatNow[idUpd2, ] <- thetaMatNew[idUpd2, ]
      sexFemNow[idUpd2] <- sexFemNew[idUpd2]
      sexPostNow[idUpd2] <- sexPostNew[idUpd2]
      dispStateNow[idUpd2] <- dispStateNew[idUpd2]
    }
    
    parPostNow <- sum(likeMortNow) + 
      sum(dtnorm(c(thetaNow), rep(defPars$priorMean, each = ncovs),
                 rep(defPars$priorSd, each = ncovs), 
                 low = rep(defPars$low, each = ncovs), log = TRUE))
    idEMnow <- which(sexFemNow == 0 & unknownFate == 1 & ageToLast >= minDispAge) 
    idSTnow <- (1:n)[!((1:n) %in% idEMnow)]  # not emigrator
    likeDispAgeNow <- CalcLikeDispAge(lambdaNow, ageToFirst, ageToLast, 
                                      dispStateNow)
    dispPostAgeNow <- likeDispAgeNow + dispStateNow * 
      log(dispStatePrior) + (1 - dispStateNow) * log(1 - dispStatePrior)
    likeDispParNow <- CalcLikeDispPars(lambdaNow, ageToFirst, ageToLast, idIM, 
                                       likeDispAgeNow)
    dispPostParNow <- sum(likeDispParNow) + priorLamNow
    
    # 5. Dispersal state:
    dispStateNew <- dispStateNow
    dispStateNew[idEMnow] <- rbinom(length(idEMnow), 1, 0.5)    
    likeDispAgeNew <- CalcLikeDispAge(lambdaNow, ageToFirst, ageToLast, 
                                      dispStateNew)
    likeDispParNew <- CalcLikeDispPars(lambdaNew, ageToFirst, ageToLast, idIM, 
                                       likeDispAgeNew)
    dispPostAgeNew <- likeDispAgeNew + dispStateNew * 
      log(dispStatePrior) + (1 - dispStateNew) * log(1 - dispStatePrior)
    
    r <- exp(dispPostAgeNew - dispPostAgeNow)[idEMnow]
    z <- runif(length(idEMnow))
    
    idUpd2 <- idEMnow[r > z]
    idUpd2 <-idUpd2[!(is.na(idUpd2))]     

    if (length(idUpd2) > 0) {
      dispStateNeow[idUpd2] <- dispStateNew[idUpd2]
      likeDispAgeNeow[idUpd2] <- likeDispAgeNew[idUpd2]
      likeDispParNow[idUpd2] <- likeDispParNew[idUpd2]
      dispPostAgeNow[idUpd2] <- dispPostAgeNew[idUpd2]
    }
    likeAgesNow <- CalcLikeAges(xNow, thetaMatNow, logPdfNow, dispStateNow)
    agePostNow <- likeAgesNow + CalcPriorAgeDist(xNow, thetaMatNow, exPrior)
    
    # 2. Dynamic Metropolis to update jumps:
    if (UpdJumps) {
      if (is.element(iter/iterUpd,c(1:50))) {
        updRate <- apply(updMat[iter - ((iterUpd - 1):0), ], 2, sum) / iterUpd
        updRate[updRate == 0] <- 1e-2
        jumpMat <- jumpMat * matrix(updRate, ncovs, defPars$le, byrow = TRUE) / 
          updTarg
        jumpLargeMat <- rbind(jumpLargeMat, c(t(jumpMat)))
      }
    }
    
    

    # 5. store results:
    parMat[iter, ] <- c(t(thetaNow))
    parPostVec[iter] <- parPostNow
    agePostMat[iter, ] <- agePostNow
    sexMat[iter, ] <- sexFemNow[idNoSex]
  }
  if (UpdJumps) {
    aveJumps <- matrix(apply(jumpLargeMat[20:50, ], 2, mean), ncovs, 
                       defPars$len, dimnames = dimnames(jumpMat))
  } else {
    aveJumps <- jumpMat
  }
  return(list(pars = parMat, parPost = parPostVec, agePost = agePostMat,
              sexEst = sexMat, jumps = aveJumps))
}
