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

# Probability that an individual that is not seen anymore has left the study area:
LikeMigr <- function(par) {
  -sum(log(par) - par * (ageLastMigr - minDispAge)) # f(x) = 1 - exp(-alpha * x), F(x) = alpha * exp(-alpha * x)
}

# Likelihoods
# Define the mortality likelihood:
CalcLikeMort <- function(x, th) {
  log(CalcPdf(th, x) / CalcSurv(th, ageTrunc))#  denominator becomes 1 for non-truncated individuals
}

# Define the emigration likelihood:
CalcLikeEmigr <- function(x, idM = idMigr, idNM = idNonMigr) {
  like <- x * 0 + 1
  like[idM] <- 1 - exp(-lamMigr * (x[idM] - minDispAge))  # f(x) = 1 - exp(-alpha * x)
  like[idNM] <- exp(-lamNonMigr * (x[idNM] - ageToLast[idNM]))  # f(x) = exp(-alpha * x)
  return(like)
}

# Define the ages at death likelihood:
CalcLikeAges <- function(x, th, idM = idMigr, idNM = idNonMigr) {
  like <- log(CalcPdf(th, x)) + # some Inf in the log of pdf because the pdf is 0
    log(CalcLikeEmigr(x, idM = idM, idNM = idNM)) 
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

CalcDispLike <- function(xd, xdm, sdxd, disp, resid) {
  dispPdf <- xd * 0
  
}
SampleDispState <- function() {
  
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
  idMnow <- idMstart
  idNMnow <- idNMstart
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
  likeMortNow <- CalcLikeMort(xNow, thetaMatNow)
  likeAgesNow <- CalcLikeAges(xStart, thetaMatNow, idM = idMnow, 
                              idNM = idNMnow)
  parPostNow <- sum(likeMortNow) + 
                sum(dtnorm(c(thetaNow), 
                rep(defPars$priorMean, each = ncovs),
                rep(defPars$priorSd, each = ncovs), 
                low = rep(defPars$low, each = ncovs), log = TRUE))
  agePostNow <- likeAgesNow + CalcPriorAgeDist(xNow, thetaMatNow, exPrior)
  sexPostNow <- likeAgesNow + (sexFemNow) * log(probFem) +
    (1 - sexFemNow) * log(1 - probFem)
  
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
      likeMortNew <- CalcLikeMort(xNow, thetaMatNew)
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
          likeMortNow <- likeMortNew
          likeAgesNow <- CalcLikeAges(xNow, thetaMatNow, idM = idMnow,
                                      idNM = idNMnow)
          agePostNow <- likeAgesNow + CalcPriorAgeDist(xNow, thetaMatNow, exPrior)
          sexPostNow <- likeAgesNow + (sexFemNow) * log(probFem) +
            (1 - sexFemNow) * log(1 - probFem)
          if (UpdJumps) updMat[iter, namesMat[pp]] <- 1
        }
      }
    }
    
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
    
    # 3. Update ages at death:
    xNew <- xNow
    xNew[idNoDeath] <- rtnorm(length(idNoDeath), xNow[idNoDeath], 0.1, 
                              low = ageToLast[idNoDeath])
    likeAgesNew <- CalcLikeAges(xNew, thetaMatNow, idM = idMnow, 
                                idNM = idNMnow)
    agePostNew <- likeAgesNew + CalcPriorAgeDist(xNew, thetaMatNow, exPrior) # inf in likeAgesNew
    
    likeMortNew <- CalcLikeMort(xNew, thetaMatNow)
    sexPostNew <- likeAgesNew + (sexFemNow) * log(probFem) +
      (1 - sexFemNow) * log(1 - probFem)
    
    r <- exp(agePostNew - agePostNow)[idNoDeath] # infinities in agePostNew & Now
    z <- runif(length(idNoDeath))
    idUpd <- idNoDeath[r > z]
    idUpd <-idUpd[!(is.na(idUpd))] 
    if (length(idUpd) > 0) {
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
    idMnew <- which(sexFemNew == 0 & unknownFate == 1 &
                      ageToLast >= minDispAge & ageToLast <= maxDispAge)  
    idNMnew <- (1:n)[!(1:n %in% idMnew)] 
    thetaMatNew <- CalcCovTheta(thetaNow, covarsNew)
    likeMortNew <- CalcLikeMort(xNow, thetaMatNew)
    likeAgesNew <- CalcLikeAges(xNow, thetaMatNew, idM = idMnew, 
                                idNM = idNMnew)
    agePostNew <- likeAgesNew + CalcPriorAgeDist(xNow, thetaMatNew, exPrior)
    
    sexPostNew <- likeAgesNew + (sexFemNew) * log(probFem) +
      (1 - sexFemNew) * log(1 - probFem)
    
    r <- exp(sexPostNew - sexPostNow)[idNoSex]
    z <- runif(length(idNoSex))
    
    idUpd2 <- idNoSex[r > z]
    idUpd2 <-idUpd2[!(is.na(idUpd2))]  # removes the NA's from idUpd (from when pdf is 0, and log like -inf)
    if (length(idUpd2) > 0) {
      likeMortNow[idUpd2] <- likeMortNew[idUpd2]
      likeAgesNow[idUpd2] <- likeAgesNew[idUpd2]
      agePostNow[idUpd2] <- agePostNew[idUpd2]
      covarsNow[idUpd2, ] <- covarsNew[idUpd2, ]  # was covarsNow, Jul, 02 Nov
      thetaMatNow[idUpd2, ] <- thetaMatNew[idUpd2, ]
      sexFemNow[idUpd2] <- sexFemNew[idUpd2]
      sexPostNow[idUpd2] <- sexPostNew[idUpd2]
    }
    
    idMnow <- which(sexFemNow == 0 & unknownFate == 1 &
                      ageToLast >= minDispAge & ageToLast <= maxDispAge)  
    idNMnow <- (1:n)[!(1:n %in% idMnow)]
    
    parPostNow <- sum(likeMortNow) + 
      sum(dtnorm(c(thetaNow), rep(defPars$priorMean, each = ncovs),
                 rep(defPars$priorSd, each = ncovs), 
                 low = rep(defPars$low, each = ncovs), log = TRUE))
    
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
