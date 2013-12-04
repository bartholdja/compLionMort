# sourced by LionMortXXX.R
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
    lowTh <- c(-Inf, 0, 0, lowTh)  # all parameters bound to 0 (a0 and b0 are in the exp(a0) form in the model, limit exp(-Inf) -> 0)
    jitter <- c(0.5, 0.2, 0.2, jitter) 
  }
  defaultTheta  <- list(length = nTh, start = startTh, jump = jumpTh, 
                        priorMean = priorMean, priorSd = priorSd, name = nameTh, 
                        low = lowTh, jitter = jitter)
  return(defaultTheta)
}

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

# Survival by shape
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
  thb <- matrix(th[, -c(1:3)], nrow(th))     # avoids function failure if nrow(th) = 1
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

### Priors:
# Mortality parameters:
CalcPriorMortPar <- function(th) {
  sum(dtnorm(c(th), rep(defPars$priorMean, each = ncovs),
             rep(defPars$priorSd, each = ncovs), 
             low = rep(defPars$low, each = ncovs), log = TRUE))
}

# Age distribution
CalcPriorAgeDist <- function(x, th) {   # Prior age distr.: S(x, prior theta) / ex(prior theta), standardised so that sums to 1
  log(CalcSurv(th, x) / exPrior)
}
# Sexes
CalcPriorSex <- function(sexFemNow) {
  (sexFemNow) * log(probFem) + (1 - sexFemNow) * log(1 - probFem)
}
# Dispersal pars
CalcLambdaPrior <- function(lambda) {
  dnorm(lambda[1], priorLamMean[1], priorLamMean[2], log = TRUE) +
    1 / dgamma(lambda[2]^2, priorLamSd[1], priorLamSd[2], log = TRUE)
}
# Dipersal state
#CalcPriorDispAge <- function (dispState) {
#  dispState * log(dispStatePrior) + (1 - dispState) * log(1 - dispStatePrior)
#}
# Dipersal state
CalcPriorDispAge <- function (dispState, sexFemNow) {
  prior <- dispState * log(dispStatePrior) + (1 - dispState) * log(1 - dispStatePrior)
  prior[sexFemNow == 1] <- 0
  return(prior)
}

### Likelihoods:
# Mortality
CalcLogPdf <- function(th, x) log(CalcMort(th, x) * CalcSurv(th, x))
CalcLikeMort <- function(x, th, logPdf) {
  logPdf - log(CalcSurv(th, ageTrunc))#  denominator becomes 1 for non-truncated individuals
}

# Dispersal state (dispersal age pdf for emigrants)
#CalcLikeDispAge <- function(lambda, dispState) {
#  like <- dispState * 
#    dlnorm(ageToLast - minDispAge, lambda[1], lambda[2], log = TRUE)
#  return(like)
#}
CalcLikeDispAge <- function(lambda, dispState) {
  like <-  dlnorm(ageToLast - minDispAge, lambda[1], lambda[2], log = FALSE)
  like[dispState  == 1] <- log(like[dispState  == 1])
  like[dispState == 0] <- 0
  return(like)
}
# Dispersal parameters (dispersal age pdf for emigrants + dispersal age cdf for immigrants)
CalcLikeDispPars <- function(lambda, idIM, likeDispAge) {
  like <- likeDispAge
  like[idIM] <- plnorm(ageToFirst[idIM] - minDispAge, lambda[1], lambda[2], 
                       log = TRUE)
  return(like)
}
# Ages at death
CalcLikeAges <- function(x, th, logPdf, dispState) {
  like <- logPdf - (1 - dispState) * (detectPar * (x - ageToLast)) 
  return(like)
}

# Construct theta by covariate matrix:
CalcCovTheta <- function(th, covars = NA) {
  if (is.na(covars[1])) {   # this might not work if there are no covariates
    theta <- matrix(th, n, defPars$length, byrow = TRUE, 
                    dimnames = list(NULL, defPars$name))  # what's the NULL for?
  } else {
    theta <- covars %*% th
  }
  class(theta) <- class(th)
  return(theta)
}

# MCMC:
RunMCMC <- function(sim) {
  startTime <- Sys.time()
  ### Now objects:
  # Mortality pars
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
  # Dispersal pars
  idEMnow <- idEMstart
  idIMnow <- idIMstart
  # idSTnow <- idSTstart
  lambdaNow <- lambdaStart
  # Ages
  xNow <- xStart
  # Sexes
  sexFemNow <- sexFemStart
  # Dispersal state
  dispStateNow <- dispStateStart
  
  ### Likelihoods:
  # Mortality pars
  logPdfNow <- CalcLogPdf(thetaMatNow, xNow)
  likeMortNow <- CalcLikeMort(xNow, thetaMatNow, logPdfNow)
  # Dispersal state
  likeDispAgeNow <- CalcLikeDispAge(lambdaNow, dispStateNow)
  # Dispersal pars
  likeDispParNow <- CalcLikeDispPars(lambdaNow, idIMnow, likeDispAgeNow)
  # Ages
  likeAgeNow <- CalcLikeAges(xNow, thetaMatNow, logPdfNow, dispStateNow)
  
  ### Priors:
  # Dispersal pars
  priorLamNow <- CalcLambdaPrior(lambdaNow)
  # Ages
  priorAgeNow <- CalcPriorAgeDist(xNow, thetaMatNow)
  # Sexes
  priorSexNow <- CalcPriorSex(sexFemNow)
  # Prior disp state
  #priorDispAgeNow <- CalcPriorDispAge(dispStateNow)
  priorDispAgeNow <- CalcPriorDispAge(dispStateNow, sexFemNow)
  
  ### Posteriors:
  # Mortality
  parPostNow <- sum(likeMortNow) + CalcPriorMortPar(thetaNow)
  # Dispersal pars
  dispPostParNow <- sum(likeDispParNow[idEMnow]) + priorLamNow
  # Ages
  agePostNow <- likeAgeNow + priorAgeNow
  # Sexes
  sexPostNow <- logPdfNow + priorSexNow
  # Dispersal state
  dispPostAgeNow <- likeDispAgeNow + priorDispAgeNow
  
  ### Output matrices and vectors:
  # Mortality pars
  parMat <- matrix(0, niter, npars, dimnames = list(NULL, thetaNames))
  namesMat <- matrix(thetaNames, ncovs, defPars$len, byrow = TRUE)
  parPostVec <- rep(0, niter)
  # Dispersal pars
  parDispMat <- matrix(0, niter, nDispPars, dimnames = list(NULL, c("Mean", "SD")))
  parDispPostVec <- rep(0, niter)
  # Ages
  agePostMat <- matrix(0, 0, n)
  # Sexes
  sexMat <- matrix(NA, 0, length(idNoSex))
  # Dispersal state
  dispStateMat <- agePostMat
  ageMat <- agePostMat
  # Objects for Dynamic Metropolis:
  jumpMat <- jumpMatStart
  if (UpdJumps) {
    updMat <- parMat
    iterUpd <- 100
    updTarg <- 0.25
    jumpLargeMat <- parMat[0, ]
    updMatLam <- parDispMat
    jumpLargeMatLam <- parDispMat[0, ]
  }
  
  # Individual runs:
  for (iter in 1:niter) {
    
    # 1. Propose mortality parameters:
    for (pp in 1:npars) {
      thetaNew <- thetaNow
      thetaNew[pp] <- rtnorm(1, thetaNow[pp], jumpMat[pp], # draws new par from truncated normal with mean = previous par value and sd = jump
                             low = rep(defPars$low, each = ncovs)[pp])
      thetaMatNew <- CalcCovTheta(thetaNew, covarsNow)
      
      logPdfNew <- CalcLogPdf(thetaMatNew, xNow)
      likeMortNew <- CalcLikeMort(xNow, thetaMatNew, logPdfNew)
      parPostNew <- sum(likeMortNew) + CalcPriorMortPar(thetaNew)
      
      r <- exp(parPostNew - parPostNow)
      if (!is.na(r)) {
        z <- runif(1)
        if (r > z) {
          thetaNow <- thetaNew
          thetaMatNow <- thetaMatNew
          logPdfNow <- logPdfNew
          likeMortNow <- likeMortNew
          parPostNow <- parPostNew
          likeAgeNow <- CalcLikeAges(xNow, thetaMatNow, logPdfNow, dispStateNow)
          priorAgeNow <- CalcPriorAgeDist(xNow, thetaMatNow)
          agePostNow <- likeAgeNow + priorAgeNow
          sexPostNow <- logPdfNow + priorSexNow
          if (UpdJumps) updMat[iter, namesMat[pp]] <- 1
        }
      }
    }
    
    # 2. Propose dispersal parameters:
    for (pp in 1:2) {
      lambdaNew <- lambdaNow
      lambdaNew[pp] <- rtnorm(1, lambdaNow[pp], lambdaJump[pp], 
                              lower = c(-Inf, 0)[pp])
      likeDispAgeNew <- CalcLikeDispAge(lambdaNew, dispStateNow)
      likeDispParNew <- CalcLikeDispPars(lambdaNew, idIMnow, likeDispAgeNew)
      priorLamNew <- CalcLambdaPrior(lambdaNew)
      dispPostParNew <- sum(likeDispParNew[idEMnow]) + priorLamNew # don't understand the limitation to idEMnow here
      r <- exp(dispPostParNew - dispPostParNow)
      if (!is.na(r)) {
        if (r > runif(1)) {
          lambdaNow <- lambdaNew
          likeDispAgeNow <- likeDispAgeNew
          likeDispParNow <- likeDispParNew
          priorLambNow <- priorLamNew
          dispPostParNow <- dispPostParNew
          dispPostAgeNow <- likeDispAgeNow + priorDispAgeNow
          if (UpdJumps) updMatLam[iter, pp] <- 1
        }
      }
    }
    
    # 3. Update ages at death:
    xNew <- xNow
    xNew[idNoDeath] <- rtnorm(length(idNoDeath), xNow[idNoDeath], 0.1, 
                              low = ageToLast[idNoDeath])
    logPdfNew <- CalcLogPdf(thetaMatNow, xNew)
    likeAgeNew <- CalcLikeAges(xNew, thetaMatNow, logPdfNew, dispStateNow)
    likeMortNew <- CalcLikeMort(xNew, thetaMatNow, logPdfNew)
    
    priorAgeNew <- CalcPriorAgeDist(xNew, thetaMatNow) 
    
    agePostNew <- likeAgeNew + priorAgeNew
    sexPostNew <- logPdfNew + priorSexNow
    
    r <- exp(agePostNew - agePostNow)[idNoDeath]
    z <- runif(length(idNoDeath))
    idUpd <- idNoDeath[r > z]
    idUpd <-idUpd[!(is.na(idUpd))] 
    if (length(idUpd) > 0) {
      xNow[idUpd] <- xNew[idUpd]
      logPdfNow[idUpd] <- logPdfNew[idUpd]
      likeMortNow[idUpd] <- likeMortNew[idUpd]
      likeAgeNow[idUpd] <- likeAgeNew[idUpd]
      priorAgeNow[idUpd] <- priorAgeNew[idUpd]
      agePostNow[idUpd] <- agePostNew[idUpd]
      sexPostNow[idUpd] <- sexPostNew[idUpd]
    }
    parPostNow <- sum(likeMortNow) + CalcPriorMortPar(thetaNow)
    
    # 4. Update unknown sexes:
    sexFemNew <- sexFemNow
    sexFemNew[idNoSex] <- rbinom(length(idNoSex), 1, 0.5)
    covarsNew <- cbind(sexFemNew, 1 - sexFemNew)
    colnames(covarsNew) <- names
    thetaMatNew <- CalcCovTheta(thetaNow, covarsNew)
    idEMnew <- which(sexFemNew == 0 & unknownFate == 1 & 
                       ageToLast >= minDispAge) 
    dispStateNew <- dispStateNow  # in my opinion, these 2 lines leave individuals of unknown sex that were
    dispStateNew[sexFemNew == 1] <- 0 # estimated as females with but now as males with a 0 for dispStateNew     
       
    logPdfNew <- CalcLogPdf(thetaMatNew, xNow)   
    priorSexNew <- CalcPriorSex(sexFemNew)
    sexPostNew <- logPdfNew + priorSexNew
    
    likeMortNew <- CalcLikeMort(xNow, thetaMatNew, logPdfNew)
    likeAgeNew <- CalcLikeAges(xNow, thetaMatNew, logPdfNew, dispStateNew)    
    priorAgeNew <- CalcPriorAgeDist(xNow, thetaMatNew)
    agePostNew <- likeAgeNew + priorAgeNew 

    
    
    r <- exp(sexPostNew - sexPostNow)[idNoSex]  
    z <- runif(length(idNoSex))
    
    idUpd2 <- idNoSex[r > z]
    idUpd2 <-idUpd2[!(is.na(idUpd2))]
    if (length(idUpd2) > 0) {
      sexFemNow[idUpd2] <- sexFemNew[idUpd2]
      covarsNow[idUpd2, ] <- covarsNew[idUpd2, ]
      thetaMatNow[idUpd2, ] <- thetaMatNew[idUpd2, ]
      dispStateNow[idUpd2] <- dispStateNew[idUpd2] 
      logPdfNow[idUpd2] <- logPdfNew[idUpd2]
      priorSexNow[idUpd2] <- priorSexNew[idUpd2]
      sexPostNow[idUpd2] <- sexPostNew[idUpd2]
      likeMortNow[idUpd2] <- likeMortNew[idUpd2]
      likeAgeNow[idUpd2] <- likeAgeNew[idUpd2]
      priorAgeNow[idUpd2] <- priorAgeNew[idUpd2]
      agePostNow[idUpd2] <- agePostNew[idUpd2]

    }
    
    parPostNow <- sum(likeMortNow) + CalcPriorMortPar(thetaNow)
    idEMnow <- which(sexFemNow == 0 & unknownFate == 1 & 
                       ageToLast >= minDispAge)
    idIMnow <-which(sexFemNow == 0 & !is.na(first) & ageToFirst >= minDispAge)
    likeDispAgeNow <- CalcLikeDispAge(lambdaNow, dispStateNow)
    #priorDispAgeNow <- CalcPriorDispAge(dispStateNow)
    priorDispAgeNow <- CalcPriorDispAge(dispStateNow, sexFemNow)
    dispPostAgeNow <- likeDispAgeNow + priorDispAgeNow
    likeDispParNow <- CalcLikeDispPars(lambdaNow, idIMnow, likeDispAgeNow)
    dispPostParNow <- sum(likeDispParNow[idEMnow]) + priorLamNow
    
    # 5. Update dispersal state:
    dispStateNew <- dispStateNow
    dispStateNew[idEMnow] <- rbinom(length(idEMnow), 1, 0.5)   
    
    likeDispAgeNew <- CalcLikeDispAge(lambdaNow, dispStateNew)
    #priorDispAgeNew <- CalcPriorDispAge(dispStateNew)
    priorDispAgeNew <- CalcPriorDispAge(dispStateNew, sexFemNow)
    dispPostAgeNew <- likeDispAgeNew + priorDispAgeNew
    
    likeDispParNew <- CalcLikeDispPars(lambdaNow, idIMnow, likeDispAgeNew)
    
    r <- exp(dispPostAgeNew - dispPostAgeNow)[idEMnow] # changed this to idEMnow, Julia 26 Nov
    z <- runif(length(idEMnow))  # same here
    
    idUpd2 <- idEMnow[r > z]
    idUpd2 <-idUpd2[!(is.na(idUpd2))]     
    
    if (length(idUpd2) > 0) {
      dispStateNow[idUpd2] <- dispStateNew[idUpd2]
      likeDispAgeNow[idUpd2] <- likeDispAgeNew[idUpd2]
      likeDispParNow[idUpd2] <- likeDispParNew[idUpd2]
      priorDispAgeNow[idUpd2] <- priorDispAgeNew[idUpd2]
      dispPostAgeNow[idUpd2] <- dispPostAgeNew[idUpd2]
    }
    likeAgeNow <- CalcLikeAges(xNow, thetaMatNow, logPdfNow, dispStateNow)
    agePostNow <- likeAgeNow + CalcPriorAgeDist(xNow, thetaMatNow)
    dispPostParNow <- sum(likeDispParNow[idEMnow]) + priorLamNow
    
    # 6. Dynamic Metropolis to update jumps:
    if (UpdJumps) {
      if (is.element(iter/iterUpd,c(1:100))) {
        # update jumps for Mort pars.:
        updRate <- apply(updMat[iter - ((iterUpd - 1):0), ], 2, sum) / iterUpd  # value between 0 and 1
        updRate[updRate == 0] <- 1e-2
        jumpMat <- jumpMat * 
          matrix(updRate, ncovs, defPars$le, byrow = TRUE) / updTarg
        jumpLargeMat <- rbind(jumpLargeMat, c(t(jumpMat)))
        
        # update jump for disp. pars.:
        updRate <- apply(updMatLam[iter - ((iterUpd - 1):0), ], 2, sum) /
          iterUpd
        updRate[updRate == 0] <- 1e-2
        lambdaJump <- lambdaJump * updRate / updTarg
        jumpLargeMatLam <- rbind(jumpLargeMatLam, lambdaJump)
      }
    }
    
    # 5. store results:
    parMat[iter, ] <- c(t(thetaNow))
    parPostVec[iter] <- parPostNow
    parDispMat[iter, ] <- lambdaNow
    if (iter %in% keep) {
      sexMat <- rbind(sexMat, sexFemNow[idNoSex])
      dispStateMat <- rbind(dispStateMat, dispStateNow)
      agePostMat <- rbind(agePostMat, agePostNow)
      ageMat <- rbind(ageMat, xNow)
    
  }
}
    if (UpdJumps) {
    aveJumps <- matrix(apply(jumpLargeMat[50:100, ], 2, mean), ncovs, 
                       defPars$len, dimnames = dimnames(jumpMat))
    aveJumpsLam <- apply(jumpLargeMatLam[50:100, ], 2, mean)
  } else {
    aveJumps <- jumpMat
    aveJumpsLam <- lambdaJump
  }
(Sys.time() - startTime)
  return(list(pars = parMat, parPost = parPostVec, agePost = agePostMat,  # why don't we return the lambda related objects?
              sexEst = sexMat, dispState = dispStateMat, 
              jumpsMort = aveJumps, jumpsDisp = aveJumpsLam, parsDisp = parDispMat))

}
