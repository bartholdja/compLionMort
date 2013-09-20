# sourced by LionMortXXX.R

# Functions:
# Basic mortality:
CalcBaseMort <- function(th, ...) UseMethod("CalcBaseMort")

CalcBaseMort.ct <- function(th, x) {
  rep(th, length(x))
}

CalcBaseMort.we <- function(th, x) {
  th[, 1] * th[, 2]^th[, 1] * x^(th[, 1] - 1)
}

CalcBaseMort.go <- function(th, x) {
  exp(th[, 1] + th[, 2] * x)
}

CalcBaseMort.lo <- function(th, x) {
  (exp(th[, 1] + th[, 2] * x)) / 
    (1 + th[, 3] * exp(th[, 1]) / th[, 2] * (exp(th[, 2] * x) - 1))
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
  thb <- th[, -c(1:3)]
  class(thb) <- class(th)
  th[, 1] * exp(-th[, 2] * x) + th[, 3] + CalcBaseMort(thb, x)
}

CalcBaseSurv <- function(th, ...) UseMethod("CalcBaseSurv")

CalcBaseSurv.ct <- function(th, x) {
  exp(-th * x)
}

CalcBaseSurv.we <- function(th, x) {
  exp(-(th[, 2] * x)^th[, 1])
}

CalcBaseSurv.go <- function(th, x) {
  Sx <- exp(exp(th[, 1]) / th[, 2] * (1 - exp(th[, 2] * x)))
  idNoTh12 <- which(th[, 1] == 0)
  Sx[idNoTh12] <- 1
  return(Sx)
}

CalcBaseSurv.lo <- function(th, x) {
  (1 + th[, 3] * exp(th[, 1]) / th[, 2] * (exp(th[, 2] * x) - 1))^(-1 / th[, 3])
}


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
  thb <- th[, -c(1:3)]
  class(thb) <- class(th)
  Sx0 <- Sx0 <- exp(exp(th[, 1]) / th[, 2] * (exp(-th[, 2] * x) - 1) - th[, 3] * x)
  idNoTh12 <- which((th[, 1] != 0 | th[, 1] == 0)) 
  return(Sx0 * CalcBaseSurv(thb, x))
}

# Ages at death density:
CalcPdf <- function(th, x) CalcMort(th, x) * CalcSurv(th, x)

# Life expectancy:
CalcEx <- function(th, minAge = 0, dx = 0.1, xMax = 1000) {
  sum(CalcSurv(th, seq(minAge, xMax, dx)) * dx)
}

# Prior for age distribution:
CalcPriorAgeDist <- function(x, th, exPrior) {   # Prior age distr.: S(x, prior theta) / ex(prior theta)
  CalcSurv(th, x) / exPrior
}

# Define the mortality likelihood:
CalcLikeMort <- function(x, th) {
  CalcPdf(th, x) / CalcSurv(th, ageTrunc) #  denominator becomes 1 for non-truncated individuals
}

# Define the emigration likelihood:
CalcLikeEmigr <- function(x, idM = idMigr, idNM = idNonMigr) {
  like <- x * 0 + 1
  like[idM] <- 1 - exp(-lamMigr * (x[idM] - 1.75))
  like[idNM] <- exp(-lamNonMigr * (x[idNM] - ageToLast[idNM]))  # f(x) = exp(-alpha * x)
  return(like)
}

# Define the full likelihood:
CalcFullLike <- function(x, th, idM = idMigr, idNM = idNonMigr) {
  like <- log(CalcLikeMort(x, th)) + 
    log(CalcLikeEmigr(x, idM = idM, idNM = idNM)) 
  return(like)
}

CalcCovTheta <- function(th, covars = NA) {
  if (is.na(covars[1])) {   # this might not work now if there are no covariates
    theta <- matrix(th, n, defPars$length, byrow = TRUE, 
                    dimnames = list(NULL, defPars$name)) 
  } else {
    theta <- covars %*% th
  }
  class(theta) <- class(th)
  return(theta)
}

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