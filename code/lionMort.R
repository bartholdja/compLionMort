rm(list = ls())

library(msm)
library(RColorBrewer)
#setwd("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/")
setwd("/Users/Viktualia/Documents/GitHub/compLionMort")
load("/Users/Viktualia/Dropbox/JuliaLions/data/hwangeMortAnal.03Sep.rdata")

# Extract variables:
study <- julian(as.Date(c("1999-06-01", "2013-06-26")), 
                origin = as.POSIXct("1970-01-01"))
n <- nrow(hwang)
birth <- julian(as.Date(hwang[, "birth.date"]), 
                origin = as.POSIXct("1970-01-01"))
last <- julian(as.Date(hwang[, "death.lastseen.date"]), 
               origin = as.POSIXct("1970-01-01"))
death <- last
death[hwang$alive == 1 | hwang$missing == 1 | hwang$presum.dead == 1] <- NA
idNoDeath <- which(is.na(death))
first <- rep(NA, n)
sex <- as.character(hwang[, 'sex'])
idLeftTr <- which(birth < study[1])
ageTrunc <- apply(cbind(study[1] - birth, 0), 1, max) / 365.25
ageToLast <- (last - birth) / 365.25
# in Pusey & Packer 1987 all males dispersed by the age of 4.2, minimum age 1.8 (1 ind out of 12)
idMigr <- which(sex == 'm' & (hwang$missing == 1 | hwang$presum.dead == 1) &
                  ageToLast >= 1.75 & ageToLast <= 4.25)  
idNonMigr <- which(hwang$alive == 1 | hwang$missing == 1 | hwang$presum.dead == 1)
idNonMigr <- idNonMigr[!(idNonMigr %in% idMigr)] # n = sum(!is.na(death)) + length(idMigr) + length(idNonMigr), everyone accounted for
idNoSex <- which(sex == "u")


# Source functions:
source("code/functions.R")

# Emigration probability of male lions aged 1.75 to 4.25, increasing with time since last seen
# Calculate probability of lions being aged 1.75 to 4.25 years of age
# to have left the area after t time steps not being seen
ageLastMigr <- ageToLast[idMigr]
LikeMigr <- function(par) {
  -sum(log(par) - par * (ageLastMigr - 1.75)) # f(x) = 1 - exp(-alpha * x), F(x) = alpha * exp(-alpha * x)
}
out <- optimise(LikeMigr, c(0, 10))
lamMigr <- out$minimum

# Non-resighting probability conditionon being alive and in the study area:
# for everyone other than male lions aged 1.75 to 4.25
lamNonMigr <- -log(0.00005) / 2

# Propose initial parameter values:
model <- "go"; shape <- "bt"
ncovs <- 2
defPars <- SetDefaultTheta()
npars <- defPars$length * ncovs
names <- c("f", "m")
thetaNames <- paste(rep(defPars$name, ncovs), 
                    rep(names, each = defPars$len), sep = '.')
thetaStart <- matrix(defPars$start, nrow = ncovs, ncol = defPars$length, 
                     byrow = TRUE, dimnames = list(names, 
                                                   defPars$name))
class(thetaStart) <- c(model, shape)

# Output storage objects:
niter <- 10000
niterRun <- 5
for (iterRun in 1:niterRun) {  # one list element per run
  if (iterRun == 1) parList <- NULL  
  parList[[iterRun]] <- matrix(0, niter, npars, dimnames = list(NULL, thetaNames))
  
}
parPostMat <- matrix(0, niterRun, niter, dimnames = list(sprintf("%d%s", 1:5, ".run"), NULL))
for (iterRun in 1:niterRun) {  # one list element per run
  if (iterRun == 1) agePostList <- NULL  
  agePostList[[iterRun]] <- matrix(0, niter, n)
}
for (iterRun in 1:niterRun) {
  if(iterRun == 1) sexList <- NULL
  sexList[[iterRun]] <- matrix(NA, niter, length(idNoSex))
}

# Propose initial ages:
xStart <- c(last - birth) / 365.25
xStart[idNoDeath] <- xStart[idNoDeath] + sample(1:10, length(idNoDeath), 
                                                replace = TRUE)

# Propose intial sexes:
sexFemStart <- rep(1, n)
sexFemStart[sex == "m"] <- 0
sexFemStart[idNoSex] <- rbinom(length(idNoSex), 1, 0.5)
covarsStart  <- matrix(c(sexFemStart, 1 - sexFemStart), n, ncovs, dimnames = list(NULL, names))

# Initial non-/migrators based on initial sexes:
idMstart <- which(sexFemStart == 0 & (hwang$missing == 1 | hwang$presum.dead == 1) &
                  ageToLast >= 1.75 & ageToLast <= 4.25)  
idNMstart <- which(hwang$alive == 1 | hwang$missing == 1 | hwang$presum.dead == 1)
idNMstart <- idNMstart[!(idNMstart %in% idMstart)] # n = sum(!is.na(death)) + length(idMigr) + length

# Calculate priors:
xv <- seq(0, 100, 0.1)
thPrior <- matrix(defPars$priorMean, length(xv), defPars$length, byrow = TRUE)
class(thPrior) <- c(model, shape)
exPrior <- sum(CalcSurv(thPrior, xv) * 0.1)

# MCMC:
for (iterRun in 1:niterRun) {
  # 1. setup "now" values:
  if(iterRun == 1) thetaNow <- thetaStart
  if (iterRun != 1) {thetaNow<- matrix(rtnorm(npars, t(thetaStart), rep(2 * defPars$jump, each = ncovs), 
                                             low = rep(defPars$low, ncovs)), 2, 5, byrow = TRUE, 
                                      dimnames = dimnames(thetaStart)); class(thetaNow) <- class(thetaStart)}
  xNow <- xStart
  sexFemNow <- sexFemStart
  idMnow <- idMstart
  idNMnow <- idNMstart
  covarsNow <- covarsStart
  
  thetaMatNow <- CalcCovTheta(thetaNow, covarsNow)
  likeMortNow <- CalcLikeMort(xStart, thetaMatNow)
  fullLikeNow <- CalcFullLike(xStart, thetaMatNow, idM = idMnow, 
                              idNM = idNMnow)
  parPostNow <- sum(fullLikeNow) + sum(dtnorm(c(thetaNow),  # I don't get this step.
                                            rep(defPars$priorMean, each = ncovs),
                                            rep(defPars$priorSd, each = ncovs), 
                                            low = rep(defPars$low, each = ncovs), 
                                            log = TRUE))
  agePostNow <- fullLikeNow + CalcPriorAgeDist(xNow, thetaMatNow, exPrior)

# Individual runs:
  for (iter in 1:niter) {
    
    # 2. Propose parameters:
    for (pp in 1:npars) {
      thetaNew <- thetaNow
      thetaNew[pp] <- rtnorm(1, thetaNow[pp], rep(defPars$jump, each = ncovs)[pp], 
                             low = rep(defPars$low, each = ncovs)[pp])
      thetaMatNew <- CalcCovTheta(thetaNew, covarsNow)
      likeMortNew <- CalcLikeMort(xNow, thetaMatNew)
      fullLikeNew <- CalcFullLike(xNow, thetaMatNew, idM = idMnow, 
                                  idNM = idNMnow)
      
      parPostNew <- sum(fullLikeNew) + sum(dtnorm(c(thetaNew), 
                                                  rep(defPars$priorMean, each = ncovs),
                                                  rep(defPars$priorSd, each = ncovs), 
                                                  low = rep(defPars$low, each = ncovs), 
                                                  log = TRUE))
      
      r <- exp(parPostNew - parPostNow)
      z <- runif(1)
      if (r > z) {
        thetaNow <- thetaNew
        thetaMatNow <- thetaMatNew
        likeMortNow <- likeMortNew
        fullLikeNow <- fullLikeNew
        parPostNow <- parPostNew
      }
    }
    
    # 3. Update ages at death:
    xNew <- xNow
    xNew[idNoDeath] <- rtnorm(length(idNoDeath), xNow[idNoDeath], 0.1, 
                              low = ageToLast[idNoDeath])
    likeMortNew <- CalcLikeMort(xNew, thetaMatNow)
    fullLikeNew <- CalcFullLike(xNew, thetaMatNow, idM = idMnow, 
                                idNM = idNMnow)
    
    agePostNew <- fullLikeNew + CalcPriorAgeDist(xNew, thetaMatNow, exPrior)
    
    r <- exp(agePostNew - agePostNow)[idNoDeath]
    z <- runif(length(idNoDeath))
    idUpd <- idNoDeath[r > z]
    if (length(idUpd) > 0) {
      likeMortNow[idUpd] <- likeMortNew[idUpd]
      fullLikeNow[idUpd] <- fullLikeNew[idUpd]
      agePostNow[idUpd] <- agePostNew[idUpd]
      xNow[idUpd] <- xNew[idUpd]
    }
    parPostNow <- sum(fullLikeNow) + sum(dtnorm(c(thetaNow), 
                                                rep(defPars$priorMean, each = ncovs),
                                                rep(defPars$priorSd, each = ncovs), 
                                                low = rep(defPars$low, each = ncovs), 
                                                log = TRUE))
    # 4. Update unknown sexes:
    sexFemNew <- sexFemNow
    sexFemNew[idNoSex] <- rbinom(length(idNoSex), 1, 0.5)
    covarsNew <- cbind(sexFemNew, 1 - sexFemNew)
    colnames(covarsNew) <- names
    idMnew <- which(sexFemNew == 0 & (hwang$missing == 1 | hwang$presum.dead == 1) &
                      ageToLast >= 1.75 & ageToLast <= 4.25)  
    idNMnew <- which(hwang$alive == 1 | hwang$missing == 1 | hwang$presum.dead == 1)
    idNMnew <- idNMnew[!(idNMnew %in% idMnew)] # n = sum(!is.na(death)) + length(idMigr) + length
    
    thetaMatNew <- CalcCovTheta(thetaNow, covarsNew)  ## was thetaStart, on purpose?
    likeMortNew <- CalcLikeMort(xNow, thetaMatNow)
    fullLikeNew <- CalcFullLike(xNow, thetaMatNow, idM = idMnew, 
                                idNM = idNMnew)
    
    agePostNew <- fullLikeNew + CalcPriorAgeDist(xNow, thetaMatNew, exPrior)
    
    r <- exp(agePostNew - agePostNow)[idNoSex]
    z <- runif(length(idNoSex))
    idUpd <- idNoSex[r > z]
    if (length(idUpd) > 0) {
      likeMortNow[idUpd] <- likeMortNew[idUpd]
      fullLikeNow[idUpd] <- fullLikeNew[idUpd]
      agePostNow[idUpd] <- agePostNew[idUpd]
      covarsNow[idUpd, ] <- covarsNow[idUpd, ]
      thetaMatNow[idUpd, ] <- thetaMatNew[idUpd, ]
      sexFemNow[idUpd] <- sexFemNew[idUpd]
    }
    idMnow <- which(sexFemNow == 0 & (hwang$missing == 1 | hwang$presum.dead == 1) &
                      ageToLast >= 1.75 & ageToLast <= 4.25)  
    idNMnow <- which(hwang$alive == 1 | hwang$missing == 1 | hwang$presum.dead == 1)
    idNMnow <- idNMnow[!(idNMnow %in% idMnow)] # n = sum(!is.na(death)) + length(idMigr) + length
    
    parPostNow <- sum(fullLikeNow) + sum(dtnorm(c(thetaNow), 
                                                rep(defPars$priorMean, each = ncovs),
                                                rep(defPars$priorSd, each = ncovs), 
                                                low = rep(defPars$low, each = ncovs), 
                                                log = TRUE))
    # 5. store results:
    parList[[iterRun]][iter, ] <- c(t(thetaNow))
    parPostMat[iterRun, iter ] <- parPostNow
    agePostList[[iterRun]][iter, ] <- agePostNow
    sexList[[iterRun]][iter, ] <- sexFemNow[idNoSex]
  }
}

pdf("results/trace008.pdf", width = 15, height = 10)
par(mfrow = c(ncovs, defPars$length))
for (i in 1:npars) {
  for (iterRun in 1:niterRun) {
  if (iterRun == 1) plot(parList[[iterRun]][ ,i], type = 'l', main = thetaNames[i])
  if (iterRun != 1) lines(parList[[iterRun]][ ,i], col = brewer.pal(npars-1, "Set1")[iterRun])
  }
}
dev.off()

# plot densities male and female separate plots:
#pdf("results/parDens005.pdf", width = 10, height = 10)
#par(mfrow = c(2, 5))
#for (i in 1:10) plot(density(parList[-c(1:1000), i]), type = 'l', 
#                    main = thetaNames[i], lwd = 3)
#dev.off()

pdf("results/parDens005.pdf", width = 12, height = 5)
par(mfrow = c(1, defPars$length))
for (i in 1:defPars$length) {
  plot(density(parList[-c(1:1000), i]), type = 'l', 
       main = c("a0", "a1", "c", "b0", "b1") [i], 
       lwd = 3, xlim = range(parList[-c(1:1000), c(i, i+2)]),
       ylim = c(0, max (c(max(density(parList[-c(1:1000), i])[[2]]), 
                          max(density(parList[-c(1:1000), 
                          i+defPars$length])[[2]])))))
  lines(density(parList[-c(1:1000), i + 5]), col = 'dark green', lwd = 3)
if (i == defPars$length) legend("topright", legend = c("female", "male"), 
                                col = c('black', 'dark green'), lwd = c(3, 3))
}
dev.off()


