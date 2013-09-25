rm(list = ls())

library(msm)
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

# Propose initial values:
model <- "go"; shape <- "bt"
ncovs <- 2
xStart <- c(last - birth) / 365.25
xStart[idNoDeath] <- xStart[idNoDeath] + sample(1:10, length(idNoDeath), 
                                                replace = TRUE)   # start ages

# Sex 
sexFem <- rep(1, n)
sexFem[sex == "m"] <- 0
sexFemNow <- sexFem
sexFemNow[idNoSex] <- rbinom(length(idNoSex), 1, 0.5)
covars <- cbind(sexFemNow, 1 - sexFemNow)
colnames(covars) <- c("f", "m")
idMnow <- which(sexFemNow == 0 & (hwang$missing == 1 | hwang$presum.dead == 1) &
                  ageToLast >= 1.75 & ageToLast <= 4.25)  
idNMnow <- which(hwang$alive == 1 | hwang$missing == 1 | hwang$presum.dead == 1)
idNMnow <- idNMnow[!(idNMnow %in% idMnow)] # n = sum(!is.na(death)) + length(idMigr) + length

defPars <- SetDefaultTheta()
thetaStart <- matrix(defPars$start, nrow = ncovs, ncol = defPars$length, 
                     byrow = TRUE, dimnames = list(colnames(covars), 
                                                   defPars$name))
class(thetaStart) <- c(model, shape)
thetaMatStart <- CalcCovTheta(thetaStart, covars)
thetaNames <- paste(rep(defPars$name,  ncovs), 
                    rep(colnames(covars), each = defPars$len), sep = '.')
npars <- length(thetaNames)
likeMortStart <- CalcLikeMort(xStart, thetaMatStart)
fullLikeStart <- CalcFullLike(xStart, thetaMatStart, idM = idMnow, 
                              idNM = idNMnow)

# Calculate priors:
xv <- seq(0, 100, 0.1)
thPrior <- matrix(defPars$priorMean, length(xv), defPars$length, byrow = TRUE)
class(thPrior) <- c(model, shape)
exPrior <- sum(CalcSurv(thPrior, xv) * 0.1)

# setup "now" values:
thetaNow <- thetaStart
thetaMatNow <- thetaMatStart
likeMortNow <- likeMortStart
fullLikeNow <- fullLikeStart
parPostNow <- sum(fullLikeNow) + sum(dtnorm(c(thetaNow),  # I don't get this step.
                                            rep(defPars$priorMean, each = ncovs),
                                            rep(defPars$priorSd, each = ncovs), 
                                            low = rep(defPars$low, each = ncovs), 
                                            log = TRUE))
xNow <- xStart
agePostNow <- fullLikeNow + CalcPriorAgeDist(xNow, thetaMatNow, exPrior)
dNow <- birth + xNow

# Output matrices:
niter <- 10000
parMat <- matrix(0, niter, npars, dimnames = list(NULL, thetaNames))
parPostVec <- rep(0, niter)
agePostMat <- matrix(0, niter, n)
sexMat <- matrix(NA, niter, length(idNoSex))

for (iter in 1:niter) {
  
  # 1. Propose parameters:
  for (pp in 1:npars) {
    thetaNew <- thetaNow
    thetaNew[pp] <- rtnorm(1, thetaNow[pp], rep(defPars$jump, each = ncovs)[pp], 
                           low = rep(defPars$low, each = ncovs)[pp])
    thetaMatNew <- CalcCovTheta(thetaNew, covars)
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
  
  # 2. Update ages at death:
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
    xNow[idUpd] <- xNew[idUpd] # added by Julia
  }
  parPostNow <- sum(fullLikeNow) + sum(dtnorm(c(thetaNow), 
                                              rep(defPars$priorMean, each = ncovs),
                                              rep(defPars$priorSd, each = ncovs), 
                                              low = rep(defPars$low, each = ncovs), 
                                              log = TRUE))
  # 3. Update unknown sexes:
  sexFemNew <- sexFemNow
  sexFemNew[idNoSex] <- rbinom(length(idNoSex), 1, 0.5)
  covarsNew <- cbind(sexFemNew, 1 - sexFemNew)
  colnames(covarsNew) <- c("f", "m")
  idMnew <- which(sexFemNew == 0 & (hwang$missing == 1 | hwang$presum.dead == 1) &
                    ageToLast >= 1.75 & ageToLast <= 4.25)  
  idNMnew <- which(hwang$alive == 1 | hwang$missing == 1 | hwang$presum.dead == 1)
  idNMnew <- idNMnew[!(idNMnew %in% idMnew)] # n = sum(!is.na(death)) + length(idMigr) + length
  
  thetaMatNew <- CalcCovTheta(thetaStart, covarsNew)  
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
    covars[idUpd, ] <- covars[idUpd, ]
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
  # 3. store results:
  parMat[iter, ] <- c(t(thetaNow))
  parPostVec[iter] <- parPostNow
  agePostMat[iter, ] <- agePostNow
  sexMat[iter, ] <- sexFemNow[idNoSex]
}

pdf("results/trace005.pdf", width = 13, height = 10)
par(mfrow = c(ncovs, defPars$length))
for (i in 1:npars) plot(parMat[, i], type = 'l', main = thetaNames[i])
dev.off()

# plot densities male and female separate plots:
#pdf("results/parDens005.pdf", width = 10, height = 10)
#par(mfrow = c(2, 5))
#for (i in 1:10) plot(density(parMat[-c(1:1000), i]), type = 'l', 
#                    main = thetaNames[i], lwd = 3)
#dev.off()

pdf("results/parDens005.pdf", width = 10, height = 5)
par(mfrow = c(1, defPars$length))
for (i in 1:defPars$length) {
  plot(density(parMat[-c(1:1000), i]), type = 'l', 
       main = c("a0", "a1", "c", "b0", "b1") [i], 
       lwd = 3, xlim = range(parMat[-c(1:1000), c(i, i+2)]),
       ylim = c(0, max (c(max(density(parMat[-c(1:1000), i])[[2]]), 
                          max(density(parMat[-c(1:1000), 
                          i+defPars$length])[[2]])))))
  lines(density(parMat[-c(1:1000), i + 5]), col = 'dark green', lwd = 3)
if (i == defPars$length) legend("topright", legend = c("female", "male"), 
                                col = c('black', 'dark green'), lwd = c(3, 3))
}
dev.off()


