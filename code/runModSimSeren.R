# Serengeti lions:
rm(list = ls())

library(msm)
if ("fernando" %in% list.files("/Users/")) {
  setwd("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLionsGithub/compLionMort/")
  load("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/data/serengeti/simDatSeren.Rdata")
} else {
  setwd("/Users/Viktualia/Documents/GitHub/compLionMort")
  load("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/data/serengeti/simDatSeren.Rdata")
  #load("N:/Barthold/simDatSeren.Rdata")
  #setwd("N:/Barthold")
}


# Source functions:
source("code/functions.R")
#source("functions.R")

# want diagnostic plots?:
plotInd <- FALSE

### Extract variables:
# number of observations
n <- nrow(dat)
# Age to first and last
ageToFirst <- dat$fsAgeYrs
first <- ageToFirst
ageToLast <- dat$ageLastSeen
# minimum dispersal age
# in Pusey & Packer 1987 all m disp by age 4.2, min 1.8 (1 of 12)
# Elliot et al (subm.) all m disp by 3.75, min 1.66 (no m survived younger than 2.6)
minDispAge <- 1.75
# sex
sex <- as.character(dat$sexNew)
# indicator and index for open fate after last seen ( 0: observed death | alive, 1: open fate)
unknownFate <- dat$noDeath
idNoDeath <- which(unknownFate == 1)
# left truncation (only immigrants for simulated data)
idLeftTr <- which(ageToFirst != 0)
ageTrunc <- rep(0, n)
ageTrunc[idLeftTr] <- ageToFirst[idLeftTr]  # born before study start
# index for unknown sex
idNoSex <- which(sex == "u" )
# probability newborns female
probFem <- 0.45

### Propose (initial) parameter values:
# Mortality
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
# Dispersal pars and state indices
lambdaStart <- c(log(2), 1) 
lambdaJump <- c(0.2, 0.2)
nDispPars <- length(lambdaStart)
idIMstart <- which(sex == "m" & first != 0 & ageToFirst >= minDispAge)
idEMstart <- which(sex == "m" & unknownFate == 1 & ageToLast >= minDispAge)
#idSTstart <- (1:n)[!((1:n) %in% idEMstart)]  # not emigrator
# Ages
xStart <- ageToLast
xStart[idNoDeath] <- xStart[idNoDeath] + sample(1:10, length(idNoDeath), 
                                                replace = TRUE)
# non-dectection probability (if alive and in study area)
detectPar <- -log(0.00005) / 2
# Initial sexes
sexFemStart <- rep(1, n)
sexFemStart[sex == "m"] <- 0
sexFemStart[idNoSex] <- rbinom(length(idNoSex), 1, probFem)
covarsStart  <- matrix(c(sexFemStart, 1 - sexFemStart), n, ncovs, dimnames = list(NULL, names))
# Dispersal state indicator
dispStateStart <- rep(0, n)
dispStateStart[idEMstart] <- 1

### For calculating priors:
xv <- seq(0, 100, 0.1)
thPrior <- matrix(defPars$priorMean, length(xv), defPars$length, byrow = TRUE)
class(thPrior) <- c(model, shape)
priorLamMean <- c(log(3), 2)
priorLamSd <- c(1, 1)
exPrior <- sum(CalcSurv(thPrior, xv) * 0.1)
dispStatePrior <- 0.7

### Build jumps matrix:
jumpMatStart <- matrix(defPars$jump, ncovs, defPars$length, byrow = TRUE,
                       dimnames = dimnames(thetaStart))

### Number of iterations and runs
niterRun <- 5

### Run dynamic Metropolis to find jumps
UpdJumps <- TRUE
niter <- 10000
burnin <- 5001
thin <- 20
keep <- seq(burnin, niter, thin)
outJump <- RunMCMC(1)

par(mfcol = c(5, 2), mar = c(3, 3, 1, 1)) 
for (i in 1:10) plot(outJump$pars[, i], type = 'l', 
                     main = colnames(outJump$pars)[i])

# Run MCMC:
UpdJumps <- FALSE
jumpMatStart <- outJump$jumpsMort
lambdaJump <- outJump$jumpsDisp
rm("outJump")
niter <- 100000
burnin <- 10001
thin <- 50
keep <- seq(burnin, niter, thin)

nsim <- 4
ncpus <- 4
require(snowfall)
sfInit(parallel = TRUE, cpus = ncpus)
sfExport(list = c(ls(), ".Random.seed"))
sfLibrary(msm, warn.conflicts = FALSE)
startTime <- Sys.time()
out <- sfClusterApplyLB(1:nsim, RunMCMC)
sfStop()
(Sys.time() - startTime)

rm(list = setdiff(ls(), c("out", "nsim", "niter", "model", "shape", "ncovs", 
                          "names", "npars", 
                          "defPars", "thetaNames")))

if ("fernando" %in% list.files("/Users/")) {
  save.image("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/results/outputHwang04Nov2.Rdata")
} else {
  save.image("N:/Barthold/simSerenOut15.Rdata")
  #save.image("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/results/simSerenOut4.Rdata")
}