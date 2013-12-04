rm(list = ls())

library(msm)
library(RColorBrewer)
if ("fernando" %in% list.files("/Users/")) {
  setwd("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLionsGithub/compLionMort/")
  load("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/data/hwange/simDatHwang.rdata")
} else {
  setwd("/Users/Viktualia/Documents/GitHub/compLionMort")
  load("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/data/hwange/simDatHwang.rdata")
}

# Source functions:
source("code/functions.R")

# Extract variables:
n <- nrow(dat)
ageTrunc <- rep(0, n)  # no truncation in simulation study
unknownFate <- dat$noDeath
idNoDeath <- which(unknownFate == 1)
ageToLast <- dat$ageLastSeen
sex <- as.character(dat[, 'sexNew'])
minDispAge <- 1.75
maxDispAge <- 4.25

# in Pusey & Packer 1987 all males dispersed by the age of 4.2, minimum age 1.8 (1 ind out of 12)
# Elliot et al (submitted) all males dispersed by the age of 3.75, minimum age 1.66 (no male survived younger than 2.6)
idMigr <- which((sex == "m") & unknownFate == 1 & (
  ageToLast >= minDispAge & ageToLast <= maxDispAge))  
# this is to be deleted later
test <- F
if (test == T) {
  test <- data.frame(ageToLast[ageToLast >= minDispAge & ageToLast <= maxDispAge & unknownFate == 1 & sex == "m"], dat$ageYrs[ageToLast >= minDispAge & ageToLast <= maxDispAge & unknownFate == 1 & sex == "m"])
  test$diff <- test$diff <- test[, 1] - test[ ,2]
  idMigr <- idMigr[-which(test$diff == 0)]
}
idNonMigr <- (1:length(sex))[!(1:length(sex) %in% idMigr)]
idNoSex <- which(sex == "u")
probFem <- 0.45 # proportion females at birth

# Emigration probability of male lions aged 2.5 to 4.25
ageLastMigr <- ageToLast[idMigr]
out <- optimise(LikeMigr, c(0, 10))
lamMigr <- out$minimum

# Non-resighting probability conditioned on being alive and in the study area:
# for everyone other than male lions aged 1.75 to 4.25
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

# Propose initial ages:
xStart <- ageToLast
xStart[idNoDeath] <- xStart[idNoDeath] + sample(1:10, length(idNoDeath), 
                                                replace = TRUE)

# Propose intial sexes:
sexFemStart <- rep(1, n)
sexFemStart[sex == "m"] <- 0
sexFemStart[idNoSex] <- rbinom(length(idNoSex), 1, probFem)
covarsStart  <- matrix(c(sexFemStart, 1 - sexFemStart), n, ncovs, dimnames = list(NULL, names))

# Initial non-/migrators based on initial sexes:
idMstart <- which(sexFemStart == 0 & unknownFate == 1 &
                    ageToLast >= minDispAge & ageToLast <= maxDispAge)  
idNMstart <- (1:n)[!(1:n %in% idMstart)]

# Calculate priors:
xv <- seq(0, 100, 0.1)
thPrior <- matrix(defPars$priorMean, length(xv), defPars$length, byrow = TRUE)
class(thPrior) <- c(model, shape)
exPrior <- sum(CalcSurv(thPrior, xv) * 0.1)

# Build jumps matrix:
jumpMatStart <- matrix(defPars$jump, ncovs, defPars$length, byrow = TRUE,
                       dimnames = dimnames(thetaStart))

# Run dynamic Metropolis to find jumps
UpdJumps <- TRUE
niter <- 5000
outJump <- RunMCMC(1)

# Run MCMC:
UpdJumps <- FALSE
jumpMatStart <- outJump$jumps
niter <- 10000
nsim <- 4
ncpus <- 4
require(snowfall)
sfInit(parallel = TRUE, cpus = ncpus)
sfExport(list = c(ls(), ".Random.seed"))
sfLibrary(msm, warn.conflicts = FALSE)
out <- sfClusterApplyLB(1:nsim, RunMCMC)
sfStop()

rm(list = setdiff(ls(), c("out", "nsim", "niter", "model", "shape", "ncovs", "lamMigr", "lamNonMigr",
                          "names", "defPars", "npars", "dat", "thetaNames", "thetaFemOr", "thetaMalOr")))

if ("fernando" %in% list.files("/Users/")) {
  save.image("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/results/simOut.Rdata")
} else {
  save.image("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/results/simHwangOut1.Rdata")
}

 
