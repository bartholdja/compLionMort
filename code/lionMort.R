rm(list = ls())

library(msm)
library(RColorBrewer)
if ("fernando" %in% list.files("/Users/")) {
  setwd("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLionsGithub/compLionMort/")
  load("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/data/hwangeMortAnal.03Sep.rdata")
} else {
  setwd("/Users/Viktualia/Documents/GitHub/compLionMort")
  load("/Users/Viktualia/Dropbox/JuliaLions/data/hwangeMortAnal.03Sep.rdata")
}
# Source functions:
source("code/functions.R")


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

# Build jumps matrix:
jumpMatStart <- matrix(defPars$jump, ncovs, defPars$length, byrow = TRUE,
                       dimnames = dimnames(thetaStart))

# Run dynamic Metropolis to find jumps
UpdJumps <- TRUE
niter <- 5000
outjump <- RunMCMC(1)

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






pdf("results/trace009.pdf", width = 15, height = 10)
par(mfrow = c(ncovs, defPars$length))
for (i in 1:npars) {
  for (sim in 1:nsim) {
  if (sim == 1) plot(out[[sim]]$par[ ,i], type = 'l', main = thetaNames[i])
  if (sim != 1) lines(out[[sim]]$par[ ,i], col = brewer.pal(npars-1, "Set1")[sim])
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


