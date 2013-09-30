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



# Emigration probability of male lions aged 1.75 to 4.25
ageLastMigr <- ageToLast[idMigr]
LikeMigr <- function(par) {
  -sum(log(par) - par * (ageLastMigr - 1.75)) # f(x) = 1 - exp(-alpha * x), F(x) = alpha * exp(-alpha * x)
}
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
  if (sim == 1) plot(out[[sim]]$pars[ ,i], type = 'l', main = thetaNames[i])
  if (sim != 1) lines(out[[sim]]$pars[ ,i], col = brewer.pal(npars-1, "Set1")[sim])
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

# Calculate mean and quantiles for parameter estimates:
keep <- seq(1000, niter, 20)
parMat <- out[[1]]$pars[keep, ]
for (i in 2:nsim) {
  parMat <- rbind(parMat, out[[i]]$pars[keep, ])
}

parQuant <- cbind(apply(parMat, 2, mean), apply(parMat, 2, sd),
                  t(apply(parMat, 2, quantile, c(0.025, 0.975))))
colnames(parQuant) <- c("Mean", "SE", "2.5%", "97.5%")

mortList <- list()
survList <- list()
xv <- seq(0, 25, 0.1)
parnames <- colnames(parMat)
for (nc in 1:ncovs) {
  idpars <- which(substr(parnames, nchar(parnames), nchar(parnames)) == names[nc])
  thetamat <- parMat[, idpars]
  nPar <- ncol(thetamat)
  mort <- apply(thetamat, 1, function(tht) {
    tht2 <- matrix(tht, length(xv), nPar, byrow = TRUE)
    class(tht2) <- c(model, shape)
    CalcMort(tht2, xv)})
  mortList[[names[nc]]] <- cbind(apply(mort, 1, mean), 
                    t(apply(mort, 1, quantile, c(0.025, 0.975))))
  surv <- apply(thetamat, 1, function(tht) {
    tht2 <- matrix(tht, length(xv), nPar, byrow = TRUE)
    class(tht2) <- c(model, shape)
    CalcSurv(tht2, xv)})
  survList[[names[nc]]] <- cbind(apply(surv, 1, mean), 
                                 t(apply(surv, 1, quantile, c(0.025, 0.975))))
  
}

rangesSurv <- sapply(names, function(nn) which(survList[[nn]][, 1] < 0.01)[1])

rangesMort <- sapply(names, function(nn) range(mortList[[nn]][1:rangesSurv[nn], ]))

par(mfrow = c(2, 1), mar = c(4, 5, 1, 1), family = 'serif')
plot(range(xv[1:max(rangesSurv)]), c(0, 1), col = NA, frame.plot = FALSE, 
     xlab = "", ylab = expression(paste("Survival, ", italic(S(x)))))
for (i in 1:2) {
  polygon(c(xv[1:rangesSurv[i]], rev(xv[1:rangesSurv[i]])), 
          c(survList[[i]][1:rangesSurv[i], 2], 
            rev(survList[[i]][1:rangesSurv[i], 3])), 
          col = c('dark red', 'dark green')[i], 
          border =  c('dark red', 'dark green')[i])
  lines(xv[1:rangesSurv[i]], survList[[i]][1:rangesSurv[i], 1], col = 'white',
        lwd = 2)
}
legend('topright', names, pch = 15, cex = 2, col = c('dark red', 'dark green'),
       bty = 'n')

plot(range(xv[1:max(rangesSurv)]), c(0, max(rangesMort)), col = NA, 
     frame.plot = FALSE, xlab = expression(paste("Age, ", italic(x))), 
     ylab = expression(paste("Mortality, ", mu(italic(x)))))
for (i in 1:2) {
  polygon(c(xv[1:rangesSurv[i]], rev(xv[1:rangesSurv[i]])), 
          c(mortList[[i]][1:rangesSurv[i], 2], 
            rev(mortList[[i]][1:rangesSurv[i], 3])), 
          col = c('dark red', 'dark green')[i], 
          border =  c('dark red', 'dark green')[i])
  lines(xv[1:rangesSurv[i]], mortList[[i]][1:rangesSurv[i], 1], col = 'white',
        lwd = 2)
}





