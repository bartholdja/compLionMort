# Simulation study to check Bayesian framework to infer
# age-specific survival in species with sex-biased dispersal

rm(list = ls())

library(msm)
library(RColorBrewer)
if ("fernando" %in% list.files("/Users/")) {
  setwd("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLionsGithub/compLionMort/")
  load("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/data/hwangeMortAnal.03Sep.rdata")
  load("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/results/output.rdata")
} else {
  setwd("/Users/Viktualia/Documents/GitHub/compLionMort")
  load("/Users/Viktualia/Dropbox/JuliaLions/data/hwangeMortAnal.03Sep.rdata")
  load("/Users/Viktualia/Dropbox/JuliaLions/results/output.rdata")
}

# Source functions:
source("code/functions.R")

# starting objects:
s <- 0.45 # prop. females among 1 year old in Serengeti (Barthold et al in review)
nMal <- 250
nFem <- s/0.5 * nMal

cdfYfem <- runif(nFem)
cdfYmal <- runif(nMal)
model <- "go"
shape <- "bt"
ageInt<- seq(0, 30, 0.1)
ageIntM <- seq(0,8, 0.1)
npars <- 5
ncovs <- 2
niter <- 10000

theta <- matrix(round(out[[1]]$pars[niter, ], 2), 1, npars * ncovs)
class(theta) <- c(model, shape)

thetaFem <- matrix(theta[, 1:npars], 1, npars)
class(thetaFem) <- class(theta)
thetaMal <- matrix(theta[,(npars + 1):ncol(theta)], 1, npars)
class(thetaMal) <- class(theta)

cdfYintFem <- CalcCdf(thetaFem, ageInt)
cdfYintMal <- CalcCdf(thetaMal, ageInt)
intFem <- findInterval(cdfYfem, cdfYintFem)
intMal <- findInterval(cdfYmal, cdfYintMal)

lamMigr <- 1.25 # app. lamMigr from Hwange data
lamNonMigr <- -log(0.00005) / 2

# ages at death
xDfem <- ageInt[intFem]
xDmal <- ageInt[intMal]

# last seen ages
idM <- which(xDmal >= 1.75)
idNMmal <- which(xDmal < 1.75)

# migrators:
cdfYintM <- 1-exp(-lamMigr * ageIntM)
cdfYM <- runif(length(idM))
intM <- findInterval(cdfYM, cdfYintM)
xM <- ageIntM[intM] + 1.75

# non-migrators:
cdfYintNM <- exp(-lamNonMigr * ageIntM)
cdfYNMmal <- runif(length(idNMMal))
intNMmal <- findInterval(cdfYNMmal, sort(cdfYintNM))
intNMmal <- length(cdfYintNM) - intNMmal
xLSmal <- xDmal
xLSmal[idMig] <- xMig