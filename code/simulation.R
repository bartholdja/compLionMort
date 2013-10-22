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
s <- 0.45 # prop. females among 1 year old in Serengeti (Barthold et al in review): 0.45
          # prop. females among fetuses Kruger (Smuts 1978): ~0.43
nMal <- 300
nFem <- s/(1-s) * nMal

cdfYfem <- runif(nFem)
cdfYmal <- runif(nMal)
model <- "go"
shape <- "bt"
ageInt<- seq(0, 30, 0.1)
ageIntM <- seq(0,8, 0.1)
npars <- 5
ncovs <- 2
niter <- 10000
alphaProbUnsex <- -1/2 * log(0.3)  # f(x) = exp(-alpha * x) - beta, f(0) = 0.9, f(2) = 0.2
betaProbUnsex <- 0.1
lamMigr <- -1.83

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

lamMigr <- 1.83 # app. lamMigr from Hwange data
lamNonMigr <- -log(0.00005) / 2

# ages at death
xDfem <- ageInt[intFem]
xDmal <- ageInt[intMal]

# unsex some individuals
dat <- data.frame(1:(nFem + nMal), c(xDfem, xDmal), c(rep("f",nFem), rep("m", nMal)))
names(dat) <- c("id", "ageYrs", "sex")
dat$sexNew <- factor(dat$sex, levels = c("f", "m", "u"))

temp <- NULL
for (i in 1:nrow(dat)) {
  if(dat$ageYrs[i] <= 1) {
    probUnsex <- 0.9 - 0.45 * dat$ageYrs[i] # prob to die unsexed 0.9 age 0, 0.2 age 2, linear decline
    temp <- rbinom(1, 1, probUnsex) } else {
      temp <- 0 }
  if (temp == 1) dat$sexNew[i] <- "u"
}
rm(temp)


# last seen ages (in Hwange 70% of males born in study are lost for data collection)
idM <- which(dat$sex == "m" & rbinom(length(dat$sex == "m"), 1, 0.7) == 1)
idNM <- (1:nrow(dat))[!(1:nrow(dat) %in% idM)]

ageIntM <- seq(2.5, 4.5, 0.1)

# migrators:
cdfYintM <- 1 - (exp(-lamMigr * (ageIntM - 2.5)))
cdfYM <- runif(length(idM))
intM <- findInterval(cdfYM, cdfYintM)
xM <- ageIntM[intM]

dat$ageLastSeen <- dat$ageYrs
for (i in 1 : length(xM)) {
if(dat$ageYrs[idM][i] > xM[i]) dat$ageLastSeen[idM][i] <- xM[i]
}


# non-migrators:
#cdfYintNM <- exp(-lamNonMigr * ageIntM)
#cdfYNMmal <- runif(length(idNMMal))
#intNMmal <- findInterval(cdfYNMmal, sort(cdfYintNM))
#intNMmal <- length(cdfYintNM) - intNMmal
#xLSmal <- xDmal
#xLSmal[idMig] <- xMig

# immigrating males
dat$immigration <- rep(0, nrow(dat))

nImmig <- 300
cdfYimmig <- runif(nImmig)
cdfYintImmig <- cdfYintMal
intImmig <- findInterval(cdfYimmig, cdfYintImmig)
# ages at death for immigrants
xDImmig <- ageInt[intImmig]

potImmig <- xDImmig[xDImmig >= 3]
ageImmig <- sample(potImmig, 33, replace = FALSE)

ageImmigID <- (dat$id[nrow(dat)]+1):(dat$id[nrow(dat)]+length(ageImmig))
sexImmig <- rep("m", length(ageImmig))
sexNewImmig <- rep("m", length(ageImmig))
immigrationImmig <- rep(1, length(ageImmig))

addImmig <- data.frame(ageImmigID, ageImmig, sexImmig, sexImmig, ageImmig,
                       immigrationImmig)
names(addImmig) <- names(dat)
datNew <- rbind(dat, addImmig)

seren <- read.csv("/Users/Viktualia/Dropbox/JuliaLions/data/SerenMalesMort01Aug13.csv",
                  na.strings = "")

serenIm <- subset(seren, seren$bornElsewher == 1)

for (i in 3:8) {
serenIm[ ,i] <- as.Date(serenIm[ ,i])
}

serenIm$datIm <- NULL
serenIm$datIm[1] <- "1900-01-01"
serenIm$datIm <- as.Date(serenIm$datIm)

for (i in 1:nrow(serenIm)) {
  for (j in 4:7){
  if (is.na(serenIm[i, j])) {} else {
    serenIm$datIm[i] <- serenIm[i, j]
  }
  }
}

serenIm$ageIm <- round((as.numeric(serenIm$datIm) -
                          as.numeric(serenIm$birthDate))/365.25, 2)
serenIm <- serenIm[-which(serenIm$ageIm == max(serenIm$ageIm)), ]
median(serenIm$ageIm[serenIm$ageIm >= 2.5])
quantile(serenIm$ageIm[serenIm$ageIm >= 2.5], 0.95)
quantile(serenIm$ageIm[serenIm$ageIm >= 2.5], 0.05)

rm(list = setdiff(ls(), c("dat")))

save.image("/Users//Viktualia/Documents/GitHub/compLionMort/data/simDatHwang.rdata")

