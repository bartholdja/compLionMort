# Simulation study to check Bayesian framework to infer
# age-specific survival in species with sex-biased dispersal

rm(list = ls())

library(msm)
library(RColorBrewer)
if ("fernando" %in% list.files("/Users/")) {
  setwd("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLionsGithub/compLionMort/")
  load("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/data/hwange/hwangeMortAnal.03Sep.rdata")
  } else {
  setwd("/Users/Viktualia/Documents/GitHub/compLionMort")
  load("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/data/hwange/hwangeMortAnal.03Sep.rdata")  
}

# Source functions:
source("code/functions.R")

# starting objects:
s <- 0.45 # prop. females among 1 year old in Serengeti (Barthold et al in review): 0.45
# prop. females among fetuses Kruger (Smuts 1978): ~0.43
nMal <- 300
nFem <- round(s/(1-s) * nMal)

n <- nMal + nFem
cdfYfem <- runif(nFem)
cdfYmal <- runif(nMal)
model <- "go"
shape <- "bt"
ageInt<- seq(0, 30, 0.1)
minDispAge <- 1.75
maxDispAge <- 4.25
lamNonMigr <- -log(0.00005) / 2
plotInd <- TRUE


thFem <- t(matrix(c(-1.4, 0.65, 0.07, -3.8, 0.2)))
thMal <- t(matrix(c(-1.2, 0.7, 0.16, -3.5, 0.23 )))
xv <- seq(0, 25, 0.1)
class(thFem) <- c(model, shape)
class(thMal) <- class(thFem)

yMortFem <- CalcMort(thFem, xv)
yMortMal <- CalcMort(thMal, xv)

if (plotInd) {
  plot(xv, yMortFem, col = 2, type = "l", lwd = 2, ylab = "Mortality hazard", 
       xlab = "Age (yrs)")
  lines(xv, yMortMal, col = 4, lwd = 2)
  legend("topleft", c("females", "males"), lwd = c(2,2), col = c(2,4))
}

ySurvFem <- CalcSurv(thFem, xv)
ySurvMal <- CalcSurv(thMal, xv)

if(plotInd) {
  plot(xv, ySurvFem, col = 2, type = "l", lwd = 2, 
       ylab = "Survival probability", xlab = "Age (yrs)")
  lines(xv, ySurvMal, col = 4, lwd = 2)
  legend("topright", c("females", "males"), lwd = c(2,2), col = c(2,4))
}

# ages at death
cdfYintFem <- CalcCdf(thFem, ageInt)
cdfYintMal <- CalcCdf(thMal, ageInt)
intFem <- findInterval(cdfYfem, cdfYintFem)
intMal <- findInterval(cdfYmal, cdfYintMal)

xDfem <- ageInt[intFem]
xDmal <- ageInt[intMal]

# create data frame
dat <- data.frame(1:(nFem + nMal), c(xDfem, xDmal), c(rep("f",nFem), 
                                                      rep("m", nMal)))
names(dat) <- c("id", "ageYrs", "sex")


# unsex some individuals
dat$sexNew <- factor(dat$sex, levels = c("f", "m", "u"))


indDiedUnsex <- rep(0,n)
probUnsex <- NULL
for (i in 1:n) {
  if(dat$ageYrs[i] <= 2) {
    probUnsex[i]<- 0.9 - 0.35 * dat$ageYrs[i] # prob to die unsexed 0.9 age 0, 0.2 age 2, linear decline
    indDiedUnsex[i] <- rbinom(1, 1, probUnsex[i])
  }
  if (indDiedUnsex[i] == 1) dat$sexNew[i] <- "u"
}
rm(indDiedUnsex, probUnsex)


# last seen ages 
# this indicates the proportion of males that are despite dispersal not 
# lost for data collection
# (in Hwange ~70% of males born in study are lost for data collection)
propLost <- 1  # set to 1 for because of how model is set-up at the moment
idM <- which(dat$sex == "m" & rbinom(length(dat$sex == "m"), 1, propLost) == 1)
idNM <- (1:nrow(dat))[!(1:nrow(dat) %in% idM)]

ageIntM <- seq(minDispAge, maxDispAge, 0.1)

# migrators:
cdfYintM <- 1 - (exp(-lamMigr * (ageIntM - minDispAge)))
if(plotInd == TRUE) {
  plot(ageIntM - minDispAge, cdfYintM, type = "l")
}
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
dat$immigration <- rep(0, n)

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
immigrationInd <- rep(1, length(ageImmig))
addImmig <- data.frame(ageImmigID, ageImmig, sexImmig, sexImmig, ageImmig,
                      immigrationInd)
names(addImmig) <- names(dat)
datNew <- rbind(dat, addImmig)
dat <- datNew

# all individ last seen < 1.75, 10% of females, 10 % of immigrating males 
# get an observed death
# of the migrators 10 % of those with death ages between min and max dispersal age
# get an observed death
dat$noDeath <- rep(NA, nrow(dat))  # 0 means  observed death, 1 no observed death
dat$noDeath[dat$ageYrs <= 1.75] <- 0 # everyone died younger than 1.75 yrs assumed dead

dat$noDeath[dat$sex == "f"] <- 0 # 10 % of females have an observed death
dat$noDeath[dat$sex == "f" & dat$ageYrs > 1.75][rbinom(length(dat$noDeath[dat$sex == "f" & dat$ageYrs > 1.75]), 1, 0.9) == 1] <- 1 # 10 % of females by chance observed death

dat$noDeath[dat$immigration == 1] <- 0 # 10 % of the immigrating males have an observed death
dat$noDeath[dat$immigration == 1][rbinom(length(dat$noDeath[dat$immigration == 1]), 1, 0.9) == 1] <- 1  # 10 % of females by chance observed death

dat$noDeath[dat$ageYrs > minDispAge & dat$ageYrs <= maxDispAge & dat$sex == "m" & dat$immigration == 0] <- 0 # 10 % of males with death ages between min and max dispersal age have observed ddeath
dat$noDeath[dat$ageYrs > minDispAge & dat$ageYrs <= maxDispAge & dat$sex == "m" & dat$immigration == 0][rbinom(length(dat$noDeath[dat$ageYrs > minDispAge & dat$ageYrs <= maxDispAge & dat$sex == "m" & dat$immigration == 0]), 1, 0.9) == 1] <- 1

dat$noDeath[is.na(dat$noDeath)] <- 1 # everyone else no observed death date

thetaFemOr <- thFem
thetaMalOr <- thMal

rm(list = setdiff(ls(), c("dat", "thetaFemOr", "thetaMalOr")))

if ("fernando" %in% list.files("/Users/")) {
  save.image("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/results/simOut.Rdata")
} else {
  save.image("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/results/simOut.Rdata")
}
save.image("/Users//Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/data/hwange/simDatHwang.rdata")
