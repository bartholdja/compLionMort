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
nMal <- 1000
nFem <- round(s/(1-s) * nMal)

n <- nMal + nFem
model <- "go"
shape <- "bt"
ageInt<- seq(0, 30, 0.1)
minDispAge <- 1.75
lambdaStart <- c(log(1.25), 1)

plotInd <- TRUE

thFem <- t(matrix(c(-1.4, 0.65, 0.07, -3.8, 0.2)))
thMal <- t(matrix(c(-1.2, 0.7, 0.26, -3.5, 0.23 )))
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
cdfYfem <- runif(nFem)
cdfYmal <- runif(nMal)
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
propLost <- 0.7 # proportion emigrating males lost for data collection
# draw migrators (from sex == "m" before unsexing)
idEM <- which(dat$sex == "m" & dat$ageYrs >= minDispAge)
idEM <- idEM[rbinom(length(idEM), 1, propLost) == 1]                
                
ageIntM <- seq(minDispAge, 25, 0.01)

# migrators:
cdfYintM <- plnorm(ageIntM - minDispAge, lambdaStart[1], lambdaStart[2], 
       log = TRUE)
if(plotInd == TRUE) {
  plot(ageIntM, exp(cdfYintM), type = "l")
}
cdfYM <- runif(length(idEM))
intM <- findInterval(cdfYM, exp(cdfYintM))
xEM <- ageIntM[intM]

dat$ageLastSeen <- dat$ageYrs
dat$ageLastSeen[idEM][dat$ageYrs[idEM] > xEM] <- xEM[dat$ageYrs[idEM] > xEM]
idEMnew <- idEM[dat$ageYrs[idEM] > xEM]

# immigrating males
dat$immigration <- rep(0, n)

# number of immigrating males
nImmig <- round((table(dat$sexNew)[2] * 1.55) - table(dat$sexNew)[2]) #prop. taken from original Serengeti data 2264 males, of that 804 immigrants
# ages at first seen
dat$fsAgeYrs <- rep(0, n)
cdfYfsIM <- runif(2 * nMal)
intM <- findInterval(cdfYfsIM, exp(cdfYintM))
xIM <- ageIntM[intM]

# ages at death for immigrants
cdfYimmig <- runif(2 * nMal)
intImmig <- findInterval(cdfYimmig, cdfYintMal)
xDImmig <- ageInt[intImmig]

datImmig <- data.frame(xIM[which(xDImmig > xIM)], xDImmig[which(xDImmig > xIM)])
names(datImmig) <- c("xIM", "xDIM")
datImmig$index <- 1:nrow(datImmig)

indexImmig <- sample(datImmig$index, nImmig, replace = FALSE)
datImmig <- datImmig[which(datImmig$index %in% indexImmig), 1:2]

datImmig$id <- (nrow(dat) + 1): (nrow(dat) + nrow(datImmig))
datImmig$sex <- rep("m", nrow(datImmig))
datImmig$sexNew <- rep("m", nrow(datImmig))
names(datImmig)[which(names(datImmig) == "xIM")] <- "fsAgeYrs"
names(datImmig)[which(names(datImmig) == "xDIM")] <- "ageYrs"
datImmig$ageLastSeen <- datImmig$ageYrs

datImmig$immigration <- rep(1, nrow(datImmig))
datNew <- rbind(dat, datImmig)
dat <- datNew

# 10 % of all individuals except for emigrators get an observed death
idST <- dat$id[!(dat$id %in% idEM)]

dat$noDeath <- rep(1, nrow(dat))  # 0 means  observed death, 1 no observed death
dat$noDeath[idST][rbinom(length(dat$noDeath[idST]), 1, 0.1) == 1] <- 0

thetaFemOr <- thFem
thetaMalOr <- thMal

dat$dispState <- rep(0, nrow(dat))
dat$dispState[idEM] <- 1 # all potential migrators
dat$dispState[idEMnew] <- 2 # actual migrators

rm(list = setdiff(ls(), c("dat", "thetaFemOr", "thetaMalOr")))

if ("fernando" %in% list.files("/Users/")) {
  save.image("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/results/simOut.Rdata")
} else {
  save.image("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/data/serengeti/simDatSeren.Rdata")
}
