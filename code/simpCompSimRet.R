rm(list = ls())

if ("fernando" %in% list.files("/Users/")) {
  setwd("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLionsGithub/compLionMort/")
  load("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/results/simOut3.Rdata")
} else {
  setwd("/Users/Viktualia/Documents/GitHub/compLionMort")
  load("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/results/simSerenOut10.Rdata")
  load("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/data/serengeti/simDatSeren.Rdata")
}

source("code/functionsNew.R")

plotInd <- TRUE
nsim <- 4
xv <- seq(0, 25, 0.1)
ySurvFem <- CalcSurv(thetaFemOr, xv)
ySurvMal <- CalcSurv(thetaMalOr, xv)
yMortFem <- CalcMort(thetaFemOr, xv)
yMortMal <- CalcMort(thetaMalOr, xv)


keep <- seq(1000, niter, 20)
parMat <- out[[1]]$pars[keep, ]
for (i in 2:nsim) {
  parMat <- rbind(parMat, out[[i]]$pars[keep, ])
}

parQuant <- cbind(apply(parMat, 2, mean), apply(parMat, 2, sd),
                  t(apply(parMat, 2, quantile, c(0.025, 0.975))))
colnames(parQuant) <- c("Mean", "SE", "2.5%", "97.5%")

thetaFemNew <-  matrix(parQuant[1:5, 1], 1, 5)
colnames(thetaFemNew) <- thetaNames[1:5]
colnames(thetaFemOr) <- thetaNames[1:5]
thetaMalNew <-  matrix(parQuant[6:10, 1], 1, 5)
colnames(thetaMalNew) <- thetaNames[6:10]
colnames(thetaMalOr) <- thetaNames[1:5]

class(thetaMalNew) <- class(thetaMalOr)
class(thetaFemNew) <- class(thetaFemOr)

ySurvFemNew <- CalcSurv(thetaFemNew, xv)
ySurvMalNew <- CalcSurv(thetaMalNew, xv)
yMortFemNew <- CalcMort(thetaFemNew, xv)
yMortMalNew <- CalcMort(thetaMalNew, xv)

if ("fernando" %in% list.files("/Users/")) {
  pdf("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/plots/hwange04Nov.pdf")
} else {
  pdf("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/plots/wtDeads.pdf")
}

if(plotInd){
  plot(xv, ySurvFem, col = 2, type = "l", lwd = 2, ylab = "Survival probability", xlab = "Age (yrs)")
  lines(xv, ySurvMal, col = 4, lwd = 2)

  lines(xv, ySurvFemNew, col = 2, lwd = 2, lty = 2)
  lines(xv, ySurvMalNew, col = 4, lwd = 2, lty = 2)
  legend("topright", c("females sim", "males sim", "females retrieved", "males retrieved"),
         lwd = c(2,2, 2, 2), lty = c(1, 1, 2, 2), col = c(2,4, 2, 4))
} 
dev.off()  


if(plotInd){
  plot(xv, yMortFem, col = 2, type = "l", lwd = 2, ylab = "Mortality rate", xlab = "Age (yrs)")
  lines(xv, yMortMal, col = 4, lwd = 2)
  
  lines(xv, yMortFemNew, col = 2, lwd = 2, lty = 2)
  lines(xv, yMortMalNew, col = 4, lwd = 2, lty = 2)
  legend("topright", c("females sim", "males sim", "females retrieved", "males retrieved"),
         lwd = c(2,2, 2, 2), lty = c(1, 1, 2, 2), col = c(2,4, 2, 4))
} 

