rm(list = ls())

library(RColorBrewer)

if ("fernando" %in% list.files("/Users/")) {
  setwd("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLionsGithub/compLionMort/")
  load("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/results/simOut3.Rdata")
  load("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/data/hwange/simDatHwang.rdata")
} else {
  setwd("/Users/Viktualia/Documents/GitHub/compLionMort")
  load("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/data/hwange/hwangeMortAnal.03Sep.rdata")
  load("/Users/Viktualia/Desktop/outputHwang.Rdata")
}

# Source functions:
source("code/functions.R")


# Calculate mean and quantiles for parameter estimates:
#nsim <- 4
#model <- "go"  
#niter <- 10000
#shape <- "bt"
#ncovs <- 2
#names <- c("f", "m")

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

pdf("/Users/Viktualia/Desktop/hwangI.pdf", height = 10)
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
     ylab = expression(paste("Mortality, ", mu(italic(x)))), ylim = c(0, 1.5))
for (i in 1:2) {
  polygon(c(xv[1:rangesSurv[i]], rev(xv[1:rangesSurv[i]])), 
          c(mortList[[i]][1:rangesSurv[i], 2], 
            rev(mortList[[i]][1:rangesSurv[i], 3])), 
          col = c('dark red', 'dark green')[i], 
          border =  c('dark red', 'dark green')[i])
  lines(xv[1:rangesSurv[i]], mortList[[i]][1:rangesSurv[i], 1], col = 'white',
        lwd = 2)
}

dev.off()
