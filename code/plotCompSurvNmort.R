rm(list = ls())

if ("fernando" %in% list.files("/Users/")) {
  setwd("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLionsGithub/compLionMort/")
  load("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/results/outputSeren1.Rdata")
}else {
  setwd("/Users/Viktualia/Documents/GitHub/compLionMort")
  load("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/results/simSerenOut14.Rdata")
}

library(RColorBrewer)
# Source functions:
source("code/functions.R")

pdf("results/traces.pdf", width = 15, height = 10)
par(mfrow = c(ncovs, defPars$length))
for (i in 1:npars) {
  for (sim in 1:nsim) {
    if (sim == 1) plot(out[[sim]]$pars[ ,i], type = 'l', main = thetaNames[i])
    if (sim != 1) lines(out[[sim]]$pars[ ,i], col = brewer.pal(nsim + 2, "Dark2")[sim])
  }
}
dev.off()

pdf("results/parDens005.pdf", width = 12, height = 5)
par(mfrow = c(1, defPars$length))
for(j in 1:nsim){
for (i in 1:defPars$length) {
  plot(density(out[[j]]$pars[-c(1:1000), i]), type = 'l', 
       main = c("a0", "a1", "c", "b0", "b1") [i], 
       col = brewer.pal(nsim + 2, "Dark2")[1],
       lwd = 3, xlim = range(out[[j]]$pars[-c(1:1000), i]),
       ylim = c(0, max (c(max(density(out[[j]]$pars[-c(1:1000), i])[[2]]), 
                          max(density(out[[j]]$pars[-c(1:1000), 
                                              i+defPars$length])[[2]])))))
  lines(density(out[[j]]$pars[-c(1:1000), i + 5]), col = brewer.pal(nsim + 2, "Dark2")[2], lwd = 3)
  if (i == defPars$length) legend("topright", legend = c("female", "male"), 
                                  col = c(brewer.pal(nsim + 2, "Dark2")[1],brewer.pal(nsim + 2, "Dark2")[2]), lwd = c(3, 3))
}
}
dev.off()


# Calculate mean and quantiles for parameter estimates:

thFem <- t(matrix(c(-1.4, 0.65, 0.07, -3.8, 0.2)))
thMal <- t(matrix(c(-1.2, 0.7, 0.26, -3.5, 0.23 )))
class(thFem) <- c(model, shape)
class(thMal) <- class(thFem)


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

# original values
ySurvFem <- CalcSurv(thFem, xv)
ySurvMal <- CalcSurv(thMal, xv)
yMortFem <- CalcMort(thFem, xv)
yMortMal <- CalcMort(thMal, xv)


pdf("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/plots/simSeren7.pdf", height = 10)
par(mfrow = c(2, 1), mar = c(4, 5, 1, 1), family = 'serif')
plot(range(xv[1:max(rangesSurv)]), c(0, 1), col = NA, frame.plot = FALSE, 
     xlab = "", ylab = expression(paste("Survival, ", italic(S(x)))))
for (i in 1:2) {
  polygon(c(xv[1:rangesSurv[i]], rev(xv[1:rangesSurv[i]])), 
          c(survList[[i]][1:rangesSurv[i], 2], 
            rev(survList[[i]][1:rangesSurv[i], 3])), 
          col = c("#66000070", '#164E0170')[i], 
          border =  c("#66000070", "#164E0170")[i])
  lines(xv[1:rangesSurv[i]], survList[[i]][1:rangesSurv[i], 1], col = 'white',
        lwd = 2)
}
legend("topright", names, pch = 15, cex = 2, col = c("#66000070", "#164E0170"),
       bty = 'n')
lines(xv, ySurvFem, col = "#660000", lwd = 2, lty = 3)
lines(xv, ySurvMal, col = "#164E01", lwd = 2, lty = 3)
plot(range(xv[1:max(rangesSurv)]), c(0, max(rangesMort)), col = NA, 
     frame.plot = FALSE, xlab = expression(paste("Age, ", italic(x))), 
     ylab = expression(paste("Mortality, ", mu(italic(x)))), ylim = c(0, 1.5))
for (i in 1:2) {
  polygon(c(xv[1:rangesSurv[i]], rev(xv[1:rangesSurv[i]])), 
          c(mortList[[i]][1:rangesSurv[i], 2], 
            rev(mortList[[i]][1:rangesSurv[i], 3])), 
          col = c("#66000070", "#164E0170")[i], 
          border =  c("#66000070", "#164E0170")[i])
  lines(xv[1:rangesSurv[i]], mortList[[i]][1:rangesSurv[i], 1], col = 'white',
        lwd = 2)
}
lines(xv, yMortFem, col = "#660000", lwd = 2, lty = 3)
lines(xv, yMortMal, col = "#164E01", lwd = 2, lty = 3)

dev.off()
