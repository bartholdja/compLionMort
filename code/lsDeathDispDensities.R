rm(list = ls())
if ("fernando" %in% list.files("/Users/")) {
  setwd("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLionsGithub/compLionMort/")
  load("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/data/hwange/hwangeMortAnal.03Sep.rdata")
} else {
  setwd("/Users/Viktualia/Documents/GitHub/compLionMort")
  load("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/data/hwange/hwangeMortAnal.03Sep.rdata")  
}

# Source functions:
source("code/functions.R")

thMal <- t(matrix(c(-1.2, 0.7, 0.16, -3.5, 0.23 )))
xv <- seq(0, 25, 0.1)
class(thMal) <- c("go", "bt")
sigma <- 0.4
mu <- 3
dispDens <- (1/(sigma * sqrt(2 * pi))) * exp(-1 * ((xv-mu)^2/(2 * sigma^2)))
plot(xv, dispDens, type = "l")
plot(xv, dispDens/20, type = "l")

pdf <- CalcPdf(thMal, xv)
pdf <- pdf * xv
plot(xv,pdf, type = "l", ylim = c(0, 0.5))
lines(xv, dispDens/10)
lines(xv, pdf+dispDens/10)