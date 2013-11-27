if ("fernando" %in% list.files("/Users/")) {
  setwd("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLionsGithub/compLionMort/")
} else {
  setwd("/Users/Viktualia/Documents/GitHub/compLionMort")
}

source("code/functions.R")

# My data simulation:
# Simulate ages at death and sexes:
n <- 1000
study <- c(1970, 2000)
minbirth <- 1965
model <- "go"
shape <- "bt"
nb <- 100
thf <- matrix(c(-1.4, 0.65, 0.07, -3.8, 0.2), 1, 5, 
              dimnames = list("f", c("a0", "a1", "c", "b0", "b1")))
thm <- matrix(c(-1.2, 0.7, 0.16, -3.5, 0.23), 1, 5, 
              dimnames = list("m", c("a0", "a1", "c", "b0", "b1")))
class(thf) <- c(model, shape)
class(thm) <- c(model, shape)
thetaReal <- matrix(c(-1.4, 0.65, 0.07, -3.8, 0.2, -1.2, 0.7, 0.16, -3.5, 0.23),
                    2, 5, byrow = TRUE, 
                    dimnames = list(c("f", "m"), 
                                    c("a0", "a1", "c", "b0", "b1")))

xv <- seq(0, 25, 0.1)
class(thetaReal) <- c(model, shape)
Fxf <- 1 - CalcSurv(thf, xv)
Fxm <- 1 - CalcSurv(thm, xv)
sr <- 0.45
sexf <- rbinom(n, 1, sr)
covar <- cbind(sexf, 1 - sexf)
colnames(covar) <- c("f", "m")
covTheReal <- CalcCovTheta(thetaReal, covar)
nf <- sum(sexf)
nm <- n - nf
bi <- sort(sample(minbirth:study[2], n, replace = TRUE))
xi <- bi*0
xi[sexf == 1] <- xv[findInterval(runif(nf), Fxf)]
xi[sexf == 0] <- xv[findInterval(runif(nm), Fxm)]
di <- bi + xi





