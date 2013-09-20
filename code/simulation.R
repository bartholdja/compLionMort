# Simulation study to check Bayesian framework to infer
# age-specific survival in species with sex-biased dispersal
library(LambertW)

# Set values
nFem <- 250
nMal <- 250
idFem <- 1:nFem
idMal <- 1:nMal
parsFem <- c(exp(-1.97), 0.1)  # b0.f, b1.f
parsMal <- c(exp(-2.12), 0.16)  # b0.f, b1.f

# Functions:
# gompertz quantile function
g.quant <- function(pars,x) {
  a <- pars[1]
  b <- pars[2]
  inside <- 1-log(x)*b/a
  out <- 1/b*log(inside)
  return(out)
}

# gompertz-makeham quantile function (Jodra 2009)
gm.quant <- function(pars, x) {
  a <- pars[1]
  b <- pars[2]
  m <- pars[3]  # makeham term
  in.lambertW <- a/m*exp(a/m-b/m*log(x))
  out <- a/(b*m)-1/m*log(x)-1/b*W(in.lambertW)
  return(out)
}


#This a function to simulate lifetimes from a Gompertz distribution
sim.g <- function(pars, n) {
  ns <- runif(n) 
  out <- g.quant(pars,x=ns)
  return(out)
}

#This is a function to simulate lifetimes from a Gompertz-Makeham distribution
sim.gm <- function(pars, n) {
  a <- pars[1]
  b <- pars[2]
  m <- pars[3]  # makeham term
  ns <- runif(n) 
  out <- gm.quant(a=a,b=b,m=m,x=ns)
  return(out)
}

# simulate 
ageDeathFem <- sim.g(n = nFem, pars = parsFem)
ageDeathMal <- sim.g(n = nMal, pars = parsMal)

# 