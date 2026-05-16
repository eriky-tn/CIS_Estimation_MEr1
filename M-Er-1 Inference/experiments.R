################################################################################
# Title: Performance Measure Estimations in Erlang Single-Server Queues
#
# Reference:
# Gomes, E. S., Cruz, F. R. B., Quinino, R. C., & Singh, S. K. (2025).
# Performance measure estimations in Erlang single-server queues.
# Communications in Statistics - Simulation and Computation, 1–14.
# https://doi.org/10.1080/03610918.2025.2501750
#
# Authors:
# Eriky S. Gomes
# Frederico R. B. Cruz
#
# Affiliation:
# Department of Statistics
# Universidade Federal de Minas Gerais (UFMG)
#
# Contact:
# eriky-tn@ufmg.br
# fcruz@est.ufmg.br
#
# Copyright (c) 2025 Gomes & Cruz
# Version: v2026
#
# Description:
# R script developed for computational experiments of the article
################################################################################

rm(list = ls())

# Format to Excel
FrExcel <- function(data){
  write.table(format(data, decimal.mark = ','),
              'clipboard', sep='\t', na = '*')
}
################################################################################
# Likelihood Function and Estimator
################################################################################

Likelihood <- function(p, x, r) {
  n <- length(x)
  y <- sum(x)
  const = 1
  for (i in 1:n) {
    const = const * factorial(x[i] + r - 1) / factorial(r - 1) / factorial(x[i])
  }
  return(const * (p / r) ^ y * (1 + p / r) ^ (-y - n * r))
}

MLERho <- function(x, r, p0 = 1e-3) {
  tol = 1e-3
  n <- length(x)
  y <- sum(x)
  if (y / n < 1)
    return(y / n)
  else
    return(1 - tol)
}

MLELq <- function(x, r, p0 = 0.95) {
  rhoMLE <- min(MLERho(x, r), p0)
  return((r + 1) / (2 * r) * rhoMLE ^ 2 / (1 - rhoMLE))
}

MLELs <- function(x, r, p0 = 0.95) {
  rhoMLE <- min(MLERho(x, r), p0)
  lqMLE <- MLELq(x, r)
  return(lqMLE + rhoMLE)
}

################################################################################
# Auxiliary functions
################################################################################

RegIB <- function(a, b, csi) {
  dRegIB <- function(u, a, b) {
    return(u ^ (a - 1) * (1 - u) ^ (b - 1))
  }
  return(integrate(dRegIB, 0, csi, a, b)[[1]] / beta(a, b))
}

GauHyp <- function(a, b, c, z) {
  stopifnot((c > b) & (b > 0))
  tol = 1e-10
  dGauHyp <- function(u, a, b, c, z) {
    return(u ^ (b - 1) * (1 - u) ^ (c - b - 1) * (1 - z * u) ^ (-a))
  }
  return(integrate(dGauHyp, 0 + tol, 1 - tol, a, b, c, z)[[1]] / beta(b, c - b))
}

PostSelfRho <- function(post, ...) {
  tol = 1e-10
  aux <- function(p, ...) {
    return(p * post(p, ...))
  }
  return(integrate(aux, 0 + tol, 1 - tol, ...)[[1]])
}

PostSelfLq <- function(post, r, ...) {
  tol = 1e-10
  aux <- function(p, ...) {
    return((r + 1) / (2 * r) * p ^ 2 / (1 - p) * post(p, ...))
  }
  return(integrate(aux, 0 + tol, 1 - tol, ...)[[1]])
}

PostSelfLs <- function(post, r, ...) {
  return(PostSelfRho(post, ...) + PostSelfLq(post, r, ...))[[1]]
}

PostRho <- function(p, x, r, lik, prior, ...) {
  n <- length(x)
  y <- sum(x)
  numerator <- function(p, x, r) {
    return(lik(p, x, r) * prior(p, r))
  }
  denominator <- function(x, r) {
    return(integrate(numerator, 0, 1, x, r)[[1]])
  }
  return(numerator(p, x, r) / denominator(x, r))
}

LqToRhoConverter <- function(Lq, r) {
  c = (r + 1) / (2 * r)
  p <- (sqrt(Lq ^ 2 + 4 * c * Lq) - Lq) / (2 * c)

  return(p)
}

BisectionRF <- function(f, lInf, lSup, tol = 1e-10, ...){
  k <- 0
  lInf <- lInf + tol
  lSup <- lSup - tol
  while(abs(lSup - lInf) > tol){
    f1 <- f(lInf, ...)
    f2 <- f(lSup, ...)
    x <- (lInf + lSup) / 2
    fx <- f(x, ...)
    ifelse(f1 * fx <= 0, lSup <- x, lInf <- x)
    k <- k + 1
  }
  root <- (lInf + lSup) / 2
  #cat(k,'iterations\n')
  return(root)
}

LsToRhoConverter <- function(Ls, r) {
  fAux <- function(p) {
    return(Ls - p * (1 + (r + 1) / (2 * r) * p / (1 - p)))
  }
  p <- BisectionRF(fAux, 0, 1)
  return(p)
}

################################################################################
# Jeffreys estimators
################################################################################

JefPrior <- function(p, r) {
  return(p ^ (-1 / 2) * (1 + p / r) ^ (-1 / 2))
}

JefPost <- function(p, x, r) {
  n <- length(x)
  y <- sum(x)
  return(
    r ^ (-1) * (p / r) ^ (y - 1 / 2) * (1 + p / r) ^ (-y - n * r - 1 / 2) /
      beta(y + 1 / 2, n * r) / RegIB(y + 1 / 2, n * r, 1 / (1 + r))
  )
}

JefSelfRhoEst <- function(x, r) {
  n <- length(x)
  y <- sum(x)
  rst <-
    r * (y + 1 / 2) / (n * r - 1) * RegIB(y + 3 / 2, n * r - 1, 1 / (1 + r)) /
    RegIB(y + 1 / 2, n * r, 1 / (1 + r))
  return(rst)
}


BayLqSelfEst <- function(x, r, post, ...) {
  p0 = 0.95
  constNorm <- integrate(post, 0, p0, x, r, ...)[[1]] ^ (-1)
  fAux <- function(p) {
    return((r + 1) / (2 * r) * p ^ 2 / (1 - p) * post(p, x, r, ...))
  }
  return(constNorm * integrate(fAux, 0, p0)[[1]])
}

BayLsSelfEst <- function(x, r, pEstimator, post, ...) {
  pEst <- pEstimator(x,r,...)
  lqEst <- BayLqSelfEst(x, r, post, ...)
  return(pEst + lqEst)
}


################################################################################
# Gaussian hypergeometric estimators
################################################################################

GauPrior <- function(p, r, a, b, c) {
  const <- 1 / beta(b, c - b) / GauHyp(a, b , c, -1 / r)
  return(const * p ^ (b - 1) * (1 - p) ^ (c - b - 1) * (1 + p / r) ^ (-a))
}

GauPost <- function(p, x, r, a, b, c) {
  n <- length(x)
  y <- sum(x)
  return(
    p ^ (y + b - 1) * (1 + p / r) ^ (-y - n * r - a) * (1 - p) ^ (c - b - 1) /
      beta(y + b, c - b) / GauHyp(y + n * r + a, y + b, c + y, -1 / r)
  )
}

GauSelfRhoEst <- function(x, r, a, b, c) {
  n <- length(x)
  y <- sum(x)
  rst <-
    (y + b) / (y + c) * GauHyp(y + n * r + a, y + b + 1, y + c + 1, -1 / r) /
    GauHyp(y + n * r + a, y + b, c + y, -1 / r)
  return(rst)
}


################################################################################
# Monte Carlo simulation
################################################################################

MonteCarlo <- function(p, size, r, fEst, ...) {
  set.seed(2024)
  rep <- 1000
  x <- numeric(size)
  estimates <- numeric(rep)
  for (i in 1:rep) {
    x <- rnbinom(size, r, 1 / (1 + p / r))
    estimates[i] <- fEst(x, r, ...)
  }
  return(estimates)
  #c(mean(estimates), sd(estimates))
}

# the number of stages 'r' cannot be a vector
MonteCarloTab <- function(p, size, r, fEst, ...) {
  tab <- matrix(
    nrow = length(p) * length(size),
    ncol = 4,
    dimnames = list(NULL, c('rho', 'size', 'mean', 'sd'))
  )
  for (i in 1:length(p)) {
    for (j in 1:length(size)) {
      est <- c(p[i], size[j], MonteCarlo(p[i], size[j], r, fEst, ...))
      tab[(i - 1) * length(size) + j,] <- est
    }
  }
  return(tab)
}

MonteCarloSimulation <- function(p, size, r, fEst, ...) {
  simulation <- vector(mode = 'list', length = length(r))
  for (i in 1:length(r)) {
    simulation[[i]] <- MonteCarloTab(p, size, r[i], fEst, ...)
    cat(i / length(r) * 100, '%\n')
  }
  return(simulation)
}

################################################################################
# A numerical example
################################################################################

set.seed(2025)
p0 = 0.95
size = 50
r = 5
p = 0.75
Lq <- (r + 1) / (2 * r) * p^2 / (1 - p)
Ls <- Lq + p
# an uniform gaussian hypergeometric prior
a = 0
b = 1
c = 2 

x <- rnbinom(size, r, 1 / (1 + p / r))
p_mle <- MLERho(x, r)
p_gh <- GauSelfRhoEst(x, r, a, b, c)
p_jef <- JefSelfRhoEst(x, r)

lq_mle <- MLELq(x, r)
lq_gh <- BayLqSelfEst(x, r, GauPost, a, b, c)
lq_jef <- BayLqSelfEst(x, r, JefPost)

ls_mle <- MLELs(x, r)
ls_gh <- BayLsSelfEst(x, r, GauSelfRhoEst, GauPost, a, b, c)
ls_jef <- BayLsSelfEst(x, r, JefSelfRhoEst, JefPost)


cat(p, Lq, Ls)
cat(p_mle, p_gh, p_jef)
cat(lq_mle, lq_gh, lq_jef)
cat(ls_mle, ls_gh, ls_jef)
################################################################################
# Parameters of the simulation study
################################################################################

p = c(0.01, 0.1, 0.2, 0.5, 0.7, 0.9, 0.99)
Lq = c(0.5, 1, 2, 3, 4, 5)
Ls = c(0.5, 1, 2, 3, 4, 5)
size = c(10, 20, 50, 100, 200)
pLq1 = sapply(Lq, LqToRhoConverter, r = 1)
pLq2 = sapply(Lq, LqToRhoConverter, r = 2)
pLs1 = sapply(Ls, LsToRhoConverter, r = 1)
pLs2 = sapply(Ls, LsToRhoConverter, r = 2)

################################################################################
# maximum likelihood estimates
################################################################################

mlRhoEst1 <- MonteCarloSimulation(p, size, r = 1, MLERho)
mlRhoEst2 <- MonteCarloSimulation(p, size, r = 2, MLERho)

mlLqEst1 <- MonteCarloSimulation(pLq1, size, r = 1, MLELq, p0 = 0.95)
mlLqEst2 <- MonteCarloSimulation(pLq2, size, r = 2, MLELq, p0 = 0.95)

mlLsEst1 <- MonteCarloSimulation(pLs1, size, r = 1, MLELs, p0 = 0.95)
mlLsEst2 <- MonteCarloSimulation(pLs2, size, r = 2, MLELs, p0 = 0.95)

mlRhoEst1
mlRhoEst2
mlLqEst1
mlLqEst2
mlLsEst1
mlLsEst2

save(mlRhoEst1,
     mlRhoEst2,
     mlLqEst1,
     mlLqEst2,
     mlLsEst1,
     mlLsEst2,
     file = 'data/MLE.rdata')

#load(file='data/MLE.rdata')
#rm(mlRhoEst, mlLqEst1, mlLqEst2, mlLsEst1, mlLsEst2)

################################################################################
# jeffreys estimates
################################################################################

jefRhoEst1 <- MonteCarloSimulation(p, size, r = 1, JefSelfRhoEst)
jefRhoEst2 <- MonteCarloSimulation(p, size, r = 1, JefSelfRhoEst)

jefLqEst1 <- MonteCarloSimulation(
  pLq1, size, r = 1, BayLqSelfEst, JefPost)

jefLqEst2 <- MonteCarloSimulation(
  pLq2, size, r = 2, BayLqSelfEst, JefPost)

jefLsEst1 <- MonteCarloSimulation(
  pLs1, size, r = 1, BayLqSelfEst, JefPost)

jefLsEst2 <- MonteCarloSimulation(
  pLs2, size, r = 2, BayLsSelfEst, JefSelfRhoEst, JefPost)

jefRhoEst1
jefRhoEst2
jefLqEst1
jefLqEst2
jefLsEst1
jefLsEst2

save(jefRhoEst1,
     jefRhoEst2,
     jefLqEst1,
     jefLqEst2,
     jefLsEst1,
     jefLsEst2,
     file = 'data/Jef.rdata')

#load(file='data/Jef.rdata')
#rm(jefRhoEst,jefLsEst)

################################################################################
# gaussian estimates
################################################################################

a <- 2
b <- 3
c <- 5

gauRhoEst1 <- MonteCarloSimulation(p, size, r = 1, GauSelfRhoEst, a, b, c)
gauRhoEst2 <- MonteCarloSimulation(p, size, r = 2, GauSelfRhoEst, a, b, c)

gauLqEst1 <- MonteCarloSimulation(
  pLq1, size, r = 1, BayLqSelfEst, GauPost, a, b, c)

gauLqEst2 <- MonteCarloSimulation(
  pLq2, size, r = 2, BayLqSelfEst, GauPost, a, b, c)

gauLsEst1 <- MonteCarloSimulation(
  pLs1, size, r = 1, BayLsSelfEst, GauSelfRhoEst, GauPost, a, b, c)

gauLsEst2 <- MonteCarloSimulation(
  pLs2, size, r = 2, BayLsSelfEst, GauSelfRhoEst, GauPost, a, b, c)

gauRhoEst1
gauRhoEst2
gauLqEst1
gauLqEst2
gauLsEst1
gauLsEst2

save(gauRhoEst1,
     gauRhoEst2,
     gauLqEst1,
     gauLqEst2,
     gauLsEst1,
     gauLsEst2,
     file = 'data/GH.rdata')

#load(file = 'data/GH.rdata')
#rm(gauRhoEst,gauLqEst,gauLsEst)

#system('shutdown -s')

################################################################################
# Tests
################################################################################

# testing MLE
p <- c(0.9)#(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
size <- c(100)#(10,20,50,100,200)
r <- c(4)#(1,4)
set.seed(2022)
MCrep <- 1000
estimates <- numeric(MCrep)
for (i in 1:MCrep) {
  x <- rnbinom(size, r, 1 / (1 + p / r))
  (estimates[i] <- MLERho(x, r))
}
estimates <- estimates[!is.na(estimates)]
c(mean(estimates), var(estimates))
MonteCarlo(p, size, r, MLERho)

#
(Lq <- p ^ 2 / (1 - p))
size <- c(100)#(10,20,50,100,200)
r <- c(4)#(1,4)
pLq <- LqToRhoConverter(Lq)
set.seed(2022)
MCrep <- 1000
estimates <- numeric(MCrep)
for (i in 1:MCrep) {
  x <- rnbinom(size, r, 1 / (1 + p / r))
  (estimates[i] <- MLELq(x, r))
}
estimates <- estimates[!is.na(estimates)]
c(mean(estimates), sd(estimates))
MonteCarlo(pLq, size, r, MLELq)

#
(Ls <- p / (1 - p))
size <- c(100)#(10,20,50,100,200)
r <- c(4)#(1,4)
pLs <- LsToRhoConverter(Ls)
set.seed(2022)
MCrep <- 1000
estimates <- numeric(MCrep)
for (i in 1:MCrep) {
  x <- rnbinom(size, r, 1 / (1 + p / r))
  (estimates[i] <- MLELs(x, r))
}
estimates <- estimates[!is.na(estimates)]
c(mean(estimates), var(estimates))
MonteCarlo(pLs, size, r, MLELs)

# testing Jeffreys
p <- c(0.5)#(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
size <- c(100)#(10,20,50,100,200)
r <- c(1)#(1,4)
set.seed(2022)
MCrep <- 1000
estimates <- numeric(MCrep)
for (i in 1:MCrep) {
  x <- rnbinom(size, r, 1 / (1 + p / r))
  (estimates[i] <- JefSelfRhoEst(x, r))
}
c(mean(estimates), var(estimates))
par(mfrow = c(1, 1))
MonteCarlo(p, size, r, JefSelfRhoEst)

#
(Lq <- p ^ 2 / (1 - p))
size <- c(100)#(10,20,50,100,200)
r <- c(1)#(1,4)
pLq <- LqToRhoConverter(Lq)
set.seed(2022)
MCrep <- 1000
estimates <- numeric(MCrep)
for (i in 1:MCrep) {
  x <- rnbinom(size, r, 1 / (1 + pLq / r))
  (estimates[i] <- JefSelfLqEst(x, r))
}
c(mean(estimates), var(estimates))
par(mfrow = c(1, 1))
MonteCarlo(pLq, size, r, JefSelfLqEst)

#
(Ls <- p / (1 - p))
size <- c(100)#(10,20,50,100,200)
r <- c(4)#(1,4)
pLs <- LsToRhoConverter(Ls)
set.seed(2022)
MCrep <- 1000
estimates <- numeric(MCrep)
for (i in 1:MCrep) {
  x <- rnbinom(size, r, 1 / (1 + pLs / r))
  estimates[i] <- JefSelfLsEst(x, r)
}
c(mean(estimates), var(estimates))
par(mfrow = c(1, 1))
MonteCarlo(pLs, size, r, JefSelfLsEst)

# testing Gauss hypergeometric
p <- c(0.5)#(0.01,0.10,0.20,0.50,0.70,0.90,0.99)
size <- c(100)#(10,20,50,100,200)
r <- c(1)#(1,4)
a <- 3
b <- 2
c <- 7
set.seed(2022)
MCrep <- 1000
estimates <- numeric(MCrep)
for (i in 1:MCrep) {
  x <- rnbinom(size, r, 1 / (1 + p / r))
  estimates[i] <- GauSelfRhoEst(x, r, a, b, c)
}
c(mean(estimates), var(estimates))
par(mfrow = c(1, 1))
MonteCarlo(p, size, r, GauSelfRhoEst, a, b, c)

#
(Lq <- p ^ 2 / (1 - p))
size <- c(100)#(10,20,50,100,200)
r <- c(1)#(1,4)
pLq <- LqToRhoConverter(Lq)
set.seed(2022)
MCrep <- 1000
estimates <- numeric(MCrep)
for (i in 1:MCrep) {
  x <- rnbinom(size, r, 1 / (1 + pLq / r))
  estimates[i] <- GauSelfLqEst(x, r, a, b, c)
}
c(mean(estimates), var(estimates))
par(mfrow = c(1, 1))
MonteCarlo(pLq, size, r, GauSelfLqEst, a, b, c)

#
(Ls <- p / (1 - p))
size <- c(100)#(10,20,50,100,200)
r <- c(4)#(1,4)
pLs <- LsToRhoConverter(Ls)
set.seed(2022)
MCrep <- 1000
estimates <- numeric(MCrep)
for (i in 1:MCrep) {
  x <- rnbinom(size, r, 1 / (1 + pLs / r))
  estimates[i] <- GauSelfLsEst(x, r, a, b, c)
}
c(mean(estimates), var(estimates))
par(mfrow = c(1, 1))
MonteCarlo(pLs, size, r, GauSelfLsEst, a, b, c)

# testing jeffreys and Gauss hypergeometric posterior
set.seed(2022)
n <- 50
r <- c(1, 2, 3, 4, 5)
a <- 1
b <- 2
c <- 3
p <- c(0.01, 0.1, 0.2, 0.5, 0.7, 0.9, 0.99)
tabJefErrors <-
  matrix(nrow = length(p),
         ncol = length(r),
         dimnames = list(p, r))
tabGauErrors <-
  matrix(nrow = length(p),
         ncol = length(r),
         dimnames = list(p, r))
for (i in 1:length(p)) {
  for (j in 1:length(r)) {
    x <- rnbinom(n, r[j], 1 / (1 + p[i] / r[j]))
    tabJefErrors[i, j] <-
      JefPost(p[i], x, r[j]) - JefPost.Def(p[i], x, r[j])
    tabGauErrors[i, j] <- GauPost(p[i], x, r[j], a, b, a + b) -
      GauPostDef(p[i], x, r[j], a, b, a + b)
  }
}
tabJefErrors
max(tabJefErrors)
tabGauErrors
max(tabGauErrors)
rm(n, r, p, tabJefErrors, tabGauErrors, x, a, b, c, i, j)

# plotting posteriors
{
  n <- 50
  r <- 1
  a <- 1
  b <- 2
  c <- 3
  p_real <- 0.8
  p <- seq(0, 1, 1e-5)
  x <- rnbinom(n, r, 1 / (1 + p_real / r))
  plot(p, JefPost(p, x, r), 'l', col = 'black') # jeffreys
  lines(p, GauPost(p, x, r, a, b, c), 'l', col = 'red') # gaussian
  # right tail behaviour near upper limit
  max(JefPost(p, x, r)) / max(GauPost(p, x, r, a, b, c))
}
rm(n, r, a, b, c, p_real, p, x)


# testing Lq convergence
fAuxGau <- function(p, x, r, a, b, c) {
  return(p ^ 2 / (1 - p) * GauPost(p, x, r, a, b, c))
}
fAuxJef <- function(p, x, r) {
  return(p ^ 2 / (1 - p) * JefPost(p, x, r))
}
f <- function(p) {
  return(1 / (1 - p) ^ 0.999)
}
integrate(f, 0, 1)[[1]]

n <- 100
a <- 1
b <- 2
c <- 3
p = 0.9
r = 2
Lq = (r + 1) / (2 * r) * p / (1 - p)
Lq

x <- rnbinom(n, r, 1 / (1 + p / r))
integrate(fAuxJef, 0, 0.9999999, x, r)[[1]]
integrate(fAuxGau, 0, 1, x, r, a, b, c)[[1]]
rm(a, b, c, p, r, x, fAuxJef, fAuxGau)


# plotting jeffreys prior for different r stages
{
  #setEPS()
  #postscript(file='FiJF.eps',width=10.5*0.75,height=8*0.75)
  par(mfrow = c(1, 1))
  p <- seq(0.01, 0.99, 0.05)
  plot (p, JefPrior(p, r = 1), 'n', ylab = 'Jeffreys density')
  lines(p,
        JefPrior(p, r = 1),
        'b',
        col = 'black',
        pch = 1)
  lines(p, JefPrior(p, r = 2), 'b', col = 'red', pch = 2)
  lines(p,
        JefPrior(p, r = 4),
        'b',
        col = 'darkgreen',
        pch = 3)
  lines(p,
        JefPrior(p, r = 10),
        'b',
        col = 'darkblue',
        pch = 4)
  legend(
    'top',
    col = c('black', 'red', 'darkgreen', 'darkblue'),
    lwd = 2,
    legend = c(
      expression("r=1"),
      expression("r=2"),
      expression("r=5"),
      expression("r=10")
    ),
    pch = c(1, 2, 3, 4),
    cex = 0.9
  )
}
graphics.off()
rm(p)

# plotting gaussian hypergeometric prior for different r stages
{
  setEPS()
  postscript(file = 'FiGH.eps',
             width = 0.75 * 10.5,
             height = 0.75 * 8)
  par(mfrow = c(1, 1))
  p <- seq(0.01, 0.99, 0.05)
  plot (p, GauPrior(
    p,
    r = 1,
    a = 3,
    b = 2,
    c = 7
  ), 'n', ylab = 'GH density')
  lines(p,
        GauPrior(
          p,
          r = 1,
          a = 3,
          b = 2,
          c = 7
        ),
        'b',
        col = 'black',
        pch = 1)
  lines(p,
        GauPrior(
          p,
          r = 4,
          a = 3,
          b = 2,
          c = 7
        ),
        'b',
        col = 'red',
        pch = 2)
  lines(p,
        GauPrior(
          p,
          r = 1,
          a = 3,
          b = 2,
          c = 5
        ),
        'b',
        col = 'darkgreen',
        pch = 3)
  lines(p,
        GauPrior(
          p,
          r = 4,
          a = 3,
          b = 2,
          c = 5
        ),
        'b',
        col = 'darkblue',
        pch = 4)
  legend(
    'top',
    col = c('black', 'red', 'darkgreen', 'darkblue'),
    lwd = 2,
    legend = c(
      expression("GH(a=3,b=2,c=7,r=1)"),
      expression("GH(a=3,b=2,c=7,r=4)"),
      expression("GH(a=3,b=2,c=5,r=1)"),
      expression("GH(a=3,b=2,c=5,r=4)")
    ),
    pch = c(1, 2, 3, 4),
    cex = 0.9
  )
}
graphics.off()
rm(p)

# testing estimators of rho
n <- 20
p <- c(0.01, 0.1, 0.2, 0.5, 0.7, 0.9, 0.99)
r <- c(1, 2, 3, 4, 5)
a <- 1
b <- 2
c <- 5
tabJefErrors <-
  matrix(nrow = length(p),
         ncol = length(r),
         dimnames = list(p, r))
tabGauErrors <-
  matrix(nrow = length(p),
         ncol = length(r),
         dimnames = list(p, r))
for (i in 1:length(p)) {
  for (j in 1:length(r)) {
    x <- rnbinom(n, r[j], 1 / (1 + p[i] / r[j]))
    tabJefErrors[i, j] <-
      JefSelfRhoEst(x, r[j]) - PostSelfRho(JefPost, x, r[j])
    tabGauErrors[i, j] <- GauSelfRhoEst(x, r[j], a, b, c) -
      PostSelfRho(GauPost, x, r[j], a, b, c)
  }
}
tabJefErrors
tabGauErrors
rm(n, p, r, a, b, c, tabJefErrors, tabGauErrors, x, i, j)

# testing estimators of Lq
set.seed(2022)
n <- 20
p <- c(0.01, 0.1, 0.2, 0.5, 0.7, 0.9, 0.99)
r <- c(1, 2, 3, 4, 5)
a <- 1
b <- 1.5
c <- 4
tabJefErrors <-
  matrix(nrow = length(p),
         ncol = length(r),
         dimnames = list(p, r))
tabGauErrors <-
  matrix(nrow = length(p),
         ncol = length(r),
         dimnames = list(p, r))
for (i in 1:length(p)) {
  for (j in 1:length(r)) {
    x <- rnbinom(n, r[j], 1 / (1 + p[i] / r[j]))
    tabJefErrors[i, j] <-
      JefSelfLqEst(x, r[j]) - PostSelfLq(JefPost, x, r[j])
    tabGauErrors[i, j] <- GauSelfLqEst(x, r[j], a, b, c) -
      PostSelfLq(GauPost, x, r[j], a, b, c)
  }
}
tabJefErrors
tabGauErrors
rm(n, p, r, a, b, c, tabJefErrors, tabGauErrors, x, i, j)

# testing estimators of Ls
{
  set.seed(2022)
  n <- 20
  p <- c(0.01, 0.1, 0.2, 0.5, 0.7, 0.9, 0.99)
  r <- c(1, 2, 3, 4, 5)
  a <- 2
  b <- 3
  c <- 4
  tabJefErrors <-
    matrix(nrow = length(p),
           ncol = length(r),
           dimnames = list(p, r))
  tabGauErrors <-
    matrix(nrow = length(p),
           ncol = length(r),
           dimnames = list(p, r))
  for (i in 1:length(p)) {
    for (j in 1:length(r)) {
      x <- rnbinom(n, r[j], 1 / (1 + p[i] / r[j]))
      tabJefErrors[i, j] <-
        JefSelfLsEst(x, r[j])#-PostSelfLs(JefPost,x,r[j])
      tabGauErrors[i, j] <- GauSelfLsEst(x, r[j], a, b, c)# -
      #PostSelfLs(GauPost,x,r[j],a,b,c)
    }
  }
}
tabJefErrors
tabGauErrors
rm(n, p, r, a, b, c, tabJefErrors, tabGauErrors, x, i, j)

###################################
# testing jeffreys
PostSelfLqJef <- function(x, r) {
  return(PostSelfLq(JefPost, x, r))
}
####################################