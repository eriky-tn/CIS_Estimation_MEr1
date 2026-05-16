################################################################################
# Tables for experiments related to the paper:
# Performance Measure Estimations in Erlang Single Server Queues
#
# Programmed by:
# Eriky S. Gomes & Frederico R. B. Cruz
# Universidade Federal de Minas Gerais
# E-mail: eriky-tn@ufmg.br
# E-mail: fcruz@est.ufmg.br
# (c) 2024 Gomes & Cruz
# v.2024.05.11
################################################################################

library(dplyr)
library(reshape2) # for melt
library(tidyr) # convert long to wide tables format


rm(list = ls())

################################################################################
# Auxiliary functions
################################################################################

# Format to Excel
FrExcel <- function(data){
  write.table(format(data, decimal.mark = ','),
              'clipboard', sep='\t', na = '*')
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

size = c(10, 20, 50, 100, 200)
r = c(1, 2)
a <- 2
b <- 3
c <- 5

################################################################################
# organizing data
################################################################################

# for r = 1

r = 1

load('MLE.rdata')
load('Jef.rdata')
load('GH.rdata')

nameMLE <- paste0('MLE r = ', r)
nameGau <- paste0('GH(2, 3, 5, -1/', r,') r = ', r)
nameJef <- paste0('Jeffreys  r = ', r)

# to data-frame
mle_rho <- cbind(name = nameMLE, mlRhoEst1[[1]])
mle_lq <- cbind(name = nameMLE, mlLqEst1[[1]])
mle_ls <- cbind(name = nameMLE, mlLsEst1[[1]])

gau_rho <- cbind(name = nameGau, gauRhoEst1[[1]])
gau_lq <- cbind(name = nameGau, gauLqEst1[[1]])
gau_ls <- cbind(name = nameGau, gauLsEst1[[1]])

jef_rho <- cbind(name = nameJef, jefRhoEst1[[1]])
jef_lq <- cbind(name = nameJef, jefLqEst1[[1]])
jef_ls <- cbind(name = nameJef, jefLsEst1[[1]])

df_rho <- rbind(mle_rho, gau_rho, jef_rho) %>% as.data.frame()
df_lq <- rbind(mle_lq, gau_lq, jef_lq) %>% as.data.frame()
df_ls <- rbind(mle_ls, gau_ls, jef_ls) %>% as.data.frame()


################################################################################
# generating tables
################################################################################

col_numeric = 2:5
col_order = c(1, 2, 3, 8, 4, 9, 5, 10, 6, 11, 7, 12) # the order is mixed on pivot_wider

tabRho <- 
  df_rho %>%
  mutate_at(col_numeric, as.numeric) %>%
  filter(size != 300 & size != 400) %>%
  mutate(var = sd^2) %>%
  select(name, rho, size, mean, var) %>%
  pivot_wider(names_from = size, values_from =  c(mean, var)) %>%
  select(all_of(col_order)) %>%
  as.data.frame

tabLq <- 
  df_lq %>%
  filter(size != 300 & size != 400) %>%
  mutate_at(col_numeric, as.numeric)%>%
  mutate(lq = round((r + 1) / (2 * r)  * rho^2 / (1 - rho), 1)) %>%
  mutate(var = sd^2) %>%
  select(name, lq, size, mean, var) %>%
  pivot_wider(names_from = size, values_from =  c(mean, var)) %>%
  select(all_of(col_order)) %>%
  as.data.frame

tabLs <- 
  df_ls %>%
  filter(size != 300 & size != 400) %>%
  mutate_at(col_numeric, as.numeric) %>%
  mutate(ls = round(rho + (r + 1) / (2 * r)  * rho^2 / (1 - rho), 1)) %>%
  mutate(var = sd^2) %>%
  select(name, ls, size, mean, var) %>%
  pivot_wider(names_from = size, values_from =  c(mean, var)) %>%
  select(all_of(col_order)) %>%
  as.data.frame

################################################################################
# tables from clipboard to excel
################################################################################

FrExcel(tabRho)
FrExcel(tabLq)
FrExcel(tabLs)

################################################################################
# r = 2
################################################################################
# organizing data
################################################################################

# clear all variables that is not a function
{
  variables <- ls()
  functions <- variables[sapply(variables, function(x) is.function(get(x)))]
  rm(list = setdiff(variables, functions))
  rm(functions, variables)
}


# for r = 2

r = 2

load('MLE.rdata')
load('Jef.rdata')
load('GH.rdata')

nameMLE <- paste0('MLE r = ', r)
nameGau <- paste0('GH(2, 3, 5, -1/', r,') r = ', r)
nameJef <- paste0('Jeffreys  r = ', r)

# to data-frame
mle_rho <- cbind(name = nameMLE, mlRhoEst2[[1]])
mle_lq <- cbind(name = nameMLE, mlLqEst2[[1]])
mle_ls <- cbind(name = nameMLE, mlLsEst2[[1]])

gau_rho <- cbind(name = nameGau, gauRhoEst2[[1]])
gau_lq <- cbind(name = nameGau, gauLqEst2[[1]])
gau_ls <- cbind(name = nameGau, gauLsEst2[[1]])

jef_rho <- cbind(name = nameJef, jefRhoEst2[[1]])
jef_lq <- cbind(name = nameJef, jefLqEst2[[1]])
jef_ls <- cbind(name = nameJef, jefLsEst2[[1]])

df_rho <- rbind(mle_rho, gau_rho, jef_rho) %>% as.data.frame()
df_lq <- rbind(mle_lq, gau_lq, jef_lq) %>% as.data.frame()
df_ls <- rbind(mle_ls, gau_ls, jef_ls) %>% as.data.frame()


################################################################################
# generating tables
################################################################################

col_numeric = 2:5
col_order = c(1, 2, 3, 8, 4, 9, 5, 10, 6, 11, 7, 12) # the order is mixed on pivot_wider

tabRho <- 
  df_rho %>%
  mutate_at(col_numeric, as.numeric) %>%
  filter(size != 300 & size != 400) %>%
  mutate(var = sd^2) %>%
  select(name, rho, size, mean, var) %>%
  pivot_wider(names_from = size, values_from =  c(mean, var)) %>%
  select(all_of(col_order)) %>%
  as.data.frame

tabLq <- 
  df_lq %>%
  filter(size != 300 & size != 400) %>%
  mutate_at(col_numeric, as.numeric)%>%
  mutate(lq = round((r + 1) / (2 * r)  * rho^2 / (1 - rho), 1)) %>%
  mutate(var = sd^2) %>%
  select(name, lq, size, mean, var) %>%
  pivot_wider(names_from = size, values_from =  c(mean, var)) %>%
  select(all_of(col_order)) %>%
  as.data.frame

tabLs <- 
  df_ls %>%
  filter(size != 300 & size != 400) %>%
  mutate_at(col_numeric, as.numeric) %>%
  mutate(ls = round(rho + (r + 1) / (2 * r)  * rho^2 / (1 - rho), 1)) %>%
  mutate(var = sd^2) %>%
  select(name, ls, size, mean, var) %>%
  pivot_wider(names_from = size, values_from =  c(mean, var)) %>%
  select(all_of(col_order)) %>%
  as.data.frame

################################################################################
# tables from clipboard to excel
################################################################################

FrExcel(tabRho)
FrExcel(tabLq)
FrExcel(tabLs)
