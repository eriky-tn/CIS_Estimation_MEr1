################################################################################
# Plotting experiments related to the paper:
# Performance Measure Estimations in Erlang Single Server Queues
#
# Programmed by:
# Eriky S. Gomes & Frederico R. B. Cruz
# Universidade Federal de Minas Gerais
# E-mail: eriky-tn@ufmg.br
# E-mail: fcruz@est.ufmg.br
# (c) 2024 Gomes & Cruz
# v.2024.05.14
################################################################################

library(ggplot2)
library(dplyr)
library(reshape2) # for melt

rm(list = ls())

################################################################################
# Auxiliary functions
################################################################################

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
# basic plotting layout --------------------------------------------------------
################################################################################

# some manual customization variables
{
#cpManual<-c('black','blue','red2','green4','orange3','purple2','navyblue','cyan3',
#            'salmon4','gold','violet','limegreen','springgreen','slateblue4')
cpManual<-c('black','black', 'blue', 'blue','red2','red2')
#ltManual<-c('dotdash','longdash','dashed','dotted','solid')
ltManual<-c('solid','solid','dashed','dashed','dotted','dotted')
shapeManual<-c(15,16,17,18,4,8)

# basic plot
ggBase <- ggplot() + theme_minimal() + labs(color='',linetype='',shape='') +
  theme(
    panel.background = element_rect(fill='white',color='black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.grid.major = element_line(linewidth = 0.5,linetype='solid',color='gray90'),
    #panel.grid.minor = element_line(linewidth = 0.5, linetype='solid',color='gray90'),
    axis.ticks = element_line(),
    axis.ticks.length = unit(c(-0,1),'mm'),
    #legend.box.spacing = unit(0,'pt'),
    legend.margin = margin(t=1,r=1,b=1,l=1,unit='mm')
  ) +
  scale_linetype_manual(values=ltManual) +
  scale_color_manual(values=cpManual) +
  scale_shape_manual(values=shapeManual)
}

plotEst <- function(myData, var_x, var_y, label_x, label_y, file_name,
                    file_ext = '.pdf', x_lim = c(0, 1), y_lim = c(0, 1),
                    leg_position = c(0.97, 0.05), leg_justification = c(0, 1)){
  myPlot <- ggBase +
    xlab(label_x) + ylab(label_y) +
    xlim(x_lim) + ylim(y_lim) +
    #scale_x_continuous(breaks=seq(0,1,0.2)) +
    #scale_y_continuous(breaks=seq(-1,1,0.2)) +
    geom_line(aes_string(x = var_x, y = var_y, color = 'name', linetype = 'name'), myData,
              lwd = 0.3) +
    geom_point(aes_string(x = var_x, y = var_y, color = 'name', shape = 'name'), size = 1.5, myData) +
    # just sets a manual legend
    # annotate('text',
    #          x = 0.03 * (x_lim[2] - x_lim[1]) + x_lim[1],
    #          y = 0.99 * (y_lim[2] - y_lim[1]) + y_lim[1],
    #          label = paste0('r = ', r)) +
    theme(
      legend.position = 'inside',
      legend.position.inside = leg_position,
      legend.justification = leg_justification
    ) +
    guides(color = guide_legend(ncol = 3), shape = guide_legend(ncol = 3),
           linetype = guide_legend(ncol = 3))
  
  # for bias only
  if(var_y == 'bias'){
    myPlot <- myPlot +
      geom_abline(intercept = 0, slope = 0, color = 'black', lwd = 0.2)
  }
  
  #myPlot
  
  # saving to file
    ggsave(
      file = paste0(file_name, file_ext),
      plot = myPlot,
      width = 5,
      heigh = 3
    )
  
  rm(myPlot)
}


################################################################################
# Parameters of the simulation study
################################################################################

Rho_MxM1 <- function(lambda, mu, mean_bat){
  Rho <- lambda * mean_bat / mu
}

Lq_MxM1 <- function(lambda, mu, mean_bat, var_bat){
  mean_bat_squared <- var_bat + mean_bat^2
  Lq <- (lambda * (mean_bat_squared - mean_bat)) / (2 * (mu - lambda * mean_bat))
  return(Lq)
}

Wq_MxM1 <- function(lambda, mu, mean_bat, var_bat){
  Lq <- Lq_MxM1(lambda, mu, mean_bat, var_bat)
  Wq <- Lq / (lambda * mean_bat)
  return(Wq)
}


{
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
}


################################################################################
# organizing data for plotting
################################################################################

{

load('MLE.rdata')
load('Jef.rdata')
load('GH.rdata')

nameMLE1 <- 'MLE r = 1'
nameMLE2 <- 'MLE r = 2'
nameGau1 <- 'GH(2, 3, 5, -1) r = 1'
nameGau2 <- 'GH(2, 3, 5, -1/2) r = 2'
nameJef1 <- 'Jeffreys r = 1'
nameJef2 <- 'Jeffreys r = 2'

# to data-frame r = 1
mle_rho <- cbind(name = nameMLE1, r = 1, mlRhoEst1[[1]])
mle_lq <- cbind(name = nameMLE1, r = 1, mlLqEst1[[1]])
mle_ls <- cbind(name = nameMLE1, r = 1, mlLsEst1[[1]])

gau_rho <- cbind(name = nameGau1, r = 1, gauRhoEst1[[1]])
gau_lq <- cbind(name = nameGau1, r = 1, gauLqEst1[[1]])
gau_ls <- cbind(name = nameGau1, r = 1, gauLsEst1[[1]])

jef_rho <- cbind(name = nameJef1, r = 1, jefRhoEst1[[1]])
jef_lq <- cbind(name = nameJef1, r = 1, jefLqEst1[[1]])
jef_ls <- cbind(name = nameJef1, r = 1, jefLsEst1[[1]])


# to data-frame r = 2
mle_rho <- rbind(mle_rho, cbind(name = nameMLE2, r = 2, mlRhoEst2[[1]]))
mle_lq <- rbind(mle_lq, cbind(name = nameMLE2, r = 2, mlLqEst2[[1]]))
mle_ls <- rbind(mle_ls, cbind(name = nameMLE2, r = 2, mlLsEst2[[1]]))

gau_rho <- rbind(gau_rho, cbind(name = nameGau2, r = 2, gauRhoEst2[[1]]))
gau_lq <- rbind(gau_lq, cbind(name = nameGau2, r = 2, gauLqEst2[[1]]))
gau_ls <- rbind(gau_ls, cbind(name = nameGau2, r = 2, gauLsEst2[[1]]))

jef_rho <- rbind(jef_rho, cbind(name = nameJef2, r = 2, jefRhoEst2[[1]]))
jef_lq <- rbind(jef_lq, cbind(name = nameJef2, r = 2, jefLqEst2[[1]]))
jef_ls <- rbind(jef_ls, cbind(name = nameJef2, r = 2, jefLsEst2[[1]]))


df_rho <- as.data.frame(rbind(mle_rho, gau_rho, jef_rho))
df_lq <- as.data.frame(rbind(mle_lq, gau_lq, jef_lq))
df_ls <- as.data.frame(rbind(mle_ls, gau_ls, jef_ls))

col_numeric <- 2:6

}

################################################################################
# plotting vs parameters
################################################################################

# grouping by parameter
{
myData_rho <-
  df_rho %>%
  mutate_at(col_numeric, as.numeric) %>%
  mutate(bias = mean - rho) %>%
  mutate(rmse = sqrt(bias^2 + sd^2)) %>%
  filter(size != 300 & size != 400) %>%
  group_by(name, rho) %>%
  summarise_at(vars(mean, bias, sd, rmse), mean)

myData_lq <-
  df_lq %>%
  mutate_at(col_numeric, as.numeric) %>%
  mutate(lq = round((r + 1) / (2 * r)  * rho^2 / (1 - rho), 1)) %>%
  mutate(bias = mean - lq) %>%
  mutate(rmse = sqrt(bias^2 + sd^2)) %>%
  filter(size != 300 & size != 400) %>%
  group_by(name, lq) %>%
  summarise_at(vars(mean, bias, sd, rmse), mean)

myData_ls <-
  df_ls %>%
  mutate_at(col_numeric, as.numeric) %>%
  mutate(ls = round(rho + (r + 1) / (2 * r)  * rho^2 / (1 - rho), 1)) %>%
  mutate(bias = mean - ls) %>%
  mutate(rmse = sqrt(bias^2 + sd^2)) %>%
  filter(size != 300 | size != 400) %>%
  group_by(name, ls) %>%
  summarise_at(vars(mean, bias, sd, rmse), mean)
}

# plotting vs parameter
{

variable_plotting <- 'rmse'
label_plotting <- 'RMSE'
leg_position <- c(0, 1.1)
leg_justification <- c(0, 1)

# rho
plotEst(myData_rho, var_x = 'rho', var_y = variable_plotting, label_x = expression(rho),
        label_y = label_plotting, file_name =  paste0(variable_plotting,'_rho'), 
        x_lim = c(0, 1), y_lim = c(0, 0.3), 
        leg_position = leg_position, leg_justification = leg_justification, 
        file_ext = '.eps')

# lq
plotEst(myData_lq, var_x = 'lq', var_y = variable_plotting, label_x = 'Lq',
        label_y = label_plotting, file_name = paste0(variable_plotting,'_lq'), 
        x_lim = c(0.5, 5), y_lim = c(0, 10),
        leg_position = leg_position, leg_justification = leg_justification,
        file_ext = '.eps')
# ls
plotEst(myData_ls, var_x = 'ls', var_y = variable_plotting, label_x = 'Ls',
        label_y = label_plotting, file_name = paste0(variable_plotting,'_ls'), 
        x_lim = c(0.5, 5), y_lim = c(0, 10),
        leg_position = leg_position, leg_justification = leg_justification, 
        file_ext = '.pdf')
}

################################################################################
# plotting vs size
################################################################################


# grouping by size
myData_rho <-
  df_rho %>%
  mutate_at(col_numeric, as.numeric) %>%
  mutate(bias = mean - rho) %>%
  mutate(rmse = sqrt(bias^2 + sd^2)) %>%
  filter(size != 300 & size != 400) %>%
  group_by(name, size) %>%
  summarise_at(vars(mean, bias, sd, rmse), mean)

myData_lq <-
  df_lq %>%
  mutate_at(col_numeric, as.numeric) %>%
  mutate(lq = round((r + 1) / (2 * r)  * rho^2 / (1 - rho), 1)) %>%
  mutate(bias = mean - lq) %>%
  mutate(rmse = sqrt(bias^2 + sd^2)) %>%
  filter(size != 300 & size != 400) %>%
  group_by(name, size) %>%
  summarise_at(vars(mean, bias, sd, rmse), mean)

myData_ls <-
  df_ls %>%
  mutate_at(col_numeric, as.numeric) %>%
  mutate(ls = round(rho + (r + 1) / (2 * r)  * rho^2 / (1 - rho), 1)) %>%
  mutate(bias = mean - ls) %>%
  mutate(rmse = sqrt(bias^2 + sd^2)) %>%
  filter(size != 300 | size != 400) %>%
  group_by(name, size) %>%
  summarise_at(vars(mean, bias, sd, rmse), mean)

# plotting vs size
{

variable_plotting <- 'rmse'
label_plotting <- 'RMSE'
leg_position <- c(0, 1.1)
leg_justification <- c(0, 1)

# rho
plotEst(myData_rho, var_x = 'size', var_y = variable_plotting, label_x = 'n',
        label_y = label_plotting, file_name =  paste0(variable_plotting,'_rho_size'), 
        x_lim = c(0, 200), y_lim = c(0, 0.3), 
        leg_position = leg_position, leg_justification = leg_justification, 
        file_ext = '.eps')
# lq
plotEst(myData_ls, var_x = 'size', var_y = variable_plotting, label_x = 'n',
        label_y = label_plotting, file_name = paste0(variable_plotting,'_lq_size'), 
        x_lim = c(0, 200), y_lim = c(0, 10),
        leg_position = leg_position, leg_justification = leg_justification, 
        file_ext = '.eps')
# ls
plotEst(myData_ls, var_x = 'size', var_y = variable_plotting, label_x = 'n',
        label_y = label_plotting, file_name = paste0(variable_plotting,'_ls_size'), 
        x_lim = c(0, 200), y_lim = c(0, 10),
        leg_position = leg_position, leg_justification = leg_justification,
        file_ext = '.pdf')
}

################################################################################
# plotting prior distributions
################################################################################

# some customization for prior plotting
{
  #cpManual<-c('black','blue','red2','green4','orange3','purple2','navyblue','cyan3',
  #            'salmon4','gold','violet','limegreen','springgreen','slateblue4')
  cpManual<-c('black', 'blue','red2','black','blue', 'red2')
  #ltManual<-c('dotdash','longdash','dashed','dotted','solid')
  ltManual<-c('solid','solid','dashed','dashed','dotted','dotted')
  shapeManual<-c(15,16,17,18,4,8)
  
  # basic plot
  ggBase <- ggplot() + theme_minimal() + labs(color='',linetype='',shape='') +
    theme(
      panel.background = element_rect(fill='white',color='black'),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #panel.grid.major = element_line(linewidth = 0.5,linetype='solid',color='gray90'),
      #panel.grid.minor = element_line(linewidth = 0.5, linetype='solid',color='gray90'),
      axis.ticks = element_line(),
      axis.ticks.length = unit(c(-0,1),'mm'),
      #legend.box.spacing = unit(0,'pt'),
      legend.margin = margin(t=1,r=1,b=1,l=1,unit='mm')
    ) +
    scale_linetype_manual(values=ltManual) +
    scale_color_manual(values=cpManual) +
    scale_shape_manual(values=shapeManual)
}

PlotDensityRho <- function(tabPlot, namePlot, 
                           legend_pos = c(0, 0), legend_justif = c(0, 0),
                           nrowLegend = 2, y_lim = c(0, 0)){
  myPlot <- ggBase + 
    xlab(expression(rho)) + 
    ylab('density') +
    expand_limits(x = c(0, 1)) +
    #scale_x_continuous(breaks=seq(0, 1.1,0.2)) +
    #scale_y_continuous(breaks=seq(-1,1,0.2)) +
    ylim(ifelse(y_lim == c(0, 0), c(min(tabPlot$value), max(tabPlot$value)), y_lim)) +
    geom_line(data = tabPlot, 
              aes(x = p, y = value,
                  color = variable, linetype = variable), lwd = 0.2) +
    geom_point(data = filter(tabPlot, row_number() %% 3 == 1), # gap between points 
               aes(x = p, y = value, 
                   color = variable, shape = variable), size = 1.5) +
    #geom_vline(xintercept = 0.95, lwd = 0.2) +
    theme(legend.position = legend_pos, 
          legend.justification = legend_justif,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(linetype = guide_legend(nrow = nrowLegend))
    
  
  #myPlot
  ggsave(
    file = paste0(namePlot),
    plot = myPlot,
    width = 6,
    height = 4,
    device = cairo_pdf
  )
}


# priors

JefPrior <- function(p, r) {
  return(p ^ (-1 / 2) * (1 + p / r) ^ (-1 / 2))
}

GauHyp <- function(a, b, c, z) {
  stopifnot((c > b) & (b > 0))
  tol = 1e-10
  dGauHyp <- function(u, a, b, c, z) {
    return(u ^ (b - 1) * (1 - u) ^ (c - b - 1) * (1 - z * u) ^ (-a))
  }
  return(integrate(dGauHyp, 0 + tol, 1 - tol, a, b, c, z)[[1]] / beta(b, c - b))
}

GauPrior <- function(p, r, a, b, c) {
  const <- 1 / beta(b, c - b) / GauHyp(a, b , c, -1 / r)
  return(const * p ^ (b - 1) * (1 - p) ^ (c - b - 1) * (1 + p / r) ^ (-a))
}

# plotting priors
p <- seq(1e-3, 1, 1e-2)
a <- c(2, 3, 3)#, 1, -1 / 2, 0)
b <- c(3, 2, 10)#, 1, 2, 5)
c <- c(5, 7, 13)#, 4, 5, 6)

gh1r1 <- GauPrior(p, r = 1, a[1], b[1], c[1])
gh2r1 <- GauPrior(p, r = 1, a[2], b[2], c[2])
gh3r1 <- GauPrior(p, r = 1, a[3], b[3], c[3])
#gh4r1 <- GauPrior(p, r = 1, a[4], b[4], c[4])
#gh5r1 <- GauPrior(p, r = 1, a[5], b[5], c[5])
#gh6r1 <- GauPrior(p, r = 1, a[6], b[6], c[6])

gh1r2 <- GauPrior(p, r = 2, a[1], b[1], c[1])
gh2r2 <- GauPrior(p, r = 2, a[2], b[2], c[2])
gh3r2 <- GauPrior(p, r = 2, a[3], b[3], c[3])
#gh4r2 <- GauPrior(p, r = 2, a[4], b[4], c[4])
#gh5r2 <- GauPrior(p, r = 2, a[5], b[5], c[5])
#gh6r2 <- GauPrior(p, r = 2, a[6], b[6], c[6])

jr1 <- JefPrior(p, r = 1)
jr2 <- JefPrior(p, r = 2)

ghp <- cbind(p, 
  gh1r1, gh2r1, gh3r1,# gh4r1, gh5r1, gh6r1,
  gh1r2, gh2r2, gh3r2#, gh4r2, gh5r2, gh6r2
  ) 
jp <- cbind(p, jr1, jr2)


names_ghp <- c(paste0('GH(', a, ', ', b, ', ', c, ', ', '-1)'),
                   paste0('GH(', a, ', ', b, ', ', c, ', ', '-1/2)'))
names_jp <- c('jeffreys r = 1', 'jeffreys r = 2')
colnames(ghp) <- c('p', names_ghp)
colnames(jp) <- c('p', names_jp)

ghp <- as.data.frame(ghp)
jp <- as.data.frame(jp)
ghp <- melt(ghp, id.vars = 'p')
jp <- melt(jp, id.vars = 'p')

PlotDensityRho(ghp, 'gh_prior.eps', y_lim = c(0, 4),
               legend_pos = c(0, 1.05), legend_justif = c(0, 1))
PlotDensityRho(jp, 'jef_prior.eps',
               legend_pos = c(0.97, 1.05), legend_justif = c(1, 1))

