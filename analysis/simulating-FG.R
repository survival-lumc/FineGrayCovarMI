# Load libraries
library(data.table)
library(crr)
library(ggplot2)
library(magrittr)
library(mstate) # for the crprep approach
library(prodlim)
library(riskRegression)

# https://github.com/erickawaguchi/fastcmprsk/blob/master/R/simulateTwoCauseFineGrayModel.R
# https://www.uniklinik-freiburg.de/imbi/stud-le/competing-risks-and-multistate-models-with-r.html

# See Zotero and Geskus slides + beyersmann book; 'simulation mass' at infinity



# Indirect method ---------------------------------------------------------


# Clean data.table way
n <- 5000
p <- 0.4

dat <- data.table(X = rbinom(n, size = 1, 0.5), U = runif(n))
dat[, D := 1 + rbinom(n, size = 1, prob = (1 - p)^exp(0.5 * X))]
dat[D == 2, time := rexp(.N, rate = exp(0.25 * X))]
dat[D == 1, time := -log(
  1 - (1 - (1 - U * (1 - (1 - p)^exp(0.5 * X)))^(1 / exp(0.5 * X))) / p
)]

mod <- FGR(Hist(time, D) ~ X, cause = 1, data = dat)
mod$crrFit$coef


# Using cause-specific hazards --------------------------------------------



#... make plots; write-up in rmarkdown
