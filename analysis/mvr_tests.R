set.seed(488641)

library(survival)
library(mice)
library(MASS)
library(mvtnorm)
library(ggforestplot)
library(ggplot2)

source("R/data-generation.R")

generate_complete_data <- function(n = 10000,
                                   p = 5,
                                   rho = 0.5,
                                   mu = rep(0, p),
                                   sds = diag(rep(1, p)),
                                   betas = rep(0.5, p),
                                   base_rate = 0.1,
                                   max_unif_cens = 10,
                                   dichot = FALSE) {


    # Generate covariate matrix
  cormat <- diag(p)
  cormat[cormat != 1] <- rho
  covmat <- sds %*% cormat %*% sds
  X <- MASS::mvrnorm(n = n, mu = mu, Sigma = covmat)
  colnames(X) <- paste0("X", seq_len(p))

  if (dichot) X <- matrix(as.numeric(X > 0), ncol = p)

  # Now we generate event times
  lp <- drop(betas %*% t(X))
  time_tilde <- rexp(n, rate = base_rate * exp(lp))
  cens <- runif(n, min = 0.001, max = max_unif_cens)
  D <- as.numeric(time_tilde < cens)
  time <- ifelse(D == 1, time_tilde, cens)

  # Lets compute f(T, D | X1, ... Xp)
  outcome_dens <- (D == 1) * (base_rate * exp(lp)) + (D == 0) * exp(-time * exp(lp))

  # Return complete data
  dat <- cbind.data.frame(time, D, X)

  # Now for densities
  densities_X <- densities_X <- matrix(
    data = 0,
    nrow = n,
    ncol = p,
    dimnames = list(seq_len(n), colnames(X))
  )

  # Compute joint covariate density, and also cond. dens of any given X given the others
  # (same for all, since correlations are the same)
  joint_dens_X <- mvtnorm::dmvnorm(x = X, mean = mu, sigma = covmat)

  # https://en.wikipedia.org/wiki/Multivariate_normal_distribution for marginal properties
  # Example with f(X_j | X_-j), and j = 1
  #j <- 1
  #jointmarg_dens_Xminj <- dmvnorm(x = X[, -j], mean = mu[-j], sigma = covmat[-j, -j])
  #cond_dens_X <- mv_dens_X / jointmarg_dens_Xminj

  for (j in seq_len(p)) {
    jointmarg_dens_Xminj <- dmvnorm(x = X[, -j], mean = mu[-j], sigma = covmat[-j, -j])
    densities_X[, j] <- joint_dens_X / jointmarg_dens_Xminj
  }
  dat_densities <- cbind.data.frame(outcome_dens, joint_dens_X, densities_X)

  # Ignore the densities if binary
  list("dat" = dat, "densities" = dat_densities)
}

add_missings <- function(dat,
                         which_miss = NULL,
                         prop_missing = 0.4) {

  if (is.null(which_miss)) colnames(dat)[grepl(colnames(dat), pattern = "^X")]

  # Ampute missing values
  amp_zero <- ampute(dat)

  # Do not impute outcome
  patts <- amp_zero$patterns
  expr_patt <- paste(paste0(c("time", "D", which_miss), " == 1"), collapse = " & ")
  inds <- with(patts, eval(parse(text = expr_patt)))
  patts <- patts[inds, ]

  # Missingness indepent of outcome
  wts <- amp_zero$weights[inds, ]
  wts[, c("time", "D")] <- 0

  ampute(dat, prop = prop_missing, patterns = patts, weights = wts)$amp
}

run_imps <- function(dat,
                     sm_form, # substantive model formula
                     m = 5,
                     cycles = 5) {

  # Add cumulative hazard
  dat$cumhaz <- compute_marginal_cumhaz(
    dat$time,
    dat$D,
    cause = 1,
    type = "cause_spec"
  )

  # Substantive model predictors
  sm_predictors <- all.vars(update(sm_form, 1 ~ .))

  # Prepare methods and predictor matrix
  meths_mice <- make.method(dat, defaultMethod = c("norm", "logreg", "polyreg", "polr"))

  # Two approaches: include or omit outcome
  predmat_incl <- predmat_excl <- make.predictorMatrix(dat)
  predmat_incl[] <- predmat_excl[] <- 0
  predmat_incl[, c(sm_predictors, "D", "cumhaz")] <- 1
  predmat_excl[, c(sm_predictors)] <- 1

  # Run imps
  imps_incl <- mice(
    data = dat,
    m = m,
    maxit = cycles,
    method = meths_mice,
    predictorMatrix = predmat_incl
  )

  imps_excl <- mice(
    data = dat,
    m = m,
    maxit = cycles,
    method = meths_mice,
    predictorMatrix = predmat_excl
  )

  # Pool and summarise already
  mods_incl <- lapply(
    complete(imps_incl, "all"),
    function(impdat) {
      mod <- do.call(what = survival::coxph, args = list(formula = sm_form, data = impdat))
    }
  )

  mods_excl <- lapply(
    complete(imps_excl, "all"),
    function(impdat) {
      mod <- do.call(what = survival::coxph, args = list(formula = sm_form, data = impdat))
    }
  )

  # Return, including CCA
  mod_cca <- coxph(sm_form, data = dat)

  dplyr::bind_rows(
    cbind(broom::tidy(mod_cca, conf.int = TRUE), method = "CCA"),
    cbind(broom::tidy(pool(mods_incl), conf.int = TRUE), method = "X, D, H(T)"),
    cbind(broom::tidy(pool(mods_excl), conf.int = TRUE), method = "Just X")
  )
}




#mice::md.pattern(dat_amp)
# make an arg on which/how many of variables have missings

# Check mechanism
#glm(is.na(X1) ~ ., data = dat_amp, family = "binomial")

p <- 5
sm_form <- reformulate(response = "Surv(time, D)", termlabels = paste0("X", seq_len(p)))
betas <- rep(0.3, p) #seq(0.25, 0.5, length.out = p) #rep(1, p)
dichot <- TRUE
dat_gen <- generate_complete_data(p = p, rho = 0.1, betas = betas, dichot = dichot)
dat_full <- dat_gen$dat
dat_gen$densities |> head()

dat_amp <- add_missings(dat_full)

if (dichot) {
  X_inds <- grepl(colnames(dat_amp), pattern = "^X")
  dat_amp[, X_inds] <- lapply(dat_amp[, X_inds], as.factor)
}

res <- run_imps(dat_amp, sm_form = sm_form, m = 5, cycles = 10)


forestplot(
  df = res,
  name = term,
  estimate = estimate,
  se = std.error,
  xtickbreaks = c(0.5, betas), #c(0.25, 0.5, 0.75, betas, 1.25),
  colour = method,
  xlim = c(0, 1.5)
)



library(ggridges)

dat_gen$densities |>
  tidyr::pivot_longer(
    cols = paste0("X", seq_len(5)),
    names_to = "term",
    values_to = "cond_dens"
  ) |>
  ggplot(aes(y = term, x = log(outcome_dens / cond_dens))) +
  #facet_wrap(~ term) +
  geom_density_ridges(alpha=0.6, bandwidth=4)


hist(log(dat_gen$densities$outcome_dens))
hist(log(dat_gen$densities$joint_dens_X))

dat_gen$densities |>
  ggplot(aes(log(joint_dens_X))) +
  geom_density() +
  geom_density(aes(log(outcome_dens)), col = "blue")

