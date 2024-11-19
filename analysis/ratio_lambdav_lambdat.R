dat[, V := (D == 1) * time + (D != 1) * cens_time]
dat[, D_star := as.numeric(D == 1)]
dat[, Lambda_V := compute_marginal_cumhaz(
  timevar = V,
  statusvar = D_star,
  cause = 1,
  type = "cause_spec"
)]
dat[, Lambda_T := compute_marginal_cumhaz(
  timevar = time,
  statusvar = D,
  cause = 1,
  type = "subdist"
)]

dat_sub <- dat[D != 1]

plot(dat_sub$time, log(dat_sub$Lambda_V / dat_sub$Lambda_T), type = "p")
