tar_load(preds_main)

# Get the coefs per imps dataset
coefs_per_imp <- preds_main[!(method %in% c("CCA", "full")), .(
  term = c("X", "Z"),
  coefs = unlist(coefs[[1]])
), by = c("method", "rep_id", "imp", "prob_space", "failure_time_model", "censoring_type")]


# Add MCSE for this and make a ribbon!!
test <- coefs_per_imp[term == "X", .(
  pooled = sapply(seq_len(.N), function(i) mean(coefs[1:i])),
  imp = seq_len(.N)
), by = c(
  "method",
  "rep_id",
  "prob_space",
  "failure_time_model",
  "censoring_type"
)][, .(
  est = mean(pooled),
  emp_se = sd(pooled),
  n = .N,
  mcse_emp_se = sd(pooled) / sqrt(2 * (.N - 1))
), by = c(
  "method",
  "imp",
  "prob_space",
  "failure_time_model",
  "censoring_type"
)]

test[, ':=' (
  emp_lower = emp_se - qnorm(0.975) * mcse_emp_se,
  emp_upper = emp_se + qnorm(0.975) * mcse_emp_se
)]

qnorm(0.975)

0 + c(-1, 1) * qnorm(0.975)

test |>
  ggplot(aes(imp, emp_se, fill = method)) +
  geom_ribbon(aes(ymin = emp_lower, ymax = emp_upper), alpha = 0.15, col = NA) +
  geom_line(linewidth = 1.5, aes(col = method, linetype = method)) +
  facet_grid(prob_space ~ failure_time_model * censoring_type, scales = "free")

test |>
  ggplot(aes(imp, mcse_emp_se, fill = method)) +
  geom_line(linewidth = 1.5, aes(col = method, linetype = method)) +
  facet_grid(prob_space ~ failure_time_model * censoring_type, scales = "free")



# Same with bias ----------------------------------------------------------


coefs_per_imp


test2 <- coefs_per_imp[term == "X" & failure_time_model == "correct_FG", .(
  pooled = sapply(seq_len(.N), function(i) mean(coefs[1:i])),
  imp = seq_len(.N)
), by = c(
  "method",
  "rep_id",
  "prob_space",
  "failure_time_model",
  "censoring_type"
)][, .(
  bias = mean(pooled - 0.75),
  mcse_bias = sd(pooled - 0.75) / sqrt(.N),
  n = .N
), by = c(
  "method",
  "imp",
  "prob_space",
  "failure_time_model",
  "censoring_type"
)]

test2[, ':=' (
  bias_lower = bias - qnorm(0.975) * mcse_bias,
  bias_upper = bias + qnorm(0.975) * mcse_bias
)]

test2 |>
  ggplot(aes(imp, bias, fill = method)) +
  geom_ribbon(aes(ymin = bias_lower, ymax = bias_upper), alpha = 0.15, col = NA) +
  geom_line(linewidth = 1.5, aes(col = method, linetype = method)) +
  facet_grid(prob_space ~ censoring_type, scales = "free")

test2 |>
  ggplot(aes(imp, mcse_bias, fill = method)) +
  geom_line(linewidth = 1.5, aes(col = method, linetype = method)) +
  facet_grid(prob_space ~  censoring_type, scales = "free")

dat <- applied_dat$dat
dat[status_ci_adm == 1][["time_ci_adm"]] |> hist(breaks = 100)

dat[status_ci_adm == 1][order(time_ci_adm)] |> View()

cont_prog_time <- sort(
  table(dat[status_ci_adm == 1]$time_ci_adm),
  decreasing = TRUE
)[1] |> names() |> as.numeric()

dat[, cont_prog := ifelse(
  time_ci_adm == cont_prog_time, "yes", "no"
)]

dat[cont_prog == "yes"]

# Continuous progression time
12 * 28 / 365.25

applied_dat_raw[time_ci_adm == 12 * 28 / 365.25 &
                  time_os_adm <= 12 * 28 / 365.25] |>
  View()

days_fourteen <- 12 * 14 / 365.25
# jitter the tingggg
#12 * 28 / 365.25 + runif()

365.12/12


dat[status_ci_adm == 1, .(max(.N)), by = time_ci_adm]

library(gtsummary)

dat[, !c("id")] |>
  tbl_summary(by = "cont_prog") |>
  add_p()

