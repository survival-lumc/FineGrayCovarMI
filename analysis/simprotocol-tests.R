# Tests:
# - Get actual censoring rates (for unif and exp)
# - Test out different pred times
# - Check nested is still not defeated by memory, otherwise always return one
# ..long data.table
# - Adjust methods for time when you do know the eventual censoring time
# - Figure out how to store/call true things after?
# - Continue with prediction tings
# - specific censoring rate for specific p?


args_event_times <- list(
  mechanism = "correct_FG",
  params = tar_read(true_params_correct_FG_0.15),
  censoring_type = "curvy_uniform",
  censoring_params = list("exponential" = 0.45, "curvy_uniform" = c(0.5, 5))
)
args_missingness <- list(mech_params = list("prob_missing" = 0.4, "mechanism_expr" = "Z"))
args_imputations <- list(m = 2, iters = 2, rjlimit = 1000)
args_predictions <- list(timepoints = seq(1, 5, by = 1))


true_betas <- c(0.75, 0.5)

params <- tar_read(true_params_correct_FG_0.25)
tar_load(pred_timepoints)
tar_load(all_sims)
all_sims[1, "preds_summary"][[1]]

reference_patients

compute_true_cuminc(
  t = pred_timepoints,
  newdat = reference_patients[4, ],
  params = tar_read(params_weibull_lfps_0.75), #tar_read(true_params_correct_FG_0.25), # tar_read(params_weibull_lfps_0.25)
  model_type = "misspec_FG" #"correct_FG" #"misspec_FG"
)


# Censoring rates test ----------------------------------------------------


dat <- generate_dataset(
  n = 1e3,
  list(
    mechanism = "correct_FG",
    params = tar_read(true_params_correct_FG_0.25),
    censoring_type = "exponential",
    censoring_params = list("exponential" = 0.45, "uniform" = c(0, 3))
  ),
  list(mech_params = list("prob_missing" = 0, "mechanism_expr" = "1"))
)

# Uniform might need to be much shorter..
# .. and need to change code so it work in case of admin censoring
dat
prop.table(table(dat$D))
plot(survfit(Surv(time, D) ~  1, data = dat))
prodlim(Hist(time, factor(D)) ~ 1, data = dat) |> plot()

FGR(Hist(time, D) ~ X +Z, data = dat, cause = 1)

dat[, ':=' (
  time_new = ifelse(D == 2, cens_time, time),
  D_new = 1 + as.numeric(D != 1)
)]

FGR(Hist(time, D) ~ X +Z, data = dat, cause = 1)$crrFit$coef
FGR(Hist(time_new, D_new) ~ X +Z, data = dat, cause = 1)$crrFit$coef
coxph(Surv(time_new, D_new == 1)~ X + Z, data = dat) |> coef()

# Visualise missing mechanism with jitter plots!!

# Pooling tests -----------------------------------------------------------


tar_load(all_sims)
pryr::object_size(all_sims) * 500
# https://stackoverflow.com/questions/70878796/how-to-store-nested-data-efficiently-in-r ??

# Bind coefs togeths
df_coefs <- rbindlist(
  with(
    all_sims,
    Map(
      cbind,
      method = method,
      coefs_summary,
      prob_space = prob_space,
      failure_time_model = failure_time_model,
      censoring_type = censoring_type#,
      #scen_summary[[1]]$params$cause1$betas
    )
  ),
  fill = TRUE
)
pryr::object_size(df_coefs)

# Continue here...

# Otherwise: don't nest, and store everything in one big df?
# Fill NAs with the rest?

# Bind predictions
df_preds <- rbindlist(
  with(all_sims, Map(cbind, method = method, sim_rep = sim_rep, preds_summary)),
  fill = TRUE
)
df_preds[is.na(imp), imp := 0]
df_preds

# Try test pooling
new_pat <- c("X" = 0.5, "Z" = 1)
new_pats <- list(
  "A" = c("X" = 0.5, "Z" = 1),
  "B" = c("X" = 0.5, "Z" = 0.5)
)
testo <- df_preds[, .(
  pred = 1 - (1 - base_cuminc)^exp(drop(unlist(coefs) %*% new_pat))
), by = c("method", "sim_rep", "time", "imp")]

df_preds[, .(
  lapply(new_pats, function(x) {
    1 - (1 - base_cuminc)^exp(drop(unlist(coefs) %*% x))
  })
), by = c("method", "sim_rep", "time", "imp")]

#for (j in cols) set(dt, j = j, value = -dt[[j]])
# for predictions?

# See also https://stackoverflow.com/questions/16846380/apply-a-function-to-every-specified-column-in-a-data-table-and-update-by-referen
# or just rbindlist the very big ones..
testo



testo[, .(
  pooled_pred = inv_cloglog(mean(cloglog(pred)))
), by = c("method", "sim_rep", "time")] |>
  ggplot(aes(time, pooled_pred)) +
  geom_line(aes(col = method, linetype = method), size = 1) +
  facet_wrap(~ sim_rep) +
  theme_bw()

essentials[, .(drop(unlist(coefs) %*% x)), by = c("time", "base_cuminc")]

#df[, unlist(summary), by = c("method", "sim_rep")]

# file:///C:/Users/efbonneville/Downloads/Manuscript%20(1).pdf
df_coefs <- df[, .(summary_coefs = list(summary[[1]][[2]])), by = c("method", "sim_rep")]
df_preds <- df[, .(summary_preds = list(summary[[1]][[1]])), by = c("method", "sim_rep")]

# Save the preds as a df!


df_coefs[, cbind(method, rbindlist(summary_coefs, fill = TRUE)), by = method]
rbindlist(df_coefs[, cbind(method, sim_rep, summary_coefs)])
df_coefs[, .(list(cbind(method, sim_rep, summary_coefs))), by = c("method", "sim_rep")]

test <- df_coefs[1, ]
test[, .(list(cbind(method, sim_rep, summary_coefs[[1]])))]

rbindlist(Map(cbind, timestamp = df_coefs$method, df_coefs$summary_coefs), fill = TRUE)

# Thank the lawddddddd https://stackoverflow.com/questions/58563899/dataframe-with-nested-dataframes-how-to-set-id-column-with-data-tablerbindli
rbindlist(df_coefs[, .(Map(cbind, method, sim_rep, summary_coefs))], fill = TRUE)

df_coefs[, .(
  list(cbind(method, sim_rep, summary_coefs[[1]]))
), by = c("method", "sim_rep")]



tidyr::unnest(df_coefs, summary_coefs)
tidyr::unnest(df_preds, summary_preds)



df_coefs[, unlist(summary_coefs, recursive = FALSE), c("method", "sim_rep")]


df[, .(dasummary[[1]][[2]]), by = c("method", "sim_rep")]
rbindlist()
df[, .(pooled_info = summary[[1]][[2]]), by = c("method", "sim_rep")]
pryr::object_size(df)
pryr::mem_used()
# https://win-vector.com/2014/05/30/trimming-the-fat-from-glm-models-in-r/


# Use to get baseline cumulative incidence at grid of monthly timepoints
1 - (1 - predict(mod_full, newdata = list(X_obs = 0, Z = 0),
                 times = c(3, 4, 5)))^exp(mod_full$crrFit$coef[1])


# Something else ----------------------------------------------------------






