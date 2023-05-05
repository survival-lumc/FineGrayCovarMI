tar_load(all_simulations)

# Example for multiple
#lapply(list(c(1, 1), c(1, 2), c(0.5, 1)), function(x) {
#  essentials[, .(drop(unlist(coefs) %*% x)), by = c("time", "base_cuminc")]
#})

# Coefficients
df_coefs <- rbindlist(
  with(
    all_simulations,
    Map(
      cbind,
      method = method,
      coefs_summary,
      prob_space = prob_space,
      failure_time_model = failure_time_model,
      censoring_type = censoring_type
    )
  )
)


df_coefs[, term := ifelse(grepl(pattern = "^X", term), "X", as.character(term))]

sim_summ_X <- simsum(
  data = df_coefs[term == "X"],
  estvarname = "estimate",
  se = "std.error",
  true = "true",
  methodvar = "method",
  ref = "full",
  by = c("censoring_type", "failure_time_model", "prob_space")
)

summary(sim_summ_X)
summ <- data.table(sim_summ_X$summ)
summ[stat == "bias"][["mcse"]] |> hist(breaks = 20, xlim = c(0, 0.02))
summ[stat == "empse"][["est"]] |> hist(breaks = 20, xlim  = c(0, 0.2))

df_coefs[term == "X"] |>
  ggplot(aes(method, estimate - true)) +
  geom_jitter(aes(col = method), size = 2.5, width = 0.25, alpha = 0.5, shape = 16) +
  facet_grid(failure_time_model * prob_space ~ censoring_type) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  geom_hline(aes(yintercept = 0), #true),
             linetype = "dashed", size = 1) +
  stat_summary(
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "crossbar",
    width = 0.75,
    aes(col = method)
  ) +
  scale_color_manual(values = Manu::get_pal("Hoiho"))

# Then Z
df_coefs[term == "Z"] |>
  ggplot(aes(method, estimate, col = method)) +
  geom_jitter(size = 2.5, width = 0.25, alpha = 0.8) +
  facet_grid(failure_time_model * prob_space ~ censoring_type) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  geom_hline(aes(yintercept = true), linetype = "dotted") +
  stat_summary(
    fun = mean,
    fun.min = mean,
    fun.max = mean,
    geom = "crossbar",
    width = 0.75,
    color = "black"
  )

# pspace thing can be edited later..


# Predictions -------------------------------------------------------------


# For now use tar_seed as rep_id
df_preds <- rbindlist(
  with(
    all_simulations,
    Map(
      cbind,
      method = method,
      preds_summary,
      prob_space = prob_space,
      failure_time_model = failure_time_model,
      censoring_type = censoring_type,
      rep_id = tar_rep
    )
  ),
  fill = TRUE
)
df_preds[is.na(imp), imp := 0]
df_preds

tar_load(reference_patients)
new_pat <- reference_patients[1, ]
new_pat
rm(all_simulations)

#anyNA(df_preds$coefs) # no na's

predict_FG <- function(base_cuminc, coefs, new_pat) {
  1 - (1 - base_cuminc)^exp(drop(unlist(coefs) %*% new_pat))
}

#predict_FG(0.1, list(c(0, 0)), as.numeric(new_pat))

#df_sub <- df_preds[imp %in% c(0, 1, 2)]
#rm(df_preds)

# Make predictions for one reference patient
df_preds[, pred := predict_FG(
  base_cuminc, coefs, as.numeric(new_pat)
), by = c(
  "method",
  "rep_id",
  "time",
  "imp",
  "prob_space",
  "failure_time_model",
  "censoring_type"
)]

# Pool
pooled_preds <- df_preds[, .(
  pooled_pred = inv_cloglog(mean(cloglog(pred)))
), by = c("method", "rep_id", "time", "prob_space", "failure_time_model", "censoring_type")]

rm(df_preds)

tar_load(all_true_cuminc)
setDT(all_true_cuminc)
pooled_preds

# Now we need to average predictions a la monte carlo
avg_preds <- pooled_preds[, .(
  avg_pred = mean(pooled_pred)
), by = c("method", "time", "prob_space", "failure_time_model", "censoring_type")]

# Try a tentative merge
new_pat
true_sub <- all_true_cuminc[X == 0 & Z == 0]

merge(avg_preds, true_sub) |> #[time %in% seq(0, 5, by = 0.5)] |>
  ggplot(aes(time, avg_pred - cuminc)) +
  geom_line(aes(col = method, group = method, linetype = method), size = 2, alpha = 0.5) +
  #geom_line(aes(y = cuminc), col = "black", size = 1.25) +
  facet_grid(failure_time_model * prob_space  ~ censoring_type, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 16) #+
 # coord_cartesian(ylim = c(0, 0.5))

# Let's zoom in one method
# Compare to CCA!!
merge(pooled_preds, true_sub)[method %in% "smcfcs_finegray"] |>
  ggplot(aes(time, pooled_pred)) +
  geom_line(aes(group = rep_id), size = 1, col = Manu::get_pal("Hoiho")[2], alpha = 0.75) +
  facet_grid(failure_time_model * prob_space  ~ censoring_type, scales = "free") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")
