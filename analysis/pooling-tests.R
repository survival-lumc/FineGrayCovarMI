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

df_coefs |>
  ggplot(aes(method, estimate, col = method, shape = term)) +
  geom_point(size = 3) +
  facet_grid(failure_time_model * prob_space ~ censoring_type) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  geom_hline(aes(yintercept = true), linetype = "dotted")
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
      rep_id = tar_seed
    )
  ),
  fill = TRUE
)
df_preds[is.na(imp), imp := 0]
df_preds

tar_load(reference_patients)
new_pat <- reference_patients[3, ]

predict_FG <- function(base_cuminc, coefs, new_pat) {
  1 - (1 - base_cuminc)^exp(drop(unlist(coefs) %*% new_pat))
}

predict_FG(0.1, list(c(0, 0)), as.numeric(new_pat))

# Make predictions for one reference patient
df_preds[, pred := predict_FG(
  base_cuminc, coefs, as.numeric(new_pat)
), by = c("method", "rep_id", "time", "imp")]

# Pool
pooled_preds <- df_preds[, .(
  pooled_pred = inv_cloglog(mean(cloglog(pred)))
), by = c("method", "rep_id", "time", "prob_space", "failure_time_model", "censoring_type")]


tar_load(all_true_cuminc)
setDT(all_true_cuminc)
pooled_preds

# Try a tentative merge
new_pat
true_sub <- all_true_cuminc[X == 0 & Z == 1]

merge(pooled_preds, true_sub)[time %in% seq(0, 5, by = 0.5)] |>
  ggplot(aes(time, pooled_pred)) +
  geom_line(aes(col = method, group = method, linetype = method), size = 2, alpha = 0.5) +
  facet_grid(prob_space * failure_time_model  ~ censoring_type) +
  geom_line(aes(y = cuminc), col = "black")
  theme_bw(base_size = 16)
