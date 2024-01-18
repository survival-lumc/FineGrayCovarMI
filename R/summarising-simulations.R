# Later extent properly with the (...),for now keep it simple
pool_nested_predictions <- function(preds_main,
                                    new_pat, # should be numeric!
                                    ...) {

  bycols <- c(
    "method",
    "rep_id",
    "imp",
    "prob_space",
    "failure_time_model",
    "censoring_type"
  )
  setkeyv(preds_main, bycols)
  df_lp <- preds_main[, .(lp = drop(new_pat %*% coefs[[1]])), by = bycols]

  df_pooled <- preds_main[df_lp][, .(
    pooled_pred = inv_cloglog(mean(cloglog(1 - (1 - base_cuminc)^exp(lp))))
  ), by = c(
    "rep_id",
    "method",
    "prob_space",
    "failure_time_model",
    "censoring_type",
    "time"
  )]
  cbind(df_pooled, "X" = new_pat[1], "Z" = new_pat[2])
}


# Ok do a jackknife ting?
rmse_mcse <- function(estimates, true, K) {

  # Keep first true value
  true_param <- true[1]

  # Calculate elements of rmse mcse
  t_bar <- mean(estimates)
  var_t <- stats::var(estimates)
  t_bar_j <- (1 / (K - 1)) * (K * t_bar - estimates)
  bias_j_sq <- (t_bar_j - true_param)^2
  s_sq_t_j <- (1 / (K - 2)) * ((K - 1) * var_t - (K / (K - 1)) * (estimates - t_bar)^2)
  rmse_j <- sqrt(bias_j_sq + s_sq_t_j)
  mse <- mean((estimates - true_param)^2)

  # Calculate rmse and mcse
  rmse <- sqrt(mse)
  mcse <- sqrt(((K - 1) / (K)) * sum((rmse_j - rmse)^2))

  return(mcse)
}
