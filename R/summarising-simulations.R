# Later extent properly with the (...),for now keep it simple
pool_nested_predictions <- function(nested_df,
                                    new_pat, # should be numeric!
                                    ...) {

  # Collapse the nested data frame with necessary identifier
  unnested_df <- rbindlist(
    with(
      tar_read(simulations_main),
      #nested_df,
      Map(
        f = cbind,
        method = method,
        preds_summary,
        prob_space = prob_space,
        failure_time_model = failure_time_model,
        censoring_type = censoring_type,
        rep_id = tar_rep,
        batch_id = tar_batch#,
        #...
      )
    ),
    fill = TRUE
  )

  # Give full and CCA methods a 'zero' id for the imputation
  unnested_df[is.na(imp), imp := 0]

  # To make prediction
  predict_nested_FG <- function(base_cuminc, coefs, new_pat) {
    1 - (1 - base_cuminc)^exp(coefs[[1]] %*% new_pat)
  }

  # Predict for a given patient
  unnested_df[, pred := predict_nested_FG(base_cuminc, coefs, new_pat), by = c(
    "method",
    "rep_id",
    "batch_id",
    "time",
    "imp",
    "prob_space",
    "failure_time_model",
    "censoring_type"
  )]

  # Now we pool the predictions after cloglog transformations
  pooled_preds <- unnested_df[, .(pooled_pred = inv_cloglog(mean(cloglog(pred)))), by = c(
    "method",
    "rep_id",
    "batch_id",
    "time",
    "prob_space",
    "failure_time_model",
    "censoring_type"
  )]

  return(pooled_preds)
}


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
