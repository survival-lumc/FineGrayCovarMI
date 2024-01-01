# Possible to-do/ideas:
# - Add proposed smcfcs.finegray to smcfcs() package; add options to kmi settings?
# - Visualise missing mechanism with jitter plots! (+ visualise the settings, together with base cumincs)
# - See {mets} for cifreg() and doubleFGR()

coefs_main <- rbindlist(
  with(
    tar_read(simulations_main),
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

coefs_main[, term := ifelse(grepl(pattern = "^X", term), "X", as.character(term))]
coefs_main[, method := factor(
  method,
  levels = c("full", "CCA", "mice_comp", "mice_subdist", "smcfcs_comp", "smcfcs_finegray"),
  labels = c(
    "Full data",
    "Compl. cases",
    "MICE cause-spec",
    "MICE subdist",
    "SMC-FCS cause-spec",
    "SMC-FCS Fine-Gray"
  )
)]
coefs_main[, censoring_type := factor(
  censoring_type,
  levels = c("none", "exponential", "curvy_uniform")
)]

coefs_main[failure_time_model == "misspec_FG"] |>
  ggplot(aes(method, estimate - true)) +
  geom_boxplot(aes(fill = method)) +
  #geom_jitter(aes(col = method), size = 2.5, width = 0.25, alpha = 0.25, shape = 16) +
  facet_grid(prob_space ~ term * censoring_type) +
  geom_hline(yintercept = 0)

coefs_main[failure_time_model == "correct_FG"] |>
  ggplot(aes(method, std.error)) +
  geom_boxplot(aes(fill = method)) +
  #geom_jitter(aes(col = method), size = 2.5, width = 0.25, alpha = 0.25, shape = 16) +
  facet_grid(prob_space ~ term * censoring_type) +
  geom_hline(yintercept = 0)
