dat_p3 <- dat_true[, .(
  HR_cs = cs_haz[X == 1] / cs_haz[X == 0],
  HR_subdist_cum = log(1 - cuminc[X == 1]) / log(1 - cuminc[X == 0]),
  HR_subdist = subdist_haz[X == 1] / subdist_haz[X == 0]
  # Use Grambauer eq. (7) instead?
), by = c("time", "cause", "prob_space", "failure_time_model")]

dat_TFP <- dat_true[, .(
  TFP = sum(cuminc),
  p_not_k = 1 - cuminc,
  cause = cause
), by = c("X", "time", "prob_space", "failure_time_model")][, .(
  ratio_Fk = p_not_k[X == 0] / p_not_k[X == 1],
  ratio_S = (1 - TFP[X == 1]) / (1 - TFP[X == 0])
), by = c("time", "prob_space", "failure_time_model", "cause")]

dat_merged <- merge(dat_TFP[cause == 1], dat_p3[cause == 1])
dat_merged
dat_merged[, HR_grambauer := HR_cs * ratio_Fk * ratio_S]

dat_merged |>
  melt(
    variable.name = "HR_type",
    value.name = "HR",
    measure.vars = c("HR_cs", "HR_subdist_cum", "HR_subdist", "HR_grambauer")
  ) |>
  ggplot(aes(time, log(HR), col = HR_type, linetype = HR_type)) +
  geom_line(linewidth = 1, aes(group = HR_type)) +
  facet_wrap(~ failure_time_model * prob_space) +
  scale_color_manual(values = 1:4)

## Grambauer corresponds to subdist ratios!!
