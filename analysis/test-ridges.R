library(ggridges)

df_relbias |>
  ggplot(
    aes(
      x = 100 * (estimate - true) / true,
      y = method,
      fill = method,
      col = method#,
      #height = stat(density)
    )
  )  +
  geom_vline(
    aes(xintercept = 0),
    linetype = "solid",
    linewidth = 1,
    col = cols[3],
    alpha = 1
  ) +
  stat_density_ridges(
    scale = 1,
    quantile_lines = TRUE,
    quantile_fun = function(value, ...) mean(value, na.rm = TRUE),
    alpha = 0.5,
    color = "black",
    linewidth = 0.75
  ) +
  #geom_density_ridges(stat = "density", scale = 1, col = NA,
  #                    alpha = 0.5, ) +
  stat_density_ridges(
    scale = 1,
    #quantile_lines = TRUE,
    #quantile_fun = function(value, ...) mean(value, na.rm = TRUE),
    alpha = 0.5,
    aes(col = method, fill = NA),
    linewidth = 0.75
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  facet_grid(
   failure_time_model * censoring_type ~ prob_space,
   labeller = all_labels
  ) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = cols[c(1, 2, 6, 4, 5)]) +
  scale_colour_manual(values = cols[c(1, 2, 6, 4, 5)]) +
  scale_x_continuous(minor_breaks = NULL,
                     breaks = c(0, 10, 25, 50, 75, -10, -25, -50, -75)) +
  labs(y = "Method", x = "100 * (Estimate - True) / True (%)")

#colorspace::lighten()

ggsave(
  filename = "analysis/figures/bias_X_ridges.pdf",
  width = 7,
  scale = 1.25,
  height = 10,
  device = cairo_pdf
)



test <- copy(coefs_main[term == "X" & !(method %in% c( "Full data"))])
test[, rel_bias := 100 * (estimate - true) / true]

test[, .(
  dens = {
    dx <- density(rel_bias)
    approx(dx$x,dx$y,xout=rel_bias)
  }
), by = c(
  "method",
  "censoring_type",
  "failure_time_model",
  "prob_space"
)]


coefs_main[term == "X" & !(method %in% c( "Full data"))] |>
  ggplot(aes(x = method, y = estimate)) +
  geom_jitter(
    aes(alpha = abs(estimate - true), col = method),
    size = 2.5, width = 0.25,
    #alpha = 0.1,
    shape = 16
  ) +
  scale_alpha_continuous(range = c(0.01, 0.5))
  facet_grid(
    failure_time_model * censoring_type ~ prob_space,
    labeller = all_labels
  ) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  scale_y_continuous(minor_breaks = NULL, breaks = c(0, 10, 25, 50, 75, -10, -25, -50, -75)) +
  scale_x_discrete(limits = rev) +
  geom_hline(
    aes(yintercept = 0),
    linetype = "solid",
    linewidth = 1,
    col = cols[3],
    alpha = 1
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    fatten = 1.5,
    #linewidth = 0.5,
    col = "black",
    alpha = 0.75,
    lineend = "round"
  ) +
  scale_color_manual(values = cols[c(1, 2, 6, 4, 5)]) +
  labs(x = "Method", y = "100 * (Estimate - True) / True (%)") +
