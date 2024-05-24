# Use example simulations data from https://ellessenne.github.io/rsimsum/articles/B-relhaz.html
library(rsimsum)
library(ggplot2)
library(data.table)

dat_sims <- data.table(rsimsum::relhaz)
head(dat_sims)

# Check number of replications per scenario
dat_sims[, .(n_replications = .N), by = c(
  "model", # method variable
  "baseline", # Scenario dimension n.1
  "n" # Scenario dimension n.2
)]

# Some pre-processing
dat_sims[, true := -0.5]
dat_sims[, relbias := 100 * (theta - true) / true]
head(dat_sims)

# To get Monte-Carlo confidence intervals(95%)
crit <- qnorm((1 - 0.95) / 2, lower.tail = FALSE)

# Set-up facet labels
facet_labs <- labeller(
  n = c("50" = "n = 50", "250" = "n = 250"),
  baseline = c("Exponential" = "Exponential baseline", "Weibull" = "Weibull baseline")
)

# Set-up global ggplot theme
theme_set(
  theme_light(base_size = 16) + #, base_family = "Roboto Condensed") # if custom font
    theme(
      strip.background = element_rect(fill = cols[2], colour = "white"),
      strip.text = element_text(colour = 'white')
    )
)

# Colours from https://github.com/G-Thomson/Manu
cols <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#BC8E7D", "#2F3D70")


# THE PLOT
dat_sims |>
  ggplot(aes(model, relbias)) +
  geom_jitter(
    aes(col = model),
    size = 2.5,
    width = 0.25,
    shape = 16,
    alpha = 0.5
  ) +
  facet_grid(baseline ~ n, labeller = facet_labs) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  ) +
  scale_y_continuous(
    minor_breaks = NULL,
    breaks = c(0, 25, 50, 75, 100, -25, -50, -75, -100)
  ) +
  scale_x_discrete(limits = rev) +
  geom_hline(
    aes(yintercept = 0),
    linetype = "solid",
    linewidth = 1,
    col = "black",
    alpha = 1
  ) +
  # Estimated relative bias + Monte-Carlo CI
  # If you comment out fun.min and fun.max, you just get the average (more aesthetic, but less info)
  stat_summary(
    fun = function(x) mean(x),
    fun.min = function(x) mean(x) - crit * (sd(x) / sqrt(length(x))),
    fun.max = function(x) mean(x) + crit * (sd(x) / sqrt(length(x))),
    geom = "crossbar",
    alpha = 0.25,
    fill = "black",
    fatten = 1,
    linewidth = 0.25
  ) +
  scale_color_manual(values = cols[c(1, 2, 4)]) +
  labs(x = "Method", y = "100 * (Estimate - True) / True (%)")
