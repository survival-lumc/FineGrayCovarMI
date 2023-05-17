source("R/compute-true-cumincs.R")
library(tidyverse)
library(patchwork)
library(extrafont)
library(Manu)

theme_set(
  theme_light(base_size = 16, base_family = "Roboto Condensed") +
    theme(
      strip.background = element_rect(fill = Manu::get_pal("Hoiho")[[2]], colour = "white"),
      strip.text = element_text(colour = 'white')
    )
)

# Timepoints
times <- c(seq(0.001, 2.5, length.out = 40),seq(2.6, 10, length.out = 10))
combs <- tidyr::crossing(
  mechanism = c("correct_FG", "misspec_FG"),
  p = c(0.15, 0.65)
)

pspace_labels <- c("0.15" = "p = 0.15", "0.65" = "p = 0.65")
failure_model_labels <- c("correct_FG" = "Well specified FG", "misspec_FG" = "Misspecified FG")
all_labels <- labeller(
  p = pspace_labels,
  mech = failure_model_labels
)

hazs <- map2(
  .x = combs$mechanism,
  .y = combs$p,
  .f = ~ {
    params_slug <- switch(
      .x,
      "misspec_FG" = paste0("params_weibull_lfps_", .y),
      "correct_FG" = paste0("true_params_correct_FG_", .y)
    )
    params <- do.call(tar_read, list(as.symbol(params_slug)))
    compute_true(
      t = times,
      newdat = cbind.data.frame(X = 0, Z = 0),
      params = params,
      model_type = .x
    ) |> cbind(p = .y, mech = .x)
  }
)

HRs <- map2(
  .x = combs$mechanism,
  .y = combs$p,
  .f = ~ {
    params_slug <- switch(
      .x,
      "misspec_FG" = paste0("params_weibull_lfps_", .y),
      "correct_FG" = paste0("true_params_correct_FG_", .y)
    )
    params <- do.call(tar_read, list(as.symbol(params_slug)))

    base <- compute_true(
      t = times,
      newdat = cbind.data.frame(X = 0, Z = 0),
      params = params,
      model_type = .x
    )

    x1 <- compute_true(
      t = times,
      newdat = cbind.data.frame(X = 1, Z = 0),
      params = params,
      model_type = .x
    )

    df_hrs <- data.table(
      time = base$time,
      HR_subdist = x1$haz_subdist1 / base$haz_subdist1,
      HR_cs1 = x1$haz_cs1 / base$haz_cs1,
      HR_cs2 = x1$haz_cs2 / base$haz_cs2
    ) |> cbind(p = .y, mech = .x)
  }
)

  #coord_cartesian(ylim = c(0, 2))



# Try all in one ----------------------------------------------------------



p1 <- rbindlist(hazs) |>
  pivot_longer(
    cols = starts_with("haz"),
    names_to = "hazard_type",
    values_to = "hazard"
  ) |>
  ggplot(aes(time, hazard)) +
  geom_line(
    aes(group = hazard_type, col = hazard_type, linetype = hazard_type),
    size = 1.5, alpha = 0.8
  ) +
  facet_wrap(p * mech ~ ., ncol = 4, labeller = all_labels) +
  #theme_minimal() +
  scale_color_manual(
    values = Manu::get_pal("Hoiho")[1:3],
    labels = c("Cause-spec cause 1", "Cause-spec cause 2", "Subdist. cause 1")
  ) +
  scale_linetype_manual(
    values = 1:3,
    labels = c("Cause-spec cause 1", "Cause-spec cause 2", "Subdist. cause 1")
  ) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(y = "Baseline hazards", x = "Time",
       linetype = NULL,
       col = NULL)

p2 <- rbindlist(hazs) |>
  pivot_longer(
    cols = starts_with("cuminc"),
    names_to = "event",
    values_to = "cuminc"
  ) |>
  ggplot(aes(time, cuminc)) +
  geom_line(
    aes(group = event, col = event, linetype = event),
    size = 1.5, alpha = 0.8
  ) +
  facet_wrap(p * mech ~ ., ncol = 4) +
  #theme_minimal() +
  scale_color_manual(
    values = Manu::get_pal("Hoiho")[1:2],
    labels = c("Cause 1", "Cause 2")
  ) +
  scale_linetype_manual(
    values = 1:2,
    labels = c("Cause 1", "Cause 2")
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(y = "Baseline cumulative incidence", x = "Time",
       linetype = NULL,
       col = NULL)


p3 <- rbindlist(HRs) |>
  pivot_longer(
    cols = starts_with("HR"),
    names_to = "type",
    values_to = "HR"
  ) |>
  ggplot(aes(time, HR)) +
  geom_line(
    aes(group = type, col = type, linetype = type),
    size = 1.5, alpha = 0.8
  ) +
  facet_wrap(p * mech ~ ., ncol = 4) +
  #theme_minimal() +
  scale_color_manual(
    values = Manu::get_pal("Hoiho")[1:3],
    labels = c("Cause-spec cause 1", "Cause-spec cause 2", "Subdist. cause 1")
  ) +
  scale_linetype_manual(
    values = 1:3,
    labels = c("Cause-spec cause 1", "Cause-spec cause 2", "Subdist. cause 1")
  ) +
  geom_hline(yintercept = 1, col = "black", linetype = "dotted") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  labs(
    y = "Hazard ratio X = 1 vs X = 0",
    x = "Time",
    linetype = NULL,
    col = NULL
  )


p3


ggsave(
  plot = p1 / p2 / p3,
  "scenarios_vis.png",
  dpi = 300,
  width = 13,
  height = 9
)



library("geomtextpath")

test[[2]] |>
  pivot_longer(
    cols = starts_with("HR"),
    names_to = "hazard_type",
    values_to = "hazard"
  ) |>
  ggplot(aes(time, hazard, col = hazard_type, linetype = hazard_type)) +
  geom_textline(aes(label = hazard_type), hjust = 0.75) +
  #geom_line(linewidth = 2) +
  labs(y = "HR") +
  coord_cartesian(ylim = c(0, 7.5)) +
  theme(legend.position = "none")


