# Figures lunchmeeting talk

# Libraries

#remotes::install_github("vankesteren/firatheme")
#https://github.com/vankesteren/firatheme/issues/6
library(firatheme)
library(data.table)
library(prodlim)
library(riskRegression)
library(survival)
library(mice)
library(ggplot2)
library(smcfcs)


# Overall
#n <- 2000
predictor_formulas <- list("cause1" = ~ X + Z, "cause2" = ~ X + Z)
#X_type = "binary"
source("R/data-generation.R")
source("R/compute-true-cumincs.R")

t_seq <- seq(0, 12, by = 0.1)

# Define some reference patients
newdat_base <- data.frame("X" = 0, "Z" = 0)
newdat_X <- data.frame("X" = 1, "Z" = 0)
newdat_interm <- data.frame("X" = 1, "Z" = 1.5)
newdat_high <- data.frame("X" = 1, "Z" = 3)
newdat_low <- data.frame("X" = 0, "Z" = -3)

# Ggplot
theme_set(theme_fira())


# Indirect ----------------------------------------------------------------


params_indirect <- list(
  "cause1" = list(
    "betas" = c(0.5, 0.25),
    "base_rate" = 0.15,
    "base_cuminc" = 0.2
  ),
  "cause2" = list(
    "betas" = c(0.5, 0.25),
    "base_rate" = 0.5, #1.2, #0.1,
    "base_shape" = 0.5
  )
)

df_base_indirect <- compute_true_cuminc(
  t = t_seq,
  newdat = newdat_X,
  params = params_indirect,
  model_type = "indirect",
  predictor_formulas = predictor_formulas
)

df_base_indirect_haz <- cbind.data.frame(
  "time" = t_seq,
  "sd_cause1" = haz_sd_cause1(t_seq, as.matrix(newdat_base), params_indirect$cause1),
  "cs_cause2" = haz_cs_cause2(t_seq, as.matrix(newdat_base), params_indirect$cause2),
  "cs_cause1" = haz_cs_cause1(t_seq, as.matrix(newdat_base), params_indirect$cause1, params_indirect$cause2)
)

df_base_indirect_hazX <- cbind.data.frame(
  "time" = t_seq,
  "sd_cause1" = haz_sd_cause1(t_seq, as.matrix(newdat_X), params_indirect$cause1),
  "cs_cause2" = haz_cs_cause2(t_seq, as.matrix(newdat_X), params_indirect$cause2),
  "cs_cause1" = haz_cs_cause1(t_seq, as.matrix(newdat_X), params_indirect$cause1, params_indirect$cause2)
)

setDT(df_base_indirect)
setDT(df_base_indirect_haz)
setDT(df_base_indirect_hazX)

# First plot baseline hazards
melt.data.table(
  data = df_base_indirect_haz,
  id.vars = "time",
  variable.name = "hazard_type",
  value.name = "hazard"
) |>
  ggplot(aes(time, hazard)) +
  geom_line(aes(linetype = hazard_type, col = hazard_type), size = 1.25) +
  scale_color_fira() +
  coord_cartesian(expand = 0, xlim = c(0, 12)) +
  labs(x = "Time", y = "Baseline hazard") +
  annotate(
    "text", x = 3, y = 0.2, label = "CSH 2",
    family =  theme_get()$text[["family"]],
    col = firaPalette()[2],
    size = 4
  ) +
  annotate(
    "text", x = 0.1, y = 0.025, label = "SDH 1",
    family =  theme_get()$text[["family"]],
    col = firaPalette()[1],
    size = 4,
    hjust = 0
  ) +
  annotate(
    "text", x = 2.5, y = 0.1, label = "CSH 1",
    family =  theme_get()$text[["family"]],
    col = firaPalette()[5],
    size = 4,
    hjust = 0
  ) +
  theme(
    legend.position = "none"
  )

firaSave("basehaz_indirect.png", device = "png", pointsize = 14, width = 8, height = 6)


# Make hazard ratios
p_cause1_indirect <- cbind.data.frame(
  "time" = t_seq,
  "haz_base" =  log(df_base_indirect_haz$cs_cause1),
  "haz_X" = log(df_base_indirect_hazX$cs_cause1),
  "HR" = log(df_base_indirect_hazX$cs_cause1 / df_base_indirect_haz$cs_cause1)
) |>
  setDT() |>
  melt.data.table(
    id.vars = "time",
    variable.name = "hazard_type",
    value.name = "hazard"
  ) |>
  ggplot(aes(time, hazard)) +
  geom_line(aes(col = hazard_type, linetype = hazard_type), size = 1.25) +
  scale_linetype_manual(values = c("dashed", "dashed", "solid")) +
  scale_color_fira() +
  coord_cartesian(expand = 0, xlim = c(0, 12), ylim = c(-7.5, 3)) +
  labs(x = "Time", y = "Log Hazard/log HR") +
  theme(
    legend.position = "none"
  ) +
  annotate(
    "text", x = 0.5, y = -0.25, label = "CSH1 | X = 1, Z = 0",
    family =  theme_get()$text[["family"]],
    col = firaPalette()[2],
    size = 4,
    hjust = 0
  ) +
  annotate(
    "text", x = 0.5, y = -4.5, label = "CSH1 | X = 0, Z = 0",
    family =  theme_get()$text[["family"]],
    col = firaPalette()[1],
    size = 4,
    hjust = 0
  ) +
  annotate(
    "text", x = 0.5, y = 1.75, label = "log HR",
    family =  theme_get()$text[["family"]],
    col = firaPalette()[5],
    size = 4,
    hjust = 0
  ) +
  ggtitle("Cause 1")


# Cause 2
p_cause2_indirect <- cbind.data.frame(
  "time" = t_seq,
  "haz_base" =  log(df_base_indirect_haz$cs_cause2),
  "haz_X" = log(df_base_indirect_hazX$cs_cause2),
  "HR" = log(df_base_indirect_hazX$cs_cause2 / df_base_indirect_haz$cs_cause2)
) |>
  setDT() |>
  melt.data.table(
    id.vars = "time",
    variable.name = "hazard_type",
    value.name = "hazard"
  ) |>
  ggplot(aes(time, hazard)) +
  geom_line(aes(col = hazard_type, linetype = hazard_type), size = 1.25) +
  scale_linetype_manual(values = c("dashed", "dashed", "solid")) +
  scale_color_fira() +
  coord_cartesian(expand = 0, xlim = c(0, 12), ylim = c(-7.5, 3)) +
  labs(x = "Time", y = "Log Hazard/log HR") +
  theme(
    legend.position = "none"
  ) +
  annotate(
    "text", x = 6.5, y = -1.25, label = "CSH2 | X = 1, Z = 0",
    family =  theme_get()$text[["family"]],
    col = firaPalette()[2],
    size = 4,
    hjust = 0
  ) +
  annotate(
    "text", x = 6.5, y = -3, label = "CSH2 | X = 0, Z = 0",
    family =  theme_get()$text[["family"]],
    col = firaPalette()[1],
    size = 4,
    hjust = 0
  ) +
  annotate(
    "text", x = 6.5, y = 1, label = "log HR",
    family =  theme_get()$text[["family"]],
    col = firaPalette()[5],
    size = 4,
    hjust = 0
  ) +
  ggtitle("Cause 2")


cowplot::plot_grid(p_cause1_indirect, p_cause2_indirect)
firaSave("HRs_indirect.png", device = "png", pointsize = 14, width = 8, height = 6)







xmelt.data.table(
  data = df_base_indirect,
  id.vars = "time",
  variable.name = "cause",
  value.name = "cuminc"
) |>
  ggplot(aes(time, cuminc)) +
  geom_area(aes(fill = cause)) +
  scale_fill_fira() +
  coord_cartesian(expand = 0, xlim = c(0, 12), ylim = c(0, 1.1))

melt.data.table(
  data = df_base_indirect_haz,
  id.vars = "time",
  variable.name = "hazard_type",
  value.name = "hazard"
) |>
  ggplot(aes(time, hazard)) +
  geom_line(aes(linetype = hazard_type, col = hazard_type), size = 1.25) +
  scale_color_fira() +
  coord_cartesian(expand = 0, xlim = c(0, 12))





# Emprirically new data
df_indirect <- generate_complete_dataset(
  n = 10000,
  params = params_indirect,
  model_type = "indirect",
  predictor_formulas,
  X_type = "binary"
)

mod <- FGR(Hist(time, D) ~ X + Z, data = df_indirect, cause = 2)

uts <- mod$crrFit$uftime
cumhaz_two <- mod$crrFit$bfitj
#plot(uts, diff(c(0, cumhaz_two)) / diff(c(0, uts)))
plot(uts, cumhaz_two, type = "s")
mod2 <- FGR(Hist(time, D) ~ X + Z, data = df_indirect, cause = 1)
lines(mod2$crrFit$uftime, mod2$crrFit$bfitj, type = "s")
table(df_indirect$D)


plot(uts[-1], diff(cumhaz_two) / diff(uts), type = "l")


#purrr:: map dfr, put plots next to

# Direct ------------------------------------------------------------------


df_base_direct <- compute_true_cuminc(
  t = t_seq,
  newdat = base,
  params = params_fgs,
  model_type = "both_fg",
  predictor_formulas = predictor_formulas
)


# Two FGs -----------------------------------------------------------------







# both_fgs
params_fgs <- list(
  "cause1" = list(
    "betas" = c(0.5, 0.25),
    "p" = 0.15,
    "base_rate" = 1,
    "base_shape" = 0.6
  ),
  "cause2" = list(
    "betas" = c(0.5, 0.25),
    "p" = 0.25,
    "base_rate" = 1,
    "base_shape" = 0.6
  )
)

df_true_both <- compute_true_cuminc(
  t = t_seq,
  newdat = newdat_base,
  params = params_fgs,
  model_type = "both_fg",
  predictor_formulas = predictor_formulas
)
setDT(df_true_both)

melt.data.table(
  data = df_true_both,
  id.vars = "time",
  variable.name = "cause",
  value.name = "cuminc"
) |>
  ggplot(aes(time, cuminc, fill = cause)) +
  geom_area() +
  coord_cartesian(expand = 0, ylim = c(0, 1.1)) +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dashed", size = 1) +
  scale_fill_viridis_d()

df_true_both |>
  ggplot(aes(time, cause1)) +
  geom_line() +
  geom_line(aes(y = cause1 + cause2)) +
  geom_hline(yintercept = params_fgs$cause1$p, linetype = "dotted")



# Here
df_true_both |>
  ggplot(aes(time, cause1)) +
  geom_line() +
  geom_line(aes(y = params_fgs$cause1$p + cause2)) +
  geom_hline(yintercept = c(
    params_fgs$cause1$p,
    params_fgs$cause1$p + params_fgs$cause2$p,
    1
  ), linetype = "dashed") +
  geom_ribbon(aes(ymin = 0, ymax = cause1),
              alpha = 0.5, fill = firaPalette()[1]) +
  geom_ribbon(aes(ymin = params_fgs$cause1$p,
                  ymax = cause2 + params_fgs$cause1$p),
              alpha = 0.5, fill = firaPalette()[2]) +
  geom_ribbon(aes(ymin = params_fgs$cause1$p + params_fgs$cause2$p,
                  ymax = 1),
              alpha = 0.5, fill = "lightgray") +
  coord_cartesian(expand = 0, xlim = c(0, 12), ylim = c(0, 1.1)) +
  labs(y = "Cumulative incidence", x = "Time") +
  annotate(
    "text", x = 9.25, y = 0.1, label = "Cause 1",
    family =  theme_get()$text[["family"]],
    #4col = firaPalette()[5],
    size = 4,
    hjust = 0
  ) +
  annotate(
    "text", x = 9.25, y = 0.3, label = "Cause 2",
    family =  theme_get()$text[["family"]],
    #4col = firaPalette()[5],
    size = 4,
    hjust = 0
  ) +
  annotate(
    "text", x = 9.25, y = 0.75, label = "Cured fraction",
    family =  theme_get()$text[["family"]],
    #4col = firaPalette()[5],
    size = 4,
    hjust = 0
  )

firaSave("basecuminc_bothFGs.png", device = "png", pointsize = 16,
         width = 8, height = 6)




# Try direct both ---------------------------------------------------------



params_direct <- list(
  "cause1" = list(
    "betas" = c(0.5, 0.25),
    "p" = 0.2,
    "base_rate" = 0.5,
    "base_shape" = 1.2
  ),
  "cause2" = list(
    "betas" = c(0.5, 0.25),
    "base_rate" = 0.15,
    "base_shape" = 0.3
  )
)


df_true_direct <- compute_true_cuminc(
  t = t_seq,
  newdat = newdat_base,
  params = params_direct,
  model_type = "direct",
  predictor_formulas = predictor_formulas
)

df_true_direct |>
  ggplot(aes(time, cause1)) +
  geom_line() +
  geom_line(aes(y = params_direct$cause1$p + cause2)) +
  geom_hline(yintercept = c(
    params_direct$cause1$p,
    params_direct$cause1$p + params_direct$cause2$p,
    1
  ), linetype = "dashed") +
  geom_ribbon(aes(ymin = 0, ymax = cause1),
              alpha = 0.5, fill = firaPalette()[1]) +
  geom_ribbon(aes(ymin = params_direct$cause1$p,
                  ymax = cause2 + params_direct$cause1$p),
              alpha = 0.5, fill = firaPalette()[2]) +
  coord_cartesian(expand = 0, xlim = c(0, 12), ylim = c(0, 1.1)) +
  labs(y = "Cumulative incidence", x = "Time") +
  annotate(
    "text", x = 9.25, y = 0.1, label = "Cause 1",
    family =  theme_get()$text[["family"]],
    #4col = firaPalette()[5],
    size = 4,
    hjust = 0
  ) +
  annotate(
    "text", x = 9.25, y = 0.3, label = "Cause 2",
    family =  theme_get()$text[["family"]],
    #4col = firaPalette()[5],
    size = 4,
    hjust = 0
  )

firaSave("basecuminc_direct.png", device = "png", pointsize = 16,
         width = 8, height = 6)


plot(df_true_direct$time, df_true_direct$cause1, ylim = c(0, 1), type = "s",
     xlim = c(0, 12000))
lines(df_true_direct$time, df_true_direct$cause2 + params_direct$cause1$p,
      ylim = c(0, 1), type = "s")
