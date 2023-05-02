library(tidyverse)
library(extrafont)
library(Manu)

# Need to add 'true' estimates using dataset pre-missings
data <- #readRDS("data-raw/sims_nocens.rds") |>
  readRDS("data-raw/sims_cens.rds") |>
  #mutate(term = fct_recode(term, "X" = "X2")) |>
  mutate(true = case_when(term == "X" ~ 0.75, term == "Z" ~ 0.5))

theme_set(
  theme_light(base_size = 16, base_family = "Roboto Condensed") +
    theme(
      strip.background = element_rect(fill = Manu::get_pal("Hoiho")[[2]], colour = "transparent"),
      strip.text = element_text(colour = 'white')
    )
)

# USE JITTERS LATER!!
data |>
  ggplot(aes(method, estimate, fill = method)) +
  #geom_jitter() +
  geom_boxplot() +
  geom_hline(aes(yintercept = true), linetype = "dashed", size = 1) +
  facet_wrap(~  term) +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(
    labels = c(
      "Compl. cases",
      "Full data",
      "MICE cause-spec",
      "MICE subdist",
      "SMC-FCS cause-spec",
      "SMC-FCS + kmi"
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = "Method", y = "Estimate")

ggsave("interim-sims.png", dpi = 300, units = "in", width = 9, height = 6)


library(rsimsum)
sims_summ <- multisimsum(
  data = data,
  par = "term",
  estvarname = "estimate",
  se = "std.error",
  true = "true",
  methodvar = "method",
  ref = "CCA"
)

sim_summ_X <- simsum(
  data = data |> filter(term == "X"),
  estvarname = "estimate",
  se = "std.error",
  true = "true",
  methodvar = "method",
  ref = "CCA"
)

summary(sim_summ_X)


# Want plots of subdist and cause-spec hazards..
# for latter; differentiate cuminc, and diving by event free surv?
