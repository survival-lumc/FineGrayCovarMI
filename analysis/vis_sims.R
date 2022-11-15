
# Need to add 'true' estimates using dataset pre-missings
data <- #readRDS("data-raw/sims_nocens.rds") |>
  readRDS("data-raw/sims_cens.rds") |>
  #mutate(term = fct_recode(term, "X" = "X2")) |>
  mutate(true = case_when(term == "X" ~ 0.75, term == "Z" ~ 0.5))

data |>
  ggplot(aes(method, estimate, fill = method)) +
  geom_boxplot() +
  theme_minimal(base_size = 22) +
  geom_hline(aes(yintercept = true), linetype = "dashed", size = 1) +
  facet_wrap(~  term) +
  scale_fill_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 45), legend.position = "none") +
  labs(x = "Method", y = "Estimate")


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
