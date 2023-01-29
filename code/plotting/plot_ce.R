library(tidyverse)
library(ggplot2)

setwd("~/Documents/Github/lds-mc-julia/")
ce_df <- read.csv("output/ce_df.csv")

ce_df |>
  select(-cvar_mstd) |>
  pivot_longer(c('ev', 'cvar', 'mstd'), names_to = 'name', values_to = 'value') |>
  pivot_wider(names_from = var, values_from = value) |>
  ggplot() +
  geom_ribbon(aes(ymin = lb, ymax = ub, x = alpha, fill = name), alpha = 0.3) +
  geom_line(aes(x = alpha, y = median, color = name)) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()

ce_df |>
  select(c(-ev, -cvar, -mstd)) |>
  rename(value = cvar_mstd) |>
  pivot_wider(names_from = var, values_from = value) |>
  ggplot() +
  geom_ribbon(aes(ymin = lb, ymax = ub, x = alpha), alpha = 0.3) +
  geom_line(aes(x = alpha, y = median)) +
  theme_minimal()
