library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggsci)
library(ggpubr)

setwd("~/Documents/Github/lds-mc-julia/")
ce_df <- read.csv("output/ce_df.csv")



plot11 <- ce_df |>
  select(-cvar_mstd) |>
  pivot_longer(c('ev', 'cvar', 'mstd'), names_to = 'name', values_to = 'value') |>
  mutate(name = factor(name, c('ev', 'cvar', 'mstd'), c("EV", "M-CVaR", "M-SD"))) |>
  pivot_wider(names_from = var, values_from = value) |>
  ggplot() +
  geom_hline(yintercept = 0, color = 'gray50') +
  geom_ribbon(aes(ymin = lb, ymax = ub, x = alpha, fill = name), alpha = 0.3) +
  geom_line(aes(x = alpha, y = median, color = name)) +
  scale_y_continuous("Change from baseline", labels = scales::percent) +
  scale_x_continuous("α") +
  ggsci::scale_color_nejm() +
  ggsci::scale_fill_nejm() +
  ggpubr::theme_pubr() +
  theme(legend.title = element_blank())

plot12 <- ce_df |>
  select(c(-ev, -cvar, -mstd)) |>
  rename(value = cvar_mstd) |>
  pivot_wider(names_from = var, values_from = value) |>
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(ymin = lb, ymax = ub, x = alpha), alpha = 0.3) +
  geom_line(aes(x = alpha, y = median)) +
  geom_hline(yintercept = 0, color = 'gray50') +
  scale_y_continuous("Change from M-SD", labels = scales::percent) +
  scale_x_continuous("α") +
  ggpubr::theme_pubr()

plot11 + plot12 + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
