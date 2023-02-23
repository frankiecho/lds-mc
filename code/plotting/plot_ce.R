library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggsci)
library(ggpubr)

setwd("~/Documents/Github/lds-mc-julia/")

## Plot heatmaps
sp_weight <- read.csv("output/spatial_weights.csv")

p11 <- sp_weight %>%
  #pivot_longer(c(x3, x4, x5)) %>%
  ggplot(aes(x = x, y = y, fill = x3)) +
  geom_tile(color = 'gray50') +
  theme_void() +
  scale_fill_viridis_c(option = 'mako', direction = -1) +
  theme(legend.position = 'bottom', legend.title = element_blank())

p12 <- sp_weight %>%
  #pivot_longer(c(x3, x4, x5)) %>%
  ggplot(aes(x = x, y = y, fill = x4)) +
  geom_tile(color = 'gray50') +
  theme_void() +
  scale_fill_viridis_c(option = 'magma', direction = -1) +
  theme(legend.position = 'bottom', legend.title = element_blank())

p13 <- sp_weight %>%
  #pivot_longer(c(x3, x4, x5)) %>%
  ggplot(aes(x = x, y = y, fill = x5)) +
  geom_tile(color = 'gray50') +
  theme_void() +
  scale_fill_viridis_c(option = 'magma', direction = -1) +
  theme(legend.position = 'bottom', legend.title = element_blank())

ggsave("plots/spatial_weight_heatmap.png", p11 + p12 + p13 + plot_annotation(tag_levels = 'a'))

## Plot change in CE
ce_df <- read.csv("output/ce_df_try.csv")

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
