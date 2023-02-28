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

heatmap <- p11 + p12 + p13 + plot_annotation(tag_levels = 'a')
heatmap
ggsave("plots/spatial_weight_heatmap.png", heatmap)

## Plot change in CE
ce_df <- read.csv("output/ce_df_baseline_500.csv") %>%
  filter(alpha <= 30)

plot11 <- ce_df |>
  dplyr::select(-cvar_mstd) |>
  pivot_longer(c('ev', 'cvar', 'mstd'), names_to = 'name', values_to = 'value') |>
  mutate(name = factor(name, c('ev', 'cvar', 'mstd'), c("EV", "M-CVaR", "M-SD"))) |>
  pivot_wider(names_from = var, values_from = value) |>
  ggplot() +
  geom_hline(yintercept = 0, color = 'gray50') +
  #geom_ribbon(aes(ymin = min, ymax = max, x = alpha, fill = name), alpha = 0.15) +
  geom_ribbon(aes(ymin = lb, ymax = ub, x = alpha, fill = name), alpha = 0.2) +
  geom_line(aes(x = alpha, y = median, color = name)) +
  scale_y_continuous("Change from baseline", labels = scales::percent) +
  scale_x_continuous("θ") +
  ggsci::scale_color_nejm() +
  ggsci::scale_fill_nejm() +
  ggpubr::theme_pubr() +
  theme(legend.title = element_blank())

plot12 <- ce_df |>
  dplyr::select(c(-ev, -cvar, -mstd)) |>
  rename(value = cvar_mstd) |>
  pivot_wider(names_from = var, values_from = value) |>
  ggplot() +
  geom_hline(yintercept = 0) +
  #geom_ribbon(aes(ymin = min, ymax = max, x = alpha), alpha = 0.15) +
  geom_ribbon(aes(ymin = lb, ymax = ub, x = alpha), alpha = 0.2) +
  geom_line(aes(x = alpha, y = median)) +
  geom_hline(yintercept = 0, color = 'gray50') +
  scale_y_continuous("Change from M-SD", labels = scales::percent) +
  scale_x_continuous("θ") +
  ggpubr::theme_pubr()

(delta_plot <- plot11 + plot12 + plot_layout(guides = "collect") & theme(legend.position = 'bottom') &  plot_annotation(tag_levels = 'a') )
ggsave("plots/delta_plot.png", delta_plot, units = 'cm', width = 20, height = 12)
