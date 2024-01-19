library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggsci)
library(ggpubr)

setwd("~/Documents/Github/lds-mc/")

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

# Decisions
p14 <- sp_weight %>%
  arrange(-x4) %>%
  mutate(rank = 1:nrow(sp_weight)) %>%
  mutate(x6 = factor(rank < 100, labels=c('No', 'Yes'))) %>%
  ggplot(aes(x = x, y = y, fill = x6)) +
  geom_tile(color = 'gray50') +
  scale_fill_manual(values = c('#dddddd', '#009E73')) +
  theme_void() +
  theme(legend.position = 'bottom', legend.title = element_blank())

heatmap <- p11 + p12 + p13 + p14 + plot_annotation(tag_levels = 'a') & plot_layout(nrow = 1)
heatmap
ggsave("plots/spatial_weight_heatmap.png", heatmap, width = 1500, height = 500, units = 'px', scale = 2)

## Plot change in CE
nsims = 5
max_alpha = 30
for (i in 1:25) {
ce_df <- read_csv(paste0("output/ce_df_param_search_", i, "_", nsims, ".csv"), show_col_types = F) %>%
  filter(alpha <= max_alpha)

plot11 <- ce_df |>
  dplyr::select(-cvar_mstd, -cvar, -mstd, -ev) |>
  dplyr::rename(cvar = cvar_ev, mstd = mstd_ev) |>
  pivot_longer(c('cvar', 'mstd'), names_to = 'name', values_to = 'value') |>
  mutate(name = factor(name, c('mstd', 'cvar', 'ev'), c("M-SD", "M-CVaR", "EV"))) |>
  pivot_wider(names_from = var, values_from = value) |>
  ggplot() +
  geom_hline(yintercept = 0, color = 'gray50') +
  #geom_ribbon(aes(ymin = min, ymax = max, x = alpha, fill = name), alpha = 0.15) +
  geom_ribbon(aes(ymin = lb, ymax = ub, x = alpha, fill = name), alpha = 0.2) +
  geom_line(aes(x = alpha, y = median, color = name)) +
  scale_y_continuous("CE: change from EV solution", labels = scales::percent) +
  scale_x_continuous("θ") +
  ggsci::scale_color_nejm() +
  ggsci::scale_fill_nejm() +
  ggpubr::theme_pubr() +
  #coord_cartesian(expand = F) +
  theme(legend.title = element_blank())

plot12 <- ce_df |>
  dplyr::select(c(-ev, -cvar, -mstd, -cvar_ev, -mstd_ev)) |>
  rename(value = cvar_mstd) |>
  pivot_wider(names_from = var, values_from = value) |>
  ggplot() +
  geom_hline(yintercept = 0) +
  #geom_ribbon(aes(ymin = min, ymax = max, x = alpha), alpha = 0.15) +
  geom_ribbon(aes(ymin = lb, ymax = ub, x = alpha), alpha = 0.2) +
  geom_line(aes(x = alpha, y = median)) +
  geom_hline(yintercept = 0, color = 'gray50') +
  scale_y_continuous("CE: change from M-SD", labels = scales::percent) +
  scale_x_continuous("θ") +
  #coord_cartesian(expand = F) +
  ggpubr::theme_pubr()

(delta_plot <- plot11 + plot12 + plot_layout(guides = "collect") & theme(legend.position = 'bottom') &  plot_annotation(tag_levels = 'a') )
ggsave(paste0("plots/delta_plot_", i, ".png"), delta_plot, units = 'cm', width = 20, height = 12)

## Plot change in contiguity
contiguity <- read_csv(paste0('output/distance_param_search_',i,'_', nsims, '.csv'), show_col_types = F) |>
  filter(alpha <= max_alpha)
contiguity |>
  pivot_longer(c('mstd','cvar')) |>
  mutate(name = factor(name, c('mstd', 'cvar', 'ev'), c("M-SD", "M-CVaR", "EV"))) |>
  pivot_wider(names_from = var, values_from = value) |>
  ggplot() +
  #geom_hline(yintercept = 0, color = 'gray50') +
  #geom_ribbon(aes(ymin = min, ymax = max, x = alpha, fill = name), alpha = 0.15) +
  geom_hline(yintercept = 0, color = 'gray50') +
  geom_ribbon(aes(ymin = lb, ymax = ub, x = alpha, fill = name), alpha = 0.2) +
  geom_line(aes(x = alpha, y = median, color = name)) +
  scale_x_continuous("θ") +
  scale_y_continuous("Distance index (change from baseline)", labels = scales::percent) +
  ggsci::scale_color_nejm() +
  ggsci::scale_fill_nejm() +
  ggpubr::theme_pubr() +
  coord_cartesian(expand = F) +
  theme(legend.title = element_blank())

ggsave(paste0("plots/distance_plot_", i, ".png"), units = 'cm', width = 15, height = 12)

## Plot downsides
downside <- read_csv(paste0('output/downside_param_search_', i, '_',nsims,'.csv'), show_col_types = F)
upside <- read_csv(paste0('output/downside_param_search_', i, '_',nsims,'.csv'), show_col_types = F)
risk_side <- bind_rows(list(downside=downside, upside=upside), .id = "side")

risk_side |>
  mutate(side = factor(side, c('downside', 'upside'), c('Downside risk', 'Upside gain'))) |>
  pivot_longer(c('mstd','cvar')) |>
  mutate(name = factor(name, c('mstd', 'cvar', 'ev'), c("M-SD", "M-CVaR", "EV"))) |>
  pivot_wider(names_from = var, values_from = value) |>
  ggplot() +
  #geom_hline(yintercept = 0, color = 'gray50') +
  #geom_ribbon(aes(ymin = min, ymax = max, x = alpha, fill = name), alpha = 0.15) +
  geom_hline(yintercept = 0.05, color = 'gray50') +
  geom_ribbon(aes(ymin = lb, ymax = ub, x = alpha, fill = name), alpha = 0.2) +
  geom_line(aes(x = alpha, y = median, color = name)) +
  scale_x_continuous("θ") +
  scale_y_continuous("Probability", labels = scales::percent) +
  ggsci::scale_color_nejm() +
  ggsci::scale_fill_nejm() +
  ggpubr::theme_pubr() +
  facet_wrap(~side) +
  coord_cartesian(expand = F) +
  theme(legend.title = element_blank())
ggsave(paste0("plots/downside_risk_", i, ".png"), units = 'cm', width = 20, height = 12)
}

## Theta-lambda matching
cvar_ce <- read_csv("output/cvar_ce.csv")
lambda_vec <- seq(0, 1, 0.05)
theta_vec <- seq(0, 50, 0.1)
cvar_ce_norm <- purrr::map_dfr(cvar_ce, function(x) (x - max(x)) / (max(x) - min(x)))
colnames(cvar_ce_norm) <- theta_vec
lambda_max <- data.frame(
  value = apply(cvar_ce_norm, 2, function(x) max(x)),
  lambda = apply(cvar_ce_norm, 2, function(x) lambda_vec[which.max(x)]),
  theta = as.character(theta_vec)
)
cvar_ce_norm$lambda <- lambda_vec

plot_theta <- c(0, 2.5, 5, 7.5, 10)
lambda_max_filter <- filter(lambda_max, theta %in% plot_theta) %>%
  mutate(theta = factor(theta, plot_theta))
cvar_ce_norm %>%
  pivot_longer(as.character(seq(0, 50, 0.1)), names_to = "theta", values_to = "value") %>%
  filter(theta %in% plot_theta) %>%
  mutate(theta = factor(theta, plot_theta)) %>%
  ggplot(aes(x = lambda, y = value, color = theta)) +
  geom_hline(yintercept = 0, linetype = 2, color = '#aaaaaa') +
  geom_point(data = lambda_max_filter) +
  geom_line(linewidth = 0.8) +
  see::scale_color_okabeito() +
  scale_y_continuous("ΔCE") +
  scale_x_continuous("𝜆")+
  labs(color = "θ") +
  coord_cartesian(ylim = c(-1, 0.2)) +
  theme_pubr() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
ggsave("plots/theta_lambda_matching.png", units = 'cm', width = 15, height = 12)
