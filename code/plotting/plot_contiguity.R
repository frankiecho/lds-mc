## Plot change in contiguity
contiguity <- read_csv('output/contiguity_baseline_100.csv') |>
  filter(alpha <= 30)
contiguity |>
  pivot_longer(c('mstd','cvar')) |>
  mutate(name = factor(name, c('mstd', 'cvar', 'ev'), c("M-SD", "M-CVaR", "EV"))) |>
  pivot_wider(names_from = var, values_from = value) |>
  ggplot() +
  #geom_hline(yintercept = 0, color = 'gray50') +
  #geom_ribbon(aes(ymin = min, ymax = max, x = alpha, fill = name), alpha = 0.151) +
  geom_hline(yintercept = 0, color = 'gray50') +
  geom_ribbon(aes(ymin = lb, ymax = ub, x = alpha, fill = name), alpha = 0.2) +
  geom_line(aes(x = alpha, y = median, color = name)) +
  scale_x_continuous("θ") +
  scale_y_continuous("Connectivity index", labels = scales::percent) +
  ggsci::scale_color_nejm() +
  ggsci::scale_fill_nejm() +
  ggpubr::theme_pubr() +
  coord_cartesian(expand = F) +
  theme(legend.title = element_blank())

ggsave("plots/contiguity_plot.png", units = 'cm', width = 12, height = 12)