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