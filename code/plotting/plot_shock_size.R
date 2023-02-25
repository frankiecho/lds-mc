library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(evmix)

shock_df <-read_csv("output/shock_df.csv")
nss <- read_csv("output/nss.csv")

s <- max(nss$s)
nss_table <- nss %>%
  group_by(nss) %>%
  summarise(prop = n())
nss_table$prop <- nss_table$prop / s

x <- seq(0,6,1)
binom_pdf <- data.frame(
  x = x,
  y = dbinom(x, 1600, 1/1600)
)

pal <- ggsci::pal_nejm()(8)

p1 <- ggplot(nss_table) +
  geom_bar(aes(x = nss, y = prop), stat = "identity", fill = 'gray60') +
  geom_line(data = binom_pdf, aes(x = x, y = y), color = pal[1]) +
  geom_point(data = binom_pdf, aes(x = x, y = y), color = pal[1]) +
  coord_cartesian(expand = F) +
  scale_x_continuous('Number of shocks') +
  scale_y_continuous('Probability') +
  theme_pubr()

gpd_fit <- fgpd(shock_df$shock_size)
gpd_x <- seq(1, max(shock_df$shock_size), 0.01)
gpd_y <- dgpd(gpd_x, gpd_fit$u, gpd_fit$sigmau, gpd_fit$xi, gpd_fit$phiu)

p2 <- ggplot(shock_df) +
  geom_histogram(aes(x = shock_size, y = ..density..), fill = 'gray60')+
  geom_line(data = data.frame(x = gpd_x, y = gpd_y), aes(x=x, y=y), color = pal[2])+
  coord_cartesian(expand = F) +
  scale_x_continuous('Size of shocks') +
  scale_y_continuous('Probability') +
  theme_pubr()

p1+p2+ plot_annotation(tag_levels = 'a')
ggsave("plots/stat_dist_shocks.png", width = 20, height = 10, units = 'cm')


## Plot heatmaps of spatial shocks
shock_heatmap_df <-read_csv("output/shock_heatmap_df.csv")

heatmap1 <- shock_heatmap_df |>
  pivot_longer(c('rv', 'r')) |>
  mutate(name = factor(name, levels = c('rv', 'r'), labels = c('Without shocks', 'With shocks'))) |>
  ggplot(aes(x = x, y = y)) +
  geom_tile(aes(fill = value), color = 'gray50') +
  theme_void()+
  facet_wrap(~name) +
  scale_fill_viridis_c(option = 'magma', direction = -1) +
  theme(legend.position = 'bottom', legend.title = element_blank())
heatmap1

heatmap2 <- shock_heatmap_df |>
  mutate(sl = factor(as.character(sl), levels = c('0','1'), labels = c('No shock', 'Shock'))) |>
  ggplot(aes(x = x, y = y)) +
  geom_tile(aes(fill = sl), color = 'gray50') +
  geom_point(data = shock_heatmap_df |> filter(so == 1), color = 'white')+
  theme_void()+
  scale_fill_manual("",values = c(pal[7],pal[1]))
  #theme(legend.position = 'none', legend.title = element_blank())
heatmap2

heatmap3 <- ggplot(shock_heatmap_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = r), color = 'gray50') +
  theme_void()+
  scale_fill_viridis_c(option = 'magma', direction = -1) +
  theme(legend.position = 'bottom', legend.title = element_blank())
heatmap3

heatmap_shock <- heatmap1 + heatmap2 + plot_layout(guides = 'collect', widths = c(2,1)) & theme(legend.position = 'bottom')
heatmap_shock
ggsave("plots/spatial_shock_heatmap.png", heatmap_shock, units = 'cm', width = 20, height = 9)
