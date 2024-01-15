library(ggplot2)
library(ggpubr)
library(evmix)

T = 100

set.seed(200)

zeta = rt(T, 3)
t = 1:T
eta = qt(1-0.1,3)
df <- data.frame(t=t, zeta=zeta)
exceed <- df[df$zeta > eta,]

p1 <- ggplot(df, aes(x = t, y = zeta)) +
  geom_hline(yintercept = eta, linetype = 2) +
  geom_line() +
  geom_point(data= exceed, color = 'red')+
  theme_pubr()

# Plot correspondence to student T distribution
shock_sizes <- c()
nshocks <- c()
nsims = 1600
for (i in 1:nsims) {
  zeta = rt(T, 3)
  t = 1:T
  eta = qt(1-0.01,3)
  df <- data.frame(t=t, zeta=zeta)
  exceed <- df[df$zeta > eta,]
  nshocks = c(nshocks, sum(df$zeta > eta))
  shock_sizes <- c(shock_sizes, exceed$zeta - eta)
}
binom_pdf <- data.frame(
  x = 0:max(nshocks),
  y = dbinom(0:max(nshocks), 100, 1/100)
)

pal <- ggsci::pal_nejm()(8)
p2 <- ggplot(data.frame(nshocks = nshocks)) +
  geom_histogram(aes(x = nshocks, y = after_stat(count)/sum(after_stat(count))), binwidth = 0.5) +
  geom_line(data = binom_pdf, aes(x = x, y = y), color = pal[1]) +
  geom_point(data = binom_pdf, aes(x = x, y = y), color = pal[1]) +
  scale_y_continuous("Proportion") +
  scale_x_continuous("Number of shocks")+
  theme_pubr()

gpd_fit <- fgpd(shock_sizes)
gpd_x <- seq(1, max(shock_sizes), 0.01)
gpd_y <- dgpd(gpd_x, gpd_fit$u, gpd_fit$sigmau, gpd_fit$xi, gpd_fit$phiu)
p3 <- ggplot(data.frame(shock_sizes = shock_sizes), aes(x = shock_sizes)) +
  geom_histogram(aes(y = after_stat(density))) +
  geom_line(data = data.frame(x = gpd_x, y = gpd_y), aes(x=x, y=y), color = pal[6])+
  theme_pubr()+
  scale_x_continuous('Size of shocks') +
  scale_y_continuous('Probability') +
  theme_pubr()

ggsave("plots/exceedance_plot.png", p2 + p3, width = 2000, height = 1000, units = 'px')
