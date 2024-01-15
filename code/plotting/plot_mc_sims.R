library(tidyverse)
library(patchwork)

nsims = 5

options(scipen = 999)
param_table <- read_csv(sprintf("output/param_df_%s.csv", nsims)) %>%
  mutate(η = as.character(η))
param_table <- param_table[1:(nrow(param_table)),] %>%
  mutate(N = ifelse(dims=='(10, 10)', 100, ifelse(dims=='(30, 30)', 900, 400)))
param_table <- param_table[1:26,]
param_table$run <- 1:nrow(param_table)

runs = 1:nrow(param_table)
alpha_vec <- c(0, 5, 10)
alpha_lab <- c(`0` = 'θ=0', `10` = 'θ=10', `20` = 'θ=20')

ce_df <- lapply(runs, \(i) read_csv(paste0("output/ce_df_param_search_", i, "_", nsims, ".csv"), show_col_types = F)) %>%
  lapply(\(x) filter(x, alpha %in% alpha_vec)) %>%
  lapply(\(x) dplyr::select(x, 'cvar_ev', 'mstd_ev', 'var', 'alpha')) %>%
  lapply(\(x) pivot_longer(x, cols = c(cvar_ev, mstd_ev))) %>%
  lapply(\(x) pivot_wider(x, names_from = 'var', values_from = 'value'))

names(ce_df) <- runs

ce_df_comb <- ce_df %>%
  bind_rows(.id = 'run') %>%
  mutate(run = as.integer(run)) %>%
  left_join(param_table, by = 'run') %>%
  mutate(N = ifelse(dims=='(10, 10)', '100', ifelse(dims=='(30, 30)', '900', '400')))

fcn_generate_mc_lineplots <- function(variable, axis_title = NULL, hide_y = TRUE) {
  selected_runs <- param_table %>%
    filter(run == 1 | !!as.symbol(variable) != param_table[1,variable][[1]])
  if (is.null(axis_title)) axis_title <- variable
  p <- ce_df_comb %>%
    mutate(name = factor(name, c('mstd_ev', 'cvar_ev'), c('M-SD', 'M-CVaR'))) %>%
    filter(run %in% selected_runs$run) %>%
    mutate(sig_diff = ifelse(lb > 0, 'cvar', ifelse(ub < 0, 'mstd', 'tie'))) %>%
    mutate(sig_diff = factor(sig_diff, levels = c('tie', 'cvar', 'mstd'))) %>%
    ggplot(aes(x = as.factor(!!as.symbol(variable)), color = name)) +
    geom_hline(yintercept = 0) +
    geom_errorbar(aes(ymin = lb, ymax = ub), width = 0, linewidth = 0.5, position = position_dodge(width=.5)) +
    geom_point(aes(y = median), position = position_dodge(width=.5)) +
    scale_y_continuous("CE change from EV", labels = scales::percent) +
    ggsci::scale_color_nejm() +
    #scale_color_manual(values = c('#333333',"#0072B5FF",'#BC3C29FF')) +
    facet_wrap(vars(alpha), strip.position="top", ncol = 1, labeller = labeller(alpha = alpha_lab)) +
    theme_bw() +
    labs(x = axis_title) +
    theme(panel.grid = element_blank(),legend.title = element_blank())
  #if (variable == 'η') p <- p + scale_x_continuous(trans=scales::pseudo_log_trans(base = 10))
  if (hide_y) p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
  if (variable=='η') p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
}

plt <- fcn_generate_mc_lineplots('ρ', hide_y = F) + fcn_generate_mc_lineplots('σ') +fcn_generate_mc_lineplots('η', 'π') +fcn_generate_mc_lineplots('yy', 'υ') +fcn_generate_mc_lineplots('budget', 'B') + fcn_generate_mc_lineplots('N') + fcn_generate_mc_lineplots('β') +
  plot_layout(nrow = 1, guides = 'collect')

ggsave('plots/mc_sim_plot_ev.png', plt, scale = 2, width = 2000, height = 1000, units = 'px')

