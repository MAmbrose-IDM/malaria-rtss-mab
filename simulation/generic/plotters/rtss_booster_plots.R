library(data.table)
library(dplyr)
library(tidyr)
library(scales)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(scales)
# wdir at rtss-scenarios
device_format <- c('pdf', 'png')
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path('simulation', 'generic', 'plotters', 'plot_bars_and_lines.R'))
source(file.path("simulation", "load_paths.R"))

age_label_values = seq(0, 8, length.out = 5)
age_labels = paste0(age_label_values, '-', (age_label_values + 1))

outcomes <- c('protective_efficacy', 'cases_averted_per100000', 'rtss_protective_efficacy', 'rtss_cases_averted_per100000')
outcomes_severe <- c('protective_efficacy_severe', 'severe_cases_averted_per100000', 'rtss_protective_efficacy_severe', 'rtss_severe_cases_averted_per100000')


f_custom_lineplot <- function(dat, xvar = 'year', yvar = 'protective_efficacy') {
  gg1 = ggplot(dat, aes_string(x = xvar, y = yvar)) +
    geom_line(aes(col = as.factor(rtss_scenario), linetype = rtss_cov, group = interaction(rtss_mode, rtss_scenario, rtss_cov)), size = 1) +
    facet_wrap(rtss_type ~ rtss_nboosters, scales = 'free_x', nrow = 1) +
    scale_color_brewer(palette = 'Dark2') +
    scale_linetype_manual(values = c('dashed', 'solid')) +
    f_getCustomTheme() +
    geom_hline(yintercept = 0) +
    scale_x_continuous(breaks = age_label_values, labels = age_labels) +
    scale_y_continuous(labels = comma) +
    geom_hline(yintercept = Inf) +
    geom_vline(xintercept = Inf) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(gg1)
}

f_custom_barplot <- function(dat, xvar = 'age_group', yvar = 'cases_averted_per100000') {
  gg1 = ggplot(dat, aes_string(x = xvar, y = yvar)) +
    geom_bar(aes(fill = as.factor(rtss_scenario),
                 group = interaction(rtss_mode, rtss_scenario, rtss_cov)), size = 1,
             stat = 'identity', position = position_dodge(width = 1)) +
    facet_wrap(rtss_type ~ rtss_nboosters, scales = 'free_x', nrow = 1) +
    scale_fill_brewer(palette = 'Dark2') +
    f_getCustomTheme() +
    geom_hline(yintercept = 0) +
    scale_y_continuous(labels = comma) +
    geom_hline(yintercept = Inf) +
    geom_vline(xintercept = Inf) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(gg1)
}

f_add_rtss_cols <- function(dat) {
  dat$rtss_type = 'EPI'
  dat$rtss_type[(grep('campboost', dat$rtss_mode))] = 'Campaign'

  dat$rtss_3rddose = '9m'
  dat$rtss_3rddose[(grep('7m', dat$rtss_mode))] = '7m'

  dat$rtss_cov = 'strict'
  dat$rtss_cov[(grep('flexible', dat$rtss_mode))] = 'flexible'

  dat$rtss_nboosters <- 1
  dat$rtss_nboosters[(grep('noboost', dat$rtss_mode))] = 0
  dat$rtss_nboosters[(grep('constant', dat$rtss_mode))] = 1
  dat$rtss_nboosters[(grep('campboost2', dat$rtss_mode))] = 2
  dat$rtss_nboosters[(grep('campboost3', dat$rtss_mode))] = 3

  dat$rtss_scenario <- paste0(dat$rtss_3rddose, '_', dat$minBoostAge)
  dat$rtss_schedule <- dat$rtss_scenario

  return(dat)
}

f_add_rtss_scen <- function(dat) {
  dat$rtss_scen <- paste0(dat$rtss_scenario, dat$rtss_type, '_', dat$rtss_nboosters)
  dat$rtss_nboosters <- as.factor(dat$rtss_nboosters)
  dat$minBoostAge <- as.factor(dat$minBoostAge)
  dat$rtss_scen <- factor(dat$rtss_scen,
                          levels = c("7m_274Campaign_1", "7m_274Campaign_2", "7m_274Campaign_3",
                                     "9m_547Campaign_1", "9m_547Campaign_2", "9m_547Campaign_3",
                                     "9m_730Campaign_1", "9m_730Campaign_2", "9m_730Campaign_3",
                                     "7m_274EPI_0", "7m_274EPI_1", "9m_730EPI_0", "9m_730EPI_1"),
                          labels = c("Campaign 7m 1-booster (min age 9m)", "Campaign 7m 2-booster (min age 9m)", "Campaign 7m 3-booster (min age 9m)",
                                     "Campaign 9m 1-booster (min age 18m)", "Campaign 9m 2-booster (min age 18m)", "Campaign 9m 3-booster (min age 18m)",
                                     "Campaign 9m 1-booster (min age 24m)", "Campaign 9m 2-booster (min age 24m)", "Campaign 9m 3-booster (min age 24m)",
                                     "EPI 7m 0-booster", "EPI 7m 1-booster", "EPI 9m 0-booster", "EPI 9m 1-booster"))
  dat$rtss_scenario[dat$rtss_type == 'EPI'] <- 'EPI'
  dat$rtss_scenario <- factor(dat$rtss_scenario,
                              levels = c('7m_274', '9m_547', '9m_730', 'EPI'),
                              labels = c('Campaign 7m\n(booster min age 9m)',
                                         'Campaign 9m\n(booster min age 18m)',
                                         'Campaign 9m\n(booster min age 24m)', 'EPI'))

  return(dat)
}

simout_dir <- file.path(projectpath, 'simulation_output', 'generic')
if (!dir.exists(paste0(simout_dir, '/_plots'))) dir.create(paste0(simout_dir, '/_plots'))

exp_sweeps <- c('Scenario_id', 'Run_Number', 'Annual EIR', 'seasonality',
                'cm_coverage', 'rtss_coverage', 'rtss_mode', 'smc_coverage', 'Cohort_birth_month',
                'cm_target_group', 'smc_target_group', 'rtss_target_group', 'minBoostAge')

max_years = c(1, 2, 5, 10)
age_groups = paste0('U', max_years)
seasonality_levels = c('constant', 'moderate_unimodal', 'high_unimodal')

theme_set(theme_bw())
customTheme <- f_getCustomTheme()
plot_dir <- file.path(simout_dir, '_plots_booster')
if (!dir.exists(plot_dir))dir.create(plot_dir)


exp_names <- c('generic_campboost_sweep_highseasonal_SMC', 'generic_campboost_flexible_highseasonal_SMC',
               'generic_campboost_7m_sweep_highseasonal_SMC', 'generic_campboost_sweep_highseasonal_SMC_minboostage18')
exp_name_counterfactual <- 'generic_counterfactual_noRTSS_SMC'

datList_each_year <- list()
datList_by_agegroup <- list()
for (exp_name in exp_names) {
  sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                      exp_name_counterfactual = exp_name_counterfactual, add_PE_perAge = TRUE,
                                      max_years = max_years, keep_birth_month = FALSE)
  datList_each_year[[length(datList_each_year) + 1]] = sim_output[[4]] %>%
    mutate(exp_name = exp_name) %>%
    as.data.frame()

  datList_by_agegroup[[length(datList_by_agegroup) + 1]] = sim_output[[3]] %>%
    mutate(exp_name = exp_name) %>%
    as.data.frame()
}


cases_scenarios_each_year <- datList_each_year %>%
  bind_rows() %>%
  f_add_rtss_cols()

cases_scenarios_by_agegroup <- datList_by_agegroup %>%
  bind_rows() %>%
  f_add_rtss_cols()

cases_scenarios_by_agegroup$age_group = factor(cases_scenarios_by_agegroup$age_group,
                                               levels = c("U1", "U2", "U5", "U10"),
                                               labels = c("U1", "U2", "U5", "U10"))

table(cases_scenarios_each_year$rtss_mode)
table(cases_scenarios_each_year$rtss_scenario, cases_scenarios_each_year$rtss_mode, exclude = NULL)
table(cases_scenarios_each_year$minBoostAge, cases_scenarios_each_year$rtss_mode, exclude = NULL)
table(cases_scenarios_each_year$rtss_nboosters, cases_scenarios_each_year$rtss_mode, exclude = NULL)
table(cases_scenarios_each_year$rtss_type, cases_scenarios_each_year$rtss_mode, exclude = NULL)


cur_df <- cases_scenarios_each_year %>%
  filter(Annual_EIR == 30 &
           rtss_coverage == 0.8 &
           smc_coverage == 0)

cur_df_agegrp <- cases_scenarios_by_agegroup %>%
  filter(Annual_EIR == 30 &
           rtss_coverage == 0.8 &
           smc_coverage == 0)

# create plot
gg0A <- cases_scenarios_each_year %>%
  filter(Annual_EIR == 30 &
           rtss_coverage == 0.8 &
           smc_coverage == 0) %>%
  mutate(severe_cases_per100000 = severe_cases * 100) %>%
  f_custom_lineplot(yvar = 'severe_cases_per100000') +
  labs(title = 'Severe cases', y = 'cases (per 100,000)',
       x = 'age of child', linetype = '', color = 'RTS,S schedule', caption = '')

gg0B <- cases_scenarios_each_year %>%
  filter(Annual_EIR == 30 &
           rtss_coverage == 0.8 &
           smc_coverage == 0) %>%
  mutate(clinical_cases_per100000 = clinical_cases * 100) %>%
  f_custom_lineplot(yvar = 'clinical_cases_per100000') +
  labs(title = 'Clinical cases', y = 'cases (per 100,000)',
       x = 'age of child', linetype = '', color = 'RTS,S schedule', caption = '')

ll <- get_legend(gg0A)
gg0A <- gg0A + theme(legend.position = 'None')
gg0B <- gg0B + theme(legend.position = 'None')
gg <- plot_grid(gg0A, gg0B, ncol = 1, labels = c("A", "B"))
gg <- plot_grid(gg, ll, rel_widths = c(1, 0.25))
gg

f_save_plot(gg, paste0('cases_by_age'), file.path(plot_dir),
            width = 13, height = 6, units = 'in', device_format = device_format)


### Loop through outcome channels
for (i in c(1:4)) {
  outcome = outcomes[i]
  outcome_severe = outcomes_severe[i]

  gg0A <- cur_df %>%
    f_custom_lineplot(yvar = outcome) +
    labs(title = 'Clinical cases', y = outcome,
         x = 'age of child', linetype = '', color = 'RTS,S schedule', caption = '')

  gg0B <- cur_df %>%
    f_custom_lineplot(yvar = outcome_severe) +
    labs(title = 'Severe cases', y = outcome_severe,
         x = 'age of child', linetype = '', color = 'RTS,S schedule', caption = '')

  ll <- get_legend(gg0A)
  gg0A <- gg0A + theme(legend.position = 'None')
  gg0B <- gg0B + theme(legend.position = 'None')
  gg <- plot_grid(gg0A, gg0B, ncol = 1, labels = c("A", "B"))
  gg <- plot_grid(gg, ll, rel_widths = c(1, 0.25))
  gg

  f_save_plot(gg, paste0(outcome, '_rtssmodes_EIR30'), file.path(plot_dir),
              width = 13, height = 6, units = 'in', device_format = device_format)

  ### Barplots
  gg0A <- cur_df_agegrp %>%
    filter(age_group == 'U5' & rtss_cov == 'strict') %>%
    f_custom_barplot(xvar = 'age_group', yvar = outcome) +
    labs(title = 'Clinical cases', y = outcome,
         x = 'age of child', linetype = '', fill = 'RTS,S schedule', caption = '')

  gg0B <- cur_df_agegrp %>%
    filter(age_group == 'U5' & rtss_cov == 'strict') %>%
    f_custom_barplot(xvar = 'age_group', yvar = outcome_severe) +
    labs(title = 'Severe cases', y = outcome_severe,
         x = 'age of child', linetype = '', fill = 'RTS,S schedule', caption = '')

  ll <- get_legend(gg0A)
  gg0A <- gg0A + theme(legend.position = 'None')
  gg0B <- gg0B + theme(legend.position = 'None')
  gg <- plot_grid(gg0A, gg0B, ncol = 1, labels = c("A", "B"))
  gg <- plot_grid(gg, ll, rel_widths = c(1, 0.25))
  gg

  f_save_plot(gg, paste0(outcome, '_bar_rtssmodes_EIR30'), file.path(plot_dir),
              width = 13, height = 6, units = 'in', device_format = device_format)
}

cur_df <- f_add_rtss_scen(cur_df)
cur_df_agegrp <- f_add_rtss_scen(cur_df_agegrp)


## Summary plot
# '#fdd870', '#fdbf11', '#e88e2d', '#ca5800'
p1 <- ggplot(data = subset(cur_df, rtss_nboosters != '0' & rtss_cov == 'strict')) +
  geom_line(aes_string(x = 'year', y = outcomes[3], col = 'rtss_scen')) +
  scale_color_manual(values = c('#e46aa7', '#ec008b', '#af1f6b',
                                '#73bfe2', '#1696d2', '#12719e',
                                '#78c26d', '#408941', '#2c5c2d',
                                '#fdbf11', '#ca5800')) +
  facet_wrap(~rtss_scenario, nrow = 1) +
  f_getCustomTheme() +
  geom_hline(yintercept = 0) +
  scale_x_continuous(breaks = age_label_values, labels = age_labels) +
  f_getCustomTheme() +
  scale_y_continuous(labels = comma) +
  geom_hline(yintercept = Inf) +
  geom_vline(xintercept = Inf) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p2 <- ggplot(data = subset(cur_df, rtss_nboosters != '0' & rtss_cov == 'strict')) +
  geom_line(aes_string(x = 'year', y = outcomes_severe[3], col = 'rtss_scen')) +
  scale_color_manual(values = c('#e46aa7', '#ec008b', '#af1f6b',
                                '#73bfe2', '#1696d2', '#12719e',
                                '#78c26d', '#408941', '#2c5c2d',
                                '#fdbf11', '#ca5800')) +
  facet_wrap(~rtss_scenario, nrow = 1) +
  f_getCustomTheme() +
  geom_hline(yintercept = 0) +  # geom_vline(xintercept =  4.5) +
  scale_x_continuous(breaks = age_label_values, labels = age_labels) +
  scale_y_continuous(labels = comma) +
  geom_hline(yintercept = Inf) +
  geom_vline(xintercept = Inf) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


b1 <- ggplot(data = subset(cur_df_agegrp, age_group %in% c('U5') &
  rtss_nboosters != '0' &
  rtss_cov == 'strict')) +
  geom_bar(aes_string(x = 'rtss_scen', y = outcomes[3], fill = 'rtss_scen'), stat = 'identity', position = position_dodge(width = 1)) +
  scale_fill_manual(values = c('#e46aa7', '#ec008b', '#af1f6b',
                               '#73bfe2', '#1696d2', '#12719e',
                               '#78c26d', '#408941', '#2c5c2d',
                               '#fdbf11', '#ca5800')) +
  facet_wrap(~age_group, nrow = 2) +
  f_getCustomTheme() +
  geom_hline(yintercept = 0) +
  scale_y_continuous(labels = comma) +
  geom_hline(yintercept = Inf) +
  geom_vline(xintercept = Inf) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

b2 <- ggplot(data = subset(cur_df_agegrp, age_group %in% c('U5') &
  rtss_nboosters != '0' &
  rtss_cov == 'strict')) +
  geom_bar(aes_string(x = 'rtss_scen', y = outcomes_severe[3], fill = 'rtss_scen'), stat = 'identity', position = position_dodge(width = 1)) +
  scale_fill_manual(values = c('#e46aa7', '#ec008b', '#af1f6b',
                               '#73bfe2', '#1696d2', '#12719e',
                               '#78c26d', '#408941', '#2c5c2d',
                               '#fdbf11', '#ca5800')) +
  facet_wrap(~age_group, nrow = 2) +
  f_getCustomTheme() +
  geom_hline(yintercept = 0) +
  scale_y_continuous(labels = comma) +
  geom_hline(yintercept = Inf) +
  geom_vline(xintercept = Inf) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


gg_legend <- get_legend(b2)
p1 <- p1 + theme(legend.position = 'None')
p2 <- p2 + theme(legend.position = 'None')
b1 <- b1 + theme(legend.position = 'None')
b2 <- b2 + theme(legend.position = 'None')

pp <- plot_grid(p1, p2, ncol = 1)
bb <- plot_grid(b1, b2, ncol = 1)


gg <- plot_grid(pp, bb, ncol = 2, rel_widths = c(1, 0.5))
gg <- plot_grid(gg, gg_legend, ncol = 2, rel_widths = c(1, 0.3))
gg

f_save_plot(gg, paste0('rtss_booster_summary_noSMC_highseason_EIR30'), file.path(plot_dir),
            width = 15, height = 6, units = 'in', device_format = device_format)

