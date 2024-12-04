library(data.table)
library(dplyr)
library(tidyr)
library(scales)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(scales)

device_format <- c('pdf', 'png')
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path('simulation', 'generic', 'plotters', 'plot_bars_and_lines.R'))
source(file.path("simulation", "load_paths.R"))
theme_set(theme_bw())
customTheme <- f_getCustomTheme()

##-------------------------------------
## Custom functions
##-------------------------------------
f_custom_lineplot <- function(dat, yvar = 'protective_efficacy') {
  gg1 = ggplot(dat, aes_string(x = 'year', y = yvar)) +
    geom_point(aes(col = rtss_scen), size = 2) +
    geom_line(aes(col = rtss_scen), size = 1) +
    scale_color_manual(values = c('#fdb863', '#e66101', '#b2abd2', '#5e3c99', 'white')) +
    customTheme +
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
    geom_bar(aes(fill = rtss_scen), stat = 'identity', position = position_dodge(width = 1)) +
    scale_fill_manual(values = c('#fdb863', '#e66101', '#b2abd2', '#5e3c99', 'white')) +
    customTheme +
    geom_hline(yintercept = 0) +
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


##-------------------------------------
## Custom objects and experiment setup
##-------------------------------------
simout_dir <- file.path(projectpath, 'simulation_output', 'generic')
if (!dir.exists(paste0(simout_dir, '/_plots'))) dir.create(paste0(simout_dir, '/_plots'))

plot_dir <- file.path(simout_dir, '_plots_booster')
if (!dir.exists(plot_dir))dir.create(plot_dir)
plot_dir <- file.path(plot_dir, '_fig_for_20211104_RTSSdose_timing')
if (!dir.exists(plot_dir))dir.create(plot_dir)

age_label_values = seq(0, 8, length.out = 5)
age_labels = paste0(age_label_values, '-', (age_label_values + 1))

exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'cm_coverage', 'rtss_coverage', 'rtss_mode', 'smc_coverage',
                'intervention_correlation', 'Cohort_birth_month')
max_years = c(1, 2, 5, 10)
age_groups = paste0('U', max_years)

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
  f_add_rtss_cols() %>%
  f_add_rtss_scen()

cases_scenarios_by_agegroup <- datList_by_agegroup %>%
  bind_rows() %>%
  f_add_rtss_cols() %>%
  f_add_rtss_scen()

cases_scenarios_by_agegroup$age_group = factor(cases_scenarios_by_agegroup$age_group,
                                               levels = c("U1", "U2", "U5", "U10"),
                                               labels = c("U1", "U2", "U5", "U10"))

table(cases_scenarios_each_year$rtss_scen)
selected_scenarios <- c("Campaign 7m 2-booster (min age 9m)",
                        "Campaign 7m 3-booster (min age 9m)",
                        "Campaign 9m 2-booster (min age 24m)",
                        "Campaign 9m 3-booster (min age 24m)")

cases_scenarios_each_year = cases_scenarios_each_year %>%
  filter(Annual_EIR == 30,
         cm_coverage == 0.6,
         smc_coverage == 0,
         rtss_cov == 'strict',
         rtss_scen %in% selected_scenarios)

cur_df = cases_scenarios_each_year %>%
  filter(rtss_coverage == 0.8)
cur_df_counter = cases_scenarios_each_year %>%
  filter(rtss_coverage == 0)


cur_df_agegrp = cases_scenarios_by_agegroup %>%
  filter(Annual_EIR == 30,
         cm_coverage == 0.6,
         smc_coverage == 0,
         rtss_cov == 'strict',
         rtss_coverage == 0.8,
         rtss_scen %in% selected_scenarios)

##-------------------------------------
## Figures
##-------------------------------------
gg0A <- ggplot(data = subset(cur_df, rtss_scen %in% selected_scenarios), aes(x = year, y = clinical_cases * 100)) +
  geom_point(aes(col = as.factor(rtss_scen), group = rtss_scen), size = 2) +
  geom_line(aes(col = as.factor(rtss_scen), group = rtss_scen), size = 1) +
  geom_point(data = cur_df_counter, aes(x = year, y = clinical_cases * 100), col = 'black', size = 2) +
  geom_line(data = cur_df_counter, aes(x = year, y = clinical_cases * 100), col = 'black', size = 1) +
  scale_color_manual(values = c('#fdb863', '#e66101', '#b2abd2', '#5e3c99', 'black')) +
  scale_x_continuous(breaks = age_label_values, labels = age_labels) +
  scale_y_continuous(labels = comma) +
  labs(title = 'Clinical cases', y = 'clinical cases (per 100,000)',
       x = 'age of child', linetype = '', color = 'RTS,S schedule', caption = '') +
  customTheme +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = Inf) +
  geom_vline(xintercept = Inf) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gg0B <- ggplot(data = subset(cur_df, rtss_scen %in% selected_scenarios), aes(x = year, y = severe_cases * 100)) +
  geom_point(aes(col = as.factor(rtss_scen), group = rtss_scen), size = 2) +
  geom_line(aes(col = as.factor(rtss_scen), group = rtss_scen), size = 1) +
  geom_point(data = cur_df_counter, aes(x = year, y = severe_cases * 100), col = 'black', size = 2) +
  geom_line(data = cur_df_counter, aes(x = year, y = severe_cases * 100), col = 'black', size = 1) +
  scale_color_manual(values = c('#fdb863', '#e66101', '#b2abd2', '#5e3c99', 'black')) +
  scale_x_continuous(breaks = age_label_values, labels = age_labels) +
  scale_y_continuous(labels = comma) +
  labs(title = 'Severe cases', y = 'severe cases (per 100,000)',
       x = 'age of child', linetype = '', color = 'RTS,S schedule', caption = '') +
  customTheme +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = Inf) +
  geom_vline(xintercept = Inf) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ll <- get_legend(gg0A)
gg0A <- gg0A + theme(legend.position = 'None')
gg0B <- gg0B + theme(legend.position = 'None')
gg <- plot_grid(gg0A, gg0B, labels = c("A", "B"))
gg <- plot_grid(gg, ll, rel_widths = c(1, 0.4))

f_save_plot(gg, paste0('cases_by_age'), file.path(plot_dir),
            width = 13, height = 4, units = 'in', device_format = device_format)


gg1 <- f_custom_lineplot(dat = cur_df, yvar = 'cases_averted_per100000') +
  labs(title = 'Clinical cases by RTS,S', y = 'clinical cases averted (per 100,000)',
       x = 'age of child',    color = 'RTS,S schedule',  caption = '')

gg2 = f_custom_lineplot(dat = cur_df, yvar = 'severe_cases_averted_per100000') +
  labs(title = 'Severe cases averted by RTS,S', y = 'severe cases averted (per 100,000)',
       x = 'age of child',  color = 'RTS,S schedule', caption = 'EIR=30, no SMC, CM=60%, high-unimodal seasonality')

ll <- get_legend(gg1)
gg1 <- gg1 + theme(legend.position = 'None')
gg2 <- gg2 + theme(legend.position = 'None')
gg <- plot_grid(gg1, gg2)
gg_line <- plot_grid(gg, ll, rel_widths = c(1, 0.4))

f_save_plot(gg_line, paste0('cases_averted_by_age'), file.path(plot_dir),
            width = 11, height = 4, units = 'in', device_format = device_format)

##### PE
gg1 = f_custom_lineplot(dat = cur_df, yvar = 'protective_efficacy') +
  labs(title = 'PE for clinical cases by RTS,S',
       y = 'protective_efficacy', x = 'age of child', color = 'RTS,S schedule', caption = '')

gg2 = f_custom_lineplot(dat = cur_df, yvar = 'protective_efficacy_severe') +
  labs(title = 'PE for severe cases by RTS,S',
       y = 'protective_efficacy',  x = 'age of child',  color = 'RTS,S schedule',
       caption = 'EIR=30, no SMC, CM=60%, high-unimodal seasonality')

ll <- get_legend(gg1)
gg1 <- gg1 + theme(legend.position = 'None')
gg2 <- gg2 + theme(legend.position = 'None')
gg <- plot_grid(gg1, gg2)
gg <- plot_grid(gg, ll, rel_widths = c(1, 0.4))

f_save_plot(gg, paste0('PE_by_age'), file.path(plot_dir),
            width = 11, height = 4, units = 'in', device_format = device_format)


######################################
### CUMULATIVE
gg1 = f_custom_barplot(cur_df_agegrp, yvar = 'cases_averted_per100000') +
  scale_y_continuous(lim = c(0, 100000), expand = c(0, 0), labels = comma) +
  labs(title = 'Clinical cases averted by RTS,S',
       y = 'clinical cases averted (per 100,000)',
       x = 'age of child', fill = 'RTS,S schedule',  caption = '')

gg2 = f_custom_barplot(cur_df_agegrp, yvar = 'severe_cases_averted_per100000') +
  scale_y_continuous(lim = c(0, 2500), expand = c(0, 0), labels = comma) +
  labs(title = 'Severe cases averted by RTS,S',
       y = 'severe cases averted (per 100,000)',   x = 'age of child',
       fill = 'RTS,S schedule',  caption = 'EIR=30, no SMC, CM=60%, high-unimodal seasonality')


ll <- get_legend(gg1)
gg1 <- gg1 + theme(legend.position = 'None')
gg2 <- gg2 + theme(legend.position = 'None')
gg <- plot_grid(gg1, gg2)
ggbar <- plot_grid(gg, ll, rel_widths = c(1, 0.4))
f_save_plot(ggbar, paste0('cases_averted_by_agegroup'), file.path(plot_dir),
            width = 11, height = 4, units = 'in', device_format = device_format)


ggcombo <- plot_grid(gg_line, ggbar, ncol = 1, labels = c("A", "B"))
f_save_plot(ggcombo, paste0('cases_averted_by_age_and_agegroup'), file.path(plot_dir),
            width = 14, height = 11, units = 'in', device_format = device_format)


analyze_EPI = T
if (analyze_EPI) {

  cases_scenarios_each_year <- datList_each_year %>%
    bind_rows() %>%
    f_add_rtss_cols() %>%
    f_add_rtss_scen() %>%
    filter(rtss_type == 'EPI',
           Annual_EIR == 30,
           rtss_coverage == 0.8,
           cm_coverage == 0.6,
           smc_coverage == 0)

  cases_scenarios_by_agegroup <- datList_by_agegroup %>%
    bind_rows() %>%
    f_add_rtss_cols() %>%
    f_add_rtss_scen() %>%
    filter(rtss_type == 'EPI',
           rtss_coverage == 0.8,
           cm_coverage == 0.6,
           smc_coverage == 0)

  cases_scenarios_by_agegroup$age_group = factor(cases_scenarios_by_agegroup$age_group,
                                                 levels = c("U1", "U2", "U5", "U10"),
                                                 labels = c("U1", "U2", "U5", "U10"))

  pp1 <- f_custom_lineplot(dat = cases_scenarios_each_year, yvar = 'cases_averted_per100000') +
    labs(title = 'Clinical cases averted by RTS,S',
         y = 'clinical cases averted\n(per 100,000)',
         x = ' ', fill = 'RTS,S schedule',   caption = '')

  pp2 <- f_custom_lineplot(dat = cases_scenarios_each_year, yvar = 'severe_cases_averted_per100000') +
    labs(title = 'Severe cases averted by RTS,S',
         y = 'severe cases averted\n(per 100,000)',
         x = ' ', fill = 'RTS,S schedule', caption = 'EIR=30, no SMC, CM=60%, high-unimodal seasonality')

  ll <- get_legend(pp1)
  pp1 <- pp1 + theme(legend.position = 'None')
  pp2 <- pp2 + theme(legend.position = 'None')
  pp <- plot_grid(pp1, pp2, ncol = 2, labels = c("A", "B"))
  pp <- plot_grid(pp, ll, rel_widths = c(1, 0.4))

  f_save_plot(pp, paste0('cases_averted_lineplot_byage_EIR30'), file.path(plot_dir),
              width = 12, height = 4, units = 'in', device_format = device_format)


  plot_dat <- cases_scenarios_by_agegroup %>% filter(age_group == "U5")
  plot_dat$Annual_EIR <- as.factor(plot_dat$Annual_EIR)
  pp1 <- f_custom_barplot(dat = plot_dat, xvar = 'Annual_EIR', yvar = 'cases_averted_per100000') +
    labs(title = 'Clinical cases averted by RTS,S\nin children under the age of 5 years',
         y = 'clinical cases averted\n(per 100,000)',  x = 'Annual EIR',
         fill = 'RTS,S schedule', caption = '')

  pp2 <- f_custom_barplot(dat = plot_dat, xvar = 'Annual_EIR', yvar = 'severe_cases_averted_per100000') +
    labs(title = 'Severe cases averted by RTS,S\nin children under the age of 5 years',
         y = 'severe cases averted\n(per 100,000)',
         x = 'Annual EIR',  fill = 'RTS,S schedule',
         caption = 'Uage=5, no SMC, CM=60%, high-unimodal seasonality')

  ll <- get_legend(pp1)
  pp1 <- pp1 + theme(legend.position = 'None')
  pp2 <- pp2 + theme(legend.position = 'None')
  pp <- plot_grid(pp1, pp2, ncol = 2, labels = c("A", "B"))
  pp <- plot_grid(pp, ll, rel_widths = c(1, 0.4))

  f_save_plot(pp, paste0('cases_averted_barplot_byEIR_U5'), file.path(plot_dir),
              width = 12, height = 6, units = 'in', device_format = device_format)


  plot_dat <- cases_scenarios_by_agegroup %>% filter(Annual_EIR == 30)
  plot_dat$Annual_EIR <- as.factor(plot_dat$Annual_EIR)
  pp1 <- f_custom_barplot(dat = plot_dat, xvar = 'age_group', yvar = 'cases_averted_per100000') +
    labs(title = 'Clinical cases averted by RTS,S\nin children under the age of 5 years',
         y = 'clinical cases averted\n(per 100,000)',
         x = 'age_group', fill = 'RTS,S schedule', caption = '')

  pp2 <- f_custom_barplot(dat = plot_dat, xvar = 'age_group', yvar = 'severe_cases_averted_per100000') +
    labs(title = 'Severe cases averted by RTS,S\nin children under the age of 5 years',
         y = 'severe cases averted\n(per 100,000)',
         x = 'age_group', fill = 'RTS,S schedule',  caption = '')

  ll <- get_legend(pp1)
  pp1 <- pp1 + theme(legend.position = 'None')
  pp2 <- pp2 + theme(legend.position = 'None')
  pp <- plot_grid(pp1, pp2, ncol = 2, labels = c("A", "B"))
  pp <- plot_grid(pp, ll, rel_widths = c(1, 0.4))

  f_save_plot(pp, paste0('cases_averted_barplot_byAge_EIR30'), file.path(plot_dir),
              width = 12, height = 6, units = 'in', device_format = device_format)

}