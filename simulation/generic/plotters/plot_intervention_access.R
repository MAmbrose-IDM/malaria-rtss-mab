library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(tidyr)

source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path('simulation', 'generic', 'plotters', 'plot_bars_and_lines.R'))
source(file.path("simulation", "load_paths.R"))

simout_dir <- file.path(projectpath, 'simulation_output', 'generic')
theme_set(theme_bw())
customTheme <- f_getCustomTheme()
myColors <- f_getColors()
max_years <- c(1, 2, 5, 10)


### New exp_name for saving plots
exp_name <- 'generic_accesscorrelation05_v3'
exp_filepath <- file.path(simout_dir, exp_name)
if (!dir.exists(exp_filepath))dir.create(exp_filepath)
plot_dir <- file.path(exp_filepath)
if (!dir.exists(plot_dir))dir.create(plot_dir)

##---------------------------------------------
## Load custom experiment scenario to combine
##---------------------------------------------
exp_sweeps <- c('Scenario_id', 'Run_Number', 'Annual EIR', 'seasonality',
                'cm_coverage', 'rtss_coverage', 'rtss_mode', 'smc_coverage', 'Cohort_birth_month',
                'cm_target_group', 'smc_target_group', 'rtss_target_group')

sim_output <- load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                     add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)

dat <- sim_output[[4]]
dat <- dat %>% filter(smc_coverage!=0)
dat <- dat %>% filter(rtss_coverage!=0)
#dat <- dat %>% filter(seasonality=='constant')
dat$access_scen <- paste0('cm_', dat$cm_target_group, ' rtss_', dat$rtss_target_group, ' smc_', dat$smc_target_group)
table(dat$access_scen)

##---------------------------------------------
## Data checking and variable formatting
##---------------------------------------------
table(dat$access_scen, dat$Annual_EIR)
table(dat$access_scen, dat$cm_coverage)
table(dat$access_scen, dat$rtss_coverage)
table(dat$access_scen, dat$smc_coverage)
table(dat$access_scen, dat$seasonality)

dat$rtss_coverage_fct <- factor(dat$rtss_coverage, levels = unique(dat$rtss_coverage),
                                labels = paste0('RTS,S ', unique(dat$rtss_coverage) * 100, '%'))
dat$smc_coverage_fct <- factor(dat$smc_coverage, levels = unique(dat$smc_coverage),
                               labels = paste0('SMC ', unique(dat$smc_coverage) * 100, '%'))
dat$scenario_name = paste0(dat$rtss_coverage_fct, '; ', dat$smc_coverage_fct)
dat$cm_coverage <- as.numeric(dat$cm_coverage)

dat$seasonality = factor(dat$seasonality, levels = c('constant', 'moderate_unimodal', 'high_unimodal'))
dat$Annual_EIR_fct <- factor(dat$Annual_EIR, levels = unique(dat$Annual_EIR),
                             labels = paste0('annual EIR ', unique(dat$Annual_EIR)))
dat$access_scen <- factor(dat$access_scen,
                          levels = c('cm_random rtss_random smc_high', 'cm_random rtss_high smc_high', 'cm_random rtss_low smc_high'),
                          labels = c('uncorrelated', 'correlated', 'anticorrelated'))
table(dat$access_scen)
dat <- dat %>% filter(access_scen!='rem')

age_label_values <- seq(0, 8, length.out = 5)
age_labels <- paste0(age_label_values, '-', (age_label_values + 1))

##---------------------------------------------
## Generate lineplot over age
##---------------------------------------------
f_custom_plot <- function(dat, yvar, facet_var, yvar_grey = '') {

  dat <- as.data.frame(dat)
  dat$facet_var <- dat[, facet_var]
  gg <- ggplot(data = dat) +
    geom_line(aes(x = year, y = get(yvar), col = access_scen, group = access_scen), size = 1) +
    theme_bw() +
    f_getCustomTheme() +
    #scale_color_manual(values = c("black", "#ce6e17", "#e29b54", "#009367", "#76c050","blue")) +
    scale_color_brewer(palette='Dark2') +
    scale_x_continuous(breaks = age_label_values, labels = age_labels) +
    geom_hline(yintercept = 0) +
    labs(y = yvar,
         x = 'age of child',
         col = 'Access correlation scenario') +
    facet_wrap(~facet_var, nrow = 3, scales = 'free_y')

  if (nchar(yvar_grey) > 1) {
    gg <- gg + geom_line(aes(x = year, y = get(yvar_grey), group = access_scen), col = 'grey', size = 1)
  }

  return(gg)

}


##-------
## Total case population
##---------
gg_severe <- f_custom_plot(dat, yvar = 'severe_cases', facet_var = 'Annual_EIR', yvar_grey = 'severe_cases_ref') +
  labs(title = 'Severe cases by EIR')

gg_clinical <- f_custom_plot(dat, yvar = 'clinical_cases', facet_var = 'Annual_EIR', yvar_grey = 'clinical_cases_ref') +
  labs(title = 'Clinical cases by EIR')

gg_legend <- get_legend(gg_clinical)
gg_clinical <- gg_clinical + theme(legend.position = 'None')
gg_severe <- gg_severe + theme(legend.position = 'None')
gg <- plot_grid(gg_clinical, gg_severe, nrow = 1)
gg <- plot_grid(gg, gg_legend, nrow = 1, rel_widths = c(1, 0.25))
gg
f_save_plot(gg, plot_name = "N_cases_byEIR", width = 12, height = 10, plot_dir = file.path(plot_dir))


##-------
## PE SMC+RTSS
##---------
gg_severe <- f_custom_plot(dat, yvar = 'protective_efficacy_severe', facet_var = 'Annual_EIR') +
  labs(title = 'Protective efficacy in severe cases by EIR')

gg_clinical <- f_custom_plot(dat, yvar = 'protective_efficacy', facet_var = 'Annual_EIR') +
  labs(title = 'Protective efficacy in clinical cases by EIR')

gg_legend <- get_legend(gg_clinical)
gg_clinical <- gg_clinical + theme(legend.position = 'None')
gg_severe <- gg_severe + theme(legend.position = 'None')
gg <- plot_grid(gg_clinical, gg_severe, nrow = 1)
gg <- plot_grid(gg, gg_legend, nrow = 1, rel_widths = c(1, 0.25))
gg
f_save_plot(gg, plot_name = "PE_SMCRTSS_byEIR", width = 12, height = 10, plot_dir = file.path(plot_dir))


##-------
## PE RTSS
##---------
gg_severe <- f_custom_plot(dat, yvar = 'rtss_protective_efficacy_severe', facet_var = 'Annual_EIR') +
  labs(title = 'Protective efficacy in severe cases by EIR')

gg_clinical <- f_custom_plot(dat, yvar = 'rtss_protective_efficacy', facet_var = 'Annual_EIR') +
  labs(title = 'Protective efficacy in clinical cases by EIR')

gg_legend <- get_legend(gg_clinical)
gg_clinical <- gg_clinical + theme(legend.position = 'None')
gg_severe <- gg_severe + theme(legend.position = 'None')
gg <- plot_grid(gg_clinical, gg_severe, nrow = 1)
gg <- plot_grid(gg, gg_legend, nrow = 1, rel_widths = c(1, 0.25))
gg
f_save_plot(gg, plot_name = "PE_RTSS_byEIR", width = 12, height = 10, plot_dir = file.path(plot_dir))

##-------
## for report
##---------
gg_severe <- dat %>%
  filter(Annual_EIR == 10) %>%
  select_at(c('year', 'scenario_name', 'access_scen', 'protective_efficacy_severe', 'rtss_protective_efficacy_severe')) %>%
  pivot_longer(cols = -c('year', 'scenario_name', 'access_scen')) %>%
  as.data.frame() %>%
  f_custom_plot(yvar = 'value', facet_var = 'name') +
  labs(title = 'Severe cases')

gg_clinical <- dat %>%
  filter(Annual_EIR == 10) %>%
  select_at(c('year', 'scenario_name', 'access_scen', 'protective_efficacy', 'rtss_protective_efficacy')) %>%
  pivot_longer(cols = -c('year', 'scenario_name', 'access_scen')) %>%
  as.data.frame() %>%
  f_custom_plot(yvar = 'value', facet_var = 'name') +
  labs(title = 'Clinical cases')

gg_legend <- get_legend(gg_clinical)
gg_clinical <- gg_clinical + theme(legend.position = 'None')
gg_severe <- gg_severe + theme(legend.position = 'None')
gg <- plot_grid(gg_clinical, gg_severe, nrow = 1)
gg <- plot_grid(gg, gg_legend, nrow = 1, rel_widths = c(1, 0.25))
gg
f_save_plot(gg, plot_name = "PE_combined_EIR10", width = 12, height = 10, plot_dir = file.path(plot_dir))


##-------
## for report
##---------
gg_severe <- dat %>%
  filter(Annual_EIR == 10) %>%
  select_at(c('year', 'scenario_name', 'access_scen','severe_cases_averted_per100000', 'rtss_severe_cases_averted_per100000')) %>%
  pivot_longer(cols = -c('year', 'scenario_name', 'access_scen')) %>%
  as.data.frame() %>%
  f_custom_plot(yvar = 'value', facet_var = 'name') +
  labs(title = 'Severe cases')

gg_clinical <- dat %>%
  filter(Annual_EIR == 10) %>%
  select_at(c('year', 'scenario_name', 'access_scen','cases_averted_per100000','rtss_cases_averted_per100000')) %>%
  pivot_longer(cols = -c('year', 'scenario_name', 'access_scen')) %>%
  as.data.frame() %>%
  f_custom_plot(yvar = 'value', facet_var = 'name') +
  labs(title = 'Clinical cases')

gg_legend <- get_legend(gg_clinical)
gg_clinical <- gg_clinical + theme(legend.position = 'None')
gg_severe <- gg_severe + theme(legend.position = 'None')
gg <- plot_grid(gg_clinical, gg_severe, nrow = 1)
gg <- plot_grid(gg, gg_legend, nrow = 1, rel_widths = c(1, 0.25))

f_save_plot(gg, plot_name = "cases_averted_combined_EIR10", width = 12, height = 10, plot_dir = file.path(plot_dir))

