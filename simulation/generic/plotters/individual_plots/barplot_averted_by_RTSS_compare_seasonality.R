# barplot_averted_by_RTSS_compare_seasonality.R


library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)


#######################################################
# setup
#######################################################
# wdir at rtss-scenarios
USER <- Sys.getenv('USERNAME')
setwd(paste0('C:/Users/', USER, '/Documents/rtss-scenarios'))

source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path('simulation', 'generic', 'plotters', 'plot_bars_and_lines.R'))
source(file.path("simulation", "load_paths.R"))

simout_dir <- file.path(projectpath, 'simulation_output', 'generic')
if (!dir.exists(paste0(simout_dir, '/_plots'))) dir.create(paste0(simout_dir, '/_plots'))
if (!dir.exists(paste0(simout_dir, '/_plots/pdf'))) dir.create(paste0(simout_dir, '/_plots/pdf'))

exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'cm_coverage', 'rtss_coverage', 'rtss_mode', 'smc_coverage',
                'intervention_correlation', 'Cohort_birth_month', 'minBoostAge')
max_years = c(1, 2, 5, 10)
age_groups = paste0('U', max_years)
seasonality_levels = c('constant', 'moderate', 'high')
rtss_mode_levels = c('constant', 'campboost', 'campboost2')
rtss_mode_levels_plot = c('standard', 'enhanced (1b)', 'enhanced (2b)')
rtss_mode_levels_lookup = rtss_mode_levels_plot
names(rtss_mode_levels_lookup) = rtss_mode_levels

theme_set(theme_bw())
customTheme <- f_getCustomTheme()
device_format <- c('pdf', 'png')



##########################################################################################################
# barplots of cases averted at different seasonality patterns (for two SMC coverages)
##########################################################################################################
exp_name = 'generic_season_sweep_EIR_10_30'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]


# ====== percent reduction (facted by RTS,S type and smc coverage) ====== #
cur_cm = 0.6
cur_corr = 0
cur_rtss = 0.8
cur_eir = 10
cur_age = 'U5'
cur_df = filter(pe_df,
                age_group == cur_age,
                intervention_correlation == cur_corr,
                cm_coverage == cur_cm,
                rtss_coverage == cur_rtss,
                Annual_EIR == cur_eir,
)
cur_df$seasonality = gsub('_unimodal', '', cur_df$seasonality)
cur_df$seasonality = factor(cur_df$seasonality, levels = seasonality_levels)
cur_df$rtss_percent_reduction = cur_df$rtss_relative_burden * -100
cur_df$rtss_percent_reduction_severe = cur_df$rtss_relative_burden_severe * -100
cur_df$rtss_mode = as.character(rtss_mode_levels_lookup[cur_df$rtss_mode])
cur_df$rtss_mode = factor(cur_df$rtss_mode, levels = rtss_mode_levels_plot)


# percent reduction - clinical cases
gg = plot_barplots(cur_df, xvar = 'seasonality', yvar = 'rtss_percent_reduction',
                   fillvar = 'seasonality', facet1 = 'rtss_mode', facet2 = 'smc_coverage',
                   SAVE = FALSE, ylab = paste0('percent reduction in ', cur_age, ' clinical cases relative to no-RTS,S'),
                   scales = 'fixed')
gg = gg + f_getColorMapping(name = 'seasonality', type = 'fill')
f_save_plot(gg, paste0('percent_reduction_cases_RTSS_seasonality_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_rtss * 100, 'RTSS_',
                       cur_corr, 'corr_',
                       cur_eir, 'EIR'),
            file.path(simout_dir, '_plots'), width = 6, height = 7, units = 'in', device_format = device_format)


# percent reduction - severe cases
gg = plot_barplots(cur_df, xvar = 'seasonality', yvar = 'rtss_percent_reduction_severe',
                   fillvar = 'seasonality', facet1 = 'rtss_mode', facet2 = 'smc_coverage',
                   SAVE = FALSE, ylab = paste0('percent reduction in ', cur_age, ' severe cases relative to no-RTS,S'),
                   scales = 'fixed')
gg = gg + f_getColorMapping(name = 'seasonality', type = 'fill')
f_save_plot(gg, paste0('percent_reduction_severe_cases_RTSS_seasonality_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_rtss * 100, 'RTSS_',
                       cur_corr, 'corr_',
                       cur_eir, 'EIR'),
            file.path(simout_dir, '_plots'), width = 6, height = 7, units = 'in', device_format = device_format)

# ====== percent reduction (facted by severe versus clinical and by SMC) ====== #
cur_cm = 0.6
cur_corr = 0
cur_rtss = 0.8
cur_eir = 10
cur_age = 'U5'
cur_mode = 'constant'
cur_df = filter(pe_df,
                age_group == cur_age,
                intervention_correlation == cur_corr,
                cm_coverage == cur_cm,
                rtss_coverage == cur_rtss,
                Annual_EIR == cur_eir,
                rtss_mode == cur_mode
)
cur_df$seasonality = gsub('_unimodal', '', cur_df$seasonality)
cur_df$seasonality = factor(cur_df$seasonality, levels = seasonality_levels)
cur_df$rtss_percent_reduction = cur_df$rtss_relative_burden * -100
cur_df$rtss_percent_reduction_severe = cur_df$rtss_relative_burden_severe * -100

# percent reduction in burden (relative to no-RTSS scenario)
gg1 = plot_barplots(cur_df, xvar = 'seasonality', yvar = 'rtss_percent_reduction',
                    fillvar = 'seasonality', facet1 = 'smc_coverage', facet2 = NULL,
                    SAVE = FALSE, ylab = paste0('percent reduction in ', cur_age, ' clinical \n cases relative to no-RTS,S'),
                    scales = 'fixed')
gg1 = gg1 + f_getColorMapping(name = 'seasonality', type = 'fill')
gg1 = gg1 + theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))

gg2 = plot_barplots(cur_df, xvar = 'seasonality', yvar = 'rtss_percent_reduction_severe',
                    fillvar = 'seasonality', facet1 = 'smc_coverage', facet2 = NULL,
                    SAVE = FALSE, ylab = paste0('percent reduction in ', cur_age, ' severe \n cases relative to no-RTS,S'),
                    scales = 'fixed')
gg2 = gg2 + f_getColorMapping(name = 'seasonality', type = 'fill')
gg2 = gg2 + theme(plot.title = element_blank(), plot.subtitle = element_blank(),
                  plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"))

tt = plot_grid(gg1, gg2, align = "v", nrow = 2, rel_heights = c(1, 1.325))
f_save_plot(tt, paste0('percent_reduction_RTSS_seasonality_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_rtss * 100, 'RTSS_',
                       cur_corr, 'corr_',
                       cur_eir, 'EIR'),
            file.path(simout_dir, '_plots'), width = 8, height = 7, units = 'in', device_format = device_format)
