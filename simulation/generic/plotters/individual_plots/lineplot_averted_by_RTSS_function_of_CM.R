# lineplot_averted_by_RTSS_function_of_CM.R


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
# lineplot of severe cases averted at different CM rates (for different age groups)
##########################################################################################################
exp_name = 'generic_cm_sweep_SMC'  #  'generic_cm_sweep_SMC' # 'generic_cm_correlated_sweep'  # 'generic_cm_sweep_SMC'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'scenario_name')
pe_df = pe_df[!is.na(pe_df$scenario_name),]

# subset to plotted scenario
cur_corr = 0
cur_rtss = 0.8
cur_smc = 0
cur_season = 'moderate_unimodal'
cur_eir = 10
cur_df = filter(pe_df,
                Annual_EIR == cur_eir,
                intervention_correlation == cur_corr,
                smc_coverage == cur_smc,
                rtss_coverage == cur_rtss,
                seasonality == cur_season
)
cur_df = cur_df[cur_df$age_group != 'U1',]
cur_df$age_group = factor(cur_df$age_group, levels = age_groups)

# NOTE: when this is plotted for clinical cases averted it will likely have an unintuitive pattern due to how EMOD counts 'new clinical cases'
#   --> higher treatment/cure rates leads to more opportunities for new infections

# severe cases
gg = plot_lines(dat = cur_df, xvar = 'cm_coverage', yvar = 'rtss_severe_cases_averted_per100000',
                colvar = 'age_group', facet1 = NULL, facet2 = NULL,
                SAVE = FALSE, xlab = 'effective case management rate', ylab = paste0('severe cases averted \n by adding RTS,S per 100,000'))
gg = gg + f_getColorMapping(name = 'age group', type = 'color')

f_save_plot(gg, paste0('severe_cases_averted_RTSS_CM_',
                       cur_season, '_',
                       cur_eir, 'EIR_',
                       cur_smc * 100, 'SMC_',
                       cur_rtss * 100, 'RTSS_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 6, height = 5 * 6 / 8, units = 'in', device_format = device_format)

