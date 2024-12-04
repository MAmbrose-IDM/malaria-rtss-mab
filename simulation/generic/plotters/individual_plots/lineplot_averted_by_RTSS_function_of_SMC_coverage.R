# plot_averted_by_RTSS_function_of_SMC_coverage.R


library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)

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




#####################################################################################################
# impact of RTS,S as a function of SMC coverage (lineplot with lines for different RTS,S schedules)
#####################################################################################################
exp_name = 'generic_rtss_sweep'  # <- only has 2 SMC coverages (0, 0.8). Or can use generic_factorial_v1, which has SMC=0, 0.6, 0.8, but only standard RTS,S
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'rtss_mode'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))

# subset simulations
cur_cm = 0.6
cur_corr = 0
cur_age = 'U5'
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   intervention_correlation == cur_corr,
                   cm_coverage == cur_cm,
                   rtss_coverage > 0
)

# clinical cases averted
gg = ggplot(pe_df_cur, aes(x = smc_coverage, y = rtss_cases_averted_per100000)) +
  geom_line(aes(col = as.factor(rtss_scenario)), size = 2) +
  ylab('cases averted (per 100,000) when adding RTS,S') +
  xlab('SMC coverage') +
  theme_bw() +
  f_getCustomTheme() +
  f_getColorMapping(name = 'RTS,S scenario', type = 'color') +
  facet_wrap(seasonality ~ Annual_EIR)
f_save_plot(gg, paste0('cases_averted_RTSS_as_function_of_SMC_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 9, height = 5, units = 'in', device_format = device_format)


# severe cases averted
gg = ggplot(pe_df_cur, aes(x = smc_coverage, y = rtss_severe_cases_averted_per100000)) +
  geom_line(aes(col = as.factor(rtss_scenario)), size = 2) +
  ylab('severe cases averted (per 100000) when adding RTS,S') +
  xlab('SMC coverage') +
  theme_bw() +
  f_getCustomTheme() +
  f_getColorMapping(name = 'RTS,S scenario', type = 'color') +
  facet_wrap(seasonality ~ Annual_EIR)

f_save_plot(gg, paste0('severe_cases_averted_RTSS_as_function_of_SMC_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 9, height = 5, units = 'in', device_format = device_format)

