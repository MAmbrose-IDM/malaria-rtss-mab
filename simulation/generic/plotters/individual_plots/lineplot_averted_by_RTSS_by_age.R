# plot_cases_averted_by_age.R


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


#######################################################################################
# plot cases averted in each age for different scenarios (panel of lineplots)
#######################################################################################
exp_name = 'generic_rtss_sweep'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
cases_scenarios_each_year = sim_output[[4]]
cases_scenarios_each_year = f_add_scenario_name(df = cases_scenarios_each_year, scenario_type = 'scenario_name2')
cases_scenarios_each_year = cases_scenarios_each_year[!is.na(cases_scenarios_each_year$scenario_name),]

# subset plotted
cur_cm = 0.6
cur_corr = 0


# =======  compare different RTS,S schedules, without SMC  ====== #
cur_smc = 0
# subset simulation to current input
cur_df = filter(cases_scenarios_each_year,
                intervention_correlation == cur_corr,
                cm_coverage == cur_cm,
                smc_coverage == 0)
# create plot
age_label_values = seq(0, 8, length.out = 5)
age_labels = paste0(age_label_values, '-', (age_label_values + 1))
# clinical cases
gg = ggplot(cur_df, aes(x=year, y=cases_averted_per100000))+
  geom_point(aes(col=as.factor(rtss_scenario), group=rtss_scenario), size=2)+
  geom_line(aes(col=as.factor(rtss_scenario), group=rtss_scenario), size=1)+
  f_getColorMapping(name='RTS,S scenario', type='color') +
  theme_bw()+
  f_getCustomTheme() +
  geom_hline(yintercept=0)+
  ylab('cases averted (per 100,000) when adding RTS,S and/or SMC') + 
  xlab('age of child') +
  scale_x_continuous(breaks=age_label_values, labels=age_labels) +
  facet_wrap(seasonality~Annual_EIR)
f_save_plot(gg, paste0('cases_averted_by_age_RTSS_and_SMC_',
                       cur_cm * 100, 'CM_',
                       cur_smc * 100, 'SMC_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 9, height = 5, units = 'in', device_format = device_format)

# severe cases
gg = ggplot(cur_df, aes(x = year, y = severe_cases_averted_per100000)) +
  geom_point(aes(col = as.factor(rtss_scenario), group = rtss_scenario), size = 2) +
  geom_line(aes(col = as.factor(rtss_scenario), group = rtss_scenario), size = 1) +
  f_getColorMapping(name = 'RTS,S scenario', type = 'color') +
  theme_bw() +
  f_getCustomTheme() +
  geom_hline(yintercept = 0) +
  ylab('severe cases averted (per 100,000) when adding RTS,S and/or SMC') +
  xlab('age of child') +
  scale_x_continuous(breaks = age_label_values, labels = age_labels) +
  facet_wrap(seasonality ~ Annual_EIR)
f_save_plot(gg, paste0('severe_cases_averted_by_age_RTSS_and_SMC_',
                       cur_cm * 100, 'CM_',
                       cur_smc * 100, 'SMC_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 9, height = 5, units = 'in', device_format = device_format)

# =======  compare RTS,S-only (standard) at 80%, SMC-only at 80%, and SMC + RTSS standard at 80%  =======  #
cases_scenarios_each_year = sim_output[[4]]
cases_scenarios_each_year = f_add_scenario_name(df = cases_scenarios_each_year, scenario_type = 'scenario_name')
cases_scenarios_each_year = cases_scenarios_each_year[!is.na(cases_scenarios_each_year$scenario_name),]
# subset simulation to current input
cur_df1 = filter(cases_scenarios_each_year,
                 intervention_correlation == cur_corr,
                 cm_coverage == cur_cm,
                 smc_coverage == 0,
                 rtss_coverage == 0.8,
                 rtss_mode == 'constant')
cur_df1$scenario_name = 'RTS,S only'
cur_df2 = filter(cases_scenarios_each_year,
                 intervention_correlation == cur_corr,
                 cm_coverage == cur_cm,
                 smc_coverage == 0.8,
                 rtss_coverage == 0)
cur_df2$scenario_name = 'SMC only'
cur_df3 = filter(cases_scenarios_each_year,
                 intervention_correlation == cur_corr,
                 cm_coverage == cur_cm,
                 smc_coverage == 0.8,
                 rtss_coverage == 0.8,
                 rtss_mode == 'constant')
cur_df3$scenario_name = 'both RTS,S and SMC'
cur_df = rbind(cur_df1, cur_df2, cur_df3)
cur_df = cur_df[!is.na(cur_df$scenario_name),]
cur_df$scenario_name = factor(cur_df$scenario_name, levels = c("RTS,S only", "SMC only", "both RTS,S and SMC"))

# create plot
age_label_values = seq(0, 8, length.out = 5)
age_labels = paste0(age_label_values, '-', (age_label_values + 1))
gg = ggplot(cur_df, aes(x = year, y = cases_averted_per100000)) +
  geom_point(aes(col = as.factor(scenario_name), group = scenario_name), size = 2) +
  geom_line(aes(col = as.factor(scenario_name), group = scenario_name), size = 1) +
  f_getColorMapping(name = 'scenario', type = 'color') +
  theme_bw() +
  f_getCustomTheme() +
  geom_hline(yintercept = 0) +
  ylab('cases averted (per 100,000) when adding RTS,S and/or SMC') +
  xlab('age of child') +
  scale_x_continuous(breaks = age_label_values, labels = age_labels) +
  facet_wrap(seasonality ~ Annual_EIR)

f_save_plot(gg, paste0('cases_averted_by_age_80RTSS_and_80SMC_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)




#######################################################################################
# plot cases averted in each age for different RTSS-SMC correlation scenarios (panel of lineplots)
#######################################################################################

exp_name = 'generic_accesscorrelation_highSMC'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
cases_scenarios_each_year = sim_output[[4]]
cases_scenarios_each_year = f_add_scenario_name(df = cases_scenarios_each_year, scenario_type = 'scenario_RTSS_SMC_withCorr2')
cases_scenarios_each_year = cases_scenarios_each_year[!is.na(cases_scenarios_each_year$scenario_name),]

# subset plotted
cur_cm = 0.6
cur_cm_target = 'random'
smc_rtss_cov = 0.8
cur_season = 'high_unimodal'
cur_scenarios = c("80% RTS,S; no SMC",
                  "no RTS,S; 80% SMC",
                  "80% RTS,S; 80% SMC; independent",
                  "80% RTS,S; 80% SMC; correlated",
                  "80% RTS,S; 80% SMC; anti-correlated")
cur_scenarios = gsub(pattern='80%', replacement=paste0(round(smc_rtss_cov*100), '%'), x=cur_scenarios)

# subset simulation to current input
cur_df = filter(cases_scenarios_each_year,
                cm_coverage == cur_cm,
                seasonality == cur_season,
                scenario_name %in% cur_scenarios,
                cm_target_group == cur_cm_target
)

cur_df = cur_df[!is.na(cur_df$scenario_name),]
# cur_df$scenario_name = factor(cur_df$scenario_name, levels = c("RTS,S only", "SMC only", "both RTS,S and SMC"))

# create plot
age_label_values = seq(0, 8, length.out = 5)
age_labels = paste0(age_label_values, '-', (age_label_values + 1))
gg = ggplot(cur_df, aes(x = year, y = cases_averted_per100000)) +
  # geom_point(aes(col = as.factor(scenario_name), group = scenario_name), size = 2) +
  geom_line(aes(col = as.factor(scenario_name), group = interaction(scenario_name, rtss_target_group)), size = 1.2) +
  # geom_line(aes(col = as.factor(scenario_name), group = interaction(scenario_name, rtss_target_group), linetype=as.factor(rtss_target_group)), size = 1.2) +
  f_getColorMapping(name = 'scenario', type = 'color') +
  theme_bw() +
  f_getCustomTheme() +
  geom_hline(yintercept = 0) +
  ylab('annual cases averted per 100,000 \n by RTS,S and SMC') +
  xlab('age of child') +
  scale_x_continuous(breaks = age_label_values, labels = age_labels) +
  facet_wrap( ~ Annual_EIR)

f_save_plot(gg, paste0('cases_averted_by_age_',100*smc_rtss_cov, 'RTSSSMCcov_accessCorr_',
                       cur_cm * 100, cur_cm_target, 'CM'),
            file.path(simout_dir, '_plots'), width = 10, height = 2.85, units = 'in', device_format = device_format)


# severe cases averted
gg = ggplot(cur_df, aes(x = year, y = severe_cases_averted_per100000)) +
  # geom_point(aes(col = as.factor(scenario_name), group = scenario_name), size = 2) +
  geom_line(aes(col = as.factor(scenario_name), group = interaction(scenario_name, rtss_target_group)), size = 1.2) +
  # geom_line(aes(col = as.factor(scenario_name), group = interaction(scenario_name, rtss_target_group), linetype=as.factor(rtss_target_group)), size = 1.2) +
  f_getColorMapping(name = 'scenario', type = 'color') +
  theme_bw() +
  f_getCustomTheme() +
  geom_hline(yintercept = 0) +
  ylab('annual severe cases averted \n per 100,000 by RTS,S and SMC') +
  xlab('age of child') +
  scale_x_continuous(breaks = age_label_values, labels = age_labels) +
  facet_wrap( ~ Annual_EIR)

f_save_plot(gg, paste0('severe_averted_by_age_',100*smc_rtss_cov, 'RTSSSMCcov_accessCorr_',
                       cur_cm * 100, cur_cm_target, 'CM'),
            file.path(simout_dir, '_plots'), width = 10, height = 2.85, units = 'in', device_format = device_format)



#######################################################################################
# plot cases averted in each age for different RTSS-SMC correlation scenarios (panel of lineplots) - one EIR, two coverages, cases and severe
#######################################################################################

exp_name = 'generic_accesscorrelation_highSMC'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
cases_scenarios_each_year = sim_output[[4]]
cases_scenarios_each_year = f_add_scenario_name(df = cases_scenarios_each_year, scenario_type = 'scenario_RTSS_SMC_withCorr3')
cases_scenarios_each_year = cases_scenarios_each_year[!is.na(cases_scenarios_each_year$scenario_name),]

# remove scenarios that are not plotted (multiple RTSS targeting when there isn't SMC; no RTSS and no SMC)
cases_scenarios_each_year = cases_scenarios_each_year[-intersect(which(cases_scenarios_each_year$smc_coverage<0.001), which(cases_scenarios_each_year$rtss_target_group !='random')),]
cases_scenarios_each_year = cases_scenarios_each_year[-intersect(which(cases_scenarios_each_year$smc_coverage<0.001), which(cases_scenarios_each_year$rtss_coverage<0.001)),]
cases_scenarios_each_year = cases_scenarios_each_year[!is.na(cases_scenarios_each_year$scenario_name),]

age_label_values = seq(0, 8, length.out = 5)
age_labels = paste0(age_label_values, '-', (age_label_values + 1))

# specify plotted scenarios
cur_eir = 10
cur_cm = 0.6
cur_cm_target = 'random'
smc_rtss_cov = c(0.6, 0.8)
cur_season = 'high_unimodal'
cur_scenarios = c("no SMC",
                  "no RTS,S",
                  "independent",
                  "correlated",
                  "anti-correlated")
cur_df = filter(cases_scenarios_each_year,
                Annual_EIR == cur_eir,
                cm_coverage == cur_cm,
                seasonality == cur_season,
                # scenario_name %in% cur_scenarios,
                cm_target_group == cur_cm_target

)

# add tag for coverage subplotting (60% coverage versus 80% coverage)
cur_df$cov_plot = NA
cur_df$cov_plot[which(cur_df$rtss_coverage == 0.6 | cur_df$smc_coverage == 0.6)] = '60% coverage'
cur_df$cov_plot[which(cur_df$rtss_coverage == 0.8 | cur_df$smc_coverage == 0.8)] = '80% coverage'




gg_case = ggplot(cur_df, aes(x = year, y = cases_averted_per100000)) +
  # geom_point(aes(col = as.factor(scenario_name), group = scenario_name), size = 2) +
  geom_line(aes(col = as.factor(scenario_name), group = interaction(scenario_name, rtss_target_group)), size = 0.8) +
  # geom_line(aes(col = as.factor(scenario_name), group = interaction(scenario_name, rtss_target_group), linetype=as.factor(rtss_target_group)), size = 1.2) +
  f_getColorMapping(name = 'scenario', type = 'color') +
  theme_bw() +
  f_getCustomTheme() +
  geom_hline(yintercept = 0) +
  ylab('annual cases averted per 100,000 \n by RTS,S and SMC') +
  xlab('age of child') +
  scale_x_continuous(breaks = age_label_values, labels = age_labels) +
  facet_wrap( ~ cov_plot)


# severe cases averted
gg_severe = ggplot(cur_df, aes(x = year, y = severe_cases_averted_per100000)) +
  # geom_point(aes(col = as.factor(scenario_name), group = scenario_name), size = 2) +
  geom_line(aes(col = as.factor(scenario_name), group = interaction(scenario_name, rtss_target_group)), size = 0.8) +
  # geom_line(aes(col = as.factor(scenario_name), group = interaction(scenario_name, rtss_target_group), linetype=as.factor(rtss_target_group)), size = 1.2) +
  f_getColorMapping(name = 'scenario', type = 'color') +
  theme_bw() +
  f_getCustomTheme() +
  geom_hline(yintercept = 0) +
  ylab('annual severe cases averted \n per 100,000 by RTS,S and SMC') +
  xlab('age of child') +
  scale_x_continuous(breaks = age_label_values, labels = age_labels) +
  facet_wrap( ~ cov_plot)


# combine plots into panel
gg_legend <- get_legend(gg_case)
gg_case <- gg_case  + theme(legend.position = 'None')
gg_severe <- gg_severe  + theme(legend.position = 'None')
gg <- plot_grid(gg_case,gg_severe, align=c('hv'),nrow=2)
gg <- plot_grid(gg, gg_legend, rel_widths = c(1, 0.5))



f_save_plot(gg, paste0('clinical_severe_averted_by_age_accessCorr_',cur_eir, 'EIR_',
                       cur_cm * 100, cur_cm_target, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 6, units = 'in', device_format = device_format)





