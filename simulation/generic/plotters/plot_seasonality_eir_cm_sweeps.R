# plot_seasonality_eir_cm_sweeps.R


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


#####################################################################################################
# grouped barplots showing cases averted when adding various coverages/schedules of RTS,S against different backgrounds (seasonality, EIR, and SMC)
#####################################################################################################
exp_name = 'generic_rtss_age_sweep'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))


# ===  show bars with all RTS,S schedules and coverages  === #

# subset simulations
cur_cm = 0.6
cur_corr = 0
cur_age = 'U5'
cur_smcs = c(0, 0.8)
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   intervention_correlation == cur_corr,
                   cm_coverage == cur_cm,
                   rtss_coverage > 0,
                   smc_coverage %in% cur_smcs,
                   # minBoostAge>20,
)

gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_cases_averted_per100000', bargroup_var = 'rtss_scenario',
                           fillvar = 'rtss_scenario', facet1 = 'seasonality', facet2 = 'smc_coverage',
                           SAVE = FALSE, ylab = 'cases averted by RTS,S per 100,000')
gg = gg + f_getColorMapping(name = 'RTS,S scenario', type = 'fill')

f_save_plot(gg, paste0('cases_averted_by_RTSS_schedules_for_smc_eir_seasonality_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            # '_minAge24m'),
            file.path(simout_dir, '_plots'), width = 8.5, height = 3.8, units = 'in', device_format = device_format)

gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_severe_cases_averted_per100000', bargroup_var = 'rtss_scenario',
                           fillvar = 'rtss_scenario', facet1 = 'seasonality', facet2 = 'smc_coverage',
                           SAVE = FALSE, ylab = 'severe cases averted by RTS,S per 100,000')
gg = gg + f_getColorMapping(name = 'RTS,S scenario', type = 'fill')


f_save_plot(gg, paste0('severe_cases_averted_by_RTSS_schedules_for_smc_eir_seasonality_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            # '_minAge24m'),
            file.path(simout_dir, '_plots'), width = 8.5, height = 3.8, units = 'in', device_format = device_format)


######################################
add_flexible_booster = F
######################################
if (add_flexible_booster) {

  #####################################################################################################
  # grouped barplots showing cases averted when adding various schedules of RTS,S
  #####################################################################################################
  exp_names = c('generic_rtss_age_sweep', 'generic_campboost3_SMC_ext')
  exp_filepath = file.path(simout_dir, exp_names[1])
  sim_output1 = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_names[1], exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                       add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
  pe_df1 = sim_output1[[3]]
  pe_df1 = f_add_scenario_name(df = pe_df1, scenario_type = 'rtss_scenario')
  pe_df1 = pe_df1[!is.na(pe_df1$rtss_scenario),]
  pe_df1$Annual_EIR = factor(pe_df1$Annual_EIR, levels = sort(unique(pe_df1$Annual_EIR)))

  simout_dir2 <- gsub('generic', 'generic_forward', simout_dir)
  exp_filepath = file.path(simout_dir2, exp_names[2])
  sim_output2 = load_Age_monthly_Cases(simout_dir = simout_dir2, exp_name = exp_names[2], exp_sweeps = c(exp_sweeps, 'rtss_mode'),
                                       add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
  pe_df2 = sim_output2[[3]]
  pe_df2$rtss_scenario = factor(pe_df2$rtss_mode,
                                levels = c('constant', 'campboost3', 'campboost3B'),
                                labels = c('1 booster all at age 24-35 mths',
                                           '1 booster spread out between age 24-59 mths',
                                           '1 booster spread out between age 24-59 mths\n+2nd booster between 48-59 mths'))

  pe_df2 = pe_df2[!is.na(pe_df2$rtss_scenario),]
  pe_df2$Annual_EIR = factor(pe_df2$Annual_EIR, levels = sort(unique(pe_df2$Annual_EIR)))

  (drop_cols2 <- colnames(pe_df2)[!(colnames(pe_df2) %in% colnames(pe_df1))])
  (drop_cols1 <- colnames(pe_df1)[!(colnames(pe_df1) %in% colnames(pe_df2))])
  pe_df1 <- pe_df1 %>% dplyr::select(-drop_cols1)
  pe_df2 <- pe_df2 %>% dplyr::select(-drop_cols2)
  table(pe_df2$smc_coverage, pe_df2$rtss_scenario)
  table(pe_df2$Annual_EIR, pe_df2$rtss_scenario)

  pe_df2 <- pe_df2 %>% dplyr::select(colnames(pe_df1))
  pe_df <- rbind(pe_df1, pe_df2)

  # ===  show bars with all RTS,S schedules and coverages  === #

  # subset simulations
  cur_cm = 0.6
  cur_corr = 0
  cur_age = 'U5'
  cur_smcs = c(0, 0.8)
  cur_seasonality = 'high_unimodal'
  pe_df_cur = pe_df %>%
    filter(age_group == cur_age,
           seasonality == cur_seasonality,
           intervention_correlation == cur_corr,
           cm_coverage == cur_cm,
           rtss_coverage > 0,
           smc_coverage %in% cur_smcs,
           # minBoostAge>20,
    )


  gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_cases_averted_per100000', bargroup_var = 'rtss_scenario',
                             fillvar = 'rtss_scenario', facet1 = 'seasonality', facet2 = 'smc_coverage',
                             SAVE = FALSE, ylab = 'cases averted by RTS,S per 100,000')

  f_save_plot(gg, paste0('cases_averted_by_RTSS_booster_scen_',
                         cur_age, '_',
                         cur_cm * 100, 'CM_',
                         cur_corr, 'corr'),
              file.path(simout_dir2, '_plots'), width = 8.5, height = 3.8, units = 'in', device_format = device_format)

  gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_severe_cases_averted_per100000', bargroup_var = 'rtss_scenario',
                             fillvar = 'rtss_scenario', facet1 = 'seasonality', facet2 = 'smc_coverage',
                             SAVE = FALSE, ylab = 'severe cases averted by RTS,S per 100,000')

  f_save_plot(gg, paste0('severe_cases_averted_by_RTSS_booster_scen_',
                         cur_age, '_',
                         cur_cm * 100, 'CM_',
                         cur_corr, 'corr'),
              # '_minAge24m'),
              file.path(simout_dir2, '_plots'), width = 8.5, height = 3.8, units = 'in', device_format = device_format)
}

# ===  simple with versus without SMC  === #
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'scenario_name2')
# subset simulations
cur_cm = 0.6
cur_corr = 0
cur_age = 'U5'
cur_smcs = c(0, 0.8)
cur_rtss_scenarios = c('80% enhanced RTS,S (1 booster, min 24m)')
cur_season = 'high_unimodal'
pe_df_cur = filter(pe_df,
                   seasonality == cur_season,
                   age_group == cur_age,
                   intervention_correlation == cur_corr,
                   cm_coverage == cur_cm,
                   rtss_coverage > 0,
                   smc_coverage %in% cur_smcs,
                   rtss_scenario %in% cur_rtss_scenarios
)
pe_df_cur$plotted_scen = as.character(pe_df_cur$scenario_name)
pe_df_cur$plotted_scen[pe_df_cur$scenario_name == '80% enhanced RTS,S (1 booster, min 24m); 0% SMC'] = 'without SMC'
pe_df_cur$plotted_scen[pe_df_cur$scenario_name == '80% enhanced RTS,S (1 booster, min 24m); 80% SMC'] = 'with SMC'
pe_df_cur$plotted_scen = factor(pe_df_cur$plotted_scen, levels = c('without SMC', 'with SMC'))

gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_cases_averted_per100000', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = 'seasonality', facet2 = NULL,
                           SAVE = FALSE, ylab = 'cases averted by\nRTS,S per 100,000')
gg = gg + f_getColorMapping(name = 'scenario', type = 'fill')

f_save_plot(gg, paste0('simple_cases_averted_by_80RTSS_with_without_80SMC_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 4.7, height = 2.8, units = 'in', device_format = device_format)
gg1 <- gg

gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_severe_cases_averted_per100000', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = 'seasonality', facet2 = NULL,
                           SAVE = FALSE, ylab = 'severe cases averted by\nRTS,S per 100,000')
gg = gg + f_getColorMapping(name = 'scenario', type = 'fill')

f_save_plot(gg, paste0('simple_severe_cases_averted_by_80RTSS_with_without_80SMC_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 4.4, height = 2.8, units = 'in', device_format = device_format)
gg2 <- gg

# ==== percent reduction for simple comparison ==== #
pe_df_cur$rtss_percent_reduction = pe_df_cur$rtss_relative_burden * -100
pe_df_cur$rtss_percent_reduction_severe = pe_df_cur$rtss_relative_burden_severe * -100
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_percent_reduction', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = 'seasonality', facet2 = NULL,  #'smc_coverage',
                           SAVE = FALSE, ylab = 'percent reduction in cases U5\nfrom adding RTS,S')
gg = gg + f_getColorMapping(name = 'scenario', type = 'fill')

f_save_plot(gg, paste0('simple_percent_reduction_clinical_RTSS_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr',
                       '_minAge24m'),
            file.path(simout_dir, '_plots'), width = 4.4, height = 2.8, units = 'in', device_format = device_format)
gg3 <- gg

pe_df_cur$rtss_percent_reduction = pe_df_cur$rtss_relative_burden * -100
pe_df_cur$rtss_percent_reduction_severe = pe_df_cur$rtss_relative_burden_severe * -100
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_percent_reduction_severe', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = 'seasonality', facet2 = NULL, # 'smc_coverage',
                           SAVE = FALSE, ylab = 'percent reduction in severe cases U5\nfrom adding RTS,S')
gg = gg + f_getColorMapping(name = 'scenario', type = 'fill')


f_save_plot(gg, paste0('simple_percent_reduction_severe_RTSS_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr',
                       '_minAge24m'),
            file.path(simout_dir, '_plots'), width = 4.4, height = 2.8, units = 'in', device_format = device_format)
gg4 <- gg

gg_legend <- get_legend(gg1)
gg1 <- gg1  + theme(legend.position = 'None')
gg2 <- gg2  + theme(legend.position = 'None')
gg3 <- gg3  + theme(legend.position = 'None')
gg4 <- gg4  + theme(legend.position = 'None')
gg <- plot_grid(gg1,gg2,gg3,gg4)
gg <- plot_grid(gg, gg_legend, rel_widths = c(1, 0.25))
f_save_plot(gg, paste0('simple_percent_reduction_severe_RTSS_combined', '_minAge24m'),
            file.path(simout_dir, '_plots'), width = 7.5, height =5, units = 'in', device_format = device_format)

rm(gg, gg1, gg2, gg3,gg4)
# ===  two RTS,S modes with versus without SMC  === #
exp_name = 'generic_rtss_sweep'  # 'generic_rtss_age_sweep'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'scenario_name2')
# subset simulations
cur_cm = 0.6
cur_corr = 0
cur_age = 'U5'
cur_smcs = c(0, 0.8)
cur_rtss_scenarios = c('no RTS,S', '80% standard RTS,S (1 booster)', '80% enhanced RTS,S (1 booster)')
cur_season = 'high_unimodal'
pe_df_cur = filter(pe_df,
                   seasonality == cur_season,
                   age_group == cur_age,
                   intervention_correlation == cur_corr,
                   cm_coverage == cur_cm,
                   rtss_coverage > 0,
                   smc_coverage %in% cur_smcs,
                   rtss_scenario %in% cur_rtss_scenarios
)
pe_df_cur$plotted_scen = as.character(pe_df_cur$scenario_name)
pe_df_cur$plotted_scen[pe_df_cur$scenario_name == '80% standard RTS,S (1 booster); 0% SMC'] = 'standard RTS,S - without SMC'
pe_df_cur$plotted_scen[pe_df_cur$scenario_name == '80% standard RTS,S (1 booster); 80% SMC'] = 'standard RTS,S - with SMC'
pe_df_cur$plotted_scen[pe_df_cur$scenario_name == '80% enhanced RTS,S (1 booster); 0% SMC'] = 'enhanced RTS,S - without SMC'
pe_df_cur$plotted_scen[pe_df_cur$scenario_name == '80% enhanced RTS,S (1 booster); 80% SMC'] = 'enhanced RTS,S - with SMC'
pe_df_cur$plotted_scen = factor(pe_df_cur$plotted_scen, levels = c('standard RTS,S - without SMC', 'standard RTS,S - with SMC', 'enhanced RTS,S - without SMC', 'enhanced RTS,S - with SMC'))

gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_cases_averted_per100000', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = 'seasonality', facet2 = NULL,
                           SAVE = FALSE, ylab = 'cases averted by \n RTS,S per 100,000')
gg = gg + f_getColorMapping(name = 'scenario', type = 'fill')

f_save_plot(gg, paste0('simple_cases_averted_by_80RTSS_stand_enh_with_without_80SMC_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 4.7, height = 2.8, units = 'in', device_format = device_format)


gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_severe_cases_averted_per100000', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = 'seasonality', facet2 = NULL,
                           SAVE = FALSE, ylab = 'severe cases averted by RTS,S per 100,000')
gg = gg + f_getColorMapping(name = 'scenario', type = 'fill')

f_save_plot(gg, paste0('simple_severe_cases_averted_by_80RTSS_stand_enh_with_without_80SMC_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 4.4, height = 2.8, units = 'in', device_format = device_format)

##########################################################################################################
# barplots of cases averted at different seasonality patterns (for two SMC coverages)
##########################################################################################################
exp_name = 'generic_season_sweep_EIR_10_30'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]


# ====== cases averted (facted by EIR and smc coverage) ====== #
# subset to plotted scenario
cur_cm = 0.6
cur_corr = 0
cur_rtss = 0.8
cur_age = 'U5'
cur_mode = 'campboost'
cur_df = filter(pe_df,
                age_group == cur_age,
                intervention_correlation == cur_corr,
                cm_coverage == cur_cm,
                rtss_coverage == cur_rtss,
                rtss_mode == cur_mode
)
cur_df$seasonality = gsub('_unimodal', '', cur_df$seasonality)
cur_df$seasonality = factor(cur_df$seasonality, levels = seasonality_levels)

gg = plot_barplots(cur_df, xvar = 'seasonality', yvar = 'rtss_cases_averted_per100000',
                   fillvar = 'seasonality', facet1 = 'Annual_EIR', facet2 = 'smc_coverage',
                   SAVE = FALSE, ylab = paste0('cases averted by adding RTS,S (per 100,000 ', cur_age, ')'))
gg = gg + f_getColorMapping(name = 'seasonality', type = 'fill')

f_save_plot(gg, paste0('cases_averted_RTSS_seasonality_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_rtss * 100, 'RTSS_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 8, height = 6, units = 'in', device_format = device_format)

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
gg1 = gg + theme(axis.title.x = element_blank(),
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

f_save_plot(gg, paste0('percent_reduction_RTSS_seasonality_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_rtss * 100, 'RTSS_',
                       cur_corr, 'corr_',
                       cur_eir, 'EIR'),
            file.path(simout_dir, '_plots'), width = 8, height = 7, units = 'in', device_format = device_format)


##########################################################################################################
# lineplot of cases averted at different EIRs (for different age groups)
##########################################################################################################
exp_name = 'generic_eir_sweep_SMC'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'scenario_name')
pe_df = pe_df[!is.na(pe_df$scenario_name),]

# subset to plotted scenario
cur_cm = 0.6
cur_corr = 0
cur_rtss = 0.8
cur_smc = 0
cur_season = 'moderate_unimodal'
cur_df = filter(pe_df,
                intervention_correlation == cur_corr,
                cm_coverage == cur_cm,
                smc_coverage == cur_smc,
                rtss_coverage == cur_rtss,
                seasonality == cur_season
)
cur_df = cur_df[cur_df$age_group != 'U1',]
cur_df$age_group = factor(cur_df$age_group, levels = age_groups)

# clinical cases averted
gg = plot_lines(dat = cur_df, xvar = 'Annual_EIR', yvar = 'rtss_cases_averted_per100000',
                colvar = 'age_group', facet1 = NULL, facet2 = NULL,
                SAVE = FALSE, xlab = 'annual EIR', ylab = paste0('new clinical cases averted \n by adding RTS,S per 100,000'))
gg = gg + f_getColorMapping(name = 'age group', type = 'color')
gg1 = gg + theme(legend.position = 'none')

# severe cases averted
gg = plot_lines(dat = cur_df, xvar = 'Annual_EIR', yvar = 'rtss_severe_cases_averted_per100000',
                colvar = 'age_group', facet1 = NULL, facet2 = NULL,
                SAVE = FALSE, xlab = 'annual EIR', ylab = paste0('severe cases averted \n by adding RTS,S per 100,000'))
gg = gg + f_getColorMapping(name = 'age group', type = 'color')
gg2 = gg


# save figure with both clinical and severe cases plots
tt = plot_grid(gg1, gg2, align = "h", nrow = 1, rel_widths = c(1, 1.19))

f_save_plot(tt, paste0('case_and_severe_averted_RTSS_EIR_',
                       cur_season, '_',
                       cur_cm * 100, 'CM_',
                       cur_smc * 100, 'SMC_',
                       cur_rtss * 100, 'RTSS_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 16 * 2 / 3, height = 5 * 2 / 3, units = 'in', device_format = device_format)


##########################################################################################################
# lineplot of severe cases averted at different CM rates (for different age groups)
##########################################################################################################
exp_name = 'generic_cm_correlated_sweep'  #  'generic_cm_sweep_SMC' # 'generic_cm_correlated_sweep'  # 'generic_cm_sweep_SMC'
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
cur_df = filter(pe_df,
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
                       cur_cm * 100, 'CM_',
                       cur_smc * 100, 'SMC_',
                       cur_rtss * 100, 'RTSS_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 6, height = 5 * 6 / 8, units = 'in', device_format = device_format)


####################################################################################################
# line plot of total cases in setting and cases averted when RTS,S added for different scenarios
##########################################################################################################
cur_age = 'U5'
cur_rtss_scenario = '80% standard RTS,S (1 booster)'
cur_cm = 0.6
cur_corr = 0

# === with different seasonalities ===#
exp_name = 'generic_eir_season_sweep'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))
pe_df$seasonality = gsub('_unimodal', '', pe_df$seasonality)
pe_df$seasonality = factor(pe_df$seasonality, levels = seasonality_levels)

# subset simulations
cur_smc = 0
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   rtss_coverage > 0.001,
                   rtss_scenario == cur_rtss_scenario,
                   smc_coverage == cur_smc,
                   cm_coverage == cur_cm,
                   intervention_correlation == cur_corr
)

tt = plot_2x2grid_lineplots(dat = pe_df_cur, xvars = c('clinical_cases_no_rtss', 'pfpr_2_10_no_rtss'),
                            yvars = c('rtss_cases_averted_per100000', 'rtss_severe_cases_averted_per100000'),
                            colvar = 'seasonality')

f_save_plot(tt, paste0('case_and_severe_averted_RTSS_seasonality_EIR_',
                       cur_cm * 100, 'CM_',
                       cur_smc * 100, 'SMC_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 9, height = 6, units = 'in', device_format = device_format)


# percent reduction patterns
pe_df_cur$rtss_percent_reduction = pe_df_cur$rtss_relative_burden * -100
pe_df_cur$rtss_percent_reduction_severe = pe_df_cur$rtss_relative_burden_severe * -100
tt = plot_2x2grid_lineplots(dat = pe_df_cur, xvars = c('clinical_cases_no_rtss', 'pfpr_2_10_no_rtss'),
                            yvars = c('rtss_percent_reduction', 'rtss_percent_reduction_severe'),
                            colvar = 'seasonality')


# === across EIRs and SMCs === #
exp_name = 'generic_eir_smc_sweep2_newRDT'  # 'generic_eir_smc_sweep'  # 'generic_eir_smc_sweep2_newRDT'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]

# subset to plotted scenario
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   rtss_coverage > 0.001,
                   rtss_scenario == cur_rtss_scenario,
                   cm_coverage == cur_cm,
                   intervention_correlation == cur_corr
)

tt = plot_2x2grid_lineplots(dat = pe_df_cur, xvars = c('clinical_cases_no_rtss', 'pfpr_2_10_no_rtss'),
                            yvars = c('rtss_cases_averted_per100000', 'rtss_severe_cases_averted_per100000'),
                            colvar = 'smc_coverage', crop_to_second_largest = TRUE)

f_save_plot(tt, paste0('case_and_severe_averted_RTSS_SMC_EIR_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 9, height = 6, units = 'in', device_format = device_format)


# gg1 = ggplot(pe_df_cur, aes(x=clinical_cases_no_rtss, y=rtss_cases_averted_per100000)) +
#   geom_point(aes(col=Annual_EIR, shape=as.factor(smc_coverage))) +
#   geom_line(aes(col=Annual_EIR, group=as.factor(smc_coverage)))+
#   ylab(paste0('annual cases ', cur_age, ' averted \n by RTS,S per 100,000')) + 
#   xlab('annual incidence without RTS,S (per 1,000)') +
#   ylim(0, NA)+
#   theme_bw()+
#   theme(legend.position='none') + 
#   f_getCustomTheme() 
# gg2 = ggplot(pe_df_cur, aes(x=clinical_cases_no_rtss, y=rtss_severe_cases_averted_per100000)) +
#   geom_point(aes(col=Annual_EIR, shape=as.factor(smc_coverage))) +
#   geom_line(aes(col=Annual_EIR, group=as.factor(smc_coverage)))+
#   ylab(paste0('annual severe cases ', cur_age, ' averted \n by RTS,S per 100,000')) + 
#   xlab('annual incidence without RTS,S (per 1,000)') +
#   ylim(0, NA)+
#   theme_bw()+
#   theme(legend.position='none') + 
#   f_getCustomTheme() 
# gg3 = ggplot(pe_df_cur, aes(x=pfpr_2_10_no_rtss, y=rtss_cases_averted_per100000)) +
#   geom_point(aes(col=Annual_EIR, shape=as.factor(smc_coverage))) +  
#   geom_line(aes(col=Annual_EIR, group=as.factor(smc_coverage))) +  
#   ylab(paste0('annual cases ', cur_age, ' averted \n by RTS,S per 100,000')) + 
#   xlab(expression(paste(italic(Pf), PR[2-10]~without~RTS, ',', S))) +
#   ylim(0, NA)+
#   theme_bw()+
#   theme(legend.position='none') + 
#   f_getCustomTheme() 
# gg4 = ggplot(pe_df_cur, aes(x=pfpr_2_10_no_rtss, y=rtss_severe_cases_averted_per100000)) +
#   geom_point(aes(col=Annual_EIR, shape=as.factor(smc_coverage))) +  
#   geom_line(aes(col=Annual_EIR, group=as.factor(smc_coverage))) +  
#   ylab(paste0('annual severe cases ', cur_age, ' averted \n by RTS,S per 100,000')) + 
#   xlab(expression(paste(italic(Pf), PR[2-10]~without~RTS, ',', S))) +
#   ylim(0, NA)+
#   theme_bw()+
#   theme(legend.position='none') + 
#   f_getCustomTheme() 
# grid.arrange(gg1,gg3,gg2,gg4, nrow=2)


########################################################################################
# check whether incidence as a function of EIR relationship is okay
########################################################################################
exp_name = 'generic_eir_sweep_SMC'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]

# subset to plotted scenario
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   rtss_coverage > 0.001,
                   rtss_scenario == cur_rtss_scenario
)
ggplot(pe_df_cur, aes(x = Annual_EIR, y = clinical_cases_ref)) +
  geom_point(aes(col = Annual_EIR, shape = as.factor(smc_coverage)), size = 3) +
  ylab('annual incidence (per 1,000 U5)') +
  xlab('annual EIR') +
  theme_bw() +
  f_getCustomTheme()


exp_name = 'generic_season_sweep_EIR_10_30'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]
# subset to plotted scenario
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   rtss_coverage > 0.001,
                   rtss_scenario == cur_rtss_scenario
)
ggplot(pe_df_cur, aes(x = Annual_EIR, y = clinical_cases_ref)) +
  geom_point(aes(col = seasonality, shape = as.factor(smc_coverage)), size = 3) +
  ylab('annual incidence (per 1,000 U5)') +
  xlab('annual EIR') +
  theme_bw() +
  f_getCustomTheme()
# xvar = 'Annual_EIR', yvar = 'rtss_cases_averted_per100000', bargroup_var='rtss_scenario',
# fillvar = 'rtss_scenario', facet1 = 'seasonality', facet2 = 'smc_coverage',
# SAVE = FALSE, ylab='cases averted by RTS,S per 100,000')

age_df = sim_output[[4]]
age_df = f_add_scenario_name(df = age_df, scenario_type = 'rtss_scenario')
age_df = age_df[!is.na(age_df$rtss_scenario),]
# subset to plotted scenario
pe_df_cur = filter(age_df,
                   age_group == '5 - 6',
                   rtss_coverage > 0.001,
                   rtss_scenario == cur_rtss_scenario
)
ggplot(pe_df_cur, aes(x = Annual_EIR, y = clinical_cases_ref)) +
  geom_point(aes(col = seasonality, shape = as.factor(smc_coverage)), size = 3) +
  ylab('annual incidence (per 1,000 aged 5-6)') +
  xlab('annual EIR') +
  theme_bw() +
  f_getCustomTheme()


####################################################################################################
# plot the age when children in each cohort first receives the RTS,S campaign booster
####################################################################################################

rtss_booster1_min_age = 2 * 365
start_day0 = round(365 + 6 * 30.4 - 7) # one week before SMC
cohort_month_shift = 0:11
start_days = start_day0 - round(30.4 * cohort_month_shift)
start_days_adjusted = c()
for (ss in 1:length(start_days)) {
  start_day = start_days[ss]
  while (start_day < rtss_booster1_min_age) {
    start_day = start_day + 365
  }
  start_days_adjusted[ss] = start_day
}
png(filename = paste0(simout_dir, '/age_at_first_booster.png'), width = 4, height = 4, units = 'in', res = 900)
plot((cohort_month_shift + 1), start_days_adjusted, ylim = c(365 * 2 - 30, 365 * 3), xlab = 'birth month', ylab = c('age (in days) when a child', 'receives their first RTS,S booster'), bty = 'L', pch = 20, cex = 1.8)
abline(v = (start_day0 %% 365) / 30.4 + 1, col = 'lightblue', lwd = 4)
points((cohort_month_shift + 1), rep(365 * 2, length(cohort_month_shift)), pch = 20, cex = 1.8, col = 'grey')
legend('topleft', c('standard', 'enhanced'), pch = 20, pt.cex = 1.8, col = c('grey', 'black'), bty = 'n')
dev.off()


gg <- ggplot() +
  geom_point(aes(x = (cohort_month_shift + 1), y = start_days_adjusted), size = 2.3) +
  geom_point(aes(x = (cohort_month_shift + 1), y = rep(365 * 2, length(cohort_month_shift))), size = 2.3, col = 'grey') +
  scale_y_continuous(lim = c(365 * 2 - 30, 365 * 3)) +
  scale_x_continuous(breaks = cohort_month_shift + 1) +
  geom_vline(xintercept = (start_day0 %% 365) / 30.4 + 1, col = 'lightblue', lwd = 2) +
  labs(x = 'birth month', y = 'age (in days) when a child\nreceives their first RTS,S booster') +
  theme_cowplot() +
  customTheme

f_save_plot(gg, paste0('age_at_first_booster_v2'), file.path(simout_dir, '_plots'), width = 5, height = 4, units = 'in', device_format = device_format)