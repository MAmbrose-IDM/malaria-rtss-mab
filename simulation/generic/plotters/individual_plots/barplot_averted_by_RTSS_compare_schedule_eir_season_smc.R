# barplot_averted_by_RTSS_compare_schedule_eir_season_smc.R


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
                'Cohort_birth_month', 'minBoostAge', 'cm_target_group', 'smc_target_group', 'rtss_target_group')
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
# grouped barplots showing cases averted when adding various coverages/schedules of RTS,S 
#    against different backgrounds (seasonality, EIR, and SMC)
#####################################################################################################
exp_name = 'generic_rtss_age_sweep'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))


# ===  show bars with all RTS,S schedules and coverages; compare across EIRs, seasonalities, and SMC coverages  === #
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
# clinical cases
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

# severe cases
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





# ===  simple barplot with versus without SMC (create panel with clinical/severe and number/percent)  === #

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

# number of clinical cases
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

# number of severe cases
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

# percent reduction clinical case
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

# percent reduction severe case
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

# combine plots into panel
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

# clinical cases
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_cases_averted_per100000', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = 'seasonality', facet2 = NULL,
                           SAVE = FALSE, ylab = 'cases averted by \n RTS,S per 100,000')
gg = gg + f_getColorMapping(name = 'scenario', type = 'fill')
f_save_plot(gg, paste0('simple_cases_averted_by_80RTSS_stand_enh_with_without_80SMC_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 6, height = 2.8, units = 'in', device_format = device_format)

# severe cases
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_severe_cases_averted_per100000', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = 'seasonality', facet2 = NULL,
                           SAVE = FALSE, ylab = 'severe cases averted by \n RTS,S per 100,000')
gg = gg + f_getColorMapping(name = 'scenario', type = 'fill')
f_save_plot(gg, paste0('simple_severe_cases_averted_by_80RTSS_stand_enh_with_without_80SMC_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 5.9, height = 2.8, units = 'in', device_format = device_format)





# ===  barplot comparing three RTS,S schedules at high versus constant seasonality, no SMC, severe and clinical cases  === #
exp_name = 'generic_standard_enhanced_season'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))

# subset simulations
cur_cm = 0.6
# cur_corr = 0
cur_age = 'U5'
cur_smcs = c(0)
cur_eirs = c(5, 10, 30)
cur_seasonalities = c('constant', 'high_unimodal')
cur_rtss_scenarios = c('80% enhanced RTS,S (1 booster, min 18m)', '80% enhanced RTS,S (1 booster, min 24m)', '80% standard RTS,S (1 booster)')
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   # intervention_correlation == cur_corr,
                   # seasonality %in% cur_seasonalities,
                   Annual_EIR %in% cur_eirs,
                   cm_coverage == cur_cm,
                   rtss_coverage > 0,
                   smc_coverage %in% cur_smcs,
                   rtss_scenario %in% cur_rtss_scenarios
)
pe_df_cur$rtss_scenario = factor(pe_df_cur$rtss_scenario, levels=cur_rtss_scenarios)
# percent reduction for simple comparison
pe_df_cur$rtss_percent_reduction = pe_df_cur$rtss_relative_burden * -100
pe_df_cur$rtss_percent_reduction_severe = pe_df_cur$rtss_relative_burden_severe * -100

# clinical cases
gg1 = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_percent_reduction', bargroup_var = 'rtss_scenario',
                           fillvar = 'rtss_scenario', facet1 = 'seasonality', facet2 = NULL,
                           SAVE = FALSE, ylab = 'percent reduction in clinical cases from adding RTS,S', scales='fixed', bar_width=0.85) + 
  scale_fill_manual(values = c('#73bfe2',
                               '#78c26d',
                               '#ca5800'))

# severe cases
gg2 = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_percent_reduction_severe', bargroup_var = 'rtss_scenario',
                           fillvar = 'rtss_scenario', facet1 = 'seasonality', facet2 = NULL,
                           SAVE = FALSE, ylab = 'percent reduction in severe cases from adding RTS,S', scales='fixed', bar_width=0.85) + 
  scale_fill_manual(values = c('#73bfe2',
                               '#78c26d',
                               '#ca5800'))

# combine plots into panel
gg_legend <- get_legend(gg1)
gg1 <- gg1  + theme(legend.position = 'None')
gg2 <- gg2  + theme(legend.position = 'None')
gg <- plot_grid(gg1,gg2, align='hv', nrow=2)
gg <- plot_grid(gg, gg_legend, rel_widths = c(1, 0.25))
f_save_plot(gg, paste0('cases_averted_by_RTSS_schedules_for_eir_seasonality_0SMC',
                       cur_age, '_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            # '_minAge24m'),
            file.path(simout_dir, '_plots'), width = 8.5, height = 3.8, units = 'in', device_format = device_format)






# ===  LINEPLOT ACROSS EIRS comparing three RTS,S schedules at high versus constant seasonality, no SMC, severe and clinical cases  === #
exp_name = 'generic_standard_enhanced_season'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))

# subset simulations
cur_cm = 0.6
# cur_corr = 0
cur_age = 'U5'
cur_smcs = c(0)
cur_seasonalities = c('constant', 'high_unimodal')
cur_rtss_scenarios = c('80% enhanced RTS,S (1 booster, min 18m)', '80% enhanced RTS,S (1 booster, min 24m)', '80% standard RTS,S (1 booster)')
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   # intervention_correlation == cur_corr,
                   seasonality %in% cur_seasonalities,
                   cm_coverage == cur_cm,
                   rtss_coverage > 0,
                   smc_coverage %in% cur_smcs,
                   rtss_scenario %in% cur_rtss_scenarios
)
pe_df_cur$rtss_scenario = factor(pe_df_cur$rtss_scenario, levels=cur_rtss_scenarios)
# percent reduction for simple comparison
pe_df_cur$rtss_percent_reduction = pe_df_cur$rtss_relative_burden * -100
pe_df_cur$rtss_percent_reduction_severe = pe_df_cur$rtss_relative_burden_severe * -100


xvar='Annual_EIR'
xlab='annual EIR'
xvar='clinical_cases_no_rtss'
xlab='pre-RTS,S incidence'
colvar='rtss_scenario'
facet1='seasonality'
# clinical cases - percent reduction
yvar='rtss_percent_reduction'
ylab = 'percent reduction in clinical cases from adding RTS,S'
gg1_p = plot_lines(dat=pe_df_cur, xvar=xvar, yvar=yvar,
                   colvar=colvar, facet1=facet1, facet2=NULL,
                   SAVE = FALSE, xlab=xlab, ylab=ylab,
                   add_watermark=FALSE) + 
  scale_color_manual(values = c('#73bfe2',
                               '#78c26d',
                               '#ca5800'))

# clinical cases - number reduction
yvar='rtss_cases_averted_per100000'
ylab = 'clinical cases averted (per 100,000 U5)\n from adding RTS,S'
gg1_n = plot_lines(dat=pe_df_cur, xvar=xvar, yvar=yvar,
                   colvar=colvar, facet1=facet1, facet2=NULL,
                   SAVE = FALSE, xlab=xlab, ylab=ylab,
                   add_watermark=FALSE) + 
  scale_color_manual(values = c('#73bfe2',
                                '#78c26d',
                                '#ca5800'))
  
  # ggplot(data = pe_df_cur, aes(x=get(xvar), y=get(yvar), group=get(colvar))) + 
  # geom_line(aes(color=get(colvar)), size=2.5) +
  # xlab(xlab)+
  # ylab(ylab) +
  # facet_wrap(~get(facetvar))+
  # scale_color_manual(values = c('#73bfe2',
  #                               '#78c26d',
  #                               '#ca5800'))



# severe cases - percent reduction
yvar='rtss_percent_reduction_severe'
ylab = 'percent reduction in severe cases from adding RTS,S'
gg2_p = plot_lines(dat=pe_df_cur, xvar=xvar, yvar=yvar,
                   colvar=colvar, facet1=facet1, facet2=NULL,
                   SAVE = FALSE, xlab=xlab, ylab=ylab,
                   add_watermark=FALSE) + 
  scale_color_manual(values = c('#73bfe2',
                                '#78c26d',
                                '#ca5800'))

# severe cases - number reduction
yvar='rtss_severe_cases_averted_per100000'
ylab = 'severe cases averted (per 100,000 U5)\n from adding RTS,S'
gg2_n = plot_lines(dat=pe_df_cur, xvar=xvar, yvar=yvar,
                   colvar=colvar, facet1=facet1, facet2=NULL,
                   SAVE = FALSE, xlab=xlab, ylab=ylab,
                   add_watermark=FALSE) + 
  scale_color_manual(values = c('#73bfe2',
                                '#78c26d',
                                '#ca5800'))


# combine plots into panel
gg_legend <- get_legend(gg1_n)
gg1_n <- gg1_n  + theme(legend.position = 'None')
gg2_n <- gg2_n  + theme(legend.position = 'None')
gg <- plot_grid(gg1_n, gg2_n, align='hv', nrow=2)
gg <- plot_grid(gg, gg_legend, rel_widths = c(1, 0.25))

f_save_plot(gg, paste0('lineplot_cases_averted_by_RTSS_schedules_for_eir_seasonality_0SMC',
                       cur_age, '_',
                       cur_cm * 100, 'CM_'),
            # '_minAge24m'),
            file.path(simout_dir, '_plots'), width = 8, height = 6, units = 'in', device_format = device_format)




# ===  with two SMC coverages versus without SMC (panel with cases/severe and number/percent reduction)  === #
exp_name = 'generic_eir_smc_sweep3'  # 'generic_rtss_age_sweep'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
# subset simulations
smc_rtss_covs = c(0.5, 0.8)
cur_cm = 0.6
cur_cm_targets = 'random'
cur_smc_targets = 'high'
cur_rtss_targets = 'random'
cur_smcs = c(0, smc_rtss_covs)
cur_rtss = c(smc_rtss_covs)
cur_season = 'high_unimodal'
cur_eirs = c(1, 10, 30)
cur_age = 'U5'
pe_df_cur = filter(pe_df,
                   seasonality == cur_season,
                   age_group == cur_age,
                   cm_coverage == cur_cm,
                   rtss_coverage %in% cur_rtss,
                   smc_coverage %in% cur_smcs,
                   Annual_EIR %in% cur_eirs,
                   # cm_target_group %in% cur_cm_targets,
                   # smc_target_group %in% cur_smc_targets,
                   # rtss_target_group %in% cur_rtss_targets
)

# percent reduction for simple comparison
pe_df_cur$rtss_percent_reduction = pe_df_cur$rtss_relative_burden * -100
pe_df_cur$rtss_percent_reduction_severe = pe_df_cur$rtss_relative_burden_severe * -100

levels = c('no SMC', '50% SMC coverage', '80% SMC coverage')
pe_df_cur$plotted_scen = paste0(round(100*(pe_df_cur$smc_coverage)),'% SMC coverage')
pe_df_cur$plotted_scen[pe_df_cur$plotted_scen == '0% SMC coverage'] = 'no SMC'
pe_df_cur$plotted_scen = factor(pe_df_cur$plotted_scen, levels = levels)
# specify desired colors for each level
myColors = c(rgb(0, 0, 0),
             rgb(0.1, 0.6, 0.3),
             rgb(0.5, 1, 0.7))
names(myColors) = levels
colorMapping = scale_fill_manual(name = 'plotted_scen', values = myColors)



# number of clinical cases
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_cases_averted_per100000', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = NULL, facet2 = NULL,
                           SAVE = FALSE, ylab = 'cases averted by\nRTS,S per 100,000 U5')
gg = gg + colorMapping
gg1 <- gg

# number of severe cases
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_severe_cases_averted_per100000', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = NULL, facet2 = NULL,
                           SAVE = FALSE, ylab = 'severe cases averted by\nRTS,S per 100,000 U5')
gg = gg + colorMapping
gg2 <- gg



# percent reduction clinical case
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_percent_reduction', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = NULL, facet2 = NULL,  #'smc_coverage',
                           SAVE = FALSE, ylab = 'percent of cases U5\n averted by RTS,S')
gg = gg + colorMapping
gg3 <- gg

# percent reduction severe case
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_percent_reduction_severe', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = NULL, facet2 = NULL, # 'smc_coverage',
                           SAVE = FALSE, ylab = 'percent of severe cases U5\n averted by RTS,S')
gg = gg + colorMapping
gg4 <- gg

# combine plots into panel
gg_legend <- get_legend(gg1)
gg1 <- gg1  + theme(legend.position = 'None')
gg2 <- gg2  + theme(legend.position = 'None')
gg3 <- gg3  + theme(legend.position = 'None')
gg4 <- gg4  + theme(legend.position = 'None')
gg <- plot_grid(gg1,gg2,gg3,gg4, align='hv')
gg <- plot_grid(gg, gg_legend, rel_widths = c(1, 0.25))

f_save_plot(gg, paste0('cases_severe_averted_by_80RTSSrandom_with_without_SMChigh_highSeasonality_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_'),
            file.path(simout_dir, '_plots'), width = 7.5, height = 5, units = 'in', device_format = device_format)
rm(gg, gg1, gg2, gg3,gg4)







# ===  with three CM coverages (without SMC) (panel with cases/severe and number/percent reduction)  === #
exp_name = 'generic_rtss_eir_cm_sweep'  # 'generic_rtss_age_sweep'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
# subset simulations
cur_cm_targets = 'random'
cur_smc_targets = 'high'
cur_rtss_targets = 'random'
cur_smcs = c(0, smc_rtss_covs)
cur_rtss = 0.8
cur_season = 'high_unimodal'
cur_eirs = c(1, 10, 30)
cur_age = 'U5'
pe_df_cur = filter(pe_df,
                   seasonality == cur_season,
                   age_group == cur_age,
                   rtss_coverage %in% cur_rtss,
                   smc_coverage < 0.001,
                   Annual_EIR %in% cur_eirs,
                   # cm_target_group %in% cur_cm_targets,
                   # smc_target_group %in% cur_smc_targets,
                   # rtss_target_group %in% cur_rtss_targets
)

# percent reduction for simple comparison
pe_df_cur$rtss_percent_reduction = pe_df_cur$rtss_relative_burden * -100
pe_df_cur$rtss_percent_reduction_severe = pe_df_cur$rtss_relative_burden_severe * -100

levels = c('30%', '60%', '90%')
pe_df_cur$plotted_scen = paste0(round(100*(pe_df_cur$cm_coverage)),'%')
pe_df_cur$plotted_scen = factor(pe_df_cur$plotted_scen, levels = levels)
# specify desired colors for each level
myColors = c(rgb(0.6, 0.1, 0.3),
             rgb(0, 0, 0),
             rgb(1, 0.5, 0.7))
names(myColors) = levels
colorMapping = scale_fill_manual(name = 'plotted_scen', values = myColors)



# number of clinical cases
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_cases_averted_per100000', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = NULL, facet2 = NULL,
                           SAVE = FALSE, ylab = 'cases averted by\nRTS,S per 100,000 U5')
gg = gg + colorMapping
gg1 <- gg

# number of severe cases
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_severe_cases_averted_per100000', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = NULL, facet2 = NULL,
                           SAVE = FALSE, ylab = 'severe cases averted by\nRTS,S per 100,000 U5')
gg = gg + colorMapping
gg2 <- gg



# percent reduction clinical case
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_percent_reduction', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = NULL, facet2 = NULL,  #'smc_coverage',
                           SAVE = FALSE, ylab = 'percent of cases U5\n averted by RTS,S')
gg = gg + colorMapping
gg3 <- gg

# percent reduction severe case
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_percent_reduction_severe', bargroup_var = 'plotted_scen',
                           fillvar = 'plotted_scen', facet1 = NULL, facet2 = NULL, # 'smc_coverage',
                           SAVE = FALSE, ylab = 'percent of severe cases U5\n averted by RTS,S')
gg = gg + colorMapping
gg4 <- gg

# combine plots into panel
gg_legend <- get_legend(gg1)
gg1 <- gg1  + theme(legend.position = 'None')
gg2 <- gg2  + theme(legend.position = 'None')
gg3 <- gg3  + theme(legend.position = 'None')
gg4 <- gg4  + theme(legend.position = 'None')
gg <- plot_grid(gg1,gg2,gg3,gg4, align='hv')
gg <- plot_grid(gg, gg_legend, rel_widths = c(1, 0.25))

f_save_plot(gg, paste0('cases_severe_averted_by_80RTSSrandom_withWithoutCMrandom_highSeasonality_',
                       cur_age, '_',
                       cur_cm * 100, 'CM_'),
            file.path(simout_dir, '_plots'), width = 7.5, height = 5, units = 'in', device_format = device_format)
rm(gg, gg1, gg2, gg3,gg4)






##########################################################################################################################
# ===  barplot with versus without SMC at different correlations (create panel with clinical/severe and number/percent)  === #
##########################################################################################################################
exp_name = 'generic_accesscorrelation_highSMC' # 'generic_accesscorrelation_highSMC'  #_highCM_cleaned' # 'generic_accesscorrelation05_v3'  # 'generic_accesscorrelation05_independentCM_v2'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'scenario_RTSS_SMC_withCorr2')
# subset simulations
smc_rtss_cov = 0.8
cur_cm = 0.6
cur_age = 'U5'
cur_smcs = c(0, smc_rtss_cov)
cur_rtss = c(smc_rtss_cov)
cur_season = 'high_unimodal'
cur_cm_target = 'high'
pe_df_cur = filter(pe_df,
                   seasonality == cur_season,
                   age_group == cur_age,
                   cm_coverage == cur_cm,
                   rtss_coverage %in% cur_rtss,
                   smc_coverage %in% cur_smcs,
                   cm_target_group == cur_cm_target
                   
)

# number of clinical cases
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_cases_averted_per100000', bargroup_var = 'scenario_name',
                           fillvar = 'scenario_name', facet1 = NULL, facet2 = NULL,
                           SAVE = FALSE, ylab = 'cases averted by\nRTS,S per 100,000')
gg = gg + f_getColorMapping(name = 'scenario', type = 'fill')
gg1 <- gg


# plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'cases_averted_per100000', bargroup_var = 'scenario_name',
#                       fillvar = 'scenario_name', facet1 = NULL, facet2 = NULL,
#                       SAVE = FALSE, ylab = 'cases averted by\nRTS,S and SMC per 100,000') +
#   f_getColorMapping(name = 'scenario', type = 'fill')


# number of severe cases
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_severe_cases_averted_per100000', bargroup_var = 'scenario_name',
                           fillvar = 'scenario_name', facet1 = NULL, facet2 = NULL,
                           SAVE = FALSE, ylab = 'severe cases averted by\nRTS,S per 100,000')
gg = gg + f_getColorMapping(name = 'scenario', type = 'fill')
gg2 <- gg

# ==== percent reduction for simple comparison ==== #
pe_df_cur$rtss_percent_reduction = pe_df_cur$rtss_relative_burden * -100
pe_df_cur$rtss_percent_reduction_severe = pe_df_cur$rtss_relative_burden_severe * -100

# percent reduction clinical case
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_percent_reduction', bargroup_var = 'scenario_name',
                           fillvar = 'scenario_name', facet1 = NULL, facet2 = NULL,  #'smc_coverage',
                           SAVE = FALSE, ylab = 'percent reduction in cases U5\nfrom adding RTS,S')
gg = gg + f_getColorMapping(name = 'scenario', type = 'fill')
gg3 <- gg

# percent reduction severe case
gg = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'rtss_percent_reduction_severe', bargroup_var = 'scenario_name',
                           fillvar = 'scenario_name', facet1 = NULL, facet2 = NULL, # 'smc_coverage',
                           SAVE = FALSE, ylab = 'percent reduction in severe cases U5\nfrom adding RTS,S')
gg = gg + f_getColorMapping(name = 'scenario', type = 'fill')
gg4 <- gg

# combine plots into panel
gg_legend <- get_legend(gg1)
gg1 <- gg1  + theme(legend.position = 'None')
gg2 <- gg2  + theme(legend.position = 'None')
gg3 <- gg3  + theme(legend.position = 'None')
gg4 <- gg4  + theme(legend.position = 'None')
gg <- plot_grid(gg1,gg2,gg3,gg4, align=c('hv'))
gg <- plot_grid(gg, gg_legend, rel_widths = c(1, 0.5))
f_save_plot(gg, paste0('cases_percent_reduction_RTSS_combined_accessCorr', '_minAge24m', '_', smc_rtss_cov, 'RTSSSMCcov', '_', cur_cm_target, 'CM'),
            file.path(simout_dir, '_plots'), width = 9.5, height =5, units = 'in', device_format = device_format)

rm(gg, gg1, gg2, gg3,gg4)







##########################################################################################################################
# ===  barplot with versus without SMC - different coverages, RTSS correlations ( panel with severe and number/percent)  === #
##########################################################################################################################
# each group of bars is a different RTSS/SMC coverage and CM correlation, rather than a different EIR


exp_name = 'generic_accesscorrelation_highSMC' # 'generic_accesscorrelation_highSMC'  #_highCM_cleaned' # 'generic_accesscorrelation05_v3'  # 'generic_accesscorrelation05_independentCM_v2'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'scenario_RTSS_SMC_withCorr3')
# subset simulations
smc_rtss_covs = c(0.6, 0.8)
cur_cm = 0.6
cur_cm_targets = 'random'  # c('random', 'high')
cur_smcs = c(0, smc_rtss_covs)
cur_rtss = c(smc_rtss_covs)
cur_season = 'high_unimodal'
cur_eir = 10
cur_age = 'U5'
pe_df_cur = filter(pe_df,
                   seasonality == cur_season,
                   Annual_EIR == cur_eir,
                   age_group == cur_age,
                   cm_coverage == cur_cm,
                   rtss_coverage %in% cur_rtss,
                   smc_coverage %in% cur_smcs,
                   cm_target_group %in% cur_cm_targets
                   
)

# create grouping variable along x-axis - out of SMC/RTSS coverage and CM targeting
if(length(cur_cm_targets)>1){
  pe_df_cur$x_group = paste0(pe_df_cur$rtss_coverage*100, '% RTSS & SMC; CM target: ', pe_df_cur$cm_target_group)
  pe_df_cur$x_group = factor(pe_df_cur$x_group, levels=c("60% RTSS & SMC; CM target: random", "60% RTSS & SMC; CM target: high", 
                                                         "80% RTSS & SMC; CM target: random", "80% RTSS & SMC; CM target: high" ))
} else{
  pe_df_cur$x_group = paste0(pe_df_cur$rtss_coverage*100, '% RTSS & SMC')
  pe_df_cur$x_group = factor(pe_df_cur$x_group, levels=c("60% RTSS & SMC", "80% RTSS & SMC" ))
}


# add the fraction of children not previously covered by an intervention that are now covered by RTSS
pe_df_cur$frac_newly_covered = NA
pe_df_cur$frac_newly_covered_no_smc = NA
pe_df_cur$frac_rtss_or_smc = NA
for (rr in 1:nrow(pe_df_cur)){
  inter_names = c('rtss','smc','cm')
  coverages = c(pe_df_cur$rtss_coverage[rr], pe_df_cur$smc_coverage[rr], pe_df_cur$cm_coverage[rr])
  targeting = c(pe_df_cur$rtss_target_group[rr], pe_df_cur$smc_target_group[rr], pe_df_cur$cm_target_group[rr])
  # pe_df_cur$frac_newly_covered[rr] = calc_frac_no_inter(inter_names, coverages, targeting, frac_in_high=0.5, createVenn = FALSE)[2]
  # pe_df_cur$frac_newly_covered_no_smc[rr] = calc_frac_no_inter(inter_names, coverages, targeting, frac_in_high=0.5, createVenn = FALSE)[3]
  pe_df_cur$frac_rtss_or_smc[rr] = calc_frac_no_inter(inter_names, coverages, targeting, frac_in_high=0.5, createVenn = FALSE)[4]
}
# ==== percent reduction for simple comparison ==== #
pe_df_cur$rtss_percent_reduction = pe_df_cur$rtss_relative_burden * -100
pe_df_cur$rtss_percent_reduction_severe = pe_df_cur$rtss_relative_burden_severe * -100



xvar='scenario_name'
fillvar='scenario_name'
bargroup_var='scenario_name'
yvar1='rtss_severe_cases_averted_per100000'
ylab1='severe cases U5 averted \n by RTS,S per 100,000'
yvar2='rtss_percent_reduction_severe'
ylab2='percent reduction in severe \n cases U5 from adding RTS,S'
# line_yvar='frac_newly_covered_no_smc'
# line_ylab='percent of those without \n SMC covered by RTS,S'
# line_ylab='percent of those without SMC \n or CM covered by RTS,S'
line_yvar='frac_rtss_or_smc'
line_ylab='percent of eligible children \n who receive SMC or RTS,S'
facet1='x_group'
facet2=NA
dat=pe_df_cur


# plot fraction of people who weren't covered by SMC who are now covered by RTS,S
gg0 = ggplot(data = dat, aes(x =as.factor(get(xvar)), y = get(line_yvar)*100)) +
  xlab(element_blank())+
  ylab(line_ylab)+
  ylim(0,100)+
  # customTheme +
  theme(panel.spacing = unit(0, "lines"), 
        # legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA, size = 0.5),
        strip.background =element_rect(fill="white")
  )+
  geom_line(aes(x =as.factor(get(xvar)), y = get(line_yvar)*100, group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))), col=rgb(0.5,0.5,0.5,0.1))+
  geom_point(aes(x =as.factor(get(xvar)), color=as.factor(get(fillvar)), y = get(line_yvar)*100, group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))), size=3.5)+
  f_getColorMapping(name = 'scenario', type = 'color') +
  facet_wrap(~get(facet1), nrow=1, strip.position = 'bottom')
if(line_yvar=='frac_rtss_or_smc'){  # add the fractions who get SMC without RTSS
  gg0 = gg0 +
    geom_line(aes(x =as.factor(get(xvar)), y = get('smc_coverage')*100, group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))), col=rgb(0.5,0.5,0.5,0.1))+
    geom_point(aes(x =as.factor(get(xvar)), y = get('smc_coverage')*100, group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))), color=rgb(100/255, 100/255, 100/255), size=3.5)
    
}

gg1 = ggplot(data = dat, aes(x =as.factor(get(xvar)), y = get(yvar1))) +
  geom_bar(stat = "identity", aes(fill = as.factor(get(fillvar)), group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))),
           position = position_dodge(width = 1), width=1) +   #'dodge') +
  # geom_point(aes(x =as.factor(get(xvar)), y = get(yvar1), group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))), color=rgb(0.5,0.5,0.5,0.5))+
  xlab(element_blank())+
  ylab(ylab1) +
  # customTheme +
  theme(panel.spacing = unit(0, "lines"), 
        # legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA, size = 0.5),
        strip.background =element_rect(fill="white")
  )+
  f_getColorMapping(name = 'scenario', type = 'fill') +
  facet_wrap(~get(facet1), nrow=1, strip.position = 'bottom')


gg2 = ggplot(data = dat, aes(x =as.factor(get(xvar)), y = get(yvar2))) +
  geom_bar(stat = "identity", aes(fill = as.factor(get(fillvar)), group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))),
           position = position_dodge(width = 1), width=1) +   #'dodge') +
  # geom_point(aes(x =as.factor(get(xvar)), y = get(yvar2), group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))), color=rgb(0.5,0.5,0.5,0.5))+
  xlab(element_blank())+
  ylab(ylab2) +
  # customTheme +
  theme(panel.spacing = unit(0, "lines"), 
        # legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA, size = 0.5),
        strip.background =element_rect(fill="white")
  )+
  f_getColorMapping(name = 'scenario', type = 'fill') +
  facet_wrap(~get(facet1), nrow=1, strip.position = 'bottom')



# combine plots into panel
gg_legend <- get_legend(gg1)
gg0 <- gg0  + theme(legend.position = 'None')
gg1 <- gg1  + theme(legend.position = 'None')
gg2 <- gg2  + theme(legend.position = 'None')
gg <- plot_grid(gg0, gg2,gg1, align=c('hv'), ncol=1, rel_heights = c(0.9,1,1))
gg <- plot_grid(gg, gg_legend, rel_widths = c(1, 0.25))

f_save_plot(gg, paste0('severe_percent_reduction_v4_RTSS_combined_accessCorr', '_minAge24m_', cur_eir, 'EIR', '_', 'multipleRTSSSMCcov', '_', cur_cm_targets, 'CM'),
            file.path(simout_dir, '_plots'), width = 4, height =7, units = 'in', device_format = device_format)





##########################################################################################################################
# ===  barplot with different correlations with CM and different RTS,S coverages ( panel with severe and number/percent)  === #
##########################################################################################################################
# each group of bars is a different RTSS/SMC coverage and CM correlation, rather than a different EIR


exp_name = 'generic_accesscorrelation_highSMC' # 'generic_accesscorrelation_highSMC'  #_highCM_cleaned' # 'generic_accesscorrelation05_v3'  # 'generic_accesscorrelation05_independentCM_v2'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'scenario_RTSS_CM_withCorr3')
# subset simulations
cur_cm = 0.6
cur_cm_targets = c('high')
cur_smcs = 0
cur_rtss = c(smc_rtss_covs)
cur_season = 'high_unimodal'
cur_eir = 10
cur_age = 'U5'
pe_df_cur = filter(pe_df,
                   seasonality == cur_season,
                   Annual_EIR == cur_eir,
                   age_group == cur_age,
                   cm_coverage == cur_cm,
                   rtss_coverage > 0.001,
                   smc_coverage %in% cur_smcs,
                   cm_target_group %in% cur_cm_targets
)

# create grouping variable along x-axis - out of SMC/RTSS coverage and CM targeting
if(length(cur_cm_targets)>1){
  pe_df_cur$x_group = paste0(pe_df_cur$rtss_coverage*100, '% RTSS; CM target: ', pe_df_cur$cm_target_group)
  pe_df_cur$x_group = factor(pe_df_cur$x_group, levels=c("60% RTSS; CM target: random", "60% RTSS; CM target: high", 
                                                         "80% RTSS; CM target: random", "80% RTSS; CM target: high" ))
} else{
  pe_df_cur$x_group = paste0(pe_df_cur$rtss_coverage*100, '% RTSS')
  pe_df_cur$x_group = factor(pe_df_cur$x_group, levels=c("60% RTSS", "80% RTSS" ))
}


# add the fraction of children not previously covered by an intervention that are now covered by RTSS
pe_df_cur$frac_rtss_or_cm = NA
for (rr in 1:nrow(pe_df_cur)){
  inter_names = c('rtss','smc','cm')
  coverages = c(pe_df_cur$rtss_coverage[rr], pe_df_cur$smc_coverage[rr], pe_df_cur$cm_coverage[rr])
  targeting = c(pe_df_cur$rtss_target_group[rr], pe_df_cur$smc_target_group[rr], pe_df_cur$cm_target_group[rr])
  pe_df_cur$frac_rtss_or_cm[rr] = calc_frac_no_inter(inter_names, coverages, targeting, frac_in_high=0.5, createVenn = FALSE)[5]
}
# ==== percent reduction for simple comparison ==== #
pe_df_cur$rtss_percent_reduction = pe_df_cur$rtss_relative_burden * -100
pe_df_cur$rtss_percent_reduction_severe = pe_df_cur$rtss_relative_burden_severe * -100



xvar='scenario_name'
fillvar='scenario_name'
bargroup_var='scenario_name'
yvar1='rtss_severe_cases_averted_per100000'
ylab1='severe cases U5 averted \n by RTS,S per 100,000'
yvar2='rtss_percent_reduction_severe'
ylab2='percent reduction in severe \n cases U5 from adding RTS,S'
# line_yvar='frac_newly_covered_no_smc'
# line_ylab='percent of those without \n SMC covered by RTS,S'
# line_ylab='percent of those without SMC \n or CM covered by RTS,S'
line_yvar='frac_rtss_or_cm'
line_ylab='percent of eligible children \n who receive CM or RTS,S'
facet1='x_group'
facet2=NA
dat=pe_df_cur


# plot fraction of people who weren't covered by SMC who are now covered by RTS,S
gg0 = ggplot(data = dat, aes(x =as.factor(get(xvar)), y = get(line_yvar)*100)) +
  xlab(element_blank())+
  ylab(line_ylab)+
  ylim(0,100)+
  # customTheme +
  theme(panel.spacing = unit(0, "lines"), 
        # legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA, size = 0.5),
        strip.background =element_rect(fill="white")
  )+
  geom_line(aes(x =as.factor(get(xvar)), y = get(line_yvar)*100, group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))), col=rgb(0.5,0.5,0.5,0.1))+
  geom_point(aes(x =as.factor(get(xvar)), color=as.factor(get(fillvar)), y = get(line_yvar)*100, group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))), size=3.5)+
  f_getColorMapping(name = 'scenario', type = 'color') +
  facet_wrap(~get(facet1), nrow=1, strip.position = 'bottom')
if(line_yvar=='frac_rtss_or_cm'){  # add the fractions who get SMC without RTSS
  gg0 = gg0 +
    geom_line(aes(x =as.factor(get(xvar)), y = get('cm_coverage')*100, group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))), col=rgb(0.5,0.5,0.5,0.1))+
    geom_point(aes(x =as.factor(get(xvar)), y = get('cm_coverage')*100, group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))), color=rgb(100/255, 100/255, 100/255), size=3.5)
}

gg1 = ggplot(data = dat, aes(x =as.factor(get(xvar)), y = get(yvar1))) +
  geom_bar(stat = "identity", aes(fill = as.factor(get(fillvar)), group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))),
           position = position_dodge(width = 1), width=1) +   #'dodge') +
  # geom_point(aes(x =as.factor(get(xvar)), y = get(yvar1), group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))), color=rgb(0.5,0.5,0.5,0.5))+
  xlab(element_blank())+
  ylab(ylab1) +
  # customTheme +
  theme(panel.spacing = unit(0, "lines"), 
        # legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA, size = 0.5),
        strip.background =element_rect(fill="white")
  )+
  f_getColorMapping(name = 'scenario', type = 'fill') +
  facet_wrap(~get(facet1), nrow=1, strip.position = 'bottom')


gg2 = ggplot(data = dat, aes(x =as.factor(get(xvar)), y = get(yvar2))) +
  geom_bar(stat = "identity", aes(fill = as.factor(get(fillvar)), group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))),
           position = position_dodge(width = 1), width=1) +   #'dodge') +
  # geom_point(aes(x =as.factor(get(xvar)), y = get(yvar2), group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))), color=rgb(0.5,0.5,0.5,0.5))+
  xlab(element_blank())+
  ylab(ylab2) +
  # customTheme +
  theme(panel.spacing = unit(0, "lines"), 
        # legend.position="none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = NA, fill = NA, size = 0.5),
        strip.background =element_rect(fill="white")
  )+
  f_getColorMapping(name = 'scenario', type = 'fill') +
  facet_wrap(~get(facet1), nrow=1, strip.position = 'bottom')



# combine plots into panel
gg_legend <- get_legend(gg1)
gg0 <- gg0  + theme(legend.position = 'None')
gg1 <- gg1  + theme(legend.position = 'None')
gg2 <- gg2  + theme(legend.position = 'None')
gg <- plot_grid(gg0, gg2,gg1, align=c('hv'), ncol=1, rel_heights = c(0.9,1,1))
gg <- plot_grid(gg, gg_legend, rel_widths = c(1, 0.25))

f_save_plot(gg, paste0('severe_percent_reduction_v4_RTSS_CM_accessCorr', '_minAge24m_', cur_eir, 'EIR', '_', '0SMC', '_', cur_cm_targets, 'CM'),
            file.path(simout_dir, '_plots'), width = 4, height =7, units = 'in', device_format = device_format)


