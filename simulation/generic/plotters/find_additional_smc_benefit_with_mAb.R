# compare_sim_mAb_SMC.R
# August 2022
# create plots showing performance of mAbs relative to RTS,S at averting burden

library(lubridate)
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(RColorBrewer)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Setup filepaths and parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setwd("C:/Users/moniqueam/Documents/malaria-rtss-mab")

source(file.path("simulation", "load_paths.R"))
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path('simulation', 'generic', 'plotters', 'plot_bars_and_lines.R'))
device_format = c('pdf', 'png')
paths = get_project_paths()
datapath = paths[1]
projectpath = paths[2]

exp_name_comp1 = 'sweep4_seeds1'
exp_name_comp2 = 'sweep5_seeds1'

simout_dir=file.path(projectpath, 'simulation_output')
if (!dir.exists(paste0(simout_dir, '/_plots'))) dir.create(paste0(simout_dir, '/_plots'))
if (!dir.exists(paste0(simout_dir, '/_plots/pdf'))) dir.create(paste0(simout_dir, '/_plots/pdf'))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Process simulation output and create heat plots comparing cases 
#   averted with SMC versus mAbs across mAb params
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
sim_output = load_Age_monthly_Cases(simout_dir=simout_dir, exp_name=c(exp_name_comp1, exp_name_comp2), add_PE_perAge=TRUE,
                                    max_years=c(5, 8), keep_birth_month=FALSE, fname='All_Age_monthly_Cases.csv')

pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'vacc_info')
pe_df = pe_df[!is.na(pe_df$vacc_info),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))

# compare mAbs with SMC against matched mAb-alone impact. 'Reference' is mAb alone.
# subset to simulations with SMC
pe_df_smc = pe_df[pe_df$smc_coverage > 0,]
pe_df_smc = pe_df_smc[pe_df_smc$vacc_coverage > 0,]
match_cols = c('Annual_EIR', 'seasonality', 'age_group', 'frac_high_access', 'vacc_char', 'vacc_coverage',
               'cm_coverage', 'frac_high_access',
               'cm_target_group', 'smc_target_group', 'vacc_target_group')
keep_cols = c('cases_averted_per100000', 'severe_cases_averted_per100000')
df_smc = pe_df_smc %>%
  dplyr::select(c(match_cols, keep_cols)) %>%
  rename(cases_averted_per100000_w_smc = cases_averted_per100000,
         severe_cases_averted_per100000_w_smc = severe_cases_averted_per100000)
# subset to mAbs-alone (the 'reference')
pe_df_mab_only = pe_df[pe_df$vacc_type %in% c('mab', 'mAb'),]
pe_df_mab_only = pe_df_mab_only[pe_df_mab_only$vacc_coverage > 0,]
pe_df_mab_only = pe_df_mab_only[pe_df_mab_only$smc_coverage < 0.001,]

# merge the reference values into the data frame
pe_df_compare <- pe_df_mab_only %>% left_join(df_smc)
pe_df_compare$additional_cases_averted_with_smc = pe_df_compare$cases_averted_per100000_w_smc - pe_df_compare$cases_averted_per100000
pe_df_compare$additional_severe_cases_averted_with_smc = pe_df_compare$severe_cases_averted_per100000_w_smc - pe_df_compare$severe_cases_averted_per100000


# subset simulations
cur_cm = 0.6
cur_age = 'U5'
cur_eirs = c(5, 10, 30)
cur_seasonalities = c('moderate_unimodal', 'high_unimodal')
pe_df_cur = filter(pe_df_compare,
                   age_group == cur_age,
                   seasonality %in% cur_seasonalities,
                   Annual_EIR %in% cur_eirs,
                   cm_coverage == cur_cm
)
pe_df_cur$seasonality = factor(pe_df_cur$seasonality, levels=cur_seasonalities)
# 
# # clinical cases
# gg1 = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'cases_averted_compared_to_smc', bargroup_var = 'vacc_info',
#                             fillvar = 'vacc_info', facet1 = 'seasonality', facet2 = NULL,
#                             SAVE = FALSE, ylab = 'cases averted (per 100k U5) \nwith mAb compared to SMC', scales='fixed', bar_width=0.85) 
# # severe cases
# gg2 = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'severe_cases_averted_compared_to_smc', bargroup_var = 'vacc_info',
#                             fillvar = 'vacc_info', facet1 = 'seasonality', facet2 = NULL,
#                             SAVE = FALSE, ylab = 'severe cases averted (per 100k U5) \nwith mAb compared to SMC', scales='fixed', bar_width=0.85)
# 
# # combine plots into panel
# gg_legend = get_legend(gg1)
# gg1 = gg1  + theme(legend.position = 'None')
# gg2 = gg2  + theme(legend.position = 'None')
# gg = plot_grid(gg1,gg2, align='hv', nrow=2)
# gg = plot_grid(gg, gg_legend, rel_widths = c(1, 0.25))
# f_save_plot(gg, paste0('cases_averted_mab_compared_smc_for_eir_seasonality',
#                        cur_age, '_',
#                        cur_cm * 100, 'CM'),
#             file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)



# create grid showing additional cases averted with SMC compared to mAbs alone
#  which set of parameters do better and worse than SMC?
# keep_cols = c('Annual_EIR', 'seasonality', 'vacc_char', 'vacc_coverage', 'vacc_target_group', 'smc_target_group', 'hh', 'max_efficacy', 'vacc_info', 'mAb_better_clinical', 'mAb_better_severe', 'cases_averted_compared_to_smc', 'severe_cases_averted_compared_to_smc')
compare_mab_smc = pe_df_cur #%>% dplyr::select(keep_cols)
compare_mab_smc = compare_mab_smc[!is.na(compare_mab_smc$additional_cases_averted_with_smc),]

extreme_cases = max(compare_mab_smc$additional_cases_averted_with_smc)
extreme_cases = ceiling(extreme_cases/10000)*10000
breaks_cases = seq(0, extreme_cases, length.out=8)
colors_cases = brewer.pal(n=length(breaks_cases), name='PRGn')
gg3 = ggplot(compare_mab_smc, aes(x=hh, y=max_efficacy, z=additional_cases_averted_with_smc)) +
  geom_contour_filled(na.rm=TRUE, breaks=breaks_cases) +
  scale_fill_manual(
    values=colors_cases, 
    # breaks = cases_breaks,
    name=paste0('additional cases averted \n(per 100k children ', cur_age, ')')
  ) +
  theme_classic() +
  facet_grid(Annual_EIR~seasonality)
f_save_plot(gg3, paste0('contour_additional_cases_averted_smc_compared_mab_only_for_eir_seasonality',
                       cur_age, '_',
                       cur_cm * 100, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)


extreme_severe = max(compare_mab_smc$additional_severe_cases_averted_with_smc)
extreme_severe = ceiling(extreme_severe/100)*100
breaks_severe = seq(0, extreme_severe, length.out=8)
colors_severe = brewer.pal(n=length(breaks_severe), name='PRGn')
gg4 = ggplot(compare_mab_smc, aes(x=hh, y=max_efficacy, z=additional_severe_cases_averted_with_smc)) +
  geom_contour_filled(na.rm=TRUE, breaks=breaks_severe) +
  scale_fill_manual(
    values=colors_severe, 
    # breaks = cases_breaks,
    name=paste0('additional severe cases averted \n(per 100k children ', cur_age, ')')
  ) +
  theme_classic()+
  facet_grid(Annual_EIR~seasonality)
f_save_plot(gg4, paste0('contour_severe_cases_averted_mab_compared_smc_for_eir_seasonality',
                       cur_age, '_',
                       cur_cm * 100, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)




