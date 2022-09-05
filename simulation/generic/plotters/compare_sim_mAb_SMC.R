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
exp_name_comp2 = 'sweep4b_seeds1'
exp_name_comp3 = 'sweep4c_seeds1'

simout_dir=file.path(projectpath, 'simulation_output')
if (!dir.exists(paste0(simout_dir, '/_plots'))) dir.create(paste0(simout_dir, '/_plots'))
if (!dir.exists(paste0(simout_dir, '/_plots/pdf'))) dir.create(paste0(simout_dir, '/_plots/pdf'))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Process simulation output and create heat plots comparing cases 
#   averted with SMC versus mAbs across mAb params
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
sim_output = load_Age_monthly_Cases(simout_dir=simout_dir, exp_name=c(exp_name_comp1, exp_name_comp2, exp_name_comp3), add_PE_perAge=TRUE,
                                    max_years=c(5, 8), keep_birth_month=FALSE, fname='All_Age_monthly_Cases.csv')
# sim_output = load_Age_monthly_Cases(simout_dir=simout_dir, exp_name=exp_name_comp1, add_PE_perAge=TRUE,
#                                     max_years=c(5, 8), keep_birth_month=FALSE, fname='All_Age_monthly_Cases.csv')
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'vacc_info')
pe_df = pe_df[!is.na(pe_df$vacc_info),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))

# compare mAbs against matched SMC impact. Start by getting 'reference' of SMC
# subset to SMC
pe_df_smc = pe_df[pe_df$smc_coverage > 0,]
match_cols = c('Annual_EIR', 'seasonality', 'age_group', 'cm_coverage', 'frac_high_access', 'cm_target_group')
keep_cols = c('cases_averted_per100000', 'severe_cases_averted_per100000')
df_smc = pe_df_smc %>%
  dplyr::select(c(match_cols, keep_cols)) %>%
  rename(cases_averted_per100000_smc = cases_averted_per100000,
         severe_cases_averted_per100000_smc = severe_cases_averted_per100000)
# subset to mAbs
pe_df_mab = pe_df[pe_df$vacc_type %in% c('mab', 'mAb'),]
pe_df_mab = pe_df_mab[pe_df_mab$vacc_coverage > 0,]

# merge the reference values into the data frame
pe_df_compare <- pe_df_mab %>% left_join(df_smc)
pe_df_compare$cases_averted_compared_to_smc = pe_df_compare$cases_averted_per100000 - pe_df_compare$cases_averted_per100000_smc
pe_df_compare$severe_cases_averted_compared_to_smc = pe_df_compare$severe_cases_averted_per100000 - pe_df_compare$severe_cases_averted_per100000_smc


# subset simulations
cur_cm = 0.6
cur_age = 'U5'
cur_eirs = c(5, 10, 30)
cur_seasonalities = c('constant', 'moderate_unimodal', 'higher_unimodal')
pe_df_cur = filter(pe_df_compare,
                   age_group == cur_age,
                   seasonality %in% cur_seasonalities,
                   Annual_EIR %in% cur_eirs,
                   cm_coverage == cur_cm,
                   vacc_coverage > 0,
)
pe_df_cur$seasonality = factor(pe_df_cur$seasonality, levels=cur_seasonalities)

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



# create grid showing where mAbs outcompete SMC
#  which set of parameters do better and worse than SMC?
pe_df_cur$mAb_better_clinical = pe_df_cur$cases_averted_compared_to_smc > 0
pe_df_cur$mAb_better_severe = pe_df_cur$severe_cases_averted_compared_to_smc > 0
keep_cols = c('Annual_EIR', 'seasonality', 'vacc_char', 'vacc_coverage', 'vacc_target_group', 'smc_target_group', 'hh', 'max_efficacy', 'vacc_info', 'mAb_better_clinical', 'mAb_better_severe', 'cases_averted_compared_to_smc', 'severe_cases_averted_compared_to_smc')
compare_mab_smc = pe_df_cur %>% dplyr::select(keep_cols)

extreme_cases = min(abs(min(compare_mab_smc$cases_averted_compared_to_smc)), max(compare_mab_smc$cases_averted_compared_to_smc))
extreme_cases = floor(extreme_cases/10000)*10000
breaks_cases = c(-Inf, seq(-1*extreme_cases, extreme_cases, length.out=6), Inf)
colors_cases = brewer.pal(n=length(breaks_cases), name='BrBG')
gg3 = ggplot(compare_mab_smc, aes(x=hh, y=max_efficacy, z=cases_averted_compared_to_smc)) +
  geom_contour_filled(na.rm=TRUE, breaks=breaks_cases) +
  scale_fill_manual(
    values=colors_cases, 
    # breaks = cases_breaks,
    name='additional cases averted'
  ) +
  theme_classic() +
  facet_grid(Annual_EIR~seasonality)
f_save_plot(gg3, paste0('contour_cases_averted_mab_compared_smc_for_eir_seasonality',
                       cur_age, '_',
                       cur_cm * 100, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)


extreme_severe = min(abs(min(compare_mab_smc$severe_cases_averted_compared_to_smc)), max(compare_mab_smc$severe_cases_averted_compared_to_smc))
extreme_severe = floor(extreme_severe/10)*10
breaks_severe = c(-Inf, seq(-1*extreme_severe, extreme_severe, length.out=6), Inf)
colors_severe = brewer.pal(n=length(breaks_severe), name='BrBG')
gg4 = ggplot(compare_mab_smc, aes(x=hh, y=max_efficacy, z=severe_cases_averted_compared_to_smc)) +
  geom_contour_filled(na.rm=TRUE, breaks=breaks_severe) +
  scale_fill_manual(
    values=colors_severe, 
    # breaks = cases_breaks,
    name='additional severe cases averted'
  ) +
  theme_classic()+
  facet_grid(Annual_EIR~seasonality)
f_save_plot(gg4, paste0('contour_severe_cases_averted_mab_compared_smc_for_eir_seasonality',
                       cur_age, '_',
                       cur_cm * 100, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)






# Cases averted relative to SMC
pe_df_cur$mAb_averted_rel_smc = (pe_df_cur$cases_averted_per100000 - pe_df_cur$cases_averted_per100000_smc) / pe_df_cur$cases_averted_per100000_smc
pe_df_cur$mAb_severe_averted_rel_smc = (pe_df_cur$severe_cases_averted_per100000 - pe_df_cur$severe_cases_averted_per100000_smc) / pe_df_cur$severe_cases_averted_per100000_smc
keep_cols = c('Annual_EIR', 'seasonality', 'vacc_char', 'vacc_coverage', 'vacc_target_group', 'hh', 'max_efficacy', 'vacc_info', 'mAb_averted_rel_smc', 'mAb_severe_averted_rel_smc')
compare_mab_smc = pe_df_cur %>% dplyr::select(keep_cols)

extreme_cases = min(abs(min(compare_mab_smc$mAb_averted_rel_smc)), max(compare_mab_smc$mAb_averted_rel_smc))
extreme_cases = floor(extreme_cases*10)/10
breaks_cases = c(-3, seq(-1*extreme_cases, extreme_cases, length.out=6), 3)
colors_cases = brewer.pal(n=length(breaks_cases), name='BrBG')
gg5 = ggplot(compare_mab_smc, aes(x=hh, y=max_efficacy, z=mAb_averted_rel_smc)) +
  geom_contour_filled(na.rm=TRUE, breaks=breaks_cases) +
  scale_fill_manual(
    values=colors_cases, 
    name='relative cases averted'
  ) +
  theme_classic() +
  facet_grid(Annual_EIR~seasonality)
f_save_plot(gg5, paste0('contour_relative_cases_averted_mab_compared_smc_for_eir_seasonality',
                        cur_age, '_',
                        cur_cm * 100, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)


extreme_severe = min(abs(min(compare_mab_smc$mAb_severe_averted_rel_smc)), max(compare_mab_smc$mAb_severe_averted_rel_smc))
extreme_severe = floor(extreme_severe*10)/10
breaks_severe = c(-3, seq(-1*extreme_severe, extreme_severe, length.out=6), 3)
colors_severe = brewer.pal(n=length(breaks_severe), name='BrBG')
gg6 = ggplot(compare_mab_smc, aes(x=hh, y=max_efficacy, z=mAb_severe_averted_rel_smc)) +
  geom_contour_filled(na.rm=TRUE, breaks=breaks_severe) +
  scale_fill_manual(
    values=colors_severe, 
    name='relative severe cases averted'
  ) +
  theme_classic()+
  facet_grid(Annual_EIR~seasonality)
f_save_plot(gg6, paste0('contour_relative_severe_cases_averted_mab_compared_smc_for_eir_seasonality',
                        cur_age, '_',
                        cur_cm * 100, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)




