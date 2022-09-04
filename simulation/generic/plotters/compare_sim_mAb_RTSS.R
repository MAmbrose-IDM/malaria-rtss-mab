# compare_sim_mAb_RTSS.R
# August 2022
# create plots showing performance of mAbs relative to RTS,S at averting burden

library(lubridate)
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(viridis)

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
# exp_name = 'sweep4_seeds1'

simout_dir=file.path(projectpath, 'simulation_output')
if (!dir.exists(paste0(simout_dir, '/_plots'))) dir.create(paste0(simout_dir, '/_plots'))
if (!dir.exists(paste0(simout_dir, '/_plots/pdf'))) dir.create(paste0(simout_dir, '/_plots/pdf'))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Process simulation output and create line plots comparing cases 
#   averted with RTS,S versus mAbs across mAb params
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# exp_name = exp_name_comp1
sim_output = load_Age_monthly_Cases(simout_dir=simout_dir, exp_name=c(exp_name_comp1, exp_name_comp2), add_PE_perAge=TRUE,
                                    max_years=c(1, 3, 5, 8), keep_birth_month=FALSE, fname='All_Age_monthly_Cases.csv')

pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'vacc_info')
pe_df = pe_df[!is.na(pe_df$vacc_info),]
pe_df$hh[pe_df$vacc_type %in% c('rtss', 'RTS,S')]=NA

# subset simulations
cur_cm = 0.6
cur_age = 'U8'
cur_smcs = c(0)
cur_eirs = c(30)
cur_seasonalities = c('constant', 'moderate_unimodal', 'high_unimodal', 'higher_unimodal')
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   seasonality %in% cur_seasonalities,
                   Annual_EIR %in% cur_eirs,
                   cm_coverage == cur_cm,
                   vacc_coverage > 0,
                   smc_coverage %in% cur_smcs,
)
pe_df_cur$seasonality = factor(pe_df_cur$seasonality, levels=cur_seasonalities)


# create plots with colors matching efficacy-through-time plot
hhs = sort(unique(pe_df_cur$hh))
hh_colors = viridis(length(hhs))
gg1 = ggplot(pe_df_cur, aes(x=hh, y=protective_efficacy, group=interaction(vacc_type, max_efficacy))) +
  geom_line(aes(linetype=as.factor(max_efficacy)))+
  geom_point(aes(col=as.factor(hh)), size=3)+
  scale_colour_manual(name='hh', breaks=hhs, values=hh_colors) +
  scale_linetype_manual(name="max efficacy",values=c("80"=3,"90"=2, "95"=1)) +
  geom_hline(data=pe_df_cur[pe_df_cur$vacc_type %in% c('rtss', 'RTS,S'),], aes(yintercept=protective_efficacy), col=rgb(0.5,0.5,0.5,0.5), linetype=1, size=1.2)+
  ylim(0,NA)+
  ylab(paste0('fraction of cases averted in children ', cur_age))+
  xlab('speed of protection decline (parameter hh)') +
  facet_wrap(facets='seasonality', nrow=1) + 
  theme_bw()
gg2 = ggplot(pe_df_cur, aes(x=hh, y=cases_averted_per100000, group=interaction(vacc_type, max_efficacy))) +
  geom_line(aes(linetype=as.factor(max_efficacy)))+
  geom_point(aes(col=as.factor(hh)), size=3)+
  scale_colour_manual(name='hh', breaks=hhs, values=hh_colors) +
  scale_linetype_manual(name="max efficacy",values=c("80"=3,"90"=2, "95"=1)) +
  geom_hline(data=pe_df_cur[pe_df_cur$vacc_type %in% c('rtss', 'RTS,S'),], aes(yintercept=cases_averted_per100000), col=rgb(0.5,0.5,0.5,0.5), linetype=1, size=1.2)+
  ylim(0,NA)+
  ylab(paste0('cases averted per 100,000 children ', cur_age))+
  xlab('speed of protection decline (parameter hh)') +
  facet_wrap(facets='seasonality', nrow=1) + 
  theme_bw()
f_save_plot(gg1, paste0('compare_rtss_mAb_frac_cases_', cur_age,'_averted_by_hh_eff.png'),
            file.path(simout_dir, '_plots'), width = 6, height = 4, units = 'in', device_format = device_format)
f_save_plot(gg2, paste0('compare_rtss_mAb_num_cases_', cur_age,'_averted_by_hh_eff.png'),
            file.path(simout_dir, '_plots'), width = 6, height = 4, units = 'in', device_format = device_format)





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Process simulation output and create bar plots comparing cases 
#   averted with RTS,S versus mAbs across EIRs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# exp_name = exp_name_comp2
simout_dir=file.path(projectpath, 'simulation_output')
sim_output = load_Age_monthly_Cases(simout_dir=simout_dir, exp_name=exp_name, add_PE_perAge=TRUE,
                                    max_years=c(5, 8), keep_birth_month=FALSE, fname='All_Age_monthly_Cases.csv')
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'vacc_info')
pe_df = pe_df[!is.na(pe_df$vacc_info),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))

# subset simulations
cur_cm = 0.6
cur_age = 'U5'
cur_smcs = c(0)
cur_eirs = c(5, 10, 30)
cur_seasonalities = c('constant', 'moderate_unimodal', 'high_unimodal')
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   seasonality %in% cur_seasonalities,
                   Annual_EIR %in% cur_eirs,
                   cm_coverage == cur_cm,
                   vacc_coverage > 0,
                   smc_coverage %in% cur_smcs,
)
pe_df_cur$seasonality = factor(pe_df_cur$seasonality, levels=cur_seasonalities)

# percent reduction for simple comparison
pe_df_cur$vacc_percent_reduction = pe_df_cur$vacc_relative_burden * -100
pe_df_cur$vacc_percent_reduction_severe = pe_df_cur$vacc_relative_burden_severe * -100

# clinical cases
gg1 = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'vacc_percent_reduction', bargroup_var = 'vacc_info',
                            fillvar = 'vacc_info', facet1 = 'seasonality', facet2 = NULL,
                            SAVE = FALSE, ylab = 'percent reduction in clinical cases \nfrom adding vaccine', scales='fixed', bar_width=0.85) 
# severe cases
gg2 = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'vacc_percent_reduction_severe', bargroup_var = 'vacc_info',
                            fillvar = 'vacc_info', facet1 = 'seasonality', facet2 = NULL,
                            SAVE = FALSE, ylab = 'percent reduction in severe cases \nfrom adding vaccine', scales='fixed', bar_width=0.85)

# combine plots into panel
gg_legend = get_legend(gg1)
gg1 = gg1  + theme(legend.position = 'None')
gg2 = gg2  + theme(legend.position = 'None')
gg = plot_grid(gg1,gg2, align='hv', nrow=2)
gg = plot_grid(gg, gg_legend, rel_widths = c(1, 0.25))
f_save_plot(gg, paste0('burden_reduction_by_vacc_types_for_eir_seasonality_0SMC',
                       cur_age, '_',
                       cur_cm * 100, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)






# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Process simulation output and create bar plots comparing 
#   performance of mAbs relative to RTSS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# exp_name = exp_name_comp2
simout_dir=file.path(projectpath, 'simulation_output')
sim_output = load_Age_monthly_Cases(simout_dir=simout_dir, exp_name=exp_name, add_PE_perAge=TRUE,
                                    max_years=c(5, 8), keep_birth_month=FALSE, fname='All_Age_monthly_Cases.csv')
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'vacc_info')
pe_df = pe_df[!is.na(pe_df$vacc_info),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))

# compare mAbs against matched RTSS impact. Start by getting 'reference' of RTSS
# subset to RTSS
pe_df_rtss = pe_df[pe_df$vacc_type %in% c('rtss', 'RTS,S'),]
pe_df_rtss = pe_df_rtss[pe_df_rtss$vacc_coverage > 0,]
match_cols = c('Annual_EIR', 'seasonality', 'age_group', 'cm_coverage', 'smc_coverage', 'frac_high_access',
              'cm_target_group', 'smc_target_group', 'vacc_target_group')
keep_cols = c('cases_averted_per100000', 'severe_cases_averted_per100000')
df_rtss = pe_df_rtss %>%
  dplyr::select(c(match_cols, keep_cols)) %>%
  rename(cases_averted_per100000_rtss = cases_averted_per100000,
         severe_cases_averted_per100000_rtss = severe_cases_averted_per100000)
# subset to mAbs
pe_df_mab = pe_df[pe_df$vacc_type %in% c('mab', 'mAb'),]
pe_df_mab = pe_df_mab[pe_df_mab$vacc_coverage > 0,]


# merge the reference values into the data frame
pe_df_compare <- pe_df_mab %>% left_join(df_rtss)
pe_df_compare$cases_averted_compared_to_rtss = pe_df_compare$cases_averted_per100000 - pe_df_compare$cases_averted_per100000_rtss
pe_df_compare$severe_cases_averted_compared_to_rtss = pe_df_compare$severe_cases_averted_per100000 - pe_df_compare$severe_cases_averted_per100000_rtss


# subset simulations
cur_cm = 0.6
cur_age = 'U5'
cur_smcs = c(0)
cur_eirs = c(5, 10, 30)
cur_seasonalities = c('constant', 'moderate_unimodal', 'high_unimodal')
pe_df_cur = filter(pe_df_compare,
                   age_group == cur_age,
                   seasonality %in% cur_seasonalities,
                   Annual_EIR %in% cur_eirs,
                   cm_coverage == cur_cm,
                   vacc_coverage > 0,
                   smc_coverage %in% cur_smcs,
)
pe_df_cur$seasonality = factor(pe_df_cur$seasonality, levels=cur_seasonalities)

# clinical cases
gg1 = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'cases_averted_compared_to_rtss', bargroup_var = 'vacc_info',
                            fillvar = 'vacc_info', facet1 = 'seasonality', facet2 = NULL,
                            SAVE = FALSE, ylab = 'cases averted (per 100k U5) \nwith mAb compared to RTS,S', scales='fixed', bar_width=0.85) 
# severe cases
gg2 = plot_grouped_barplots(dat = pe_df_cur, xvar = 'Annual_EIR', yvar = 'severe_cases_averted_compared_to_rtss', bargroup_var = 'vacc_info',
                            fillvar = 'vacc_info', facet1 = 'seasonality', facet2 = NULL,
                            SAVE = FALSE, ylab = 'severe cases averted (per 100k U5) \nwith mAb compared to RTS,S', scales='fixed', bar_width=0.85)

# combine plots into panel
gg_legend = get_legend(gg1)
gg1 = gg1  + theme(legend.position = 'None')
gg2 = gg2  + theme(legend.position = 'None')
gg = plot_grid(gg1,gg2, align='hv', nrow=2)
gg = plot_grid(gg, gg_legend, rel_widths = c(1, 0.25))
f_save_plot(gg, paste0('cases_averted_mab_compared_rtss_for_eir_seasonality_0SMC',
                       cur_age, '_',
                       cur_cm * 100, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)


########################################################
# heat/contour maps
########################################################
# Cases averted
# create grid showing where mAbs outcompete RTSS
#  which set of parameters do better and worse than RTSS?
pe_df_cur$mAb_better_clinical = pe_df_cur$cases_averted_compared_to_rtss > 0
pe_df_cur$mAb_better_severe = pe_df_cur$severe_cases_averted_compared_to_rtss > 0
keep_cols = c('Annual_EIR', 'seasonality', 'vacc_char', 'vacc_coverage', 'vacc_target_group', 'hh', 'max_efficacy', 'vacc_info', 'mAb_better_clinical', 'mAb_better_severe', 'cases_averted_compared_to_rtss', 'severe_cases_averted_compared_to_rtss')
compare_mab_rtss = pe_df_cur %>% dplyr::select(keep_cols)

extreme_cases = min(abs(min(compare_mab_rtss$cases_averted_compared_to_rtss)), max(compare_mab_rtss$cases_averted_compared_to_rtss))
extreme_cases = 30000#floor(extreme_cases/10000)*10000
breaks_cases = c(-Inf, seq(-1*extreme_cases, extreme_cases, length.out=6), Inf)
colors_cases = brewer.pal(n=length(breaks_cases), name='BrBG')
gg3 = ggplot(compare_mab_rtss, aes(x=hh, y=max_efficacy, z=cases_averted_compared_to_rtss)) +
  geom_contour_filled(na.rm=TRUE, breaks=breaks_cases) +
  scale_fill_manual(
    values=colors_cases, 
    # breaks = cases_breaks,
    name='additional cases averted'
  ) +
  theme_classic() +
  facet_grid(Annual_EIR~seasonality)
f_save_plot(gg3, paste0('contour_cases_averted_mab_compared_rtss_for_eir_seasonality',
                        cur_age, '_',
                        cur_cm * 100, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)


extreme_severe = min(abs(min(compare_mab_rtss$severe_cases_averted_compared_to_rtss)), max(compare_mab_rtss$severe_cases_averted_compared_to_rtss))
extreme_severe = floor(extreme_severe/10)*10
breaks_severe = c(-Inf, seq(-1*extreme_severe, extreme_severe, length.out=6), Inf)
colors_severe = brewer.pal(n=length(breaks_severe), name='BrBG')
gg4 = ggplot(compare_mab_rtss, aes(x=hh, y=max_efficacy, z=severe_cases_averted_compared_to_rtss)) +
  geom_contour_filled(na.rm=TRUE, breaks=breaks_severe) +
  scale_fill_manual(
    values=colors_severe, 
    # breaks = cases_breaks,
    name='additional severe cases averted'
  ) +
  theme_classic()+
  facet_grid(Annual_EIR~seasonality)
f_save_plot(gg4, paste0('contour_severe_cases_averted_mab_compared_rtss_for_eir_seasonality',
                        cur_age, '_',
                        cur_cm * 100, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)



# Cases averted relative to RTSS
pe_df_cur$mAb_averted_rel_rtss = (pe_df_cur$cases_averted_per100000 - pe_df_cur$cases_averted_per100000_rtss) / pe_df_cur$cases_averted_per100000_rtss
pe_df_cur$mAb_severe_averted_rel_rtss = (pe_df_cur$severe_cases_averted_per100000 - pe_df_cur$severe_cases_averted_per100000_rtss) / pe_df_cur$severe_cases_averted_per100000_rtss
keep_cols = c('Annual_EIR', 'seasonality', 'vacc_char', 'vacc_coverage', 'vacc_target_group', 'hh', 'max_efficacy', 'vacc_info', 'mAb_averted_rel_rtss', 'mAb_severe_averted_rel_rtss')
compare_mab_rtss = pe_df_cur %>% dplyr::select(keep_cols)

extreme_cases = min(abs(min(compare_mab_rtss$mAb_averted_rel_rtss)), max(compare_mab_rtss$mAb_averted_rel_rtss))
extreme_cases = 0.3#floor(extreme_cases*10)/10
breaks_cases = c(-3, seq(-1*extreme_cases, extreme_cases, length.out=6), 3)
colors_cases = brewer.pal(n=length(breaks_cases), name='BrBG')
gg5 = ggplot(compare_mab_rtss, aes(x=hh, y=max_efficacy, z=mAb_averted_rel_rtss)) +
  geom_contour_filled(na.rm=TRUE, breaks=breaks_cases) +
  scale_fill_manual(
    values=colors_cases, 
    name='relative cases averted'
  ) +
  theme_classic() +
  facet_grid(Annual_EIR~seasonality)
f_save_plot(gg5, paste0('contour_relative_cases_averted_mab_compared_rtss_for_eir_seasonality',
                        cur_age, '_',
                        cur_cm * 100, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)


extreme_severe = min(abs(min(compare_mab_rtss$mAb_severe_averted_rel_rtss)), max(compare_mab_rtss$mAb_severe_averted_rel_rtss))
extreme_severe = floor(extreme_severe*10)/10
breaks_severe = c(-3, seq(-1*extreme_severe, extreme_severe, length.out=6), 3)
colors_severe = brewer.pal(n=length(breaks_severe), name='BrBG')
gg6 = ggplot(compare_mab_rtss, aes(x=hh, y=max_efficacy, z=mAb_severe_averted_rel_rtss)) +
  geom_contour_filled(na.rm=TRUE, breaks=breaks_severe) +
  scale_fill_manual(
    values=colors_severe, 
    name='relative severe cases averted'
  ) +
  theme_classic()+
  facet_grid(Annual_EIR~seasonality)
f_save_plot(gg6, paste0('contour_relative_severe_cases_averted_mab_compared_rtss_for_eir_seasonality',
                        cur_age, '_',
                        cur_cm * 100, 'CM'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)


















