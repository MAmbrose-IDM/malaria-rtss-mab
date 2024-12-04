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

exp_name_comp1 = 'mAb_sweep4_seeds1'
exp_name_comp2 = 'mAb_sweep4b_seeds1'
exp_name_comp3 = 'mAb_sweep4d_seeds1'
exp_name_comp4 = 'mAb_sweep4e_seeds1'

simout_dir=file.path(projectpath, 'simulation_output')
if (!dir.exists(paste0(simout_dir, '/_plots'))) dir.create(paste0(simout_dir, '/_plots'))
if (!dir.exists(paste0(simout_dir, '/_plots/pdf'))) dir.create(paste0(simout_dir, '/_plots/pdf'))



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Process simulation output and create line plots comparing cases 
#   averted with SMC versus mAbs across mAb params
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# exp_name = exp_name_comp1
sim_output = load_Age_monthly_Cases(simout_dir=simout_dir, exp_name=c(exp_name_comp1, exp_name_comp2, exp_name_comp3, exp_name_comp4), add_PE_perAge=TRUE,
                                    max_years=c(1, 3, 5, 8), keep_birth_month=FALSE, fname='All_Age_monthly_Cases.csv')

pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'vacc_info')
pe_df = pe_df[!is.na(pe_df$vacc_info),]

# subset simulations
cur_cm = 0.6
cur_age = 'U5'
cur_eirs = c(20)
cur_seasonalities = c('constant', 'moderate_unimodal', 'higher_unimodal')#, 'high_unimodal')
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   seasonality %in% cur_seasonalities,
                   Annual_EIR %in% cur_eirs,
                   cm_coverage == cur_cm,
)
pe_df_cur$seasonality = factor(pe_df_cur$seasonality, levels=cur_seasonalities)
# get rid of RTS,S (while leaving SMC)
pe_df_cur = pe_df_cur[-intersect(which(pe_df_cur$vacc_type %in% c('rtss', 'RTS,S')), which(pe_df_cur$smc_coverage <0.01)),]

# remove the highest max efficacy mAbs for plot
if(length(intersect(which(pe_df_cur$max_efficacy >98), which(pe_df_cur$vacc_type %in% c('mab', 'mAb'))))>0) pe_df_cur = pe_df_cur[-intersect(which(pe_df_cur$max_efficacy >98), which(pe_df_cur$vacc_type %in% c('mab', 'mAb'))),]

# create plots with colors matching efficacy-through-time plot
clist = get_intervention_colors(pe_df_cur, level_name = 'hh')
intervention_names = clist[[1]]
intervention_colors = clist[[2]]

gg1 = ggplot(pe_df_cur[pe_df_cur$smc_coverage < 0.01,], aes(x=hh, y=protective_efficacy, group=interaction(vacc_type, max_efficacy))) +
  geom_line(aes(linetype=as.factor(max_efficacy)))+
  geom_point(aes(col=as.factor(hh)), size=3)+
  scale_colour_manual(name='hh', breaks=intervention_names, values=intervention_colors) +
  scale_linetype_manual(name="max efficacy",values=c("80"=3,"90"=2, "95"=1)) +
  geom_hline(data=pe_df_cur[pe_df_cur$smc_coverage > 0.01,], aes(yintercept=protective_efficacy), col=rgb(0.5,0.5,0.5,0.5), linetype=1, size=1.2)+
  ylim(0,NA)+
  ylab(paste0('fraction of cases averted in children ', cur_age))+
  xlab('speed of protection decline (parameter hh)') +
  facet_wrap(facets='seasonality', nrow=1) + 
  theme_bw()
gg2 = ggplot(pe_df_cur[pe_df_cur$smc_coverage < 0.01,], aes(x=hh, y=cases_averted_per100000, group=interaction(vacc_type, max_efficacy))) +
  geom_line(aes(linetype=as.factor(max_efficacy)))+
  geom_point(aes(col=as.factor(hh)), size=3)+
  scale_colour_manual(name='hh', breaks=intervention_names, values=intervention_colors) +
  scale_linetype_manual(name="max efficacy",values=c("80"=3,"90"=2, "95"=1)) +
  geom_hline(data=pe_df_cur[pe_df_cur$smc_coverage > 0.01,], aes(yintercept=cases_averted_per100000), col=rgb(0.5,0.5,0.5,0.5), linetype=1, size=1.2)+
  ylim(0,NA)+
  ylab(paste0('cases averted per 100,000 children ', cur_age))+
  xlab('speed of protection decline (parameter hh)') +
  facet_wrap(facets='seasonality', nrow=1) + 
  theme_bw()
f_save_plot(gg1, paste0('compare_smc_mAb_frac_cases_', cur_age,'_averted_by_hh_eff_EIR', cur_eirs[1], '.png'),
            file.path(simout_dir, '_plots'), width = 7.5, height = 4.5, units = 'in', device_format = device_format)
f_save_plot(gg2, paste0('compare_smc_mAb_num_cases_', cur_age,'_averted_by_hh_eff_EIR', cur_eirs[1], '.png'),
            file.path(simout_dir, '_plots'), width = 7.5, height = 4.5, units = 'in', device_format = device_format)





# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Process simulation output and create heat plots comparing cases 
#   averted with SMC versus mAbs across mAb params
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
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
pe_df_compare$cases_averted_compared_to_smc = (pe_df_compare$cases_averted_per100000 - pe_df_compare$cases_averted_per100000_smc)/100
pe_df_compare$severe_cases_averted_compared_to_smc = (pe_df_compare$severe_cases_averted_per100000 - pe_df_compare$severe_cases_averted_per100000_smc)/100


# subset simulations
cur_cm = 0.6
cur_age = 'U5'
cur_eirs = c(5, 10, 30)
cur_seasonalities = c('constant', 'moderate_unimodal', 'high_unimodal', 'higher_unimodal')
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
extreme_cases = 600  # floor(extreme_cases/100)*100
breaks_cases = c(-Inf, seq(-1*extreme_cases, extreme_cases, length.out=8), Inf)
colors_cases = brewer.pal(n=length(breaks_cases), name='BrBG')
colors_cases[length(breaks_cases)/2] = '#EFF2DA'
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


# build plot for GR with example point for mAb product
example_product = data.frame(hh=10, max_efficacy=90, seasonality=factor(cur_seasonalities, levels=cur_seasonalities))
example_product = merge(example_product, data.frame(Annual_EIR=factor(cur_eirs, levels=sort(cur_eirs))), all=TRUE)
gg3b = ggplot() +
  geom_contour_filled(data=compare_mab_smc, aes(x=hh, y=max_efficacy, z=cases_averted_compared_to_smc), na.rm=TRUE, breaks=breaks_cases) +
  scale_fill_manual(
    values=colors_cases, 
    # breaks = cases_breaks,
    name='additional cases averted'
  ) +
  geom_point(data=example_product, aes(x=hh, y=max_efficacy), shape=20, size=3, color='red') +
  theme_classic() +
  facet_grid(Annual_EIR~seasonality)
f_save_plot(gg3b, paste0('contour_cases_averted_mab_compared_smc_for_eir_seasonality',
                        cur_age, '_',
                        cur_cm * 100, 'CM_GRexample'),
            file.path(simout_dir, '_plots'), width = 8, height = 5, units = 'in', device_format = device_format)



extreme_severe = min(abs(min(compare_mab_smc$severe_cases_averted_compared_to_smc)), max(compare_mab_smc$severe_cases_averted_compared_to_smc))
extreme_severe = 0.5  # floor(extreme_severe)
breaks_severe = c(-Inf, seq(-1*extreme_severe, extreme_severe, length.out=6), Inf)
colors_severe = brewer.pal(n=length(breaks_severe), name='BrBG')
colors_severe[length(breaks_severe)/2] = '#EFF2DA'
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
colors_cases[length(breaks_cases)/2] = '#EFF2DA'
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
colors_severe[length(breaks_severe)/2] = '#EFF2DA'
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




