# compare_sim_mAb_SMC.R
# August 2022
# create plots showing performance of mAbs relative to RTS,S at averting burden

library(lubridate)
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
library(viridis)
library(ggpattern)

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

exp_name_comp1 = 'sweep6_seeds1'
exp_name_comp2 = 'sweep7_seeds1'

simout_dir=file.path(projectpath, 'simulation_output')
if (!dir.exists(paste0(simout_dir, '/_plots'))) dir.create(paste0(simout_dir, '/_plots'))
if (!dir.exists(paste0(simout_dir, '/_plots/pdf'))) dir.create(paste0(simout_dir, '/_plots/pdf'))




###############################################################
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# cases and cases averted by age
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
###############################################################
# age labels for plots
age_label_values = seq(0, 8, length.out = 5)
age_labels = paste0(age_label_values, '-', (age_label_values + 1))

########################################################################################
# Process simulation output for cases and cases averted by age
########################################################################################
sim_output = load_Age_monthly_Cases(simout_dir=simout_dir, exp_name=c(exp_name_comp1, exp_name_comp2), add_PE_perAge=TRUE,
                                    max_years=c(5, 8), keep_birth_month=FALSE, fname='All_Age_monthly_Cases.csv')
# sim_output = load_Age_monthly_Cases(simout_dir=simout_dir, exp_name=exp_name_comp1, add_PE_perAge=TRUE,
#                                     max_years=c(2, 5, 8), keep_birth_month=FALSE, fname='All_Age_monthly_Cases.csv')
pe_df = sim_output[[4]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'smc_rtss_mab_age')
pe_df = pe_df[!is.na(pe_df$vacc_info),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))


# subset simulations
cur_cm = 0.6
cur_eirs = unique(pe_df$Annual_EIR)
cur_seasonalities = c('constant', 'moderate_unimodal', 'high_unimodal')
pe_df_cur = filter(pe_df,
                   seasonality %in% cur_seasonalities,
                   Annual_EIR %in% cur_eirs,
                   cm_coverage == cur_cm,
)
pe_df_cur$seasonality = factor(pe_df_cur$seasonality, levels=cur_seasonalities)
# remove some of the new vaccines
pe_df_cur = pe_df_cur[-which(pe_df_cur$hh<5),]

# only look at U2 and U5 targets and plot to age 6 (concerns about incidence-by-age and recalibration needs)
pe_df_cur = pe_df_cur[which(pe_df_cur$max_target_age %in% c(2, 5)),]
pe_df_cur = pe_df_cur[which(pe_df_cur$year < 8),]

# remove duplicate no-intervention simulations
pe_df_cur$max_target_age[which(as.character(pe_df_cur$smc_rtss_mab) == 'none')] = 2


########################################################################################
# mAbs, SMC, and RTSS on same age plots
########################################################################################
clist = get_intervention_colors(level_name = 'smc_rtss_mab')
intervention_names = clist[[1]]
intervention_colors = clist[[2]]
# plot clinical cases by age
gg1 = ggplot(pe_df_cur, aes(x=year, y=clinical_cases))+
  geom_point(aes(col=as.factor(smc_rtss_mab), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=1.5)+
  geom_line(aes(col=as.factor(smc_rtss_mab), linetype=as.factor(max_target_age), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=1)+
  geom_line(data=pe_df_cur[which(as.character(pe_df_cur$smc_rtss_mab) == 'none'),], aes(x=year, y=clinical_cases), color='black', size=2)+
  scale_alpha_discrete(range=c(1,0.3))+
  scale_colour_manual(name='smc_rtss_mab', breaks=intervention_names, values=intervention_colors) +
  theme_bw()+
  f_getCustomTheme() +
  geom_hline(yintercept=0)+
  ylab('clinical cases (per 100,000)') + 
  xlab('age of child') +
  scale_x_continuous(breaks=age_label_values, labels=age_labels) +
  facet_grid(cols=vars(seasonality), rows=vars(Annual_EIR))

f_save_plot(gg1, paste0('cases_by_age_mAb_RTSS_SMC.png'),
            file.path(simout_dir, '_plots'), width = 12, height = 8, units = 'in', device_format = device_format)



# remove the no-intervention rows
pe_df_cur = pe_df_cur[as.character(pe_df_cur$smc_rtss_mab) != 'none',]
# create plot of clinical cases averted
gg2 = ggplot(pe_df_cur, aes(x=year, y=cases_averted_per100000))+
  geom_point(aes(col=as.factor(smc_rtss_mab), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=2)+
  # geom_line(aes(col=as.factor(smc_rtss_mab), linetype=as.factor(max_target_age), group=smc_rtss_mab_age, alpha=1/(max_target_age*2)+0.75), size=1)+
  # geom_line(aes(col=as.factor(smc_rtss_mab), linetype=as.factor(max_target_age), group=smc_rtss_mab_age, alpha=factor(max_target_age, levels=sort(unique(max_target_age), decreasing=TRUE))), size=1)+
  geom_line(aes(col=as.factor(smc_rtss_mab), linetype=as.factor(max_target_age), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=1)+
  scale_alpha_discrete(range=c(1,0.3))+
  scale_colour_manual(name='smc_rtss_mab', breaks=intervention_names, values=intervention_colors) +
  theme_bw()+
  f_getCustomTheme() +
  geom_hline(yintercept=0)+
  ylab('cases averted (per 100,000)') + 
  xlab('age of child') +
  scale_x_continuous(breaks=age_label_values, labels=age_labels) +
  facet_grid(cols=vars(seasonality), rows=vars(Annual_EIR))
f_save_plot(gg2, paste0('averted_cases_by_age_mAb_RTSS_SMC.png'),
            file.path(simout_dir, '_plots'), width = 12, height = 8, units = 'in', device_format = device_format)


########################################################################################
# mAbs alone on age plots
########################################################################################

# subset simulations to mAbs and no-intervention only
pe_df_cur = filter(pe_df,
                   seasonality %in% cur_seasonalities,
                   Annual_EIR %in% cur_eirs,
                   cm_coverage == cur_cm,
                   smc_coverage < 0.01
)
# remove RTS,S
pe_df_cur = pe_df_cur[-intersect(which(pe_df_cur$vacc_type %in% c('rtss', 'RTS,S')), which(pe_df_cur$vacc_coverage >0)),]

# add factors for seasonality
pe_df_cur$seasonality = factor(pe_df_cur$seasonality, levels=cur_seasonalities)
# remove some of the new vaccines
pe_df_cur = pe_df_cur[-which(pe_df_cur$hh<5),]

# only look at U2 and U5 targets and plot to age 6 (concerns about incidence-by-age and recalibration needs)
pe_df_cur = pe_df_cur[which(pe_df_cur$max_target_age %in% c(2, 5)),]
pe_df_cur = pe_df_cur[which(pe_df_cur$year < 8),]

# remove duplicate no-intervention simulations
pe_df_cur$max_target_age[which(as.character(pe_df_cur$smc_rtss_mab) == 'none')] = 2

# plot clinical cases by age
gg3 = ggplot(pe_df_cur, aes(x=year, y=clinical_cases))+
  geom_point(aes(col=as.factor(smc_rtss_mab), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=1.5)+
  geom_line(aes(col=as.factor(smc_rtss_mab), linetype=as.factor(max_target_age), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=1)+
  geom_line(data=pe_df_cur[which(as.character(pe_df_cur$smc_rtss_mab) == 'none'),], aes(x=year, y=clinical_cases), color='black', size=2)+
  scale_alpha_discrete(range=c(1,0.3))+
  scale_colour_manual(name='smc_rtss_mab', breaks=intervention_names, values=intervention_colors) +
  theme_bw()+
  f_getCustomTheme() +
  geom_hline(yintercept=0)+
  ylab('clinical cases (per 1000)') + 
  xlab('age of child') +
  scale_x_continuous(breaks=age_label_values, labels=age_labels) +
  facet_grid(cols=vars(seasonality), rows=vars(Annual_EIR))

f_save_plot(gg3, paste0('cases_by_age_mAb.png'),
            file.path(simout_dir, '_plots'), width = 12, height = 8, units = 'in', device_format = device_format)



# remove the no-intervention rows
pe_df_cur = pe_df_cur[as.character(pe_df_cur$smc_rtss_mab) != 'none',]
# create plot of clinical cases averted
gg4 = ggplot(pe_df_cur, aes(x=year, y=cases_averted_per100000))+
  geom_point(aes(col=as.factor(smc_rtss_mab), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=2)+
  # geom_line(aes(col=as.factor(smc_rtss_mab), linetype=as.factor(max_target_age), group=smc_rtss_mab_age, alpha=1/(max_target_age*2)+0.75), size=1)+
  # geom_line(aes(col=as.factor(smc_rtss_mab), linetype=as.factor(max_target_age), group=smc_rtss_mab_age, alpha=factor(max_target_age, levels=sort(unique(max_target_age), decreasing=TRUE))), size=1)+
  geom_line(aes(col=as.factor(smc_rtss_mab), linetype=as.factor(max_target_age), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=1)+
  scale_alpha_discrete(range=c(1,0.3))+
  scale_colour_manual(name='smc_rtss_mab', breaks=intervention_names, values=intervention_colors) +
  theme_bw()+
  f_getCustomTheme() +
  geom_hline(yintercept=0)+
  ylab('cases averted (per 100,000)') + 
  xlab('age of child') +
  scale_x_continuous(breaks=age_label_values, labels=age_labels) +
  facet_grid(cols=vars(seasonality), rows=vars(Annual_EIR))
f_save_plot(gg4, paste0('averted_cases_by_age_mAb.png'),
            file.path(simout_dir, '_plots'), width = 12, height = 8, units = 'in', device_format = device_format)


######################################## 
# GR plots build sequence
######################################## 
gg3_GR_width = 7
gg3_GR_height = 5
# subset simulations to certain seasonalities and EIRs and to mAbs and no-intervention only
pe_df_cur = filter(pe_df,
                   seasonality %in% c('moderate_unimodal'),
                   Annual_EIR %in% c(30, 60, 120),
                   cm_coverage == cur_cm,
                   smc_coverage < 0.01
)
# remove RTS,S
pe_df_cur = pe_df_cur[-intersect(which(pe_df_cur$vacc_type %in% c('rtss', 'RTS,S')), which(pe_df_cur$vacc_coverage >0)),]
# remove some of the new vaccines
pe_df_cur = pe_df_cur[-which(pe_df_cur$hh<5),]

# only look at U2 and U5 targets and plot to age 6 (concerns about incidence-by-age and recalibration needs)
pe_df_cur = pe_df_cur[which(pe_df_cur$max_target_age %in% c(2, 5)),]
pe_df_cur = pe_df_cur[which(pe_df_cur$year < 8),]

# remove duplicate no-intervention simulations
pe_df_cur$max_target_age[which(as.character(pe_df_cur$smc_rtss_mab) == 'none')] = 2

# plot clinical cases by age
gg3_GR3 = ggplot(pe_df_cur, aes(x=year, y=clinical_cases))+
  geom_point(aes(col=as.factor(smc_rtss_mab), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=1.5)+
  geom_line(aes(col=as.factor(smc_rtss_mab), linetype=as.factor(max_target_age), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=1)+
  geom_line(data=pe_df_cur[which(as.character(pe_df_cur$smc_rtss_mab) == 'none'),], aes(x=year, y=clinical_cases), color='black', size=2)+
  scale_alpha_discrete(range=c(1,0.3))+
  scale_colour_manual(name='smc_rtss_mab', breaks=intervention_names, values=intervention_colors) +
  theme_bw()+
  f_getCustomTheme() +
  geom_hline(yintercept=0)+
  ylab('clinical cases (per 1000)') + 
  xlab('age of child') +
  scale_x_continuous(breaks=age_label_values, labels=age_labels) +
  facet_grid( rows=vars(Annual_EIR)) + 
  theme_classic()
# set legend spacing separately from main plot for consistency across builds
gg3_GR_legend = get_legend(gg3_GR3)
p3 = grid.arrange(gg3_GR3 + theme(legend.position="none"),
                   gg3_GR_legend, nrow=1,widths=c(2, 1))
f_save_plot(p3, paste0('cases_by_age_mAb_GR3.png'),
            file.path(simout_dir, '_plots'), width = gg3_GR_width, height = gg3_GR_height, units = 'in', device_format = device_format)


# remove older age group targeting
pe_df_cur = pe_df_cur[which(pe_df_cur$max_target_age %in% c(2)),]
gg3_GR2 = ggplot(pe_df_cur, aes(x=year, y=clinical_cases))+
  geom_point(aes(col=as.factor(smc_rtss_mab), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=1.5)+
  geom_line(aes(col=as.factor(smc_rtss_mab), linetype=as.factor(max_target_age), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=1)+
  geom_line(data=pe_df_cur[which(as.character(pe_df_cur$smc_rtss_mab) == 'none'),], aes(x=year, y=clinical_cases), color='black', size=2)+
  scale_alpha_discrete(range=c(1,0.3))+
  scale_colour_manual(name='smc_rtss_mab', breaks=intervention_names, values=intervention_colors) +
  theme_bw()+
  f_getCustomTheme() +
  geom_hline(yintercept=0)+
  ylab('clinical cases (per 1000)') + 
  xlab('age of child') +
  scale_x_continuous(breaks=age_label_values, labels=age_labels) +
  facet_grid(rows=vars(Annual_EIR)) + 
  theme_classic()
# set legend spacing separately from main plot for consistency across builds
p2 = grid.arrange(gg3_GR2 + theme(legend.position="none"),
                  gg3_GR_legend, nrow=1,widths=c(2, 1))
f_save_plot(p2, paste0('cases_by_age_mAb_GR2.png'),
            file.path(simout_dir, '_plots'), width = gg3_GR_width, height = gg3_GR_height, units = 'in', device_format = device_format)


# only show no-intervention scenario
pe_df_cur = pe_df_cur[which(pe_df_cur$vacc_coverage < 0.01),]
gg3_GR1 = ggplot(pe_df_cur, aes(x=year, y=clinical_cases))+
  geom_point(aes(col=as.factor(smc_rtss_mab), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=1.5)+
  geom_line(aes(col=as.factor(smc_rtss_mab), linetype=as.factor(max_target_age), group=smc_rtss_mab_age, alpha=as.factor(max_target_age)), size=1)+
  geom_line(data=pe_df_cur[which(as.character(pe_df_cur$smc_rtss_mab) == 'none'),], aes(x=year, y=clinical_cases), color='black', size=2)+
  scale_alpha_discrete(range=c(1,0.3))+
  scale_colour_manual(name='smc_rtss_mab', breaks=intervention_names, values=intervention_colors) +
  theme_bw()+
  f_getCustomTheme() +
  geom_hline(yintercept=0)+
  ylab('clinical cases (per 1000)') + 
  xlab('age of child') +
  scale_x_continuous(breaks=age_label_values, labels=age_labels) +
  facet_grid(rows=vars(Annual_EIR)) + 
  theme_classic()
# set legend spacing separately from main plot for consistency across builds
p1 = grid.arrange(gg3_GR1 + theme(legend.position="none"),
                  gg3_GR_legend, nrow=1,widths=c(2, 1))
f_save_plot(p1, paste0('cases_by_age_mAb_GR1.png'),
            file.path(simout_dir, '_plots'), width = gg3_GR_width, height = gg3_GR_height, units = 'in', device_format = device_format)




###############################################################
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# cases averted per 'dose' or 'dollar' of mAbs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
###############################################################



###########################################################################
# plot cases averted per 'dose' of mAbs
###########################################################################
coverage = 0.8
cost_per_visit = 0.75  # value should include everything aside from cost of drug itself (this value does not depend on dose, same for all ages)
price_per_g = 30
# how many healthcare interactions were in the U2 versus U5 vaccine versions?
# average number of doses received by each age-year. Starts at 4 months, but given seasonally and we simulate with cohorts
# for the U2 target, individuals between 4 and 23.99 months receive the vaccine seasonally. 
num_months_u2 = 24-4
num_months_u5 = 12*5-4

# given bodyweight scaling, how many units of product were used for U5 versus U2?
bodyweights = c(6.65, 9.1, 11.05, 12.95, 14.7, 15.8, 17.2, 20.05, 21.85, 22.35)
total_weight_dosed_per_child_u2 = sum(bodyweights[1:2]*c(8, 12))/12
total_weight_dosed_per_child_u5 = sum(bodyweights[1:5]*c(8,rep(12,4)))/12
# assume 20mg/kg administered to a given person (of weight Ykg) with price of $X per gram of mAb. Price for dose is X/1000*20*Ykg
total_dose_price_per_child_u2 = price_per_g / 1000 * 20 * total_weight_dosed_per_child_u2
total_dose_price_per_child_u5 = price_per_g / 1000 * 20 * total_weight_dosed_per_child_u5


pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'vacc_info')


cur_cm = 0.6
cur_eirs = unique(pe_df$Annual_EIR)
cur_age_numeric = 8
cur_age = paste0('U', cur_age_numeric)
cur_seasonalities = c('constant', 'moderate_unimodal', 'high_unimodal')
pe_df_cur = filter(pe_df,
                   seasonality %in% cur_seasonalities,
                   Annual_EIR %in% cur_eirs,
                   age_group == cur_age,
                   cm_coverage == cur_cm,
                   smc_coverage == 0,
                   vacc_coverage > 0.01,
                   vacc_type %in% c('mab', 'mAb')
)
pe_df_cur$seasonality = factor(pe_df_cur$seasonality, levels=cur_seasonalities)
# remove some of the new vaccines
pe_df_cur = pe_df_cur[-which(pe_df_cur$hh<5),]

# subset to U2 results for those receiving vaccine to age U2 and U5 results for those receving the vaccine to U5
# pe_df_cur = pe_df_cur[-intersect(which(pe_df_cur$age_group == 'U5'), which(pe_df_cur$max_target_age == '2')),]
# pe_df_cur = pe_df_cur[-intersect(which(pe_df_cur$age_group == 'U2'), which(pe_df_cur$max_target_age == '5')),]

# remove the U8 targeting - not plotted here
pe_df_cur = pe_df_cur[-which(pe_df_cur$max_target_age == '8'),]

# get total number of cases averted (from average averted in each year)
pe_df_cur$total_cases_averted_across_years = cur_age_numeric * pe_df_cur$cases_averted_per100000
pe_df_cur$total_cases_averted_across_years_per_child = pe_df_cur$total_cases_averted_across_years /100000
# pe_df_cur$total_cases_averted_across_years[pe_df_cur$age_group == 'U2'] = 5 * pe_df_cur$cases_averted_per100000[pe_df_cur$age_group == 'U2']
# pe_df_cur$total_cases_averted_across_years[pe_df_cur$age_group == 'U5'] = 5 * pe_df_cur$cases_averted_per100000[pe_df_cur$age_group == 'U5']

# cases averted per unit of the product
pe_df_cur$cases_averted_per_product_unit = NA
pe_df_cur$cases_averted_per_product_unit[pe_df_cur$max_target_age == '2'] = pe_df_cur$total_cases_averted_across_years[pe_df_cur$max_target_age == '2'] / (total_dose_price_per_child_u2 * pe_df_cur$pop * coverage)
pe_df_cur$cases_averted_per_product_unit[pe_df_cur$max_target_age == '5'] = pe_df_cur$total_cases_averted_across_years[pe_df_cur$max_target_age == '5'] / (total_dose_price_per_child_u5 * pe_df_cur$pop * coverage)

# cases averted per visit/healthcare interaction
pe_df_cur$cases_averted_per_visit = NA
pe_df_cur$cases_averted_per_visit[pe_df_cur$max_target_age == '2'] = pe_df_cur$total_cases_averted_across_years[pe_df_cur$max_target_age == '2'] / (num_months_u2/12 * pe_df_cur$pop * coverage)
pe_df_cur$cases_averted_per_visit[pe_df_cur$max_target_age == '5'] = pe_df_cur$total_cases_averted_across_years[pe_df_cur$max_target_age == '5'] / (num_months_u5/12 * pe_df_cur$pop * coverage)

# cases averted per dollar with current price estimates
pe_df_cur$cases_averted_per_dollar = NA
pe_df_cur$cases_averted_per_dollar[pe_df_cur$max_target_age == '2'] = pe_df_cur$total_cases_averted_across_years[pe_df_cur$max_target_age == '2'] / (num_months_u2/12*cost_per_visit * pe_df_cur$pop * coverage + total_dose_price_per_child_u2 * pe_df_cur$pop * coverage)
pe_df_cur$cases_averted_per_dollar[pe_df_cur$max_target_age == '5'] = pe_df_cur$total_cases_averted_across_years[pe_df_cur$max_target_age == '5'] / (num_months_u5/12*cost_per_visit * pe_df_cur$pop * coverage + total_dose_price_per_child_u5 * pe_df_cur$pop * coverage)


# ggplot(pe_df_cur, aes(x=vacc_info, y=cases_averted_per_product_unit, fill=as.factor(max_target_age)))+
#   geom_bar(stat='identity', position=position_dodge()) +
#   facet_grid(cols=vars(seasonality), rows=vars(Annual_EIR))
# 
# ggplot(pe_df_cur, aes(x=vacc_info, y=total_cases_averted_across_years, fill=as.factor(max_target_age)))+
#   geom_bar(stat='identity', position=position_dodge()) +
#   facet_grid(cols=vars(seasonality), rows=vars(Annual_EIR))

# get plotting colors
clist = get_intervention_colors(level_name = 'vacc_info')
intervention_names = clist[[1]]
intervention_colors = clist[[2]]
gg5a = ggplot(pe_df_cur, aes(x=vacc_info, y=total_cases_averted_across_years_per_child, pattern=as.factor(max_target_age), fill=as.factor(vacc_info)))+
  geom_col_pattern(position = position_dodge(preserve = "single"),
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.3,
                   pattern_spacing = 0.03,
                   pattern_key_scale_factor = 0.9
  ) +
  scale_fill_manual(name='vacc_info', breaks=intervention_names, values=intervention_colors) +
  scale_pattern_manual(values = c('2' = "none", '5' = "stripe")) +
  facet_grid(cols=vars(seasonality), rows=vars(Annual_EIR)) +
  theme_classic()

f_save_plot(gg5a, paste0('total_cases_averted_by_mAbs_across_years_byMaxAge.png'),
            file.path(simout_dir, '_plots'), width = 12, height = 8, units = 'in', device_format = device_format)


# ggplot(pe_df_cur, aes(x=vacc_info, y=cases_averted_per_product_unit, pattern=as.factor(max_target_age), fill=as.factor(vacc_info), color=as.factor(vacc_info)))+
#   geom_col_pattern(position = position_dodge(preserve = "single"),
#                  pattern_fill = "black",
#                  pattern_angle = 45,
#                  pattern_density = 0.3,
#                  pattern_spacing = 0.03,
#                  pattern_key_scale_factor = 0.9
#                  ) +
#   scale_fill_manual(name='vacc_info', breaks=intervention_names, values=intervention_colors) +
#   scale_pattern_manual(values = c('2' = "none", '5' = "stripe")) +
#   facet_grid(cols=vars(seasonality), rows=vars(Annual_EIR))
# 
# 
# ggplot(pe_df_cur, aes(x=vacc_info, y=cases_averted_per_visit, pattern=as.factor(max_target_age), fill=as.factor(vacc_info), color=as.factor(vacc_info)))+
#   geom_col_pattern(position = position_dodge(preserve = "single"),
#                    pattern_fill = "black",
#                    pattern_angle = 45,
#                    pattern_density = 0.3,
#                    pattern_spacing = 0.03,
#                    pattern_key_scale_factor = 0.9
#   ) +
#   scale_fill_manual(name='vacc_info', breaks=intervention_names, values=intervention_colors) +
#   scale_pattern_manual(values = c('2' = "none", '5' = "stripe")) +
#   facet_grid(cols=vars(seasonality), rows=vars(Annual_EIR))

gg5b = ggplot(pe_df_cur, aes(x=vacc_info, y=cases_averted_per_dollar, pattern=as.factor(max_target_age), fill=as.factor(vacc_info)))+
  geom_col_pattern(position = position_dodge(preserve = "single"),
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.3,
                   pattern_spacing = 0.03,
                   pattern_key_scale_factor = 0.9
  ) +
  scale_fill_manual(name='vacc_info', breaks=intervention_names, values=intervention_colors) +
  scale_pattern_manual(values = c('2' = "none", '5' = "stripe")) +
  facet_grid(cols=vars(seasonality), rows=vars(Annual_EIR)) +
  theme_classic()

f_save_plot(gg5b, paste0('cases_averted_by_mAbs_per_dollar.png'),
            file.path(simout_dir, '_plots'), width = 12, height = 8, units = 'in', device_format = device_format)



######################################## 
# GR plots
######################################## 
gg3_GR_width = 7 *2/5
gg3_GR_height = 5
# additional subsetting of simulations to certain seasonalities and EIRs
pe_df_cur = filter(pe_df_cur,
                   seasonality %in% c('moderate_unimodal'),
                   Annual_EIR %in% c(30, 60, 120),
)

gg5a_GR = ggplot(pe_df_cur, aes(x=vacc_info, y=total_cases_averted_across_years_per_child, pattern=as.factor(max_target_age), fill=as.factor(vacc_info)))+
  geom_col_pattern(position = position_dodge(preserve = "single"),
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.3,
                   pattern_spacing = 0.03,
                   pattern_key_scale_factor = 0.9
  ) +
  scale_fill_manual(name='vacc_info', breaks=intervention_names, values=intervention_colors) +
  scale_pattern_manual(values = c('2' = "none", '5' = "stripe")) +
  facet_grid(rows=vars(Annual_EIR)) +
  theme_classic() +
  theme(legend.position = 'none')
# # set legend spacing separately from main plot for consistency across builds
# gg5_GR_legend = get_legend(gg5a_GR)
# p5a = grid.arrange(gg5a_GR + theme(legend.position="none"),
#                   gg5_GR_legend, nrow=1,widths=c(2, 1))
f_save_plot(gg5a_GR, paste0('total_cases_averted_by_mAbs_across_years_byMaxAge_GR.png'),
            file.path(simout_dir, '_plots'), width = gg3_GR_width, height = gg3_GR_height, units = 'in', device_format = device_format)


gg5b_GR = ggplot(pe_df_cur, aes(x=vacc_info, y=cases_averted_per_dollar, pattern=as.factor(max_target_age), fill=as.factor(vacc_info)))+
  geom_col_pattern(position = position_dodge(preserve = "single"),
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.3,
                   pattern_spacing = 0.03,
                   pattern_key_scale_factor = 0.9
  ) +
  scale_fill_manual(name='vacc_info', breaks=intervention_names, values=intervention_colors) +
  scale_pattern_manual(values = c('2' = "none", '5' = "stripe")) +
  facet_grid(rows=vars(Annual_EIR)) +
  theme_classic() +
  theme(legend.position = 'none')
# # set legend spacing separately from main plot for consistency across builds
# p5b = grid.arrange(gg5b_GR + theme(legend.position="none"),
#                    gg5_GR_legend, nrow=1,widths=c(2, 1))
f_save_plot(gg5b_GR, paste0('cases_averted_by_mAbs_per_dollar_GR.png'),
            file.path(simout_dir, '_plots'), width = gg3_GR_width, height = gg3_GR_height, units = 'in', device_format = device_format)





###############################################################
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# cases averted across EIRs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
###############################################################


pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'vacc_info')
pe_df = pe_df[!is.na(pe_df$vacc_info),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))

# subset simulations
cur_cm = 0.6
cur_eirs = unique(pe_df$Annual_EIR)
cur_seasonalities = c('constant', 'moderate_unimodal', 'high_unimodal')
cur_max_target_age = 5
cur_age_numeric = 5
cur_age_group = paste0('U', cur_age_numeric)
pe_df_cur = filter(pe_df,
                   seasonality %in% cur_seasonalities,
                   Annual_EIR %in% cur_eirs,
                   cm_coverage == cur_cm,
                   smc_coverage == 0,
                   vacc_coverage > 0.01,
                   vacc_type %in% c('mab', 'mAb'),  #, 'rtss', 'RTS,S'),
                   max_target_age == cur_max_target_age,
                   age_group == cur_age_group
)
pe_df_cur$seasonality = factor(pe_df_cur$seasonality, levels=cur_seasonalities)
# remove some of the new vaccines
pe_df_cur = pe_df_cur[-which(pe_df_cur$hh<5),]


########################################################################################
# mAbs, SMC, and RTSS on same age plots
########################################################################################
pe_df_cur$total_cases_averted_per_child = pe_df_cur$cases_averted_per100000 * cur_age_numeric / 100000
# plot clinical cases by age
gg6 = ggplot(pe_df_cur, aes(x=Annual_EIR, y=total_cases_averted_per_child))+
  geom_point(aes(col=as.factor(vacc_info), group=vacc_info), size=1.5)+
  geom_line(aes(col=as.factor(vacc_info), group=vacc_info), size=1)+
  scale_colour_manual(name='smc_rtss_mab', breaks=intervention_names, values=intervention_colors) +
  theme_bw()+
  f_getCustomTheme() +
  geom_hline(yintercept=0)+
  ylab(paste0('total clinical cases averted over first  ', cur_age_numeric,' years per child')) + 
  xlab('annual EIR') +
  facet_wrap(.~seasonality)+
  theme_classic()

f_save_plot(gg6, paste0('cases_averted_by_mAbs_across_EIRs.png'),
            file.path(simout_dir, '_plots'), width = 10, height = 4, units = 'in', device_format = device_format)



