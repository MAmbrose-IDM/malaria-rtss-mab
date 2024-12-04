# RTSS_validation_Phase3.R
# July 2022
# create plots of RTS,S simluation output to compare against Figure S2.4 from Penny et al. 

library(lubridate)
library(dplyr)
library(ggplot2)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Setup filepaths and parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setwd("C:/Users/moniqueam/Documents/malaria-rtss-mab")

source(file.path("simulation", "load_paths.R"))
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
device_format = c('pdf', 'png')
paths = get_project_paths()
datapath = paths[1]
projectpath = paths[2]

booster_string = 'boost'  #'no.boost' or 'boost'
if(booster_string == 'boost'){
  exp_name1 = 'TEST_validation_phase3_wBooster_firstParamSet'; bb1='_bb80' # NA
  exp_name2 = 'validation_phase3_wBooster'; bb2 = '_bb60'  # (was previously mislabeled _bb80
  exp_name3 = 'validation_Phase3_wBooster_p1_cleaned_4seeds'; bb3=''
  exp_name4 = 'validation_Phase3_wBooster_p2_cleaned_4seeds'; bb4=''
  exp_names = c(exp_name1, exp_name2, exp_name3, exp_name4)
  bbs = c(bb1, bb2, bb3, bb4)
  # exp_names = c(exp_name2, exp_name3, exp_name4)
  # bbs = c(bb2, bb3, bb4)
} else {
  exp_name1 = 'TEST_validation_phase3_noBooster_firstParamSet'; bb2=''  # NA
  exp_name2 = 'validation_phase3_noBooster'; bb1=''
  exp_names = c(exp_name1, exp_name2)
  bbs = c(bb1, bb2)
}

# reference_filepath = file.path(datapath, 'rtss_phase3/kintampo_trial_summary_3month.csv')
reference_filepath = file.path(datapath, 'rtss_phase3/RTSS_32month_allCases.csv')


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Process simulation output
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# combine separately-run experiments with different vaccine parameter sets
for (ee in 1:length(exp_names)){
  exp_filepath = file.path(projectpath, 'simulation_output', exp_names[ee])
  cases_filepath = file.path(exp_filepath, 'All_Age_monthly_Cases.csv')
  cases_df_cur = read.csv(cases_filepath)
  cases_df_cur$vacc_char = paste0(cases_df_cur$vacc_char, bbs[ee])
  if (ee==1){
    cases_df = cases_df_cur
  } else{
    cases_df = merge(cases_df, cases_df_cur, all=TRUE)
  }
}

# subset to time after vaccination
vaccine_date = as.Date('2021-01-01')  # date vaccination given in simulations
cases_df$date = as.Date(cases_df$date)
cases_df = cases_df[cases_df$date >= vaccine_date,]


# get months since vaccination
elapsed_months <- function(start_date, end_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}
cases_df$elapsed_months = elapsed_months(vaccine_date, cases_df$date)

# indicate time group (three months per group)
cases_df$time.group = floor(cases_df$elapsed_months / 3) + 1

# get total within each three-month period
cases_df_grouped = cases_df %>%
  dplyr::group_by(Scenario_id, Run_Number, time.group, Annual.EIR, vacc_coverage, vacc_char, seasonality) %>%
  dplyr::summarise(total_cases = sum(New.Clinical.Cases),
                   total_severe_cases = sum(New.Severe.Cases),
                   average_pfpr = mean(PfHRP2.Prevalence),
                   average_pop = mean(Statistical.Population),
                   middle_month = mean(date)
  ) %>%
  dplyr::ungroup()

# set the simulations without vaccination to have the same vacc_char=NA
cases_df_grouped$vacc_char[cases_df_grouped$vacc_coverage <0.001] = NA

# get averages across runs
cases_df_ave = cases_df_grouped %>%
  dplyr::group_by(time.group, Annual.EIR, vacc_coverage, vacc_char, seasonality) %>%  # Scenario_id
  dplyr::summarise(clinical_cases = mean(total_cases),
                   severe_cases = mean(total_severe_cases),
                   pfpr = mean(average_pfpr),
                   pop = mean(average_pop),
                   middle_month = mean(middle_month)
  ) %>%
  dplyr::ungroup()

# calculate the PE
# merge pairs of simulations with versus without RTSS
with_vacc = cases_df_ave[cases_df_ave$vacc_coverage > 0, which(colnames(cases_df_ave) %in% c('time.group', 'Annual.EIR', 'seasonality', 'clinical_cases', 'middle_month', 'vacc_char'))]
colnames(with_vacc)[which(colnames(with_vacc) == 'clinical_cases')] = 'clinical_cases_vacc'
without_vacc = cases_df_ave[cases_df_ave$vacc_coverage < 0.00001, which(colnames(cases_df_ave) %in% c('time.group', 'Annual.EIR', 'seasonality', 'clinical_cases', 'middle_month'))]
colnames(without_vacc)[which(colnames(without_vacc) == 'clinical_cases')] = 'clinical_cases_no_vacc'

comparison_df = merge(with_vacc, without_vacc, by = c('time.group', 'Annual.EIR', 'middle_month', 'seasonality'))
comparison_df$PE_sim = 1 - comparison_df$clinical_cases_vacc / comparison_df$clinical_cases_no_vacc
comparison_df$time_since_3rd_dose = as.numeric((comparison_df$middle_month - vaccine_date) / 365)
colnames(comparison_df)[colnames(comparison_df)=='seasonality'] = 'site'
comparison_df_sim_only = comparison_df

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Process reference dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# also include reference data PE
ref_df = read.csv(reference_filepath)
ref_df_subset = ref_df
ref_df_subset = ref_df_subset[!grepl("M", ref_df_subset$range),]
ref_df_subset = ref_df_subset[ref_df_subset$age.group == 'children',]
ref_df_subset = ref_df_subset[ref_df_subset$site != 'OVERALL',]
ref_df_subset = ref_df_subset[ref_df_subset$measure == 'clinical',]
ref_df_subset = ref_df_subset[ref_df_subset$case.defn == 'primary',]
ref_df_subset$n.1 = as.numeric(ref_df_subset$n.1)
ref_df_subset$N.1 = as.numeric(ref_df_subset$N.1)
ref_df_subset$n.2 = as.numeric(ref_df_subset$n.2)
ref_df_subset$N.2 = as.numeric(ref_df_subset$N.2)
ref_df_subset$T.1 = as.numeric(ref_df_subset$T.1)
ref_df_subset$T.2 = as.numeric(ref_df_subset$T.2)
ref_df_subset$nT.1 = as.numeric(ref_df_subset$nT.1)
ref_df_subset$nT.2 = as.numeric(ref_df_subset$nT.2)
ref_df_subset$nT.LL.1 = as.numeric(ref_df_subset$nT.LL.1)
ref_df_subset$nT.UL.1 = as.numeric(ref_df_subset$nT.UL.1)
ref_df_subset$nT.LL.2 = as.numeric(ref_df_subset$nT.LL.2)
ref_df_subset$nT.UL.2 = as.numeric(ref_df_subset$nT.UL.2)
ref_df_subset$PE = as.numeric(ref_df_subset$PE)

ref_df_subset$PE2 = 1 - (ref_df_subset$nT.1/ref_df_subset$T.1) / (ref_df_subset$nT.2/ref_df_subset$T.2)
# second version using raw counts - not sure what adjustment was, since data was undocumented
#    ref_df_subset$PE2 = 1 - (ref_df_subset$n.1/ref_df_subset$N.1) / (ref_df_subset$n.2/ref_df_subset$N.2)
ref_df_subset$PE2l = 1 - (ref_df_subset$nT.LL.1/ref_df_subset$T.1) / (ref_df_subset$nT.LL.2/ref_df_subset$T.2)
ref_df_subset$PE2u = 1 - (ref_df_subset$nT.UL.1/ref_df_subset$T.1) / (ref_df_subset$nT.UL.2/ref_df_subset$T.2)

# get rid of values where sample sizes are extremely small
size_limit = 10
ref_df_subset$size_too_small = ref_df_subset$n.1<size_limit & ref_df_subset$n.2<size_limit
ref_df_subset$PE[ref_df_subset$size_too_small] = NA
ref_df_subset$PE2[ref_df_subset$size_too_small] = NA
ref_df_subset$PE2l[ref_df_subset$size_too_small] = NA
ref_df_subset$PE2u[ref_df_subset$size_too_small] = NA
ref_df_subset$N.1[ref_df_subset$size_too_small] = NA
ref_df_subset$N.2[ref_df_subset$size_too_small] = NA

ref_df_subset = ref_df_subset[ref_df_subset$comparison == booster_string, c('time.group', 'PE2', 'PE2l', 'PE2u', 'site', 'nT.1', 'T.1', 'n.1', 'N.1', 'nT.2', 'T.2', 'n.2', 'N.2')]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Plot comparison between simulation and reference
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
comparison_df = merge(comparison_df_sim_only, ref_df_subset, by = c('time.group', 'site'))


# plot similar to S2.4
gg = ggplot(data = comparison_df) +
  geom_line(aes(x = time_since_3rd_dose, y = PE_sim, col=factor(vacc_char)), lwd = 1) +
  geom_point(aes(x = time_since_3rd_dose, y = PE2)) +
  geom_errorbar(aes(x = time_since_3rd_dose, ymin=PE2l, ymax=PE2u)) +
  theme_light() +
  ylab('efficacy against clinical malaria') +
  xlab('time since third dose (years)') +
  ylim(-0.5, 1) + 
  facet_wrap(facets='site',nrow=3)
f_save_plot(gg, paste0('PE2_comparison'),
            file.path(exp_filepath), width = 12, height = 8, units = 'in', device_format = device_format)








#################################################
# which vaccine parametres did best?
#################################################

# calculate which vacc_char generates best log-likelihood with data
# ... nevermind, will take longer to do it correctly. instead, do a very rough and not-correct approach (probably okay for our purposes, though)
comparison_df2 = comparison_df
# given PE_sim, N1, N2, and either n1 or n2, what is (very, very rough) probability of observing the other n?
comparison_df2$p2 = (comparison_df2$n.1/comparison_df2$N.1) / (1- comparison_df2$PE_sim)
comparison_df2$nT2_prob = dbinom(x=comparison_df2$n.2, size=comparison_df2$N.2, prob=comparison_df2$p2)
comparison_df2$p1 = (comparison_df2$n.2/comparison_df2$N.2) * (1- comparison_df2$PE_sim)
comparison_df2$nT1_prob = dbinom(x=comparison_df2$n.1, size=comparison_df2$N.1, prob=comparison_df2$p1)
comparison_df2$ave_prob = (comparison_df2$nT1_prob + comparison_df2$nT2_prob)/2

vacc_char_sets = unique(comparison_df2$vacc_char)
rough_probs = rep(NA, length(vacc_char_sets))
for(vc in 1:length(vacc_char_sets)){
  comparison_cur = comparison_df2[comparison_df2$vacc_char==vacc_char_sets[vc],]
  print(nrow(comparison_cur))
  rough_probs[vc] = sum(log(comparison_cur$ave_prob), na.rm=TRUE)
}
print(paste0('Our vaccine winner is: ', vacc_char_sets[rough_probs == max(rough_probs)]))
df = data.frame('vacc_char'=vacc_char_sets, 'rough_loglik'=rough_probs)
write.csv(df, paste0(exp_filepath, '/compare_vacc_char_sets_sizeLimit', size_limit,'.csv'))



# plot the selected values from Phase 3 and Chandramohan
if(booster_string == 'boost'){
  comparison_df3 = comparison_df[as.character(comparison_df$vacc_char) %in% c('hh40_nn200_bb60', 'hh40_nn200_bb70', 'hh40_nn200_bb80'),]
} else{
  comparison_df3 = comparison_df[as.character(comparison_df$vacc_char) %in% c('hh40_nn200'),]
}
gg2 = ggplot(data = comparison_df3) +
  geom_line(aes(x = time_since_3rd_dose, y = PE_sim, col=factor(vacc_char)), lwd = 1) +
  geom_point(aes(x = time_since_3rd_dose, y = PE2)) +
  geom_errorbar(aes(x = time_since_3rd_dose, ymin=PE2l, ymax=PE2u)) +
  theme_light() +
  ylab('efficacy against clinical malaria') +
  xlab('time since third dose (years)') +
  ylim(-0.5, 1) + 
  facet_wrap(facets='site',nrow=3)
f_save_plot(gg2, paste0('PE2_comparison_SelectParams_sizeLimit', size_limit),
            file.path(exp_filepath), width = 12, height = 8, units = 'in', device_format = device_format)









# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Compare PE across reference sites - is it different for sites with different EIRs?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# add rough EIR associated with each site
site_eirs = data.frame(site=c("Agogo", "Bagamoyo", "Kilifi", "Kintampo",
                              "Kombewa","Korogwe","Lambarene", "Lilongwe" ,
                              "Manhica","Nanoro" ,"Siaya"), 
                       eir=c(4.9, 1.4, 0.4, 11.1,
                             9.7, 0.6, 1, 2.1,
                             1.1, 33.3, 27.8))
ref_df_eir = merge(ref_df_subset, site_eirs, by='site', all=TRUE)
# add EIR categories: under 5, 5-15, >15 
ref_df_eir$eir_category = '<5'
ref_df_eir$eir_category[ref_df_eir$eir>5] = '5-15'
ref_df_eir$eir_category[ref_df_eir$eir>15] = '>15'
ref_df_eir$eir_category = factor(ref_df_eir$eir_category, levels=c('<5', '5-15', '>15'))
# plot, with point shape determined by EIR
ggplot(ref_df_eir, aes(x=time.group, y=PE2, group=site, color=eir_category)) + 
  geom_point() +
  geom_smooth(data=ref_df_eir, aes(x=time.group, y=PE2, group=eir_category), method='loess') +
  theme_classic()



# aggregated from week 0-32 ; look at faceted by age group and boost/no boost
ref_df_agg = ref_df[intersect(which(ref_df$range == '[M0-32]'), which(ref_df$measure =='clinical')),]
ref_df_agg = ref_df_agg[ref_df_agg$site != 'OVERALL',]
ref_df_agg = ref_df_agg[ref_df_agg$case.defn == 'primary',]
ref_df_agg$PE = as.numeric(ref_df_agg$PE)
# merge with EIR groupings
ref_df_agg = merge(ref_df_agg, site_eirs, by='site', all=TRUE)
# add EIR categories: under 5, 5-15, >15 
ref_df_agg$eir_category = '<5'
ref_df_agg$eir_category[ref_df_agg$eir>5] = '5-15'
ref_df_agg$eir_category[ref_df_agg$eir>15] = '>15'
ref_df_agg$eir_category = factor(ref_df_agg$eir_category, levels=c('<5', '5-15', '>15'))

ggplot(ref_df_agg, aes(x=eir_category, y=PE))+#, color=age.group, shape=comparison))+
  geom_violin() +
  geom_point() +
  theme_classic() +
  facet_grid(rows=vars(comparison), cols=vars(age.group))

ggplot(ref_df_agg, aes(x=eir_category, y=PE))+#, color=age.group, shape=comparison))+
  geom_violin() +
  geom_point() +
  theme_classic() +
  facet_grid(rows=vars(comparison))

ggplot(ref_df_agg, aes(x=eir_category, y=PE))+#, color=age.group, shape=comparison))+
  geom_violin() +
  geom_point() +
  theme_classic() 










# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Compare the residual between reference and simulation points across reference sites - is it different for sites with different EIRs?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# plot the selected values from Phase 3 and Chandramohan
if(booster_string == 'boost'){
  comparison_df3 = comparison_df[as.character(comparison_df$vacc_char) %in% c('hh40_nn200_bb80'),]
} else{
  comparison_df3 = comparison_df[as.character(comparison_df$vacc_char) %in% c('hh40_nn200'),]
}

comparison_df3$eir_category = '<5'
comparison_df3$eir_category[comparison_df3$Annual.EIR>5] = '5-15'
comparison_df3$eir_category[comparison_df3$Annual.EIR>15] = '>15'
comparison_df3$eir_category = factor(comparison_df3$eir_category, levels=c('<5', '5-15', '>15'))

# plot faceted by EIR category, with site as color
ggplot(data = comparison_df3) +
  geom_line(aes(x = time_since_3rd_dose, y = PE_sim, col=factor(site)), lwd = 1) +
  geom_point(aes(x = time_since_3rd_dose, y = PE2, col=factor(site))) +
  geom_errorbar(aes(x = time_since_3rd_dose, ymin=PE2l, ymax=PE2u)) +
  theme_light() +
  ylab('efficacy against clinical malaria') +
  xlab('time since third dose (years)') +
  ylim(-0.5, 1) + 
  facet_wrap(facets='eir_category',nrow=3)


# plot all on same space, with EIR category as color
ggplot(data = comparison_df3) +
  geom_line(aes(x = time_since_3rd_dose, y = PE_sim, col=factor(eir_category), group_by=site), lwd = 1) +
  geom_line(aes(x = time_since_3rd_dose, y = PE2, col=factor(eir_category), group_by=site), lwd = 2, alpha=0.2) +
  geom_point(aes(x = time_since_3rd_dose, y = PE2, col=factor(eir_category), group_by=site), size=2) +
  geom_errorbar(aes(x = time_since_3rd_dose, ymin=PE2l, ymax=PE2u)) +
  theme_light() +
  ylab('efficacy against clinical malaria') +
  xlab('time since third dose (years)') +
  ylim(-0.5, 1) 

ggplot(data = comparison_df3) +
  geom_line(aes(x = time_since_3rd_dose, y = PE_sim, col=factor(eir_category), group_by=site), lwd = 1) +
  # geom_line(aes(x = time_since_3rd_dose, y = PE2, col=factor(eir_category), group_by=site), lwd = 2, alpha=0.2) +
  geom_point(aes(x = time_since_3rd_dose, y = PE2, col=factor(eir_category), group_by=site), size=4) +
  # geom_errorbar(aes(x = time_since_3rd_dose, ymin=PE2l, ymax=PE2u)) +
  theme_light() +
  ylab('efficacy against clinical malaria') +
  xlab('time since third dose (years)') +
  ylim(-0.5, 1) 

# get difference between reference point and matched simulation value
comparison_df3$diff = comparison_df3$PE2 - comparison_df3$PE_sim

