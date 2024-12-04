# 01_process_sim_outputs.R
# August 2021
# from analyzer outputs, aggregate and format simulation results into format that is read in by Shiny app
#    calculate protective efficacy of various intervention scenarios relative to CM-only for U1, U2, and U5

library(lubridate)
library(dplyr)
library(ggplot2)
library(data.table)


boxpath = 'C:/Users/mambrose/Box/rtss_scenarios'
script_path = 'C:/Users/mambrose/Documents/rtss-scenarios'
exp_filepath = paste0(boxpath, '/simulation_output/generic/generic_standard_enhanced_season')  #/generic_cm_sweep_SMC') #generic_season_sweep_SMC')  #/generic_factorial_v1')  #rtss_generic_factorial')
cases_filepath = paste0(exp_filepath, '/All_Age_monthly_Cases.csv')

source(paste0(script_path, '/simulation/generic/plotters/plotter_helpers.R'))

if(!file.exists(cases_filepath)){
  # combine output from separate seeds or experiments
  # seed_endings = paste0('_r', 0:3)
  seed_endings = c('_r1', '_r2')
  for(rr in seed_endings){
    cur_run_df = fread(paste0(exp_filepath, '/All_Age_monthly_Cases', rr, '.csv'))
    if(!('minBoostAge' %in% colnames(cur_run_df))) cur_run_df$minBoostAge = 24
    if(rr==seed_endings[1]){
      combined_df = cur_run_df
    } else{
      combined_df = rbind(combined_df, cur_run_df)
    }
  }
  write.csv(combined_df, cases_filepath)
}


sim_output = load_Age_monthly_Cases(simout_dir=simout_dir, exp_name=exp_name, exp_sweeps=exp_sweeps,
                                    add_PE_perAge=TRUE, max_years=max_years, keep_birth_month=FALSE)
pe_df = sim_output[[3]]
pe_df$scenario_name = paste0(round(100*pe_df$rtss_coverage), '% RTS,S; ', round(100*pe_df$smc_coverage), '% SMC')
pe_df$scenario_name = factor(pe_df$scenario_name, levels=scenario_levels)

write.csv(pe_df, paste0(exp_filepath, '/protective_efficacy_estimates.csv'), row.names=FALSE)




##### combine multiple simulation outputs with added tags
tag_name = 'minBoostAge'
input_dir_1 = paste0(boxpath, '/simulation_output/generic/generic_enhanced_season_minSeasonBoostAge24m')
tag1 = 24

input_dir_2 = paste0(boxpath, '/simulation_output/generic/generic_enhanced_season_minSeasonBoostAge18m')
tag2=18

output_dir = paste0(boxpath, '/simulation_output/generic/generic_enhanced_season')


input1 =  fread(paste0(input_dir_1, '/All_Age_monthly_Cases.csv'))
input1[,tag_name] = tag1
input2 =  fread(paste0(input_dir_2, '/All_Age_monthly_Cases.csv'))
input2[,tag_name] = tag2
output = rbind(input1, input2)
if(!dir.exists(output_dir)) dir.create(output_dir)
write.csv(output, paste0(output_dir, '/All_Age_monthly_Cases.csv'), row.names=FALSE)






# plot simulation output's incidence as a function of age by eir
exp_filepath = paste0(boxpath, '/simulation_output/generic/generic_eir_sweep_SMC')   #generic_season_sweep_EIR_10_30') #generic_eir_sweep_SMC') 
sim_df = fread(paste0(exp_filepath, '/All_Age_monthly_Cases.csv'))
sim_df_ref = sim_df[sim_df$rtss_coverage == 0 & sim_df$smc_coverage==0 & sim_df$Cohort_birth_month==0,]
colnames(sim_df_ref) = gsub(' ', '_', colnames(sim_df_ref))
sim_df_ref$clinical_cases_per_person_per_year = sim_df_ref$New_Clinical_Cases/sim_df_ref$Statistical_Population * 12
sim_df_ref$date = as.Date(sim_df_ref$date)
sim_df_ref$age_months = as.numeric(round((sim_df_ref$date - min(sim_df_ref$date))/30.4))



ggplot(sim_df_ref, aes(x=age_months, y=clinical_cases_per_person_per_year)) + 
  geom_point() + 
  facet_wrap(seasonality~Annual_EIR)

ggplot(sim_df_ref, aes(x=age_months, y=PfHRP2_Prevalence)) + 
  geom_point() + 
  facet_wrap(seasonality~Annual_EIR)

ggplot(sim_df_ref, aes(x=PfHRP2_Prevalence, y=clinical_cases_per_person_per_year)) + 
  geom_point(aes(col=age_months)) + 
  facet_wrap(seasonality~Annual_EIR)





process_pfpr_sim_ouptut = function(filepath, exp_sweeps, rdt_sensitivity_multiplier=1){
  
  prev_df = fread(filepath) %>%
    rename_with(~gsub(" ", "_", .x)) %>%
    mutate(date = as.Date(date),
           year = lubridate::year(date),
           age = year - min(year))
  # get means across runs, cohorts
  if(any(!(exp_sweeps %in% colnames(prev_df)))) exp_sweeps <- exp_sweeps[-which(!(exp_sweeps %in% colnames(prev_df)))]
  if(('Cohort_birth_month' %in% exp_sweeps)) exp_sweeps <- exp_sweeps[-which(exp_sweeps=='Cohort_birth_month')]
  prev_df_ave = prev_df %>%
    dplyr::group_by_at(c(exp_sweeps, 'date')) %>%
    dplyr::summarise(New_Clinical_Cases = mean(New_Clinical_Cases),
                     True_Prevalence = mean(True_Prevalence),
                     Blood_Smear_Parasite_Prevalence = mean(Blood_Smear_Parasite_Prevalence),
                     PfHRP2_Prevalence = mean(PfHRP2_Prevalence),
                     PCR_Parasite_Prevalence = mean(PCR_Parasite_Prevalence),
                     Statistical_Population = mean(Statistical_Population)) %>%
    dplyr::ungroup() %>%
    mutate(incidence = New_Clinical_Cases / Statistical_Population * 1000) %>%
    mutate(date = as.Date(date),
           year = lubridate::year(date),
           age = year - min(year))
  prev_df_ave$PfHRP2_Prevalence = prev_df_ave$PfHRP2_Prevalence * rdt_sensitivity_multiplier
  return(prev_df_ave)
}

add_PfPR_2_10 = function(prev_df_ave, pfpr_colname, exp_sweeps){
  # cut out ages under 2 and over 10
  df = prev_df_ave[prev_df_ave$age>=2 & prev_df_ave$age<10, ]
  
  # aggregate across ages within each scenario
  if(any(!(exp_sweeps %in% colnames(df)))) exp_sweeps <- exp_sweeps[-which(!(exp_sweeps %in% colnames(df)))]
  df_ave = df %>%
    dplyr::group_by_at(c(exp_sweeps)) %>%
    dplyr::summarise(PfPR_2_10 = mean(get(pfpr_colname))) %>%
    dplyr::ungroup()

  # merge the 2-10 PfPR values into the main dataframe
  merged_df = merge(prev_df_ave, df_ave, by=exp_sweeps)
  return(merged_df)
}


# general parameters
output_path = 'C:/Users/mambrose/OneDrive - Institute for Disease Modeling/MalariaDiagnostics/prevalence_comparisons/simulation_output'
exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'intervention_correlation', 'cm_coverage', 'ipti_coverage',
                'rtss_coverage', 'smc_coverage', 'frac_high_access', 'Cohort_birth_month','rtss_mode')
rdt_sensitivity_multiplier = 1

# plot timeseries of different measures of prevalence for various EIRs
exp_name = 'generic_eir_sweep_noRTSSSMC'
filepath = paste0(output_path, '/', exp_name, '/All_Age_monthly_prevalence.csv')
prev_df_ave = process_pfpr_sim_ouptut(filepath=filepath, exp_sweeps=exp_sweeps, rdt_sensitivity_multiplier=rdt_sensitivity_multiplier)
prev_df_ave = add_PfPR_2_10(prev_df_ave=prev_df_ave, pfpr_colname='PfHRP2_Prevalence', exp_sweeps=exp_sweeps)
prev_df_ave1 = prev_df_ave


### with SMC, only one seasonality
exp_name = 'generic_eir_smc_sweep'
filepath = paste0(output_path, '/', exp_name, '/All_Age_monthly_prevalence.csv')
exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'intervention_correlation', 'cm_coverage', 'ipti_coverage',
                'rtss_coverage', 'smc_coverage', 'frac_high_access', 'Cohort_birth_month','rtss_mode')

prev_df_ave = process_pfpr_sim_ouptut(filepath=filepath, exp_sweeps=exp_sweeps, rdt_sensitivity_multiplier=rdt_sensitivity_multiplier)
prev_df_ave = add_PfPR_2_10(prev_df_ave=prev_df_ave, pfpr_colname='PfHRP2_Prevalence', exp_sweeps=exp_sweeps)
prev_df_ave2 = prev_df_ave

ggplot() + 
  geom_point(data=prev_df_ave1, aes(x=PfHRP2_Prevalence, y=Blood_Smear_Parasite_Prevalence, col=as.factor(Annual_EIR), shape=as.factor(smc_coverage))) +
  geom_point(data=prev_df_ave2, aes(x=PfHRP2_Prevalence, y=Blood_Smear_Parasite_Prevalence, col=as.factor(Annual_EIR), shape=as.factor(smc_coverage))) +
  geom_abline(intercept=0, slope=1) +
  ylim(0,1) + 
  xlim(0,1) + 
  ylab('microscopy prevalence') + 
  xlab('HRP2 RTD prevalence')


# ggplot(prev_df_ave, aes(x=date, y=True_Prevalence)) + 
#   geom_line(aes(col=as.factor(Annual_EIR), linetype=as.factor(smc_coverage))) + 
#   geom_line(aes(y=PCR_Parasite_Prevalence, col=as.factor(Annual_EIR)))
# ggplot(prev_df_ave, aes(x=PCR_Parasite_Prevalence, y=PfHRP2_Prevalence)) + 
#   geom_line(aes(col=as.factor(Annual_EIR), linetype=as.factor(smc_coverage))) +
#   geom_abline(intercept=0, slope=1)
# ggplot(prev_df_ave, aes(x=Blood_Smear_Parasite_Prevalence, y=PfHRP2_Prevalence)) + 
#   geom_point(aes(col=as.factor(Annual_EIR), shape=as.factor(smc_coverage))) +
#   geom_abline(intercept=0, slope=1)






# incidence against PfPR2-10
ggplot() + 
  geom_point(data=prev_df_ave1, aes(x=PfPR_2_10, y=New_Clinical_Cases, col=as.factor(Annual_EIR), shape=as.factor(smc_coverage))) +
  geom_point(data=prev_df_ave2, aes(x=PfPR_2_10, y=New_Clinical_Cases, col=as.factor(Annual_EIR), shape=as.factor(smc_coverage))) +
  xlim(0,1) + 
  xlab('HRP2 RDT prevalence 2-10') +
  ylab('incidence')





# timeseries comparison between old and new RDTs

# plot timeseries of different measures of prevalence for new and old RDT version
exp_name = 'test_RDT_updates_oldRDT'
filepath = paste0(output_path, '/', exp_name, '/All_Age_monthly_prevalence.csv')
prev_df_ave = process_pfpr_sim_ouptut(filepath=filepath, exp_sweeps=exp_sweeps, rdt_sensitivity_multiplier=rdt_sensitivity_multiplier)
prev_df_ave = add_PfPR_2_10(prev_df_ave=prev_df_ave, pfpr_colname='PfHRP2_Prevalence', exp_sweeps=exp_sweeps)
prev_df_ave1 = prev_df_ave
prev_df_ave1$RDT_age = 'old'

exp_name = 'test_RDT_updates_newRDT'
filepath = paste0(output_path, '/', exp_name, '/All_Age_monthly_prevalence.csv')
exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'intervention_correlation', 'cm_coverage', 'ipti_coverage',
                'rtss_coverage', 'smc_coverage', 'frac_high_access', 'Cohort_birth_month','rtss_mode')

prev_df_ave = process_pfpr_sim_ouptut(filepath=filepath, exp_sweeps=exp_sweeps, rdt_sensitivity_multiplier=rdt_sensitivity_multiplier)
prev_df_ave = add_PfPR_2_10(prev_df_ave=prev_df_ave, pfpr_colname='PfHRP2_Prevalence', exp_sweeps=exp_sweeps)
prev_df_ave2 = prev_df_ave
prev_df_ave2$RDT_age = 'new'

prev_df_ave = rbind(prev_df_ave1, prev_df_ave2)




ggplot(data=prev_df_ave, aes(x=date, y=PfHRP2_Prevalence,  group=as.factor(interaction(RDT_age, seasonality)), color=as.factor(RDT_age)))  +
  geom_line(aes(linetype=as.factor(seasonality))) +
  facet_grid(Annual_EIR ~ cm_coverage)

ggplot(data=prev_df_ave, aes(x=date, y=Blood_Smear_Parasite_Prevalence,  group=as.factor(interaction(RDT_age, seasonality)), color=as.factor(RDT_age)))  +
  geom_line(aes(linetype=as.factor(seasonality))) +
  facet_grid(Annual_EIR ~ cm_coverage)

ggplot(data=prev_df_ave, aes(x=date, y=PCR_Parasite_Prevalence,  group=as.factor(interaction(RDT_age, seasonality)), color=as.factor(RDT_age)))  +
  geom_line(aes(linetype=as.factor(seasonality))) +
  facet_grid(Annual_EIR ~ cm_coverage)

ggplot(data=prev_df_ave, aes(x=date, y=New_Clinical_Cases,  group=as.factor(interaction(RDT_age, seasonality)), color=as.factor(RDT_age)))  +
  geom_line(aes(linetype=as.factor(seasonality))) +
  facet_grid(Annual_EIR ~ cm_coverage)



# 
# #######################################
# #############   ARCHIVE   #############
# #######################################
# 
# # the remainder of this script has been moved to  plot_seasonality_eir_cm_sweeps.R, along with helper functions from the /simulation/generic/plotters directory
# 
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# #  calculate average clinical cases per 1000 per year for each age
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# 
# cases_df = fread(cases_filepath)
# colnames(cases_df) = gsub(' ', '_', colnames(cases_df))
# 
# # get total within each year
# cases_df$date = as.Date(cases_df$date)
# cases_df$year = lubridate::year(cases_df$date)
# cases_df$year = cases_df$year - min(cases_df$year) 
# 
# cases_df_annual = cases_df %>% dplyr::group_by(Scenario_id, Run_Number, Cohort_birth_month, year, Annual_EIR, seasonality, intervention_correlation, cm_coverage, rtss_coverage, smc_coverage) %>% 
#   dplyr::summarise(total_cases = sum(New.Clinical.Cases),
#                    total_severe_cases = sum(New.Severe.Cases),
#                    average_pfpr = mean(PfHRP2.Prevalence),
#                    average_pop = mean(Statistical.Population),
#   ) %>% 
#   dplyr::ungroup()
# 
# # get averages across runs and across birth cohorts
# cases_df_annual_ave = cases_df_annual %>% dplyr::group_by(Scenario_id, year, Annual_EIR, seasonality, intervention_correlation, cm_coverage, rtss_coverage, smc_coverage) %>% 
#   dplyr::summarise(clinical_cases = mean(total_cases),
#                    severe_cases = mean(total_severe_cases),
#                    pfpr = mean(average_pfpr),
#                    pop = mean(average_pop)
#   ) %>% 
#   dplyr::ungroup()
# 
# # clinical cases per 1000 people per year
# cases_df_annual_ave$clinical_cases = cases_df_annual_ave$clinical_cases / cases_df_annual_ave$pop * 1000
# cases_df_annual_ave$severe_cases = cases_df_annual_ave$severe_cases / cases_df_annual_ave$pop * 1000
# 
# 
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# #  calculate CM-only burdens for each scenario group
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# 
# # subset to scenarios with no RTS,S, no SMC
# no_cm_references = cases_df_annual_ave[intersect(which(cases_df_annual_ave$rtss_coverage==0),
#                                                  which(cases_df_annual_ave$smc_coverage==0)),]
# 
# # rename column to cm_only_clinical_cases
# no_cm_references = no_cm_references[, c('year', 'Annual_EIR', 'seasonality', 'cm_coverage', 'clinical_cases', 'severe_cases', 'intervention_correlation')]
# colnames(no_cm_references)[which(colnames(no_cm_references) == 'clinical_cases')] = 'cm_only_clinical_cases'
# colnames(no_cm_references)[which(colnames(no_cm_references) == 'severe_cases')] = 'cm_only_severe_cases'
# 
# 
# 
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# #  calculate no-RTS,S burdens as references for each scenario group
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# 
# # subset to scenarios with no RTS,S
# no_rtss_references = cases_df_annual_ave[which(cases_df_annual_ave$rtss_coverage==0),]
# 
# # rename column to cm_only_clinical_cases
# no_rtss_references = no_rtss_references[, c('year', 'Annual_EIR', 'seasonality', 'cm_coverage', 'smc_coverage', 'clinical_cases', 'severe_cases', 'intervention_correlation')]
# colnames(no_rtss_references)[which(colnames(no_rtss_references) == 'clinical_cases')] = 'no_rtss_clinical_cases'
# colnames(no_rtss_references)[which(colnames(no_rtss_references) == 'severe_cases')] = 'no_rtss_severe_cases'
# 
# 
# 
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# #  merge the CM-only and no-RTS,S baselines with scenarios
# # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# # remove the CM-only scenarios from the data frame
# # cases_scenarios = cases_df_annual_ave[-intersect(which(cases_df_annual_ave$rtss_coverage==0), which(cases_df_annual_ave$smc_coverage==0)),]
# cases_scenarios = cases_df_annual_ave
# 
# # merge dataframes
# cases_scenarios_references = merge(cases_scenarios, no_cm_references, by=c('year', 'Annual_EIR', 'seasonality', 'intervention_correlation', 'cm_coverage'), all.x=TRUE)
# cases_scenarios_references = merge(cases_scenarios_references, no_rtss_references, by=c('year', 'Annual_EIR', 'seasonality', 'intervention_correlation', 'cm_coverage', 'smc_coverage'), all.x=TRUE)
# 
# # aggregate dataframes to U1, U2, U5, and U10
# max_years = c(1, 2, 5, 10)
# for(i_max in 1:length(max_years)){
#   cases_scenarios_references_Um = cases_scenarios_references[cases_scenarios_references$year <= max_years[i_max],]
#   # get cases per 1000 across all included ages
#   cases_age_aggregated = cases_scenarios_references_Um %>% dplyr::group_by(Scenario_id, Annual_EIR, seasonality, intervention_correlation, cm_coverage, rtss_coverage, smc_coverage) %>% 
#     dplyr::summarise(clinical_cases = sum(clinical_cases),
#                      cm_only_clinical_cases = sum(cm_only_clinical_cases),
#                      no_rtss_clinical_cases = sum(no_rtss_clinical_cases),
#                      severe_cases = sum(severe_cases),
#                      cm_only_severe_cases = sum(cm_only_severe_cases),
#                      no_rtss_severe_cases = sum(no_rtss_severe_cases)
#     ) %>% 
#     dplyr::ungroup()
#   # convert to cases per 1000 per year
#   cases_age_aggregated$clinical_cases = cases_age_aggregated$clinical_cases / max_years[i_max]
#   cases_age_aggregated$cm_only_clinical_cases = cases_age_aggregated$cm_only_clinical_cases / max_years[i_max]
#   cases_age_aggregated$no_rtss_clinical_cases = cases_age_aggregated$no_rtss_clinical_cases / max_years[i_max]
#   cases_age_aggregated$severe_cases = cases_age_aggregated$severe_cases / max_years[i_max]
#   cases_age_aggregated$cm_only_severe_cases = cases_age_aggregated$cm_only_severe_cases / max_years[i_max]
#   cases_age_aggregated$no_rtss_severe_cases = cases_age_aggregated$no_rtss_severe_cases / max_years[i_max]
#   
#   # calculate protective efficacy and save dataframes for each age group
#   cases_age_aggregated$protective_efficacy = 1 - cases_age_aggregated$clinical_cases / cases_age_aggregated$cm_only_clinical_cases
#   cases_age_aggregated$protective_efficacy_severe = 1 - cases_age_aggregated$severe_cases / cases_age_aggregated$cm_only_severe_cases
#   cases_age_aggregated$rtss_protective_efficacy = 1 - cases_age_aggregated$clinical_cases / cases_age_aggregated$no_rtss_clinical_cases
#   cases_age_aggregated$rtss_protective_efficacy_severe = 1 - cases_age_aggregated$severe_cases / cases_age_aggregated$no_rtss_severe_cases
#   
#   # calculate burden relative to no RTS,S/SMC for each age group
#   cases_age_aggregated$relative_burden = (cases_age_aggregated$clinical_cases - cases_age_aggregated$cm_only_clinical_cases) / cases_age_aggregated$cm_only_clinical_cases
#   cases_age_aggregated$relative_burden_severe = (cases_age_aggregated$severe_cases - cases_age_aggregated$cm_only_severe_cases) / cases_age_aggregated$cm_only_severe_cases
#   
#   # cases averted per 100,000 children
#   cases_age_aggregated$cases_averted_per100000 = (cases_age_aggregated$cm_only_clinical_cases - cases_age_aggregated$clinical_cases) * 100
#   cases_age_aggregated$severe_cases_averted_per100000 = (cases_age_aggregated$cm_only_severe_cases - cases_age_aggregated$severe_cases) * 100
#   cases_age_aggregated$rtss_cases_averted_per100000 = (cases_age_aggregated$no_rtss_clinical_cases - cases_age_aggregated$clinical_cases) * 100
#   cases_age_aggregated$rtss_severe_cases_averted_per100000 = (cases_age_aggregated$no_rtss_severe_cases - cases_age_aggregated$severe_cases) * 100
#   
#   # label age group
#   cases_age_aggregated$age_group = paste0('U', max_years[i_max])
#   
#   
#   
#   if(i_max ==1){
#     pe_df = cases_age_aggregated
#   } else{
#     pe_df = rbind(pe_df, cases_age_aggregated)
#   }
# }
# 
# write.csv(pe_df, paste0(exp_filepath, '/protective_efficacy_estimates.csv'), row.names=FALSE)
# 
# pe_df = read.csv(paste0(exp_filepath, '/protective_efficacy_estimates.csv'))
# 
# # use function to get data frames
# sim_output = load_Age_monthly_Cases(simout_dir=paste0(boxpath, '/simulation_output/generic'), exp_name='generic_factorial_v1', exp_sweeps = c('Scenario_id', 'Annual_EIR', 'seasonality', 'cm_coverage', 'ipti_coverage', 'rtss_coverage', 'smc_coverage',
#                                                                                                                                               'intervention_correlation', 'frac_high_access', 'Cohort_birth_month', 'rtss_mode'),
#                                     add_PE_perAge=TRUE, max_years = c(1, 2, 5, 10), keep_birth_month=FALSE)
# cases_df_annual_ave = sim_output[[2]]
# pe_df = sim_output[[3]]
# cases_scenarios_each_year = sim_output[[4]]
# 
# ########################################################################
# # plot cases averted as a function of age for different scenarios
# 
# 
# cases_scenarios_each_year = cases_scenarios_references
# 
# # calculate protective efficacy and cases averted for each year
# 
# 
# # calculate protective efficacy and save dataframes for each age group
# cases_scenarios_each_year$protective_efficacy = 1 - cases_scenarios_each_year$clinical_cases / cases_scenarios_each_year$cm_only_clinical_cases
# cases_scenarios_each_year$protective_efficacy_severe = 1 - cases_scenarios_each_year$severe_cases / cases_scenarios_each_year$cm_only_severe_cases
# # calculate burden relative to no RTS,S/SMC for each age group
# cases_scenarios_each_year$relative_burden = (cases_scenarios_each_year$clinical_cases - cases_scenarios_each_year$cm_only_clinical_cases) / cases_scenarios_each_year$cm_only_clinical_cases
# cases_scenarios_each_year$relative_burden_severe = (cases_scenarios_each_year$severe_cases - cases_scenarios_each_year$cm_only_severe_cases) / cases_scenarios_each_year$cm_only_severe_cases
# # cases averted per 100,000 children
# cases_scenarios_each_year$cases_averted_per100000 = (cases_scenarios_each_year$cm_only_clinical_cases - cases_scenarios_each_year$clinical_cases) * 100
# cases_scenarios_each_year$severe_cases_averted_per100000 = (cases_scenarios_each_year$cm_only_severe_cases - cases_scenarios_each_year$severe_cases) * 100
# # label age group
# cases_scenarios_each_year$age_group = paste0(cases_scenarios_each_year$year, '-', (cases_scenarios_each_year$year+1))
# 
# 
# 
# 
# # myColors = c(rep(rgb(0.4,0.7,1),2), rep(rgb(0.4,1,0.7),2), rep(rgb(0.2,0.95,0.9),4))#, rep(rgb(1,0.4,0.7),2), rep(rgb(0.9,0.2,1),2))
# myColors = c(rgb(0.4,0.7,1), rgb(0.2,0.3,1), 
#              rep(rgb(0.4,1,0.7),2), 
#              rgb(0.2,0.95,0.9), rgb(0.01,0.85,0.8),
#              rgb(0,0.6,0.55), rgb(0,0.6,0.55))#, rep(rgb(1,0.4,0.7),2), rep(rgb(0.9,0.2,1),2))
# # plot(1:length(myColors), col=myColors, pch=20, cex=3)
# scenario_levels = c("60% RTS,S; 0% SMC", "80% RTS,S; 0% SMC", 
#                     "0% RTS,S; 60% SMC",  "0% RTS,S; 80% SMC",   
#                     "60% RTS,S; 60% SMC", "60% RTS,S; 80% SMC", 
#                     "80% RTS,S; 60% SMC", "80% RTS,S; 80% SMC")
# names(myColors) = scenario_levels
# fillScale = scale_fill_manual(name = "scenario_name", values = myColors)
# colScale = scale_color_manual(name = "scenario", values = myColors)
# 
# # subset
# cur_cm = 0.6
# cur_corr = 0
# cur_smc=0.8
# # subset simulation to current input
# cur_df = filter(cases_scenarios_each_year,
#                 intervention_correlation == cur_corr,
#                 cm_coverage == cur_cm,
#                 smc_coverage == cur_smc)
# cur_df$scenario_name = paste0(round(100*cur_df$rtss_coverage), '% RTS,S; ', round(100*cur_df$smc_coverage), '% SMC')
# cur_df$scenario_name = factor(cur_df$scenario_name, levels=scenario_levels)
# cur_df = cur_df[!is.na(cur_df$scenario_name),]
# 
# # create plot
# age_label_values = seq(0, 8, length.out=5)
# age_labels = paste0(age_label_values, '-', (age_label_values + 1))
# gg = ggplot(cur_df, aes(x=year, y=cases_averted_per100000))+
#   geom_point(aes(col=as.factor(scenario_name), group=scenario_name), size=2)+
#   geom_line(aes(col=as.factor(scenario_name), group=scenario_name), size=1)+
#   f_getColorMapping(name='scenario', type='color') +
#   theme_bw()+
#   f_getCustomTheme() +
#   geom_hline(yintercept=0)+
#   ylab('cases averted (per 100000) when adding RTS,S and/or SMC') + 
#   xlab('age of child') +
#   scale_x_continuous(breaks=age_label_values, labels=age_labels) +
#   facet_wrap(seasonality~Annual_EIR)
# if(!dir.exists(paste0(exp_filepath, '/_plots'))) dir.create(paste0(exp_filepath, '/_plots'))
# ggsave(filename=paste0(exp_filepath, '/_plots/cases_averted_by_age_RTSS_and_SMC_', 
#                        cur_cm*100, 'CM_',
#                        cur_smc*100, 'SMC_',
#                        cur_corr, 'corr',
#                        '.png'),
#        gg, width=8, height=5, units='in')
# 
# 
# # create plot
# ggplot(cur_df)+
#   geom_bar(stat = "identity", position='dodge', aes(x=as.factor(year), y=cases_averted_per100000, fill=as.factor(scenario_name), group=interaction(year, scenario_name)))+
#   theme_bw()+
#   f_getCustomTheme() +
#   geom_hline(yintercept=0)+
#   ylab('Cases averted (per 100000) when adding RTS,S and/or SMC') + 
#   xlab('Age of child') +
#   facet_wrap(seasonality~Annual_EIR)
# 
# 
# ######################################################################################
# # in an area with smc_coverage = X, how many cases averted by adding RTS,S?
# cur_age = 'U5'
# pe_df_cur = filter(pe_df,
#                    age_group == cur_age,
#                    intervention_correlation == cur_corr,
#                    cm_coverage == cur_cm
#                    )
# gg=ggplot(pe_df_cur, aes(x=smc_coverage, y=rtss_cases_averted_per100000)) +
#   geom_line(aes(col=as.factor(rtss_coverage)), size=2)+
#   ylab('Cases averted (per 100000) when adding RTS,S') + 
#   xlab('SMC coverage') +
#   theme_bw()+
#   f_getCustomTheme() +
#   facet_wrap(seasonality~Annual_EIR)
# if(!dir.exists(paste0(exp_filepath, '/_plots'))) dir.create(paste0(exp_filepath, '/_plots'))
# ggsave(filename=paste0(exp_filepath, '/_plots/cases_averted_RTSS_as_function_of_SMC_', 
#                        cur_age, '_',
#                        cur_cm*100, 'CM_',
#                        cur_corr, 'corr',
#                        '.png'),
#        gg, width=8, height=5, units='in')
# 
# # test_result = load_Age_monthly_Cases(simout_dir=paste0(boxpath, '/simulation_output/generic'),
# #                                      exp_name='rtss_generic_factorial',
# #                                      exp_sweeps =  c('Annual_EIR', 'seasonality', 'cm_coverage', 'rtss_coverage', 'smc_coverage', 'ipti_coverage',
# #                                                      'intervention_correlation',  'Cohort_birth_month'))
# 
# 
# 
# 
# ###########################################################################################
# # barplots of cases averted at different seasonality patterns (for two SMC coverages)
# 
# # subset to plotted scenario
# cur_cm = 0.6
# cur_corr = 0
# cur_rtss=0.8
# cur_age = 'U2'
# cur_df = filter(pe_df,
#                 age_group == cur_age,
#                 intervention_correlation == cur_corr,
#                 cm_coverage == cur_cm,
#                 rtss_coverage == cur_rtss
#                )
# cur_df$scenario_name = paste0(round(100*cur_df$rtss_coverage), '% RTS,S; ', round(100*cur_df$smc_coverage), '% SMC')
# cur_df$scenario_name = factor(cur_df$scenario_name, levels=scenario_levels)
# cur_df$seasonality = factor(cur_df$seasonality, levels=c('constant', 'moderate_unimodal', 'high_unimodal'))
# 
# 
# gg=plot_barplots(cur_df, xvar = 'seasonality', yvar = 'rtss_cases_averted_per100000',
#               fillvar = 'seasonality', facet1 = 'smc_coverage', facet2 = NULL,
#               SAVE = FALSE, ylab=paste0('cases averted by adding RTS,S per 100,000 ', cur_age))
# if(!dir.exists(paste0(exp_filepath, '/_plots'))) dir.create(paste0(exp_filepath, '/_plots'))
# ggsave(filename=paste0(exp_filepath, '/_plots/cases_averted_RTSS_seasonality_', 
#                        cur_age, '_',
#                        cur_cm*100, 'CM_',
#                        cur_corr, 'corr',
#                        '.png'),
#        gg, width=8, height=6, units='in')
# 
# 
# 
# 
# 
# 
# ###########################################################################################
# # lineplot of cases averted at different EIRs (for different age groups)
# 
# # subset to plotted scenario
# cur_cm = 0.6
# cur_corr = 0
# cur_rtss=0.8
# cur_smc=0
# cur_season = 'moderate_unimodal'
# cur_df = filter(pe_df,
#                 intervention_correlation == cur_corr,
#                 cm_coverage == cur_cm,
#                 smc_coverage == cur_smc,
#                 rtss_coverage == cur_rtss,
#                 seasonality == cur_season
# )
# cur_df$scenario_name = paste0(round(100*cur_df$rtss_coverage), '% RTS,S; ', round(100*cur_df$smc_coverage), '% SMC')
# cur_df$scenario_name = factor(cur_df$scenario_name, levels=scenario_levels)
# cur_df$age_group = factor(cur_df$age_group, levels=c('U1', 'U2', 'U5', 'U10'))
# 
# 
# gg=plot_lines(dat=cur_df, xvar = 'Annual_EIR', yvar = 'rtss_severe_cases_averted_per100000',
#                  colvar = 'age_group', facet1 = NULL, facet2 = NULL,
#                  SAVE = FALSE, xlab='Annual EIR', ylab=paste0('severe cases averted by adding RTS,S per 100,000'))
# 
# if(!dir.exists(paste0(exp_filepath, '/_plots'))) dir.create(paste0(exp_filepath, '/_plots'))
# ggsave(filename=paste0(exp_filepath, '/_plots/severe_cases_averted_RTSS_EIR_', 
#                        cur_season, '_',
#                        cur_cm*100, 'CM_',
#                        cur_smc*100, 'SMC_',
#                        cur_rtss*100, 'RTSS_',
#                        cur_corr, 'corr',
#                        '.png'),
#        gg, width=8, height=5, units='in')
# 
# 
# 
# 
# 
# ###########################################################################################
# # lineplot of severe cases averted at different CM rates (for different age groups)
# 
# # subset to plotted scenario
# cur_corr = 0
# cur_rtss=0.8
# cur_smc=0
# cur_season = 'moderate_unimodal'
# cur_df = filter(pe_df,
#                 intervention_correlation == cur_corr,
#                 smc_coverage == cur_smc,
#                 rtss_coverage == cur_rtss,
#                 seasonality == cur_season
# )
# cur_df$scenario_name = paste0(round(100*cur_df$rtss_coverage), '% RTS,S; ', round(100*cur_df$smc_coverage), '% SMC')
# cur_df$scenario_name = factor(cur_df$scenario_name, levels=scenario_levels)
# cur_df$age_group = factor(cur_df$age_group, levels=c('U1', 'U2', 'U5', 'U10'))
# 
# 
# gg=plot_lines(dat=cur_df, xvar = 'cm_coverage', yvar = 'rtss_severe_cases_averted_per100000',
#               colvar = 'age_group', facet1 = NULL, facet2 = NULL,
#               SAVE = FALSE, xlab='Effective case management rate', ylab=paste0('severe cases averted by adding RTS,S per 100,000'))
# 
# if(!dir.exists(paste0(exp_filepath, '/_plots'))) dir.create(paste0(exp_filepath, '/_plots'))
# ggsave(filename=paste0(exp_filepath, '/_plots/severe_cases_averted_RTSS_CM_', 
#                        cur_season, '_',
#                        cur_cm*100, 'CM_',
#                        cur_smc*100, 'SMC_',
#                        cur_rtss*100, 'RTSS_',
#                        cur_corr, 'corr',
#                        '.png'),
#        gg, width=8, height=5, units='in')
