# RTSS_validation_KintampoPE.R
# August 2021
# create plots of RTS,S simluation output to compare against Figure S2.4 from Penny et al. 

library(lubridate)
library(dplyr)
library(ggplot2)

setwd("C:/Users/moniqueam/Documents/malaria-rtss-mab")

source(file.path("simulation", "load_paths.R"))
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
device_format = c('pdf', 'png')
paths = get_project_paths()
datapath = paths[1]
projectpath = paths[2]

exp_name = 'TEST_validation_phase3Kintampo' # 'rtss_validation_KintampoPhase3_wBooster'
exp_filepath = file.path(projectpath, 'simulation_output', exp_name)
cases_filepath = file.path(exp_filepath, 'All_Age_monthly_Cases.csv')
reference_filepath = file.path(datapath, 'rtss_phase3/kintampo_trial_summary_3month.csv')
booster_string = 'boost'  #'no.boost' or 'boost'

ref_df = read.csv(reference_filepath)
cases_df = read.csv(cases_filepath)

# date vaccination given in simulations
vaccine_date = as.Date('2021-01-01')


# subset to time after vaccination
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

# get averages across runs
cases_df_ave = cases_df_grouped %>%
  dplyr::group_by(Scenario_id, time.group, Annual.EIR, vacc_coverage, vacc_char, seasonality) %>%
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

# also include reference data PE
ref_df_subset = ref_df[ref_df$comparison == booster_string, c('time.group', 'PE')]
comparison_df = merge(comparison_df, ref_df_subset, by = 'time.group')


# plot similar to S2.4
gg = ggplot(data = comparison_df) +
  geom_line(aes(x = time_since_3rd_dose, y = PE_sim, col=factor(vacc_char)), lwd = 2) +
  geom_point(aes(x = time_since_3rd_dose, y = PE)) +
  theme_light() +
  ylab('efficacy against clinical malaria') +
  xlab('time since third dose (years)') +
  ylim(-0.5, 1) + 
  facet_wrap(facets='seasonality')
f_save_plot(gg, paste0('PE_comparison'),
            file.path(exp_filepath), width = 3.25, height = 2.8, units = 'in', device_format = device_format)


# ggplot(data=comparison_df) + 
#   geom_line(aes(x=time_since_3rd_dose, y=PE_sim), col=rgb(0.4, 0.6, 0.3), lwd=2)+ 
#   geom_point(aes(x=time_since_3rd_dose, y=PE)) +
#   theme_light()+
#   ylab('efficacy against clinical malaria') + 
#   xlab('time since third dose (years)') +
#   ylim(-0.5,1) + 
#   theme_classic()



