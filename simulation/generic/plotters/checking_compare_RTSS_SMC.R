
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

# exp_name_comp1 = 'sweep1_seeds1'
# exp_name_comp2 = 'sweep2_seeds1'
exp_name = 'sweep4_seeds1'

simout_dir=file.path(projectpath, 'simulation_output')
if (!dir.exists(paste0(simout_dir, '/_plots'))) dir.create(paste0(simout_dir, '/_plots'))
if (!dir.exists(paste0(simout_dir, '/_plots/pdf'))) dir.create(paste0(simout_dir, '/_plots/pdf'))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Process simulation output and create line plots comparing cases 
#   averted with RTS,S versus mAbs across mAb params
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# exp_name = exp_name_comp1
sim_output = load_Age_monthly_Cases(simout_dir=simout_dir, exp_name=exp_name, add_PE_perAge=TRUE,
                                    max_years=c(5, 8), keep_birth_month=FALSE, fname='All_Age_monthly_Cases.csv')

pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'vacc_info')
pe_df = pe_df[!is.na(pe_df$vacc_info),]
pe_df$hh[pe_df$vacc_type %in% c('rtss', 'RTS,S')]=NA
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))

# subset simulations
cur_cm = 0.6
cur_age = 'U5'
cur_eirs = c(5, 10, 30)
cur_seasonalities = c('constant', 'moderate_unimodal', 'high_unimodal')
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   seasonality %in% cur_seasonalities,
                   Annual_EIR %in% cur_eirs,
                   cm_coverage == cur_cm,
)
pe_df_cur$seasonality = factor(pe_df_cur$seasonality, levels=cur_seasonalities)
# remove mAb rows
pe_df_cur = pe_df_cur[-intersect(which(pe_df_cur$vacc_type %in% c('mab', 'mAb')), which(pe_df_cur$vacc_coverage>0)),]




# compare RTSS against matched SMC impact. Start by getting 'reference' of SMC
match_cols = c('Annual_EIR', 'seasonality', 'age_group', 'cm_coverage', 'frac_high_access', 'cm_target_group')
keep_cols = c('cases_averted_per100000', 'severe_cases_averted_per100000')
# subset to SMC
pe_df_smc = pe_df_cur[pe_df_cur$smc_coverage > 0,]
df_smc = pe_df_smc %>%
  dplyr::select(c(match_cols, keep_cols)) %>%
  rename(cases_averted_per100000_smc = cases_averted_per100000,
         severe_cases_averted_per100000_smc = severe_cases_averted_per100000)
# subset to RTSS
pe_df_rtss = pe_df_cur[pe_df_cur$vacc_coverage > 0,]
df_rtss = pe_df_rtss %>%
  dplyr::select(c(match_cols, keep_cols)) %>%
  rename(cases_averted_per100000_rtss = cases_averted_per100000,
         severe_cases_averted_per100000_rtss = severe_cases_averted_per100000)

# merge the reference values into the data frame
pe_df_compare <- df_rtss %>% left_join(df_smc)
ggplot(pe_df_compare, aes(x=cases_averted_per100000_rtss, y=cases_averted_per100000_smc, col=seasonality, shape=Annual_EIR), size=2)+
  geom_point()+
  ylim(15000, 120000)+
  xlim(15000, 120000)+
  geom_abline(intercept=0, slope=1) +
  theme_bw()
plot(pe_df_compare$cases_averted_per100000_rtss, pe_df_compare$cases_averted_per100000_smc)
abline(a=0, b=1)




##########
# see whether SMC is mis-timed
cases_df = fread(file.path(simout_dir, exp_name, 'All_Age_monthly_Cases.csv')) %>% rename_with(~gsub(" ", "_", .x))


exp_sweeps = c('Scenario_id', 'Annual_EIR', 'seasonality', 'Cohort_birth_month', 'vacc_char',
                              'vacc_coverage', 'cm_coverage', 'smc_coverage', 'frac_high_access',
                              'cm_target_group', 'smc_target_group', 'vacc_target_group')
if (any(!(exp_sweeps %in% colnames(cases_df)))) exp_sweeps <- exp_sweeps[-which(!(exp_sweeps %in% colnames(cases_df)))]
# if ('Cohort_birth_month' %in% exp_sweeps) exp_sweeps <- exp_sweeps[-which(exp_sweeps == 'Cohort_birth_month')]
# get averages across runs and optional across birth cohorts; calculate rates per 1000 people per year
cases_df_ave = cases_df %>%
  dplyr::group_by_at(c(exp_sweeps, 'date')) %>%
  dplyr::summarise(clinical_cases = mean(New_Clinical_Cases),
                   severe_cases = mean(New_Severe_Cases),
                   pfpr = mean(PfHRP2_Prevalence),
                   pop = mean(Statistical_Population)
  ) %>%
  dplyr::ungroup() %>%
  mutate(clinical_cases = clinical_cases / pop * 1000,
         severe_cases = severe_cases / pop * 1000)

# subset to relevant RTSS and SMC scenarios
cases_df_ave = f_add_scenario_name(df = cases_df_ave, scenario_type = 'vacc_info')
cases_df_ave$Annual_EIR = factor(cases_df_ave$Annual_EIR, levels = sort(unique(cases_df_ave$Annual_EIR)))
cur_seasonalities = c('constant', 'moderate_unimodal', 'high_unimodal')
cases_df_ave$seasonality = factor(cases_df_ave$seasonality, levels=cur_seasonalities)
cases_df_ave = cases_df_ave[-intersect(which(cases_df_ave$vacc_type %in% c('mab', 'mAb')), which(cases_df_ave$vacc_coverage>0)),]
# cases_df_ave = cases_df_ave[-intersect(which(cases_df_ave$smc_coverage<0.001), which(cases_df_ave$vacc_coverage<0.001)),]
cases_df_ave$label='none'
cases_df_ave$label[cases_df_ave$smc_coverage>0.001] = 'SMC'
cases_df_ave$label[cases_df_ave$vacc_coverage>0.001] = 'RTS,S'
# subset to dates
cases_df_ave = cases_df_ave[cases_df_ave$date<as.Date("2024-01-01"),]
cases_df_ave = cases_df_ave[cases_df_ave$date>=as.Date("2023-01-01"),]
# subset to single cohort birth month (otherwise will look like flat seasonality because we're averaging over cohorts months)
cases_df_ave = cases_df_ave[cases_df_ave$Cohort_birth_month == 0,]

ggplot(cases_df_ave, aes(x=date, y=clinical_cases, color=label))+
  geom_line() +
  facet_grid(Annual_EIR~seasonality)
