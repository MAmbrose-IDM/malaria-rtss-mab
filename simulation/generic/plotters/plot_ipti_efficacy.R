library(data.table)
library(dplyr)
library(ggplot2)

#wdir at rtss-scenarios
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path("simulation", "load_paths.R"))

simout_dir <- file.path(projectpath, 'simulation_output', 'generic_forward')
theme_set(theme_bw())
customTheme <- f_getCustomTheme()


exp_name <- "generic_IPTi_scl_only_IPTi_scl_v0"
exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'cm_coverage', 'ipti_coverage', 'rtss_coverage', 'smc_coverage',
                'intervention_correlation', 'frac_high_access', 'Cohort_birth_month', 'rtss_mode') #ipti_mode

plot_dir <- file.path(simout_dir, exp_name, '_plots')
if (!dir.exists(plot_dir))dir.create(plot_dir)

if (file.exists(file.path(simout_dir, exp_name, 'All_Age_Annual_Cases.csv'))) {
  fname = 'All_Age_Annual_Cases.csv'
}else {
  fname = 'All_Age_monthly_Cases.csv'
}

use_summaryreport = FALSE
if (!use_summaryreport) {
  plotname = 'IPTi_effect_sizes'
  cases_df_list <- load_Age_monthly_Cases(simout_dir, exp_name, exp_sweeps, max_years = c(0), fname = fname)

  ##Aggregated Runs and Cohort_birth_month
  dat <- cases_df_list[[4]] %>%
    filter(year==0) %>%
    dplyr::select(Scenario_id, cm_coverage, ipti_coverage, Annual_EIR, year, clinical_cases, pfpr, pop)

  grp_channels = c('cm_coverage',  'Annual_EIR', 'year') #ipti_mode
  dat <- data.table(dat, key =grp_channels)
  dat[, clinical_casesAverted := clinical_cases[ipti_coverage == 0.0] - clinical_cases, by = grp_channels]
  dat[, protective_efficacy := clinical_casesAverted / clinical_cases[ipti_coverage == 0.0], by =grp_channels]

}else {
  plotname = 'IPTi_effect_sizes_summaryReport'
  dat <- fread(file.path(simout_dir, exp_name, 'U1_PfPR_ClinicalIncidence.csv')) %>%
    rename_with(~gsub(" ", "_", .x)) %>%
    mutate(cases_per1000 = Cases_U1 / Pop_U1 * 1000) %>%
    select(year, Annual_EIR, cm_coverage, ipti_coverage, Cases_U1, cases_per1000, Pop_U1)

  dat <- data.table(dat, key = c('cm_coverage', 'Annual_EIR', 'year'))
  dat[, clinical_casesAverted := cases_per1000[ipti_coverage == 0.0] - cases_per1000, by = c('cm_coverage', 'Annual_EIR', 'year')]  #ipti_mode
  dat[, protective_efficacy := clinical_casesAverted / cases_per1000[ipti_coverage == 0.0], by = c('cm_coverage', 'Annual_EIR', 'year')]  #ipti_mode
}

dat$ipti_coverage_fct <- as.factor(dat$ipti_coverage)
dat$cm_coverage_fct <- as.factor(dat$cm_coverage)
dat$Annual_EIR_fct <- factor(dat$Annual_EIR, levels = unique(dat$Annual_EIR),
                             labels = paste0("annual EIR: ", unique(dat$Annual_EIR)))

pplot <- ggplot(data = subset(dat, ipti_coverage > 0)) +
  geom_bar(aes(x = cm_coverage_fct, y = protective_efficacy * 100, fill = ipti_coverage_fct,
               group = interaction(cm_coverage_fct, ipti_coverage_fct)), size = 0.2, col = 'lightgrey', stat = 'identity',
           position = position_dodge(1)) +
  scale_y_continuous(lim = c(0, 40)) +
  facet_wrap(~Annual_EIR_fct, nrow = 1, scales = 'free') +
  labs(col = 'ipti_coverage', fill = 'ipti_coverage', x = 'Treatment coverage', y = 'Protective efficacy') +
  scale_fill_brewer(palette = 'YlGnBu') +
  customTheme


print(pplot)
f_save_plot(pplot, plot_name = plotname,
            width = 10, height = 4,
            plot_dir = file.path(simout_dir, exp_name, '_plots'))


lineplot_over_time = TRUE
if (fname == 'All_Age_Annual_Cases.csv')lineplot_over_time = FALSE
if (lineplot_over_time) {
  ##_____ efficacy over age
  datU5 <- fread(file.path(simout_dir, exp_name, 'U5_PfPR_ClinicalIncidence.csv')) %>%
    rename_with(~gsub(" ", "_", .x)) %>%
    filter(ipti_coverage == 0) %>%
    mutate(PfPR_U5 = round(PfPR_U5 * 100, 0)) %>%
    select(Annual_EIR, cm_coverage, PfPR_U5)


  dat <- fread(file.path(simout_dir, exp_name, fname)) %>%
    rename_with(~gsub(" ", "_", .x)) %>%
    mutate(date = as.Date(date, format = '%Y-%m-%d'),
           timestep = as.numeric(date - as.Date('2020-01-01')),
           month = timestep / 30,
           New_Clinical_Cases = New_Clinical_Cases / Statistical_Population * 1000) %>%
    filter(timestep <= 730) %>%
    select(timestep, month, cm_coverage, ipti_coverage, Annual_EIR, New_Clinical_Cases)

  dat <- dat %>%
    left_join(datU5) %>%
    mutate(PfPR_U5_fct = paste0('(PfPR U5: ', PfPR_U5, '%)'),
           Annual_EIR_fct = paste0(Annual_EIR, '\n(PfPR U5: ', PfPR_U5, '%)'))

  dat$ipti_coverage_fct <- as.factor(dat$ipti_coverage)
  dat$cm_coverage_fct <- as.factor(dat$cm_coverage)
  dat$Annual_EIR_fct <- factor(dat$Annual_EIR, levels = unique(dat$Annual_EIR),
                               labels = paste0("annual EIR: ", unique(dat$Annual_EIR)))

  pplot <- ggplot() +
    geom_line(data = subset(dat, ipti_coverage == 0), aes(x = as.numeric(month), y = New_Clinical_Cases,
                                                          group = ipti_coverage_fct), col = 'black', size = 0.7) +
    geom_line(data = subset(dat, ipti_coverage > 0), aes(x = as.numeric(month), y = New_Clinical_Cases,
                                                         col = ipti_coverage_fct, group = ipti_coverage_fct), size = 1.1) +
    facet_wrap(cm_coverage_fct ~ Annual_EIR_fct + PfPR_U5_fct, scales = 'free') +
    labs(col = 'ipti coverage', x = 'Age (month)', y = 'New_Clinical_Cases per 1000') +
    #scale_fill_brewer(palette = 'YlGnBu') +
    customTheme

  print(pplot)
  f_save_plot(pplot, plot_name = 'IPTi_new_clinical_cases_by_EIR_CM',
              width = 14, height = 12,
              plot_dir = file.path(simout_dir, exp_name, '_plots'))


  ### first x days after treatment
  pplot <- ggplot() +
    geom_line(data = subset(dat, timestep <= 140 & ipti_coverage == 0),
              aes(x = as.numeric(month), y = New_Clinical_Cases,
                  group = ipti_coverage_fct), col = 'black', size = 0.7) +
    geom_line(data = subset(dat, timestep <= 140 & ipti_coverage > 0),
              aes(x = as.numeric(month), y = New_Clinical_Cases,
                  col = ipti_coverage_fct, group = ipti_coverage_fct), size = 1.1) +
    facet_wrap(cm_coverage_fct ~ Annual_EIR + PfPR_U5_fct, scales = 'free') +
    labs(col = 'ipti coverage', x = 'Age (months)', y = 'New_Clinical_Cases per 1000') +
    #scale_fill_brewer(palette = 'YlGnBu') +
    customTheme

  print(pplot)
  f_save_plot(pplot, plot_name = 'IPTi_new_clinical_cases_by_EIR_CM_zoom',
              width = 14, height = 12,
              plot_dir = file.path(simout_dir, exp_name, '_plots'))
}