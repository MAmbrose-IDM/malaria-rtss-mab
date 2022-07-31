# RTSS_validation_ChandramohanTrial.R
# August 2021
# create plots of RTS,S/SMC simluation output to compare against Figure 2 from Chandramohan et al. 2021

library(lubridate)
library(dplyr)
library(ggplot2)
library(scales)
library(lemon)
library(reshape2)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Setup filepaths and parameters
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setwd("C:/Users/moniqueam/Documents/malaria-rtss-mab")

source(file.path("simulation", "load_paths.R"))
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path('simulation', 'generic', 'plotters', 'plot_bars_and_lines.R'))
paths = get_project_paths()
datapath = paths[1]
projectpath = paths[2]

exp_filepath = file.path(projectpath, 'simulation_output/validation_ChandramohanTrial_p4_cleaned_3seeds')  # validation_ChandramohanTrial_cleaned_1seed
cases_filepath = file.path(exp_filepath, 'All_Age_monthly_Cases.csv')
cases_df = read.csv(cases_filepath)
eir_cur = 6  # 6, 10, 15
cm_cur = 0.8  # 0.9, 0.8
device_format <- c('pdf','png')

TIMESERIES_WITH_SCATTER = F
if (TIMESERIES_WITH_SCATTER) {
  library(tidyr)
  library(cowplot)
}


reference_filepath = file.path(datapath, 'reference_Thompson/plotGrabbed_FigA1.csv')
ref_df = read.csv(reference_filepath)


# shift all dates so that start date is June 1, 2016 (one year before 3rd vaccine)
sim_start_date = as.Date('2016-06-01')
cases_df$date = as.Date(cases_df$date)
cases_df$date0 = cases_df$date
sim_shift = min(cases_df$date) - sim_start_date
cases_df$date = round_date(cases_df$date0 - sim_shift, unit='month')  # round to first of month
#cases_df$date = cases_df$date0 - sim_shift  # round to first of month ## if Rccp error occurs for round_date
#cases_df$date = make_date(year(cases_df$date), month(cases_df$date))

#######################################################################################
# timeseries: comparison with Chandramohan et al. 2021 Figure 2
#######################################################################################


# get average within each reported timepoint (month) across runs
cases_df_ave = cases_df %>%
  dplyr::group_by(Scenario_id, date, Annual.EIR, seasonality, Cohort_birth_month, cm_coverage, smc_coverage, vacc_coverage, vacc_char, vacc_mode, cm_target_group, smc_target_group, vacc_target_group) %>%
  dplyr::summarise(clinical_cases = mean(Received_Treatment),
                   clinical_cases_total = mean(New.Clinical.Cases),
                   severe_cases = mean(New.Severe.Cases),
                   pfpr = mean(PfHRP2.Prevalence),
                   pop = mean(Statistical.Population)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(clinical_cases = clinical_cases / pop * 1000,
                clinical_cases_total = clinical_cases_total / pop * 1000,
                severe_cases = severe_cases / pop * 1000)

# label intervention arm
cases_df_ave$intervention_arm = 'No intervention'
cases_df_ave$intervention_arm[intersect(which(cases_df_ave$smc_coverage > 0), which(cases_df_ave$vacc_coverage > 0))] = 'Combination'
cases_df_ave$intervention_arm[intersect(which(cases_df_ave$smc_coverage > 0), which(cases_df_ave$vacc_coverage < 0.0001))] = 'Chemoprevention Alone'
cases_df_ave$intervention_arm[intersect(which(cases_df_ave$smc_coverage < 0.0001), which(cases_df_ave$vacc_coverage > 0))] = 'Vaccine Alone'
cases_df_ave$intervention_arm = factor(cases_df_ave$intervention_arm, levels = c('Chemoprevention Alone', 'Vaccine Alone', 'Combination', 'No intervention'))

# create colors for each intervention arm
myColors = c('No intervention' = rgb(0, 0, 0),
             'Chemoprevention Alone' = rgb(110 / 255, 147 / 255, 201 / 255),
             'Vaccine Alone' = rgb(207 / 255, 106 / 255, 91 / 255),
             'Combination' = rgb(160 / 255, 113 / 255, 178 / 255))
colScale = scale_colour_manual(name = "intervention_arm", values = myColors)

# subset to plotted timeperiod
date_start_plot = as.Date('2017-04-01')
date_end_plot = as.Date('2020-03-05')
cases_df_subset = cases_df_ave[intersect(which(cases_df_ave$date >= date_start_plot), which(cases_df_ave$date < date_end_plot)),]
# plotted_dates = seq(date_start_plot, date_start_plot+365*3-30, length.out=36)
plotted_dates = c(as.Date(paste0('2017-', seq(4, 12, 2), '-01')), as.Date(paste0('2018-', seq(2, 12, 2), '-01')), as.Date(paste0('2019-', seq(2, 12, 2), '-01')), as.Date(paste0('2020-02-01')))

# subset to scenario
seasonality_options = unique(cases_df_subset$seasonality)
for (ss in 1:length(seasonality_options)) {
  cases_df_cur = cases_df_subset[intersect(intersect(which(cases_df_subset$seasonality == seasonality_options[ss]),
                                                     which(cases_df_subset$Annual.EIR == eir_cur)),
                                           which(cases_df_subset$cm_coverage == cm_cur)),]
  cases_df_cur = cases_df_cur[as.character(cases_df_cur$intervention_arm) != 'No intervention',]
  gg = ggplot(cases_df_cur, aes(x = date, y = clinical_cases, col = intervention_arm)) +
    geom_point() +
    geom_line(lwd = 0.5, alpha = 0.2) +
    colScale +
    theme_classic() +
    ylab('Clinical cases per 1000 per month') +
    # ylim(0,175)+
    # scale_x_continuous(breaks=18) +
    scale_x_date(labels = date_format("%B %Y"), breaks = plotted_dates) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank(), axis.line = element_line(), axis.title.x = element_blank(), legend.position = 'none') +
    facet_rep_wrap(facets = vars(intervention_arm), ncol = 1) #, scales='free_y')
  # facet_wrap(facets=vars(intervention_arm), ncol=1, scales = "free_x")
  f_save_plot(gg, paste0('timeseries_reportedClinicalCases_EIR',eir_cur, '_', seasonality_options[ss], '_cm', cm_cur * 100),
              file.path(exp_filepath), width = 6.4, height = 9.5, units = 'in', device_format = device_format)

  print(gg)

  if (TIMESERIES_WITH_SCATTER) {
    cases_df_cur <- cases_df_subset[intersect(which(cases_df_subset$seasonality == seasonality_options[ss]), which(cases_df_subset$Annual.EIR == eir_cur)),]
    cases_df_cur <- cases_df_cur[as.character(cases_df_cur$intervention_arm) != "No intervention",]
    gg <- ggplot(cases_df_cur, aes(x = date, y = clinical_cases, col = intervention_arm)) +
      geom_point() +
      geom_line(lwd = 0.5, alpha = 0.2) +
      colScale +
      theme_classic() +
      ylab("Clinical cases per 1000 per month") +
      scale_x_date(labels = date_format("%B %Y"), breaks = plotted_dates) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank(), axis.line = element_line(), axis.title.x = element_blank(), legend.position = "none") +
      facet_rep_wrap(facets = vars(intervention_arm), ncol = 1)


    sim_dat <- cases_df_cur %>%
      dplyr::select(date, intervention_arm, clinical_cases) %>%
      rename(clinical_cases_sim = clinical_cases) # %>% mutate(df_type='simulation')

    ref_dat <- fread(file.path(datapath,"data", "Chandramohan_2021", "chandramohan_2021_Fig1.csv")) %>%
      dplyr::select(-timepoint) %>%
      pivot_longer(cols = -c(date)) %>%
      rename(intervention_arm = name, clinical_cases_ref = value) %>%
      mutate(
        clinical_cases_ref = clinical_cases_ref,
        date = as.Date(date, format = "%m/%d/%Y")
      ) # df_type='reference'

    sim_dat$intervention_arm <- factor(sim_dat$intervention_arm,
                                       levels = c("Chemoprevention Alone", "Vaccine Alone", "Combination"),
                                       labels = c("Chemoprevention Alone", "Vaccine Alone", "Combination")
    )
    ref_dat$intervention_arm <- factor(ref_dat$intervention_arm,
                                       levels = c("Chemoprevention Alone", "Vaccine Alone", "Combination"),
                                       labels = c("Chemoprevention Alone", "Vaccine Alone", "Combination")
    )
    sim_dat$date <- as.character(sim_dat$date)
    ref_dat$date <- as.character(ref_dat$date)
    plot_df <- left_join(sim_dat, ref_dat)

    plot_df <- plot_df %>%
      group_by(date, intervention_arm) %>%
      summarize(clinical_cases_sim = mean(clinical_cases_sim), clinical_cases_ref = mean(clinical_cases_ref))

    sim_dat <- sim_dat %>%
      group_by(date, intervention_arm) %>%
      summarize(clinical_cases_sim = mean(clinical_cases_sim))

    pp2 <- ggplot(data = plot_df) +
      geom_abline(intercept = 0, slope = 1) +
      geom_smooth(aes(x = clinical_cases_ref, y = clinical_cases_sim, col = intervention_arm, fill = intervention_arm), size = 1.3, method = 'lm') +
      geom_point(aes(x = clinical_cases_ref, y = clinical_cases_sim, fill = intervention_arm), shape = 21, size = 2) +
      scale_y_continuous(lim = c(-5, 200)) +
      scale_x_continuous(lim = c(-5, 200)) +
      theme_classic() +
      scale_color_manual(values = myColors) +
      scale_fill_manual(values = myColors) +
      labs(y = "Clinical cases per 1000 per month") +
      theme(strip.background = element_blank(), axis.line = element_line(), axis.title.x = element_blank(), legend.position = "none") +
      facet_rep_wrap(facets = vars(intervention_arm), ncol = 1, scales = "free")


    pp1 <- ggplot(data = sim_dat) +
      geom_line(data = sim_dat, aes(x = as.Date(date), y = clinical_cases_sim, col = intervention_arm), lwd = 0.5) +
      geom_line(data = sim_dat, aes(x = as.Date(date), y = clinical_cases_sim, col = intervention_arm), lwd = 0.5) +
      geom_point(data = ref_dat, aes(x = as.Date(date), y = clinical_cases_ref), col = "black", shape = 22) +
      theme_classic() +
      scale_color_manual(values = myColors) +
      ylab("Clinical cases per 1000 per month") +
      # scale_x_continuous(breaks=18) +
      scale_x_date(labels = date_format("%B %Y"), breaks = plotted_dates) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), strip.background = element_blank(), axis.line = element_line(), axis.title.x = element_blank(), legend.position = "none") +
      facet_rep_wrap(facets = vars(intervention_arm), ncol = 1)

    gg <- plot_grid(pp1, pp2, ncol = 2, rel_widths = c(1, 0.75), rel_heights = c(1, 0.8))
    gg

    f_save_plot(gg, paste0('validation_plot_scatter_', eir_cur, '_', seasonality_options[ss], '_cm', cm_cur * 100),
                file.path(exp_filepath), width = 7, height = 7, units = 'in', device_format = device_format)
  }
}


##################################################################################################################
# calculate incidence rate ratios between intervention arms (compare to Figure A1 from seasonal_use_case_2_pager)
##################################################################################################################


# get total within each Trial Year (from April to March each year) - only include  the first three years

# calculate year of trial
cases_df$trial_year_diff = time_length(difftime(cases_df$date, as.Date("2017-04-01")), "years")
cases_df_trial_years = cases_df[cases_df$trial_year_diff > -0.001,]
cases_df_trial_years$trial_year_diff = round(cases_df_trial_years$trial_year_diff, 2)
cases_df_trial_years$trial_year_diff = floor(cases_df_trial_years$trial_year_diff)
cases_df_trial_years = cases_df_trial_years[cases_df_trial_years$trial_year_diff <= 2,]
cases_df_trial_years$trial_year_diff = cases_df_trial_years$trial_year_diff + 1

# aggregate within a trial year
cases_df_annual = cases_df_trial_years %>%
  dplyr::group_by(Scenario_id, Run_Number, trial_year_diff, Annual.EIR, seasonality, Cohort_birth_month, cm_coverage, smc_coverage, vacc_coverage, vacc_char, vacc_mode, cm_target_group, smc_target_group, vacc_target_group) %>%
  dplyr::summarise(total_treated_cases = sum(Received_Treatment),
                   total_cases = sum(New.Clinical.Cases),
                   total_severe_cases = sum(New.Severe.Cases),
                   average_pfpr = mean(PfHRP2.Prevalence),
                   average_pop = mean(Statistical.Population),
  ) %>%
  dplyr::ungroup()
# get averages across runs
cases_df_annual_ave = cases_df_annual %>%
  dplyr::group_by(Scenario_id, trial_year_diff, Annual.EIR, seasonality, cm_coverage, smc_coverage, vacc_coverage, vacc_char, vacc_mode, cm_target_group, smc_target_group, vacc_target_group) %>%
  dplyr::summarise(treated_clinical_cases = mean(total_treated_cases),
                   clinical_cases = mean(total_cases),
                   severe_cases = mean(total_severe_cases),
                   pfpr = mean(average_pfpr),
                   pop = mean(average_pop)
  ) %>%
  dplyr::ungroup()


# aggregate across all trial years
cases_df_all_years = cases_df_trial_years %>%
  dplyr::group_by(Scenario_id, Run_Number, Annual.EIR, seasonality, Cohort_birth_month, cm_coverage, smc_coverage, vacc_coverage, vacc_char, vacc_mode, cm_target_group, smc_target_group, vacc_target_group) %>%
  dplyr::summarise(total_treated_cases = sum(Received_Treatment),
                   total_cases = sum(New.Clinical.Cases),
                   total_severe_cases = sum(New.Severe.Cases),
                   average_pfpr = mean(PfHRP2.Prevalence),
                   average_pop = mean(Statistical.Population),
  ) %>%
  dplyr::ungroup()
# get averages across runs
cases_df_all_years_ave = cases_df_all_years %>%
  dplyr::group_by(Scenario_id, Annual.EIR, seasonality, cm_coverage, smc_coverage, vacc_coverage, vacc_char, vacc_mode, cm_target_group, smc_target_group, vacc_target_group) %>%
  dplyr::summarise(treated_clinical_cases = mean(total_treated_cases),
                   clinical_cases = mean(total_cases),
                   severe_cases = mean(total_severe_cases),
                   pfpr = mean(average_pfpr),
                   pop = mean(average_pop)
  ) %>%
  dplyr::ungroup()

# label intervention arm
cases_df_annual_ave$intervention_arm = 'No intervention'
cases_df_annual_ave$intervention_arm[intersect(which(cases_df_annual_ave$smc_coverage > 0), which(cases_df_annual_ave$vacc_coverage > 0))] = 'Combination'
cases_df_annual_ave$intervention_arm[intersect(which(cases_df_annual_ave$smc_coverage > 0), which(cases_df_annual_ave$vacc_coverage < 0.0001))] = 'Chemoprevention Alone'
cases_df_annual_ave$intervention_arm[intersect(which(cases_df_annual_ave$smc_coverage < 0.0001), which(cases_df_annual_ave$vacc_coverage > 0))] = 'Vaccine Alone'
cases_df_annual_ave$intervention_arm = factor(cases_df_annual_ave$intervention_arm, levels = c('Chemoprevention Alone', 'Vaccine Alone', 'Combination', 'No intervention'))
cases_df_all_years_ave$intervention_arm = 'No intervention'
cases_df_all_years_ave$intervention_arm[intersect(which(cases_df_all_years_ave$smc_coverage > 0), which(cases_df_all_years_ave$vacc_coverage > 0))] = 'Combination'
cases_df_all_years_ave$intervention_arm[intersect(which(cases_df_all_years_ave$smc_coverage > 0), which(cases_df_all_years_ave$vacc_coverage < 0.0001))] = 'Chemoprevention Alone'
cases_df_all_years_ave$intervention_arm[intersect(which(cases_df_all_years_ave$smc_coverage < 0.0001), which(cases_df_all_years_ave$vacc_coverage > 0))] = 'Vaccine Alone'
cases_df_all_years_ave$intervention_arm = factor(cases_df_all_years_ave$intervention_arm, levels = c('Chemoprevention Alone', 'Vaccine Alone', 'Combination', 'No intervention'))


# subset to scenario
seasonality_options = unique(cases_df_annual_ave$seasonality)
for (ss in 1:length(seasonality_options)) {
  # subset to relevant part of data frame - separate years
  cases_df_cur = cases_df_annual_ave[intersect(intersect(which(cases_df_annual_ave$seasonality == seasonality_options[ss]),
                                                         which(cases_df_annual_ave$Annual.EIR == eir_cur)),
                                               which(cases_df_annual_ave$cm_coverage == cm_cur)),]

  cases_df_cur = cases_df_cur[as.character(cases_df_cur$intervention_arm) != 'No intervention',]
  cases_df_cur = cases_df_cur[c('intervention_arm', 'trial_year_diff', 'treated_clinical_cases')]
  # subset to relevant part of data frame - aggregated years
  cases_df_cur_all = cases_df_all_years_ave[intersect(intersect(which(cases_df_all_years_ave$seasonality == seasonality_options[ss]),
                                                                which(cases_df_all_years_ave$Annual.EIR == eir_cur)),
                                                      which(cases_df_all_years_ave$cm_coverage == cm_cur)),]
  cases_df_cur_all = cases_df_cur_all[as.character(cases_df_cur_all$intervention_arm) != 'No intervention',]
  cases_df_cur_all = cases_df_cur_all[c('intervention_arm', 'treated_clinical_cases')]
  cases_df_cur_all$trial_year_diff = 'Aggregated'

  # each intervention arm in a column
  cases_df_wide = reshape2::dcast(cases_df_cur, trial_year_diff ~ intervention_arm, value.var = 'treated_clinical_cases')
  cases_df_all_wide = reshape2::dcast(cases_df_cur_all, trial_year_diff ~ intervention_arm, value.var = 'treated_clinical_cases')
  cases_df_wide_2 = merge(cases_df_wide, cases_df_all_wide, by = c('trial_year_diff', 'Chemoprevention Alone', 'Vaccine Alone', 'Combination'), all = TRUE)

  # calculate rate ratios
  cases_df_wide_2$rtss_over_smc = cases_df_wide_2$`Vaccine Alone` / cases_df_wide_2$`Chemoprevention Alone`
  cases_df_wide_2$combination_over_rtss = cases_df_wide_2$Combination / cases_df_wide_2$`Vaccine Alone`
  cases_df_wide_2$combination_over_smc = cases_df_wide_2$Combination / cases_df_wide_2$`Chemoprevention Alone`

  rate_ratios = reshape2::melt(cases_df_wide_2[c('trial_year_diff', 'rtss_over_smc', 'combination_over_rtss', 'combination_over_smc')], id.vars = 'trial_year_diff')
  rate_ratios$trial_year = rate_ratios$trial_year_diff
  rate_ratios$trial_year[rate_ratios$trial_year %in% paste0(1:10)] = paste0('Year ', rate_ratios$trial_year[rate_ratios$trial_year %in% paste0(1:10)])

  # merge with reference dataset
  rate_ratios = merge(rate_ratios, ref_df, by = c('trial_year', 'variable'))

  # plot rate ratios
  gg = ggplot(rate_ratios, aes(x = variable, y = value)) +
    geom_point(color = rgb(0.6, 0.8, 1), pch = 17, size = 3) +
    geom_point(aes(x = variable, y = ref_incidence_rate_ratio), pch = 20, cex = 2) +
    theme_light() +
    geom_hline(yintercept = 1) +
    ylim(0.05, max(max(rate_ratios$value), 1.35)) +
    ylab('Incidence rate ratio') +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title.x = element_blank()) + #+, strip.background = element_blank(), axis.line=element_line(), legend.position='none') +
    facet_wrap(facets = vars(trial_year), nrow = 1)
  # facet_wrap(facets=vars(intervention_arm), ncol=1, scales = "free_x")
  f_save_plot(gg, paste0('incidenceRateRatios_EIR',eir_cur, '_', seasonality_options[ss], '_cm', cm_cur * 100),
              file.path(exp_filepath), width = 7.5, height = 2.75, units = 'in', device_format = device_format)

  print(gg)
}
