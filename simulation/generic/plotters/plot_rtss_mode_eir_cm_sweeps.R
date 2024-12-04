library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)

source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path('simulation', 'generic', 'plotters', 'plot_bars_and_lines.R'))
source(file.path("simulation", "load_paths.R"))

simout_dir <- file.path(projectpath, 'simulation_output', 'generic_forward')
theme_set(theme_bw())
customTheme <- f_getCustomTheme()
myColors <- f_getColors()
max_years = c(1, 2, 5, 10)
age_groups = paste0('U', max_years)
add_watermark = FALSE

# subset for plot
cur_cm = 0.6
cur_corr = 0
cur_rtss = c(0.6, 0.8)
cur_smc = c(0.6)

## Age labels
age_label_values = seq(0, 8, length.out = 5)
age_labels = paste0(age_label_values, '-', (age_label_values + 1))

## Simulation experiment
exp_name <- 'generic_campboost3_SMC'
exp_filepath <- file.path(simout_dir, exp_name)
plot_dir <- file.path(exp_filepath, '_plots')
if (!dir.exists(plot_dir))dir.create(plot_dir)

exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'cm_coverage', 'rtss_coverage', 'smc_coverage',
                'intervention_correlation', 'Cohort_birth_month', 'rtss_mode')
sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)

generate_lineplot_over_age = TRUE
generate_cumulative_barplot = TRUE

##__________________________________
### Lineplot over age
##__________________________________
if (generate_lineplot_over_age) {
  dat = sim_output[[4]]

  dat$cm_coverage <- as.numeric(dat$cm_coverage)
  dat$smc_coverage <- as.numeric(dat$smc_coverage)
  dat$intervention_correlation <- as.numeric(dat$intervention_correlation)
  dat$age_group = factor(dat$age_group, levels = c("U1", "U2", "U5", "U10"))
  dat$seasonality = factor(dat$seasonality, levels = c('constant', 'moderate_unimodal', 'high_unimodal'))
  dat$rtss_mode_fct = factor(dat$rtss_mode,
                             levels = c('constant', 'campboost3', 'campboost3B'),
                             labels = c('1 booster all at age 24-35 mths',
                                        '1 booster spread out between age 24-59 mths',
                                        '1 booster spread out between age 24-59 mths\n+2nd booster between 48-59 mths'))

  dat$rtss_coverage_fct <- factor(dat$rtss_coverage,
                                  levels = unique(dat$rtss_coverage),
                                  labels = paste0('RTS,S ', unique(dat$rtss_coverage) * 100, '%'))
  dat$smc_coverage_fct <- factor(dat$smc_coverage,
                                 levels = unique(dat$smc_coverage),
                                 labels = paste0('SMC ', unique(dat$smc_coverage) * 100, '%'))
  dat$Annual_EIR_fct <- factor(dat$Annual_EIR,
                               levels = unique(dat$Annual_EIR),
                               labels = paste0('annual EIR ', unique(dat$Annual_EIR)))


  outcomes <- colnames(dat)[grep("rtss", colnames(dat))]
  outcomes <- outcomes[!(outcomes %in% c("rtss_coverage", "severe_cases_no_rtss", "rtss_mode_fct",
                                         "rtss_coverage_fct", "clinical_cases_no_rtss", "rtss_mode"))]

  for (outcome in outcomes) {
    #outcome=outcomes[1]
    if ((outcome) %in% colnames(dat)) {
      print(paste0('generate plot for ', outcome))
      gg <- dat %>%
        filter(
          intervention_correlation == cur_corr,
          cm_coverage %in% cur_cm,
          rtss_coverage %in% cur_rtss,
          smc_coverage %in% cur_smc) %>%
        ggplot(aes(x = year, y = get(outcome))) +
        geom_point(aes(col = rtss_mode_fct, group = rtss_mode_fct), size = 2) +
        geom_line(aes(col = rtss_mode_fct, group = rtss_mode_fct), size = 1) +
        theme_bw() +
        f_getCustomTheme() +
        scale_color_manual(values = myColors) +
        scale_x_continuous(breaks = age_label_values, labels = age_labels) +
        geom_hline(yintercept = 0) +
        labs(y = gsub("_", " ", gsub("rtss_", "", outcome)),
             x = 'age of child',
             col = 'enhanced RTS,S') +
        facet_wrap(Annual_EIR_fct ~ rtss_coverage_fct)

      if (add_watermark)gg <- watermark(gg, hjust = 1.5, vjust = -3)
      if (!dir.exists(file.path(plot_dir, 'lineplots')))dir.create(file.path(plot_dir, 'lineplots'))
      f_save_plot(gg, plot_name = paste0(outcome, ''),
                  width = 12, height = 6, plot_dir = file.path(plot_dir, 'lineplots'))
    }else {
      print(paste0(outcome, ' not found in dataframe'))
    }
  }
}

##__________________________________
### Cumulative barplot
##__________________________________
if (generate_cumulative_barplot) {
  dat <- sim_output[[3]]

  dat$scenario_name = paste0(round(100 * dat$rtss_coverage), '% RTS,S; ', round(100 * dat$smc_coverage), '% SMC')
  dat$cm_coverage <- as.numeric(dat$cm_coverage)
  dat$smc_coverage <- as.numeric(dat$smc_coverage)
  dat$intervention_correlation <- as.numeric(dat$intervention_correlation)

  dat$seasonality = factor(dat$seasonality, levels = c('constant','moderate_unimodal','high_unimodal'))
  dat$rtss_mode_fct = factor(dat$rtss_mode, levels = c('constant', 'campboost3', 'campboost3B'),
                             labels = c('1 booster all at age 24-35 mths',
                                        '1 booster spread out between age 24-59 mths',
                                        '1 booster spread out between age 24-59 mths\n+2nd booster between 48-59 mths'))

  dat$rtss_coverage_fct <- factor(dat$rtss_coverage,
                                  levels = unique(dat$rtss_coverage),
                                  labels = paste0('RTS,S ', unique(dat$rtss_coverage) * 100, '%'))
  dat$smc_coverage_fct <- factor(dat$smc_coverage, levels = unique(dat$smc_coverage),
                                 labels = paste0('SMC ', unique(dat$smc_coverage) * 100, '%'))


  dat$Annual_EIR_fct <- factor(dat$Annual_EIR, levels = unique(dat$Annual_EIR),
                               labels = paste0('annual EIR ', unique(dat$Annual_EIR)))

  dat$age_group = factor(dat$age_group, levels = c("U1", "U2", "U5", "U10"),
                         labels = c("U1", "U2", "U5", "U10"))

  for (Uage in c("U1", "U2", "U5", "U10")) {
    # Uage = "U5"
    for (outcome in outcomes) {
      #outcome = outcomes[1]
      gg <- dat %>%
        filter(age_group == Uage) %>%
        filter(
          intervention_correlation == cur_corr,
          cm_coverage %in% cur_cm,
          rtss_coverage %in% cur_rtss,
          smc_coverage %in% cur_smc) %>%
        ggplot() +
        geom_bar(aes(x = as.factor(as.numeric(rtss_mode_fct)),
                     y = get(outcome), fill = rtss_mode_fct), stat = 'identity') +
        labs(y = gsub("_", " ",
                      gsub("rtss_", "", outcome)),
             x = 'enhanced RTS,S',
             fill = 'enhanced RTS,S') +
        facet_wrap(Annual_EIR_fct ~ rtss_coverage_fct) +
        theme_bw() +
        f_getCustomTheme() +
        scale_fill_manual(values = myColors) +
        customTheme

      if (!dir.exists(file.path(plot_dir, 'barplots')))dir.create(file.path(plot_dir, 'barplots'))
      f_save_plot(gg, plot_name = paste0(Uage, '_', outcome, ''),
                  width = 12, height = 6, plot_dir = file.path(plot_dir, 'barplots'))
    }
  }
}

