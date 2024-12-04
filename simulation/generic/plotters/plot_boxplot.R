library(data.table)
library(dplyr)
library(ggplot2)

#wdir at rtss-scenarios
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path("simulation", "load_paths.R"))

simout_dir <- file.path(projectpath, 'simulation_output', 'generic_forward')
theme_set(theme_bw())
customTheme <- f_getCustomTheme()


plot_boxplots <- function(dat, xvar = 'Annual_EIR', yvar = 'cases_averted_per100000',
                          fillvar = 'rtss_coverage', facet1 = 'smc_coverage', facet2 = 'seasonality',
                          Uage = 'U5', SAVE = TRUE, caption = '', add_watermark = TRUE) {

  dat <- as.data.frame(dat)
  dat$facet1 = dat[, facet1]
  dat$facet2 = dat[, facet2]
  #dat$facet1 = paste0(facet1, '\n', dat[, facet1])
  #dat$facet2 = paste0(facet2, '\n', dat[, facet2])

  pplot <- ggplot(data = dat) +
    geom_boxplot(aes(x = as.factor(get(xvar)), y = get(yvar), fill = as.factor(get(fillvar)),
                     group = interaction(get(xvar), get(fillvar))),
                 position = position_dodge(width = 1),
                 outlier.shape = NA) +
    scale_fill_brewer(palette = "YlGnBu",
                      guide = guide_legend(reverse = TRUE)) +
    labs(x = xvar, y = yvar, fill = fillvar,
         caption = caption) +
    customTheme +
    theme(panel.spacing = unit(1, "lines"))

  if (!is.null(facet1) & !is.null(facet2)) pplot <- pplot + facet_grid(facet1 ~ facet2, scales = "free")
  if (!is.null(facet1) & is.null(facet2)) pplot <- pplot + facet_wrap(~facet1, scales = "free")

  if (add_watermark)pplot <- watermark(pplot)

  if (SAVE) {
    if (!dir.exists(file.path(simout_dir, exp_name, '_plots', Uage))) {
      dir.create(file.path(simout_dir, exp_name, '_plots', Uage))
    }
    f_save_plot(pplot, plot_name = paste('boxplot', yvar, facet1, sep = '_'),
                width = 8, height = 6,
                plot_dir = file.path(simout_dir, exp_name, '_plots', Uage))

  }
  return(pplot)
}

#_Run_plots______________________________________________________#
#exp_name <- "generic_settings_test_SMC_constrained"
exp_name <- "generic_campboost_RTSS_testrunv0"

exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'cm_coverage', 'ipti_coverage', 'rtss_coverage', 'smc_coverage',
                'intervention_correlation', 'frac_high_access', 'Cohort_birth_month', 'rtss_mode')

plot_dir <- file.path(simout_dir, exp_name, '_plots')
if (!dir.exists(plot_dir))dir.create(plot_dir)

cases_df_list <- load_Age_monthly_Cases(simout_dir, exp_name, exp_sweeps) #keep_birth_month=T
dat <- cases_df_list[[2]]

if ('age_group' %in% colnames(dat)) {
  table(dat$age_group)
  dat$age_group <- factor(dat$age_group,
                          levels = c('U1', 'U2', 'U5', 'U10'),
                          labels = c('U1', 'U2', 'U5', 'U10'))

  cases_outcomes <- c('protective_efficacy', 'relative_burden', 'cases_averted_per100000')
  severe_cases_outcomes <- c('protective_efficacy_severe', 'relative_burden_severe',
                             'severe_cases_averted_per100000')
}else {
  dat$age_group <- dat$year
  lvls <- paste0(unique(dat$year), 'to', unique(dat$year) + 1)
  dat$age_group <- factor(dat$age_group, levels = unique(dat$year), labels = lvls)
  #table(dat$age_group, dat$year, exclude = NULL)
  cases_outcomes <- c('clinical_cases')
  severe_cases_outcomes <- c('severe_cases')
}


plot_boxplots(subset(dat), xvar = 'smc_coverage', yvar = cases_outcomes[1],
              fillvar = 'rtss_coverage', facet1 = 'age_group', facet2 = NULL,
              Uage = 'all', SAVE = T)





for (ageGrp in unique(dat$age_group)) {
  for (outcome in c(cases_outcomes, severe_cases_outcomes)) {
    if (outcome %in% colnames(dat)) {
      plot_boxplots(subset(dat, age_group == ageGrp), xvar = 'smc_coverage', yvar = outcome,
                    fillvar = 'rtss_coverage', facet1 = NULL, facet2 = NULL,
                    Uage = ageGrp, SAVE = TRUE)

      plot_boxplots(subset(dat, age_group == ageGrp), xvar = 'smc_coverage', yvar = outcome,
                    fillvar = 'rtss_coverage', facet1 = 'Annual_EIR', facet2 = NULL,
                    Uage = ageGrp, SAVE = TRUE)
    }else {
      print(paste0(outcome, " not found in dataframe"))
    }
  }
}
