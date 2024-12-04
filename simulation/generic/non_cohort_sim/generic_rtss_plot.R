library(data.table)
library(dplyr)
library(ggplot2)

#wdir at rtss-scenarios
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path("simulation", "load_paths.R"))

source(file.path("simulation", "load_paths.R"))
simout_dir <- file.path(projectpath, 'simulation_output', 'generic')
theme_set(theme_bw())
customTheme <- f_getCustomTheme()


plot_timeline <- function(plotdat, Uage = "U1", SAVE = TRUE) {
  pplot <- ggplot(data = plotdat) +
    geom_line(aes(x = year, y = incidence, col = as.factor(coverage))) +
    facet_wrap(~Annual_EIR + cm_cov_U5, nrow = 4, scales = "free",
               labeller = labeller(.cols = label_both)) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = "", x = "", y = paste0("incidence ", Uage),
         color = "RTS,S coverage", caption = 'RTS,S simulated with booster') +
    customTheme

  if (SAVE) {
    f_save_plot(pplot, plot_name = paste0(Uage, '_indicence_timeline'), width = 12, height = 8,
                plot_dir = file.path(simout_dir, exp_name, '_plots'))
  }
  return(pplot)


}

plot_effect_size <- function(dat, Uage, selected_year = 2023, SAVE = TRUE) {

  dat <- dat %>%
    ungroup() %>%
    filter(year == selected_year) %>%
    dplyr::select(year, coverage, Annual_EIR, cm_cov_U5, incidence)

  dat_counterfactual <- dat %>%
    filter(coverage == 0) %>%
    filter(year == selected_year) %>%
    dplyr::select(year, Annual_EIR, cm_cov_U5, incidence) %>%
    rename(incidence_counterfactual = incidence)

  plotdat <- dat %>%
    left_join(dat_counterfactual) %>%
    group_by(coverage, Annual_EIR, cm_cov_U5) %>%
    mutate(eff = 1 - (incidence / incidence_counterfactual))

  plotdat$cm_cov_U5 <- as.factor(plotdat$cm_cov_U5)
  pplot <- ggplot(data = plotdat) +
    geom_bar(aes(x = cm_cov_U5, y = eff, fill = as.factor(coverage),
                 group = interaction(coverage, cm_cov_U5)),
             position = position_dodge(width = 1), stat = 'identity') +
    #geom_boxplot(aes(x = cm_cov_U5, y = eff, fill = as.factor(coverage),
    #                 group = interaction(coverage, cm_cov_U5)),
    #             position = position_dodge(width = 1),
    #             outlier.shape = NA) +
    #geom_jitter(aes(x = cm_cov_U5, y = eff, group = as.factor(coverage)),
    #            position = position_dodge(width = 1),
    #            col = 'white', fill = 'black', shape = 21, alpha = 0.7) +
    facet_wrap(~Annual_EIR, scales = "free", labeller = labeller(.cols = label_both)) +
    scale_fill_brewer(palette = "Dark2") +
    labs(title = paste0("Effect size of RTS,S (EPI deployment) in ", Uage, "\n one year after intervention start"),
         x = "Case management",
         y = 'Effect size\n(relative reduction to no RTS,S)',
         fill = "RTS,S coverage",
         caption = "Simulated with booster, fixed coverage 80%") +
    theme(panel.spacing = unit(1, "lines")) +
    customTheme

  if (SAVE) {
    f_save_plot(pplot, plot_name = paste0(Uage, '_indicence_effectsize'), width = 12, height = 8,
                plot_dir = file.path(simout_dir, exp_name, '_plots'))

  }
  return(pplot)
}


#_Run_plots______________________________________________________#
exp_name <- "MR_generic_rtss_sweep_nobooster"
plot_dir <- file.path(simout_dir, exp_name, '_plots')
if (!dir.exists(plot_dir))dir.create(plot_dir)

simdat <- load_sim_data(simout_dir, exp_name)
U1dat <- simdat[[1]] %>% filter(year > 2015)
U5dat <- simdat[[2]] %>% filter(year > 2015)
rm(simdat)

plot_timeline(U1dat, "U1")
plot_timeline(U5dat, "U5")

plot_effect_size(U1dat, "U1", selected_year = 2025)
plot_effect_size(U5dat, "U5", selected_year = 2025)
