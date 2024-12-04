library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
#wdir at rtss-scenarios
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path("simulation", "load_paths.R"))

simout_dir <- file.path(projectpath, 'simulation_output', 'generic_forward')
theme_set(theme_minimal())
customTheme <- f_getCustomTheme()

agegrp_labels <- c("U1", "U2", "U5", "U10")
exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'cm_coverage', 'ipti_coverage', 'rtss_coverage', 'smc_coverage',
                'intervention_correlation', 'frac_high_access', 'Cohort_birth_month', 'rtss_mode')
keepVARS <- c('age_group', 'Annual_EIR', 'seasonality', 'intervention_correlation', 'protective_efficacy',
              'cm_coverage', 'rtss_coverage', 'smc_coverage')#'rtss_mode'

cases_outcomes <- c('protective_efficacy', 'relative_burden', 'cases_averted_per100000')
severe_cases_outcomes <- c('protective_efficacy_severe', 'relative_burden_severe', 'severe_cases_averted_per100000')


##--------------------------------------
## Define custom functions
##--------------------------------------
get_interpolatedgrid <- function(dat, cov1, cov2, outcome) {
  ### Define regression model per group
  dfLM <- dat %>%
    as.data.table() %>%
    rename(cov1 = as.name(cov1),
           cov2 = as.name(cov2)) %>%
    dplyr::group_by(facet1, facet2) %>%
    do(fitlm = loess(get(outcome) ~ cov1 + cov2, data = .))

  ### Make predictions
  pred_list <- list()
  for (i in c(1:nrow(dfLM))) {
    pred_dat <- data.frame(expand.grid('cov1' = seq(0, 1, 0.05), 'cov2' = seq(0, 1, 0.05)))
    pred_dat$y_pred <- NA

    model <- dfLM[i, 'fitlm'][[1]]
    pred_dat$y_pred <- predict(model[[1]], newdata = pred_dat, type = "response")
    pred_dat$facet1 <- dfLM[i, 'facet1'][[1]]
    pred_dat$facet2 <- dfLM[i, 'facet2'][[1]]
    pred_dat[, outcome] <- pred_dat[, 'y_pred']
    pred_dat[, cov1] <- pred_dat[, 'cov1']
    pred_dat[, cov2] <- pred_dat[, 'cov2']
    pred_dat <- pred_dat %>% dplyr::select(-cov1, -cov2, -y_pred)
    pred_list[[length(pred_list) + 1]] <- pred_dat
  }

  ### Combine data-list
  pred_dat <- pred_list %>% bind_rows()
  return(pred_dat)
}

plot_heatmap <- function(dat, cov1 = 'smc_coverage', cov2 = 'rtss_coverage', outcome = 'protective_efficacy',
                         facet1 = 'Annual_EIR', facet2 = 'cm_coverage', Uage = 'U5', SAVE = TRUE, smoothgrid = FALSE,
                         add_watermark = FALSE, plot_dir = NULL, smoothgrid_scatter = F) {

  if (is.null(plot_dir))plot_dir <- file.path(simout_dir, exp_name, '_plots')
  # Facets wit labels
  dat$facet1 = dat[, facet1]
  dat$facet2 = dat[, facet2]

  if (smoothgrid) {
    dat_raw <- dat
    dat <- get_interpolatedgrid(dat_raw, cov1, cov2, outcome)
    plot_dir <- file.path(simout_dir, exp_name, '_plots_smooth')

    if (smoothgrid_scatter) {
      dat_raw_i <- dat_raw %>% select_at(.vars = c(outcome, 'facet1', 'facet2', 'smc_coverage', 'rtss_coverage'))
      colnames(dat_raw_i)[colnames(dat_raw_i) == outcome] <- paste0(outcome, '_raw')
      cordat <- dat %>%
        filter(smc_coverage %in% dat_raw$smc_coverage & rtss_coverage %in% dat_raw$rtss_coverage) %>%
        left_join(dat_raw_i)
      corval <- as.numeric(cor.test(cordat[, outcome], cordat[, paste0(outcome, '_raw')], method = "spearman")$estimate)

      pfit_scatterplot <- dat %>%
        filter(smc_coverage %in% dat_raw$smc_coverage & rtss_coverage %in% dat_raw$rtss_coverage) %>%
        left_join(dat_raw_i) %>%
        ggplot() +
        geom_point(aes(x = get(paste0(outcome, '_raw')), y = get(outcome), group = interaction(smc_coverage, rtss_coverage))) +
        geom_smooth(aes(x = get(paste0(outcome, '_raw')), y = get(outcome)), se = FALSE, col = 'blue') +
        geom_abline(intercept = 0, slope = 1) +
        facet_grid(facet1 ~ facet2)
      labs(x = paste0(outcome, '_raw'), y = outcome, caption = paste0("Spearman's rho", corval))

      f_save_plot(pfit_scatterplot, plot_name = paste('fit_scatterplot', outcome,
                                                      gsub('_coverage', '', cov1),
                                                      gsub('_coverage', '', cov2),
                                                      facet1, facet2, sep = "_"),
                  width = 8, height = 6,
                  plot_dir = file.path(plot_dir, Uage))
    }
  }

  dat <- as.data.table(dat)
  stepsize = 0.2

  pplot <- ggplot(data = dat, aes(x = get(cov1), y = get(cov2), z = get(outcome), fill = get(outcome))) +
    geom_tile() +
    labs(x = gsub("_", " ", cov1), y = gsub("_", " ", cov2),
         fill = gsub("_", " ", outcome)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, stepsize)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1, stepsize)) +
    customTheme +
    theme(panel.spacing = unit(0.9, "lines"))

  if (smoothgrid)pplot <- pplot + geom_contour(col = "black")

  if (sum(grep('protective', outcome)) > 0) {
    pplot <- pplot + scale_fill_viridis_c(option = "A", breaks = seq(0, 1, 0.2), lim = c(0, 1))
  }else {
    pplot <- pplot + scale_fill_viridis_c(option = "A")

  }

  if (!is.null(facet1)) {
    if (length(facet1) == 1) {
      if (!is.null(facet2)) pplot <- pplot + facet_grid(facet1 ~ facet2)
      if (is.null(facet2)) pplot <- pplot + facet_wrap(~facet1)
      facet1_label = facet1
    }else {
      if (!is.null(facet2)) pplot <- pplot + facet_grid(facet1a + facet1b ~ facet2)
      if (is.null(facet2)) pplot <- pplot + facet_wrap(~facet1a + facet1b)
      facet1_label = paste(facet1a, facet1b, sep = '_')
    }
  }

  if (add_watermark)pplot <- watermark(pplot)

  if (SAVE) {
    if (!dir.exists(file.path(plot_dir))) dir.create(file.path(plot_dir))
    if (!dir.exists(file.path(plot_dir, Uage))) dir.create(file.path(plot_dir, Uage))
    f_save_plot(pplot, plot_name = paste('heatmap', outcome,
                                         gsub('_coverage', '', cov1),
                                         gsub('_coverage', '', cov2),
                                         facet1_label, facet2, sep = "_"),
                width = 8, height = 6,
                plot_dir = file.path(plot_dir, Uage))

  }
  return(pplot)

}


##--------------------------------------
## Load simulation data and generate figures
##--------------------------------------
#exp_name <- "generic_heatmap_constant_cor0_SMC"
#exp_name <- "generic_heatmap_constant_cor1_SMC"
#exp_name <- "generic_heatmap_highseason_cor1_SMC"
#exp_name <- "generic_heatmap_highseason_cor0_SMC"
exp_name <- "generic_heatmap_highseason_RTSSbooster_cor0_SMC"

dat <- f_load_sim_PEestimates(simout_dir, exp_name, exp_sweeps)

#table(dat$age_group)
dat$age_group <- factor(dat$age_group, levels = agegrp_labels, labels = agegrp_labels)

plegend <- get_legend(plot_heatmap(subset(dat), Uage = 'all',
                                   outcome = cases_outcomes[1],
                                   cov1 = 'smc_coverage', cov2 = 'rtss_coverage',
                                   facet1 = 'Annual_EIR', facet2 = 'age_group',
                                   SAVE = F, smoothgrid = T))

### Heatmap across all U-ages
for (outc in c(cases_outcomes, severe_cases_outcomes)) {
  facet1 <- 'Annual_EIR'
  print(outc)
  plot_heatmap(subset(dat, cm_coverage == 0.6), Uage = 'all',
               outcome = outc,
               cov1 = 'smc_coverage', cov2 = 'rtss_coverage',
               facet1 = facet1, facet2 = 'age_group',
               SAVE = T, smoothgrid = T, smoothgrid_scatter = T)
}

### Custom plot,
dat <- dat %>% filter(age_group %in% c("U2", "U5"), cm_coverage == 0.6)

facet1 <- 'Annual_EIR' # 'Annual_EIR'
p1 <- plot_heatmap(subset(dat, cm_coverage == 0.6), Uage = 'all',
                   outcome = cases_outcomes[1],
                   cov1 = 'smc_coverage', cov2 = 'rtss_coverage',
                   facet1 = facet1, facet2 = 'age_group',
                   SAVE = F, smoothgrid = T) +
  labs(title = 'Clinical cases') +
  theme(legend.position = 'none')

p2 <- plot_heatmap(subset(dat), Uage = 'all',
                   outcome = severe_cases_outcomes[1],
                   cov1 = 'smc_coverage', cov2 = 'rtss_coverage',
                   facet1 = 'Annual_EIR', facet2 = 'age_group',
                   SAVE = F, smoothgrid = T) +
  labs(title = 'Severe cases') +
  theme(legend.position = 'none')


pplot <- plot_grid(p1, p2, nrow = 1)
pplot <- plot_grid(pplot, plegend, nrow = 1, rel_widths = c(1, 0.25))
print(pplot)

f_save_plot(pplot, plot_name = paste0('heatmap_PE_custom_smooth'),
            width = 12, height = 6,
            plot_dir = file.path(simout_dir, exp_name))
rm(plegend, p1, p2, pplot)

#### Without smoothing
plegend <- get_legend(plot_heatmap(subset(dat), Uage = 'all',
                                   outcome = cases_outcomes[1],
                                   cov1 = 'smc_coverage', cov2 = 'rtss_coverage',
                                   facet1 = 'Annual_EIR', facet2 = 'age_group',
                                   SAVE = F, smoothgrid = F))

p1 <- plot_heatmap(subset(dat), Uage = 'all',
                   outcome = cases_outcomes[1],
                   cov1 = 'smc_coverage', cov2 = 'rtss_coverage',
                   facet1 = 'Annual_EIR', facet2 = 'age_group',
                   SAVE = F, smoothgrid = F) +
  labs(title = 'Clinical cases') +
  theme(legend.position = 'none')

p2 <- plot_heatmap(subset(dat), Uage = 'all',
                   outcome = severe_cases_outcomes[1],
                   cov1 = 'smc_coverage', cov2 = 'rtss_coverage',
                   facet1 = 'Annual_EIR', facet2 = 'age_group',
                   SAVE = F, smoothgrid = F) +
  labs(title = 'Severe cases') +
  theme(legend.position = 'none')


pplot <- plot_grid(p1, p2, nrow = 1)
pplot <- plot_grid(pplot, plegend, nrow = 1, rel_widths = c(1, 0.25))

f_save_plot(pplot, plot_name = paste0('heatmap_PE_custom'),
            width = 12, height = 6,
            plot_dir = file.path(simout_dir, exp_name))
rm(plegend, p1, p2, pplot)


##---------------------------------------
## Comparison by defined parameter
##---------------------------------------
EXP_COMPARISON <- 'seasonality'
#EXP_COMPARISON <- 'intervention_correlation'
#EXP_COMPARISON <- 'rtss_mode'

if (EXP_COMPARISON == 'seasonality') {
  ## High versus constant seasonality
  exp_names <- c("generic_heatmap_highseason_cor0_SMC", "generic_heatmap_constant_cor0_SMC")
  exp_name <- exp_names[1]
}
if (EXP_COMPARISON == 'intervention_correlation') {
  ## Wit/without access correlation
  #exp_names <- c("generic_heatmap_highseason_cor0_SMC","generic_heatmap_highseason_cor1_SMC")
  exp_names <- c("generic_heatmap_constant_cor0_SMC", "generic_heatmap_constant_cor1_SMC")
  #exp_names <- c("generic_accesscorrelation_SMC")
  exp_name <- exp_names[1]
}
if (EXP_COMPARISON == 'rtss_mode') {
  ##  Standard vs enhanced RTSS
  exp_names <- c("generic_heatmap_highseason_cor0_SMC", "generic_heatmap_highseason_RTSSbooster_cor0_SMC")
  exp_name <- exp_names[1]
}

if (length(exp_names) == 2) {
  dat1 <- f_load_sim_PEestimates(simout_dir, exp_names[1], exp_sweeps) %>%
    mutate(expname = gsub("heatmap_", "", gsub("generic_", "", exp_names[1]))) %>%
    dplyr::select_at(.vars = c(keepVARS, 'expname'))
  dat2 <- f_load_sim_PEestimates(simout_dir, exp_names[2], exp_sweeps) %>%
    mutate(expname = gsub("heatmap_", "", gsub("generic_", "", exp_names[1]))) %>%
    dplyr::select_at(.vars = c(keepVARS, 'expname'))
  dat <- as.data.frame(rbind(dat1, dat2))
}
if (length(exp_names) == 1) {
  dat <- f_load_sim_PEestimates(simout_dir, exp_name, exp_sweeps) %>%
    mutate(expname = gsub("heatmap_", "", gsub("generic_", "", exp_names[1]))) %>%
    dplyr::select_at(.vars = c(keepVARS, 'expname'))
}

dat$age_group <- factor(dat$age_group, levels = agegrp_labels, labels = agegrp_labels)
dat <- subset(dat, Annual_EIR == 10 & cm_coverage == 0.6)
dat <- subset(dat, rtss_coverage %in% c(0, 0.2, 0.4, 0.8))
dat <- subset(dat, age_group %in% c("U5"))

#(grp_var <- EXP_COMPARISON)
grp_var = 'expname'
pplot <- ggplot(data = subset(dat, seasonality != 'constant')) +
  geom_point(aes(x = smc_coverage, y = get(cases_outcomes[1]), shape = get(grp_var),
                 col = get(grp_var)), size = 2) +
  geom_line(aes(x = smc_coverage, y = get(cases_outcomes[1]), linetype = get(grp_var),
                col = get(grp_var))) +
  facet_wrap(~rtss_coverage, nrow = 1) +
  scale_y_continuous(lim = c(0, 0.5)) +
  scale_linetype_manual(values = c('dashed', 'solid')) +
  scale_color_manual(values = c('grey', 'maroon4')) +
  labs(subtitle = 'Facets by RTS,S coverage', x = "SMC coverage", y = cases_outcomes[1],
       col = grp_var, shape = grp_var, linetype = grp_var) +
  customTheme

print(pplot)
f_save_plot(p2, plot_name = paste0('coverage_', EXP_COMPARISON, '_U5'),
            width = 12, height = 4,
            plot_dir = file.path(simout_dir, exp_name))
rm(pplot, dat, dat1, dat2, exp_names, exp_name)

