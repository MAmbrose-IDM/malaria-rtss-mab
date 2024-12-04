library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(cowplot)

source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path("simulation", "load_paths.R"))
cols <- colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[3:12]

theme_set(theme_bw())
customTheme <- f_getCustomTheme()

exp_name <- "generic_eir_cm_sweep_no_RTSS"
exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'cm_coverage', 'ipti_coverage', 'rtss_coverage', 'smc_coverage',
                'intervention_correlation', 'frac_high_access', 'Cohort_birth_month', 'rtss_mode')
simout_dir <- file.path(projectpath, 'simulation_output', 'eir_pfpr')

cases_df_list <- load_Age_monthly_Cases(simout_dir, exp_name, exp_sweeps)
dat <- cases_df_list[[2]] %>%
  filter(year <= 5) %>%   #Children aged U5 years, including first 3 mths
  filter(intervention_correlation == 0 &
           smc_coverage == 0 &
           rtss_coverage == 0) %>%
  group_by(Annual_EIR, seasonality, cm_coverage) %>%
  summarise(clinical_cases = sum(clinical_cases),
            PfPR_U5 = mean(pfpr) * 100) %>%
  mutate(clinical_cases = clinical_cases / 5)

dat$cm_coverage_fct <- as.factor(dat$cm_coverage)
dat$PfPR_U5_label <- paste0(round(dat$PfPR_U5, 0), '%')
dat$clinical_cases_label <- paste0(round(dat$clinical_cases, 0))
dat$seasonality <- factor(dat$seasonality,
                          levels = c("constant", "moderate_unimodal", "high_unimodal"),
                          labels = c("constant", "moderate_unimodal", "high_unimodal"))

### Add average across CM values and seasonalites, keep unaggregated data
dat <- dat %>%
  group_by(Annual_EIR) %>%
  mutate(PfPR_U5_mean = round(mean(PfPR_U5), 1),
         clinical_cases_mean = round(mean(clinical_cases), 0))
fwrite(dat, file = file.path(simout_dir, exp_name, 'eir_U5_pfpr_dat.csv'))

## Save another copy in simout_dir
dat %>%
  dplyr::select(Annual_EIR, seasonality, cm_coverage, PfPR_U5, clinical_cases, PfPR_U5_mean, clinical_cases_mean) %>%
  fwrite(file = file.path(simout_dir, 'generic_eir_U5_pfpr_table.csv'))

dat <- dat %>% filter(seasonality != "moderate_unimodal")

f_custom_plot <- function(xvar, yvar) {
  gg <- ggplot(data = dat, aes(x = get(xvar), y = get(yvar),
                               shape = seasonality, linetype = seasonality, col = cm_coverage_fct)) +
    geom_point(size = 2) +
    geom_line() +
    scale_color_manual(values = cols[3:10], guide = guide_legend(reverse = TRUE)) +
    scale_fill_manual(values = cols[3:10], guide = guide_legend(reverse = TRUE)) +
    labs(x = gsub("_", " ", xvar),
         y = gsub("_", " ", yvar),
         fill = 'Case management', col = 'Case management') +
    scale_y_continuous(labels = comma) +
    customTheme
  return(gg)
}

pplot1 <- f_custom_plot('Annual_EIR', 'PfPR_U5') + labs(title = "PfPR-EIR\n", y = 'PfPR (HRP2) U5')
pplot2 <- f_custom_plot('Annual_EIR', 'clinical_cases') + labs(title = "Incidence-EIR\n", y = 'clinical cases U5')
pplot3 <- f_custom_plot('PfPR_U5', 'clinical_cases') + labs(title = "Incidence-PfPR\n", y = 'clinical cases U5', x = 'PfPR (HRP2) U5')

### Combine and save
plegend <- get_legend(pplot1)
pplot1 <- pplot1 + theme(legend.position = "None")
pplot2 <- pplot2 + theme(legend.position = "None")
pplot3 <- pplot3 + theme(legend.position = "None")
pplot <- plot_grid(pplot1, pplot2, pplot3, nrow = 1)
pplot <- plot_grid(pplot, plegend, nrow = 1, rel_widths = c(1, 0.25))
print(pplot)
f_save_plot(pplot, plot_name = paste0('combined_EIR_PfPR_incidence_curves'),
            width = 14, height = 4,
            plot_dir = file.path(simout_dir, exp_name))
