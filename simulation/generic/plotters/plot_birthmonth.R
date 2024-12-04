library(data.table)
library(dplyr)
library(ggplot2)

#wdir at rtss-scenarios
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path("simulation", "load_paths.R"))

simout_dir <- file.path(projectpath, 'simulation_output', 'generic_forward')
theme_set(theme_bw())
customTheme <- f_getCustomTheme()


#_Run_plots______________________________________________________#
#exp_name <- "generic_settings_test_SMC_constrained"
#exp_name <- "generic_heatmap_highseason_SMCv0"
#exp_name <- "generic_campboost3_double_RTSS_test"
exp_name <- "generic_campboost3_SMC"

exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'cm_coverage', 'ipti_coverage', 'rtss_coverage', 'smc_coverage',
                'intervention_correlation', 'frac_high_access', 'Cohort_birth_month', 'rtss_mode')

plot_dir <- file.path(simout_dir, exp_name, '_plots')
if (!dir.exists(plot_dir))dir.create(plot_dir)

## Age month at dose
mth_adj <- function(x, agedose = 9) {
  mth <- rep(seq(1, 12), 5)
  xadj <- mth[x + agedose]
  return(xadj)
}


cases_df_list <- load_Age_monthly_Cases(simout_dir, exp_name, exp_sweeps, keep_birth_month = T)
dat <- cases_df_list[[4]]

dat <- dat %>%
  mutate(Cohort_RTSS_month = mth_adj(Cohort_birth_month, agedose = 18)) %>%
  as.data.table()

dat$age_group <- factor(dat$age_group,
                        levels = c("0 - 1", "1 - 2", "2 - 3", "3 - 4", "4 - 5", "5 - 6", "6 - 7", "7 - 8", "8 - 9", "9 - 10"),
                        labels = c("0 - 1", "1 - 2", "2 - 3", "3 - 4", "4 - 5", "5 - 6", "6 - 7", "7 - 8", "8 - 9", "9 - 10"))

seasonality_dir <- file.path(projectpath, "simulation_inputs/scenario_files/generic/Seasonality")
seasondat <- fread(file.path(seasonality_dir, 'seasonality_multipliers_high_unimodal.csv'))
seasondat$month <- c(1:12)

### RTSS overall
seasondat$month <- as.numeric(as.character(seasondat$month))
dat$Cohort_birth_month <- as.numeric(as.character(dat$Cohort_birth_month))
table(dat$age_group)
table(dat$rtss_coverage, dat$smc_coverage)
table(dat$Annual_EIR, dat$cm_coverage)
table(dat$rtss_coverage, dat$cm_coverage)

dat$rtss_mode_fct <- factor(dat$rtss_mode, levels = c('constant', 'campboost3', 'campboost3B'),
                            labels = c('EPI-RTS,S\n(initial dose at 9 months,\n booster at 24 months of age)',
                                       'campaign-RTS,S\n(initial dose at 9 months,\n booster at 24-46 months of age)',
                                       'campaign-RTS,S\n(initial dose at 9 months,\n booster at 24-46 months of age\n +booster at 47-59 months)'))

#outcome <- severe_cases_outcomes[1]
outcome <- cases_outcomes[1]

pplot <- ggplot(data = subset(dat, age_group %in% c('0 - 1', '1 - 2', '2 - 3', '3 - 4', '4 - 5', '5 - 6') &
  Annual_EIR == 10 &
  cm_coverage == 0.6 &
  smc_coverage %in% c(0) &
  rtss_coverage %in% c(0.8))) +
  geom_bar(data = seasondat, aes(x = month, y = monthly_EIR_multiplier / 6), alpha = 0.4, stat = "identity") +
  geom_line(aes(x = Cohort_birth_month + 1, y = get(outcome), group = interaction(rtss_coverage, age_group), col = age_group)) +
  scale_color_viridis_d(option = "B", begin = 0.2, end = 0.9, direction = -1) +
  labs(title = '',
       y = gsub("_", " ", outcome),
       x = 'Birth month of child',
       color = 'Age of child\n(years)') +
  customTheme +
  scale_x_continuous(labels = c(1:12), breaks = c(1:12)) +
  theme(legend.position = 'right', panel.spacing = unit(1, "lines")) +
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.2),
                     sec.axis = sec_axis(trans = ~. * 6, name = "Monthly EIR")) +
  facet_wrap(~rtss_mode_fct, scales = 'free_y')

f_save_plot(pplot, plot_name = paste0(outcome, '_birthmonth_plot'),
            width = 12, height = 4,
            plot_dir = file.path(simout_dir, exp_name, '_plots'))


# dat$season_born <- 'begin_dry'
# dat$season_born[dat$Cohort_birth_month %in% c(4:7)] <- 'mid_dry'
# dat$season_born[dat$Cohort_birth_month %in% c(8:11)] <- 'wet'
#
# pplot <- ggplot(data = subset(dat,
#   Annual_EIR == 10 &
#   cm_coverage==0.6&
#   smc_coverage %in% c(0) &
#   rtss_coverage %in% c(0.8))) +
#   # geom_bar(data = seasondat, aes(x = month, y = monthly_EIR_multiplier / 6), alpha = 0.4, stat = "identity") +
#   geom_smooth(aes(x = age_group, y = get(cases_outcomes[1]), group = interaction(season_born), col = as.factor(season_born)), se = FALSE) +
#   scale_color_viridis_d(option = "B", begin = 0.2, end = 0.9, direction = -1) +
#   labs(y = cases_outcomes[1], x = 'Age', col = 'Season when child was born') +
#   customTheme +
#   theme(legend.position = 'right', panel.spacing = unit(1, "lines")) +
#   scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.2))+
#   facet_wrap(~rtss_mode)
#
#
# f_save_plot(pplot, plot_name = 'custom_boxplot',
#             width = 10, height = 4,
#             plot_dir = file.path(simout_dir, exp_name, '_plots', Uage))



