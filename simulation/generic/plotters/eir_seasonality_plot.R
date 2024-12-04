require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)

source(file.path("simulation", "load_paths.R"))
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path('simulation', 'generic', 'plotters', 'plot_bars_and_lines.R'))

simout_dir <- file.path(projectpath, 'simulation_output', 'generic_forward')
theme_set(theme_bw())
customTheme <- f_getCustomTheme()

dat_eir <- fread(file.path(projectpath, 'simulation_inputs/scenario_files/generic/Seasonality', 'seasonality_eir_multipliers.csv')) %>%
  pivot_longer(cols = -c(month))

eir_values <- as.data.frame(cbind('multiplier' = c(1, 5, 10, 30)))
dat_eir <- dat_eir %>%
  left_join(eir_values, by = character()) %>%
  mutate(value_scl = value * multiplier)

selected_seasonalities <- c("constant","moderate unimodal", "high unimodal")
dat_eir$name <- factor(dat_eir$name, levels = c("constant","moderate_unimodal", "high_unimodal"),
                       labels = c("constant","moderate unimodal", "high unimodal"))

pplot <- ggplot(data = subset(dat_eir, name %in% selected_seasonalities)) +
  geom_line(aes(x = month, y = value_scl, col = as.factor(multiplier), group = multiplier)) +
  scale_y_continuous(breaks = c(1:6), labels = c(1:6), lim = c(0, 7), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1:12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~name) +
  labs(x = "", y = "Monthly EIR", col = 'Annual EIR') +
  customTheme

f_save_plot(pplot, plot_name = 'seasonality_inputs',
            width = 14, height = 4,
            plot_dir = file.path(projectpath,'simulation_inputs/scenario_files/generic/Seasonality'))