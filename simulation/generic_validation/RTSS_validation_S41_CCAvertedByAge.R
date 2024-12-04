# RTSS_validation_S41_CCAvertedByAge.R
# August 2021
# create plots of RTS,S simluation output to compare against Figure S4.1 from Penny et al. 

library(lubridate)
library(dplyr)
library(ggplot2)


source(file.path("simulation", "load_paths.R"))
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
device_format = c('pdf', 'png')


exp_filepath = file.path(projectpath, 'simulation_output/generic/rtss_validation_SweepEIR_noBooster')
cases_filepath = file.path(exp_filepath, 'All_Age_monthly_Cases.csv')

cases_df = read.csv(cases_filepath)

# get total within each year
cases_df$date = as.Date(cases_df$date)
cases_df$year = lubridate::year(cases_df$date)
cases_df$year = cases_df$year - min(cases_df$year)

cases_df_annual = cases_df %>%
  dplyr::group_by(Scenario_id, Run_Number, year, Annual.EIR, rtss_coverage) %>%
  dplyr::summarise(total_cases = sum(New.Clinical.Cases),
                   total_severe_cases = sum(New.Severe.Cases),
                   average_pfpr = mean(PfHRP2.Prevalence),
                   average_pop = mean(Statistical.Population),
  ) %>%
  dplyr::ungroup()

# get averages across runs
cases_df_annual_ave = cases_df_annual %>%
  dplyr::group_by(Scenario_id, year, Annual.EIR, rtss_coverage) %>%
  dplyr::summarise(clinical_cases = mean(total_cases),
                   severe_cases = mean(total_severe_cases),
                   pfpr = mean(average_pfpr),
                   pop = mean(average_pop)
  ) %>%
  dplyr::ungroup()


# calculate the number of clinical cases averted for each with-RTSS / without-RTSS pair
# merge pairs of simulations with versus without RTSS
with_rtss = cases_df_annual_ave[cases_df_annual_ave$rtss_coverage > 0, which(colnames(cases_df_annual_ave) %in% c('year', 'Annual.EIR', 'clinical_cases'))]
colnames(with_rtss)[which(colnames(with_rtss) == 'clinical_cases')] = 'clinical_cases_rtss'
without_rtss = cases_df_annual_ave[cases_df_annual_ave$rtss_coverage < 0.00001, which(colnames(cases_df_annual_ave) %in% c('year', 'Annual.EIR', 'clinical_cases'))]
colnames(without_rtss)[which(colnames(without_rtss) == 'clinical_cases')] = 'clinical_cases_no_rtss'

comparison_df = merge(with_rtss, without_rtss, by = c('year', 'Annual.EIR'))
comparison_df$cases_averted = comparison_df$clinical_cases_no_rtss - comparison_df$clinical_cases_rtss

# subset to EIR corresponding to Figure S4.1
comparison_df_subset = comparison_df[comparison_df$Annual.EIR %in% c(0.23, 0.91, 3.61, 9.3, 16.54),]
comparison_df_subset$cases_averted_per_100000 = comparison_df_subset$cases_averted * 100
gg = ggplot(data = comparison_df_subset, aes(x = year, y = cases_averted_per_100000)) +
  geom_bar(stat = 'identity', fill = rgb(0.1, 0.55, 0.4)) +
  theme_light() +
  ylab('clinical cases averted per 100000 fully vaccinated (no booster)') +
  xlab('age') +
  facet_wrap(facets = vars(Annual.EIR), ncol = 1, scales = 'free')
f_save_plot(gg, paste0('clinical_cases_averted_barplot'),
            file.path(exp_filepath), width = 2.25, height = 7.2, units = 'in', device_format = device_format)


plot_heights = c(5000, 15000, 40000, 100000)
plot_mins = c(-500, -1000, -4000, -10000)
eir_vals = c(0.23, 0.91, 3.61, 9.3)
for (i_eir in 1:length(eir_vals)) {
  comparison_df_subset_eir = comparison_df_subset[comparison_df_subset$Annual.EIR == eir_vals[i_eir],]
  gg = ggplot(data = comparison_df_subset_eir, aes(x = year, y = cases_averted)) +
    geom_bar(stat = 'identity', fill = rgb(0.1, 0.55, 0.4)) +
    theme_light() +
    ylab('clinical cases averted per 1000 fully vaccinated (no booster)') +
    xlab('age') +
    ylim(plot_mins[i_eir] / 100, plot_heights[i_eir] / 100)
  f_save_plot(gg, paste0('clinical_cases_averted_barplot_eir', round(eir_vals[i_eir] * 100)),
              file.path(exp_filepath), width = 2.1, height = 1.5, units = 'in', device_format = device_format)

}


# compare EIR-PfPR(2-10) relationship
cases_df_annual_2_10 = cases_df_annual[intersect(which(cases_df_annual$year > 1), which(cases_df_annual$year < 10)),]
pfpr_2_10_average = cases_df_annual_2_10 %>%
  dplyr::group_by(Annual.EIR) %>%
  dplyr::summarise(average_pfpr = mean(average_pfpr)) %>%
  dplyr::ungroup()

plot(cases_df_annual_2_10$average_pfpr[cases_df_annual_2_10$rtss_coverage < 0.001], cases_df_annual_2_10$Annual.EIR[cases_df_annual_2_10$rtss_coverage < 0.001], col = 'brown', pch = 20, cex = 0.5,
     xlab = 'PfPR (2-10)', ylab = 'annual EIR', bty = 'L')
points(pfpr_2_10_average$average_pfpr, pfpr_2_10_average$Annual.EIR, pch = 20, col = 'brown', cex = 2)
# add points for the old relationship (shared on Slack by Caitlin on August 26 2021)
previous_eir_pfpr_relationship = data.frame('Prev_2_10' = c(3, 10, 20, 30, 40, 50, 60, 70, 75) / 100,
                                            'Annual.EIR' = c(0.23, 0.91, 2.06, 3.61, 5.75, 9.3, 16.54, 26.72, 36.68))
points(previous_eir_pfpr_relationship$Prev_2_10, previous_eir_pfpr_relationship$Annual.EIR, pch = 20, col = 'blue', cex = 2)
legend('topleft', legend = c('old values', 'new values'), pch = 20, col = c('blue', 'brown'), bty = 'n')

pfpr_2_10_average$Prev_2_10 = pfpr_2_10_average$average_pfpr * 100
print(pfpr_2_10_average[, c('Prev_2_10', 'Annual.EIR')])
