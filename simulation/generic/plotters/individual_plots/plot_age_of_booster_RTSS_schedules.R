# plot_age_of_booster_RTSS_schedules.R

# plot the age when children in each cohort first receives the RTS,S campaign booster

plot_days = FALSE
rtss_booster1_min_age = 2 * 365
start_day0 = round(365 + 6 * 30.4 - 7) # one week before SMC
cohort_month_shift = 0:11
start_days = start_day0 - round(30.4 * cohort_month_shift)
start_days_adjusted = c()
for (ss in 1:length(start_days)) {
  start_day = start_days[ss]
  while (start_day < rtss_booster1_min_age) {
    start_day = start_day + 365
  }
  start_days_adjusted[ss] = start_day
}
if (plot_days){
  # plot in days
  png(filename = paste0(simout_dir, '/_plots/age_at_first_booster.png'), width = 4, height = 4, units = 'in', res = 900)
  plot((cohort_month_shift + 1), start_days_adjusted, ylim = c(365 * 2 - 30, 365 * 3), xlab = 'birth month', ylab = c('age (in days) when a child', 'receives their first RTS,S booster'), bty = 'L', pch = 20, cex = 1.8)
  abline(v = (start_day0 %% 365) / 30.4 + 1, col = 'lightblue', lwd = 4)
  points((cohort_month_shift + 1), rep(365 * 2, length(cohort_month_shift)), pch = 20, cex = 1.8, col = 'grey')
  legend('topleft', c('standard', 'enhanced'), pch = 20, pt.cex = 1.8, col = c('grey', 'black'), bty = 'n')
  dev.off()
  
  
  gg <- ggplot() +
    geom_point(aes(x = (cohort_month_shift + 1), y = start_days_adjusted), size = 2.3) +
    geom_point(aes(x = (cohort_month_shift + 1), y = rep(365 * 2, length(cohort_month_shift))), size = 2.3, col = 'grey') +
    scale_y_continuous(lim = c(365 * 2 - 30, 365 * 3)) +
    scale_x_continuous(breaks = cohort_month_shift + 1) +
    geom_vline(xintercept = (start_day0 %% 365) / 30.4 + 1, col = 'lightblue', lwd = 2) +
    labs(x = 'birth month', y = 'age (in days) when a child\nreceives their first RTS,S booster') +
    theme_cowplot() +
    customTheme
  
  f_save_plot(gg, paste0('age_at_first_booster_v2'), file.path(simout_dir, '_plots'), width = 5, height = 4, units = 'in', device_format = device_format)
}else{
  # plot in months
  png(filename = paste0(simout_dir, '/_plots/age_at_first_booster.png'), width = 4, height = 4, units = 'in', res = 900)
  plot((cohort_month_shift + 1), start_days_adjusted/30.4, ylim = c(365 * 2 - 30, 365 * 3)/30.4, xlab = 'birth month', ylab = c('age (in months) when a child', 'receives their first RTS,S booster'), bty = 'L', pch = 20, cex = 1.8)
  abline(v = (start_day0 %% 365) / 30.4 + 1, col = 'lightblue', lwd = 4)
  points((cohort_month_shift + 1), rep(365 * 2, length(cohort_month_shift))/30.4, pch = 20, cex = 1.8, col = 'grey')
  legend('topleft', c('standard', 'enhanced'), pch = 20, pt.cex = 1.8, col = c('grey', 'black'), bty = 'n')
  dev.off()
  pdf(file = paste0(simout_dir, '/_plots/pdf/age_at_first_booster.pdf'), width = 4, height = 4, useDingbats=FALSE)
  plot((cohort_month_shift + 1), start_days_adjusted/30.4, ylim = c(365 * 2 - 30, 365 * 3)/30.4, xlab = 'birth month', ylab = c('age (in months) when a child', 'receives their first RTS,S booster'), bty = 'L', pch = 20, cex = 1.8)
  abline(v = (start_day0 %% 365) / 30.4 + 1, col = 'lightblue', lwd = 4)
  points((cohort_month_shift + 1), rep(365 * 2, length(cohort_month_shift))/30.4, pch = 20, cex = 1.8, col = 'grey')
  legend('topleft', c('standard', 'enhanced'), pch = 20, pt.cex = 1.8, col = c('grey', 'black'), bty = 'n')
  dev.off()
  
  gg <- ggplot() +
    geom_point(aes(x = (cohort_month_shift + 1), y = start_days_adjusted/30.4), size = 2.3) +
    geom_point(aes(x = (cohort_month_shift + 1), y = rep(365 * 2, length(cohort_month_shift))/30.4), size = 2.3, col = 'grey') +
    scale_y_continuous(lim = c(365 * 2 - 30, 365 * 3)/30.4) +
    scale_x_continuous(breaks = cohort_month_shift + 1) +
    geom_vline(xintercept = (start_day0 %% 365) / 30.4 + 1, col = 'lightblue', lwd = 2) +
    labs(x = 'birth month', y = 'age (in months) when a child\nreceives their first RTS,S booster') +
    theme_cowplot() +
    customTheme
  
  f_save_plot(gg, paste0('age_at_first_booster_v2'), file.path(simout_dir, '_plots'), width = 5, height = 4, units = 'in', device_format = device_format)
}

