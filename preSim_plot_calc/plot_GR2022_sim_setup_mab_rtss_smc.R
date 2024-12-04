# plot_GR2022_sim_setup_mab_rtss_smc.R
# plot EIRs and seasonality patterns, along with timings for SMC/RTSS/mAb distributions

library(reshape2)
library(ggplot2)
library(ggalt)

# project filepath
projectpath = 'C:/Users/moniqueam/Dropbox (IDM)/Malaria Team Folder/projects/mAb_rtss_smc_comparison'

# set which seasonalities to plot
seasonality_types = c('constant', 'moderate_unimodal', 'high_unimodal', 'higher_unimodal')
# set which EIRs to plot
annual_eir = c(5, 10, 30, 60, 120)
# set dates for intervention distributions
# mab_month = 175/30.4
# rtss_month = 175/30.4
# smc_months = c(175, 205, 235, 265)/30.4
mab_month = 190/30.4
rtss_month = 190/30.4
smc_months = (15+c(175, 205, 235, 265))/30.4

# seasonality profiles
season_month_scale = read.csv(paste0(projectpath, '/simulation_inputs/seasonality/seasonality_eir_multipliers.csv'))
season_month_scale = season_month_scale[, which(colnames(season_month_scale) %in% c('month', seasonality_types))]
# add another point mid-way through the month with the same value
season_month_scale1 = season_month_scale
season_month_scale1$month = season_month_scale1$month + 0.35
season_month_scale = rbind(season_month_scale, season_month_scale1)

# combine seasonality profiles and EIRs into long-format data frame
season_month_scale2 = melt(season_month_scale, id.vars=c('month'))
colnames(season_month_scale2)[colnames(season_month_scale2) == 'variable'] = 'seasonality'
colnames(season_month_scale2)[colnames(season_month_scale2) == 'value'] = 'month_scalar'
monthly_eirs = merge(season_month_scale2, data.frame(annual_eir), all=TRUE)
monthly_eirs$monthly_eir = monthly_eirs$month_scalar * monthly_eirs$annual_eir
monthly_eirs = monthly_eirs[order(monthly_eirs$seasonality, monthly_eirs$annual_eir, monthly_eirs$month),]
monthly_eirs$seasonality = factor(monthly_eirs$seasonality, levels=seasonality_types)

# create dataframe of intervention dates
intervention_dates = data.frame(month=c(mab_month, rtss_month, smc_months), intervention=c('mAb', 'RTS,S', rep('SMC', length(smc_months))))
intervention_dates = merge(intervention_dates, data.frame('seasonality'=seasonality_types))
intervention_dates$seasonality = factor(intervention_dates$seasonality, levels=seasonality_types)
max_month_eir = max(monthly_eirs$monthly_eir)

monthly_eirs$daily_eir = monthly_eirs$monthly_eir/30.4

# create plot, faceted by seasonality
gg = ggplot() +
  geom_xspline(data=monthly_eirs, aes(x=month, y=monthly_eir, col=annual_eir, group=annual_eir), spline_shape=0.3, size=1)+
  geom_ribbon(data=monthly_eirs, aes(x=month, ymax=monthly_eir, fill=annual_eir, group=annual_eir),ymin=0,alpha=0.3) +
  # geom_point(data=monthly_eirs, aes(x=month, y=monthly_eir, col=annual_eir, group=annual_eir), size=2)+
  geom_point(data=intervention_dates[intervention_dates$intervention=='SMC',], aes(x=month, y=(max_month_eir+5)), size=2, shape=25) +  #
  # geom_point(data=intervention_dates, aes(x=month, fill=factor(intervention), y=(max_month_eir+20 - 3*as.numeric(factor(intervention)))), size=2, shape=25) +  #
  # scale_fill_manual(breaks=c('mAb', 'RTS,S', 'SMC'), values=c(rgb(149/255,140/255,180/255), rgb(36/255,162/255,140/255), rgb(209/255,171/255,133/255))) +
  # scale_color_manual(breaks=c('mAb', 'RTS,S', 'SMC'), values=c(rgb(149/255,140/255,180/255), rgb(36/255,162/255,140/255), rgb(209/255,171/255,133/255))) +
  scale_x_continuous(breaks=c(1, 4, 7, 10), labels=c('Jan', 'Apr', 'Jul', 'Oct')) +
  facet_wrap('seasonality', nrow=1) +
  theme_classic()

ggsave(plot=gg, filename=paste0(projectpath, '/nonSimFigures/GR_illustration_EIR_intervention_timing.png'),
       width = 8, height = 3, units = 'in')


# create plot, faceted by seasonality - with DAILY EIR instead of monthly EIR
gg2 = ggplot() +
  geom_xspline(data=monthly_eirs, aes(x=month, y=daily_eir, col=annual_eir, group=annual_eir), spline_shape=0.3, size=1)+
  geom_ribbon(data=monthly_eirs, aes(x=month, ymax=daily_eir, fill=annual_eir, group=annual_eir),ymin=0,alpha=0.3) +
  geom_xspline(data=monthly_eirs, aes(x=month, y=daily_eir, col=annual_eir, group=annual_eir), spline_shape=0.3, size=1)+
  # geom_point(data=monthly_eirs, aes(x=month, y=monthly_eir, col=annual_eir, group=annual_eir), size=2)+
  geom_point(data=intervention_dates[intervention_dates$intervention=='SMC',], aes(x=month, y=(max_month_eir/30.4+.1)), size=2, shape=25) +  #
  # geom_point(data=intervention_dates, aes(x=month, fill=factor(intervention), y=(max_month_eir+20 - 3*as.numeric(factor(intervention)))), size=2, shape=25) +  #
  # scale_fill_manual(breaks=c('mAb', 'RTS,S', 'SMC'), values=c(rgb(149/255,140/255,180/255), rgb(36/255,162/255,140/255), rgb(209/255,171/255,133/255))) +
  # scale_color_manual(breaks=c('mAb', 'RTS,S', 'SMC'), values=c(rgb(149/255,140/255,180/255), rgb(36/255,162/255,140/255), rgb(209/255,171/255,133/255))) +
  scale_x_continuous(breaks=c(1, 4, 7, 10), labels=c('Jan', 'Apr', 'Jul', 'Oct')) +
  facet_wrap('seasonality', nrow=1) +
  theme_classic()

ggsave(plot=gg2, filename=paste0(projectpath, '/nonSimFigures/GR_illustration_EIR_intervention_timing.png'),
       width = 8, height = 3, units = 'in')






# create plot of daily EIR for highest seasonality, compared to mAb efficacy-through-time - note, need to run plot_pkpd_mAb_RTSS.R first to generate.
monthly_eirs2 = monthly_eirs
monthly_eirs2$month = monthly_eirs2$month + 12
monthly_eirs3 = monthly_eirs2
monthly_eirs3$month = monthly_eirs2$month + 12
monthly_eirs2 = rbind(monthly_eirs, monthly_eirs2, monthly_eirs3)
monthly_eirs2$day = monthly_eirs2$month * 30.4
monthly_eirs2 = monthly_eirs2[monthly_eirs2$annual_eir == 60,]
monthly_eirs2 = monthly_eirs2[monthly_eirs2$seasonality == 'higher_unimodal',]
daily_eirs = monthly_eirs2
daily_eirs$day = daily_eirs$day - mab_month * 30.4
daily_eirs = daily_eirs[daily_eirs$day>=0,]
daily_eirs = daily_eirs[daily_eirs$day<=(365*2),]
# gg3 = ggplot() +
#   geom_xspline(data=monthly_eirs2, aes(x=day, y=daily_eir, group=annual_eir), col='dodgerblue', spline_shape=0.3, size=1)+
#   geom_ribbon(data=monthly_eirs2, aes(x=day, ymax=daily_eir, group=annual_eir), fill='dodgerblue', ymin=0,alpha=0.3) +
#   geom_point(data=intervention_dates[intervention_dates$intervention=='mAb',], aes(x=month*30.4, y=(max_month_eir/30.4+.1)), size=2, shape=25) +  #
#   # scale_x_continuous(breaks=c(1, 4, 7, 10), labels=c('Jan', 'Apr', 'Jul', 'Oct')) +
#   theme_classic()
gg3 = ggplot() +
  geom_xspline(data=daily_eirs, aes(x=day, y=daily_eir, group=annual_eir), col='dodgerblue', spline_shape=0.3, size=1)+
  geom_ribbon(data=daily_eirs, aes(x=day, ymax=daily_eir, group=annual_eir), fill='dodgerblue', ymin=0,alpha=0.3) +
  # geom_point(data=intervention_dates[intervention_dates$intervention=='mAb',], aes(x=month*30.4, y=(max_month_eir/30.4+.1)), size=2, shape=25) +  #
  # scale_x_continuous(breaks=c(1, 4, 7, 10), labels=c('Jan', 'Apr', 'Jul', 'Oct')) +
  theme_classic()
ggsave(plot=gg3, filename=paste0(projectpath, '/nonSimFigures/GR_illustration_seasonality_after_intervention.png'),
       width = 4, height = 3, units = 'in')


