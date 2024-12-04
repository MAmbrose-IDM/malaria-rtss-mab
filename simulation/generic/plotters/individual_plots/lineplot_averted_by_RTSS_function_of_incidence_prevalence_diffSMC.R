# lineplot_averted_by_RTSS_function_of_incidence_prevalence_diffSMC.R


library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)


#######################################################
# setup
#######################################################
# wdir at rtss-scenarios
USER <- Sys.getenv('USERNAME')
setwd(paste0('C:/Users/', USER, '/Documents/rtss-scenarios'))

source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))
source(file.path('simulation', 'generic', 'plotters', 'plot_bars_and_lines.R'))
source(file.path("simulation", "load_paths.R"))

simout_dir <- file.path(projectpath, 'simulation_output', 'generic')
if (!dir.exists(paste0(simout_dir, '/_plots'))) dir.create(paste0(simout_dir, '/_plots'))
if (!dir.exists(paste0(simout_dir, '/_plots/pdf'))) dir.create(paste0(simout_dir, '/_plots/pdf'))

exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'cm_coverage', 'rtss_coverage', 'rtss_mode', 'smc_coverage',
                'intervention_correlation', 'Cohort_birth_month', 'minBoostAge')
max_years = c(1, 2, 5, 10)
age_groups = paste0('U', max_years)
seasonality_levels = c('high', 'moderate', 'constant')
rtss_mode_levels = c('constant', 'campboost', 'campboost2')
rtss_mode_levels_plot = c('standard', 'enhanced (1b)', 'enhanced (2b)')
rtss_mode_levels_lookup = rtss_mode_levels_plot
names(rtss_mode_levels_lookup) = rtss_mode_levels

theme_set(theme_bw())
customTheme <- f_getCustomTheme()
device_format <- c('pdf', 'png')





####################################################################################################
# panel of three line plots of incidence (x-axis) and cases averted when RTS,S added (y-axis) for:
  # 1) different effective CMs
  # 2) different SMC coverages
  # 3) different seasonalities
####################################################################################################

# plot by incidence, prevalence, or EIR on x-axis?
x_axis_background_metric = 'prevalence'  # eir, incidence, prevalence
if(x_axis_background_metric=='incidence'){
  x_var ='clinical_cases_no_rtss'
  xmax=2800
} else if(x_axis_background_metric=='prevalence'){
  x_var='pfpr_2_10_no_rtss'
  xmax=0.5
} else if(x_axis_background_metric=='eir'){
  x_var='Annual_EIR'
  xmax=40
}

# what burdenreduction metric should be plotted on the y-axis? 
burden_metric='severeAverted'  # severeAverted, clinicalAverted, percentReductSevere, percentReductClinical
if(burden_metric=='severeAverted'){
  y_var='rtss_severe_cases_averted_per100000'
  ymax=700
} else if(burden_metric=='clinicalAverted'){
  y_var='rtss_cases_averted_per100000'
  ymax=50000
} else if(burden_metric=='percentReductSevere'){
  y_var='rtss_percent_reduction_severe'
  ymax=30
} else if(burden_metric=='percentReductClinical'){
  y_var='rtss_percent_reduction'
  ymax=30
}

# default characteristics of plotted scenarios
cur_age = 'U5'
cur_rtss_scenario = '80% standard RTS,S (1 booster)'
cur_cm = 0.6
cur_corr = 0
cur_season = 'high_unimodal'
cur_smc=0
include_eir_labels=TRUE
eir_labels = c(1, 10, 30)


#   =====   seasonality   =====   #
exp_name = 'generic_eir_season_sweep2'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))
pe_df$seasonality = gsub('_unimodal', '', pe_df$seasonality)
pe_df$seasonality = factor(pe_df$seasonality, levels = seasonality_levels)

if(grepl('percent', burden_metric)){
  # percent reduction
  pe_df$rtss_percent_reduction = pe_df$rtss_relative_burden * -100
  pe_df$rtss_percent_reduction_severe = pe_df$rtss_relative_burden_severe * -100
}

# subset simulations
cur_seasonalities = c('moderate', 'high')
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   rtss_coverage > 0.001,
                   rtss_scenario == cur_rtss_scenario,
                   smc_coverage == cur_smc,
                   cm_coverage == cur_cm,
                   intervention_correlation == cur_corr,
                   seasonality %in% cur_seasonalities
)

# specify desired colors for each level
levels = c('high', 'moderate', 'constant')
myColors = c(rgb(0, 0, 0),
             rgb(0.1, 0.3, 0.6),
             rgb(0.5, 0.7, 1)
)
names(myColors) = levels
colorMapping = scale_color_manual(name = 'seasonality', values = myColors)


# plot
gg_season = plot_lineplots_with_labels(dat = pe_df_cur, xvar = x_var, yvar = y_var, 
                                       colvar = 'seasonality', colorMapping=colorMapping, shapevar = NA, add_watermark = FALSE,  # , crop_to_second_largest=FALSE
                                       ymax=ymax, xmax=xmax)

if(include_eir_labels){
  pe_df_cur_subset = filter(pe_df_cur,
                            Annual_EIR %in% eir_labels)
  gg_season = gg_season + 
    geom_point(data=pe_df_cur_subset, aes(x=get(x_var), y=get(y_var)), color=rgb(0.5,0.5,0.5, 0.5), size=2)
  
  # if(x_axis_background_metric != 'eir'){
  #   pe_df_cur_subset2 = filter(pe_df_cur_subset,
  #                              cm_coverage == cur_cm)
  #   pe_df_cur_subset2$eir_label = paste0(pe_df_cur_subset2$Annual_EIR, ' ibpa')
  #   gg_cm = gg_cm + geom_text(data=pe_df_cur_subset2, aes(x=get(x_var), y=get(y_var), label=eir_label), color=rgb(0.5,0.5,0.5), vjust=-0.5, hjust=1, size=3)
  # }
}

# # seasonality - clinical cases only (create and save a separate figure)
# gg_season_clinical = plot_lineplots_with_labels(dat = pe_df_cur, xvar = 'clinical_cases_no_rtss', yvar = 'rtss_cases_averted_per100000', 
#                                        colvar = 'seasonality', colorMapping=colorMapping, shapevar = NA, add_watermark = FALSE,  # , crop_to_second_largest=FALSE
#                                        ymax=40000, xmax=xmax)
# f_save_plot(gg_season_clinical, paste0('clinical_averted_RTSS_by_incidence_seasonality_',
#                        cur_corr, 'corr'),
#             file.path(simout_dir, '_plots'), width = 3, height = 3.5, units = 'in', device_format = device_format)


#   =====   CM   =====   #
exp_name = 'generic_rtss_eir_cm_sweep'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]
if(grepl('percent', burden_metric)){
  # percent reduction
  pe_df$rtss_percent_reduction = pe_df$rtss_relative_burden * -100
  pe_df$rtss_percent_reduction_severe = pe_df$rtss_relative_burden_severe * -100
}

# subset to plotted scenario
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   rtss_coverage > 0.001,
                   rtss_scenario == cur_rtss_scenario,
                   smc_coverage == cur_smc,
                   intervention_correlation == cur_corr,
                   seasonality == cur_season
)

# specify desired colors for each level
levels = c(0.3, 0.6, 0.9)
myColors = c(rgb(0.6, 0.1, 0.3),
             rgb(0, 0, 0),
             rgb(1, 0.5, 0.7))
names(myColors) = levels

colorMapping = scale_color_manual(name = 'cm_coverage', values = myColors)
# plot
gg_cm = plot_lineplots_with_labels(dat = pe_df_cur, xvar = x_var, yvar = y_var, 
                                    colvar = 'cm_coverage', colorMapping=colorMapping, shapevar = NA, add_watermark = FALSE,  # , crop_to_second_largest=FALSE
                                   ymax=ymax, xmax=xmax)

# add labels for EIR
if(include_eir_labels){
  pe_df_cur_subset = filter(pe_df_cur,
                            Annual_EIR %in% eir_labels)
  gg_cm = gg_cm + 
    geom_point(data=pe_df_cur_subset, aes(x=get(x_var), y=get(y_var)), color=rgb(0.5,0.5,0.5, 0.5), size=2)
    
    # if(x_axis_background_metric != 'eir'){
    #   pe_df_cur_subset2 = filter(pe_df_cur_subset,
    #                              cm_coverage == cur_cm)
    #   pe_df_cur_subset2$eir_label = paste0(pe_df_cur_subset2$Annual_EIR, ' ibpa')
    #   gg_cm = gg_cm + geom_text(data=pe_df_cur_subset2, aes(x=get(x_var), y=get(y_var), label=eir_label), color=rgb(0.5,0.5,0.5), vjust=-0.5, hjust=1, size=3)
    # }
}

#   =====   SMC   =====   #
exp_name = 'generic_eir_smc_sweep3'  # 'generic_eir_smc_sweep2_newRDT'  # 'generic_eir_smc_sweep'  # 'generic_eir_smc_sweep2_newRDT'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]
if(grepl('percent', burden_metric)){
  # percent reduction
  pe_df$rtss_percent_reduction = pe_df$rtss_relative_burden * -100
  pe_df$rtss_percent_reduction_severe = pe_df$rtss_relative_burden_severe * -100
}

# subset to plotted scenario
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   rtss_coverage > 0.001,
                   rtss_scenario == cur_rtss_scenario,
                   cm_coverage == cur_cm,
                   intervention_correlation == cur_corr,
                   seasonality == cur_season
)

# specify desired colors for each level
levels = c(0, 0.5, 0.8)
myColors = c(rgb(0, 0, 0),
             rgb(0.1, 0.6, 0.3),
             rgb(0.5, 1, 0.7))
names(myColors) = levels
colorMapping = scale_color_manual(name = 'smc_coverage', values = myColors)
# plot
gg_smc = plot_lineplots_with_labels(dat = pe_df_cur, xvar = x_var, yvar = y_var, 
                                    colvar = 'smc_coverage', colorMapping=colorMapping, shapevar = NA, add_watermark = FALSE,  # , crop_to_second_largest=FALSE
                                    ymax=ymax, xmax=xmax)
# add labels for EIR
if(include_eir_labels){
  pe_df_cur_subset = filter(pe_df_cur,
                          Annual_EIR %in% eir_labels)
  gg_smc = gg_smc + 
    geom_point(data=pe_df_cur_subset, aes(x=get(x_var), y=get(y_var)), color=rgb(0.5,0.5,0.5, 0.5), size=2)
    
  if(x_axis_background_metric != 'eir'){
    pe_df_cur_subset2 = filter(pe_df_cur_subset,
                              smc_coverage < 0.001)
    pe_df_cur_subset2$eir_label = paste0(pe_df_cur_subset2$Annual_EIR, ' ibpa')
    gg_smc = gg_smc + geom_text(data=pe_df_cur_subset2, aes(x=get(x_var), y=get(y_var), label=eir_label), color=rgb(0.5,0.5,0.5), vjust=-0.5, hjust=1, size=3)
  }
}

# combine plots in grid
tt=plot_grid(gg_season, gg_cm, gg_smc, nrow=1, align = "h")
f_save_plot(tt, paste0(burden_metric, '_RTSS_by_',x_axis_background_metric,'_CM_SMC_seasonality_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 9, height = 3.5, units = 'in', device_format = device_format)












####################################################################################################
# line plot of total cases in setting and cases averted when RTS,S added for different scenarios
####################################################################################################
cur_age = 'U5'
cur_rtss_scenario = '80% standard RTS,S (1 booster)'
cur_cm = 0.6
cur_corr = 0


# === across EIRs and SMCs (a line for each SMC coverage) === #
exp_name = 'generic_eir_smc_sweep2_newRDT'  # 'generic_eir_smc_sweep'  # 'generic_eir_smc_sweep2_newRDT'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = exp_sweeps,
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]

# subset to plotted scenario
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   rtss_coverage > 0.001,
                   rtss_scenario == cur_rtss_scenario,
                   cm_coverage == cur_cm,
                   intervention_correlation == cur_corr
)

tt = plot_2x2grid_lineplots(dat = pe_df_cur, xvars = c('clinical_cases_no_rtss', 'pfpr_2_10_no_rtss'),
                            yvars = c('rtss_cases_averted_per100000', 'rtss_severe_cases_averted_per100000'),
                            colvar = 'smc_coverage', crop_to_second_largest = TRUE)

f_save_plot(tt, paste0('case_and_severe_averted_RTSS_SMC_EIR_',
                       cur_cm * 100, 'CM_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 9, height = 6, units = 'in', device_format = device_format)



# === a line for each type of seasonality ===#
exp_name = 'generic_eir_season_sweep'
exp_filepath = paste0(simout_dir, '/', exp_name)

sim_output = load_Age_monthly_Cases(simout_dir = simout_dir, exp_name = exp_name, exp_sweeps = c(exp_sweeps, 'minBoostAge'),
                                    add_PE_perAge = TRUE, max_years = max_years, keep_birth_month = FALSE)
pe_df = sim_output[[3]]
pe_df = f_add_scenario_name(df = pe_df, scenario_type = 'rtss_scenario')
pe_df = pe_df[!is.na(pe_df$rtss_scenario),]
pe_df$Annual_EIR = factor(pe_df$Annual_EIR, levels = sort(unique(pe_df$Annual_EIR)))
pe_df$seasonality = gsub('_unimodal', '', pe_df$seasonality)
pe_df$seasonality = factor(pe_df$seasonality, levels = seasonality_levels)

# subset simulations
cur_smc = 0
pe_df_cur = filter(pe_df,
                   age_group == cur_age,
                   rtss_coverage > 0.001,
                   rtss_scenario == cur_rtss_scenario,
                   smc_coverage == cur_smc,
                   cm_coverage == cur_cm,
                   intervention_correlation == cur_corr
)

tt = plot_2x2grid_lineplots(dat = pe_df_cur, xvars = c('clinical_cases_no_rtss', 'pfpr_2_10_no_rtss'),
                            yvars = c('rtss_cases_averted_per100000', 'rtss_severe_cases_averted_per100000'),
                            colvar = 'seasonality')

f_save_plot(tt, paste0('case_and_severe_averted_RTSS_seasonality_EIR_',
                       cur_cm * 100, 'CM_',
                       cur_smc * 100, 'SMC_',
                       cur_corr, 'corr'),
            file.path(simout_dir, '_plots'), width = 9, height = 6, units = 'in', device_format = device_format)


# percent reduction patterns
pe_df_cur$rtss_percent_reduction = pe_df_cur$rtss_relative_burden * -100
pe_df_cur$rtss_percent_reduction_severe = pe_df_cur$rtss_relative_burden_severe * -100
tt = plot_2x2grid_lineplots(dat = pe_df_cur, xvars = c('clinical_cases_no_rtss', 'pfpr_2_10_no_rtss'),
                            yvars = c('rtss_percent_reduction', 'rtss_percent_reduction_severe'),
                            colvar = 'seasonality')



# gg1 = ggplot(pe_df_cur, aes(x=clinical_cases_no_rtss, y=rtss_cases_averted_per100000)) +
#   geom_point(aes(col=Annual_EIR, shape=as.factor(smc_coverage))) +
#   geom_line(aes(col=Annual_EIR, group=as.factor(smc_coverage)))+
#   ylab(paste0('annual cases ', cur_age, ' averted \n by RTS,S per 100,000')) + 
#   xlab('annual incidence without RTS,S (per 1,000)') +
#   ylim(0, NA)+
#   theme_bw()+
#   theme(legend.position='none') + 
#   f_getCustomTheme() 
# gg2 = ggplot(pe_df_cur, aes(x=clinical_cases_no_rtss, y=rtss_severe_cases_averted_per100000)) +
#   geom_point(aes(col=Annual_EIR, shape=as.factor(smc_coverage))) +
#   geom_line(aes(col=Annual_EIR, group=as.factor(smc_coverage)))+
#   ylab(paste0('annual severe cases ', cur_age, ' averted \n by RTS,S per 100,000')) + 
#   xlab('annual incidence without RTS,S (per 1,000)') +
#   ylim(0, NA)+
#   theme_bw()+
#   theme(legend.position='none') + 
#   f_getCustomTheme() 
# gg3 = ggplot(pe_df_cur, aes(x=pfpr_2_10_no_rtss, y=rtss_cases_averted_per100000)) +
#   geom_point(aes(col=Annual_EIR, shape=as.factor(smc_coverage))) +  
#   geom_line(aes(col=Annual_EIR, group=as.factor(smc_coverage))) +  
#   ylab(paste0('annual cases ', cur_age, ' averted \n by RTS,S per 100,000')) + 
#   xlab(expression(paste(italic(Pf), PR[2-10]~without~RTS, ',', S))) +
#   ylim(0, NA)+
#   theme_bw()+
#   theme(legend.position='none') + 
#   f_getCustomTheme() 
# gg4 = ggplot(pe_df_cur, aes(x=pfpr_2_10_no_rtss, y=rtss_severe_cases_averted_per100000)) +
#   geom_point(aes(col=Annual_EIR, shape=as.factor(smc_coverage))) +  
#   geom_line(aes(col=Annual_EIR, group=as.factor(smc_coverage))) +  
#   ylab(paste0('annual severe cases ', cur_age, ' averted \n by RTS,S per 100,000')) + 
#   xlab(expression(paste(italic(Pf), PR[2-10]~without~RTS, ',', S))) +
#   ylim(0, NA)+
#   theme_bw()+
#   theme(legend.position='none') + 
#   f_getCustomTheme() 
# grid.arrange(gg1,gg3,gg2,gg4, nrow=2)
