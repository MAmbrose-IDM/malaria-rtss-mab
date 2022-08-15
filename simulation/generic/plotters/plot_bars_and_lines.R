library(data.table)
library(dplyr)
library(ggplot2)

#wdir at rtss-scenarios
source(file.path('simulation', 'generic', 'plotters', 'plotter_helpers.R'))


plot_barplots <- function(dat, xvar = 'scenario_name', yvar = 'cases_averted_per100000',
                          fillvar = 'scenario_name', facet1 = 'seasonality', facet2 = NULL,
                          SAVE = FALSE, ylab='Cases averted per 100000', scales='free',add_watermark = FALSE) {

  
  dat <- as.data.frame(dat)
  dat$facet1 = dat[,facet1]  #paste0(facet1, '\n', dat[, facet1])
  dat$facet2 = paste0(facet2, '\n', dat[, facet2])
  
  pplot <- ggplot(data = dat, aes(x =as.factor(get(xvar)), y = get(yvar))) +
    geom_bar(stat = "identity", aes(fill = as.factor(get(fillvar))),
             position = position_dodge(width = 1)) +
    scale_fill_brewer(palette = "YlGnBu",
                      guide = guide_legend(reverse = TRUE)) +
    xlab(element_blank())+
    ylab(ylab) +
    customTheme +
    theme(panel.spacing = unit(1, "lines"), 
          legend.position="none",
          axis.text.x = element_text(angle = 45, vjust = 0.5))
  
  if (!is.null(facet1) & !is.null(facet2)) pplot <- pplot + facet_grid(facet1 ~ facet2, scales = scales)
  if (!is.null(facet1) & is.null(facet2)) pplot <- pplot + facet_wrap(~facet1, scales = scales)
  
  if (SAVE) {
    dir.create(file.path(simout_dir, exp_name, '_plots'))
    f_save_plot(pplot, plot_name = paste('barplot', yvar, facet1, sep = '_'),
                width = 8, height = 6,
                plot_dir = file.path(simout_dir, exp_name, '_plots'))
  }
  return(pplot)
}




plot_grouped_barplots <- function(dat, xvar = 'Annual_EIR', yvar = 'cases_averted_per100000', bargroup_var='scenario_name',
                                  fillvar = 'scenario_name', facet1 = 'seasonality', facet2 = 'smc_coverage',
                                  SAVE = FALSE, ylab='Cases averted per 100000',add_watermark = FALSE, scales='free', bar_width=0.9) {
  
  dat <- as.data.frame(dat)
  dat$facet1 = paste0(facet1, '\n', dat[, facet1])
  dat$facet2 = paste0(facet2, '\n', dat[, facet2])
  
  pplot <- ggplot(data = dat, aes(x =as.factor(get(xvar)), y = get(yvar))) +
    geom_bar(stat = "identity", aes(fill = as.factor(get(fillvar)), group=interaction(as.factor(get(bargroup_var)), as.factor(facet1), as.factor(facet2))),
             position = position_dodge(width = 0.9), width=bar_width) +
    # position = position_dodge(width = 1)) +
    scale_fill_brewer(palette = "Spectral", # "YlGnBu",
                      guide = guide_legend(reverse = TRUE)) +
    xlab(element_blank())+
    ylab(ylab) +
    theme_bw() #+
    # theme(panel.spacing = unit(1, "lines"), 
    #       # legend.position="none",
    #       axis.text.x = element_text(angle = 45, vjust = 0.5))
  
  if (!is.null(facet1) & !is.null(facet2)) pplot <- pplot + facet_grid(facet1 ~ facet2, scales = scales)
  if (!is.null(facet1) & is.null(facet2)) pplot <- pplot + facet_wrap(~facet1, scales = scales)

  if (add_watermark)pplot <- watermark(pplot)

  if (SAVE) {
    dir.create(file.path(simout_dir, exp_name, '_plots'))
    f_save_plot(pplot, plot_name = paste('barplot', yvar, facet1, sep = '_'),
                width = 8, height = 6,
                plot_dir = file.path(simout_dir, exp_name, '_plots'))
  }
  return(pplot)
}



get_label_name = function(column_name){
  lookup_labels = c('rtss_cases_averted_per100000' = 'annual cases averted \n by RTS,S per 100,000', 
                    'rtss_severe_cases_averted_per100000' = 'annual severe cases averted \n by RTS,S per 100,000',
                    'clinical_cases_no_rtss' = 'annual incidence (per 1,000) \n without RTS,S',
                    'pfpr_2_10_no_rtss'='annual PfPR (2-10) \n without RTS,S',
                    'smc_coverage' = 'SMC coverage')
  label = as.character(lookup_labels[column_name])
  if(is.na(label)) label = column_name
  return(label)
}

get_max_xlim = function(dat, xvar, group_var, crop_to_second_largest){
  if(crop_to_second_largest){
    aa = aggregate(dat[[xvar]], by=list(dat[[group_var]]), max)
    max_vals = aa$x
    # get second largest
    return(max(max_vals[-which.max(max_vals)]))
  } else{
    return(max(dat[[xvar]]))
  }

}

plot_2x2grid_lineplots <- function(dat, xvars = c('clinical_cases_no_rtss', 'pfpr_2_10_no_rtss'), yvars = c('rtss_cases_averted_per100000', 'rtss_severe_cases_averted_per100000'), 
                                  colvar = 'seasonality', shapevar = NA, crop_to_second_largest=FALSE, add_watermark = FALSE) {
  
  dat <- as.data.frame(dat)
  
  gg1 = ggplot(dat, aes(x=get(xvars[1]), y=get(yvars[1]))) +
    geom_line(aes(col=as.factor(get(colvar))), size=1.12)+  
    geom_point(aes(col=as.factor(get(colvar))), size=1.2)+  
    xlab(get_label_name(xvars[1])) +
    ylab(get_label_name(yvars[1])) + 
    coord_cartesian(xlim=c(0,get_max_xlim(dat, xvar=xvars[1], group_var=colvar, crop_to_second_largest=crop_to_second_largest))) +
    theme_bw()+
    guides(col=guide_legend(get_label_name(colvar))) +
    theme(legend.position=c(0.8,0.23), 
          legend.background = element_rect(fill="transparent")) + 
    f_getCustomTheme() +
    f_getColorMapping(name=colvar, type='color')
  gg2 = ggplot(pe_df_cur, aes(x=get(xvars[2]), y=get(yvars[1]))) +
    geom_line(aes(col=as.factor(get(colvar))), size=1.12)+  
    geom_point(aes(col=as.factor(get(colvar))), size=1.2)+  
    xlab(get_label_name(xvars[2])) +
    ylab(get_label_name(yvars[1])) + 
    coord_cartesian(xlim=c(0,get_max_xlim(dat, xvar=xvars[2], group_var=colvar, crop_to_second_largest=crop_to_second_largest))) +
    theme_bw()+
    theme(legend.position='none') + 
    f_getCustomTheme() +
    f_getColorMapping(name=colvar, type='color')
  gg3 = ggplot(pe_df_cur, aes(x=get(xvars[1]), y=get(yvars[2]))) +
    geom_line(aes(col=as.factor(get(colvar))), size=1.12)+  
    geom_point(aes(col=as.factor(get(colvar))), size=1.2)+  
    xlab(get_label_name(xvars[1])) +
    ylab(get_label_name(yvars[2])) + 
    coord_cartesian(xlim=c(0,get_max_xlim(dat, xvar=xvars[1], group_var=colvar, crop_to_second_largest=crop_to_second_largest))) +
    theme_bw()+
    theme(legend.position='none') + 
    f_getCustomTheme() +
    f_getColorMapping(name=colvar, type='color')
  gg4 = ggplot(pe_df_cur, aes(x=get(xvars[2]), y=get(yvars[2]))) +
    geom_line(aes(col=as.factor(get(colvar))), size=1.12)+  
    geom_point(aes(col=as.factor(get(colvar))), size=1.2)+  
    xlab(get_label_name(xvars[2])) +
    ylab(get_label_name(yvars[2])) + 
    coord_cartesian(xlim=c(0,get_max_xlim(dat, xvar=xvars[2], group_var=colvar, crop_to_second_largest=crop_to_second_largest))) +
    theme_bw()+
    theme(legend.position='none') + 
    f_getCustomTheme() +
    f_getColorMapping(name=colvar, type='color')

    # combine plots in grid
  tt=plot_grid(gg1,gg2,gg3,gg4, nrow=2, align = "hv")
  
  if (add_watermark)tt <- watermark(tt)

  return(tt)
}


plot_lineplots_with_labels <- function(dat, xvar = 'clinical_cases_no_rtss', yvar = 'rtss_severe_cases_averted_per100000', 
                                   colvar = 'seasonality', colorMapping=NA, shapevar = NA, add_watermark = FALSE, ymax=800, xmax=3000) {  # crop_to_second_largest=FALSE, 
  
  dat <- as.data.frame(dat)
  if(is.factor(dat[[xvar]])) dat[[xvar]] = as.numeric(as.character(dat[[xvar]]))

  gg = ggplot(dat, aes(x=get(xvar), y=get(yvar))) +
    geom_line(aes(group=get(colvar), col=as.factor(get(colvar))), size=1.12)+  
    # geom_point(aes(col=as.factor(get(colvar))), size=1.2)+  
    xlab(get_label_name(xvar)) +
    ylab(get_label_name(yvar)) + 
    coord_cartesian(xlim=c(0, xmax), ylim=c(0,ymax)) +  # get_max_xlim(dat, xvar=xvar, group_var=colvar, crop_to_second_largest=crop_to_second_largest))) +
    theme_bw()+
    guides(col=guide_legend(get_label_name(colvar))) +
    theme_classic() +
    theme(legend.position=c(0.3,0.8), 
          legend.background = element_rect(fill="transparent")) + 
    colorMapping
  
  
  if (add_watermark)gg <- watermark(gg)
  
  return(gg)
}



plot_lines <- function(dat, xvar = 'Annual_EIR', yvar = 'rtss_cases_averted_per100000',
                       colvar = 'age_group', facet1 = 'smc_coverage', facet2 = NULL,
                       SAVE = FALSE, xlab='annual EIR', ylab='Cases averted per 100000',
                       add_watermark=FALSE, scales='free') {
  
  dat <- as.data.frame(dat)
  dat$facet1 = paste0(facet1, '\n', dat[, facet1])
  dat$facet2 = paste0(facet1, '\n', dat[, facet2])
  
  pplot <- ggplot(data = dat, aes(x=get(xvar), y = get(yvar), group = get(colvar))) +
    geom_line(aes(col = as.factor(get(colvar))), size=1) +
    scale_color_brewer(palette = "GnBu") +
    xlab(xlab)+
    ylab(ylab) +
    # ylim(0,400)+
    customTheme 
  # theme(panel.spacing = unit(1, "lines"))
  
  if (!is.null(facet1) & !is.null(facet2)) pplot <- pplot + facet_grid(facet1 ~ facet2, scales = scales)
  if (!is.null(facet1) & is.null(facet2)) pplot <- pplot + facet_wrap(~facet1, scales = scales)

  if (add_watermark)pplot <- watermark(pplot)

  if (SAVE) {
    dir.create(file.path(simout_dir, exp_name, '_plots'))
    f_save_plot(pplot, plot_name = paste('barplot', yvar, facet1, sep = '_'),
                width = 8, height = 6,
                plot_dir = file.path(simout_dir, exp_name, '_plots'))
  }
  return(pplot)
}


plot_burden_by_age <- function(cases_scenarios_each_year, cur_corr, cur_cm, cur_smc, add_watermark=FALSE){
  # subset simulation to current input
  cur_df = filter(cases_scenarios_each_year,
                  intervention_correlation == cur_corr,
                  cm_coverage == cur_cm,
                  smc_coverage == cur_smc)
  cur_df = cur_df[!is.na(cur_df$scenario_name),]
  
  # create plot
  age_label_values = seq(0, 8, length.out=5)
  age_labels = paste0(age_label_values, '-', (age_label_values + 1))
  gg = ggplot(cur_df, aes(x=year, y=cases_averted_per100000))+
    geom_point(aes(col=as.factor(scenario_name), group=scenario_name), size=2)+
    geom_line(aes(col=as.factor(scenario_name), group=scenario_name), size=1)+
    f_getColorMapping(name='scenario', type='color') +
    theme_bw()+
    f_getCustomTheme() +
    geom_hline(yintercept=0)+
    ylab('cases averted (per 100000) when adding RTS,S and/or SMC') + 
    xlab('age of child') +
    scale_x_continuous(breaks=age_label_values, labels=age_labels) +
    facet_wrap(seasonality~Annual_EIR)

  if (add_watermark)gg <- watermark(gg)
  return(gg)
}


