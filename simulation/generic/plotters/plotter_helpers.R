library(dplyr)
library(data.table)
library('ggVennDiagram')

#require(devtools)
#source_gist("9d50c86e21161188a6c17c7c9650ad5d", filename = "f_save_plot.R")
#source_gist("9d50c86e21161188a6c17c7c9650ad5d", filename = "f_getCustomTheme.R")

watermark <- function(p, label = "PRELIMINARY", col = "darkgrey", cex = 4.5, hjust = 1, vjust = -0.5) {
  p <- p + annotate("text", x = Inf, y = -Inf, label = label,
                    hjust = hjust, vjust = vjust, angle = -30, col = col, cex = cex,
                    fontface = "bold", alpha = 0.8)
  return(p)
}

f_save_plot <- function(pplot, plot_name, plot_dir, width = 14, height = 8, units = 'in', device_format = c('pdf','png')) {
  if ('png' %in% device_format) {
    ggsave(paste0(plot_name, ".png"), plot = pplot, path = plot_dir,
           width = width, height = height, units = units, device = "png")
  }
  if ('pdf' %in% device_format) {
    if (!dir.exists(file.path(plot_dir, "pdf"))) { dir.create(file.path(plot_dir, "pdf")) }
    ggsave(paste0(plot_name, ".pdf"), plot = pplot, path = file.path(plot_dir, "pdf"),
           width = width, height = height, units = units, device = "pdf", useDingbats = FALSE)
  }
}

f_load_sim_PEestimates <- function(simout_dir, exp_name, exp_sweeps) {
  fname <- file.path(simout_dir, exp_name, 'protective_efficacy_estimates.csv')
  if (file.exists(fname)) {
    dat <- fread(fname)
  }else {
    cases_df_list <- load_Age_monthly_Cases(simout_dir, exp_name, exp_sweeps)
    dat <- cases_df_list[[3]]
    fwrite(dat, file = fname)
  }
  return(dat)
}

f_getCustomTheme <- function(fontscl = 1) {

  customTheme <- theme(
    strip.text.x = element_text(size = 12 * fontscl, face = "bold"),
    strip.text.y = element_text(size = 12 * fontscl, face = "bold"),
    strip.background = element_blank(),
    plot.title = element_text(size = 14 * fontscl, vjust = -1, hjust = 0),
    plot.subtitle = element_text(size = 12 * fontscl),
    plot.caption = element_text(size = 9 * fontscl),
    legend.title = element_text(size = 12 * fontscl),
    legend.text = element_text(size = 12 * fontscl),
    axis.title.x = element_text(size = 12 * fontscl),
    axis.text.x = element_text(size = 12 * fontscl),
    axis.title.y = element_text(size = 12 * fontscl),
    axis.text.y = element_text(size = 12 * fontscl),
    axis.ticks = element_line(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

  return(customTheme)
}

# specify colors for plotted interventions (based on hh for mAbs and specified manually otherwise)
get_intervention_colors = function(df, level_name = 'smc_rtss_mab', rtss_col=rgb(0.5,0.5,0.5), smc_col='#924900', none_col='black', hh_vals=c(2, 5, 10, 20, 40, 60, 80)){
  # Note: in ggplot, use:   scale_colour_manual(name='smc_rtss_mab', breaks=intervention_names, values=intervention_colors) +
  if(is.factor(df[[level_name]])){
    intervention_names = levels(df[[level_name]])
  } else{
    intervention_names = sort(unique((df[[level_name]])))
  }
  intervention_colors = rep("#b6dbff", length(intervention_names))
  hh_cols = viridis(length(hh_vals))
  # pie(rep(1, length(hh_cols)), labels=hh_vals, col=hh_cols)
  for (hh in 1:length(hh_vals)){
    if(level_name == 'hh'){
      intervention_colors[intervention_names == hh_vals[hh]] = hh_cols[hh]
    } else{
      intervention_colors[grep(paste0('hh=',hh_vals[hh]), intervention_names)] = hh_cols[hh]
    }
  }
  intervention_colors[grep('RTS,S', intervention_names)] = rtss_col
  intervention_colors[grep('SMC', intervention_names)] = smc_col
  intervention_colors[intervention_names == 'none'] = none_col
  # pie(rep(1, length(intervention_colors)), labels=intervention_names, col=intervention_colors)
  return(list(intervention_names, intervention_colors))
}
# example colors
# pal <- c("#004949","#009292",
#          +          "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
#          +          "#924900")
# pie(rep(1,length(pal)), col=pal)



f_getColors <- function() {
  myColors = c(rgb(0.6, 0.6, 1), rgb(0.3, 0.3, 0.9),
               rgb(0.65, 0.95, 0.65), rgb(0.4, 0.9, 0.4),
               rgb(0.2, 0.92, 0.85), rgb(0, 0.8, 0.6),
               rgb(0, 0.6, 0.8), rgb(0, 0.5, 0.5))
  return(myColors)
}

f_getColorMapping <- function(name = 'scenario', type = 'fill') {

  if (name == 'scenario') {
    levels = c("60% RTS,S; 0% SMC", "80% RTS,S; 0% SMC",
               "0% RTS,S; 60% SMC", "0% RTS,S; 80% SMC",
               "60% RTS,S; 60% SMC", "60% RTS,S; 80% SMC",
               "80% RTS,S; 60% SMC", "80% RTS,S; 80% SMC",
               "RTS,S only", "SMC only", "both RTS,S and SMC",
               "without SMC", "with SMC",
               'standard RTS,S - without SMC', 'standard RTS,S - with SMC', 'enhanced RTS,S - without SMC', 'enhanced RTS,S - with SMC', 
               
               "no RTS,S; no SMC", 
               "80% RTS,S; no SMC",
               "no RTS,S; 80% SMC",
               "80% RTS,S; 80% SMC; independent",
               "80% RTS,S; 80% SMC; correlated",
               "80% RTS,S; 80% SMC; correlatedSMC",
               "80% RTS,S; 80% SMC; correlatedRTSS",
               "80% RTS,S; 80% SMC; anti-correlated",
               
               "60% RTS,S; no SMC",
               "no RTS,S; 60% SMC",
               "60% RTS,S; 60% SMC; independent",
               "60% RTS,S; 60% SMC; correlated",
               "60% RTS,S; 60% SMC; correlatedSMC",
               "60% RTS,S; 60% SMC; correlatedRTSS",
               "60% RTS,S; 60% SMC; anti-correlated",
               
               "no SMC",
               "no RTS,S",
               "correlated",
               "independent",
               "anti-correlated")
    myColors = c(rgb(0.6, 0.6, 1), rgb(0.3, 0.3, 0.9),
                 rgb(0.65, 0.95, 0.65), rgb(0.4, 0.9, 0.4),
                 rgb(0.2, 0.92, 0.85), rgb(0, 0.8, 0.6),
                 rgb(0, 0.6, 0.8), rgb(0, 0.5, 0.5),
                 rgb(0.3, 0.3, 0.9), rgb(0.4, 0.9, 0.4), rgb(0, 0.5, 0.5),
                 rgb(0.25, 0.25, 0.8), rgb(0.1, 0.7, 0.6),
                 rgb(0.25, 0.25, 0.8), rgb(0.3, 0.8, 0.95), rgb(0.6, 0.0, 0.6), rgb(0.9, 0.4, 0.95),
                 
                 rgb(0.5,0.5,0.5),
                 rgb(210/255, 170/255, 187/255),
                 rgb(100/255, 100/255, 100/255),
                 rgb(65/255, 196/255, 196/255),
                 rgb(110/255, 170/255, 240/255),
                 rgb(164/255, 219/255, 133/255),
                 rgb(164/255, 219/255, 133/255),
                 rgb(164/255, 219/255, 133/255),
                 
                 rgb(210/255, 170/255, 187/255),
                 rgb(100/255, 100/255, 100/255),
                 rgb(65/255, 196/255, 196/255),
                 rgb(110/255, 170/255, 240/255),
                 rgb(164/255, 219/255, 133/255),
                 rgb(164/255, 219/255, 133/255),
                 rgb(164/255, 219/255, 133/255),
                 
                 rgb(210/255, 170/255, 187/255),
                 rgb(100/255, 100/255, 100/255),
                 rgb(65/255, 196/255, 196/255),
                 rgb(110/255, 170/255, 240/255),
                 rgb(164/255, 219/255, 133/255)

                 # rgb(255/255, 170/255, 187/255),
                 # rgb(153/255, 153/255, 51/255),
                 # rgb(68/255, 187/255, 153/255),
                 # rgb(153/255, 221/255, 255/255),
                 # rgb(119/255, 170/255, 221/255),
                 # rgb(119/255, 170/255, 221/255),
                 # rgb(100/255, 155/255, 205/255)
                ) 
  } else if (name == 'age group') {
    levels = c('U1', 'U2', 'U5', 'U10')
    myColors = c(rgb(0 / 255, 21 / 255, 156 / 255),
                 rgb(120 / 255, 0 / 255, 170 / 255),
                 rgb(200 / 255, 67 / 255, 115 / 255),
                 rgb(250 / 255, 125 / 255, 125 / 255))
  } else if (name == 'RTS,S coverage') {
    levels = c(0, 0.6, 0.8)
    myColors = c(rgb(0.65, 0.95, 0.65),
                 rgb(0.2, 0.92, 0.85),
                 rgb(0, 0.6, 0.8))
  } else if (name == 'seasonality') {
    levels = c('constant', 'moderate_unimodal', 'high_unimodal',
               'moderate', 'high')
    myColors = c(rgb(207 / 255, 237 / 255, 139 / 255),
                 rgb(127 / 255, 205 / 255, 187 / 255),
                 rgb(44 / 255, 127 / 255, 184 / 255),
                 rgb(127 / 255, 205 / 255, 187 / 255),
                 rgb(44 / 255, 127 / 255, 184 / 255))
  } else if (name == 'smc_coverage') {
    levels = c(0, 0.6, 0.8)
    myColors = c(rgb(0.1, 0.6, 0.3),
                 rgb(0.3, 0.8, 0.5),
                 rgb(0.5, 1, 0.7))
  } else if (name == 'RTS,S scenario') {
    levels = c('no RTS,S', '60% standard RTS,S (1 booster)', '80% standard RTS,S (1 booster)',
               '80% enhanced RTS,S (1 booster)', '80% enhanced RTS,S (2 boosters)',
               '80% enhanced RTS,S (1 booster, min 24m)', '80% enhanced RTS,S (2 boosters, min 24m)',
               '80% enhanced RTS,S (1 booster, min 18m)', '80% enhanced RTS,S (2 boosters, min 18m)')
    myColors = c(rgb(0.65, 0.95, 0.65), rgb(0.2, 0.92, 0.85), rgb(0, 0.6, 0.8),
                 rgb(0.7, 0.5, 0.9), rgb(0.4, 0.1, 0.7),
                 rgb(0.7, 0.5, 0.9), rgb(0.4, 0.1, 0.7),
                 rgb(1, 0.5, 0.9), rgb(0.8, 0.1, 0.7))
  }
  # plot(1:length(myColors), col=myColors, pch=20, cex=3)
  names(myColors) = levels
  if (type == 'fill') {
    return(scale_fill_manual(name = name, values = myColors))
  } else { # color
    return(colScale = scale_color_manual(name = name, values = myColors))
  }
}


f_add_scenario_name <- function(df, scenario_type) {
  if (scenario_type == 'vacc_info') {
    # add vaccine characteristics from vacc_char tag
    df$vacc_type = str_extract(df$vacc_char, pattern='[a-z]*')
    df$vacc_type[df$vacc_type=='mab'] = 'mAb'
    df$vacc_type[df$vacc_type=='rtss'] = 'RTS,S'
    df$initial_conc = as.numeric(str_replace(str_extract(df$vacc_char, pattern='ic[0-9]*'), pattern='ic', replacement=''))
    # df$max_efficacy = as.numeric(str_replace(str_extract(df$vacc_char, pattern='ime[0-9]*\\.*[0-9]*'), pattern='ime', replacement=''))
    df$max_efficacy = as.numeric(str_replace(str_extract(df$vacc_char, pattern='ime[0-9]*'), pattern='ime', replacement=''))
    df$hh = str_replace(str_extract(df$vacc_char, pattern='ih[0-9]*'), pattern='ih', replacement='')
    df$hh = as.numeric(df$hh)
    df$vacc_info = paste0(df$vacc_type, ': ', df$max_efficacy, '% max, hh=', df$hh, ', initC=', df$initial_conc)
    df$vacc_info[df$vacc_coverage < 0.001] = 'none'
    # get the desired levels, sorted appropriately
    unique_combos = distinct(df[,c('vacc_type', 'max_efficacy', 'hh', 'initial_conc')])
    unique_combos$vacc_info = paste0(unique_combos$vacc_type, ': ', unique_combos$max_efficacy, '% max, hh=', unique_combos$hh, ', initC=', unique_combos$initial_conc)
    unique_combos = unique_combos[order(unique_combos$max_efficacy),]
    unique_combos = unique_combos[order(unique_combos$hh,decreasing=TRUE),]
    unique_combos = unique_combos[order(unique_combos$vacc_type),]
    df$vacc_info = factor(df$vacc_info, levels = c('none', unique_combos$vacc_info))
  } else if(scenario_type == 'smc_rtss_mab'){
    df = f_add_scenario_name(df, 'vacc_info')
    df$smc_rtss_mab = as.character(df$vacc_info)
    df$smc_rtss_mab[df$vacc_coverage <0.01] = ''
    df$smc_rtss_mab[df$smc_coverage >0] = paste0('SMC. ', df$smc_rtss_mab[df$smc_coverage >0])
    df$smc_rtss_mab[df$smc_rtss_mab==''] = 'none'
    df$smc_rtss_mab = factor(df$smc_rtss_mab, levels=c('SMC. ', levels(df$vacc_info), paste0('SMC. ', levels(df$vacc_info))))
  } else if (scenario_type == 'smc_rtss_mab_age'){
    df = f_add_scenario_name(df, 'smc_rtss_mab')
    df$smc_rtss_mab_age = as.character(df$smc_rtss_mab)
    df$smc_rtss_mab_age[df$smc_rtss_mab_age !='none'] = paste0('U', df$max_target_age[df$smc_rtss_mab_age !='none'], '. ', df$smc_rtss_mab_age[df$smc_rtss_mab_age !='none'])
    all_combos = expand.grid(levels(df$smc_rtss_mab), sort(unique(df$max_target_age)))
    df$smc_rtss_mab_age = factor(df$smc_rtss_mab_age, levels=c('none', paste0('U', all_combos$Var2, '. ', all_combos$Var1)))
    
  }else {
    warning(paste0('scenario_type: ', scenario_type, ' not recognized by the f_add_scenario_name function.'))
  }
  return(df)
}

add_counterfactual <- function(dat, exp_name_counterfactual, fname = 'All_Age_monthly_Cases.csv') {

  ## RTSS counterfactual
  if (min(dat$rtss_coverage) > 0) {
    cases_df_counterfactual <- fread(file.path(simout_dir, exp_name_counterfactual, fname)) %>%
      rename_with(~gsub(" ", "_", .x)) %>%
      mutate(Scenario_id = paste0(Scenario_id, '_counterfactual_rtss')) %>%
      filter(rtss_coverage == 0 &
               smc_coverage %in% unique(dat$smc_coverage) &
               cm_coverage %in% unique(dat$cm_coverage) &
               seasonality %in% unique(dat$seasonality) &
               Annual_EIR %in% unique(dat$Annual_EIR))

    if ('rtss_mode' %in% colnames(dat)) {
      cases_df_counterfactual <- cases_df_counterfactual %>% dplyr::select(-rtss_mode)
      ## repeat rtss 0 coverage for all rtss_modes
      rtss_mode = unique(dat$rtss_mode)
      cases_df_counterfactual <- cases_df_counterfactual %>%
        left_join(as.data.frame(rtss_mode), by = character())
    }
    if ('minBoostAge' %in% colnames(dat)) {
      cases_df_counterfactual <- cases_df_counterfactual %>% dplyr::select(-minBoostAge)
      ## repeat rtss 0 coverage for all minBoostAge
      minBoostAge = unique(dat$minBoostAge)
      cases_df_counterfactual <- cases_df_counterfactual %>%
        left_join(as.data.frame(minBoostAge), by = character())
    }

    cases_df_counterfactual <- cases_df_counterfactual %>% dplyr::select(colnames(dat))
    if (sum(grep('/', cases_df_counterfactual$date[1])) > 0) {
      cases_df_counterfactual$date <- as.Date(cases_df_counterfactual$date, format = '%d/%m/%Y')
    }
    dat$date <- as.Date(dat$date)
    dat <- rbind(dat, cases_df_counterfactual)
  }

  ## SMC counterfactual
  if (min(dat$smc_coverage) > 0) {
    cases_df_counterfactual <- fread(file.path(simout_dir, exp_name_counterfactual, fname)) %>%
      rename_with(~gsub(" ", "_", .x)) %>%
      mutate(Scenario_id = paste0(Scenario_id, '_counterfactual_smc')) %>%
      filter(smc_coverage == 0 &
               rtss_coverage %in% unique(dat$rtss_coverage) &
               cm_coverage %in% unique(dat$cm_coverage) &
               seasonality %in% unique(dat$seasonality) &
               Annual_EIR %in% unique(dat$Annual_EIR)) %>%
      dplyr::select(colnames(dat))
    if (sum(grep('/', cases_df_counterfactual$date[1])) > 0) {
      cases_df_counterfactual$date <- as.Date(cases_df_counterfactual$date, format = '%d/%m/%Y')
    }
    dat$date <- as.Date(dat$date)
    dat <- rbind(dat, cases_df_counterfactual)
  }
  return(dat)
}
# 
# load_summary_data <- function(simout_dir, exp_name, exp_sweeps = Null) {
# 
#   if (is.null(exp_sweeps)) {
#     exp_sweeps <- c('Scenario_id', 'Annual_EIR', 'seasonality', 'Cohort_birth_month', 'vacc_char',
#                     'vacc_coverage', 'cm_coverage', 'smc_coverage', 'frac_high_access',
#                     'cm_target_group', 'smc_target_group', 'rtss_target_group')
#   }
#   grpVars1 <- c('Run_Number', exp_sweeps)
# 
#   U5dat <- fread(file.path(simout_dir, exp_name, 'U5_PfPR_ClinicalIncidence.csv')) %>%
#     rename_with(~gsub(" ", "_", .x)) %>%
#     dplyr::group_by_at(grpVars1) %>%
#     dplyr::summarize(PfPR_U5 = mean(PfPR_U5),
#                      Cases_U5 = sum(Cases_U5),
#                      Pop_U5 = mean(Pop_U5)) %>%
#     dplyr::group_by_at(exp_sweeps) %>%
#     dplyr::summarize(PfPR_U5 = mean(PfPR_U5),
#                      Cases_U5 = mean(Cases_U5),
#                      Pop_U5 = mean(Pop_U5)) %>%
#     dplyr::mutate(pos_U5 = PfPR_U5 * Pop_U5,
#                   total_cases_U5 = Cases_U5 * Pop_U5 * 30 / 365) %>%
#     dplyr::mutate(incidence = total_cases_U5 / Pop_U5 * 1000)
# 
#   U1dat <- fread(file.path(simout_dir, exp_name, 'U1_PfPR_ClinicalIncidence.csv')) %>%
#     rename_with(~gsub(" ", "_", .x)) %>%
#     dplyr::group_by_at(grpVars1) %>%
#     dplyr::summarize(PfPR_U1 = mean(PfPR_U1),
#                      Cases_U1 = sum(Cases_U1),
#                      Pop_U1 = mean(Pop_U1)) %>%
#     dplyr::group_by_at(exp_sweeps) %>%
#     dplyr::summarize(PfPR_U1 = mean(PfPR_U1),
#                      Cases_U1 = mean(Cases_U1),
#                      Pop_U1 = mean(Pop_U1)) %>%
#     dplyr::mutate(pos_U1 = PfPR_U1 * Pop_U1,
#                   total_cases_U1 = Cases_U1 * Pop_U1 * 30 / 365) %>%
#     dplyr::mutate(incidence = total_cases_U1 / Pop_U1 * 1000)
# 
#   simdat_list <- list(U1dat, U5dat)
#   return(simdat_list)
# }


load_Age_monthly_Cases <- function(simout_dir, exp_name, exp_sweeps = NULL, add_PE_perAge = TRUE, max_years = c(1, 2, 5, 8),
                                   exp_name_counterfactual = NULL, keep_birth_month = FALSE, fname = 'All_Age_monthly_Cases.csv') {

  if (length(exp_name)>1){
    for(ii in 1:length(exp_name)){
      cases_df_cur <- fread(file.path(simout_dir, exp_name[ii], fname)) %>% rename_with(~gsub(" ", "_", .x))
      if(ii==1){
        cases_df = cases_df_cur
      } else{
        cases_df = merge(cases_df, cases_df_cur, all=TRUE)
      }
    }
  } else{
    cases_df <- fread(file.path(simout_dir, exp_name, fname)) %>% rename_with(~gsub(" ", "_", .x))
  }

  # if not simulated, add counterfactual scenarios, matched by exp_sweeps except RTSS coverage (and SMC coverage)
  if (!is.null(exp_name_counterfactual))cases_df <- add_counterfactual(dat = cases_df, exp_name_counterfactual)


  # When using All_Age_Annual_Cases as returned from IPTi postprocessing, year had already been generated
  if (fname != 'All_Age_Annual_Cases.csv') {
    cases_df <- cases_df %>%
      mutate(date = as.Date(date),
             year = lubridate::year(date),
             year = year - min(year, na.rm=TRUE))
  }

  exp_sweeps <- c(exp_sweeps, c('Scenario_id', 'Annual_EIR', 'seasonality', 'Cohort_birth_month', 'max_target_age', 'vacc_char',
                                'vacc_coverage', 'cm_coverage', 'smc_coverage', 'frac_high_access',
                                'cm_target_group', 'smc_target_group', 'vacc_target_group'))
  if (any(!(exp_sweeps %in% colnames(cases_df)))) exp_sweeps <- exp_sweeps[-which(!(exp_sweeps %in% colnames(cases_df)))]

  # aggregate to the year level
  cases_df_annual = cases_df %>%
    dplyr::group_by_at(c(exp_sweeps, 'Run_Number', 'year')) %>%
    dplyr::summarise(total_cases = sum(New_Clinical_Cases),
                     total_severe_cases = sum(New_Severe_Cases),
                     average_pfpr = mean(PfHRP2_Prevalence),
                     average_pop = mean(Statistical_Population)) %>%
    dplyr::ungroup()

  # get averages across runs and optional across birth cohorts; calculate rates per 1000 people per year
  if (!keep_birth_month & ('Cohort_birth_month' %in% exp_sweeps)) exp_sweeps <- exp_sweeps[-which(exp_sweeps == 'Cohort_birth_month')]
  if (('Scenario_id' %in% exp_sweeps)) exp_sweeps <- exp_sweeps[-which(exp_sweeps == 'Scenario_id')]
  cases_df_annual_ave = cases_df_annual %>%
    dplyr::group_by_at(c(exp_sweeps, 'year')) %>%
    dplyr::summarise(clinical_cases = mean(total_cases),
                     severe_cases = mean(total_severe_cases),
                     pfpr = mean(average_pfpr),
                     pop = mean(average_pop)
    ) %>%
    dplyr::ungroup() %>%
    mutate(clinical_cases = clinical_cases / pop * 1000,
           severe_cases = severe_cases / pop * 1000)

  if (add_PE_perAge) {
    cases_against_reference <- get_PE_and_averted_estimates(dat = cases_df_annual_ave, exp_sweeps = exp_sweeps, max_years = max_years)
    cases_age_aggregated = cases_against_reference[[1]]
    cases_each_year = cases_against_reference[[2]]
    simdat <- list(cases_df_annual, cases_df_annual_ave, cases_age_aggregated, cases_each_year)
  }else {
    simdat <- list(cases_df_annual, cases_df_annual_ave)
  }
  return(simdat)

}


add_PfPR_2_10 = function(prev_df_ave, pfpr_colname, exp_sweeps) {
  # cut out ages under 2 and over 10
  df = prev_df_ave[prev_df_ave$year >= 2 & prev_df_ave$year < 10,]
  # remove the year column from df (since PfPR_2_10 is for the full setting, not a specific year)
  df = df[,-which(colnames(df)=='year')]
  
  # aggregate across ages within each scenario
  if (any(!(exp_sweeps %in% colnames(df)))) exp_sweeps <- exp_sweeps[-which(!(exp_sweeps %in% colnames(df)))]
  df_ave = df %>%
    dplyr::group_by_at(c(exp_sweeps)) %>%
    dplyr::summarise(PfPR_2_10 = mean(get(pfpr_colname))) %>%
    dplyr::ungroup()

  # merge the 2-10 PfPR values into the main dataframe
  merged_df = merge(prev_df_ave, df_ave, by = exp_sweeps, all=TRUE)
  return(merged_df)
}


get_PE_and_averted_estimates <- function(dat, exp_sweeps, max_years = c(1, 2, 5, 8), reference = 'cm_coverage') {
  # calculate PE and cases averted relative to matched CM-only scenario and to the matched no-RTS,S scenario

  # Comparison 1: reference scenarios are those with one intervention only (i.e., CM-only)
  intervention_coverages <- c('cm_coverage', 'vacc_coverage', 'smc_coverage', 'ipti_coverage')
  if (any(!(intervention_coverages %in% colnames(dat)))) intervention_coverages <- intervention_coverages[-which(!(intervention_coverages %in% colnames(dat)))]
  intervention_coverages <- intervention_coverages[-which(intervention_coverages == reference)]

  # remove non-reference scenarios and associated columns
  remove_cols <- c(intervention_coverages, 'Scenario_id', 'vacc_char', 'vacc_target_group', 'smc_target_group', 'frac_high_access', 'max_target_age')
  if (any(!remove_cols %in% colnames(dat))) remove_cols <- remove_cols[-which(!(remove_cols %in% colnames(dat)))]
  df_references = dat %>%
    filter(across(intervention_coverages, ~. == 0)) %>%
    dplyr::select(c(all_of(exp_sweeps), 'year', 'clinical_cases', 'severe_cases', 'pfpr')) %>%
    dplyr::select(-all_of(remove_cols)) %>%
    rename(clinical_cases_ref = clinical_cases,
           severe_cases_ref = severe_cases)
  # add PfPR_2_10 and get mean across any duplicate simulations
  ref_exp_sweeps = c(exp_sweeps, 'year')
  if (any(!(ref_exp_sweeps %in% colnames(df_references)))) ref_exp_sweeps <- ref_exp_sweeps[-which(!(ref_exp_sweeps %in% colnames(df_references)))]
  df_references = add_PfPR_2_10(prev_df_ave = df_references, pfpr_colname = 'pfpr', exp_sweeps = ref_exp_sweeps) %>%
    rename(pfpr_2_10_ref = PfPR_2_10) %>%
    dplyr::select(-pfpr) %>%
    dplyr::group_by_at(c(ref_exp_sweeps)) %>%
    dplyr::summarise(clinical_cases_ref = mean(clinical_cases_ref),
                     severe_cases_ref = mean(severe_cases_ref),
                     pfpr_2_10_ref = mean(pfpr_2_10_ref))
  # remove any duplicate rows that remain
  df_references = dplyr::distinct(df_references)

  # Comparison 2: reference scenarios are those without RTS,S
  # remove non-reference scenarios and associated columns
  remove_cols <- c('vacc_coverage', 'Scenario_id', 'vacc_target_group', 'vacc_char')
  if (any(!remove_cols %in% colnames(dat))) remove_cols <- remove_cols[-which(!(remove_cols %in% colnames(dat)))]
  df_no_vacc_references = dat %>%
    filter(vacc_coverage == 0) %>%
    dplyr::select(c(all_of(exp_sweeps), 'year', 'clinical_cases', 'severe_cases', 'pfpr')) %>%
    dplyr::select(-all_of(remove_cols)) %>%
    rename(clinical_cases_no_vacc = clinical_cases,
           severe_cases_no_vacc = severe_cases)
  # add PfPR_2_10 and get mean across any duplicate simulations
  ref_exp_sweeps = c(exp_sweeps, 'year')
  if (any(!(ref_exp_sweeps %in% colnames(df_no_vacc_references)))) ref_exp_sweeps <- ref_exp_sweeps[-which(!(ref_exp_sweeps %in% colnames(df_no_vacc_references)))]
  df_no_vacc_references = add_PfPR_2_10(prev_df_ave = df_no_vacc_references, pfpr_colname = 'pfpr', exp_sweeps = ref_exp_sweeps) %>%
    rename(pfpr_2_10_no_vacc = PfPR_2_10) %>%
    dplyr::select(-pfpr) %>%
    dplyr::group_by_at(c(ref_exp_sweeps)) %>%
    dplyr::summarise(clinical_cases_no_vacc = mean(clinical_cases_no_vacc),
                     severe_cases_no_vacc = mean(severe_cases_no_vacc),
                     pfpr_2_10_no_vacc = mean(pfpr_2_10_no_vacc))
  # remove any duplicate rows that remain
  df_no_vacc_references = dplyr::distinct(df_no_vacc_references)
  
  # merge the reference values into the data frame
  cases_scenarios_references <- dat %>% left_join(df_references)
  cases_scenarios_references <- cases_scenarios_references %>% left_join(df_no_vacc_references)


  # = = = = = 
  # calculate protective efficacy and cases averted for each age year
  # = = = = = 
  # calculate key simulation summary results
  cases_scenarios_each_year <- cases_scenarios_references %>%
    mutate(
      # calculate protective efficacy relative to CM-only (no RTS,S/SMC/IPTi) for each age group
      protective_efficacy = 1 - clinical_cases / clinical_cases_ref,
      protective_efficacy_severe = 1 - severe_cases / severe_cases_ref,
      # calculate protective efficacy relative to no-RTS,S for each age group
      vacc_protective_efficacy = 1 - clinical_cases / clinical_cases_no_vacc,
      vacc_protective_efficacy_severe = 1 - severe_cases / severe_cases_no_vacc,
      # calculate burden relative to no RTS,S/SMC/IPTi for each age group
      relative_burden = (clinical_cases - clinical_cases_ref) / clinical_cases_ref,
      relative_burden_severe = (severe_cases - severe_cases_ref) / severe_cases_ref,
      # calculate burden relative to no RTS,S for each age group
      vacc_relative_burden = (clinical_cases - clinical_cases_no_vacc) / clinical_cases_no_vacc,
      vacc_relative_burden_severe = (severe_cases - severe_cases_no_vacc) / severe_cases_no_vacc,
      # cases and severe cases averted per 100,000 children, relative to CM-only scenario
      cases_averted_per100000 = (clinical_cases_ref - clinical_cases) * 100,
      severe_cases_averted_per100000 = (severe_cases_ref - severe_cases) * 100,
      # cases and severe cases averted per 100,000 children, relative to no-RTS,S scenario
      vacc_cases_averted_per100000 = (clinical_cases_no_vacc - clinical_cases) * 100,
      vacc_severe_cases_averted_per100000 = (severe_cases_no_vacc - severe_cases) * 100,
      age_group = paste(year, '-', (year + 1)), 
      pop = pop
    )


  # = = = = = 
  # calculate protective efficacy and cases averted for each age group
  # = = = = = 
  # aggregate dataframes to U1, U2, U5, and U10
  cases_age_aggregated_list <- list()
  for (i_max in max_years) {
    # get cases per 1000 across all included ages
    tdf = cases_scenarios_references %>%
      filter(year < i_max) %>%
      ungroup() %>%
      dplyr::group_by_at(exp_sweeps) %>%
      dplyr::summarise(clinical_cases = mean(clinical_cases),
                       clinical_cases_ref = mean(clinical_cases_ref),
                       clinical_cases_no_vacc = mean(clinical_cases_no_vacc),
                       severe_cases = mean(severe_cases),
                       severe_cases_ref = mean(severe_cases_ref),
                       severe_cases_no_vacc = mean(severe_cases_no_vacc),
                       pfpr_2_10_no_vacc = mean(pfpr_2_10_no_vacc),
                       pfpr_2_10_ref = mean(pfpr_2_10_ref),
                       pop = mean(pop)) %>%
      mutate(age_group = paste0('U', i_max)) %>%
      dplyr::ungroup()

    # calculate key simulation summary results
    tdf <- tdf %>%
      mutate(
        # calculate protective efficacy relative to CM-only (no RTS,S/SMC/IPTi) for each age group
        protective_efficacy = 1 - clinical_cases / clinical_cases_ref,
        protective_efficacy_severe = 1 - severe_cases / severe_cases_ref,
        # calculate protective efficacy relative to no-RTS,S for each age group
        vacc_protective_efficacy = 1 - clinical_cases / clinical_cases_no_vacc,
        vacc_protective_efficacy_severe = 1 - severe_cases / severe_cases_no_vacc,
        # calculate burden relative to no RTS,S/SMC/IPTi for each age group
        relative_burden = (clinical_cases - clinical_cases_ref) / clinical_cases_ref,
        relative_burden_severe = (severe_cases - severe_cases_ref) / severe_cases_ref,
        # calculate burden relative to no RTS,S for each age group
        vacc_relative_burden = (clinical_cases - clinical_cases_no_vacc) / clinical_cases_no_vacc,
        vacc_relative_burden_severe = (severe_cases - severe_cases_no_vacc) / severe_cases_no_vacc,
        # cases and severe cases averted per 100,000 children, relative to CM-only scenario
        cases_averted_per100000 = (clinical_cases_ref - clinical_cases) * 100,
        severe_cases_averted_per100000 = (severe_cases_ref - severe_cases) * 100,
        # cases and severe cases averted per 100,000 children, relative to no-RTS,S scenario
        vacc_cases_averted_per100000 = (clinical_cases_no_vacc - clinical_cases) * 100,
        vacc_severe_cases_averted_per100000 = (severe_cases_no_vacc - severe_cases) * 100, 
        pop = pop)
    cases_age_aggregated_list[[length(cases_age_aggregated_list) + 1]] <- tdf
  }

  cases_age_aggregated <- bind_rows(cases_age_aggregated_list)


  return(list(cases_age_aggregated, cases_scenarios_each_year))
}











calc_high_low_access_coverages = function(coverage_all, high_access_frac){
  if ((high_access_frac < 1) & (coverage_all >= high_access_frac)){
    coverage_high = 1
    coverage_low = (coverage_all - high_access_frac) / (1 - high_access_frac)
  } else{
    coverage_high = coverage_all / high_access_frac
    coverage_low = 0
  }
  return(c(coverage_high, coverage_low))
}


calc_frac_no_inter = function(inter_names, coverages, targeting, frac_in_high, createVenn = FALSE){
  sample_pop = 10000
  high_access_nums = 1:round(sample_pop*frac_in_high)
  low_access_nums = round(sample_pop*frac_in_high+1):sample_pop
  
  # assign each intervention to people in the population
  xx = list()
  for(ii in 1:length(inter_names)){
    if(targeting[ii]=='high'){
      high_low_cov = calc_high_low_access_coverages(coverage_all=coverages[ii], high_access_frac=frac_in_high)
      # sample from high-access
      cov_ids = sample(high_access_nums, size=high_low_cov[1]*length(high_access_nums), replace=FALSE)
      # sample from low-access
      cov_ids = c(cov_ids, sample(low_access_nums, size=high_low_cov[2]*length(low_access_nums), replace=FALSE))
    } else if (targeting[ii]=='low'){
      low_high_cov = calc_high_low_access_coverages(coverage_all=coverages[ii], high_access_frac=(1-frac_in_high))
      # sample from low-access
      cov_ids = sample(low_access_nums, size=low_high_cov[1]*length(low_access_nums), replace=FALSE)
      # sample from high-access
      cov_ids = c(cov_ids, sample(high_access_nums, size=low_high_cov[2]*length(high_access_nums), replace=FALSE))
    } else{ # random
      cov_ids = sample(1:sample_pop, size=coverages[ii]*sample_pop, replace=FALSE)
    }
    xx[[ii]] = cov_ids
  }
# add venn diagrams of children reached by each intervention
  if(createVenn){
    ggVennDiagram(xx, 
                  label = "percent")
  }
  
  # get fraction not receiving any intervention
  inc_ids = unique(unlist(xx))
  frac_no_inter = (sample_pop-length(inc_ids))/sample_pop
  
  
  # get fraction receiving either RTS,S or SMC
  inc_ids = unique(unlist(xx[which(inter_names %in%  c('smc', 'rtss'))]))
  frac_rtss_or_smc = (length(inc_ids))/sample_pop
  
  # get fraction receiving either RTS,S or SMC
  inc_ids = unique(unlist(xx[which(inter_names %in%  c('cm', 'rtss'))]))
  frac_rtss_or_cm = (length(inc_ids))/sample_pop
  
  # get fraction not previously covered by an intervention that are now covered by RTSS
  smc_cm_ids = unique(unlist(xx[which(inter_names != 'rtss')]))
  no_smc_cm_ids = which(!(1:sample_pop %in% smc_cm_ids))
  now_rtss = intersect(xx[[which(inter_names == 'rtss')]], no_smc_cm_ids)
  frac_newly_covered = length(now_rtss)/length(no_smc_cm_ids)
  
  
  # get fraction not previously covered by SMC that are now covered by RTSS
  smc_ids = unique(xx[[which(inter_names == 'smc')]])
  no_smc_ids = which(!(1:sample_pop %in% smc_ids))
  now_rtss = intersect(xx[[which(inter_names == 'rtss')]], no_smc_ids)
  frac_newly_covered_not_smc = length(now_rtss)/length(no_smc_ids)
  
  
  # get fraction not previously covered by CM that are now covered by RTSS
  cm_ids = unique(xx[[which(inter_names == 'cm')]])
  no_cm_ids = which(!(1:sample_pop %in% cm_ids))
  now_rtss = intersect(xx[[which(inter_names == 'rtss')]], no_cm_ids)
  frac_newly_covered_not_cm = length(now_rtss)/length(no_cm_ids)
  
  
  return(c(frac_no_inter, frac_newly_covered, frac_newly_covered_not_smc, frac_rtss_or_smc, frac_rtss_or_cm))
}


