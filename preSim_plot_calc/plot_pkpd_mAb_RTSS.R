# plot_mAb_efficacy_time

library(ggplot2)
library(viridis)

setwd("C:/Users/moniqueam/Documents/malaria-rtss-mab")
source(file.path('simulation', 'load_paths.R'))

paths = get_project_paths()
datapath = paths[1]
projectpath = paths[2]

##############################################################################
# current parameter eyeballed approximate guestimates
##############################################################################
# monoclonals
mab_initial_concentration = 1000  # 1100
mab_max_efficacy = 0.95
mab_fast_frac=0.722  # 0.7
mab_k1=6.5  # 8
mab_k2=95  # 100
mab_m2 = -0.1
mab_hh=6
mab_nn=1.4

# RTS,S - approximate current functions
rtss1_initial_concentration = 90
rtss1_max_efficacy = 0.8
rtss1_fast_frac=0.88
rtss1_k1=46
rtss1_k2=583
rtss1_m2 = -0.01
rtss1_hh=6
rtss1_nn=1.4

# RTS,S - match concentration decay data
rtss2_fast_frac=0.95
rtss2_k1=120
rtss2_k2=1000
rtss2_m1 = 280

# RTS,S - Imperial parameters from Penny et al
rtss3_fast_frac=0.88
rtss3_k1=46
rtss3_k2=583
rtss3_m1 = 280
# additional parameters to match previous curve along with the concentration params
rtss3_max_efficacy = 0.8
rtss3_initial_concentration = 620
rtss3_hh = 40
rtss3_nn = 1.4


# params after fitting with Phase 3 data
rtss4_fast_frac=0.88
rtss4_k1=46
rtss4_k2=583
rtss4_max_efficacy = 0.8  # for booster, was 0.6
rtss4_initial_concentration = 620
rtss4_hh = 40
rtss4_nn = 2

# time plotted
xx=1:(365*3)

##################################################################################
#--------------------------------------------------------------------------------#
# functions
#--------------------------------------------------------------------------------#
##################################################################################

######################################################################
# function using PKPD to get concentration and efficacy through time
######################################################################
calc_concentration_through_time = function(initial_concentration, fast_frac, k1, k2, xx=seq(1,365*3)){
  # concentration_through_time = c(initial_concentration*exp(-xx[xx<=time_switch]/k1), (initial_concentration*exp(-1*time_switch/k1))*exp(-(xx[xx>time_switch]-time_switch)/k2))
  concentration_through_time = initial_concentration*(fast_frac*exp(-xx / k1) + (1-fast_frac)*exp(-xx / k2))
  return(concentration_through_time)
}

calc_eff_by_conc = function(concentration_values, max_efficacy, hh=NA, nn=NA, m2=NA, hill_func=TRUE){
  if(hill_func){
    efficacy = max_efficacy / (1+(hh/concentration_values)^nn)
  } else{
    efficacy = max_efficacy * (1-exp(m2*concentration_values))
  }
  return(efficacy)
}

calc_effacy_through_time = function(initial_concentration, max_efficacy, fast_frac, k1, k2, m2=NA, hh=NA, nn=NA, hill_func=TRUE, xx=seq(1,365*3), booster_day=NA, create_plot_panel=FALSE){  # time_switch
  concentration_through_time = calc_concentration_through_time(initial_concentration, fast_frac, k1, k2, xx)

  if(!is.na(booster_day)){
    after_booster_initial_concentration = concentration_through_time[booster_day] + initial_concentration
    bb = 1:(max(xx)-booster_day)
    concentration_through_time = c(concentration_through_time[1:booster_day], 
                                   calc_concentration_through_time(initial_concentration=after_booster_initial_concentration, fast_frac=fast_frac, k1=k1, k2=k2, xx=bb))
                                   # after_booster_initial_concentration*(fast_frac*exp(-bb / k1) + (1-fast_frac)*exp(-bb / k2)))
  }
  efficacy_through_time = calc_eff_by_conc(concentration_values=concentration_through_time, max_efficacy=max_efficacy, hh=hh, nn=nn, m2=m2, hill_func=hill_func)

  if(create_plot_panel){
    par(mfrow=c(1,3))
    # concentration through time
    plot(xx, concentration_through_time, type='l', ylim=c(0,initial_concentration), bty='L', xlab='time', main='Concentration through time')
    # Efficacy by concentration
    concentration_sweep = c(seq(0.001,1,length.out=100), seq(1.1, initial_concentration, length.out=1000))
    if(hill_func){
      efficacy_by_concentration = max_efficacy / (1+(hh/concentration_sweep)^nn)
    } else{
      efficacy_by_concentration = max_efficacy * (1-exp(m2*concentration_sweep))
    }
    plot(concentration_sweep, efficacy_by_concentration, ylim=c(0,1), type='l', bty='L', xlab='concentration', log='x', main='Efficacy by concentration')
    # Efficacy through time 
    plot(xx, efficacy_through_time, type='l', ylim=c(0,1), bty='L', xlab='time', main='Efficacy through time')
    par(mfrow=c(1,1))
  }
  return(list(efficacy_through_time, concentration_through_time))
}



######################################################################
# functions for plotting
######################################################################
make_concentration_plot = function(df_mg, fast_frac, k1, k2, m1, title='', conc_col='serum_concentration_ug_ml', time_col='time_weeks', max_days=7*24, include_exp=TRUE){
  initial_concentration = df_mg[[conc_col]][1]
  if(grepl('week', time_col)){
    time_multiplier = 7
  } else if(grepl('year', time_col)){
    time_multiplier = 365
  } else{
    time_multiplier = 1
  }
  plot(df_mg[[time_col]]*time_multiplier, df_mg[[conc_col]], ylim=c(1,1800), pch=20, log='y', xlim=c(0,max_days), bty='L', xlab='time (days)', ylab='concentration (ug/mL)', main=title)
  xx=1:max_days
  if (include_exp){
    # simple exponential
    concentration_through_time = initial_concentration*exp(-xx/m1)
    lines(xx, concentration_through_time, col='blue')
  }
  # current PKPD with biphasic
  concentration_through_time2 =  calc_concentration_through_time(initial_concentration, fast_frac, k1, k2, xx)
  lines(xx, concentration_through_time2, type='l', ylim=c(0,initial_concentration), bty='L', xlab='time', main='Concentration through time', col='red')
}


######################################################################
# functions for finding best parameter fit
######################################################################
get_ref_concTime_distance = function(df_mg, conc_col='serum_concentration_ug_ml', time_col='time_weeks', fast_frac, k1, k2, get_rel_diff=TRUE){
  initial_concentration = df_mg[[conc_col]][1]
  if(grepl('week', time_col)){
    time_multiplier = 7
  } else if(grepl('year', time_col)){
    time_multiplier = 365
  } else{
    time_multiplier = 1
  }
  df_mg$model_val = calc_concentration_through_time(initial_concentration, fast_frac, k1, k2, df_mg[[time_col]]*time_multiplier)
  df_mg$abs_rel_diff = abs((df_mg[[conc_col]] - df_mg$model_val) / df_mg[[conc_col]])
  df_mg$abs_diff = abs((df_mg[[conc_col]] - df_mg$model_val))
  if (get_rel_diff){
    return(mean(df_mg$abs_rel_diff))
  } else{
    return(mean(df_mg$abs_diff))
  }
}

get_conc_infect_lik = function(concentration_df, conc_col='serum_concentration_ug_ml', max_efficacy, hh, nn, num_challenges=1){
  concentration_df$prob_protect = calc_eff_by_conc(concentration_values=concentration_df[[conc_col]], max_efficacy=max_efficacy, hh=hh, nn=nn)
  concentration_df$lik = (concentration_df$prob_protect)^num_challenges
  concentration_df$lik[concentration_df$infected] = (1-(concentration_df$prob_protect[concentration_df$infected])^num_challenges)
  loglik = sum(log(concentration_df$lik))
  return(loglik)
}


fit_and_plot_concTime = function(ref_cur, ref_csvs, ref_names, ref_weights, time_col, conc_col, max_days=7*24, overwrite_file=FALSE){
  if(!file.exists(paste0(projectpath, '/nonSimFigures/mAb_concDecay_fit_', ref_cur, '.csv')) | overwrite_file){
    fast_frac_vals = seq(0.5,1,length.out=10)
    k1_vals = seq(4, 10, length.out=25)
    k2_vals = seq(60, 130, length.out=25)
    ## -> results: best_match = c(5,11,13), corresponding to fast_frac=0.722, k1=6.5, k2=95. min(diff_outputs) = 4.9
    diff_outputs = array(0, dim=c(length(fast_frac_vals), length(k1_vals), length(k2_vals)),
                         dimnames = list(fast_frac = fast_frac_vals,  k1 = k1_vals, k2 = k2_vals))
    for(i1 in 1:length(fast_frac_vals)){
      for(i2 in 1:length(k1_vals)){
        for(i3 in 1:length(k2_vals)){
          for(rr in 1:length(ref_csvs)){
            df_mg = read.csv(ref_csvs[rr])
            cur_diff = get_ref_concTime_distance(df_mg, conc_col='serum_concentration_ug_ml', time_col=time_col, fast_frac=fast_frac_vals[i1], k1=k1_vals[i2], k2=k2_vals[i3], get_rel_diff=TRUE)
            diff_outputs[i1,i2,i3] = diff_outputs[i1,i2,i3] + ref_weights[rr] * cur_diff
          }
        }
      }
    }
    # create dataframe with param values and difference values
    diff_df = as.data.frame.table(diff_outputs)
    indx = sapply(diff_df, is.factor)
    diff_df[indx] = lapply(diff_df[indx], function(x) as.numeric(as.character(x)))
    diff_df_sorted = diff_df[order(diff_df$Freq, decreasing=FALSE),]
    colnames(diff_df_sorted)[colnames(diff_df_sorted)=='Freq'] = 'difference'
    
    write.csv(diff_df_sorted, file=paste0(projectpath, '/nonSimFigures/mAb_concDecay_fit_', ref_cur, '.csv'), row.names=FALSE)
  } else{
    diff_df_sorted = read.csv(paste0(projectpath, '/nonSimFigures/mAb_concDecay_fit_', ref_cur, '.csv'))
  }
  
  # create plot of fit to data
  if(grepl('week', time_col)){
    time_multiplier = 7
  } else if(grepl('year', time_col)){
    time_multiplier = 365
  } else{
    time_multiplier = 1
  }
  xx=1:max_days
  png(filename=paste0(projectpath, '/nonSimFigures/mAb_concDecay_', ref_cur, '.png'), width=6, height=5, res=900, units='in')
  par(mfrow=c(2,2))
  nvals=10
  colors = rainbow(nvals)
  colors[1]='black'
  for(rr in 1:length(ref_csvs)){
    df_mg = read.csv(ref_csvs[rr])
    initial_concentration = df_mg[[conc_col]][1]
    # plot reference data
    plot(df_mg[[time_col]]*time_multiplier, df_mg[[conc_col]], ylim=c(1,1800), pch=20, log='y', xlim=c(0,max_days), bty='L', xlab='time (days)', ylab='concentration (ug/mL)', main=ref_names[rr])
    # current PKPD with biphasic decay for top parameter sets
    for(ll in 1:nvals){
      concentration_through_time2 =  calc_concentration_through_time(initial_concentration, fast_frac=diff_df_sorted$fast_frac[ll], k1=diff_df_sorted$k1[ll], k2=diff_df_sorted$k2[ll], xx)
      lines(xx, concentration_through_time2, type='l', lwd=1, col=colors[ll])
    }
    if(rr==length(ref_csvs)){
      legend('left', bty='n', paste0('ff ', round(diff_df_sorted$fast_frac[1:nvals],2), '; k1 ', round(diff_df_sorted$k1[1:nvals],1), '; k2 ', round(diff_df_sorted$k2[1:nvals],2)), col=colors, lwd=1, cex=0.5)
    }
  }
  dev.off()
  par(mfrow=c(1,1))
}


##################################################################################
#--------------------------------------------------------------------------------#
# mAb fitting
#--------------------------------------------------------------------------------#
##################################################################################

##############################################################################
# estimate concentration through time (parameters fast_frac, k1, k2)
##############################################################################
# since I don't have actual data, just plot-grabbed values for means, use distance between model and mean datapoint, weighted by number of individuals

max_days = 7*24

# includes two reference datasets: Gaudinski et al with CIS43LS (for concentration-through-time only, but longer period of time) and Wu et al with L9LS (both concentration through time and efficacy by concentration)
ref_cur = 'Gaudinski'
ref_csvs = paste0(datapath, '/vacc_pkpd/', c('rough_dataGrab_5mgSC_Gaudinski_SF1.csv', 'rough_dataGrab_5mg_Gaudinski_SF1.csv', 'rough_dataGrab_20mg_Gaudinski_SF1.csv', 'rough_dataGrab_40mg_Gaudinski_SF1.csv'))
ref_names = c('5mg/kg (SC)', '5mg/kg (IV)', '20mg/kg (IV)', '40mg/kg (IV)')
ref_weights = c(3, 4, 5, 6)
time_col = 'time_weeks'
conc_col = 'serum_concentration_ug_ml'
fit_and_plot_concTime(ref_cur, ref_csvs, ref_names, ref_weights, time_col, conc_col, max_days=max_days, overwrite_file=FALSE)
ref_csvs_Gaudinski = ref_csvs
time_col_Gaudinski = time_col

ref_cur = 'Wu'
ref_csvs = paste0(datapath, '/Wu et al_2022/', c('Wu_Fig3_1iv.csv', 'Wu_Fig3_5sc.csv', 'Wu_Fig3_5iv.csv', 'Wu_Fig3_20iv.csv'))
ref_names = c('1mg/kg (IV)','5mg/kg (SC)', '5mg/kg (IV)', '20mg/kg (IV)')
ref_weights = c(5, 1, 5, 5)
time_col='time_days'
conc_col = 'serum_concentration_ug_ml'
fit_and_plot_concTime(ref_cur, ref_csvs, ref_names, ref_weights, time_col, conc_col, max_days=max_days, overwrite_file=FALSE)
ref_csvs_Wu = ref_csvs
time_col_Wu = time_col

# - - - - - - - - - - - - - - - - #
# compare Guadinski and Wu data
# - - - - - - - - - - - - - - - - #
# create plot of fit to data
png(filename=paste0(projectpath, '/nonSimFigures/mAb_concDecay_compareRefs.png'), width=6, height=5, res=900, units='in')
par(mfrow=c(1,1))
plot(NA, xlim=c(0,max_days), ylim=c(1,1800), log='y', bty='L', xlab='time (days)', ylab='concentration (ug/mL)', main='Comparison of CIS43LS and L9LS Phase I')
for(rr in 1:4){
  df_mg_g = read.csv(ref_csvs_Gaudinski[rr])
  df_mg_w = read.csv(ref_csvs_Wu[rr])
  # plot reference data
  lines(df_mg_g[[time_col_Gaudinski]]*7, df_mg_g[[conc_col]], col='blue')
  points(df_mg_g[[time_col_Gaudinski]]*7, df_mg_g[[conc_col]], col='blue', pch=20)
  lines(df_mg_w[[time_col_Wu]]*1, df_mg_w[[conc_col]], col='red')
  points(df_mg_w[[time_col_Wu]]*1, df_mg_w[[conc_col]], col='red', pch=20)
}
legend('topright', c('Gaudinski et al (CIS43LS)', 'Wu et al (L9LS)'), col=c('blue', 'red'), pch=20, bty='n')
dev.off()




##############################################################################
# estimate efficacy-by-concentration parameters (max efficacy, hh, nn)
##############################################################################
concentration_df = read.csv(paste0(datapath, '/Wu et al_2022/Wu_Fig3_concentrations_at_challenge.csv'))
overwrite_file = FALSE
num_challenges = 3
if(!file.exists(paste0(projectpath, '/nonSimFigures/mAb_effConc_fit_Wu_', num_challenges, 'challenges.csv')) | overwrite_file){
  max_eff_vals = c(0.8, 0.85, 0.9, 0.95, 0.98)
  hh_vals = seq(1,12,length.out=30)
  nn_vals = seq(0.1,10, length.out=30)
  loglik_outputs = array(0, dim=c(length(max_eff_vals), length(hh_vals), length(nn_vals)),
                         dimnames = list(max_efficacy = max_eff_vals,  hh = hh_vals, nn = nn_vals))
  for(i1 in 1:length(max_eff_vals)){
    for(i2 in 1:length(hh_vals)){
      for(i3 in 1:length(nn_vals)){
        loglik_outputs[i1, i2, i3] = get_conc_infect_lik(concentration_df, max_efficacy = max_eff_vals[i1], hh=hh_vals[i2], nn=nn_vals[i3], num_challenges=num_challenges)
      }
    }
  }
  # create dataframe with param values and  loglikelihood values
  loglik_df = as.data.frame.table(loglik_outputs)
  indx = sapply(loglik_df, is.factor)
  loglik_df[indx] = lapply(loglik_df[indx], function(x) as.numeric(as.character(x)))
  loglik_df_sorted = loglik_df[order(loglik_df$Freq, decreasing=TRUE),]
  colnames(loglik_df_sorted)[colnames(loglik_df_sorted)=='Freq'] = 'loglik'
  write.csv(loglik_df_sorted, file=paste0(projectpath, '/nonSimFigures/mAb_effConc_fit_Wu_', num_challenges, 'challenges.csv'), row.names=FALSE)
} else{
  loglik_df_sorted = read.csv(paste0(projectpath, '/nonSimFigures/mAb_effConc_fit_Wu_', num_challenges, 'challenges.csv'))
}
# create plot showing comparison of top fits against data
png(filename=paste0(projectpath, '/nonSimFigures/mAb_effConc_wData_fit_Wu_', num_challenges, 'challenges.png'), width=5, height=4, res=900, units='in')
concs = 1:ceiling(max(concentration_df$serum_concentration_ug_ml, na.rm=TRUE))
par(mfrow=c(1,1))
nvals=10
colors = rainbow(nvals)
colors[1]='black'
plot(concs, calc_eff_by_conc(concentration_values=concs, max_efficacy=loglik_df_sorted$max_efficacy[1], hh=loglik_df_sorted$hh[1], nn=loglik_df_sorted$nn[1]),
     type='l', lwd=2, col=colors[1], log='x', ylim=c(0,1), ylab = 'probability of protection', xlab='concentration (ug/ml)', bty='L')
for(ll in 2:nvals){
  lines(concs, calc_eff_by_conc(concentration_values=concs, max_efficacy=loglik_df_sorted$max_efficacy[ll], hh=loglik_df_sorted$hh[ll], nn=loglik_df_sorted$nn[ll]),
        type='l', lwd=1, col=colors[ll])
}
points(concentration_df$serum_concentration_ug_ml, abs(1-concentration_df$infected), pch=20)
legend('right', bty='n', paste0('ME ', round(loglik_df_sorted$max_efficacy[1:nvals],2), '; hh ', round(loglik_df_sorted$hh[1:nvals],1), '; nn ', round(loglik_df_sorted$nn[1:nvals],2)), col=colors, lwd=1, cex=0.85)
dev.off()
  


########################################
# plot efficacy through time for different starting concentrations using all fit parameters from Wu
########################################
diff_df_sorted_Wu = read.csv(paste0(projectpath, '/nonSimFigures/mAb_concDecay_fit_Wu.csv'))
loglik_df_sorted = read.csv(paste0(projectpath, '/nonSimFigures/mAb_effConc_fit_Wu_', num_challenges, 'challenges.csv'))
initial_concentrations=c(50, 100, 1000); ltys = c(1,2, 3)   # c(100, 1000)
png(filename=paste0(projectpath, '/nonSimFigures/effTime_mAb_topFits_Wu_', num_challenges, 'challenges.png'), width=6, height=5, res=900, units='in')
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='days', ylab='efficacy at preventing infection', main='Efficacy against infection through time')
xx=1:(365*2)
nvals=10
colors = rainbow(nvals)
colors[1]='black'
for(i1 in 1:length(initial_concentrations)){
  for(ll_index in 1:nvals){
    cc_index = 1  # ll_index  # index to use for concentration-through-time fit parameters (probably either 1 or ll_index)
    lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=loglik_df_sorted$max_efficacy[ll_index], fast_frac=diff_df_sorted_Wu$fast_frac[cc_index], k1=diff_df_sorted_Wu$k1[cc_index], k2=diff_df_sorted_Wu$k2[cc_index], hh=loglik_df_sorted$hh[ll_index], nn=loglik_df_sorted$nn[ll_index], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
          col=colors[ll_index], lty=ltys[i1])
  }
}
lines(xx, calc_effacy_through_time(initial_concentration=rtss4_initial_concentration, max_efficacy=rtss4_max_efficacy, fast_frac=rtss4_fast_frac, k1=rtss4_k1, k2=rtss4_k2, nn=rtss4_nn, hh=rtss4_hh, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
      col=rgb(0,0,0,0.2), lwd=3.5, lty=1)
legend('topright', bty='n', paste0('ME ', round(loglik_df_sorted$max_efficacy[1:nvals],2), '; hh ', round(loglik_df_sorted$hh[1:nvals],1), '; nn ', round(loglik_df_sorted$nn[1:nvals],2)), col=colors, lwd=1, cex=0.85)
legend('bottomleft', bty='n', paste0('initial concentration: ', initial_concentrations, ' ug/ml'), col='black', lty=ltys, lwd=1, cex=0.85)
dev.off()


# separate lines for concentration fits
diff_df_sorted_Wu = read.csv(paste0(projectpath, '/nonSimFigures/mAb_concDecay_fit_Wu.csv'))
loglik_df_sorted = read.csv(paste0(projectpath, '/nonSimFigures/mAb_effConc_fit_Wu_', num_challenges, 'challenges.csv'))
initial_concentrations=c(20, 100, 1000); ltys = c(1,2)
png(filename=paste0(projectpath, '/nonSimFigures/compare_GR_sweep_with_Wu_fit_', num_challenges, 'challenges.png'), width=6, height=5, res=900, units='in')
nvals=10
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='days', ylab='efficacy at preventing infection', main='Efficacy against infection through time')
colors = c('red','orange', 'yellow')
xx=1:(365*2)
for(i1 in 1:length(initial_concentrations)){
  for(ll_index in 1:nvals){
    for(cc_index in 1:nvals){  # ll_index  # index to use for concentration-through-time fit parameters (probably either 1 or ll_index)
      lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=loglik_df_sorted$max_efficacy[ll_index], fast_frac=diff_df_sorted_Wu$fast_frac[cc_index], k1=diff_df_sorted_Wu$k1[cc_index], k2=diff_df_sorted_Wu$k2[cc_index], hh=loglik_df_sorted$hh[ll_index], nn=loglik_df_sorted$nn[ll_index], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
            col=colors[i1], lty=1)
    }
  }
}
# single line from 'best' match
lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[2], max_efficacy=loglik_df_sorted$max_efficacy[1], fast_frac=diff_df_sorted_Wu$fast_frac[1], k1=diff_df_sorted_Wu$k1[1], k2=diff_df_sorted_Wu$k2[1], hh=loglik_df_sorted$hh[1], nn=loglik_df_sorted$nn[1], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
      col='black', lty=1)
# RTS,S
lines(xx, calc_effacy_through_time(initial_concentration=rtss4_initial_concentration, max_efficacy=rtss4_max_efficacy, fast_frac=rtss4_fast_frac, k1=rtss4_k1, k2=rtss4_k2, nn=rtss4_nn, hh=rtss4_hh, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
      col=rgb(0,0,0,0.2), lwd=3.5, lty=1)
legend('bottomleft', bty='n', paste0('Wu fits - initial concentration: ', initial_concentrations, ' ug/ml'),  col=colors, lwd=1, cex=0.85)

# GR values
initial_concentrations=c(1000)
lwds=c(1.5)
max_efficacies=c(0.8, 0.9, 0.95, 0.98)
ltys = c(4, 3, 2, 1)
fast_frac=0.722
k1=6.5
k2=95
nns = c(1.4)
hhs = c(2, 5, 10, 20, 40, 60)   # c(5, 10, 20, 40, 60)     c(10, 20, 40) 
colors = viridis(length(hhs)+1)
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 1:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=colors[i4], lwd=lwds[i1], lty=ltys[i2])
      }
    }
  }
}

dev.off()



##################################################################################
#--------------------------------------------------------------------------------#
# GR (Sept 2022) plots: efficacy-through-time plots
#--------------------------------------------------------------------------------#
##################################################################################

xx = 1:(365*2)
# first set of PKPD profiles for GR (9/10/2022)
initial_concentrations=c(1000)
lwds=c(1.5)
max_efficacies=c(0.8, 0.9, 0.95) #, 0.99)
ltys = c(3, 2, 1)
fast_frac=0.722
k1=6.5
k2=95
nns = c(1.4)
hhs = c(2, 5, 10, 20, 40, 60)   # c(5, 10, 20, 40, 60)     c(10, 20, 40) 
colors = viridis(length(hhs)+1)
png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_plannedParams4b.png'), width=4, height=4, res=900, units='in')
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='days', ylab='efficacy at preventing infection')
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 1:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=colors[i4], lwd=lwds[i1], lty=ltys[i2])
      }
    }
  }
}
lines(xx, calc_effacy_through_time(initial_concentration=rtss4_initial_concentration, max_efficacy=rtss4_max_efficacy, fast_frac=rtss4_fast_frac, k1=rtss4_k1, k2=rtss4_k2, nn=rtss4_nn, hh=rtss4_hh, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
      col=rgb(0,0,0,0.2), lwd=3.5, lty=1)
# legend(x=500, y=1, legend=hhs, lty=1, col=colors, title='EC50', bty='n')
# legend(x=500, y=0.5, legend=max_efficacies, lty=ltys, col='black', title='max efficacy', bty='n')
dev.off()


# build without RTSS
png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_plannedParams4_noRTSS.png'), width=4, height=4, res=900, units='in')
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='days', ylab='efficacy at preventing infection')
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 1:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=colors[i4], lwd=lwds[i1], lty=ltys[i2])
      }
    }
  }
}
dev.off()


# build with only one max efficacy not transparent
png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_plannedParams4_oneMaxEff.png'), width=4, height=4, res=900, units='in')
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='days', ylab='efficacy at preventing infection')
colors = viridis(length(hhs)+1, alpha=0.15)
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 1:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=colors[i4], lwd=lwds[i1], lty=ltys[i2])
      }
    }
  }
}
colors = viridis(length(hhs)+1)
for(i1 in 1:length(initial_concentrations)){
  for(i2 in length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 1:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=colors[i4], lwd=lwds[i1], lty=ltys[i2])
      }
    }
  }
}
dev.off()


# build with only one hh (EC50) not transparent
png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_plannedParams4_oneEC50.png'), width=4, height=4, res=900, units='in')
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='days', ylab='efficacy at preventing infection')
colors = viridis(length(hhs)+1, alpha=0.15)
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 1:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=colors[i4], lwd=lwds[i1], lty=ltys[i2])
      }
    }
  }
}
colors = viridis(length(hhs)+1)
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 3){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=colors[i4], lwd=lwds[i1], lty=ltys[i2])
      }
    }
  }
}
dev.off()




# build with only one hh (EC50) and max_efficacy to represent red dot example
png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_plannedParams4_oneME_oneEC50.png'), width=4, height=4, res=900, units='in')
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='days', ylab='efficacy at preventing infection')
for(i1 in 1:length(initial_concentrations)){
  for(i2 in length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 3){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col='red', lwd=3, lty=ltys[i2])
      }
    }
  }
}
dev.off()





# second set of PKPD profiles for Sept GR (9/4/2022)
initial_concentrations=c(1000)
lwds=c(3.5)
max_efficacies=c(0.95)
ltys = c(1)
fast_frac=0.722
k1=6.5
k2=95
nns = c(1.4)
hhs = c(5, 10, 20, 40, 60)
colors = viridis(length(hhs)+1)
png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_plannedParams5_GR3.png'), width=3.25, height=3.5, res=900, units='in')
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='days', ylab='efficacy at preventing infection')
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 2:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=colors[i4], lwd=lwds[i1], lty=ltys[i2])
      }
    }
  }
}
dev.off()




# third set of PKPD profiles for Sept GR (9/4/2022)
initial_concentrations=c(1000)
lwds=c(3.5)
max_efficacies=c(0.95)
ltys = c(1)
fast_frac=0.722
k1=6.5
k2=95
nns = c(1.4)
hhs = c(5, 10, 20, 40, 60)
colors = viridis(length(hhs)+1)
png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_plannedParams5_GR4.png'), width=3.25, height=3.5, res=900, units='in')
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='days', ylab='efficacy at preventing infection')
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 2:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=colors[i4], lwd=lwds[i1], lty=ltys[i2])
      }
    }
  }
}
lines(xx, calc_effacy_through_time(initial_concentration=1000, max_efficacy=0.995, fast_frac=fast_frac, k1=k1, k2=k2, nn=2, hh=5, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
      col="#b6dbff", lwd=lwds[1], lty=1)
# legend(x=500, y=1, legend=hhs, lty=1, col=colors, title='EC50', bty='n')
# legend(x=500, y=0.5, legend=max_efficacies, lty=ltys, col='black', title='max efficacy', bty='n')
dev.off()




# example generic efficacy-through-time plot for GR (9/7/2022)
initial_concentrations=c(1000)
lwds=c(1.5)
max_efficacies=c( 0.95)
ltys = c(1)
fast_frac=0.722
k1=6.5
k2=95
nns = c(1.4)
hhs = c(10)  
png(filename=paste0(projectpath, '/nonSimFigures/example_efficacy_through_Time.png'), width=3, height=3, res=900, units='in')
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='days', ylab='efficacy at preventing infection')
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 1:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col='black', lwd=lwds[i1], lty=ltys[i2])
      }
    }
  }
}
dev.off()








xx = 1:(365*2)
# new set of PKPD profiles for GR, with earlier decay instead of starting with long plateau (smaller nn) (9/8/2022)
initial_concentrations=c(1000)
lwds=c(1.5)
max_efficacies=c(0.8, 0.9, 0.95)
ltys = c(3, 2, 1)
fast_frac=0.722
k1=6.5
k2=95
nns = c(1.4)
hhs = c(5, 10, 20, 40, 60)   # c(5, 10, 20, 40, 60)     c(10, 20, 40) 
colors = viridis(length(hhs)+1)
png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_plannedParams4bTEST2.png'), width=4, height=4, res=900, units='in')
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='days', ylab='efficacy at preventing infection')
for(i1 in 1:length(initial_concentrations)){
for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 1:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=colors[i4], lwd=lwds[i1], lty=ltys[i2])
  }
  }
  }
}
lines(xx, calc_effacy_through_time(initial_concentration=rtss4_initial_concentration, max_efficacy=rtss4_max_efficacy, fast_frac=rtss4_fast_frac, k1=rtss4_k1, k2=rtss4_k2, nn=rtss4_nn, hh=rtss4_hh, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
      col=rgb(0,0,0,0.2), lwd=3.5, lty=1)
# legend(x=500, y=1, legend=hhs, lty=1, col=colors, title='EC50', bty='n')
# legend(x=500, y=0.5, legend=max_efficacies, lty=ltys, col='black', title='max efficacy', bty='n')
dev.off()


##################################################################################
#--------------------------------------------------------------------------------#
# other plots
#--------------------------------------------------------------------------------#
##################################################################################

###################################################################################################################################################
# RTS,S protective efficacy through time from prior params versus pkpd functions to closely match current shape
###################################################################################################################################################
png(filename=paste0(projectpath, '/nonSimFigures/emod_function_rtss.png'), width=5.5, height=4.5, res=900, units='in')
plot(NA, ylim=c(0,1), xlim=c(0,max(xx)), ylab='protective efficacy', xlab='days since dose', bty='L')
# values described in RTSS report
init=0.8; dur=410; color='grey'
lines(init*exp(xx*-1/dur), ylab='protective efficacy', xlab='days since dose', bty='L', lwd=2, col=color)
# values used in RTSS analysis
init=0.8; dur=592; color='black'
lines(init*exp(xx*-1/dur), ylab='protective efficacy', xlab='days since dose', bty='L', lwd=2, col=color)
# imitate with PKPD from Phase 3
initial_concentration = rtss3_initial_concentration
max_efficacy = rtss3_max_efficacy
fast_frac=rtss3_fast_frac
k1=rtss3_k1
k2=rtss3_k2
hh=40
nn=2
eff_conc = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, hh=hh, nn=nn, hill_func=TRUE, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]]
lines(eff_conc, lwd=2, col='salmon')
legend('topright', c('RTSS analysis (initial=0.8, duration=592)', 'RTSS report (initial=0.8, duration=410)', 'selected from Phase 3 comparisons'), 
       col=c('black', 'grey', 'salmon'), lwd=2, bty='n')
dev.off()



##################################
# PKPD examples for monoclonals
##################################
# what it should look like given concentration through time and efficacy by concentration (PKPD dynamics)
# set parameters for old concentration-efficacy and Hill-function-like concentration-efficacy curves
initial_concentration = mab_initial_concentration
fast_frac=mab_fast_frac
k1=mab_k1
k2=mab_k2
max_efficacy = mab_max_efficacy
m2 = mab_m2
m2b = -2
nn=mab_nn
hh=mab_hh
hh2=0.3


# time plotted
xx=1:(365*3)
png(filename=paste0(projectpath, '/nonSimFigures/pkpd_function_example_oldFunc.png'), width=6, height=2.5, res=900, units='in')
eff_conc = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, m2=m2, hill_func=FALSE, xx=xx, booster_day=NA, create_plot_panel=TRUE)
dev.off()
png(filename=paste0(projectpath, '/nonSimFigures/pkpd_function_example_oldFunc_m2b.png'), width=6, height=2.5, res=900, units='in')
eff_conc = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, m2=m2b, hill_func=FALSE, xx=xx, booster_day=NA, create_plot_panel=TRUE)
dev.off()
png(filename=paste0(projectpath, '/nonSimFigures/pkpd_function_example_hill.png'), width=6, height=2.5, res=900, units='in')
eff_conc = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, nn=nn, hh=hh, hill_func=TRUE, xx=xx, booster_day=NA, create_plot_panel=TRUE)
dev.off()
png(filename=paste0(projectpath, '/nonSimFigures/pkpd_function_example_hill2.png'), width=6, height=2.5, res=900, units='in')
eff_conc = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, nn=nn, hh=hh2, hill_func=TRUE, xx=xx, booster_day=NA, create_plot_panel=TRUE)
dev.off()



#################################################
# will exponential decay look too different? Yes.
#################################################
initial_concentration = mab_initial_concentration
fast_frac=mab_fast_frac
k1=mab_k1
k2=mab_k2
max_efficacy = mab_max_efficacy
m2 = mab_m2
nn=mab_nn
hh=mab_hh

png(filename=paste0(projectpath, '/nonSimFigures/emod_function_mab.png'), width=5.5, height=4.5, res=900, units='in')
plot(NA, ylim=c(0,1), xlim=c(0,max(xx)), ylab='protective efficacy', xlab='days since dose', bty='L')

# try to represent with exponential decay
init=max_efficacy; dur=365*3; color='orange'
lines(init*exp(xx*-1/dur), ylab='protective efficacy', xlab='days since dose', bty='L', lwd=2, col=color)
init=max_efficacy; dur=592; color='darkred'
lines(init*exp(xx*-1/dur), ylab='protective efficacy', xlab='days since dose', bty='L', lwd=2, col=color)
init=max_efficacy; dur=410; color='black'
lines(init*exp(xx*-1/dur), ylab='protective efficacy', xlab='days since dose', bty='L', lwd=2, col=color)
init=max_efficacy; dur=200; color='grey'
lines(init*exp(xx*-1/dur), ylab='protective efficacy', xlab='days since dose', bty='L', lwd=2, col=color)

# manually specify decay-of-efficacy vector
color='dodgerblue'
efficacy_time = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, nn=nn, hh=hh, hill_func=TRUE, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]]
lines(efficacy_time, col=color, lwd=2)

legend('topright', c('exponential (duration=1095)', 'exponential (duration=592)', 'exponential (duration=410)', 'exponential (duration=200)', 'E(T)'), 
       col=c('orange', 'darkred', 'black', 'grey', 'dodgerblue'), lwd=2, bty='n')
dev.off()



##############################################################################
# compare exponential and biphasic decay against RTSS anti-CS antibody data
##############################################################################
# RTSS params from eyeballing reference plots

fast_frac=rtss2_fast_frac
k1=rtss2_k1
k2=rtss2_k2
m1 = rtss2_m1
# compare roughly data-extracted values with simple exponential versus biphasic decay
conc_col='titre_eu_ml'
png(filename=paste0(projectpath, '/nonSimFigures/RTSS_concDecay_wData.png'), width=9, height=4.5, res=900, units='in')
par(mfrow=c(1,2))
max_days = round(365*3.5)
# df_mg = read.csv(paste0(datapath, '/vacc_pkpd/White_517_boost.csv'))
# make_concentration_plot(df_mg, fast_frac, k1, k2, m1, title='5-17 months (booster)', conc_col='titre_eu_ml', time_col='time_years', max_days=max_days)
df_mg = read.csv(paste0(datapath, '/vacc_pkpd/White_517_noBoost.csv'))
make_concentration_plot(df_mg, fast_frac=rtss2_fast_frac, k1=rtss2_k1, k2=rtss2_k2, m1=rtss2_m1, title='5-17 months (no booster)', conc_col='titre_eu_ml', time_col='time_years', max_days=max_days)
# add Imperial curves
initial_concentration = df_mg[[conc_col]][1]
concentration_through_time2 =  calc_concentration_through_time(initial_concentration, fast_frac=rtss3_fast_frac, k1=rtss3_k1, k2=rtss3_k2, xx=1:max_days)
lines(1:max_days, concentration_through_time2, type='l', col='darkred')
legend('topright', c('data','exponential', 'biphasic', 'Imperial'), pch=c(20,NA,NA,NA), lwd=c(NA,1,1,1), col=c('black','blue','red', 'darkred'), bty='n')
# df_mg = read.csv(paste0(datapath, '/vacc_pkpd/White_612_boost.csv'))
# make_concentration_plot(df_mg, fast_frac, k1, k2, m1, title='6-12 weeks (booster)', conc_col='titre_eu_ml', time_col='time_years', max_days=max_days)
df_mg = read.csv(paste0(datapath, '/vacc_pkpd/White_612_noBoost.csv'))
make_concentration_plot(df_mg, fast_frac=rtss2_fast_frac, k1=rtss2_k1, k2=rtss2_k2, m1=rtss2_m1, title='6-12 weeks (no booster)', conc_col='titre_eu_ml', time_col='time_years', max_days=max_days)
# add Imperial curves
initial_concentration = df_mg[[conc_col]][1]
concentration_through_time2 =  calc_concentration_through_time(initial_concentration, fast_frac=rtss3_fast_frac, k1=rtss3_k1, k2=rtss3_k2, xx=1:max_days)
lines(1:max_days, concentration_through_time2, type='l', col='darkred')
par(mfrow=c(1,1))
dev.off()

# RTSS params from Imperial values in Penny et al 2016

fast_frac=rtss3_fast_frac
k1=rtss3_k1
k2=rtss3_k2
m1 = rtss3_m1
# compare roughly data-extracted values with simple exponential versus biphasic decay
png(filename=paste0(projectpath, '/nonSimFigures/RTSS_concDecay_wData_ImperialParams.png'), width=6, height=3, res=900, units='in')
par(mfrow=c(1,2))
max_days = round(365*3.5)
# df_mg = read.csv(paste0(datapath, '/vacc_pkpd/White_517_boost.csv'))
# make_concentration_plot(df_mg, fast_frac, k1, k2, m1, title='5-17 months (booster)', conc_col='titre_eu_ml', time_col='time_years', max_days=max_days)
df_mg = read.csv(paste0(datapath, '/vacc_pkpd/White_517_noBoost.csv'))
make_concentration_plot(df_mg, fast_frac, k1, k2, m1, title='5-17 months (no booster)', conc_col='titre_eu_ml', time_col='time_years', max_days=max_days)
legend('topleft', c('data','exponential', 'biphasic'), pch=c(20,NA,NA), lwd=c(NA,1,1), col=c('black','blue','red'), bty='n')
# df_mg = read.csv(paste0(datapath, '/vacc_pkpd/White_612_boost.csv'))
# make_concentration_plot(df_mg, fast_frac, k1, k2, m1, title='6-12 weeks (booster)', conc_col='titre_eu_ml', time_col='time_years', max_days=max_days)
df_mg = read.csv(paste0(datapath, '/vacc_pkpd/White_612_noBoost.csv'))
make_concentration_plot(df_mg, fast_frac, k1, k2, m1, title='6-12 weeks (no booster)', conc_col='titre_eu_ml', time_col='time_years', max_days=max_days)
par(mfrow=c(1,1))
dev.off()





##############################################################################
# Efficacy through time given sets of mAb params
##############################################################################

initial_concentrations=c(500,1000,2000)
max_efficacies=c(0.8,0.95)
fast_frac=mab_fast_frac
k1=mab_k1
k2=mab_k2
nns = c(1.4,2,4)
hhs = c(1, 6,30)
color_set=list(c(rgb(1,0,0), rgb(0,0,1), rgb(0.3,0.3,0.3)),c(rgb(1,0.5,0.5), rgb(0.5,0.5,1), rgb(0.7,0.7,0.7)))
lwds = c(1,2)
ltys=c(1,2, 3)
# png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_OldParams.png'), width=6, height=5, res=900, units='in')
plot(NA, xlim=c(0,max(xx)), ylim=c(0,1), bty='L', xlab='time', ylab='efficacy', main='Alternative TPPs')
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 1:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=color_set[[i4]][i3], lwd=lwds[i2], lty=ltys[i1])
      }
    }
  }
}
# dev.off()


# parameter sets currently (7/29/2022) planned for Sept GR analyses
initial_concentrations=c(1000)
lwds=c(1)
max_efficacies=c(0.8, 0.9, 0.95)
ltys = c(3, 2, 1)
fast_frac=0.7
k1=8
k2=100
nns = c(2)
hhs = c(10, 20, 40)   # c(5, 10, 20, 40, 60)     c(10, 20, 40) 
colors = viridis(length(hhs))
png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_plannedParams2.png'), width=6, height=5, res=900, units='in')
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='time', ylab='efficacy', main='Alternative TPPs')
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 1:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=colors[i4], lwd=lwds[i1], lty=ltys[i2])
      }
    }
  }
}
lines(xx, calc_effacy_through_time(initial_concentration=rtss4_initial_concentration, max_efficacy=rtss4_max_efficacy, fast_frac=rtss4_fast_frac, k1=rtss4_k1, k2=rtss4_k2, nn=rtss4_nn, hh=rtss4_hh, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
      col=rgb(0,0,0,0.2), lwd=3.5, lty=1)
legend(x=500, y=1, legend=hhs, lty=1, col=colors, title='hh', bty='n')
legend(x=500, y=0.5, legend=max_efficacies, lty=ltys, col='black', title='max efficacy', bty='n')
dev.off()

# parameter sets currently (7/29/2022) planned for Sept GR analyses
initial_concentrations=c(1000)
max_efficacies=c(0.8, 0.9, 0.95)
ltys = c(3, 2, 1)
fast_frac=0.7
k1=8
k2=100
nns = c(0.5, 2)
lwds=c(1, 2)
hhs = c(10, 20, 40)   # c(5, 10, 20, 40, 60)     c(10, 20, 40) 
colors = viridis(length(hhs))
png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_otherPossParams.png'), width=6, height=5, res=900, units='in')
plot(NA, xlim=c(0,365*2), ylim=c(0,1), bty='L', xlab='time', ylab='efficacy', main='Alternative TPPs')
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 1:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              col=colors[i4], lwd=lwds[i3], lty=ltys[i2])
      }
    }
  }
}
lines(xx, calc_effacy_through_time(initial_concentration=rtss4_initial_concentration, max_efficacy=rtss4_max_efficacy, fast_frac=rtss4_fast_frac, k1=rtss4_k1, k2=rtss4_k2, nn=rtss4_nn, hh=rtss4_hh, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
      col=rgb(0,0,0,0.2), lwd=3.5, lty=1)
legend(x=500, y=1, legend=hhs, lty=1, col=colors, title='hh', bty='n')
legend(x=500, y=0.5, legend=max_efficacies, lty=ltys, col='black', title='max efficacy', bty='n')
dev.off()



initial_concentrations=c(100, 1000, 2000)
max_efficacies=c(0.8)
fast_frac=mab_fast_frac
k1=mab_k1
k2=mab_k2
nns = c(0.5, 2, 4, 8)
hhs = c(0.1, 10, 20, 40)
color_set=list(c(rgb(1,0,0), rgb(0,0,1), rgb(0.3,0.3,0.3)),c(rgb(1,0.5,0.5), rgb(0.5,0.5,1), rgb(0.7,0.7,0.7)))
lwds = c(1,2)
ltys=c(1,2, 3)
# png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_OldParams.png'), width=6, height=5, res=900, units='in')
plot(NA, xlim=c(0,max(xx)), ylim=c(0,1), bty='L', xlab='time', ylab='efficacy', main='Alternative TPPs')
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(nns)){
      for(i4 in 1:length(hhs)){
        lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, nn=nns[i3], hh=hhs[i4], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
              # col=color_set[[i4]][i3], lwd=lwds[i2], lty=ltys[i1])
              col='lightblue', lwd=lwds[i2], lty=ltys[i1])
      }
    }
  }
}
lines(xx, calc_effacy_through_time(initial_concentration=rtss4_initial_concentration, max_efficacy=rtss4_max_efficacy, fast_frac=rtss4_fast_frac, k1=rtss4_k1, k2=rtss4_k2, nn=rtss4_nn, hh=rtss4_hh, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
      col=rgb(0,0,0,0.2), lwd=3.5, lty=1)


##############################################################################
# booster doses
##############################################################################
# ideally, boosters would increase the concentration and thereby change efficacy. However, if we model efficacy directly, we'll need a different additive or multiplicative approach
# boosters work through increasing concentration
max_efficacy = 0.8
initial_concentration = 1100
fast_frac=0.7
k1=8
k2=100
m2 = -0.1
booster_day = 365

png(filename=paste0(projectpath, '/nonSimFigures/est_mAb_pkpd_wBoosterDay', booster_day, '.png'), width=6, height=2.5, res=900, units='in')
eff_conc_pkpd = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, m2=m2, xx=xx, booster_day=booster_day, create_plot_panel=TRUE)
dev.off()


# from RTSS
xx=1:(365*3)
fast_frac = rtss4_fast_frac
k1=rtss4_k1
k2=rtss4_k2
max_efficacy=rtss4_max_efficacy
initial_concentration=rtss4_initial_concentration
initial_concentration_booster = rtss4_initial_concentration/2
hh=rtss4_hh
nn=rtss4_nn
booster_day=547 # 365, 547
png(filename=paste0(projectpath, '/nonSimFigures/est_rtss_pkpd_wBoosterDay', booster_day, '.png'), width=6, height=2.5, res=900, units='in')
eff_conc_pkpd = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, hh=hh, nn=nn, xx=xx, booster_day=booster_day, hill_func=TRUE, create_plot_panel=TRUE)
dev.off()



# _paramSet2
# max_efficacy = 0.8
# initial_concentration = 1100
# fast_frac=0
# k1=8
# k2=250
# m2 = -.01

# compare alternative approach of just adding efficacies against 'correct' PKPD approach of adding concentrations
png(filename=paste0(projectpath, '/nonSimFigures/compare_conc_eff_booster_approaches_', '_maxEff', round(max_efficacy*100), '.png'), width=2.5*length(booster_days), height=3, res=900, units='in')
par(mfrow=c(1,3))
booster_days = c(182, 365, 365*2)
for(booster_day in booster_days){
  eff_conc_pkpd = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, m2=m2, xx=xx, booster_day=booster_day, create_plot_panel=FALSE)
  # what is an alternative way of doing efficacy through time for use in EMOD (i.e., without using concentrations)?
  # option 1) add new efficacy on top of old
  eff_conc_no_boost = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, m2=m2, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]]
  eff_conc_with_boost = eff_conc_no_boost
  eff_conc_with_boost[booster_day:length(eff_conc_with_boost)] = eff_conc_with_boost[booster_day:length(eff_conc_with_boost)] + eff_conc_no_boost[1:(length(eff_conc_with_boost) - booster_day + 1)]
  eff_conc_with_boost[eff_conc_with_boost>1]=1
  # version removing previous vaccine and applying new one
  eff_conc_with_boost2 = eff_conc_no_boost
  eff_conc_with_boost2[booster_day:length(eff_conc_with_boost)] = eff_conc_no_boost[1:(length(eff_conc_with_boost) - booster_day + 1)]
  
  plot(xx, eff_conc_pkpd[[1]], type='l', col='black', ylim=c(0,1), bty='L', xlab='time', ylab='efficacy', main=paste0('Booster day ', booster_day), lwd=2)
  lines(xx, eff_conc_with_boost, col='blue')
  lines(xx, eff_conc_with_boost2, col='red')
  
}
legend('topright', c('PKPD', 'additive', 'remove old'), lwd=c(2,1,1), col=c('black', 'blue', 'red'), bty='n')
dev.off()
par(mfrow=c(1,1))





# For RTS,S parameters: compare alternative approach of just adding efficacies against 'correct' PKPD approach of adding concentrations
png(filename=paste0(projectpath, '/nonSimFigures/compare_conc_eff_booster_approaches_RTSS_', '_maxEff', round(max_efficacy*100), '.png'), width=2.5*length(booster_days), height=3, res=900, units='in')
par(mfrow=c(1,3))
# from RTSS
xx=1:(365*3)
fast_frac = rtss4_fast_frac
k1=rtss4_k1
k2=rtss4_k2
max_efficacy=rtss4_max_efficacy
initial_concentration=rtss4_initial_concentration
initial_concentration_booster = rtss4_initial_concentration/2
hh=rtss4_hh
nn=rtss4_nn
booster_days = c(182, 365, 365*2)
for(booster_day in booster_days){
  eff_conc_pkpd = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, hh=hh, nn=nn, xx=xx, booster_day=booster_day, create_plot_panel=FALSE)
  # what is an alternative way of doing efficacy through time for use in EMOD (i.e., without using concentrations)?
  # option 1) add new efficacy on top of old
  eff_conc_no_boost = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, hh=hh, nn=nn, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]]
  eff_conc_with_boost = eff_conc_no_boost
  eff_conc_with_boost[booster_day:length(eff_conc_with_boost)] = eff_conc_with_boost[booster_day:length(eff_conc_with_boost)] + eff_conc_no_boost[1:(length(eff_conc_with_boost) - booster_day + 1)]
  eff_conc_with_boost[eff_conc_with_boost>1]=1
  # version removing previous vaccine and applying new one
  eff_conc_with_boost2 = eff_conc_no_boost
  eff_conc_with_boost2[booster_day:length(eff_conc_with_boost)] = eff_conc_no_boost[1:(length(eff_conc_with_boost) - booster_day + 1)]
  
  plot(xx, eff_conc_pkpd[[1]], type='l', col='black', ylim=c(0,1), bty='L', xlab='time', ylab='efficacy', main=paste0('Booster day ', booster_day), lwd=2)
  lines(xx, eff_conc_with_boost, col='blue')
  lines(xx, eff_conc_with_boost2, col='red')
  
}
legend('topright', c('PKPD', 'additive', 'remove old'), lwd=c(2,1,1), col=c('black', 'blue', 'red'), bty='n')
dev.off()
par(mfrow=c(1,1))

##############################################################################
# initial concentration by dose
##############################################################################
dose = c(5, 20, 40)  # mg/kg
mean_concentration = c(198, 934, 1764)  # ug/ml
lower_concentration = c(198-28, 934-293, 1764-260)  # ug/ml
higher_concentration = c(198+28, 934+293, 1764+260)  # ug/ml
convert = mean(mean_concentration/dose)
convert_lower = mean(lower_concentration/dose)
convert_higher = mean(higher_concentration/dose)
dose*convert





##############################################################################
# misc...
##############################################################################

zz=c(seq(0.0001,1,0.01), seq(2,100))
df = data.frame('zz'=zz, 'yy'=2000-zz, 'yy2'=2000-zz*2, 'yy3'=(1-exp(-0.1*zz)))
ggplot(df, aes(x=zz, y=yy3))+geom_line()+scale_x_continuous(trans='log10')



# rtss anti-CSP antibody titres against efficacy
protected=c(2+2+1+5,3+1+9+13, 1+10+6+17)
total_num=c(16+6+6+28, 8+7+14+29, 4+16+8+28)
protected/total_num
plot(c((3+78)/2, (78+183)/2,(183+1136)/2),protected/total_num)
