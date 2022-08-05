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



######################################################################
# function using PKPD to get concentration and efficacy through time
######################################################################
calc_concentration_through_time = function(initial_concentration, fast_frac, k1, k2, xx=seq(1,365*3)){
  # concentration_through_time = c(initial_concentration*exp(-xx[xx<=time_switch]/k1), (initial_concentration*exp(-1*time_switch/k1))*exp(-(xx[xx>time_switch]-time_switch)/k2))
  concentration_through_time = initial_concentration*(fast_frac*exp(-xx / k1) + (1-fast_frac)*exp(-xx / k2))
  return(concentration_through_time)
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
  if(hill_func){
    efficacy_through_time = max_efficacy / (1+(hh/concentration_through_time)^nn)
  } else{
    efficacy_through_time = max_efficacy * (1-exp(m2*concentration_through_time))
  }

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





###################################################################################################################################################
# RTS,S protective efficacy through time from current assumptions versus pkpd functions to closely match current shape
###################################################################################################################################################
# time plotted
xx=1:(365*3)

png(filename=paste0(projectpath, '/nonSimFigures/emod_function_rtss.png'), width=5.5, height=4.5, res=900, units='in')
plot(NA, ylim=c(0,1), xlim=c(0,max(xx)), ylab='protective efficacy', xlab='days since dose', bty='L')

# values described in RTSS report
init=0.8; dur=410; color='grey'
lines(init*exp(xx*-1/dur), ylab='protective efficacy', xlab='days since dose', bty='L', lwd=2, col=color)

# values used in RTSS analysis
init=0.8; dur=592; color='black'
lines(init*exp(xx*-1/dur), ylab='protective efficacy', xlab='days since dose', bty='L', lwd=2, col=color)
# 
# # imitate with PKPD
# initial_concentration = rtss1_initial_concentration
# max_efficacy = rtss1_max_efficacy
# fast_frac=rtss1_fast_frac
# k1=rtss1_k1
# k2=rtss1_k2
# hh=rtss1_hh
# nn=rtss1_nn
# eff_conc = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, hh=hh, nn=nn, hill_func=TRUE, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]]
# lines(eff_conc, lwd=2, col='green')

# imitate with PKPD
initial_concentration = rtss3_initial_concentration
max_efficacy = rtss3_max_efficacy
fast_frac=rtss3_fast_frac
k1=rtss3_k1
k2=rtss3_k2
# hh=rtss3_hh
# nn=rtss3_nn
# eff_conc = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, hh=hh, nn=nn, hill_func=TRUE, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]]
# lines(eff_conc, lwd=2, col='green')
# # a couple other parameter sets
# initial_concentration = rtss3_initial_concentration
# max_efficacy = rtss3_max_efficacy
# fast_frac=rtss3_fast_frac
# k1=rtss3_k1
# k2=rtss3_k2
# hh=20
# nn=rtss1_nn
# eff_conc = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, hh=hh, nn=nn, hill_func=TRUE, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]]
# lines(eff_conc, lwd=2, col='blue')
hh=40
nn=2
eff_conc = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, hh=hh, nn=nn, hill_func=TRUE, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]]
lines(eff_conc, lwd=2, col='salmon')
# legend('topright', c('RTSS analysis (initial=0.8, duration=592)', 'RTSS report (initial=0.8, duration=410)', 'approximate match with PKPD', 'other PKPD params (v1)', 'other PKPD params (v2)', 'selected from Phase 3 comparisons'), 
#        col=c('black', 'grey', 'green', 'blue', 'dodgerblue', 'salmon'), lwd=2, bty='n')
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
# compare exponential and biphasic decay against mAb concentration data
##############################################################################

make_concentration_plot = function(df_mg, fast_frac, k1, k2, m1, title='', conc_col='serum_concentration_ug_ml', time_col='time_weeks', max_days=7*24, include_exp=TRUE){
  initial_concentration = df_mg[[conc_col]][1]
  if(grepl('week', time_col)){
    time_multiplier = 7
  } else if(grepl('year', time_col)){
    time_multiplier = 365
  } else{
    time_multiplier = 1
  }
  plot(df_mg[[time_col]]*time_multiplier, df_mg[[conc_col]], ylim=c(1,1800), pch=20, log='y', xlim=c(0,max_days), bty='L', xlab='time (days)', ylab='concentration', main=title)
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


fast_frac=mab_fast_frac
k1=mab_k1
k2=mab_k2
m1 = 50


ref_csvs = paste0(datapath, '/vacc_pkpd/', c('rough_dataGrab_5mgSC_Gaudinski_SF1.csv', 'rough_dataGrab_5mg_Gaudinski_SF1.csv', 'rough_dataGrab_20mg_Gaudinski_SF1.csv', 'rough_dataGrab_40mg_Gaudinski_SF1.csv'))
ref_names = c('5mg/kg (SC)', '5mg/kg (IV)', '20mg/kg (IV)', '40mg/kg (IV)')
ref_weights = c(3, 4, 5, 6)
# compare roughly data-extracted values with simple exponential versus biphasic decay
png(filename=paste0(projectpath, '/nonSimFigures/mAb_concDecay_wData.png'), width=6, height=6, res=900, units='in')
par(mfrow=c(2,2))
for(rr in 1:length(ref_csvs)){
  df_mg = read.csv(ref_csvs[rr])
  make_concentration_plot(df_mg, fast_frac, k1, k2, m1, title=ref_names[rr])
  if(rr==1){
    legend('topleft', c('data','exponential', 'biphasic'), pch=c(20,NA,NA), lwd=c(NA,1,1), col=c('black','blue','red'), bty='n')
  }
}
par(mfrow=c(1,1))
dev.off()


# fit best pkpd parameters to mAb data
# since I don't have actual data, just plot-grabbed values for means, use distance between model and mean datapoint, weighted by number of individuals
get_ref_func_distance = function(df_mg, conc_col='serum_concentration_ug_ml', time_col='time_weeks', fast_frac, k1, k2, get_rel_diff=TRUE){
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
# # iteration 1
# fast_frac_vals = seq(0,1,length.out=15)
# k1_vals = 10^seq(-1, 3, length.out=20)
# k2_vals = 10^seq(-1, 3, length.out=20)
## -> results: best_match = c(5, 15, 10), corresponding to fast_frac=0.29, k1=88.59, k2=7.85. min(diff_outputs) = 5.03
# # iteration 2
# fast_frac_vals = seq(0.5,1,length.out=10)
# k1_vals = 10^seq(-1, 3, length.out=25)
# k2_vals = 10^seq(-1, 3, length.out=25)
## -> results: best_match = c(5, 12, 19), corresponding to fast_frac=0.722, k1=6.81, k2=100. min(diff_outputs) = 4.94
# iteration 3
fast_frac_vals = seq(0.5,1,length.out=10)
k1_vals = seq(4, 10, length.out=25)
k2_vals = seq(60, 130, length.out=25)
## -> results: best_match = c(5,11,13), corresponding to fast_frac=0.722, k1=6.5, k2=95. min(diff_outputs) = 4.9
diff_outputs = array(0, dim=c(length(fast_frac_vals), length(k1_vals), length(k2_vals)))
for(i1 in 1:length(fast_frac_vals)){
  for(i2 in 1:length(k1_vals)){
    for(i3 in 1:length(k2_vals)){
      for(rr in 1:length(ref_csvs)){
        df_mg = read.csv(ref_csvs[rr])
        cur_diff = get_ref_func_distance(df_mg, conc_col='serum_concentration_ug_ml', time_col='time_weeks', fast_frac=fast_frac_vals[i1], k1=k1_vals[i2], k2=k2_vals[i3], get_rel_diff=TRUE)
        diff_outputs[i1,i2,i3] = diff_outputs[i1,i2,i3] + ref_weights[rr] * cur_diff
      }
    }
  }
}



best_match = which(diff_outputs == min(diff_outputs), arr.ind=TRUE)[1,]
png(filename=paste0(projectpath, '/nonSimFigures/mAb_concDecay_wData_fit.png'), width=6, height=5, res=900, units='in')
par(mfrow=c(2,2))
for(rr in 1:length(ref_csvs)){
  df_mg = read.csv(ref_csvs[rr])
  make_concentration_plot(df_mg, fast_frac=fast_frac_vals[best_match[1]], k1=k1_vals[best_match[2]], k2=k2_vals[best_match[3]], m1=NA, title=ref_names[rr],include_exp=FALSE)
}
dev.off()
par(mfrow=c(1,1))
print(best_match)
print(min(diff_outputs))
print(fast_frac_vals[best_match[1]])
print(k1_vals[best_match[2]])
print(k2_vals[best_match[3]])



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
