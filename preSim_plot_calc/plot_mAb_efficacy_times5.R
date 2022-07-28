# plot_mAb_efficacy_time

library(ggplot2)

setwd("C:/Users/moniqueam/Documents/malaria-rtss-mab")
source(file.path('simulation', 'load_paths.R'))

paths = get_project_paths()
datapath = paths[1]
projectpath = paths[2]

##############################################################################
# current parameter eyeballed approximate guestimates
##############################################################################
# monoclonals
mab_initial_concentration = 1100
mab_max_efficacy = 0.95
mab_fast_frac=0.7
mab_k1=8
mab_k2=100
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
rtss3_initial_concentration = 620
rtss3_hh = 40
rtss3_nn = 1.4



######################################################################
# function using PKPD to get concentration and efficacy through time
######################################################################
calc_concentration_through_time = function(initial_concentration, max_efficacy, fast_frac, k1, k2, xx=seq(1,365*3)){
  # concentration_through_time = c(initial_concentration*exp(-xx[xx<=time_switch]/k1), (initial_concentration*exp(-1*time_switch/k1))*exp(-(xx[xx>time_switch]-time_switch)/k2))
  concentration_through_time = initial_concentration*(fast_frac*exp(-xx / k1) + (1-fast_frac)*exp(-xx / k2))
  return(concentration_through_time)
}

calc_effacy_through_time = function(initial_concentration, max_efficacy, fast_frac, k1, k2, m2=NA, hh=NA, nn=NA, hill_func=TRUE, xx=seq(1,365*3), booster_day=NA, create_plot_panel=FALSE){  # time_switch
  concentration_through_time = calc_concentration_through_time(initial_concentration, max_efficacy, fast_frac, k1, k2, xx)

  if(!is.na(booster_day)){
    after_booster_initial_concentration = concentration_through_time[booster_day] + initial_concentration
    bb = 1:(max(xx)-booster_day)
    concentration_through_time = c(concentration_through_time[1:booster_day], 
                                   after_booster_initial_concentration*(fast_frac*exp(-bb / k1) + (1-fast_frac)*exp(-bb / k2)))
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
hh=rtss3_hh
nn=rtss3_nn
eff_conc = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, hh=hh, nn=nn, hill_func=TRUE, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]]
lines(eff_conc, lwd=2, col='green')

# a couple other parameter sets
initial_concentration = rtss3_initial_concentration
max_efficacy = rtss3_max_efficacy
fast_frac=rtss3_fast_frac
k1=rtss3_k1
k2=rtss3_k2
hh=20
nn=rtss1_nn
eff_conc = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, hh=hh, nn=nn, hill_func=TRUE, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]]
lines(eff_conc, lwd=2, col='blue')
hh=60
nn=4#rtss1_nn
eff_conc = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=k1, k2=k2, hh=hh, nn=nn, hill_func=TRUE, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]]
lines(eff_conc, lwd=2, col='dodgerblue')


legend('topright', c('RTSS analysis (initial=0.8, duration=592)', 'RTSS report (initial=0.8, duration=410)', 'approximate match with PKPD', 'other PKPD params (v1)', 'other PKPD params (v2)'), 
       col=c('black', 'grey', 'green', 'blue', 'dodgerblue'), lwd=2, bty='n')
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
efficacy_time = calc_effacy_through_time(initial_concentration=initial_concentration, max_efficacy=max_efficacy, fast_frac=fast_frac, k1=m1, k2=m1, nn=nn, hh=hh, hill_func=TRUE, xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]]
lines(efficacy_time, col=color, lwd=2)

legend('topright', c('exponential (duration=1095)', 'exponential (duration=592)', 'exponential (duration=410)', 'exponential (duration=200)', 'E(T)'), 
       col=c('orange', 'darkred', 'black', 'grey', 'dodgerblue'), lwd=2, bty='n')
dev.off()



##############################################################################
# compare exponential and biphasic decay against mAb concentration data
##############################################################################

make_concentration_plot = function(df_mg, fast_frac, k1, k2, m1, title='', conc_col='serum_concentration_ug_ml', time_col='time_weeks', max_days=7*24){
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
  # simple exponential
  concentration_through_time = initial_concentration*exp(-xx/m1)
  lines(xx, concentration_through_time, col='blue')
  # current PKPD with biphasic
  concentration_through_time2 =  calc_concentration_through_time(initial_concentration, max_efficacy, fast_frac, k1, k2, xx)
  lines(xx, concentration_through_time2, type='l', ylim=c(0,initial_concentration), bty='L', xlab='time', main='Concentration through time', col='red')
}


fast_frac=mab_fast_frac
k1=mab_k1
k2=mab_k2
m1 = 50

# compare roughly data-extracted values with simple exponential versus biphasic decay
png(filename=paste0(projectpath, '/nonSimFigures/mAb_concDecay_wData.png'), width=6, height=6, res=900, units='in')
par(mfrow=c(2,2))
# 5mg (SC)
df_mg = read.csv(paste0(datapath, '/vacc_pkpd/rough_dataGrab_5mgSC_Gaudinski_SF1.csv'))
make_concentration_plot(df_mg, fast_frac, k1, k2, m1, title='5mg (SC)')
legend('topleft', c('data','exponential', 'biphasic'), pch=c(20,NA,NA), lwd=c(NA,1,1), col=c('black','blue','red'), bty='n')
# 5mg
df_mg = read.csv(paste0(datapath, '/vacc_pkpd/rough_dataGrab_5mg_Gaudinski_SF1.csv'))
make_concentration_plot(df_mg, fast_frac, k1, k2, m1, title='5mg (IV)')
# 20mg
df_mg = read.csv(paste0(datapath, '/vacc_pkpd/rough_dataGrab_20mg_Gaudinski_SF1.csv'))
make_concentration_plot(df_mg, fast_frac, k1, k2, m1, title='20mg (IV)')
# 40mg
df_mg = read.csv(paste0(datapath, '/vacc_pkpd/rough_dataGrab_40mg_Gaudinski_SF1.csv'))
make_concentration_plot(df_mg, fast_frac, k1, k2, m1, title='40mg (IV)')
par(mfrow=c(1,1))
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
png(filename=paste0(projectpath, '/nonSimFigures/RTSS_concDecay_wData.png'), width=6, height=6, res=900, units='in')
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

# RTSS params from Imperial values in Penny et al 2016

fast_frac=rtss3_fast_frac
k1=rtss3_k1
k2=rtss3_k2
m1 = rtss3_m1
# compare roughly data-extracted values with simple exponential versus biphasic decay
png(filename=paste0(projectpath, '/nonSimFigures/RTSS_concDecay_wData_ImperialParams.png'), width=6, height=6, res=900, units='in')
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
m2s = c(-0.1, -1, -10)
colors=c('red', 'blue', 'grey')
lwds = c(1,2)
ltys=c(1,2, 3)
png(filename=paste0(projectpath, '/nonSimFigures/sweep_mAb_OldParams.png'), width=6, height=2.5, res=900, units='in')
plot(NA, xlim=c(0,max(xx)), ylim=c(0,1), bty='L', xlab='time', ylab='efficacy', main='Alternative TPPs')
for(i1 in 1:length(initial_concentrations)){
  for(i2 in 1:length(max_efficacies)){
    for(i3 in 1:length(m2s)){
      lines(xx, calc_effacy_through_time(initial_concentration=initial_concentrations[i1], max_efficacy=max_efficacies[i2], fast_frac=fast_frac, k1=k1, k2=k2, m2=m2s[i3], xx=xx, booster_day=NA, create_plot_panel=FALSE)[[1]],
            col=colors[i1], lwd=lwds[i2], lty=ltys[i3])
    }
  }
}
dev.off()


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
