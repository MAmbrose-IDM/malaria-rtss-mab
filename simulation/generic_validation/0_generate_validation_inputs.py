import os
import sys
import numpy as np
import pandas as pd
import itertools

sys.path.append('../')
from simulation.load_paths import load_box_paths

datapath, projectpath = load_box_paths()
base_scenario_filepath = ''  # os.path.join('scenario_files', 'generic')

# set which validation scenario to generate files for
validation_scenario = 'Phase3_noBooster'  # 'SweepEIR_wBooster', 'SweepEIR_noBooster', 'KintampoPhase3_wBooster', KintampoPhase3_noBooster', 'ChandramohanTrial'

# SMC parameters
num_smc_rounds = 4
smc_start_month = 7
num_repeated_years = 0
# vacc parameters
initial_conc = 620
initial_fast_frac = 0.88
initial_k1 = 46
initial_k2 = 583
initial_max_efficacy = 0.8
initial_hhs = [20, 40, 60]
initial_nns = [1.4, 4]
booster_conc = 620
booster_fast_frac = 0.88
booster_k1 = 46
booster_k2 = 583
booster_max_efficacy = 0.8
booster_hhs = [20, 40, 60]
booster_nns = [1.4, 4]
vacc_total_time = 365*3

agemin_vacc = 0.25
agemax_vacc = 5
vacc_target_group = ['random']
cm_target_group = ['random']
smc_target_group = ['random']

# scenario-specific parameters
if 'SweepEIR' in validation_scenario:
    annual_EIR = [0.23, 0.91, 2.06, 3.61, 5.75, 9.3, 16.54, 26.72, 36.68]  # np.linspace(5,100,4)
    cm_coverage = [0.45]  # np.linspace(0.2, 1, 5)[:-1]
    seasonality = ['constant']
    intervention_correlation = [0]
    vacc_coverage = [0, 1]
    smc_coverage = [0]
    #
    # # vacc arguments
    # booster_coverage_vacc = 1 if 'wBooster' in validation_scenario else 0
    # vacc_filename_description = 'SweepEIR'
    # vacc_3rd_dose_age = 274
    # vacc_coverage_multipliers = [1, booster_coverage_vacc]
    # vacc_rounds = [3, 4]
    # vacc_types = ['simple', 'booster']
    # vacc_initial_killing = [initial_effect_vacc, booster_effect_vacc]
    # vacc_days = [vacc_3rd_dose_age, 810]

elif 'Phase3' in validation_scenario:
    # TODO: iterate through phase 3 sites
    phase3_eirs = pd.read_csv(os.path.join(datapath, 'normalizedEIRprofilesFromRTSSincidenceChildren.csv'))
    site_name = 'Kintampo'
    annual_EIR = [11.1]
    cm_coverage = [0.9]
    seasonality = [site_name]
    intervention_correlation = [0]
    vacc_coverage = [0, 1]
    smc_coverage = [0]

    # vacc arguments
    booster_coverage_vacc = 1 if 'wBooster' in validation_scenario else 0
    vacc_filename_description = 'Phase3_%s' % site_name
    vacc_3rd_dose_age = 365
    vacc_coverage_multipliers = [1, booster_coverage_vacc]
    vacc_rounds = [3, 4]
    # vacc_types = ['simple', 'booster']
    # vacc_initial_killing = [initial_effect_vacc, booster_effect_vacc]
    vacc_days = [vacc_3rd_dose_age, 912]

elif 'ChandramohanTrial' in validation_scenario:
    annual_EIR = [6, 10, 15]  # [10, 15]#, 20]
    cm_coverage = [0.8, 0.9]
    seasonality = ['higher_JuneStart']#, 'aveNanoroKintampo_JuneStart']
    intervention_correlation = [0]
    vacc_coverage = [0, 0.934]
    smc_coverage = [0, 0.89]

    # smc arguments
    smc_start_month = 14.8
    num_repeated_years = 2

    # vacc arguments
    booster_coverage_vacc = 0.95
    vacc_filename_description = 'ChandramohanTrial'
    # vacc_3rd_dose_age = 365  # At age 1 in June (simulation starts in June)
    # vacc_coverage_multipliers = [1, booster_coverage_vacc, booster_coverage_vacc]
    # vacc_rounds = [3, 4, 5]
    # vacc_types = ['simple', 'booster', 'booster']
    # vacc_initial_killing = [initial_effect_vacc, booster_effect_vacc, booster_effect_vacc]
    # vacc_days = [vacc_3rd_dose_age, 365*2, 365*3]

else:
    print('validation scenario not recognized')
    annual_EIR = [0]
    cm_coverage = [0]
    seasonality = ['constant']
    intervention_correlation = [0]
    vacc_coverage = [0]
    smc_coverage = [0]
    # vacc arguments
    booster_coverage_vacc = 0
    vacc_filename_description = 'ERROR'
    vacc_3rd_dose_age = 365
    # vacc_coverage_multipliers = [1]
    # vacc_rounds = [3]
    # vacc_types = ['simple']
    # vacc_initial_killing = [initial_effect_vacc]
    # vacc_days = [vacc_3rd_dose_age]






########################################
# functions
########################################

def remove_duplicate_scenarios(scen_csv):
    # remove rows that are effectively duplicates: these are created from taking the full factorial of parameter sets -
    #    when an intervention coverage is zero, it is not necessary to have separate simulations for different access
    #    groups or schedule types of that intervention
    scen_df = pd.read_csv(os.path.join(projectpath, 'simulation_inputs', scen_csv),  encoding='latin')
    # remove duplicate RTS,S rows when RTS,S coverage is zero
    vacc0_subset_cols = ['annual_EIR', 'seasonality', 'cm_coverage', 'ipti_coverage', 'smc_coverage',
                         'frac_high_access', 'cm_target_group', 'smc_target_group']  # columns to identify duplicates
    vacc0_subset_cols = list(set(vacc0_subset_cols).intersection(scen_df.columns))
    scen_df_vacc0 = scen_df[scen_df.vacc_coverage == 0]
    scen_df_vacc_remainder = scen_df[scen_df.vacc_coverage != 0]
    scen_df_vacc0_fixed = scen_df_vacc0[~scen_df_vacc0.duplicated(subset=vacc0_subset_cols, keep='first')]
    scen_df = pd.concat([scen_df_vacc0_fixed, scen_df_vacc_remainder])

    # remove duplicate SMC rows when SMC coverage is zero
    smc0_subset_cols = ['annual_EIR', 'seasonality', 'cm_coverage', 'ipti_coverage',
                         'vacc_coverage', 'vacc_characteristics', 'vacc_campaign', 'frac_high_access', 'vacc_mode', 'minBoostAge', 'vacc_booster1_min_age',
                         'cm_target_group', 'vacc_target_group']  # columns to identify duplicates
    smc0_subset_cols = list(set(smc0_subset_cols).intersection(scen_df.columns))
    scen_df_smc0 = scen_df[scen_df.smc_coverage == 0]
    scen_df_smc_remainder = scen_df[scen_df.smc_coverage != 0]
    scen_df_smc0_fixed = scen_df_smc0[~scen_df_smc0.duplicated(subset=smc0_subset_cols, keep='first')]
    scen_df = pd.concat([scen_df_smc0_fixed, scen_df_smc_remainder])

    # re-sort rows
    scen_df = scen_df.sort_values(by=['annual_EIR', 'seasonality', 'cm_coverage', 'smc_coverage', 'vacc_coverage', 'vacc_target_group'])

    # re-label settings
    scen_df['setting_id'] = ['HX' + str(ii) for ii in list(range(scen_df.shape[0]))]
    scen_df.set_index(keys='setting_id', inplace=True)

    scen_df.to_csv(os.path.join(projectpath, 'simulation_inputs',
                                scen_csv.replace('.csv', '_cleaned.csv')), encoding='latin')



########################################
# create intervention input files
########################################
# CM
for cm in cm_coverage:
    if cm == 0:
        df = pd.DataFrame()
    else:
        df = pd.DataFrame({'U5_coverage': [cm], 'adult_coverage': [cm], 'severe_cases': np.min([0.8, cm*2]),
                          'simday': [0], 'duration': [-1], 'start_day_override': [-1], 'run_col': 'run'})
    df.to_csv(os.path.join(projectpath, 'simulation_inputs', 'CM',
                           'CM_constant_%icoverage.csv' % (100 * float(cm))))

# RTS,S - two files: one campaign and one vaccine characteristics
for vacc in vacc_coverage:
    if vacc == 0:
        df = pd.DataFrame()
        vacc_char_files = []
    else:
        vacc_coverages = vacc_coverage_multipliers
        vacc_coverages[0] = vacc
        num_vacc = len(vacc_coverages)
        vacc_camp_df = pd.DataFrame({
            'DS_Name': 'all',
            'coverage_levels': vacc_coverages,
            'simday': vacc_days,
            'deploy_type': ['EPI_cohort'] * num_vacc,
            'agemin': [agemin_vacc] * num_vacc,
            'agemax': [agemax_vacc] * num_vacc,
            'runcol': 'run'
        })
        vacc_camp_df.to_csv(os.path.join(projectpath, 'simulation_inputs', 'vaccines',
                               'RTSS_campaign_%s_%icoverage_%ibooster.csv' % (
                               vacc_filename_description, 100 * float(vacc), 100 * float(booster_coverage_vacc))))
        vacc_char_files = []
        for ii in range(len(initial_hhs)):
            for jj in range(len(initial_nns)):
                cur_filename = 'vaccines/RTSS_characteristics_%s_hh%i_nn%i.csv' % (vacc_filename_description, initial_hhs[ii],
                                                                          round(100*initial_nns[jj]))
                vacc_char_df = pd.DataFrame({
                    'vacc_type': ['initial', 'booster'],
                    'vacc_waning_type': ['pkpd']*2,
                    'initial_concentration': [initial_conc, booster_conc],
                    'max_efficacy': [initial_max_efficacy, booster_max_efficacy],
                    'fast_frac': [initial_fast_frac, booster_fast_frac],
                    'k1': [initial_k1, booster_k1],
                    'k2': [initial_k2, booster_k2],
                    'hh': [initial_hhs[ii], booster_hhs[ii]],
                    'nn': [initial_nns[jj], booster_nns[jj]],
                    'total_time': [vacc_total_time]*2,
                })
                vacc_char_df.to_csv(os.path.join(projectpath, 'simulation_inputs', cur_filename))
                vacc_char_files = vacc_char_files + [cur_filename]


# SMC - assumes the same coverage in all rounds
for smc in smc_coverage:
    if smc == 0:
        df_all_years = pd.DataFrame()
    else:
        df_r1 = pd.DataFrame({'round': [1], 'coverage': [smc], 'max_age': [5],
                          'simday': [(smc_start_month-1)*30], 'duration': [-1], 'run_col': 'run'})
        df_base_year = df_r1.copy()
        for rr in range(num_smc_rounds-1):
            df_new_round = df_r1.copy()
            df_new_round['round'] = rr+2
            df_new_round['simday'] = df_r1['simday'] + (1+rr)*30
            df_base_year = df_base_year.append(df_new_round)

        # copy data frame for all of the included simulation years
        df_all_years = df_base_year.copy()
        for yy in range(num_repeated_years):
            df_new_year = df_base_year.copy()
            df_new_year['simday'] = df_base_year['simday'] + 365 * (yy+1)
            df_all_years = df_all_years.append(df_new_year)
    df_all_years.to_csv(os.path.join(projectpath, 'simulation_inputs', 'SMC',
                                     'SMC_%s_%icoverage.csv' % (vacc_filename_description, 100 * float(smc))))



########################################
# create scenario description csv
########################################

# Create full factorial
df_array = np.array(list(itertools.product(annual_EIR, cm_coverage, seasonality, intervention_correlation, vacc_coverage, smc_coverage, vacc_char_files, vacc_target_group, smc_target_group, cm_target_group)))
df = pd.DataFrame(df_array)

# Optionally drop certain combinations
# [...]

# Prepare and save df
df.columns = ['annual_EIR', 'cm_coverage', 'seasonality', 'intervention_correlation', 'vacc_coverage', 'smc_coverage', 'vacc_characteristics', 'vacc_target_group', 'smc_target_group', 'cm_target_group']
df.reset_index(inplace=True)
df = df.rename(columns={'index': 'setting_id'})
df['setting_id'] = 'HX' + df['setting_id'].astype(str)

# for the correlated demographics file, use the minimum of the coverages as the fraction of people in the high-access group
df['frac_high_access'] = [np.min([float(df['cm_coverage'][yy]),
                                  (float(df['vacc_coverage'][yy]) if float(df['vacc_coverage'][yy])>0 else 1),
                                  (float(df['smc_coverage'][yy]) if float(df['smc_coverage'][yy])>0 else 1)]) for yy in range(len(df))]
df.loc[df.intervention_correlation == '0', 'frac_high_access'] = 0

# add demographics, seasonality, and intervention csv filepaths
# df['demographics_filename'] = os.path.join(base_scenario_filepath, 'Demographics', 'generic_demographics_cohort_uncorrelated.json')
# df.loc[df.intervention_correlation == 1, 'demographics_filename'] = [os.path.join(base_scenario_filepath, 'Demographics', 'generic_demographics_cohort_correlated%i.json' % int(float(df['frac_high_access'][yy])*100)) for yy in range(len(df))]
df['demographics_filename'] = [os.path.join(base_scenario_filepath, 'Demographics', 'generic_demographics_cohort_uncorrelated.json') if df['intervention_correlation'][yy] == '0' else
                               os.path.join(base_scenario_filepath, 'Demographics', 'generic_demographics_cohort_correlated%i.json' % int(float(df['frac_high_access'][yy])*100)) for yy in range(len(df))]
# df['seasonality_filename'] = os.path.join(base_scenario_filepath, 'Seasonality', 'seasonality_eir_multipliers')
df['CM_filename'] = [os.path.join(base_scenario_filepath, 'CM', 'CM_constant_%icoverage.csv' % (100 * float(yy))) for yy in df['cm_coverage']]
df['SMC_filename'] = [os.path.join(base_scenario_filepath, 'SMC', 'SMC_%s_%icoverage.csv' % (vacc_filename_description, 100 * float(yy))) for yy in df['smc_coverage']]
df.loc[df.smc_coverage == '0', 'SMC_filename'] = ''
df['VACC_filename'] = [os.path.join(base_scenario_filepath, 'vaccines', 'RTSS_campaign_%s_%icoverage_%ibooster.csv' % (vacc_filename_description, 100 * float(yy), 100 * float(booster_coverage_vacc))) for yy in df['vacc_coverage']]

scen_csv = 'coordinator_files/validation_%s.csv' % validation_scenario
df.to_csv(os.path.join(projectpath, 'simulation_inputs', scen_csv), index=False)
remove_duplicate_scenarios(scen_csv)
