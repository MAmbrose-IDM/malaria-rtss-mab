import os
import numpy as np
import pandas as pd
import itertools


########################################
# functions to generate sim input files
########################################

def get_parameter_space():
    """ get dictionary of default parameter values that will be updated for specific scenarios
    """
    parameter_space = {
        # ---  parameters that are explored in a full factorial setup
        # site specifications
        'annual_EIR': [1, 10, 30, 60],
        'seasonality': ['constant', 'moderate_unimodal', 'high_unimodal'],
        # CM, SMC, and IPTi
        'cm_coverage': [0.6, 0.8],
        'smc_coverage': [0, 0.6, 0.8],
        'ipti_coverage': [0],
        'ipti_mode': ['basic'],  # options are: 'basic', 'extended3tp', 'extended'
        # Vaccine
        'vacc_coverage': [0, 0.6, 0.8],
        'vacc_deploy_type': ['EPI_cohort'],  # options are: 'EPI_cohort' (EPI main dose and booster dose),
        #              'seasonal' (all doses seasonally)
        #              'campboost' (EPI main dose and 1 campaign booster dose),
        #              'campboost2' (EPI main dose and 2 campaign booster doses)
        'initial_hhs': [30, 40],  # [20, 40, 60],
        'initial_nns': [1.4, 2, 4],  # [1.4, 4],
        'booster_hhs': [30, 40],  # [20, 40, 60],
        'booster_nns': [1.4, 2, 4],  # [1.4, 4],
        # for RTS,S, it seems that the best parameters are hh=40, nn=2
        # Intervention targeting
        'vacc_target_group': ['random'],
        # specify whether vaccine is targeted toward an access group or given at random
        'smc_target_group': ['random'],  # specify whether SMC is targeted toward an access group or given at random
        'cm_target_group': ['random'],  # specify whether CM is targeted toward an access group or given at random

        # --- additional parameters that apply to all simulations
        # CM, SMC, and IPTi
        'ipti_touchpoints': [61, 91, 274],
        'ipti_postprocess': True,  # if true IPTi wont be simulated but scenario csv used for postprocessing
        'num_smc_rounds': 4,
        'smc_start_month': 7,
        'num_repeated_years': 0,  # number of repeated years for SMC
        # Vaccine
        'agemin_vacc': 0.25,
        'agemax_vacc': 5,
        'vacc_coverage_multipliers': [1, 0],  # coverage of initial and booster doses
        'vacc_filename_description': 'default',
        'vacc_days': [365, 365*2],
        'vacc_total_time': 365 * 3,
        'initial_conc': 620,
        'initial_fast_frac': 0.88,
        'initial_k1': 46,
        'initial_k2': 583,
        'initial_max_efficacy': 0.8,
        'booster_conc': 620,
        'booster_fast_frac': 0.88,
        'booster_k1': 46,
        'booster_k2': 583,
        'booster_max_efficacy': 0.6,  # 0.8
    }
    return parameter_space


def create_intervention_inputs(param_dic, projectpath):
    """ create input csvs with the relevant intervention specifications for the scenarios"""
    cm_coverage = param_dic['cm_coverage']
    smc_coverage = param_dic['smc_coverage']
    vacc_coverage = param_dic['vacc_coverage']
    vacc_coverage_multipliers = param_dic['vacc_coverage_multipliers']
    vacc_days = param_dic['vacc_days']
    vacc_deploy_type = param_dic['vacc_deploy_type']
    agemin_vacc = param_dic['agemin_vacc']
    agemax_vacc = param_dic['agemax_vacc']
    vacc_filename_description = param_dic['vacc_filename_description']
    initial_hhs = param_dic['initial_hhs']
    booster_hhs = param_dic['booster_hhs']
    initial_nns = param_dic['initial_nns']
    booster_nns = param_dic['booster_nns']
    initial_conc = param_dic['initial_conc']
    booster_conc = param_dic['booster_conc']
    initial_max_efficacy = param_dic['initial_max_efficacy']
    booster_max_efficacy = param_dic['booster_max_efficacy']
    initial_fast_frac = param_dic['initial_fast_frac']
    booster_fast_frac = param_dic['booster_fast_frac']
    initial_k1 = param_dic['initial_k1']
    booster_k1 = param_dic['booster_k1']
    initial_k2 = param_dic['initial_k2']
    booster_k2 = param_dic['booster_k2']
    vacc_total_time = param_dic['vacc_total_time']
    smc_start_month = param_dic['smc_start_month']
    num_smc_rounds = param_dic['num_smc_rounds']
    num_repeated_years = param_dic['num_repeated_years']

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
            vacc_char_files = []
        else:
            vacc_coverages = vacc_coverage_multipliers
            booster_coverage_vacc = vacc_coverages[1]
            vacc_coverages[0] = vacc
            num_vacc = len(vacc_coverages)
            vacc_camp_df = pd.DataFrame({
                'DS_Name': ['all'] * num_vacc,
                'coverage_levels': vacc_coverages,
                'vacc_day': vacc_days,
                'deploy_type': vacc_deploy_type * num_vacc,
                'agemin': [agemin_vacc] * num_vacc,
                'agemax': [agemax_vacc] * num_vacc,
                'runcol': ['run'] * num_vacc
            })
            vacc_camp_df.to_csv(os.path.join(projectpath, 'simulation_inputs', 'vaccines',
                                   'RTSS_campaign_%s_%icoverage_%ibooster.csv' % (
                                   vacc_filename_description, 100 * float(vacc), 100 * float(booster_coverage_vacc))))
            vacc_char_files = []
            for ii in range(len(initial_hhs)):
                for jj in range(len(initial_nns)):
                    cur_filename = 'vaccines/RTSS_characteristics_%s_hh%i_nn%i_bb%i.csv' % (vacc_filename_description, initial_hhs[ii],
                                                                              round(100*initial_nns[jj]), round(100*booster_max_efficacy))
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
                                  'simday': [(smc_start_month-1)*30], 'duration': [-1], 'run_col': ['run']})
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
    return vacc_char_files


def create_coordinator_csvs(param_dic, base_scenario_filepath, vacc_char_files=[], high_access_frac=None):
    annual_EIR = param_dic['annual_EIR']
    seasonality = param_dic['seasonality']
    cm_coverage = param_dic['cm_coverage']
    smc_coverage = param_dic['smc_coverage']
    vacc_coverage = param_dic['vacc_coverage']
    vacc_filename_description = param_dic['vacc_filename_description']
    cm_target_group = param_dic['cm_target_group']
    smc_target_group = param_dic['smc_target_group']
    vacc_target_group = param_dic['vacc_target_group']
    booster_coverage_vacc = param_dic['vacc_coverage_multipliers'][1]

    # Create full factorial
    df_array = np.array(list(
        itertools.product(annual_EIR, cm_coverage, seasonality, vacc_coverage, smc_coverage,
                          vacc_char_files, vacc_target_group, smc_target_group, cm_target_group)))
    df = pd.DataFrame(df_array)

    # Prepare and save df
    df.columns = ['annual_EIR', 'cm_coverage', 'seasonality', 'vacc_coverage',
                  'smc_coverage', 'vacc_characteristics', 'vacc_target_group', 'smc_target_group', 'cm_target_group']
    df.reset_index(inplace=True)
    df = df.rename(columns={'index': 'setting_id'})
    df['setting_id'] = 'HX' + df['setting_id'].astype(str)

    # for the correlated demographics file, use the minimum of the coverages as the fraction of people in the high-access group

    # for the correlated demographics file, use the minimum of the coverages as the fraction of people in the high-access group unless otherwise specified
    if type(high_access_frac) == int or type(high_access_frac) == float:
        df['frac_high_access'] = high_access_frac
    else:
        df['frac_high_access'] = [np.min([float(df['cm_coverage'][yy]),
                                          (float(df['vacc_coverage'][yy]) if float(df['vacc_coverage'][yy]) > 0 else 1),
                                          (float(df['smc_coverage'][yy]) if float(df['smc_coverage'][yy]) > 0 else 1)])
                                  for yy in range(len(df))]

    # add demographics, seasonality, and intervention csv filepaths
    df['demographics_filename'] = [os.path.join(base_scenario_filepath, 'Demographics',
                                   'generic_demographics_cohort_correlated%i.json' % int(float(df['frac_high_access'][yy]) * 100))
                                   for yy in range(len(df))]
    df.loc[df.frac_high_access == 0, 'demographics_filename'] = os.path.join(base_scenario_filepath, 'Demographics',
                                   'generic_demographics_cohort_uncorrelated.json')

    df['CM_filename'] = [os.path.join(base_scenario_filepath, 'CM', 'CM_constant_%icoverage.csv' % (100 * float(yy)))
                         for yy in df['cm_coverage']]
    df['SMC_filename'] = [os.path.join(base_scenario_filepath, 'SMC',
                                       'SMC_%s_%icoverage.csv' % (vacc_filename_description, 100 * float(yy))) for yy in
                          df['smc_coverage']]
    df.loc[df.smc_coverage == '0', 'SMC_filename'] = ''
    df['VACC_filename'] = [os.path.join(base_scenario_filepath, 'vaccines',
                                        'RTSS_campaign_%s_%icoverage_%ibooster.csv' % (
                                        vacc_filename_description, 100 * float(yy), 100 * float(booster_coverage_vacc)))
                           for yy in df['vacc_coverage']]
    df.loc[df.vacc_coverage == '0', 'VACC_filename'] = ''

    return df


def combine_input_files(scen_df1, scen_df2):
    scen_df = pd.concat([scen_df1, scen_df2], sort=False).reset_index(drop=True)
    # re-label settings
    scen_df['setting_id'] = ['HX' + str(ii) for ii in list(range(scen_df.shape[0]))]
    scen_df.set_index(keys='setting_id', inplace=True)
    return scen_df


def combine_input_files_from_csv(scen_csv1, scen_csv2, combined_scen_csv, projectpath):
    scen_df1 = pd.read_csv(os.path.join(projectpath, 'simulation_inputs', scen_csv1),
                          encoding='latin')
    scen_df2 = pd.read_csv(os.path.join(projectpath, 'simulation_inputs', scen_csv2),
                          encoding='latin')
    scen_df = pd.concat([scen_df1, scen_df2]).reset_index(drop=True)

    # re-label settings
    scen_df['setting_id'] = ['HX' + str(ii) for ii in list(range(scen_df.shape[0]))]
    scen_df.set_index(keys='setting_id', inplace=True)

    scen_df.to_csv(os.path.join(projectpath, 'simulation_inputs', combined_scen_csv), encoding='latin')



def remove_duplicate_scenarios(scen_csv, projectpath):
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
