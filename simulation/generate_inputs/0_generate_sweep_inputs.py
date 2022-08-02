import os
import sys
import pandas as pd

from simulation.generate_inputs.helpers_input_file_generation import get_parameter_space, create_intervention_inputs, \
    create_coordinator_csvs, remove_duplicate_scenarios, combine_input_files
sys.path.append('../')
from simulation.load_paths import load_box_paths
datapath, projectpath = load_box_paths()
base_scenario_filepath = ''  # os.path.join('scenario_files', 'generic')

# set which validation scenario to generate files for
sweep1 = True  # compares RTS,S and mAbs for a larger number of mAb hh values
sweep2 = False  # compares RTS,S and mAbs for a larger number of setting EIRs

# shared default parameters for sweeps
d_seasonal_campaign_month = 6
d_vacc_dates = [round(d_seasonal_campaign_month*30.4+365*xx) for xx in range(0, 10)]
d_vaccine_coverage = 0.8
d_vacc_coverage_multipliers = [1] * len(d_vacc_dates)
d_cm_coverage = 0.6
d_smc_coverage = 0
d_vacc_target_group = ['high']
d_smc_target_group = ['high']
d_cm_target_group = ['random']
d_high_access_frac = 0.6

# RTS,S-specific vaccine parameters
r_initial_hhs = [40]
r_initial_nns = [2]
r_booster_hhs = [40]
r_booster_nns = [2]
r_vacc_filename_description = 'rtss'
r_initial_conc = 620
r_initial_fast_frac = 0.88
r_initial_k1 = 46
r_initial_k2 = 583
r_initial_max_efficacy = [0.8]
r_booster_conc = 620
r_booster_fast_frac = 0.88
r_booster_k1 = 46
r_booster_k2 = 583
r_booster_max_efficacy = [0.8]
r_vacc_total_time = 365 * 3

# mAb-specific vaccine parameters
m_initial_hhs = [40]
m_initial_nns = [2]
m_booster_hhs = [40]
m_booster_nns = [2]
m_vacc_filename_description = 'mab'
m_initial_conc = 1000
m_initial_fast_frac = 0.722
m_initial_k1 = 6.5
m_initial_k2 = 95
m_initial_max_efficacy = [0.95]
m_booster_conc = 1000
m_booster_fast_frac = 0.722
m_booster_k1 = 6.5
m_booster_k2 = 95
m_booster_max_efficacy = [0.95]
m_vacc_total_time = 600


# scenario-specific parameters
if sweep1:
    validation_scenario = 'sweep1'
    param_dic = get_parameter_space()
    # update dictionary for this scenario
    # setting
    param_dic.update({'annual_EIR': [30]})
    param_dic.update({'seasonality': ['constant', 'high_unimodal']})
    # intervention campaigns
    param_dic.update({'vacc_days': d_vacc_dates})
    param_dic.update({'vacc_deploy_type': ['seasonal'] * len(d_vacc_dates)})
    param_dic.update({'vacc_coverage': [0, d_vaccine_coverage]})
    param_dic.update({'vacc_coverage_multipliers': [1] * len(d_vacc_dates)})
    param_dic.update({'cm_coverage': [d_cm_coverage]})
    param_dic.update({'smc_coverage': [0]})
    param_dic.update({'num_repeated_years': 0})  # for SMC
    param_dic.update({'vacc_target_group': d_vacc_target_group})
    param_dic.update({'cm_target_group': d_cm_target_group})
    param_dic.update({'smc_target_group': d_smc_target_group})
    param_dic.update({'high_access_frac': d_high_access_frac})

    # === RTS,S and no-vaccine scenarios === #
    # vaccine pk/pd
    param_dic.update({'vacc_filename_description': r_vacc_filename_description})
    param_dic.update({'initial_hhs': r_initial_hhs})
    param_dic.update({'initial_nns': r_initial_nns})
    param_dic.update({'booster_hhs': r_booster_hhs})
    param_dic.update({'booster_nns': r_booster_nns})
    param_dic.update({'initial_conc': r_initial_conc})
    param_dic.update({'initial_fast_frac': r_initial_fast_frac})
    param_dic.update({'initial_k1': r_initial_k1})
    param_dic.update({'initial_k2': r_initial_k2})
    param_dic.update({'initial_max_efficacy': r_initial_max_efficacy})
    param_dic.update({'booster_conc': r_booster_conc})
    param_dic.update({'booster_fast_frac': r_booster_fast_frac})
    param_dic.update({'booster_k1': r_booster_k1})
    param_dic.update({'booster_k2': r_booster_k2})
    param_dic.update({'booster_max_efficacy': r_booster_max_efficacy})
    param_dic.update({'vacc_total_time': r_vacc_total_time})

    # create input csvs
    vacc_char_files = create_intervention_inputs(param_dic=param_dic, projectpath=projectpath)
    # create coordinator csv
    df_rtss = create_coordinator_csvs(param_dic=param_dic, base_scenario_filepath=base_scenario_filepath,
                                      vacc_char_files=vacc_char_files)

    # === mAb scenarios === #
    # vaccine pk/pd
    param_dic.update({'vacc_filename_description': m_vacc_filename_description})
    param_dic.update({'initial_max_efficacy': [0.8, 0.9, 0.95]})
    param_dic.update({'booster_max_efficacy': [0.8, 0.9, 0.95]})
    param_dic.update({'initial_hhs': [5, 10, 20, 40, 60]})
    param_dic.update({'booster_hhs': [5, 10, 20, 40, 60]})
    param_dic.update({'initial_nns': m_initial_nns})
    param_dic.update({'booster_nns': m_booster_nns})
    param_dic.update({'initial_conc': m_initial_conc})
    param_dic.update({'booster_conc': m_booster_conc})
    param_dic.update({'initial_fast_frac': m_initial_fast_frac})
    param_dic.update({'booster_fast_frac': m_booster_fast_frac})
    param_dic.update({'initial_k1': m_initial_k1})
    param_dic.update({'booster_k1': m_booster_k1})
    param_dic.update({'initial_k2': m_initial_k2})
    param_dic.update({'booster_k2': m_booster_k2})
    param_dic.update({'vacc_total_time': m_vacc_total_time})
    # create input csvs
    vacc_char_files = create_intervention_inputs(param_dic=param_dic, projectpath=projectpath)
    # create coordinator csv
    df_mab = create_coordinator_csvs(param_dic=param_dic, base_scenario_filepath=base_scenario_filepath,
                                 vacc_char_files=vacc_char_files)

    # === combine and save dataframes === #
    df = combine_input_files(scen_df1=df_rtss, scen_df2=df_mab)
    scen_csv = 'coordinator_files/%s.csv' % validation_scenario
    df.to_csv(os.path.join(projectpath, 'simulation_inputs', scen_csv), index=False)
    remove_duplicate_scenarios(scen_csv, projectpath=projectpath)

if sweep2:
    validation_scenario = 'sweep2'
    param_dic = get_parameter_space()
    # update dictionary for this scenario
    # setting
    param_dic.update({'annual_EIR': [5, 10, 30]})
    param_dic.update({'seasonality': ['constant', 'high_unimodal']})
    # intervention campaigns
    param_dic.update({'vacc_days': d_vacc_dates})
    param_dic.update({'vacc_deploy_type': ['seasonal'] * len(d_vacc_dates)})
    param_dic.update({'vacc_coverage': [0, d_vaccine_coverage]})
    param_dic.update({'vacc_coverage_multipliers': [1] * len(d_vacc_dates)})
    param_dic.update({'cm_coverage': [d_cm_coverage]})
    param_dic.update({'smc_coverage': [0]})
    param_dic.update({'num_repeated_years': 0})  # for SMC
    param_dic.update({'vacc_target_group': d_vacc_target_group})
    param_dic.update({'cm_target_group': d_cm_target_group})
    param_dic.update({'smc_target_group': d_smc_target_group})
    param_dic.update({'high_access_frac': d_high_access_frac})

    # === RTS,S and no-vaccine scenarios === #
    # vaccine pk/pd
    param_dic.update({'vacc_filename_description': r_vacc_filename_description})
    param_dic.update({'initial_hhs': r_initial_hhs})
    param_dic.update({'initial_nns': r_initial_nns})
    param_dic.update({'booster_hhs': r_booster_hhs})
    param_dic.update({'booster_nns': r_booster_nns})
    param_dic.update({'initial_conc': r_initial_conc})
    param_dic.update({'initial_fast_frac': r_initial_fast_frac})
    param_dic.update({'initial_k1': r_initial_k1})
    param_dic.update({'initial_k2': r_initial_k2})
    param_dic.update({'initial_max_efficacy': r_initial_max_efficacy})
    param_dic.update({'booster_conc': r_booster_conc})
    param_dic.update({'booster_fast_frac': r_booster_fast_frac})
    param_dic.update({'booster_k1': r_booster_k1})
    param_dic.update({'booster_k2': r_booster_k2})
    param_dic.update({'booster_max_efficacy': r_booster_max_efficacy})
    param_dic.update({'vacc_total_time': r_vacc_total_time})

    # create input csvs
    vacc_char_files = create_intervention_inputs(param_dic=param_dic, projectpath=projectpath)
    # create coordinator csv
    df_rtss = create_coordinator_csvs(param_dic=param_dic, base_scenario_filepath=base_scenario_filepath,
                                 vacc_char_files=vacc_char_files)

    # === mAb scenarios === #
    # vaccine pk/pd
    param_dic.update({'vacc_filename_description': m_vacc_filename_description})
    param_dic.update({'initial_max_efficacy': [0.8, 0.9, 0.95]})
    param_dic.update({'booster_max_efficacy': [0.8, 0.9, 0.95]})
    param_dic.update({'initial_hhs': [10, 20, 40]})
    param_dic.update({'booster_hhs': [10, 20, 40]})
    param_dic.update({'initial_nns': m_initial_nns})
    param_dic.update({'booster_nns': m_booster_nns})
    param_dic.update({'initial_conc': m_initial_conc})
    param_dic.update({'booster_conc': m_booster_conc})
    param_dic.update({'initial_fast_frac': m_initial_fast_frac})
    param_dic.update({'booster_fast_frac': m_booster_fast_frac})
    param_dic.update({'initial_k1': m_initial_k1})
    param_dic.update({'booster_k1': m_booster_k1})
    param_dic.update({'initial_k2': m_initial_k2})
    param_dic.update({'booster_k2': m_booster_k2})
    param_dic.update({'vacc_total_time': m_vacc_total_time})

    # create input csvs
    vacc_char_files = create_intervention_inputs(param_dic=param_dic, projectpath=projectpath)
    # create coordinator csv
    df_mab = create_coordinator_csvs(param_dic=param_dic, base_scenario_filepath=base_scenario_filepath,
                                 vacc_char_files=vacc_char_files)

    # === combine and save dataframes === #
    df = combine_input_files(scen_df1=df_rtss, scen_df2=df_mab)
    scen_csv = 'coordinator_files/%s.csv' % validation_scenario
    df.to_csv(os.path.join(projectpath, 'simulation_inputs', scen_csv), index=False)
    remove_duplicate_scenarios(scen_csv, projectpath=projectpath)
