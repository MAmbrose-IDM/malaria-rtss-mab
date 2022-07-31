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
validation_scenario = 'Phase3_wBooster_p2'  # 'SweepEIR_wBooster', 'SweepEIR_noBooster', 'Phase3_wBooster', Phase3_noBooster', 'ChandramohanTrial'

# scenario-specific parameters
# if 'SweepEIR' in validation_scenario:
    # annual_EIR = [0.23, 0.91, 2.06, 3.61, 5.75, 9.3, 16.54, 26.72, 36.68]  # np.linspace(5,100,4)
    # cm_coverage = [0.45]  # np.linspace(0.2, 1, 5)[:-1]
    # seasonality = ['constant']
    # intervention_correlation = [0]
    # vacc_coverage = [0, 1]
    # smc_coverage = [0]
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
    # param_dic = get_parameter_space()
    # # update dictionary for this scenario
    # # create input csvs
    # vacc_char_files = create_intervention_inputs(param_dic=param_dic, projectpath=projectpath)
    # # create coordinator csv
    # df = create_coordinator_csvs(param_dic=param_dic, base_scenario_filepath=base_scenario_filepath,
    #                              vacc_char_files=vacc_char_files)
    # scen_csv = 'coordinator_files/validation_%s.csv' % validation_scenario
    # df.to_csv(os.path.join(projectpath, 'simulation_inputs', scen_csv), index=False)
    # remove_duplicate_scenarios(scen_csv, projectpath=projectpath)

if 'Phase3' in validation_scenario:
    booster_coverage_vacc = 1 if 'wBooster' in validation_scenario else 0
    param_dic = get_parameter_space()
    # update dictionary for this scenario
    param_dic.update({'cm_coverage': [0.9]})
    param_dic.update({'vacc_coverage': [0, 1]})
    param_dic.update({'smc_coverage': [0]})
    param_dic.update({'vacc_filename_description': 'Phase3'})
    param_dic.update({'vacc_days': [365, 912]})
    param_dic.update({'vacc_coverage_multipliers': [1, booster_coverage_vacc]})
    param_dic.update({'initial_hhs': [30, 40]})  # [20, 30, 40, 60]
    param_dic.update({'initial_nns': [2]})   # [1.4, 2, 4]
    param_dic.update({'booster_hhs': [30, 40]})  # [20, 30, 40, 60]
    param_dic.update({'booster_nns': [2]})   # [1.4, 2, 4]
    param_dic.update({'initial_max_efficacy': 0.8})
    param_dic.update({'booster_max_efficacy': 0.8})  # 0.6, 0.7, 0.8
    # create input csvs
    vacc_char_files = create_intervention_inputs(param_dic=param_dic, projectpath=projectpath)

    phase3_annualEIRs = pd.read_csv(os.path.join(datapath, 'rtss_phase3/phase3_annualEIRs.csv'))
    df = pd.DataFrame()
    for r, row in phase3_annualEIRs.iterrows():
        param_dic.update({'annual_EIR': [row['annual_eir']]})
        param_dic.update({'seasonality': [row['site']]})

        # create coordinator csv
        df_cur = create_coordinator_csvs(param_dic=param_dic, base_scenario_filepath=base_scenario_filepath,
                                         vacc_char_files=vacc_char_files, high_access_frac=0)
        # combine dataframes from sites
        df = combine_input_files(scen_df1=df, scen_df2=df_cur)
    scen_csv = 'coordinator_files/validation_%s.csv' % validation_scenario
    df.to_csv(os.path.join(projectpath, 'simulation_inputs', scen_csv), index=False)
    remove_duplicate_scenarios(scen_csv, projectpath=projectpath)

elif 'ChandramohanTrial' in validation_scenario:
    param_dic = get_parameter_space()
    # update dictionary for this scenario
    param_dic.update({'annual_EIR': [6, 10, 15]})
    param_dic.update({'seasonality': ['higher_JuneStart']})
    param_dic.update({'cm_coverage': [0.8, 0.9]})
    param_dic.update({'vacc_coverage': [0, 0.934]})
    param_dic.update({'smc_coverage': [0, 0.89]})
    param_dic.update({'smc_start_month': 14.8})
    param_dic.update({'num_repeated_years': 2})
    param_dic.update({'vacc_coverage_multipliers': [1, 0.95*0.934, 0.95*0.934]})
    param_dic.update({'vacc_filename_description': 'ChandramohanTrial'})
    param_dic.update({'vacc_days': [365, 365*2, 365*3]})
    param_dic.update({'initial_hhs': [40]})
    param_dic.update({'initial_nns': [2]})
    param_dic.update({'booster_hhs': [40]})
    param_dic.update({'booster_nns': [2]})
    param_dic.update({'initial_max_efficacy': 0.8})
    param_dic.update({'booster_max_efficacy': 0.8})   # 0.6  , 0.7, 0.75, 0.8

    vacc_char_files = create_intervention_inputs(param_dic=param_dic, projectpath=projectpath)
    # create coordinator csv
    df = create_coordinator_csvs(param_dic=param_dic, base_scenario_filepath=base_scenario_filepath,
                                 vacc_char_files=vacc_char_files, high_access_frac=0)
    scen_csv = 'coordinator_files/validation_%s.csv' % validation_scenario
    df.to_csv(os.path.join(projectpath, 'simulation_inputs', scen_csv), index=False)
    remove_duplicate_scenarios(scen_csv, projectpath=projectpath)


else:
    print('validation scenario not recognized')
    # annual_EIR = [0]
    # cm_coverage = [0]
    # seasonality = ['constant']
    # intervention_correlation = [0]
    # vacc_coverage = [0]
    # smc_coverage = [0]
    # # vacc arguments
    # booster_coverage_vacc = 0
    # vacc_filename_description = 'ERROR'
    # vacc_3rd_dose_age = 365
    # vacc_coverage_multipliers = [1]
    # vacc_rounds = [3]
    # vacc_types = ['simple']
    # vacc_initial_killing = [initial_effect_vacc]
    # vacc_days = [vacc_3rd_dose_age]
    # param_dic = get_parameter_space()
    # # update dictionary for this scenario
    # # create input csvs
    # vacc_char_files = create_intervention_inputs(param_dic=param_dic, projectpath=projectpath)
    # # create coordinator csv
    # df = create_coordinator_csvs(param_dic=param_dic, base_scenario_filepath=base_scenario_filepath,
    #                              vacc_char_files=vacc_char_files)
    # scen_csv = 'coordinator_files/validation_%s.csv' % validation_scenario
    # df.to_csv(os.path.join(projectpath, 'simulation_inputs', scen_csv), index=False)
    # remove_duplicate_scenarios(scen_csv, projectpath=projectpath)



















