import os
import sys
import logging
import pandas as pd
import time
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser
from hbhi.set_up_general import set_spaq_params
from simulation.generic.setup_helper import setup_setting, setup_simulation, add_generic_interventions
from simulation.generic.update_drug_params import update_drug_config

sys.path.append('../')
from simulation.load_paths import load_box_paths

logger = logging.getLogger(__name__)
if os.name == "posix":
    SetupParser.default_block = 'NUCLUSTER'
else:
    SetupParser.default_block = 'HPC'

datapath, projectpath = load_box_paths(parser_default=SetupParser.default_block)
scenario_csv = os.path.join(projectpath, 'simulation_inputs', 'scenario_files', 'generic', 'generic_settings_test.csv')
scen_df = pd.read_csv(scenario_csv)

num_seeds = 5
years = 6  # for general runs: 10, for SweepEIR: 20, for KintampoPhase3 and ChandramohanTrial: 6
use_12_cohorts_flag = False  # True  # simulate one cohort born in each month (seasonality and SMC dates updated accordingly)
sdxpyr_kill = [1.58, 1.58, 1.2]
sdxpyr_c50 = [10.83, 10.83, 10.83]
sdxpyr_decayT1_2 = [10, 12, 15]

expname = 'rtss_validation_ChandramohanTrial_sweepSP3_sahelSeason'  # 'rtss_validation_ChandramohanTrial_ageEIR'  # 'generic_validation_noSARisk'  # 'GEN_projections_TEST'
ds_name = 'run_col'


if __name__ == "__main__":
    # setup basic config for simulations
    cb = setup_simulation(years)

    # REFERENCE CSV FOR TRANSMISSION SETTINGS
    scen_csv = 'generic_validation_ChandramohanTrial.csv'  # "generic_validation_SweepEIR_noBooster.csv"   #  "generic_settings_test.csv"  #  "generic_validation_kintampo_approximation.csv"  # "generic_validation_cases_averted_by_age.csv"  #
    scen_df = pd.read_csv(os.path.join(projectpath, 'simulation_inputs/scenario_files/generic', scen_csv),
                     encoding='latin')
    scen_df['DS_Name'] = scen_df['setting_id'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8')
    scen_df = scen_df.set_index('DS_Name')
    eir_monthly_multipliers = pd.read_csv(os.path.join(projectpath, 'simulation_inputs/scenario_files/generic', 'Seasonality',
                           'seasonality_eir_multipliers.csv'))
    if use_12_cohorts_flag:
        cohort_month_shift_values = range(12)
    else:
        cohort_month_shift_values = [0]


    # from malaria.reports.MalariaReport import add_survey_report
    # add_survey_report(cb, [100, 300, 780], reporting_interval=1)

    # BUILDER
    builder = ModBuilder.from_list([[ModFn(setup_setting,
                                           scen_df=scen_df,
                                           id=id,
                                           eir_monthly_multipliers=eir_monthly_multipliers,
                                           EIR_scale='monthly',
                                           cohort_month_shift=cohort_month_shift),
                                     ModFn(add_generic_interventions,
                                           projectpath=projectpath,
                                           scen_df=scen_df,
                                           id=id,
                                           cohort_month_shift=cohort_month_shift),
                                     ModFn(update_drug_config,  # note: this call updates all other drug params to the values specified in the update_drug_config file
                                           list_drug_param=['SulfadoxinePyrimethamine', 'Drug_PKPD_C50'],
                                           value=sdxpyr_c50[sp],
                                           list_drug_param2=['SulfadoxinePyrimethamine', 'Max_Drug_IRBC_Kill'],
                                           value2=sdxpyr_kill[sp],
                                           list_drug_param3=['SulfadoxinePyrimethamine', 'Drug_Decay_T1'],
                                           value3=sdxpyr_decayT1_2[sp],
                                           list_drug_param4=['SulfadoxinePyrimethamine', 'Drug_Decay_T2'],
                                           value4=sdxpyr_decayT1_2[sp],
                                           ),
                                     ModFn(DTKConfigBuilder.set_param, 'SP_c50', sdxpyr_c50[sp]),
                                     ModFn(DTKConfigBuilder.set_param, 'SP_kill', sdxpyr_kill[sp]),
                                     ModFn(DTKConfigBuilder.set_param, 'SP_decayT', sdxpyr_decayT1_2[sp]),
                                     ModFn(DTKConfigBuilder.set_param, 'Setting_id', id),
                                     ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
                                     ModFn(DTKConfigBuilder.set_param, 'Cohort_birth_month', cohort_month_shift),
                                     ]
                                    # for id in ['HX30', 'HX31', 'HX32']#'HX3', 'HX219', 'HX255', 'HX0', 'HX216', 'HX252']  # Arbitrary setting for testing
                                    for id in scen_df.index
                                    for cohort_month_shift in cohort_month_shift_values
                                    for sp in range(len(sdxpyr_c50))
                                    for x in range(num_seeds)
                                    ])

    run_sim_args = {
        'exp_name': expname,
        'config_builder': cb,
        'exp_builder': builder
    }

    SetupParser.init()
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)

    if os.name == "posix":
        time.sleep(20)
        exp_manager.wait_for_finished(verbose=True)
        assert (exp_manager.succeeded())
