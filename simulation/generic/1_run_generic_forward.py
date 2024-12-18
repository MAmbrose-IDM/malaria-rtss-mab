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

sys.path.append('../../')
from simulation.generic.setup_helper import setup_setting, setup_simulation, add_generic_interventions, check_colnames, fix_scen_csv_access_cols
from simulation.load_paths import load_box_paths

logger = logging.getLogger(__name__)
if os.name == "posix":
    SetupParser.default_block = 'NUCLUSTER'
else:
    SetupParser.default_block = 'HPC'

datapath, projectpath = load_box_paths(parser_default=SetupParser.default_block)

# COORDINATOR CSV FOR SETTING UP SIMULATION SCENARIOS
""" reference csv's examples
TESTING
test_mAb.csv
VALIDATION
validation_Phase3_noBooster_cleaned.csv
validation_Phase3_wBooster_cleaned.csv
validation_ChandramohanTrial_cleaned.csv
NO-INTERVENTIONS
___.csv
"""
num_seeds = 1
scen_csv = 'test.csv'
expname = 'mAb_%s' % scen_csv.replace('.csv', '_seeds%i' % num_seeds)
ds_name = 'run_col'

years = 3  #!8  # for general runs: 8, for Phase3 and ChandramohanTrial: 6
use_12_cohorts_flag = False #!! True  # for general runs: True, for validation experiments: False
                            # If True, simulates one cohort born in each month (seasonality and SMC dates updated accordingly)

if __name__ == "__main__":
    # setup basic config for simulations
    cb = setup_simulation(years)
    scen_df = pd.read_csv(os.path.join(projectpath, 'simulation_inputs/coordinator_files', scen_csv),
                          encoding='latin')

    check_colnames(scen_df)

    # update columns for old versions of scen_csv
    fix_scen_csv_access_cols(scen_df)

    # use same csv for IPTi simulation and postprocessing, if postprocessing run without ipti
    # do nothing if ipti not specified
    try:
        if str(scen_df.at[1, 'ipti_postprocess']).lower() == 'true':
            scen_df = scen_df.loc[scen_df['ipti_coverage'] == 0]
    except:
        pass

    scen_df['DS_Name'] = scen_df['setting_id'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode(
        'utf-8')
    scen_df = scen_df.set_index('DS_Name')
    eir_monthly_multipliers = pd.read_csv(
        os.path.join(projectpath, 'simulation_inputs', 'Seasonality',
                     'seasonality_eir_multipliers.csv'))

    if use_12_cohorts_flag:
        cohort_month_shift_values = range(12)
    else:
        cohort_month_shift_values = [0]

    # REPORTS
    summaryreport = False
    if summaryreport:
        from malaria.reports.MalariaReport import add_summary_report
        add_summary_report(cb, start=0, interval=365, duration_days=365, age_bins=[1, 125], description='U1')
        add_summary_report(cb, start=0, interval=365 * 5, duration_days=365 * 5, age_bins=[0.25, 5, 125], description='U5')

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
                                     # ModFn(set_spaq_params),  <-- currently the dtk-tools-malaria parameters are the most updated
                                     ModFn(DTKConfigBuilder.set_param, 'Setting_id', id),
                                     ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
                                     ModFn(DTKConfigBuilder.set_param, 'Cohort_birth_month', cohort_month_shift),
                                     ]
                                    # for id in ['HX3', 'HX7', 'HX11']  # Arbitrary setting for testing
                                    for id in scen_df.index
                                    for cohort_month_shift in cohort_month_shift_values
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
