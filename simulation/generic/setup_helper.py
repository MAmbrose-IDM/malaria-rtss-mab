import os
import re
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from dtk.interventions.input_EIR import add_InputEIR
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from malaria.site.input_EIR_by_site import mAb_vs_EIR
from malaria.reports.MalariaReport import add_filtered_report
from malaria.reports.MalariaReport import add_filtered_report
from malaria.reports.MalariaReport import add_event_counter_report
from hbhi.utils import add_monthly_parasitemia_rep_by_year, add_nmf_trt, generate_multiyr_df, read_main_dfs, tryread_df
from simulation.generic.set_up_generic_interventions import InterventionSuite, \
    add_all_interventions  # , update_smc_access_ips


def tryread_df2(input_path, scen_df, id, fname):
    try:
        df = tryread_df(os.path.join(input_path, '%s' % scen_df.at[id, fname]))
    except:
        # print(f"WARNING: {fname} not in scen_df.")
        df = pd.DataFrame()
    return df


def check_colnames(df):
    interventions = ['SMC', 'VACC', 'ITN', 'CM', 'IRS', 'IPTi']
    fnames = [fname for fname in [f'{int}_filename' for int in interventions] if fname not in df.columns]
    print(f"WARNING: {','.join(fnames)} not in scen_df.")


def fix_scen_csv_access_cols(df):
    # if any of the new columns specifying which access group each intervention should target are missing, assume no access-group targeting
    if 'vacc_target_group' not in df.columns:
        df['vacc_target_group'] = 'random'
    if 'smc_target_group' not in df.columns:
        df['smc_target_group'] = 'random'
    if 'cm_target_group' not in df.columns:
        df['cm_target_group'] = 'random'
    return df


def monthly_to_daily_EIR(monthly_EIR):
    """
    modified from
    monthly_to_daily_EIR function, by Erin Zwick
    https://github.com/numalariamodeling/smc-spaq/blob/master/helper_scripts/monthly_to_daily_EIR.py
    """
    if len(monthly_EIR) == 12:
        monthly_EIR = monthly_EIR + [(monthly_EIR[0] + monthly_EIR[-1]) / 2]
    x_monthly = np.linspace(0, 364, num=13, endpoint=True)
    x_daily = np.linspace(0, 364, num=365, endpoint=True)
    EIR = interp1d(x_monthly, monthly_EIR, kind='cubic')
    daily_EIR = EIR(x_daily)
    daily_EIR /= 30

    daily_EIR = daily_EIR.tolist()
    daily_EIR = [0 if eir <0 else eir for eir in daily_EIR]
    return daily_EIR


def setup_setting(cb, scen_df, id, eir_monthly_multipliers, EIR_scale='monthly', cohort_month_shift=0):
    if 'smc_coverage' not in scen_df.columns:
        scen_df['smc_coverage'] = 0
    if 'ipti_coverage' not in scen_df.columns:
        scen_df['ipti_coverage'] = 0
    if 'vacc_coverage' not in scen_df.columns:
        scen_df['vacc_coverage'] = 0
    if 'vacc_characteristics' not in scen_df.columns:
        scen_df['vacc_characteristics'] = 0
    if 'vacc_mode' not in scen_df.columns:
        scen_df['vacc_mode'] = 'constant'
    if 'ipti_mode' not in scen_df.columns:
        scen_df['ipti_mode'] = ''
    if 'vacc_target_group' not in scen_df.columns:
        scen_df['vacc_target_group'] = 'random'
    if 'smc_target_group' not in scen_df.columns:
        scen_df['smc_target_group'] = 'random'
    if 'cm_target_group' not in scen_df.columns:
        scen_df['cm_target_group'] = 'random'

    scen_row = scen_df[scen_df['setting_id'] == id]

    # DEMOGRAPHICS
    cb.update_params({'Demographics_Filenames': [os.path.join(scen_row['demographics_filename'][0])],
                      'Age_Initialization_Distribution_Type': 'DISTRIBUTION_SIMPLE'})

    annual_eir = float(scen_row['annual_EIR'][0])
    monthly_eir_scalers0 = eir_monthly_multipliers[scen_row['seasonality'][0]].tolist()
    # shift monthly EIRs to account for cohort born cohort_month_shift months into the year
    if cohort_month_shift > 0:
        monthly_eir_scalers = monthly_eir_scalers0[cohort_month_shift:12] + monthly_eir_scalers0[0:cohort_month_shift]
    else:
        monthly_eir_scalers = monthly_eir_scalers0
    eir_sum = sum([x for x in monthly_eir_scalers])

    if annual_eir is None:
        annual_eir = eir_sum
        monthly_eir = monthly_eir_scalers
    else:
        annual_eir = annual_eir
        monthly_eir = [(x / eir_sum) * annual_eir for x in monthly_eir_scalers]

    if EIR_scale == 'monthly':
        add_InputEIR(cb, start_day=0, monthlyEIRs=monthly_eir)  # , EIR_type='MONTHLY', , age_dependence='OFF'
    if EIR_scale == 'daily':
        daily_eir = monthly_to_daily_EIR(monthly_eir)
        daily_eir = [(x / sum(daily_eir)) * annual_eir for x in daily_eir]
        add_InputEIR(cb, start_day=0, EIR_type='DAILY',  dailyEIRs=daily_eir)  # , age_dependence='OFF'

    """Adjust Maternal_Antibody_Protection for transmission intensity"""
    Maternal_Antibody_Protection = 0.1327
    mAb = Maternal_Antibody_Protection * mAb_vs_EIR(annual_eir)
    cb.update_params({'Maternal_Antibody_Protection': mAb})

    vacc_char_sub = re.search(r'hh\d+_nn\d+', scen_row['vacc_characteristics'][0])
    if vacc_char_sub:
        vacc_char_vals = vacc_char_sub[0]
    else:
        vacc_char_vals = 'NA'
    return {'seasonality': scen_row['seasonality'][0],
            'Annual EIR': annual_eir,
            'mAb': mAb,
            'cm_coverage': scen_row['cm_coverage'][0],
            'vacc_coverage': scen_row['vacc_coverage'][0],
            'smc_coverage': scen_row['smc_coverage'][0],
            'ipti_coverage': scen_row['ipti_coverage'][0],
            'vacc_target_group': scen_row['vacc_target_group'][0],
            'smc_target_group': scen_row['smc_target_group'][0],
            'cm_target_group': scen_row['cm_target_group'][0],
            'frac_high_access': scen_row['frac_high_access'][0],
            'vacc_mode': scen_row['vacc_mode'][0],
            'ipti_mode': scen_row['ipti_mode'][0],
            'vacc_char': vacc_char_vals
            }


def setup_simulation(years):
    # BASIC SETUP
    cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')

    # Run time
    cb.update_params({'Simulation_Duration': years * 365 + 1})

    # Logging
    cb.update_params({
        'logLevel_JsonConfigurable': 'ERROR',
        'logLevel_VectorHabitat': 'ERROR',
        'logLevel_StandardEventCoordinator': 'ERROR',
        'logLevel_SusceptibilityMalaria': 'ERROR'
    })

    # Demographics
    cb.update_params({
        'Enable_Birth': 0,
        'Enable_Births': 0,
        'Enable_Demographics_Birth': 0,
        'Enable_Demographics_Risk': 0,
        'Birth_Rate_Dependence': 'FIXED_BIRTH_RATE',
        'Enable_Initial_Prevalence': 0,
        'Age_Dependent_Biting_Risk_Type': 'OFF',
        'x_Birth': 1,
        'x_Base_Population': 1,
        'Maternal_Antibodies_Type': 'SIMPLE_WANING',  # 'CONSTANT_INITIAL_IMMUNITY',

        "Age_Initialization_Distribution_Type": 'DISTRIBUTION_COMPLEX',
        'Disable_IP_Whitelist': 1
    })

    # Serialization (none)
    cb.update_params({
        'Serialization_Type': 'NONE',
        'Serialized_Population_Writing_Type': 'NONE',
        'Serialized_Population_Reading_Type': 'NONE',
    })

    # Vector and climate
    cb.update_params({
        "Vector_Species_Names": [],
        'x_temporary_Larval_Habitat': 0,
        'Climate_Model': 'CLIMATE_CONSTANT'  # "CLIMATE_BY_DATA"
    })

    # Reporting
    cb.update_params({
        'Enable_Default_Reporting': 0,
        'Enable_Property_Output': 0,
        'Enable_Vector_Species_Report': 0,
        'Report_Detection_Threshold_Blood_Smear_Parasites': 10,  # 50
        "Parasite_Smear_Sensitivity": 0.02,  # 50/uL
        'RDT_Sensitivity': 0.1
    })
    cb.update_params({
        "Report_Event_Recorder": 1,
        "Report_Event_Recorder_Individual_Properties": ['VaccineStatus', 'AccessToInterventions', 'vaccine_selected'],
        "Report_Event_Recorder_Ignore_Events_In_List": 0,
        "Report_Event_Recorder_Events": ['Births', 'PropertyChange'],
        'Custom_Individual_Events': ['Received_Treatment', 'Received_Severe_Treatment',
                                     'Received_NMF_Treatment', 'Received_Vaccine', 'Received_Campaign_Drugs',
                                     'event_add_new_vaccine']
    })

    # Filtered report (all years):
    num_year = 1
    start = 1  # 1 + (years - num_year) * 365
    end = 1 + years * 365
    add_filtered_report(cb, start=start, end=end)

    # CUSTOM REPORTS
    add_event_counter_report(cb, event_trigger_list=['Received_Treatment', 'Received_Severe_Treatment',
                                                     'Received_Vaccine', 'Received_Campaign_Drugs'],
                             duration=years * 365 + 1)
    #
    # for year in range(years):
    #     add_monthly_parasitemia_rep_by_year(cb, num_year=years, tot_year=years,
    #                                         sim_start_year=2020,
    #                                         yr_plusone=True, prefix='Monthly')
    # add_monthly_parasitemia_rep_by_year(cb, num_year=years, tot_year=years,
    #                                     sim_start_year=2020,
    #                                     yr_plusone=True,
    #                                     age_bins=[1, 5, 120],
    #                                     prefix='FineMonthly')
    return cb


def add_generic_interventions(cb, projectpath, scen_df, id, ds_name='run_col', cohort_month_shift=0):
    # INTERVENTIONS
    int_suite = InterventionSuite()

    int_suite.hs_ds_col = ds_name
    int_suite.ipti_ds_col = ds_name
    int_suite.smc_ds_col = ds_name

    input_path = os.path.join(projectpath, 'simulation_inputs')
    hs_df = tryread_df2(input_path, scen_df, id, 'CM_filename')
    ipti_df = tryread_df2(input_path, scen_df, id, 'IPTi_filename')
    smc_df = tryread_df2(input_path, scen_df, id, 'SMC_filename')
    vacc_char_df = tryread_df2(input_path, scen_df, id, 'vacc_characteristics')
    vacc_df = tryread_df2(input_path, scen_df, id, 'VACC_filename')

    add_all_interventions(cb,
                          int_suite,
                          my_ds='run',  # to not subset dataframe as each intervention csv is generic
                          high_access_ip_frac=scen_df.at[id, 'frac_high_access'],
                          vacc_target_group=scen_df.at[id, 'vacc_target_group'],
                          smc_target_group=scen_df.at[id, 'smc_target_group'],
                          cm_target_group=scen_df.at[id, 'cm_target_group'],
                          hs_df=hs_df,
                          ipti_df=ipti_df,
                          smc_df=smc_df,
                          vacc_char_df=vacc_char_df,
                          vacc_df=vacc_df,
                          cohort_month_shift=cohort_month_shift)

    return {'Scenario_id': id}
