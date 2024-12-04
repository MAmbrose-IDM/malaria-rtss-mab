import logging
from shutil import copyfile
import os
import time
import numpy as np
import pandas as pd
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.interventions.property_change import change_individual_property
from simtools.SetupParser import SetupParser
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.Utilities.Experiments import retrieve_experiment
from malaria.interventions.health_seeking import add_health_seeking
from malaria.reports.MalariaReport import add_filtered_report, add_summary_report
from dtk.interventions.input_EIR import add_InputEIR
from simtools.ModBuilder import ModBuilder, ModFn
from malaria.site.input_EIR_by_site import mAb_vs_EIR
from malaria.reports.MalariaReport import add_event_counter_report
from hbhi.set_up_general import set_spaq_params
from malaria.interventions.malaria_vaccine import add_vaccine
from hbhi.utils import add_monthly_parasitemia_rep_by_year, add_nmf_trt

logger = logging.getLogger(__name__)
if os.name == "posix":
    SetupParser.default_block = 'NUCLUSTER'
    burnin_id = ""
else:
    SetupParser.default_block = 'HPC'
    burnin_id = ""
cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')


def set_EIR(cb, annual_eir, EIR_scale):
    Maternal_Antibody_Protection = 0.1327
    mAb = Maternal_Antibody_Protection * mAb_vs_EIR(annual_eir)
    cb.update_params({'Maternal_Antibody_Protection': mAb})

    """Use constant transmission for testing"""
    if EIR_scale == 'monthly':
        add_InputEIR(cb, start_day=0, EIR_type='MONTHLY', monthlyEIRs=[annual_eir / 12] * 12)
    if EIR_scale == 'daily':
        add_InputEIR(cb, start_day=0, EIR_type='DAILY', dailyEIRs=[annual_eir / 365] * 365)

    return {'Annual EIR': annual_eir,
            'EIR_scale': EIR_scale,
            'mAb': mAb}


def setup_simulation(cb, years, sim_start_year, interval=30):
    cb.update_params({'Demographics_Filenames': [os.path.join('Demographics/Aba North_demographics_wVaxSMC.json')],
                      'Age_Initialization_Distribution_Type': 'DISTRIBUTION_COMPLEX'})
    cb.update_params({'Vector_Species_Names': [],
                      'Enable_Vital_Dynamics': 1,
                      'Enable_Births': 1,
                      'Birth_Rate_Dependence': 'FIXED_BIRTH_RATE',
                      'x_Birth': 10,
                      'x_Base_Population': 1,
                      'Disable_IP_Whitelist': 1,
                      'Maternal_Antibodies_Type': 'SIMPLE_WANING',  # 'CONSTANT_INITIAL_IMMUNITY',
                      "Incubation_Period_Distribution": "EXPONENTIAL_DISTRIBUTION",
                      'Climate_Model': 'CLIMATE_CONSTANT',
                      'Parasite_Smear_Sensitivity': 0.02,
                      'logLevel_JsonConfigurable': 'ERROR',
                      'Simulation_Duration': years * 365
                      })

    ## REPORTS
    add_filtered_report(cb, start=0, end=years * 365)
    for year in range(years):
        add_summary_report(cb, start=(365 * year), interval=interval, duration_days=365,
                           age_bins=[1, 2, 3, 4, 5, 125], description=f'FineMonthly{year + sim_start_year}',
                           parasitemia_bins=[10, 50, 1e9])

        add_summary_report(cb, start=(365 * year), interval=interval, duration_days=365,
                           age_bins=[0.25, 5, 125], description=f'Monthly{year + sim_start_year}',
                           parasitemia_bins=[10, 50, 1e9])

    add_event_counter_report(cb, event_trigger_list=['Received_Vaccine'],
                             duration=years * 365 + 1)

    if serialize:
        cb.update_params({
            'Serialization_Time_Steps': [365 * years],
            'Serialization_Type': 'TIMESTEP',
            'Serialized_Population_Writing_Type': 'TIMESTEP',
            'Serialization_Mask_Node_Write': 0,
            'Serialization_Precision': 'REDUCED'
        })
    else:
        cb.update_params({
            'Serialization_Type': 'NONE',
            'Serialized_Population_Writing_Type': 'NONE'
        })

    if pull_from_serialization:
        expt = retrieve_experiment(burnin_id)
        ser_df = pd.DataFrame([x.tags for x in expt.simulations])
        ser_df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])
        ser_path = ser_df['outpath'].values[0]

        cb.update_params({
            'Serialized_Population_Reading_Type': 'READ',
            'Serialized_Population_Path': os.path.join(ser_path, 'output'),
            'Serialized_Population_Filenames': ['state-03650.dtk'],
            'Enable_Random_Generator_From_Serialized_Population': 0,
            'Serialization_Mask_Node_Read': 0
        })
    else:
        cb.update_params({
            'Serialized_Population_Reading_Type': 'NONE'
        })


def fill_list(x_list, rtss_touchpoints):
    if len(x_list) is not len(rtss_touchpoints):
        return [x_list[0]] * len(rtss_touchpoints)
    else:
        return x_list


def validate_zip(x_list, rtss_touchpoints):
    false_lengths = [i for i, x in enumerate(x_list) if len(x) != len(rtss_touchpoints)]
    if len(false_lengths) == 0:
        print("zip valid, all same lengths")
    else:
        raise ValueError("zip not valid, items not of same length!")


def case_management(cb, cm_cov_U5=0.418, cm_cov_adults=0.418):
    add_health_seeking(cb, start_day=0,
                       targets=[{'trigger': 'NewClinicalCase', 'coverage': cm_cov_U5, 'agemin': 0, 'agemax': 5,
                                 'seek': 1, 'rate': 0.3},
                                {'trigger': 'NewClinicalCase', 'coverage': cm_cov_adults, 'agemin': 5, 'agemax': 100,
                                 'seek': 1, 'rate': 0.3}],
                       drug=['Artemether', 'Lumefantrine'], duration=-1)
    add_health_seeking(cb, start_day=0,
                       targets=[{'trigger': 'NewSevereCase', 'coverage': 0.49, 'seek': 1, 'rate': 0.5}],
                       drug=['Artemether', 'Lumefantrine'], duration=-1,
                       broadcast_event_name='Received_Severe_Treatment')

    return {'cm_cov_U5': cm_cov_U5,
            'cm_cov_adults': cm_cov_adults}


def rtss_intervention(cb, start_days, coverage_levels, rtss_types, rtss_touchpoints, deploy_type='EPI',
                      ageminmax=[0.25, 5], repetitions=5, tsteps_btwn_repetitions=365,
                      delay_distribution_name=None, Std_Dev_list=[7],
                      Decay_Time_Constant_list=[592.4066511650307], Initial_Effect_list=[0.8],
                      Waning_Class_list=["WaningEffectExponential"]):
    if not deploy_type == 'EPI':
        """Use campaign-style deployment
        Currently only supports constant efficacy and coverage for each round
        Booster not constrained to previous vaccine status
        """
        rtss_event_names = []
        vaccine_params = {"Waning_Config": {"Initial_Effect": Initial_Effect[0],
                                            "Decay_Time_Constant": Decay_Time_Constant[0], "class": Waning_Class[0]}}
        for start_day in start_days:
            add_vaccine(cb,
                        vaccine_type='RTSS',
                        vaccine_params=vaccine_params,
                        start_days=[start_day],
                        coverage=coverage_levels[0],
                        repetitions=repetitions,
                        tsteps_btwn_repetitions=tsteps_btwn_repetitions,
                        target_group={'agemin': ageminmax[0], 'agemax': ageminmax[1]})
    else:
        """
        Use birthtriggered deployment for vaccine distributed along EPI
        Allows round-specific coveraged and vaccine parameters
        Booster vaccine is only given if vaccine status is 'GotVaccine'
        """
        rtss_event_names = [f'RTSS_{x + 1}_eligible' for x in range(len(rtss_touchpoints))]

        """Per default take first item if length does not match"""
        coverage_levels = fill_list(coverage_levels, rtss_touchpoints)
        Std_Dev_list = fill_list(Std_Dev_list, rtss_touchpoints)
        Initial_Effect_list = fill_list(Initial_Effect_list, rtss_touchpoints)
        Decay_Time_Constant_list = fill_list(Decay_Time_Constant_list, rtss_touchpoints)
        Waning_Class_list = fill_list(Waning_Class_list, rtss_touchpoints)
        if len(rtss_types) != len(rtss_touchpoints):
            raise ValueError('Please provide rtss_type for each rtss_touchpoint')

        validate_zip([coverage_levels, rtss_types, rtss_event_names, Std_Dev_list,
                      Initial_Effect_list, Decay_Time_Constant_list, Waning_Class_list], rtss_touchpoints)
        for tp_time_trigger, coverage, vtype, event_name, std, init_eff, decay_t, decay_c in \
                zip(rtss_touchpoints, coverage_levels, rtss_types, rtss_event_names, Std_Dev_list,
                    Initial_Effect_list, Decay_Time_Constant_list, Waning_Class_list):

            vaccine_params = {"Waning_Config": {"Initial_Effect": init_eff,
                                                "Decay_Time_Constant": decay_t,
                                                "class": decay_c}}

            if delay_distribution_name == "GAUSSIAN_DISTRIBUTION":
                delay_distribution = {"Delay_Period_Distribution": "GAUSSIAN_DISTRIBUTION",
                                      "Delay_Period_Gaussian_Mean": tp_time_trigger,
                                      "Delay_Period_Gaussian_Std_Dev": std}
                tp_time_trigger = None
            else:
                delay_distribution = None

            if not vtype == 'booster':
                add_vaccine(cb,
                            vaccine_type='RTSS',
                            vaccine_params=vaccine_params,
                            start_days=start_days,
                            coverage=coverage,
                            delay_distribution=delay_distribution,
                            triggered_delay=tp_time_trigger,
                            trigger_condition_list=[event_name],
                            birthtriggered=True)
                change_individual_property(cb,
                                           target_property_name='VaccineStatus',
                                           target_property_value='GotVaccine',
                                           trigger_condition_list=['Received_Vaccine'],
                                           blackout_flag=False)
            else:
                """Booster vaccine, targeting 80% of those who got the vaccine"""
                add_vaccine(cb,
                            vaccine_type='RTSS',
                            vaccine_params=vaccine_params,
                            start_days=start_days,
                            coverage=0.8,
                            delay_distribution=delay_distribution,
                            triggered_delay=tp_time_trigger,
                            trigger_condition_list=[event_name],
                            ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                            birthtriggered=True)

                change_individual_property(cb,
                                           target_property_name='BoosterStatus',
                                           target_property_value='GotBooster',
                                           ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                                           trigger_condition_list=['Received_Vaccine'],
                                           blackout_flag=False)

    cb.update_params({'Report_Event_Recorder': 1,
                      'Report_Event_Recorder_Ignore_Events_In_List': 0,
                      'Report_Event_Recorder_Events': ['Births', 'PropertyChange',
                                                       'Received_Vaccine'] + rtss_event_names,
                      'Report_Event_Recorder_Individual_Properties': ['VaccineStatus', 'BoosterStatus'],
                      'Custom_Individual_Events': ['Received_Vaccine']
                      })

    return {'deploy_type': deploy_type,
            'coverage': coverage_levels,
            'rtss_touchpoints': rtss_touchpoints,
            'rtss_types': rtss_types,
            'Initial_Effect_list': Initial_Effect_list,
            'Waning_Class_list': Waning_Class_list}


if __name__ == "__main__":

    serialize = False
    pull_from_serialization = False

    if serialize:
        expname = "test_burnin_v6"
        num_seeds = 1
        burnin = 10
        years = burnin
        sim_start_year = 2010
        rtss_start = None
        report_interval = 365
    else:
        expname = "generic_rtss_sweep"
        num_seeds = 5
        if pull_from_serialization:
            burnin = 1
            sim_start_year = 2020
            rtss_start = 365 + 365 * burnin
        else:
            burnin = 10
            sim_start_year = 2020 - burnin
            rtss_start = 365 + 365 + 365 * burnin
        years = 7 + burnin
        report_interval = 30

    season_start_day = 1  # Start first day of the year
    setup_simulation(cb, years, sim_start_year, report_interval)
    builder = ModBuilder.from_list([[ModFn(set_EIR, annual_eir, 'daily'),
                                     ModFn(rtss_intervention,
                                           start_days=[rtss_start + season_start_day],
                                           coverage_levels=[cov],
                                           deploy_type='EPI',
                                           rtss_touchpoints=touchpoints,
                                           rtss_types=type,
                                           delay_distribution_name="CONSTANT_DISTRIBUTION",
                                           # Std_Dev_list= [std],
                                           # Decay_Time_Constant_list=[t_decay],
                                           Initial_Effect_list=[eff, eff_booster]),
                                     ModFn(case_management, cm_cov_U5),
                                     ModFn(set_spaq_params),
                                     ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
                                     ModFn(DTKConfigBuilder.set_param, 'coverage', cov),
                                     ]
                                    for annual_eir in [0.23, 2.06, 3.61, 5.75, 9.3, 16.54, 26.72, 36.68]
                                    for cov in [0, 0.6, 0.8, 1]  # np.linspace(0, 1, 9)
                                    for std in [7]  # [7, 14, 21]
                                    for eff in [0.8]  # np.linspace(0, 1, 9)
                                    for eff_booster in [0.4]
                                    for touchpoints in [[274, 730]]  # [[274],[274, 730]]
                                    for type in [['simple', 'booster']]  # [['simple'],['simple', 'booster']]
                                    # for t_decay in [592.4066511650307]
                                    for cm_cov_U5 in [0.45, 0.6, 0.8]  # np.linspace(0, 1, 11)
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
