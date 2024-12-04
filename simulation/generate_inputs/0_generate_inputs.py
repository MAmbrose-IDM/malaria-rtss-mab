import os
import sys
import numpy as np
import pandas as pd
import itertools

sys.path.append('../')
from simulation.load_paths import load_box_paths


def get_parameter_space():
    """all parameter until deployment_parameters are full factorial combinations
    When creating master csv either SMC or IPTi coverages will be selected and the other one set to 0
    """
    parameter_space = {
        'annual_EIR': [1, 10, 30, 60],  # np.linspace(5,100,4)
        'cm_coverage': [0.6, 0.8],  # np.linspace(0.2, 1, 5)[:-1]
        'seasonality': ['constant', 'moderate_unimodal', 'high_unimodal'],
        'rtss_target_group': ['random'],  # specify whether RTS,S is targeted toward an access group or given at random
        'smc_target_group': ['random'],  # specify whether SMC is targeted toward an access group or given at random
        'cm_target_group': ['random'],  # specify whether CM is targeted toward an access group or given at random
        'rtss_coverage': [0, 0.6, 0.8],
        'rtss_mode': ['constant'],  # options are: 'constant' (EPI main dose and booster dose),
        #              'campboost' (EPI main dose and 1 campaign booster dose),
        #              'campboost2' (EPI main dose and 2 campaign booster doses)
        'smc_coverage': [0, 0.6, 0.8],
        'ipti_coverage': [0, 0.6, 0.8],
        'deployment_parameters': '',
        'rtss_age_days': [274, 730],
        'rtss_booster1_min_age': [730],  # per default minimum age 24 months to receive 1st booster
        'booster_coverage_rtss': 0.8,
        'max_age_rtss': 5,
        'decay_time_constant_rtss': 592.4066512,
        'decay_class_rtss': 'WaningEffectExponential',
        'ipti_mode': ['basic'],  # options are: 'basic', 'extended3tp', 'extended'
        'ipti_touchpoints': [61, 91, 274],
        'ipti_postprocess': True,  # if true IPTi wont be simulated but scenario csv used for postprocessing
        'num_smc_rounds': 4,
        'smc_start_month': 7,
        'num_repeated_years': 20,
        # # old correlation parameters
        # 'intervention_correlation': [''],
        # 'intervention_cor_cm': [''], # If False, CM coverage is independent from SMC and RTSS coverage (given intervention_correlation=1)
        # 'antiaccess_high': [''], # Specify whether 'SMC' or 'RTSS' is the high coverage intervention (given intervention_correlation=1)
    }
    return parameter_space


def create_scenarios_mastercsv(param_dic, fname_out='generic_settings_test', smcipti='SMC', high_access_frac=None):
    annual_EIR = param_dic['annual_EIR']
    cm_coverage = param_dic['cm_coverage']
    seasonality = param_dic['seasonality']
    rtss_coverage = param_dic['rtss_coverage']
    rtss_mode = param_dic['rtss_mode']
    rtss_booster1_min_age = param_dic['rtss_booster1_min_age']
    smc_coverage = param_dic['smc_coverage']
    ipti_coverage = param_dic['ipti_coverage']
    ipti_mode = param_dic['ipti_mode']
    rtss_target_group = param_dic['rtss_target_group']
    smc_target_group = param_dic['smc_target_group']
    cm_target_group = param_dic['cm_target_group']

    # Create full factorial
    if smcipti == 'SMC':
        ipti_coverage = [0]
    elif smcipti == 'IPTi':
        smc_coverage = [0]
    elif smcipti == 'RTSS':
        smc_coverage = [0]
        ipti_coverage = [0]
    else:
        print(f'Warning: {smcipti} not specified, run RTSS only')
        smcipti = 'RTSS'  # RTSS only no SMC or IPTi
        smc_coverage = [0]
        ipti_coverage = [0]

    df_array = np.array(
        list(itertools.product(annual_EIR, cm_coverage, seasonality,
                               rtss_target_group, smc_target_group, cm_target_group,
                               rtss_mode, rtss_coverage, smc_coverage, ipti_coverage, ipti_mode,rtss_booster1_min_age)))
    df = pd.DataFrame(df_array)
    # Optionally drop certain combinations
    # [...]

    # Prepare and save df
    df.columns = ['annual_EIR', 'cm_coverage', 'seasonality',
                  'rtss_target_group', 'smc_target_group', 'cm_target_group',
                  'rtss_mode', 'rtss_coverage', 'smc_coverage', 'ipti_coverage', 'ipti_mode','rtss_booster1_min_age']

    df.reset_index(inplace=True)
    df = df.rename(columns={'index': 'setting_id'})
    df['setting_id'] = 'HX' + df['setting_id'].astype(str)

    # for the correlated demographics file, use the minimum of the coverages as the fraction of people in the high-access group unless otherwise specified
    if type(high_access_frac) == int or type(high_access_frac) == float:
        df['frac_high_access'] = high_access_frac
    else:
        df['frac_high_access'] = [np.min([float(df['cm_coverage'][yy]),
                                          (float(df['rtss_coverage'][yy]) if float(df['rtss_coverage'][yy]) > 0 else 1),
                                          (float(df['smc_coverage'][yy]) if float(df['smc_coverage'][yy]) > 0 else 1),
                                          (float(df['ipti_coverage'][yy]) if float(df['ipti_coverage'][yy]) > 0 else 1)])
                                  for yy in range(len(df))]

    # add demographics, seasonality, and intervention csv filepaths
    df['demographics_filename'] = [os.path.join(base_scenario_filepath, 'Demographics',
                                   'generic_demographics_cohort_correlated%i.json' % int(float(df['frac_high_access'][yy]) * 100))
                                   for yy in range(len(df))]
    df.loc[df.frac_high_access == 0, 'demographics_filename'] = os.path.join(base_scenario_filepath, 'Demographics',
                                   'generic_demographics_cohort_uncorrelated.json')
    # df['seasonality_filename'] = os.path.join(base_scenario_filepath, 'Seasonality', 'seasonality_eir_multipliers')
    df['CM_filename'] = [os.path.join(base_scenario_filepath, 'CM', 'CM_constant_%icoverage.csv' % (100 * float(yy)))
                         for yy in df['cm_coverage']]
    df['RTSS_filename'] = [
        os.path.join(base_scenario_filepath, 'RTSS', f'RTSS_{xx}_{int(100 * float(yy))}coverage.csv') for xx, yy in
        zip(df['rtss_mode'], df['rtss_coverage'])]

    df['SMC_filename'] = [os.path.join(base_scenario_filepath, 'SMC', 'SMC_constant_%icoverage.csv' % (100 * float(yy)))
                         for yy in df['smc_coverage']]

    if smcipti == 'IPTi':
        df['IPTi_filename'] = [
            os.path.join(base_scenario_filepath, 'IPTi', f'IPTi_{xx}_{int(100 * float(yy))}coverage.csv') for xx, yy in
            zip(df['ipti_mode'], df['ipti_coverage'])]
        df['ipti_postprocess'] = param_dic['ipti_postprocess']

    print(f'Saving  {fname_out}_{smcipti}.csv')
    df.to_csv(os.path.join(projectpath, 'simulation_inputs/scenario_files/generic', f'{fname_out}_{smcipti}.csv'),
              index=False)


def intervention_inputs(param_dic, CM=False, RTSS=False, SMC=False, IPTi=False):
    if CM == True:
        cm_coverage = param_dic['cm_coverage']
        for cm in cm_coverage:
            if cm == 0:
                df = pd.DataFrame()
            else:
                df = pd.DataFrame({'U5_coverage': [cm], 'adult_coverage': [cm], 'severe_cases': np.min([0.8, cm * 2]),
                                   'simday': [0], 'duration': [-1], 'start_day_override': [-1], 'run_col': 'run'})
            print(f'Writing CM_constant_{int(100 * float(cm))}coverage.csv')
            df.to_csv(os.path.join(projectpath, 'simulation_inputs/scenario_files/generic', 'CM',
                                   'CM_constant_%icoverage.csv' % (100 * float(cm))))

    if RTSS == True:
        # RTS,S parameters
        for rtss_mode in param_dic['rtss_mode']:
            rtss_age_days = param_dic['rtss_age_days']
            if 'campboost' in rtss_mode:
                if rtss_mode == 'campboost':
                    num_rtss_rounds = 2  # main EPI round and 1 campaign booster
                else:
                    if 'campboost_flexible' in rtss_mode:
                        num_rtss_rounds = 4  # main EPI round and 1 booster within 3 years with defined coverages=
                        rtss_types = ['simple', 'booster1', 'booster2',
                                      'booster3']  # 1,2,3 for years to receive booster within 2-5 years
                    elif 'campboost2' in rtss_mode:
                        if 'campboost2_flexible' in rtss_mode:
                            num_rtss_rounds = 5  # main EPI round and 1 booster within 3 years with defined coverages, +1 booster in yr 4-5
                            rtss_types = ['simple', 'booster1', 'booster2', 'booster3', 'booster3']  # last deployment will get an property_restrictions 'GotBooster1'
                        else:
                            num_rtss_rounds = 3  # main EPI round and 3 boosters
                            rtss_types = ['simple', 'booster1', 'booster2']
                    elif 'campboost3' in rtss_mode:
                        num_rtss_rounds = 4  # main EPI round and 3 boosters
                        rtss_types = ['simple', 'booster1', 'booster2', 'booster3']
                    else:
                        num_rtss_rounds = 3  # main EPI round and 2 campaign boosters

                season_start_month = param_dic['smc_start_month']
                rtss_campaign_booster = (round(( season_start_month - 1) * 30.4)) - 7 + 365  # 1 week before SMC (adjusted for min age 24/36 months in set_up_interventions)

                ## seasonal campaign to be repeated every year ('campboost')
                rtss_booster_days = [rtss_campaign_booster + 365 * x for x in range(num_rtss_rounds - 1)]
                rtss_age_days = [rtss_age_days[0]] + rtss_booster_days  # initial EPI dose as well as dates of campaign booster doses

                try:
                    rtss_types
                except:
                    rtss_types = ['simple', 'booster1'] + ['booster2'] * (len(rtss_booster_days) - 1)
                effect_rtss = [0.8, 0.4] + [0.4] * (len(rtss_booster_days) - 1)
                deploy_types = ['EPI_cohort', 'campboost'] + ['campboost'] * (len(rtss_booster_days) - 1)
                tsteps_btwn_repetitions = [-1] * num_rtss_rounds
            elif 'noboost' in rtss_mode:
                num_rtss_rounds = 1
                rtss_types = ['simple']
                effect_rtss = [0.8]
                deploy_types = ['EPI_cohort']
                tsteps_btwn_repetitions = [-1]
            else:
                num_rtss_rounds = 2
                if len(rtss_age_days) != 2:
                    raise ValueError('two dates should be given for the EPI RTS,S deployment')
                ## Single deployment at defined ages (simulation days)
                rtss_types = ['simple', 'booster']
                effect_rtss = [0.8, 0.4]
                deploy_types = ['EPI_cohort'] * 2
                tsteps_btwn_repetitions = [-1] * 2

            booster_coverage_rtss = param_dic['booster_coverage_rtss']
            max_age_rtss = param_dic['max_age_rtss']
            decay_time_constant_rtss = param_dic['decay_time_constant_rtss']
            decay_class_rtss = param_dic['decay_class_rtss']

            # iterate through RTS,S scenarios, creating csvs for each (one row for initial dose, one row for booster)
            rtss_coverage = param_dic['rtss_coverage']
            for rtss in rtss_coverage:
                if rtss == 0:
                    df = pd.DataFrame()
                else:
                    if rtss_mode == 'campboost_flexible':
                        coverage_levels = [rtss] + [0.5] * 3  # corresponding to [0.5, 0.25, 0.125]
                    elif rtss_mode == 'campboost2_flexible':
                        coverage_levels = [rtss] + [0.5] * 3 + [booster_coverage_rtss]
                    else:
                        coverage_levels = [rtss] + [booster_coverage_rtss] * (num_rtss_rounds - 1)
                    df = pd.DataFrame(
                        {'coverage_levels': coverage_levels,
                         'round': list(range(1, num_rtss_rounds + 1)),
                         'rtss_types': rtss_types,
                         'initial_killing': effect_rtss,
                         'decay_time_constant': [decay_time_constant_rtss] * num_rtss_rounds,
                         'decay_class': [decay_class_rtss] * num_rtss_rounds,
                         'rtss_touchpoints': ['NA'] * num_rtss_rounds,
                         'RTSS_day': rtss_age_days,
                         'repetitions': [1] * num_rtss_rounds,
                         'tsteps_btwn_repetitions': tsteps_btwn_repetitions,
                         'deploy_type': deploy_types,
                         'agemin': [0] * num_rtss_rounds,
                         'agemax': [max_age_rtss] * num_rtss_rounds,
                         'distribution_name': ['CONSTANT_DISTRIBUTION'] * num_rtss_rounds,
                         'distribution_std': [1] * num_rtss_rounds,
                         'run_col': ['run'] * num_rtss_rounds})
                rtss_csv = f'RTSS_{rtss_mode}_{int(rtss * 100)}coverage.csv'

                """Add individual property restrictions (optional)"""
                if rtss > 0 and 'flexible' in rtss_mode:
                    df['rtss_property_restrictions'] = ''
                    df.at[len(df) - 1, 'rtss_property_restrictions'] = 'GotBooster1'
                print(f'Writing {rtss_csv}')
                df.to_csv(os.path.join(projectpath, 'simulation_inputs/scenario_files/generic', 'RTSS', rtss_csv))

    if SMC == True:
        # SMC parameters
        num_smc_rounds = param_dic['num_smc_rounds']
        smc_start_month = param_dic['smc_start_month']
        num_repeated_years = param_dic['num_repeated_years']

        # SMC - assumes the same coverage in all rounds, with the first of four rounds starting in the 7th month;
        smc_coverage = param_dic['smc_coverage']
        for smc in smc_coverage:
            if smc == 0:
                df_all_years = pd.DataFrame()
            else:
                df_r1 = pd.DataFrame({'round': [1], 'coverage': [smc], 'max_age': [5],
                                      'simday': [(smc_start_month - 1) * 30], 'duration': [-1], 'run_col': 'run'})
                df_base_year = df_r1.copy()
                for rr in range(num_smc_rounds - 1):
                    df_new_round = df_r1.copy()
                    df_new_round['round'] = rr + 2
                    df_new_round['simday'] = df_r1['simday'] + (1 + rr) * 30
                    df_base_year = df_base_year.append(df_new_round)

                # copy data frame for all of the included simulation years
                df_all_years = df_base_year.copy()
                for yy in range(num_repeated_years):
                    df_new_year = df_base_year.copy()
                    df_new_year['simday'] = df_base_year['simday'] + 365 * (yy + 1)
                    df_all_years = df_all_years.append(df_new_year)
            print(f'Writing SMC_constant_{int(100 * float(smc))}coverage.csv ')
            df_all_years.to_csv(os.path.join(projectpath, 'simulation_inputs/scenario_files/generic', 'SMC',
                                             'SMC_constant_%icoverage.csv' % (100 * float(smc))))

    if IPTi == True:
        # IPTi - assumes the same coverage in all rounds
        # Basic IPTi 2-3-9
        # extended IPTi 3-6-9-12-15-18
        basic_ipti = [2, 3, 9]
        extended3tp_ipti = [3, 6, 9]
        extended_ipti = [3, 6, 9, 12, 15, 18]

        for ipti_mode in param_dic['ipti_mode']:
            ipti_coverage = param_dic['ipti_coverage']

            if ipti_mode == 'basic' or ipti_mode == 'constant':
                ipti_mode = 'basic'
                ipti_touchpoints = [x * 365 / 12 for x in basic_ipti]
            elif ipti_mode == 'extended3tp':
                ipti_touchpoints = [x * 365 / 12 for x in extended3tp_ipti]
            elif ipti_mode == 'extended':
                ipti_touchpoints = [x * 365 / 12 for x in extended_ipti]
            else:
                ipti_mode == 'custom'
                ipti_touchpoints = param_dic['ipti_touchpoints']
            ntp = len(ipti_touchpoints)
            # EPI_cohort using campaign style deployment, distribution_name and distribution_std not used
            for ipti in ipti_coverage:
                if ipti == 0:
                    df = pd.DataFrame()
                else:
                    df = pd.DataFrame({'round': range(1, ntp + 1),
                                       'coverage_levels': [ipti] * ntp,
                                       'prop_no_ipti': [1 - ipti] * ntp,
                                       'ipti_touchpoints': ipti_touchpoints,
                                       'IPTi_day': ipti_touchpoints,
                                       'repetitions': [1] * ntp,
                                       'tsteps_btwn_repetitions': [-1] * ntp,
                                       'agemin': [0] * ntp,
                                       'agemax': [2] * ntp,
                                       'deploy_type': ['EPI_cohort'] * ntp,
                                       'drug_code': ['SP'] * ntp,
                                       'run_col': ['run'] * ntp})

                # fname = f'IPTi_{len(ipti_touchpoints)}tp_constant_{int(100 * ipti)}coverage.csv'
                fname = f'IPTi_{ipti_mode}_{int(100 * ipti)}coverage.csv'
                print(f'Writing {fname}')
                df.to_csv(os.path.join(projectpath, 'simulation_inputs/scenario_files/generic', 'IPTi', fname),
                          index=False)


def combine_input_files(scen_csv1, scen_csv2, combined_scen_csv):
    scen_df1 = pd.read_csv(os.path.join(projectpath, 'simulation_inputs/scenario_files/generic', scen_csv1),
                          encoding='latin')
    scen_df2 = pd.read_csv(os.path.join(projectpath, 'simulation_inputs/scenario_files/generic', scen_csv2),
                          encoding='latin')
    scen_df = pd.concat([scen_df1, scen_df2]).reset_index(drop=True)

    # re-label settings
    scen_df['setting_id'] = ['HX' + str(ii) for ii in list(range(scen_df.shape[0]))]
    scen_df.set_index(keys='setting_id', inplace=True)

    scen_df.to_csv(os.path.join(projectpath, 'simulation_inputs/scenario_files/generic', combined_scen_csv), encoding='latin')


def remove_duplicate_noRTSS_scenarios(scen_csv):
    # remove rows that are effectively duplicates: these are created from taking the full factorial of parameter sets -
    #    when an intervention coverage is zero, it is not necessary to have separate simulations for different access
    #    groups or schedule types of that intervention
    scen_df = pd.read_csv(os.path.join(projectpath, 'simulation_inputs/scenario_files/generic', scen_csv),
                          encoding='latin')
    # remove duplicate RTS,S rows when RTS,S coverage is zero
    rtss0_subset_cols = ['annual_EIR', 'seasonality', 'cm_coverage', 'ipti_coverage', 'smc_coverage',
                         'frac_high_access', 'cm_target_group', 'smc_target_group']  # columns to identify duplicates
    rtss0_subset_cols = list(set(rtss0_subset_cols).intersection(scen_df.columns))
    scen_df_rtss0 = scen_df[scen_df.rtss_coverage == 0]
    scen_df_rtss_remainder = scen_df[scen_df.rtss_coverage != 0]
    scen_df_rtss0_fixed = scen_df_rtss0[~scen_df_rtss0.duplicated(subset=rtss0_subset_cols, keep='first')]
    scen_df = pd.concat([scen_df_rtss0_fixed, scen_df_rtss_remainder])

    # remove duplicate SMC rows when SMC coverage is zero
    smc0_subset_cols = ['annual_EIR', 'seasonality', 'cm_coverage', 'ipti_coverage',
                         'rtss_coverage', 'frac_high_access', 'rtss_mode', 'minBoostAge', 'rtss_booster1_min_age',
                         'cm_target_group', 'rtss_target_group']  # columns to identify duplicates
    smc0_subset_cols = list(set(smc0_subset_cols).intersection(scen_df.columns))
    scen_df_smc0 = scen_df[scen_df.smc_coverage == 0]
    scen_df_smc_remainder = scen_df[scen_df.smc_coverage != 0]
    scen_df_smc0_fixed = scen_df_smc0[~scen_df_smc0.duplicated(subset=smc0_subset_cols, keep='first')]
    scen_df = pd.concat([scen_df_smc0_fixed, scen_df_smc_remainder])

    # re-sort rows
    scen_df = scen_df.sort_values(by=['annual_EIR', 'seasonality', 'cm_coverage', 'smc_coverage', 'rtss_coverage', 'rtss_target_group'])

    # re-label settings
    scen_df['setting_id'] = ['HX' + str(ii) for ii in list(range(scen_df.shape[0]))]
    scen_df.set_index(keys='setting_id', inplace=True)

    scen_df.to_csv(os.path.join(projectpath, 'simulation_inputs/scenario_files/generic',
                                scen_csv.replace('.csv', '_cleaned.csv')), encoding='latin')


if __name__ == "__main__":
    datapath, projectpath = load_box_paths()
    base_scenario_filepath = ''  # os.path.join('scenario_files', 'generic')
    p = get_parameter_space() # View parameter defaults and options
    print(p)

    """Scenario flags"""
    COUNTERFACTUAL = False
    COVERAGE_CORRELATION = False
    RTSS_BOOSTER = False
    RTSS_TYPES_SWEEP = False
    COVERAGE_HEATMAP = False
    ACCESS_CORRELATION = False
    ACCESS_CORRELATION_SMC_RTSS_ONLY = False
    ACCESS_CORRELATION_CM_RTSS_ONLY = False
    SEASONALITY_SWEEP = False
    SEASONALITY_SWEEP_REDUCED = True
    CM_SWEEP = False
    EIR_SWEEP = False
    EIR_CM_SWEEP = False
    EIR_CM_SWEEP_WITH_RTSS = False
    EIR_SEASON_SWEEP = False
    EIR_SMC_SWEEP = False
    IPTI_SWEEPS = False

    if COUNTERFACTUAL:
        """Coverage correlation scenarios"""
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0]})
        param_dic.update({'smc_coverage': [0, 0.6, 0.8]})
        param_dic.update({'cm_coverage': [0.6, 0.8]})
        param_dic.update({'annual_EIR': [1,5,10,30, 60]})  # constrain levels
        param_dic.update({'seasonality': ['constant','high_unimodal','moderate_unimodal']})
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_counterfactual_noRTSS')

    if COVERAGE_CORRELATION:
        """Coverage correlation scenarios"""
        param_dic = get_parameter_space()
        create_scenarios_mastercsv(param_dic, smcipti='RTSS')
        create_scenarios_mastercsv(param_dic, smcipti='SMC')
        create_scenarios_mastercsv(param_dic, smcipti='IPTi')
        intervention_inputs(param_dic, CM=True, RTSS=True, SMC=True, IPTi=True)

        """Lower level CM"""
        param_dic.update({'cm_coverage': [0.45]})
        intervention_inputs(param_dic, CM=True)
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_lowCM')

    if RTSS_BOOSTER:
        """RTSS booster scenarios______________________________________________________2
        2nd booster depending on seasonality, use same as SMC (1 day apart)
        Constrained full factorial
        """
        param_dic = get_parameter_space()
        param_dic.update({'rtss_mode': ['constant', 'campboost']})
        param_dic.update({'seasonality': ['high_unimodal']})
        param_dic.update({'annual_EIR': [5,10,30]})  # constrain levels

        create_scenarios_mastercsv(param_dic, smcipti='RTSS', fname_out='generic_campboost')
        # create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_campboost')
        # create_scenarios_mastercsv(param_dic, smcipti='IPTi', fname_out='generic_campboost')
        intervention_inputs(param_dic, RTSS=True)

        """RTSS alternative booster scenarios_____________________________________________2b
        single booster during 24-59 months at defined coverages per year 
        """
        param_dic = get_parameter_space()
        param_dic.update({'rtss_mode': ['campboost_flexible']})
        intervention_inputs(param_dic, RTSS=True)

        """
        if one booster is given before 48 months, give second booster
        """
        param_dic = get_parameter_space()
        param_dic.update({'rtss_mode': ['campboost2_flexible']})
        param_dic.update({'booster_coverage_rtss': 1})  # Assume 100% of those who got booster1 before 48 mths get booster2
        intervention_inputs(param_dic, RTSS=True)

        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0.8]})  # counterfactual already exists
        param_dic.update({'smc_coverage': [0, 0.8]})
        param_dic.update({'cm_coverage': [0.6]})
        param_dic.update({'rtss_mode': ['campboost_flexible', 'campboost2_flexible']})
        param_dic.update({'intervention_correlation': [0]})
        param_dic.update({'annual_EIR': [5,10, 30]})  # constrain levels
        param_dic.update({'seasonality': ['high_unimodal']})

        create_scenarios_mastercsv(param_dic, smcipti='RTSS', fname_out='generic_campboost_flexible_highseasonal')
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_campboost_flexible_highseasonal')

    if RTSS_TYPES_SWEEP:
        param_dic = get_parameter_space()
        param_dic.update({'smc_coverage': [0, 0.8]})
        param_dic.update({'cm_coverage': [0.6]})
        param_dic.update({'ipti_postprocess': False})
        param_dic.update({'ipti_coverage': [0]})
        param_dic.update({'seasonality': ['constant', 'high_unimodal']})
        param_dic.update({'annual_EIR': [5, 10, 30]})

        # EPI RTS,S
        param_dic.update({'rtss_coverage': [0, 0.6, 0.8]})
        param_dic.update({'rtss_mode': ['constant']})
        param_dic.update({'rtss_age_days': [274, 730]})
        intervention_inputs(param_dic, RTSS=True)
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_epi_rtss')

        ##--------------------------
        ## RTS,S booster csv's
        ##--------------------------
        param_dic = get_parameter_space()
        param_dic.update({'ipti_coverage': [0, 0.6, 0.8]})

        param_dic.update({'rtss_mode': ['noboost_7m']})
        param_dic.update({'rtss_age_days': [213]})
        intervention_inputs(param_dic, RTSS=True)

        param_dic.update({'rtss_mode': ['noboost']})
        param_dic.update({'rtss_age_days': [274]})
        intervention_inputs(param_dic, RTSS=True)

        param_dic.update({'rtss_mode': ['constant_7m']})
        param_dic.update({'rtss_age_days': [213, 730]})
        intervention_inputs(param_dic, RTSS=True)

        param_dic.update({'rtss_mode': ['campboost_7m']})
        param_dic.update({'rtss_age_days': [213, 730]})
        intervention_inputs(param_dic, RTSS=True)

        param_dic.update({'rtss_mode': ['campboost', 'campboost2']})
        param_dic.update({'rtss_age_days': [274, 730, 800]})  # dates of boosters will be overwritten
        intervention_inputs(param_dic, RTSS=True)

        # Enhanced RTS,S with 2 or 3 boosters and initial sequence 7 months, boosters starting at 9 months of age
        param_dic.update({'rtss_mode': ['campboost2_7m']})
        param_dic.update({'rtss_age_days': [213, 274, 639]})  # dates of boosters will be overwritten depending on seasonal campaign time
        intervention_inputs(param_dic, RTSS=True)

        param_dic.update({'rtss_mode': ['campboost3_7m']})
        param_dic.update({'rtss_age_days': [213, 274, 639, 1004]})  # dates of boosters will be overwritten depending on seasonal campaign time
        intervention_inputs(param_dic, RTSS=True)

        param_dic.update({'rtss_mode': ['campboost3']})
        param_dic.update({'rtss_age_days': [274, 730, 1095, 1460]})  # dates of boosters will be overwritten depending on seasonal campaign time
        intervention_inputs(param_dic, RTSS=True)

        ##--------------------------
        # Enhanced RTS,S master csv's
        ##--------------------------
        param_dic.update({'rtss_mode': ['campboost', 'campboost2']})
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_camp_rtss')

        ## Master csv, seasonal,
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0.8]})  # counterfactual already exists
        param_dic.update({'smc_coverage': [0, 0.8]})
        param_dic.update({'cm_coverage': [0.6]})
        param_dic.update({'ipti_coverage': [0]})
        param_dic.update({'seasonality': ['high_unimodal']})
        param_dic.update({'annual_EIR': [5,10,30]})

        ## default boosters
        param_dic.update({'rtss_booster1_min_age': [730]})
        param_dic.update({'rtss_mode': ['noboost','constant','campboost', 'campboost2', 'campboost3']})
        create_scenarios_mastercsv(param_dic, smcipti='RTSS', fname_out='generic_campboost_sweep_highseasonal')
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_campboost_sweep_highseasonal')

        ## early boosters
        param_dic.update({'rtss_booster1_min_age': [int(round(9*365/12,0))]})
        param_dic.update({'rtss_mode': ['noboost_7m','constant_7m','campboost_7m', 'campboost2_7m', 'campboost3_7m']})
        create_scenarios_mastercsv(param_dic, smcipti='RTSS', fname_out='generic_campboost_7m_sweep_highseasonal')
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_campboost_7m_sweep_highseasonal')

        ## Master csv, constant seasonality
        param_dic.update({'seasonality': ['constant']})

        ## default boosters
        param_dic.update({'rtss_booster1_min_age': [730]})
        param_dic.update({'rtss_mode': ['noboost','constant','campboost', 'campboost2', 'campboost3']})
        create_scenarios_mastercsv(param_dic, smcipti='RTSS', fname_out='generic_campboost_sweep_constant')
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_campboost_sweep_constant')

        ## early boosters
        param_dic.update({'rtss_booster1_min_age': [int(round(9 * 365 / 12, 0))]})
        param_dic.update({'rtss_mode': ['noboost_7m','constant_7m','campboost_7m', 'campboost2_7m', 'campboost3_7m']})
        create_scenarios_mastercsv(param_dic, smcipti='RTSS', fname_out='generic_campboost_7m_sweep_constant')
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_campboost_7m_sweep_constant')



    if COVERAGE_HEATMAP:
        """Coverage heatmap scenarios______________________________________________________3
        example for SMC only for IPTi can be simpler setup and handled in postprocessing for scaling method
        """
        param_dic = get_parameter_space()
        param_dic.update({'annual_EIR': [5,10,30]})
        param_dic.update({'rtss_coverage': np.arange(0, 1.1, 0.20)})
        param_dic.update({'smc_coverage': np.arange(0, 1.1, 0.20)})
        param_dic.update({'cm_coverage':  [0.6]})
        intervention_inputs(param_dic, RTSS=True)
        intervention_inputs(param_dic, SMC=True)

        access_correlation_types = ['random', 'high']
        for i in [0, 1]:
            param_dic.update({'rtss_target_group': [access_correlation_types[i]]})
            param_dic.update({'smc_target_group': [access_correlation_types[i]]})
            param_dic.update({'cm_target_group': [access_correlation_types[i]]})
            param_dic.update({'seasonality': ['constant']})
            create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out=f'generic_heatmap_constant_cor{i}')
            param_dic.update({'seasonality': ['high_unimodal']})
            create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out=f'generic_heatmap_highseason_cor{i}')

        ### With RTS,S seasonal booster
        param_dic = get_parameter_space()
        param_dic.update({'annual_EIR': [5, 10, 30]})
        param_dic.update({'rtss_coverage': np.arange(0, 1.1, 0.20)})
        param_dic.update({'smc_coverage': np.arange(0, 1.1, 0.20)})
        param_dic.update({'cm_coverage':  [0.6]})
        param_dic.update({'rtss_mode': ['campboost_flexible']})
        intervention_inputs(param_dic, RTSS=True)

        for i in [0, 1]:
            param_dic.update({'rtss_target_group': [access_correlation_types[i]]})
            param_dic.update({'smc_target_group': [access_correlation_types[i]]})
            param_dic.update({'cm_target_group': [access_correlation_types[i]]})
            param_dic.update({'seasonality': ['high_unimodal']})
            create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out=f'generic_heatmap_highseason_RTSSbooster_cor{i}')

        param_dic = get_parameter_space()
        rtss_modes = ['constant', 'campboost', 'campboost2', 'campboost3']
        for rtss_mode in  rtss_modes:
            param_dic.update({'rtss_mode': [rtss_mode]})
            param_dic.update({'rtss_coverage': np.arange(0, 1.1, 0.10)})
            intervention_inputs(param_dic, RTSS=True)
        param_dic.update({'smc_coverage': np.arange(0, 1.1, 0.10)})
        intervention_inputs(param_dic, SMC=True)

    if ACCESS_CORRELATION:
        """Sweep access correlation____________________________________________________"""
        # EPI RTS,S
        # With and without correlated access for SMC, RTS,S and CM; with and without anti-correlated access for RTS,S and SMC
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0,0.6, 0.8]})
        param_dic.update({'smc_coverage': [0,0.6, 0.8]})
        param_dic.update({'cm_coverage': [0.6,0.8]})
        param_dic.update({'rtss_target_group': ['random', 'high', 'low']})
        param_dic.update({'smc_target_group': ['random', 'high', 'low']})
        param_dic.update({'cm_target_group': ['random', 'high']})
        param_dic.update({'seasonality': ['high_unimodal']})
        param_dic.update({'annual_EIR': [5,10,30]})
        param_dic.update({'rtss_mode': ['constant']})
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out=f'generic_accesscorrelation')

        param_dic.update({'rtss_coverage': [0,0.6, 0.8]})
        param_dic.update({'ipti_coverage': [0,0.6, 0.8]})
        param_dic.update({'seasonality': ['constant']})
        param_dic.update({'ipti_postprocess': False})
        create_scenarios_mastercsv(param_dic, smcipti='IPTi', fname_out=f'generic_accesscorrelation')

        # The below access correlation and access anti-correlation are now included in the above sweep across correlation patterns
        # ### Access correlation without CM
        # param_dic = get_parameter_space()
        # param_dic.update({'rtss_coverage': [0,0.6, 0.8]})
        # param_dic.update({'smc_coverage': [0,0.6, 0.8]})
        # param_dic.update({'cm_coverage': [0.6,0.8]})
        # param_dic.update({'intervention_correlation': [1]})
        # param_dic.update({'seasonality': ['high_unimodal']})
        # param_dic.update({'annual_EIR': [10]})  # [5,10,30]
        # param_dic.update({'rtss_mode': ['constant']})
        # param_dic.update({'intervention_cor_cm': ['False']})
        # create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out=f'generic_accesscorrelation_independentCM')
        # #IPTi equivalent
        # param_dic.update({'ipti_postprocess': False})
        # param_dic.update({'seasonality': ['constant']})
        # create_scenarios_mastercsv(param_dic, smcipti='IPTi', fname_out=f'generic_accesscorrelation_independentCM')
        #
        # ### Anti-access correlation without CM
        # param_dic = get_parameter_space()
        # param_dic.update({'rtss_coverage': [0,0.6, 0.8]})
        # param_dic.update({'smc_coverage': [0,0.6, 0.8]})
        # param_dic.update({'cm_coverage': [0.6,0.8]})
        # param_dic.update({'intervention_correlation': [1]})
        # param_dic.update({'seasonality': ['high_unimodal']})
        # param_dic.update({'annual_EIR': [10]})  # [5,10,30]
        # param_dic.update({'rtss_mode': ['constant']})
        # param_dic.update({'intervention_cor_cm': ['False']})
        # param_dic.update({'antiaccess_high': ['SMC','RTSS']})
        # create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out=f'generic_antiaccesscorrelation_independentCM')
        # #IPTi equivalent
        # param_dic.update({'ipti_postprocess': False})
        # param_dic.update({'seasonality': ['constant']})
        # param_dic.update({'antiaccess_high': ['IPTi','RTSS']})
        # create_scenarios_mastercsv(param_dic, smcipti='IPTi', fname_out=f'generic_antiaccesscorrelation_independentCM')


    if ACCESS_CORRELATION_SMC_RTSS_ONLY:
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0, 0.8]})
        param_dic.update({'smc_coverage': [0, 0.8]})
        param_dic.update({'cm_coverage': [0.6]})
        param_dic.update({'rtss_target_group': ['random', 'high', 'low']})
        param_dic.update({'smc_target_group': ['high']})
        param_dic.update({'cm_target_group': ['random']})
        param_dic.update({'seasonality': ['high_unimodal']})
        param_dic.update({'annual_EIR': [5, 10, 30]})
        param_dic.update({'rtss_mode': ['constant']})
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out=f'generic_accesscorrelation_highSMC_randomCM', high_access_frac=0.5)

        # lower RTS,S and SMC coverages
        param_dic.update({'rtss_coverage': [0, 0.6]})
        param_dic.update({'smc_coverage': [0, 0.6]})
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out=f'generic_accesscorrelation_highSMC_randomCM_60covAll', high_access_frac=0.5)

        combine_input_files(scen_csv1='generic_accesscorrelation_highSMC_randomCM_SMC.csv',
                            scen_csv2='generic_accesscorrelation_highSMC_randomCM_60covAll_SMC.csv',
                            combined_scen_csv='generic_accesscorrelation_highSMC_randomCM.csv')
        remove_duplicate_noRTSS_scenarios(scen_csv='generic_accesscorrelation_highSMC_randomCM.csv')

    if ACCESS_CORRELATION_CM_RTSS_ONLY:
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0, 0.8]})
        param_dic.update({'smc_coverage': [0]})
        param_dic.update({'cm_coverage': [0.6]})
        param_dic.update({'rtss_target_group': ['random', 'high', 'low']})
        param_dic.update({'smc_target_group': ['random']})
        param_dic.update({'cm_target_group': ['high']})
        param_dic.update({'seasonality': ['high_unimodal']})
        param_dic.update({'annual_EIR': [5, 10, 30]})
        param_dic.update({'rtss_mode': ['constant']})
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out=f'generic_accesscorrelationCM05_v3', high_access_frac=0.5)

    if SEASONALITY_SWEEP:
        """Sweep seasonality"""
        # EPI RTS,S
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0, 0.8]})
        param_dic.update({'smc_coverage': [0, 0.8]})
        param_dic.update({'cm_coverage': [0.6]})
        param_dic.update({'seasonality': ['constant', 'moderate_unimodal', 'high_unimodal']})
        param_dic.update({'annual_EIR': [10]})

        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_season_sweep')

        # Enhanced RTS,S
        param_dic.update({'rtss_coverage': [0, 0.8]})
        param_dic.update({'rtss_mode': ['campboost', 'campboost2']})
        param_dic.update({'rtss_age_days': [274, 730, 800]})  # dates of boosters will be overwritten
        intervention_inputs(param_dic, RTSS=True)
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_season_enhanced_rtss')


    if SEASONALITY_SWEEP_REDUCED:
        """Sweep seasonality"""
        # EPI RTS,S
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0, 0.8]})
        param_dic.update({'smc_coverage': [0]})
        param_dic.update({'cm_coverage': [0.6]})
        param_dic.update({'seasonality': ['constant', 'moderate_unimodal', 'high_unimodal']})
        param_dic.update({'annual_EIR': [1, 5, 10, 30]})

        create_scenarios_mastercsv(param_dic, smcipti='RTSS', fname_out='generic_standard_season')
        remove_duplicate_noRTSS_scenarios(scen_csv='generic_standard_season_RTSS.csv')

        # Enhanced RTS,S
        param_dic.update({'rtss_coverage': [0.8]})
        param_dic.update({'rtss_mode': ['campboost']})
        param_dic.update({'rtss_age_days': [274, 730]})  # dates of boosters will be overwritten
        intervention_inputs(param_dic, RTSS=True)
        create_scenarios_mastercsv(param_dic, smcipti='RTSS', fname_out='generic_enhanced_season')
        remove_duplicate_noRTSS_scenarios(scen_csv='generic_enhanced_season_RTSS.csv')

    if CM_SWEEP:
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0, 0.8]})
        param_dic.update({'smc_coverage': [0, 0.8]})
        param_dic.update({'cm_coverage': [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]})
        param_dic.update({'seasonality': ['moderate_unimodal']})
        param_dic.update({'annual_EIR': [10]})

        intervention_inputs(param_dic, CM=True)
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_cm_sweep')

        param_dic.update({'rtss_target_group': ['high']})
        param_dic.update({'smc_target_group': ['high']})
        param_dic.update({'cm_target_group': ['high']})
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_cm_sweep_cor1')

        """Sweep CM (correlated interventions)______________________________________________________
        """
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0, 0.8]})
        param_dic.update({'smc_coverage': [0]})
        param_dic.update({'cm_coverage': [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]})
        param_dic.update({'rtss_target_group': ['high']})
        param_dic.update({'smc_target_group': ['high']})
        param_dic.update({'cm_target_group': ['high']})
        param_dic.update({'seasonality': ['moderate_unimodal']})
        param_dic.update({'annual_EIR': [10]})

        intervention_inputs(param_dic, CM=True)
        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_cm_correlated_sweep')

    if EIR_SWEEP:
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0, 0.6, 0.8]})
        param_dic.update({'smc_coverage': [0, 0.8]})
        param_dic.update({'cm_coverage': [0.6]})
        param_dic.update({'seasonality': ['moderate_unimodal']})
        param_dic.update({'annual_EIR': [1, 5, 10, 15, 20, 30, 60]})

        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_eir_sweep')

    if EIR_CM_SWEEP:
        param_dic = get_parameter_space()
        param_dic.update({'cm_coverage': [0.4, 0.5, 0.6, 0.8]})
        param_dic.update({'rtss_coverage': [0]})
        param_dic.update({'smc_coverage': [0]})
        param_dic.update({'cm_coverage': [0.4, 0.5, 0.6, 0.8]})
        param_dic.update({'annual_EIR': [1, 5, 10, 15, 20, 30, 60]})
        param_dic.update({'seasonality': ['constant', 'moderate_unimodal', 'high_unimodal']})

        create_scenarios_mastercsv(param_dic, smcipti='CM')
        intervention_inputs(param_dic, CM=True)
        create_scenarios_mastercsv(param_dic, smcipti='RTSS', fname_out='generic_eir_cm_sweep_no')

    if EIR_CM_SWEEP_WITH_RTSS:
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0, 0.8]})
        param_dic.update({'smc_coverage': [0]})
        param_dic.update({'cm_coverage': [0.3, 0.6, 0.9]})
        param_dic.update({'annual_EIR': [1, 5, 10, 15, 20, 30, 60]})
        param_dic.update({'seasonality': ['high_unimodal']})

        create_scenarios_mastercsv(param_dic, smcipti='CM')
        intervention_inputs(param_dic, CM=True)
        create_scenarios_mastercsv(param_dic, smcipti='RTSS', fname_out='generic_rtss_eir_cm_sweep2')

    if EIR_SEASON_SWEEP:
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0, 0.8]})
        param_dic.update({'smc_coverage': [0]})
        param_dic.update({'cm_coverage': [0.6]})
        param_dic.update({'seasonality': ['constant', 'moderate_unimodal', 'high_unimodal']})
        param_dic.update({'annual_EIR': [1, 5, 10, 15, 20, 30, 60]})

        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_eir_season_sweep2')

    if EIR_SMC_SWEEP:
        param_dic = get_parameter_space()
        param_dic.update({'rtss_coverage': [0, 0.8]})
        param_dic.update({'smc_coverage': [0, 0.5, 0.8]})
        param_dic.update({'cm_coverage': [0.6]})
        param_dic.update({'seasonality': ['high_unimodal']})
        param_dic.update({'annual_EIR': [1, 5, 10, 15, 20, 30, 60, 70]})

        create_scenarios_mastercsv(param_dic, smcipti='SMC', fname_out='generic_eir_smc_sweep2')

    if IPTI_SWEEPS:
        # Basic IPTi 2-3-9
        # extended IPTi 3-6-9-12-15-18
        param_dic = get_parameter_space()
        param_dic.update({'ipti_postprocess': False})
        param_dic.update({'ipti_coverage': [0, 0.6, 0.8, 1]})
        param_dic.update({'ipti_mode': ['basic', 'extended3tp', 'extended']})
        intervention_inputs(param_dic, IPTi=True)

        param_dic.update({'rtss_coverage': [0]})
        param_dic.update({'seasonality': ['constant']})
        create_scenarios_mastercsv(param_dic, smcipti='IPTi', fname_out='generic_noRTSS')

        param_dic = get_parameter_space()
        param_dic.update({'ipti_postprocess': False})
        param_dic.update({'seasonality': ['constant']})
        param_dic.update({'ipti_mode': ['basic', 'extended']})
        create_scenarios_mastercsv(param_dic, smcipti='IPTi',
                                   fname_out='generic_RTSS')  # IPTi will automatically be attached


