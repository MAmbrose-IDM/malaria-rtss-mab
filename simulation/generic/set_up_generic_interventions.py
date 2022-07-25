import pandas as pd
import numpy as np
import math as math

from dtk.interventions.irs import add_IRS
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.interventions.property_change import change_individual_property
from malaria.interventions.adherent_drug import configure_adherent_drug
from malaria.interventions.health_seeking import add_health_seeking
from malaria.interventions.malaria_diagnostic import add_diagnostic_survey
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign
from malaria.interventions.malaria_vaccine import add_vaccine
from malaria.reports.MalariaReport import add_event_counter_report


def calc_high_low_access_coverages(coverage_all, high_access_frac):
    if (high_access_frac < 1) & (coverage_all >= high_access_frac):
        coverage_high = 1
        coverage_low = (coverage_all - high_access_frac) / (1 - high_access_frac)
    else:
        coverage_high = coverage_all / high_access_frac
        coverage_low = 0
    return [coverage_high, coverage_low]



def get_concentration_at_time(self, tt, initial_concentration, fast_frac, k1, k2):
    concentration_at_tt = initial_concentration * (fast_frac * math.exp(-1 * tt / k1) + (1 - fast_frac) * math.exp(-1 * tt / k2))
    return concentration_at_tt


def get_time_efficacy_values(self, initial_concentration, max_efficacy, fast_frac, k1, k2, mm, total_time, booster_interval=365):
    # get concentration and efficacy through time for a single dose
    concentration_through_time = [get_concentration_at_time(tt, initial_concentration, fast_frac, k1, k2) for tt in range(total_time)]
    efficacy_through_time = max_efficacy * (1 - math.exp(mm * concentration_through_time))

    # get concentration and efficacy through time after a second dose. calculate the added efficacy through time compared to a single dose
    if booster_interval < total_time:
        after_booster_initial_concentration = concentration_through_time[booster_interval] + initial_concentration
    else:
        after_booster_initial_concentration = 0
    concentration_through_time_after_booster = [get_concentration_at_time(tt, after_booster_initial_concentration, fast_frac, k1, k2) for tt in range(total_time)]
    efficacy_through_time_after_booster = max_efficacy * (1 - math.exp(mm * concentration_through_time_after_booster))
    # calcuate how much _extra_ protection the booster provides over original efficacy (this will be added to existing efficacy from prior vaccine)
    protection_from_prior_dose = efficacy_through_time[booster_interval:] + [0] * (total_time - len(efficacy_through_time[booster_interval:]))
    added_booster_efficacy = efficacy_through_time_after_booster - protection_from_prior_dose
    return [[i for i in range(total_time)], efficacy_through_time, added_booster_efficacy]


class InterventionSuite:
    # hs
    hs_ds_col = 'repDS'
    hs_duration = None

    hs_start_col = 'simday'
    hs_coverage_age = {  # column: [agemin, agemax]
        'U5_coverage': [0, 5],
        'adult_coverage': [5, 100]
    }
    hs_severe_coverage_age = {
        'severe_cases': [0, 100]
    }
    hs_rates = 0.3
    hs_severe_rates = 0.5

    # itn
    itn_ds_col = 'repDS'
    itn_leak_factor = 0.9
    itn_cov_cols = ['U5_ITN_use', 'six_nine_ITN_use', 'ten_eighteen_ITN_use', 'over_eighteen_ITN_use']
    itn_cov_age_bin = [0, 5, 10, 18]
    itn_seasonal_values = [0.032, 0.032, 0.0378, 0.154, 0.177, 0.105, 0.25, 0.32, 0.23, 0.18, 0.032]
    itn_seasonal_months = [0, 32, 60, 91, 121, 152, 182, 213, 244, 274, 364]
    itn_retention_in_yr = 1.51
    itn_preg_max_months = 48
    itn_preg_monthly_births_per = 0.00001

    # smc
    smc_ds_col = 'DS_Name'
    smc_adherence = True
    smc_default_adherence = 0.8
    smc_adherence_multiplier = 1
    smc_sp_resist_day1_multiply = 1
    smc_drug_code = 'SDX_PYR'
    smc_coverage_col = 'coverage'
    smc_agemin = 0.25
    smc_agemax = 5
    smc_max_age_col = 'max_age'
    smc_TAT_col = 'TAT'
    smc_leakage = True
    smc_leak_agemax = 10
    smc_leak_coverage = 0.081

    # irs
    irs_ds_col = 'DS_Name'
    irs_start_col = 'IRS_day'
    irs_coverage_col = 'coverage'
    irs_box_dur_col = 'box_duration'
    irs_decay_t_col = 'decay_time_constant'
    irs_init_eff_col = 'initial_killing'

    # vacc
    vacc_ds_col = 'DS_name'

    # rtss
    rtss_ds_col = 'DS_name'
    rtss_type_col = 'rtss_types'
    rtss_start_col = 'RTSS_day'
    rtss_coverage_col = 'coverage_levels'
    rtss_touchpoint_col = 'rtss_touchpoints'  # days since births!
    rtss_deploy_type_col = 'deploy_type'
    rtss_distribution_col = 'distribution_name'
    rtss_std_col = 'distribution_std'
    rtss_min_age_col = 'agemin'
    rtss_max_age_col = 'agemax'
    rtss_repetitions = 'repetitions'
    rtss_tsteps_btwn_repetitions = 'tsteps_btwn_repetitions'
    rtss_init_eff_col = 'initial_killing'
    rtss_decay_t_col = 'decay_time_constant'
    rtss_decay_class_col = 'decay_class'
    rtss_property_restrictions = 'rtss_property_restrictions'

    # ipti
    ipti_ds_col = 'DS_Name'
    ipti_start_col = 'IPTi_day'
    ipti_coverage_col = 'coverage_levels'
    ipti_touchpoint_col = 'ipti_touchpoints'  # days since births!
    ipti_distribution_col = 'distribution_name'
    ipti_std_col = 'distribution_std'
    ipti_agemin = 'agemin'
    ipti_agemax = 'agemax'
    ipti_repetitions = 'repetitions'
    ipti_tsteps_btwn_repetitions = 'tsteps_btwn_repetitions'


    def add_pkpd_vacc(self, cb, vacc_df, my_ds, high_access_ip_frac=0, vacc_target_group='random', cohort_month_shift=0):
        vacc_df = vacc_df[vacc_df[self.vacc_ds_col].str.upper() == my_ds.upper()]

        """Use campaign-style deployment targeted to specific ages (since births disabled in simulation)"""
        for r, row in vacc_df.iterrows():
            # calculate campaign day, given cohort month shift
            start_day0 = row['simday']
            start_day = start_day0 - round(30.4 * cohort_month_shift)
            if start_day > 0:
                start_day = start_day + 365

            # TODO: add IP-dependent efficacy parameters in case of booster
            change_booster_IPs = False  # updated in booster statements if applicable
            if 'adaptive_booster_addition' in vacc_df.columns:
                change_booster_IPs = row['adaptive_booster_addition']  # change booster IPs when specifying coverage per booster dose

            """Set vaccine properties (e.g., initial efficacy, waning)"""
            try:
                waning_type = row['vacc_waning_type']
            except:
                waning_type = 'exponential'

            if waning_type == 'pkpd':
                initial_concentration = row['initial_concentration']
                max_efficacy = row['max_efficacy']
                fast_frac = row['fast_frac']
                k1 = row['k1']
                k2 = row['k2']
                mm = row['m']
                total_time = row['total_time']
                booster_adjust = row['booster_adjust']
                if booster_adjust:
                    booster_interval = row['booster_interval']
                else:
                    booster_interval = 365 * 100
                time_efficacy_values = get_time_efficacy_values(initial_concentration, max_efficacy, fast_frac, k1, k2, mm, total_time, booster_interval)
                time_efficacy_initial = time_efficacy_values[1][0]
                time_efficacy_multipliers = [time_efficacy_values[1][yy] / time_efficacy_initial for yy in range(len(time_efficacy_values[1]))]
                time_efficacy_boost_initial = time_efficacy_values[2][0]
                time_efficacy_boost_multipliers = [time_efficacy_values[2][yy] / time_efficacy_boost_initial for yy in range(len(time_efficacy_values[2]))]
                vaccine_params_no_boost = {"Waning_Config": {"Initial_Effect": time_efficacy_initial,
                                                             "Durability_Map":{
                                                                 "Times": time_efficacy_values[1],
                                                                 "Values": time_efficacy_multipliers
                                                             },
                                                             "Reference_Timer": 1,
                                                             "Expire_At_Durability_Map_End": 1,
                                                             "class": "WaningEffectMapLinear"}}
                vaccine_params_boost = {"Waning_Config": {"Initial_Effect": time_efficacy_boost_initial,
                                                        "Durability_Map":{
                                                            "Times": time_efficacy_values[1],
                                                            "Values": time_efficacy_boost_multipliers
                                                        },
                                                        "Reference_Timer": 1,
                                                        "Expire_At_Durability_Map_End": 1,
                                                        "class": "WaningEffectMapLinear"}}
            else:
                if waning_type != 'exponential':
                    raise ValueError("Unknown vaccine decay type, assuming exponential decay")
                vaccine_params_no_boost = {"Waning_Config": {"Initial_Effect": row['vacc_init_eff'],
                                                    "Decay_Time_Constant": row['vacc_decay_t'],
                                                    "class": "WaningEffectExponential"}}
                vaccine_params_boost = vaccine_params_no_boost


            """Set group of individuals to receive vaccine and add campaign"""
            if high_access_ip_frac > 0 and vacc_target_group in ['low', 'high']:  # different coverages in different access groups
                if vacc_target_group == 'low':
                    # vaccine is preferentially given to the 'low-access' group
                    high_low_coverages = calc_high_low_access_coverages(coverage_all=row[self.rtss_coverage_col],
                                                                        high_access_frac=1 - high_access_ip_frac)
                    high_low_coverages = [high_low_coverages[1], high_low_coverages[0]]
                elif vacc_target_group == 'high':
                    high_low_coverages = calc_high_low_access_coverages(coverage_all=row[self.rtss_coverage_col],
                                                                        high_access_frac=high_access_ip_frac)

                # high-access coverage
                # TODO: add IPs for whether this should be treated as a booster
                add_vaccine(cb,
                            vaccine_type='RTSS',
                            vaccine_params=vaccine_params_no_boost,
                            start_days=[start_day],
                            coverage=high_low_coverages[0],
                            repetitions=row[self.rtss_repetitions],
                            tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                            target_group={'agemin': row[self.rtss_min_age_col],
                                          'agemax': row[self.rtss_max_age_col]},
                            ind_property_restrictions=[{'AccessToInterventions': 'higher'},
                                                       {'vaccine_selected': 'True'}],
                            disqualifying_properties=[{"remove_vaccine": "True"}])

                # low-access coverage
                add_vaccine(cb,
                            vaccine_type='RTSS',
                            vaccine_params=vaccine_params_no_boost,
                            start_days=[start_day],
                            coverage=high_low_coverages[1],
                            repetitions=row[self.rtss_repetitions],
                            tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                            target_group={'agemin': row[self.rtss_min_age_col],
                                          'agemax': row[self.rtss_max_age_col]},
                            ind_property_restrictions=[{'AccessToInterventions': 'lower'},
                                                       {'vaccine_selected': 'True'}],
                            disqualifying_properties=[{"remove_vaccine": "True"}])

            else:  # uniform probability of getting the vaccine
                if vacc_target_group != 'random':
                    print('WARNING: name for RTS,S access-group targeting not recognized, assuming random access.')
                add_vaccine(cb,
                            vaccine_type='RTSS',
                            vaccine_params=vaccine_params_no_boost,
                            start_days=[start_day],
                            coverage=row[self.rtss_coverage_col],
                            repetitions=row[self.rtss_repetitions],
                            tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                            target_group={'agemin': row[self.rtss_min_age_col],
                                          'agemax': row[self.rtss_max_age_col]},
                            ind_property_restrictions=[{'vaccine_selected': 'True'}],
                            disqualifying_properties=[{"remove_vaccine": "True"}])

            # self.change_rtss_ips(cb, change_booster_IPs=change_booster_IPs)
            cb.update_params({
                "Report_Event_Recorder_Events": ['Births', 'PropertyChange', 'Received_Vaccine', 'Received_Treatment']
            })

        return len(vacc_df)


    def add_ds_hs(self, cb, hs_df, my_ds, high_access_ip_frac=0, cm_target_group='random'):
        ds_col = self.hs_ds_col
        duration = self.hs_duration
        df = hs_df[hs_df[ds_col] == my_ds]
        for r, row in df.iterrows():
            self.add_hs_from_file(cb, row, duration=duration, high_access_ip_frac=high_access_ip_frac, cm_target_group=cm_target_group)

        return len(df)

    def add_hs_from_file(self, cb, row, duration, high_access_ip_frac=0, cm_target_group='random'):
        rates = self.hs_rates
        severe_rates = self.hs_severe_rates

        start_day = row[self.hs_start_col]  # if start_day_override < 0 else start_day_override
        if duration is None:
            duration = row['duration']

        # Uncomplicated - access depends on IP
        if (high_access_ip_frac > 0) and (cm_target_group != 'random'):
            for key, value in self.hs_coverage_age.items():
                if cm_target_group == 'low':  # CM is preferentially given to the 'low-access' group
                    high_low_coverages = calc_high_low_access_coverages(coverage_all=row[key],
                                                                        high_access_frac=1 - high_access_ip_frac)
                    high_low_coverages = [high_low_coverages[1], high_low_coverages[0]]
                elif cm_target_group == 'high':
                    high_low_coverages = calc_high_low_access_coverages(coverage_all=row[key],
                                                                        high_access_frac=high_access_ip_frac)
                else:
                    print('WARNING: name for CM access-group targeting not recognized.')

                # high-access coverage
                targets = [{'trigger': 'NewClinicalCase',
                            'coverage': high_low_coverages[0],
                            'agemin': value[0],
                            'agemax': value[1],
                            'seek': 1,
                            'rate': rates
                            }]
                add_health_seeking(cb, start_day=start_day,
                                   targets=targets,
                                   drug=['Artemether', 'Lumefantrine'], duration=duration,
                                   ind_property_restrictions=[{'AccessToInterventions': 'higher'}])
                # low-access coverage
                targets = [{'trigger': 'NewClinicalCase',
                            'coverage': high_low_coverages[1],
                            'agemin': value[0],
                            'agemax': value[1],
                            'seek': 1,
                            'rate': rates
                            }]
                add_health_seeking(cb, start_day=start_day,
                                   targets=targets,
                                   drug=['Artemether', 'Lumefantrine'], duration=duration,
                                   ind_property_restrictions=[{'AccessToInterventions': 'lower'}])
        else:
            targets = []
            for key, value in self.hs_coverage_age.items():
                targets.append({
                    'trigger': 'NewClinicalCase',
                    'coverage': row[key],
                    'agemin': value[0],
                    'agemax': value[1],
                    'seek': 1,
                    'rate': rates
                })
            add_health_seeking(cb, start_day=start_day,
                               targets=targets,
                               drug=['Artemether', 'Lumefantrine'], duration=duration)

        # Severe - same coverage for all access-IPS
        targets = []
        for key, value in self.hs_severe_coverage_age.items():
            targets.append({
                'trigger': 'NewSevereCase',
                'coverage': row[key],
                'agemin': value[0],
                'agemax': value[1],
                'seek': 1,
                'rate': severe_rates
            })
        add_health_seeking(cb, start_day=start_day,
                           targets=targets,
                           drug=['Artemether', 'Lumefantrine'], duration=duration,
                           broadcast_event_name='Received_Severe_Treatment')

    def add_ds_smc(self, cb, smc_df, my_ds, high_access_ip_frac=0, smc_target_group='random', cohort_month_shift=0):

        ds_col = self.smc_ds_col
        adherence_multiplier = self.smc_adherence_multiplier
        sp_resist_day1_multiply = self.smc_sp_resist_day1_multiply
        df = smc_df[smc_df[ds_col] == my_ds]
        drug_code = self.smc_drug_code
        if self.smc_adherence:
            if 'adherence' in smc_df.columns.values:
                adherent_drug_configs = self.smc_adherent_configuration(cb=cb,
                                                                        adherence=df['adherence'].values[
                                                                                      0] * adherence_multiplier,
                                                                        sp_resist_day1_multiply=sp_resist_day1_multiply)
            else:
                default_adherence = self.smc_default_adherence
                adherent_drug_configs = self.smc_adherent_configuration(cb=cb,
                                                                        adherence=default_adherence * adherence_multiplier,
                                                                        sp_resist_day1_multiply=sp_resist_day1_multiply)
            drug_code = None
            adherent_drug_configs = [adherent_drug_configs]

        else:
            adherent_drug_configs = None

        def det_agemax(type, dfmax, argmax):
            if type == 'fixed':
                return (argmax)
            elif type == 'df':
                return (dfmax)

        agemin = self.smc_agemin
        agemax = self.smc_agemax
        for r, row in df.iterrows():
            if self.smc_max_age_col in smc_df.columns.values:
                max_smc_age = row[self.smc_max_age_col]
            else:
                max_smc_age = 5
            if self.smc_TAT_col in smc_df.columns.values:
                TAT = row[self.smc_TAT_col]
            else:
                TAT = 0  # assume kids with fever don't get SMC at all and are also not referred

            # calculate SMC day, given cohort month shift
            start_day0 = row['simday']
            start_day = start_day0 - round(30.4 * cohort_month_shift)
            if start_day > 0:
                if high_access_ip_frac > 0 and smc_target_group != 'random':  # different coverages in different access groups
                    if smc_target_group == 'low':
                        # SMC is preferentially given to the 'low-access' group
                        high_low_coverages = calc_high_low_access_coverages(coverage_all=row[self.smc_coverage_col],
                                                                            high_access_frac=1-high_access_ip_frac)
                        high_low_coverages = [high_low_coverages[1], high_low_coverages[0]]
                    elif smc_target_group == 'high':
                        high_low_coverages = calc_high_low_access_coverages(coverage_all=row[self.smc_coverage_col],
                                                                            high_access_frac=high_access_ip_frac)
                    else:
                        print('WARNING: name for SMC access-group targeting not recognized.')
                    # higher access
                    add_drug_campaign(cb, 'SMC', drug_code, start_days=[start_day],
                                      coverage=high_low_coverages[0],
                                      target_group={'agemin': agemin,
                                                    'agemax': agemax},
                                      listening_duration=2,
                                      trigger_condition_list=['No_SMC_Fever'],
                                      ind_property_restrictions=[{'AccessToInterventions': 'higher'}],
                                      adherent_drug_configs=adherent_drug_configs.copy())
                    # lower access
                    add_drug_campaign(cb, 'SMC', drug_code, start_days=[start_day],
                                      coverage=high_low_coverages[1],
                                      target_group={'agemin': agemin,
                                                    'agemax': agemax},
                                      listening_duration=2,
                                      trigger_condition_list=['No_SMC_Fever'],
                                      ind_property_restrictions=[{'AccessToInterventions': 'lower'}],
                                      adherent_drug_configs=adherent_drug_configs.copy())
                    if TAT:
                        # higher access
                        add_drug_campaign(cb, 'MDA', drug_code='AL', start_days=[start_day],
                                          coverage=high_low_coverages[0],
                                          target_group={'agemin': agemin,
                                                        'agemax': agemax},
                                          listening_duration=2,
                                          trigger_condition_list=['Has_SMC_Fever'],
                                          ind_property_restrictions=[{'AccessToInterventions': 'higher'}],
                                          receiving_drugs_event_name='Received_TAT_Treatment')
                        # lower access
                        add_drug_campaign(cb, 'MDA', drug_code='AL', start_days=[start_day],
                                          coverage=high_low_coverages[1],
                                          target_group={'agemin': agemin,
                                                        'agemax': agemax},
                                          listening_duration=2,
                                          trigger_condition_list=['Has_SMC_Fever'],
                                          ind_property_restrictions=[{'AccessToInterventions': 'lower'}],
                                          receiving_drugs_event_name='Received_TAT_Treatment')
                else:  # uniform access
                    add_drug_campaign(cb, 'SMC', drug_code, start_days=[start_day],
                                      coverage=row[self.smc_coverage_col],
                                      target_group={'agemin': agemin,
                                                    'agemax': agemax},
                                      listening_duration=2,
                                      trigger_condition_list=['No_SMC_Fever'],
                                      adherent_drug_configs=adherent_drug_configs.copy())
                    if TAT:
                        add_drug_campaign(cb, 'MDA', drug_code='AL', start_days=[start_day],
                                          coverage=row[self.smc_coverage_col],
                                          target_group={'agemin': agemin,
                                                        'agemax': agemax},
                                          listening_duration=2,
                                          trigger_condition_list=['Has_SMC_Fever'],
                                          receiving_drugs_event_name='Received_TAT_Treatment')

                add_diagnostic_survey(cb, start_day=start_day,
                                      coverage=1,
                                      target={"agemin": agemin,
                                              "agemax": agemax},
                                      diagnostic_type='FEVER',
                                      diagnostic_threshold=0.5,
                                      negative_diagnosis_configs=[{
                                          "Broadcast_Event": "No_SMC_Fever",
                                          "class": "BroadcastEvent"}],
                                      positive_diagnosis_configs=[{
                                          "Broadcast_Event": "Has_SMC_Fever",
                                          "class": "BroadcastEvent"}]
                                      )

        if self.smc_leakage and len(df) > 0:
            # leakage to between age max of SMC group and user defined max
            leak_agemin = agemax
            # adjust start days to account for cohort birth month
            start_days0 = df['simday'].values - round(30.4 * cohort_month_shift)
            start_days = start_days0[start_days0 > 0]
            add_drug_campaign(cb, 'SMC', drug_code, start_days=start_days,
                              coverage=self.smc_leak_coverage,
                              target_group={'agemin': leak_agemin, 'agemax': self.smc_leak_agemax},
                              adherent_drug_configs=adherent_drug_configs.copy())

        return len(df)

    def smc_adherent_configuration(self, cb, adherence, sp_resist_day1_multiply):

        smc_adherent_config = configure_adherent_drug(cb,
                                                      doses=[["SulfadoxinePyrimethamine", 'Amodiaquine'],
                                                             ['Amodiaquine'],
                                                             ['Amodiaquine']],
                                                      dose_interval=1,
                                                      non_adherence_options=['Stop'],
                                                      non_adherence_distribution=[1],
                                                      adherence_config={
                                                          "class": "WaningEffectMapCount",
                                                          "Initial_Effect": 1,
                                                          "Durability_Map": {
                                                              "Times": [
                                                                  1.0,
                                                                  2.0,
                                                                  3.0
                                                              ],
                                                              "Values": [
                                                                  sp_resist_day1_multiply,  # for day 1
                                                                  adherence,  # day 2
                                                                  adherence  # day 3
                                                              ]
                                                          }
                                                      }
                                                      )
        return smc_adherent_config

    def change_rtss_ips(self, cb, change_booster_IPs=False):
        change_individual_property(cb,
                                   target_property_name='VaccineStatus',
                                   target_property_value='GotVaccine',
                                   ind_property_restrictions=[{'VaccineStatus': 'None'}],
                                   trigger_condition_list=['Received_Vaccine'],
                                   blackout_flag=False)
        # UPDATE 2021-09-15: Currently, the generic model is set up to give Booster1 in the first year and Booster2 in the second year (this differs from the Nigeria setup).
        #     Booster2 is given to anyone who received the first vaccine, regardless of whether they also received Booster1.
        if change_booster_IPs:
            change_individual_property(cb,
                                       target_property_name='VaccineStatus',
                                       target_property_value='GotBooster1',
                                       ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                                       trigger_condition_list=['Received_Vaccine'],
                                       blackout_flag=False)
            change_individual_property(cb,
                                       target_property_name='VaccineStatus',
                                       target_property_value='GotBooster2',
                                       ind_property_restrictions=[{'VaccineStatus': 'GotBooster1'}],
                                       trigger_condition_list=['Received_Vaccine'],
                                       blackout_flag=False)

    def add_ds_rtss(self, cb, rtss_df, my_ds, high_access_ip_frac=0, rtss_target_group='random',rtss_booster1_min_age=730, cohort_month_shift=0):
        rtss_df = rtss_df[rtss_df[self.rtss_ds_col].str.upper() == my_ds.upper()]
        change_booster_IPs = False # updated in booster statements if applicable
        
        rtss_booster2_min_age = rtss_booster1_min_age + 365
        rtss_booster3_min_age = rtss_booster2_min_age + 365
    
        if 'rtss_property_restrictions' in rtss_df.columns:
            change_booster_IPs = True # change booster IPs when specifying coverage per booster dose
        """Use campaign-style deployment targeted to specific ages (since births disabled in simulation)"""
        for r, row in rtss_df.iterrows():
            try:
                Waning_Class = rtss_df['rtss_decay_class_col'].unique()[0]
            except:
                Waning_Class = "WaningEffectExponential"

            vaccine_params = {"Waning_Config": {"Initial_Effect": row[self.rtss_init_eff_col],
                                                "Decay_Time_Constant": row[self.rtss_decay_t_col],
                                                "class": Waning_Class}}

            vtype = row[self.rtss_type_col]
            if vtype == 'booster' or vtype =='booster1':
                """everyone who received the original vaccine has the same probability of getting a booster, regardless of IP 
                (otherwise, there will be different booster coverages depending on the fraction of people in the high IP)
                if campboost, receive booster between 24-35 month
                """
                if row[self.rtss_deploy_type_col] == 'EPI_cohort':
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[row[self.rtss_start_col]],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                                target_group={'agemin': row[self.rtss_min_age_col], 'agemax': row[self.rtss_max_age_col]})
                elif row[self.rtss_deploy_type_col] == 'campboost':
                    # calculate RTS,S booster day, given cohort month shift
                    start_day0 = row[self.rtss_start_col]
                    start_day = start_day0 - round(30.4 * cohort_month_shift)
                    # if booster would occur before the eligible age, wait until the next year
                    while start_day < rtss_booster1_min_age:
                        start_day = start_day + 365
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[start_day],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                                target_group={'agemin': row[self.rtss_min_age_col], 'agemax': row[self.rtss_max_age_col]})

            elif vtype == 'booster2':
                """everyone who received the original vaccine has the same probability of getting a booster, regardless of IP 
                (otherwise, there will be different booster coverages depending on the fraction of people in the high IP)
                if campboost, receive booster between 36-47 month
                """
                try:
                    booster_restr = row[self.rtss_property_restrictions]
                    if booster_restr == 'GotBooster1':
                        change_booster_IPs = True
                    if len(booster_restr) < 1:
                        booster_restr = 'GotVaccine'
                except:
                    booster_restr = 'GotVaccine'
                if row[self.rtss_deploy_type_col] == 'EPI_cohort':
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[row[self.rtss_start_col]],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                ind_property_restrictions=[{'VaccineStatus': booster_restr}],  # even if someone didn't get booster1 during the 24-month EPI visit, they can still get a booster during an EPI visit the following year
                                target_group={'agemin': row[self.rtss_min_age_col], 'agemax': row[self.rtss_max_age_col]})
                elif row[self.rtss_deploy_type_col] == 'campboost':
                    # calculate RTS,S booster day, given cohort month shift
                    start_day0 = row[self.rtss_start_col]
                    start_day = start_day0 - round(30.4 * cohort_month_shift)
                    # if booster would occur before the eligible age, wait until the next year
                    while start_day < rtss_booster2_min_age:
                        start_day = start_day + 365
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[start_day],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                ind_property_restrictions=[{'VaccineStatus': booster_restr}],  # even if someone didn't get booster1 during the first campaign, they can still get a booster during the following campaign
                                target_group={'agemin': row[self.rtss_min_age_col], 'agemax': row[self.rtss_max_age_col]})
            elif vtype == 'booster3':
                """ if campboost, receive booster between 48-59 month
                """
                try:
                    booster_restr = row[self.rtss_property_restrictions]
                    if len(booster_restr)<1:
                        booster_restr = 'GotVaccine'
                except:
                    booster_restr = 'GotVaccine'
                if row[self.rtss_deploy_type_col] == 'EPI_cohort':
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[row[self.rtss_start_col]],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                ind_property_restrictions=[{'VaccineStatus': booster_restr}],
                                target_group={'agemin': row[self.rtss_min_age_col], 'agemax': row[self.rtss_max_age_col]})
                elif row[self.rtss_deploy_type_col] == 'campboost':
                    # calculate RTS,S booster day, given cohort month shift
                    start_day0 = row[self.rtss_start_col]
                    start_day = start_day0 - round(30.4 * cohort_month_shift)
                    # if booster would occur before the eligible age, wait until the next year
                    while start_day < rtss_booster3_min_age:
                        start_day = start_day + 365
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[start_day],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                ind_property_restrictions=[{'VaccineStatus': booster_restr}],  # even if someone didn't get booster1 during the previous campaigns, they can still get a booster during the following campaign
                                target_group={'agemin': row[self.rtss_min_age_col], 'agemax': row[self.rtss_max_age_col]})
            else:
                if high_access_ip_frac > 0 and rtss_target_group != 'random':  # different coverages in different access groups
                    if rtss_target_group == 'low':
                        # RTSS is preferentially given to the 'low-access' group
                        high_low_coverages = calc_high_low_access_coverages(coverage_all=row[self.rtss_coverage_col],
                                                                            high_access_frac=1 - high_access_ip_frac)
                        high_low_coverages = [high_low_coverages[1], high_low_coverages[0]]
                    elif rtss_target_group == 'high':
                        high_low_coverages = calc_high_low_access_coverages(coverage_all=row[self.rtss_coverage_col],
                                                                            high_access_frac=high_access_ip_frac)
                    else:
                        print('WARNING: name for RTS,S access-group targeting not recognized.')
                    # high-access coverage
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[row[self.rtss_start_col]],
                                coverage=high_low_coverages[0],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                target_group={'agemin': row[self.rtss_min_age_col],
                                              'agemax': row[self.rtss_max_age_col]},
                                ind_property_restrictions=[{'AccessToInterventions': 'higher'}])

                    # low-access coverage
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[row[self.rtss_start_col]],
                                coverage=high_low_coverages[1],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                target_group={'agemin': row[self.rtss_min_age_col],
                                              'agemax': row[self.rtss_max_age_col]},
                                ind_property_restrictions=[{'AccessToInterventions': 'lower'}])

                else:  # uniform probability of getting the vaccine
                    add_vaccine(cb,
                                vaccine_type='RTSS',
                                vaccine_params=vaccine_params,
                                start_days=[row[self.rtss_start_col]],
                                coverage=row[self.rtss_coverage_col],
                                repetitions=row[self.rtss_repetitions],
                                tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                                target_group={'agemin': row[self.rtss_min_age_col],
                                              'agemax': row[self.rtss_max_age_col]})

            self.change_rtss_ips(cb,change_booster_IPs=change_booster_IPs)
            cb.update_params({
                "Report_Event_Recorder_Events": ['Births', 'PropertyChange', 'Received_Vaccine', 'Received_Treatment']
            })

        return len(rtss_df)



    def add_ds_ipti(self, cb, ipti_df, my_ds, high_access_ip_frac=0):
        ipti_df = ipti_df[ipti_df[self.ipti_ds_col].str.upper() == my_ds.upper()]
        EPI = ipti_df['deploy_type'].unique()[0] == 'EPI'
        try:
            drug_code = ipti_df['drug_code'].unique()[0]
        except:
            drug_code = 'SDX_PYR' # 'SP'
        if not EPI:
            """Use campaign-style deployment"""
            for r, row in ipti_df.iterrows():
                if high_access_ip_frac > 0:  # different coverages in different access groups
                    high_low_coverages = calc_high_low_access_coverages(coverage_all=row[self.ipti_coverage_col],
                                                                        high_access_frac=high_access_ip_frac)
                    # higher access
                    add_drug_campaign(cb, campaign_type='SMC',
                                      drug_code=drug_code,
                                      coverage=high_low_coverages[0],
                                      start_days=[row[self.ipti_start_col]],
                                      target_group={'agemin': row[self.ipti_agemin], 'agemax': row[self.ipti_agemax]},
                                      repetitions=row[self.ipti_repetitions],
                                      tsteps_btwn_repetitions=row[self.ipti_tsteps_btwn_repetitions],
                                      ind_property_restrictions=[{'AccessToInterventions': 'higher'}],
                                      receiving_drugs_event_name='Received_Campaign_Drugs')

                    # lower access
                    add_drug_campaign(cb, campaign_type='SMC',
                                      drug_code=drug_code,
                                      coverage=high_low_coverages[1],
                                      start_days=[row[self.ipti_start_col]],
                                      target_group={'agemin': row[self.ipti_agemin], 'agemax': row[self.ipti_agemax]},
                                      repetitions=row[self.ipti_repetitions],
                                      tsteps_btwn_repetitions=row[self.ipti_tsteps_btwn_repetitions],
                                      ind_property_restrictions=[{'AccessToInterventions': 'lower'}],
                                      receiving_drugs_event_name='Received_Campaign_Drugs')
                else:
                    add_drug_campaign(cb, campaign_type='SMC',
                                      drug_code=drug_code,
                                      coverage=row[self.ipti_coverage_col],
                                      start_days=[row[self.ipti_start_col]],
                                      target_group={'agemin': row[self.ipti_agemin], 'agemax': row[self.ipti_agemax]},
                                      repetitions=row[self.ipti_repetitions],
                                      tsteps_btwn_repetitions=row[self.ipti_tsteps_btwn_repetitions],
                                      receiving_drugs_event_name='Received_Campaign_Drugs')

        else:
            """Use birthtriggered deployment for IPTi administered along EPI"""
            coverage_levels = list(ipti_df[self.ipti_coverage_col].values)
            ipti_touchpoints = list(ipti_df[self.ipti_touchpoint_col].values)
            start_days = list(ipti_df[self.ipti_start_col].unique())
            delay_distribution_name = list(ipti_df[self.ipti_distribution_col].values)[0]
            Std_Dev_list = list(ipti_df[self.ipti_std_col].values)
            ipti_event_names = [f'IPTi_{x + 1}' for x in range(len(ipti_touchpoints))]

            for tp_time_trigger, coverage, event_name, std in zip(ipti_touchpoints, coverage_levels,
                                                                  ipti_event_names, Std_Dev_list):
                if delay_distribution_name == "GAUSSIAN_DISTRIBUTION":
                    delay_distribution = {"Delay_Period_Distribution": "GAUSSIAN_DISTRIBUTION",
                                          "Delay_Period_Gaussian_Mean": tp_time_trigger,
                                          "Delay_Period_Gaussian_Std_Dev": std}
                    tp_time_trigger = None
                else:
                    delay_distribution = None

                add_drug_campaign(cb,
                                  campaign_type='IPTi',
                                  drug_code=drug_code,
                                  start_days=start_days,
                                  coverage=coverage,
                                  delay_distribution=delay_distribution,
                                  triggered_campaign_delay=tp_time_trigger,
                                  trigger_name=event_name)

        return len(ipti_df)


def add_all_interventions(cb, int_suite, my_ds, high_access_ip_frac=0,
                          rtss_target_group='random', smc_target_group='random', cm_target_group='random',
                          rtss_booster1_min_age=730,
                          hs_df=pd.DataFrame(),
                          smc_df=pd.DataFrame(),
                          rtss_df=pd.DataFrame(),
                          ipti_df=pd.DataFrame(),
                          addtl_smc_func=None,
                          cohort_month_shift=0):
    event_list = ['Received_NMF_Treatment']

    if not smc_df.empty:
        if addtl_smc_func:
            addtl_smc = addtl_smc_func(cb, smc_df, my_ds)
        has_smc = int_suite.add_ds_smc(cb, smc_df, my_ds, high_access_ip_frac, smc_target_group, cohort_month_shift)
        if has_smc > 0:
            event_list.append('Received_Campaign_Drugs')

    if not hs_df.empty:
        """per default CM is included in intervention access correlation, otherwise set high_access_ip_frac_cm to 0"""
        has_cm = int_suite.add_ds_hs(cb, hs_df, my_ds, high_access_ip_frac, cm_target_group)
        if has_cm:
            event_list.append('Received_Treatment')
            event_list.append('Received_Severe_Treatment')

    if not rtss_df.empty:
        has_rtss = int_suite.add_ds_rtss(cb, rtss_df, my_ds, high_access_ip_frac, rtss_target_group,rtss_booster1_min_age, cohort_month_shift)
        if has_rtss > 0:
            event_list.append('Received_Vaccine')

    if not ipti_df.empty:
        has_ipti = int_suite.add_ds_ipti(cb, ipti_df, my_ds, high_access_ip_frac)
        if has_ipti > 0:
            event_list.append('Received_Campaign_Drugs')

    event_list = list(np.unique(event_list))
    add_event_counter_report(cb, event_trigger_list=event_list)
    return {}


