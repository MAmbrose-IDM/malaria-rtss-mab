import pandas as pd
import numpy as np
import math as math

from dtk.interventions.irs import add_IRS
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.utils.Campaign.CampaignClass import CampaignEvent, NodeSetAll, StandardInterventionDistributionEventCoordinator, BroadcastEvent
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



def get_concentration_at_time(tt, initial_concentration, fast_frac, k1, k2):
    concentration_at_tt = initial_concentration * (fast_frac * math.exp(-1 * tt / k1) + (1 - fast_frac) * math.exp(-1 * tt / k2))
    return concentration_at_tt


def get_time_efficacy_values(initial_concentration, max_efficacy, fast_frac, k1, k2, hh, nn, total_time):
    # get concentration and efficacy through time
    concentration_through_time = [get_concentration_at_time(tt, initial_concentration, fast_frac, k1, k2) for tt in range(total_time)]
    # efficacy_through_time = [max_efficacy * (1 - math.exp(mm * cc)) for cc in concentration_through_time]
    efficacy_through_time = [max_efficacy / (1 + math.pow((hh / cc), nn)) for cc in concentration_through_time]
    return [[i for i in range(total_time)], efficacy_through_time]


def get_vacc_params_from_pkpd_df(row):
    initial_concentration = row['initial_concentration']
    max_efficacy = row['max_efficacy']
    fast_frac = row['fast_frac']
    k1 = row['k1']
    k2 = row['k2']
    hh = row['hh']
    nn = row['nn']
    total_time = row['total_time']
    time_efficacy_values = get_time_efficacy_values(initial_concentration, max_efficacy, fast_frac, k1, k2, hh, nn, total_time)
    return time_efficacy_values


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
    vacc_ds_col = 'DS_Name'
    filter_by_ds = False
    vacc_type_col = 'vacc_types'
    vacc_start_col = 'vacc_day'
    vacc_coverage_col = 'coverage_levels'
    vacc_touchpoint_col = 'vacc_touchpoints'  # days since births!
    vacc_deploy_type_col = 'deploy_type'
    vacc_distribution_col = 'distribution_name'
    vacc_std_col = 'distribution_std'
    vacc_min_age_col = 'agemin'
    vacc_max_age_col = 'agemax'
    vacc_repetitions = 'repetitions'
    vacc_tsteps_btwn_repetitions = 'tsteps_btwn_repetitions'
    vacc_init_eff_col = 'initial_killing'
    vacc_decay_t_col = 'decay_time_constant'
    vacc_decay_class_col = 'decay_class'
    vacc_property_restrictions = 'vacc_property_restrictions'

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

    def add_triggered_vacc(self, cb, vacc_char_df, my_ds=''):

        if self.filter_by_ds:
            vacc_char_df = vacc_char_df[vacc_char_df[self.vacc_ds_col].str.upper() == my_ds.upper()]
        if 'vacc_type' in vacc_char_df.columns:
            row_initial = vacc_char_df[vacc_char_df['vacc_type'] == 'initial'].iloc[0]
            row_boost = vacc_char_df[vacc_char_df['vacc_type'] == 'booster'].iloc[0]
        else:
            row_initial = vacc_char_df.iloc[0]
            row_boost = vacc_char_df.iloc[0]

        """Set vaccine properties (e.g., initial efficacy, waning)"""
        # TODO: currently, assumes that boost has same waning type as initial dose. Should probably change this.
        try:
            waning_type = row_initial['vacc_waning_type']
        except:
            waning_type = 'exponential'

        if waning_type == 'pkpd':
            time_efficacy_values_initial = get_vacc_params_from_pkpd_df(row_initial)
            time_efficacy_values_boost = get_vacc_params_from_pkpd_df(row_boost)
            time_efficacy_initial = time_efficacy_values_initial[1][0]
            time_efficacy_multipliers = [time_efficacy_values_initial[1][yy] / time_efficacy_initial for yy in
                                         range(len(time_efficacy_values_initial[1]))]
            time_efficacy_boost_initial = time_efficacy_values_boost[1][0]
            time_efficacy_boost_multipliers = [time_efficacy_values_boost[1][yy] / time_efficacy_boost_initial for yy in
                                         range(len(time_efficacy_values_boost[1]))]

            vaccine_params_initial = {"Waning_Config": {"Initial_Effect": time_efficacy_initial,
                                                        "Durability_Map": {
                                                            "Times": time_efficacy_values_initial[0],
                                                            "Values": time_efficacy_multipliers
                                                        },
                                                        "Reference_Timer": 1,
                                                        "Expire_At_Durability_Map_End": 1,
                                                        "class": "WaningEffectMapLinear"},
                                                         "Disqualifying_Properties": ["vaccine_selected:Yes"]}
            vaccine_params_boost = {"Waning_Config": {"Initial_Effect": time_efficacy_boost_initial,
                                                      "Durability_Map": {
                                                          "Times": time_efficacy_values_boost[0],
                                                          "Values": time_efficacy_boost_multipliers
                                                      },
                                                      "Reference_Timer": 1,
                                                      "Expire_At_Durability_Map_End": 1,
                                                      "class": "WaningEffectMapLinear"},
                                                      "Disqualifying_Properties": ["vaccine_selected:Yes"]}

            # vaccine_params_initial = {"Waning_Config": {"Initial_Effect": 1,
            #                                             "Box_Duration": 400,
            #                                             "class": "WaningEffectBox"},
            #                           "Disqualifying_Properties": ["vaccine_selected:Yes"]}
            # vaccine_params_boost = {"Waning_Config": {"Initial_Effect": 0.5,
            #                                             "Box_Duration": 400,
            #                                             "class": "WaningEffectBox"},
            #                         "Disqualifying_Properties": ["vaccine_selected:Yes"]}
        else:
            raise ValueError("Unknown vaccine decay type. Only 'pkpd' currently supported.")
            # if waning_type != 'exponential':
            #     raise ValueError("Unknown vaccine decay type, assuming exponential decay")
            # vaccine_params_no_boost = {"Waning_Config": {"Initial_Effect": row['vacc_init_eff'],
            #                                              "Decay_Time_Constant": row['vacc_decay_t'],
            #                                              "class": "WaningEffectExponential"}}
            # vaccine_params_boost = vaccine_params_no_boost

        # vaccine is added in response to broadcast event
        # initial vaccine
        add_vaccine(cb,
                    vaccine_type='RTSS',
                    vaccine_params=vaccine_params_initial,
                    start_days=[0],
                    coverage=1,
                    repetitions=1,
                    tsteps_btwn_repetitions=-1,
                    listening_duration=-1,
                    ind_property_restrictions=[{'VaccineStatus': 'None'}],
                    trigger_condition_list=['event_add_new_vaccine'])
        # booster vaccines
        add_vaccine(cb,
                    vaccine_type='RTSS',
                    vaccine_params=vaccine_params_boost,
                    start_days=[0],
                    coverage=1,
                    repetitions=1,
                    tsteps_btwn_repetitions=-1,
                    listening_duration=-1,
                    ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                    trigger_condition_list=['event_add_new_vaccine'])

    def change_vacc_ips(self, cb):
        change_individual_property(cb,
                                   target_property_name='VaccineStatus',
                                   target_property_value='GotVaccine',
                                   ind_property_restrictions=[{'VaccineStatus': 'None'}],
                                   trigger_condition_list=['Received_Vaccine'],
                                   blackout_flag=False)

    def add_pkpd_vacc(self, cb, vacc_df, my_ds='', high_access_ip_frac=0, vacc_target_group='random', cohort_month_shift=0):

        # Sequence of vaccine events:
        #  - on start_day, select people to receive vaccine (change their IP vaccine_selected to True). This will remove old vaccines
        #  - on start_day+1, create a campaign which changes IP vaccine_selected to False (allowing new vaccines to be given) and broadcasts an event that triggers a node-level intervention where the new vaccine is distributed to these same individuals. The vaccine should have Disqualifying_Properties set to {“vaccine_selected”: “True”}. Also change VaccineStatus IP to ReceivedVaccine.
        if self.filter_by_ds:
            vacc_df = vacc_df[vacc_df[self.vacc_ds_col].str.upper() == my_ds.upper()]

        """Use campaign-style deployment targeted to specific ages (since births disabled in simulation)"""
        for r, row in vacc_df.iterrows():
            # calculate vaccine delivery day, given cohort month shift. If EPI type, don't adjust for cohort month
            #    (because vaccine is given according to individual's age instead of in a mass campaign)
            start_day0 = row['vacc_day']
            if row['deploy_type'] == 'EPI_cohort':
                start_day = start_day0
            elif 'season' in row['deploy_type']:
                start_day = start_day0 - round(30.4 * cohort_month_shift)
                if start_day < 0:
                    start_day = start_day + 365
            else:
                print('WARNING: vaccine delivery name not recognized, age-based vaccination.')
                start_day = start_day0

            # Select people to receive vaccine (change their IP vaccine_selected to True)
            """Set group of individuals to receive vaccine and change IPs accordingly"""
            if high_access_ip_frac > 0 and vacc_target_group in ['low', 'high']:  # different coverages in different access groups
                if vacc_target_group == 'low':
                    # vaccine is preferentially given to the 'low-access' group
                    high_low_coverages = calc_high_low_access_coverages(coverage_all=row[self.vacc_coverage_col],
                                                                        high_access_frac=1 - high_access_ip_frac)
                    high_low_coverages = [high_low_coverages[1], high_low_coverages[0]]
                elif vacc_target_group == 'high':
                    high_low_coverages = calc_high_low_access_coverages(coverage_all=row[self.vacc_coverage_col],
                                                                        high_access_frac=high_access_ip_frac)
                # high-access coverage
                change_individual_property(cb,
                                           start_day=start_day,
                                           coverage=high_low_coverages[0],
                                           target_property_name='vaccine_selected',
                                           target_property_value='Yes',
                                           target_group={'agemin': row[self.vacc_min_age_col],
                                                         'agemax': row[self.vacc_max_age_col]},
                                           ind_property_restrictions=[{'AccessToInterventions': 'higher'}],
                                           blackout_flag=False)
                # low-access coverage
                change_individual_property(cb,
                                           start_day=start_day,
                                           coverage=high_low_coverages[1],
                                           target_property_name='vaccine_selected',
                                           target_property_value='Yes',
                                           target_group={'agemin': row[self.vacc_min_age_col],
                                                         'agemax': row[self.vacc_max_age_col]},
                                           ind_property_restrictions=[{'AccessToInterventions': 'lower'}],
                                           blackout_flag=False)
            else:  # uniform probability of getting the vaccine
                if vacc_target_group != 'random':
                    print('WARNING: name for RTS,S access-group targeting not recognized, assuming random access.')
                change_individual_property(cb,
                                           start_day=start_day,
                                           coverage=row[self.vacc_coverage_col],
                                           target_property_name='vaccine_selected',
                                           target_property_value='Yes',
                                           target_group={'agemin': row[self.vacc_min_age_col],
                                                         'agemax': row[self.vacc_max_age_col]},
                                           blackout_flag=False)

            # On start_day+1, create a campaign which changes IP vaccine_selected to No (allowing new vaccines to be given) and broadcasts an event that triggers a node-level intervention where the new vaccine is distributed to these same individuals. The vaccine should have Disqualifying_Properties set to {“vaccine_selected”: “Yes”}. Also change VaccineStatus IP to ReceivedVaccine.
            vacc_broadcast_event = CampaignEvent(Start_Day=start_day+1,
                                                 Nodeset_Config=NodeSetAll(),
                                                 Event_Coordinator_Config=StandardInterventionDistributionEventCoordinator(
                                                     Demographic_Coverage=1,
                                                     Property_Restrictions=['vaccine_selected:Yes'],
                                                     Number_Repetitions=1,
                                                     Timesteps_Between_Repetitions=-1,
                                                     Intervention_Config=BroadcastEvent(
                                                         Broadcast_Event='event_add_new_vaccine',
                                                         New_Property_Value='vaccine_selected:No'),
                                                         ))
            cb.add_event(vacc_broadcast_event)

        self.change_vacc_ips(cb)
        cb.update_params({
            "Report_Event_Recorder_Events": ['Births', 'PropertyChange', 'Received_Vaccine', 'Received_Treatment']
        })

        return len(vacc_df)


    def add_ds_hs(self, cb, hs_df, my_ds, high_access_ip_frac=0, cm_target_group='random'):
        ds_col = self.hs_ds_col
        duration = self.hs_duration
        if self.filter_by_ds:
            df = hs_df[hs_df[ds_col] == my_ds]
        else:
            df = hs_df
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
        if self.filter_by_ds:
            df = smc_df[smc_df[ds_col] == my_ds]
        else:
            df = smc_df
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


    def add_ds_ipti(self, cb, ipti_df, my_ds, high_access_ip_frac=0):
        if self.filter_by_ds:
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


def add_all_interventions(cb, int_suite, my_ds='', high_access_ip_frac=0,
                          vacc_target_group='random', smc_target_group='random', cm_target_group='random',
                          hs_df=pd.DataFrame(),
                          smc_df=pd.DataFrame(),
                          vacc_char_df=pd.DataFrame(),
                          vacc_df=pd.DataFrame(),
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

    if not vacc_df.empty:
        int_suite.add_triggered_vacc(cb, vacc_char_df, my_ds)
        has_vacc = int_suite.add_pkpd_vacc(cb, vacc_df, my_ds, high_access_ip_frac, vacc_target_group, cohort_month_shift)
        if has_vacc > 0:
            event_list.append('Received_Vaccine')
            event_list.append('event_add_new_vaccine')

    if not ipti_df.empty:
        has_ipti = int_suite.add_ds_ipti(cb, ipti_df, my_ds, high_access_ip_frac)
        if has_ipti > 0:
            event_list.append('Received_Campaign_Drugs')

    event_list = list(np.unique(event_list))
    add_event_counter_report(cb, event_trigger_list=event_list)
    return {}


