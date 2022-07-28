import pandas as pd
import numpy as np
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
import datetime
import os
import sys

sys.path.append('../')


class monthlyTreatedCasesAnalyzer(BaseAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, channels=None, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026):
        super(monthlyTreatedCasesAnalyzer, self).__init__(working_dir=working_dir,
                                                          filenames=["output/ReportEventCounter.json",
                                                                     "output/ReportMalariaFiltered.json"]
                                                          )
        self.sweep_variables = sweep_variables or ["admin_name", "Run_Number"]
        if channels is None:
            self.channels = ['Received_Treatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment']
        else:
            self.channels = channels
        self.inset_channels = ['Statistical Population', 'New Clinical Cases', 'New Severe Cases', 'PfHRP2 Prevalence']
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year

    # added to bypass failed cases
    # def filter(self, simulation):
    #     return simulation.status.name == 'Succeeded'

    def select_simulation_data(self, data, simulation):

        simdata = pd.DataFrame({x: data[self.filenames[0]]['Channels'][x]['Data'] for x in self.channels})
        simdata['Time'] = simdata.index

        d = pd.DataFrame({x: data[self.filenames[1]]['Channels'][x]['Data'] for x in self.inset_channels})
        d['Time'] = d.index

        if len(self.channels) > 0:
            simdata = pd.merge(left=simdata, right=d, on='Time')
        else:
            simdata = d
        simdata['Day'] = simdata['Time'] % 365
        simdata['month'] = simdata['Day'].apply(lambda x: self.monthparser((x + 1) % 365))
        simdata['year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))

        adf = pd.concat(selected).reset_index(drop=True)
        adf['date'] = adf.apply(lambda x: datetime.date(x['year'], x['month'], 1), axis=1)

        sum_channels = self.channels + ['New Clinical Cases', 'New Severe Cases']
        mean_channels = ['Statistical Population', 'PfHRP2 Prevalence']

        df = adf.groupby(self.sweep_variables + ['date'])[sum_channels].agg(np.sum).reset_index()
        pdf = adf.groupby(self.sweep_variables + ['date'])[mean_channels].agg(np.mean).reset_index()

        adf = pd.merge(left=pdf, right=df, on=(self.sweep_variables + ['date']))
        adf.to_csv(os.path.join(self.working_dir, self.expt_name, 'All_Age_monthly_Cases.csv'), index=False)


if __name__ == "__main__":

    from simtools.Analysis.AnalyzeManager import AnalyzeManager
    from simtools.SetupParser import SetupParser
    from simulation.load_paths import load_box_paths

    datapath, projectpath = load_box_paths()

    SetupParser.default_block = 'HPC'
    SetupParser.init()

    working_dir = os.path.join(projectpath, 'simulation_output')
    start_year = 2020  # assume simulation begins on Jan 1 of this year
    end_year = 2025 #! # general: 2029, Kintampo and Chandramohan: 2025, SweepEIR: 2039  # simulation ends on Dec 31 of this year

    expt_ids = {
        'TEST_validation_phase3_wBooster': '5a34bc52-b20e-ed11-a9fb-b88303911bc1',
        'TEST_validation_phase3_noBooster': '1ffe1a64-b40e-ed11-a9fb-b88303911bc1',
    }

    sweep_variables = ['Scenario_id', 'Run_Number', 'Annual EIR', 'seasonality', 'Cohort_birth_month',
                       'cm_coverage', 'smc_coverage', 'vacc_coverage', 'vacc_char', 'vacc_mode',
                       'cm_target_group', 'smc_target_group', 'vacc_target_group']
    report_count_channels = ['Received_Treatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment']
    for expname, expid in expt_ids.items():
        print('running expt %s' % expname)

        cur_monthlyTreatedCasesAnalyzer = monthlyTreatedCasesAnalyzer(expt_name=expname,
                                                                      channels=report_count_channels,
                                                                      sweep_variables=sweep_variables,
                                                                      working_dir=working_dir,
                                                                      start_year=start_year,
                                                                      end_year=end_year)

        am = AnalyzeManager(expid, analyzers=[cur_monthlyTreatedCasesAnalyzer], force_analyze=True)
        am.analyze()