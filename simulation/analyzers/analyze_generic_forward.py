import argparse
import os
import sys
import pandas as pd
import numpy as np
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simulation.analyzers.analyze_generic import monthlyTreatedCasesAnalyzer
from simulation.load_paths import load_box_paths


def parse_args():
    description = "Simulation specifications"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-n",
        "--expname",
        type=str,
        help="Name of simulation experiment",
        default=None
    )
    parser.add_argument(
        "-i",
        "--expid",
        type=str,
        nargs='+',
        help="Unique ID of simulation experiment"
    )
    parser.add_argument(
        "-syr",
        "--start_year",
        type=int,
        help="First year of reports to analyze",
        default=2020
    )
    parser.add_argument(
        "-eyr",
        "--end_year",
        type=int,
        help="Last year of reports to analyze",
        default=2025
    )
    return parser.parse_args()


class CohortPfPRAnalyzerU5(BaseAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir="."):
        super(CohortPfPRAnalyzerU5, self).__init__(working_dir=working_dir,
                                                   filenames=["output/MalariaSummaryReport_U5.json"]
                                                   )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name

    def select_simulation_data(self, data, simulation):

        fname = self.filenames[0]
        d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin']
        pfpr = [x[1] for x in d]
        d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin']
        clinical_cases = [x[1] for x in d]
        d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin']
        severe_cases = [x[1] for x in d]
        d = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin']
        pop = [x[1] for x in d]
        adf = pd.DataFrame({'year': 5,
                            'PfPR U5': pfpr,
                            'Cases U5': clinical_cases,
                            'Severe cases U5': severe_cases,
                            'Pop U5': pop})

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                adf[sweep_var] = simulation.tags[sweep_var]
        return adf

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, 'U5_PfPR_ClinicalIncidence.csv')), index=False)


class CohortPfPRAnalyzerU1(BaseAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir="."):
        super(CohortPfPRAnalyzerU1, self).__init__(working_dir=working_dir,
                                                   filenames=["output/MalariaSummaryReport_U1.json"]
                                                   )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name

    def select_simulation_data(self, data, simulation):

        fname = self.filenames[0]
        d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin']
        pfpr = [x[0] for x in d]
        d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin']
        clinical_cases = [x[0] for x in d]
        d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin']
        severe_cases = [x[0] for x in d]
        d = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin']
        pop = [x[0] for x in d]
        adf = pd.DataFrame({'year': 1,
                            'PfPR U1': pfpr,
                            'Cases U1': clinical_cases,
                            'Severe cases U1': severe_cases,
                            'Pop U1': pop})

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                adf[sweep_var] = simulation.tags[sweep_var]
        return adf

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv((os.path.join(self.working_dir, self.expt_name, 'U1_PfPR_ClinicalIncidence.csv')), index=False)


if __name__ == "__main__":

    if os.name == "posix":
        SetupParser.default_block = 'NUCLUSTER'
    else:
        SetupParser.default_block = 'HPC'
    datapath, projectpath = load_box_paths(parser_default=SetupParser.default_block)
    SetupParser.init()

    working_dir = os.path.join(projectpath, 'simulation_output/generic_forward')
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    """Simulation arguments"""
    args = parse_args()
    expname = args.expname
    expid = args.expid
    start_year = args.start_year
    end_year = args.end_year

    sweep_variables = ["Run_Number", "Scenario_id", "Setting_id", "intervention_correlation", "frac_high_access",
                       "Cohort_birth_month", "Annual EIR", "seasonality", "cm_coverage",
                       "ipti_coverage", "rtss_coverage", "smc_coverage",
                       "rtss_mode"]  # ,"intervention_cor_cm","antiaccess_high"  ##"ipti_mode"

    treatment_channels = ['Received_NMF_Treatment', 'Received_Severe_Treatment', 'Received_Treatment']
    analyzers_summaryReport = [
        CohortPfPRAnalyzerU5(expt_name=expname,
                             sweep_variables=sweep_variables,
                             working_dir=working_dir),
        CohortPfPRAnalyzerU1(expt_name=expname,
                             sweep_variables=sweep_variables,
                             working_dir=working_dir)
    ]
    analyzers = [
        monthlyTreatedCasesAnalyzer(expt_name=expname,
                                    sweep_variables=sweep_variables,
                                    channels=treatment_channels,
                                    working_dir=working_dir,
                                    start_year=start_year,
                                    end_year=end_year)
    ]
    am = AnalyzeManager(expid, analyzers=analyzers)  # + analyzers_summaryReport
    am.analyze()
