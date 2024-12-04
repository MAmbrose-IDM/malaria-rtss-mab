import argparse
import os
import sys
import pandas as pd
import numpy as np
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simulation.analyzers.analyzer_collections import MonthlyPfPRAnalyzer, IndividualEventsAnalyzer, \
    MonthlyU1PfPRAnalyzer, MonthlyTreatedCasesAnalyzer
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
        default=2024
    )
    return parser.parse_args()


if __name__ == "__main__":

    if os.name == "posix":
        SetupParser.default_block = 'NUCLUSTER'
    else:
        SetupParser.default_block = 'HPC'
    datapath, projectpath = load_box_paths(parser_default=SetupParser.default_block)
    SetupParser.init()

    working_dir = os.path.join(projectpath, 'simulation_output/generic')
    if not os.path.exists(working_dir):
        os.mkdir(working_dir)

    """Simulation arguments"""
    args = parse_args()
    expname = args.expname
    expid = args.expid
    start_year = args.start_year
    end_year = args.end_year

    sweep_variables = ["Run_Number", "Annual EIR", "cm_cov_U5", "coverage"]

    analyzers = [
        MonthlyPfPRAnalyzer(expt_name=expname,
                            sweep_variables=sweep_variables,
                            working_dir=working_dir,
                            start_year=start_year,
                            end_year=end_year),
        MonthlyU1PfPRAnalyzer(expt_name=expname,
                              sweep_variables=sweep_variables,
                              working_dir=working_dir,
                              start_year=start_year,
                              end_year=end_year),
        MonthlyTreatedCasesAnalyzer(expt_name=expname,
                                    sweep_variables=sweep_variables,
                                    channels=['Received_Vaccine'],
                                    working_dir=working_dir,
                                    start_year=start_year,
                                    end_year=end_year),
        # IndividualEventsAnalyzer(expt_name=expname,
        #                          sweep_variables=sweep_variables,
        #                          working_dir=working_dir,
        #                          yrs_to_keep=3,
        #                          start_year=start_year,
        #                          end_year=end_year)
    ]
    am = AnalyzeManager(expid, analyzers=analyzers)
    am.analyze()
