from simtools.Analysis.SSMTAnalysis import SSMTAnalysis
from simtools.SetupParser import SetupParser
from simulation.analyzers.analyze_generic import monthlyTreatedCasesAnalyzer
from simtools.Utilities.COMPSUtilities import get_most_recent_experiment_id_by_name

wi_name_base = "ssmt_analyzer_"
working_dir = '.'
start_year = 2020  # assume simulation begins on Jan 1 of this year
end_year = 2025  # 2039  # simulation ends on Dec 31 of this year

# can specify group of experiments by name...
experiment_name_stem = ''
experiment_name_tail = ''
experiment_numbers = list(range(5))
experiments = {}
# ... or can specify specific set of experiments by id
experiments = {
    'TEST_validation_phase3_wBooster': '5a34bc52-b20e-ed11-a9fb-b88303911bc1',
    # 'TEST_validation_phase3_noBooster': '1ffe1a64-b40e-ed11-a9fb-b88303911bc1',
}

sweep_variables = ['Scenario_id', 'Run_Number', 'Annual EIR', 'seasonality', 'Cohort_birth_month',
                   'cm_coverage', 'smc_coverage', 'vacc_coverage', 'vacc_char', 'vacc_mode',
                   'cm_target_group', 'smc_target_group', 'vacc_target_group']
report_count_channels = ['Received_Treatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment']

if __name__ == "__main__":

    SetupParser.default_block = 'HPC'
    SetupParser.init()

    if not bool(experiments):
        for ee in range(len(experiment_numbers)):
            cur_exp_name = '%s%i%s' % (experiment_name_stem, experiment_numbers[ee], experiment_name_tail)
            cur_exp_id = get_most_recent_experiment_id_by_name(cur_exp_name)
            experiments[cur_exp_name] = str(cur_exp_id)

    analyzers = [monthlyTreatedCasesAnalyzer]

    for expt_name, exp_id in experiments.items():
        wi_name = '%s_%s' % (wi_name_base, expt_name)

        args_each = {'expt_name': expt_name,
                     'channels': report_count_channels,
                     'sweep_variables': sweep_variables,
                     'working_dir': working_dir,
                     'start_year': start_year,
                     'end_year': end_year}

        analysis = SSMTAnalysis(experiment_ids=[exp_id],
                                analyzers=analyzers,
                                analyzers_args=[args_each],
                                analysis_name=wi_name)
        analysis.analyze()
