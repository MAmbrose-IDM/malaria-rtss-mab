import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import sys

sys.path.append('../../')
from simulation.load_paths import load_box_paths

mpl.rcParams['pdf.fonttype'] = 42
sns.set_style('whitegrid', {'axes.linewidth': 0.5})
palette = sns.color_palette('Dark2', 10)


def plot_vaccine_report(simoutdir, simpop=20000):
    df = pd.read_csv(os.path.join(simoutdir, "ReportEventRecorder.csv"))
    df['Age_mth'] = df['Age'] / 7 / 4  # Age in months

    grp_channels = ['Event_Name', 'Age_mth', 'Individual_ID']
    df = df.loc[df['Event_Name'] == 'Received_Vaccine']
    adf = df.groupby(grp_channels)[['Node_ID']].agg('count').reset_index().rename(
        columns={'Node_ID': 'count'})

    grp_channels = ['Event_Name', 'Individual_ID']
    adf2 = adf.groupby(grp_channels)[['count']].agg('count').reset_index().rename(
        columns={'count': 'n_vaccines'})

    grp_channels = ['n_vaccines']
    adf3 = adf2.groupby(grp_channels)[['Event_Name']].agg('count').reset_index().rename(
        columns={'Event_Name': 'n_children'})
    adf3['prop_vacc'] = adf3['n_children'] / simpop  # Sim pop hardcoded
    not_vacc = simpop - sum(adf3['n_children'])
    prop_not_vacc = not_vacc / simpop
    adf3.loc[len(adf3.index)] = [0, not_vacc, prop_not_vacc]

    """Plot1: number of vaccines received, denominator all children"""
    fig = plt.figure(figsize=(6, 3))
    fig.suptitle('', y=0.98, fontsize=12)
    fig.subplots_adjust(right=0.97, left=0.10, hspace=0.4, wspace=0.2, top=0.90, bottom=0.13)
    axes = [fig.add_subplot(1, 1, x + 1) for x in range(1)]
    ax = axes[0]
    ax.bar(adf3['n_vaccines'], adf3['prop_vacc'], align='center', width=0.9, color=palette[1])
    ax.set_ylabel(f"Proportion children")
    ax.set_xlabel('Number of vaccine doses received')
    ax.set_xticks(adf3['n_vaccines'])
    fig.savefig(os.path.join(simoutdir, f'number_of_doses.png'))
    fig.savefig(os.path.join(simoutdir,'pdf', f'number_of_doses.pdf'))

    """Plot2: number of boosters received, denominator children who got initial dose"""
    adf4 = adf.copy()
    adf4['age_at_dose'] = adf4["Age_mth"].astype(int)
    adf4 = adf4.drop(['Age_mth', 'Event_Name'], axis=1)
    adf4 = adf4.pivot(index="Individual_ID", columns="age_at_dose", values="count")
    cols = adf4.columns
    adf4['ndosestotal'] = adf4.iloc[:, 0:].sum(axis=1)
    adf4 = adf4.reset_index()
    adf4['dummy'] = 'nchildren'

    adf4 = adf4.groupby('dummy')[list(cols) ].agg('sum').T.reset_index()
    initial = adf4.at[0, 'nchildren'].item()
    adf4['perc'] = adf4['nchildren'] / initial
    adf4['age_at_dose'] = [str(x) for x in adf4['age_at_dose']]

    fig = plt.figure(figsize=(6, 3))
    fig.suptitle('Proportion of children receiving boosters after initial dose', y=0.98, fontsize=12)
    fig.subplots_adjust(right=0.97, left=0.10, hspace=0.4, wspace=0.2, top=0.90, bottom=0.13)
    axes = [fig.add_subplot(1, 1, x + 1) for x in range(1)]
    ax = axes[0]
    ax.bar(adf4['age_at_dose'], adf4['perc'], align='center', width=0.9, color=palette[1])
    ax.set_ylabel(f"Proportion children")
    ax.set_xlabel('')
    fig.savefig(os.path.join(simoutdir, f'booster_coverage.png'))
    fig.savefig(os.path.join(simoutdir,'pdf', f'booster_coverage.pdf'))

    adf4['simpop'] = simpop
    adf4.to_csv(os.path.join(simoutdir, f'booster_coverage_summary.csv'))

    return adf4

def plot_vaccine_report_flexible(simoutdir, simpop=20000):
    df = pd.read_csv(os.path.join(simoutdir, "ReportEventRecorder.csv"))
    df['Age_mth'] = df['Age'] / 7 / 4  # Age in months

    grp_channels = ['Event_Name', 'Age_mth', 'Individual_ID']
    df = df.loc[df['Event_Name'] == 'Received_Vaccine']
    adf = df.groupby(grp_channels)[['Node_ID']].agg('count').reset_index().rename(
        columns={'Node_ID': 'count'})

    grp_channels = ['Event_Name', 'Individual_ID']
    adf2 = adf.groupby(grp_channels)[['count']].agg('count').reset_index().rename(
        columns={'count': 'n_vaccines'})

    grp_channels = ['n_vaccines']
    adf3 = adf2.groupby(grp_channels)[['Event_Name']].agg('count').reset_index().rename(
        columns={'Event_Name': 'n_children'})
    adf3['prop_vacc'] = adf3['n_children'] / simpop  # Sim pop hardcoded
    not_vacc = simpop - sum(adf3['n_children'])
    prop_not_vacc = not_vacc / simpop
    adf3.loc[len(adf3.index)] = [0, not_vacc, prop_not_vacc]

    """Plot1: number of vaccines received, denominator all children"""
    fig = plt.figure(figsize=(6, 3))
    fig.suptitle('', y=0.98, fontsize=12)
    fig.subplots_adjust(right=0.97, left=0.10, hspace=0.4, wspace=0.2, top=0.90, bottom=0.13)
    axes = [fig.add_subplot(1, 1, x + 1) for x in range(1)]
    ax = axes[0]
    ax.bar(adf3['n_vaccines'], adf3['prop_vacc'], align='center', width=0.9, color=palette[1])
    ax.set_ylabel(f"Proportion children")
    ax.set_xlabel('Number of vaccine doses received')
    ax.set_xticks(adf3['n_vaccines'])
    fig.savefig(os.path.join(simoutdir, f'number_of_doses.png'))
    fig.savefig(os.path.join(simoutdir,'pdf', f'number_of_doses.pdf'))

    """Plot2: number of boosters received, denominator children who got initial dose"""
    adf4 = adf.copy()
    adf4['age_at_dose'] = adf4["Age_mth"].astype(int)
    adf4 = adf4.drop(['Age_mth', 'Event_Name'], axis=1)
    adf4 = adf4.pivot(index="Individual_ID", columns="age_at_dose", values="count")
    cols = ['initial'] + [f'booster1_{i + 1}' for i in range(len(adf4.columns) - 1)]
    adf4.columns = cols
    adf4['ndosestotal'] = adf4.iloc[:, 0:].sum(axis=1)
    adf4['booster2'] = 0
    adf4['anybooster'] = 0
    adf4.loc[adf4['ndosestotal'] == max(adf4['ndosestotal']), 'booster2'] = 1
    adf4.loc[adf4['booster2'] == 1, cols[-1]] = 0
    adf4.loc[adf4['ndosestotal'] > 1, 'anybooster'] = 1
    adf4 = adf4.reset_index()
    adf4['dummy'] = 'nchildren'

    adf4 = adf4.groupby('dummy')[cols + ['booster2', 'anybooster']].agg('sum').T.reset_index()
    initial = adf4.at[0, 'nchildren'].item()
    nobooster = initial - adf4.at[len(adf4) - 1, 'nchildren'].item()
    adf4.loc[len(adf4.index)] = ['nobooster', nobooster]
    adf4 = adf4.loc[adf4["index"] != 'anybooster']
    adf4['perc'] = adf4['nchildren'] / initial

    fig = plt.figure(figsize=(6, 3))
    fig.suptitle('Proportion of children receiving boosters after initial dose', y=0.98, fontsize=12)
    fig.subplots_adjust(right=0.97, left=0.10, hspace=0.4, wspace=0.2, top=0.90, bottom=0.13)
    axes = [fig.add_subplot(1, 1, x + 1) for x in range(1)]
    ax = axes[0]
    ax.bar(adf4['index'], adf4['perc'], align='center', width=0.9, color=palette[1])
    ax.set_ylabel(f"Proportion children")
    ax.set_xlabel('')
    fig.savefig(os.path.join(simoutdir, f'booster_coverage.png'))
    fig.savefig(os.path.join(simoutdir,'pdf', f'booster_coverage.pdf'))

    adf4['simpop'] = simpop
    adf4.to_csv(os.path.join(simoutdir, f'booster_coverage_summary.csv'))

    return adf4


if __name__ == "__main__":
    run_in_test_mode = True
    _, projectpath = load_box_paths()

    exp_name =  "generic/generic_thirddoseage_campaign_RTSS_test"
    #exp_name =  "generic_forward/generic_campboost3_double_RTSS_test"
    simoutdir = os.path.join(projectpath, 'simulation_output', exp_name)
    plot_vaccine_report(simoutdir, simpop=20000)
    #plot_vaccine_report_flexible(simoutdir, simpop=20000)
