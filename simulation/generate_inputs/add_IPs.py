import json
import os
from dtk.utils.Campaign.CampaignClass import InterventionForCurrentPartners
import pandas as pd
import numpy as np
from input_file_generation.add_properties_to_demographics import generate_demographics_properties
import sys
#sys.path.append('../') # not needed if running from base directory
from simulation.load_paths import load_box_paths

datapath, projectpath = load_box_paths(parser_default="NUCLUSTER")
#inputs_path = os.path.join(projectpath, '..', 'hbhi_nigeria', 'simulation_inputs')
inputs_path = os.path.join(projectpath, 'simulation_inputs', 'nigeria')
#inputs_path = os.path.join(projectpath, 'simulation_inputs', 'burkina')

if __name__ == '__main__' :

    master_csv = os.path.join(projectpath, 'nigeria_LGA_pop.csv')
    #master_csv = os.path.join(projectpath, 'burkina_DS_pop.csv')
    df = pd.read_csv(master_csv, encoding='latin')
    df['LGA'] = df['LGA'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8')
    #df['DS_Name'] = df['DS_Name'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8')

    HF_times = [0.5, 1, 2, 4, 7, 50]
    inpatient_times = [0.5, 1, 2, 4, 7, 10, 15, 20, 50]

    IPs = [
        #{'Property': 'AccessToCare',
        # 'Values': ['HFWithin%.1fInpatientWithin%.1f' % (x,y) for x in HF_times for y in inpatient_times],
        # 'Initial_Distribution': [1] + [0 for x in range(len(HF_times)*len(inpatient_times)-1)],
        # 'Transitions': []},
        {'Property' : 'DrugStatus',
         'Values' : ['None', 'RecentDrug'],
         'Initial_Distribution': [1, 0],
         'Transitions' : []},
        {'Property': 'SMCAccess',
         'Values': ['Low', 'High'],
         'Initial_Distribution': [0.5, 0.5],
         'Transitions': []},
        {'Property': 'AgeGroup',
         'Values': ['Under15', '15to30', '30to50', '50plus'],
         'Initial_Distribution': [1, 0, 0, 0],
         'Transitions': []},
        {'Property': 'VaccineStatus',
         'Values': ['None', 'GotVaccine', 'GotBooster1', 'GotBooster2', 'GotBooster3'],
         'Initial_Distribution': [1, 0, 0, 0, 0],
         'Transitions': []},
    ]

    for r, row in df.iterrows() :
        hfca = row['LGA']
        #hfca = row['DS_Name']
        #tt_df = pd.read_csv(os.path.join(inputs_path, hfca, 'travel_to_care_matrix.csv'))
        #tt_df.set_index('hf time', inplace=True)

        #ddf0 = pd.DataFrame( { 'Property_Value' : IPs[0]['Values'],
        #                      'Initial_Distribution' : [tt_df.get_value(x, 'inp time %.1f' % y) for x in HF_times for y in inpatient_times]
        #                      })
        #ddf0['Property'] = 'AccessToCare'

        ddf1 = pd.DataFrame( { 'Property' : ['DrugStatus']*2,
                              'Property_Value' : ['None', 'RecentDrug'],
                              'Initial_Distribution': [1, 0]})

        ddf2 = pd.DataFrame( { 'Property' : ['SMCAccess']*2,
                              'Property_Value' : ['Low', 'High'],
                              'Initial_Distribution': [0.5, 0.5]})

        ddf3 = pd.DataFrame( { 'Property' : ['AgeGroup']*4,
                              'Property_Value' : ['Under15', '15to30', '30to50', '50plus'],
                              'Initial_Distribution': [1, 0, 0, 0]})

        ddf4 = pd.DataFrame( { 'Property' : ['VaccineStatus']*5,
                              'Property_Value' : ['None', 'GotVaccine', 'GotBooster1', 'GotBooster2', 'GotBooster3'],
                              'Initial_Distribution': [1, 0, 0, 0, 0]})

        adf = pd.concat([ddf1, ddf2, ddf3, ddf4])
        adf['Property_Type'] = 'IP'
        adf['node'] = 1

        print(hfca)
        demo_fname = os.path.join(inputs_path, hfca, '%s_demographics.json' % hfca)
        IP_demo_fname = os.path.join(inputs_path, hfca, '%s_demographics_wVaxSMC3B.json' % hfca)

        generate_demographics_properties(refdemo_fname=demo_fname,
                                         output_filename=IP_demo_fname,
                                         as_overlay=False,
                                         IPs=IPs,
                                         df=adf)
        f = open(IP_demo_fname)
        js = json.load(f)
        js['Defaults']['IndividualAttributes']["RiskDistributionFlag"] = 3

        f = open(IP_demo_fname, 'w')
        json.dump(js, f, sort_keys=True, indent=4, separators=(',', ': '))
