import os
import sys
import numpy as np
import pandas as pd

sys.path.append('../')
from simulation.load_paths import load_box_paths


def create_eir_seasonality(projectpath, datapath):
    """create seasonality files with monthly EIR rescaling values. Include seasonalities for Phase III trial sites,
    Chandramohan site, and various sites used for vaccine sweep explorations."""
    phase3_eirs = pd.read_csv(os.path.join(datapath, 'normalizedEIRprofilesFromRTSSincidenceChildren.csv'))
    df_phase3 = pd.DataFrame(data={'month': [ii+1 for ii in range(12)]})
    for col in phase3_eirs.columns:
        if not col == 'month':
            site_name = col
            raw_eirs = phase3_eirs[col].values.tolist()
            rescaled_eirs = raw_eirs / np.sum(raw_eirs)
            df_phase3[site_name] = rescaled_eirs


    """
    OLD:  for highest-seasonality setting, use (roughly) Sahel scaling from Yamba et al. 2020
    https://www.mdpi.com/2306-5729/5/2/31/html
    eir_raw_high = [0.1,0.1,0.1,0.01,0.1,0.1,5,15,30,10,2,1]
    use the Nanoro monthly EIRs from Caitlin's RTS,S work
    """
    eir_raw_high = [0.025013202, 0.016779953, 0.00207845, 0.00, 0.002746537, 0.045913459, 0.166354365, 0.207947834,
                    0.208551483, 0.16492372, 0.104241929, 0.055449068]
    eir_high = eir_raw_high / np.sum(eir_raw_high)
    eir_high_JuneStart = [*eir_high[5:12], *eir_high[0:5]]
    eir_sahel = [0.1, 0.1, 0.1, 0.01, 0.1, 0.1, 5, 15, 30, 10, 2, 1]
    eir_sahel_JuneStart = [*eir_sahel[5:12], *eir_sahel[0:5]] / np.sum(eir_sahel)
    eir_higher_JuneStart = [2.37440071e-02, 1.22541062e-01, 2.22065556e-01, 3.40459020e-01,
                            1.61189619e-01, 6.78665164e-02, 3.55973099e-02, 1.32938786e-02,
                            9.17725409e-03, 1.82650259e-03, 7.87277594e-05,
                            2.16054609e-03]  # average between sahel and nanoro
    eir_higher_JuneStart = eir_higher_JuneStart / np.sum(eir_higher_JuneStart)

    """
    # OLD: for medium-seasonality setting, use rough values from West Africa zone from Yamba et al. 2020 
    # https://www.mdpi.com/2306-5729/5/2/31/html
    # (older version) East Africa zone: eir_raw_med = [5, 8, 10, 15, 22, 24, 25, 14, 13, 11, 8, 14]
    # (older version) Guinea zone:  eir_raw_med = [4, 3, 2, 7, 17, 18, 20, 13, 15, 13, 10, 4]
    # eir_raw_med = [12, 20, 18, 17, 22, 17, 16, 24, 43, 28, 19, 11]
    # however, switch values from month 5 and 7 to get smoother pattern and same SMC months as the high-transmission setting
    # eir_raw_med = [12, 20, 18, 17, 16, 17, 22, 24, 43, 28, 19, 11]
    # use the Kintampo monthly EIRs from Caitlin's RTS,S work
    """
    eir_raw_kintampo = [0.02676951, 0.06805808, 0.08325771, 0.10140654, 0.11887477, 0.146098, 0.10276769, 0.09233212,
                        0.07009982, 0.08779492, 0.0707804, 0.03176044]
    eir_kintampo = eir_raw_kintampo / np.sum(eir_raw_kintampo)
    eir_kintampo_JuneStart = [*eir_kintampo[5:12], *eir_kintampo[0:5]]

    """ shift schedule so that the sequence of the four highest months start in July 
    (to match the high-transmission and SMC delivery schedules)
    """
    eir_raw_med = [0.08779492, 0.0707804, 0.03176044, 0.02676951, 0.06805808, 0.08325771, 0.10140654, 0.11887477,
                   0.146098, 0.10276769, 0.09233212, 0.07009982]
    eir_med = eir_raw_med / np.sum(eir_raw_med)

    # take the averages between Nanoro and Kintamop with a June start
    eir_aveNanoroKintampo_JuneStart = [(a + b) / 2 for a, b in zip(eir_high_JuneStart, eir_kintampo_JuneStart)]
    eir_aveNanoroKintampo_JuneStart = eir_aveNanoroKintampo_JuneStart / np.sum(eir_aveNanoroKintampo_JuneStart)

    # for the constant seasonality
    eir_flat = [1 / 12] * 12

    df2 = pd.DataFrame(
        data={'high_unimodal': eir_high, 'moderate_unimodal': eir_med, 'constant': eir_flat,
              'kintampo': eir_kintampo,
              'high_unimodal_JuneStart': eir_high_JuneStart,
              'kintampo_JuneStart': eir_kintampo_JuneStart,
              'aveNanoroKintampo_JuneStart': eir_aveNanoroKintampo_JuneStart,
              'sahel_JuneStart': eir_sahel_JuneStart,
              'higher_JuneStart': eir_higher_JuneStart})
    df = pd.concat([df_phase3, df2], axis=1)
    df.to_csv(os.path.join(projectpath, 'simulation_inputs', 'seasonality',
                           'seasonality_eir_multipliers.csv'), index=False)


if __name__ == "__main__":
    datapath, projectpath = load_box_paths()
    create_eir_seasonality(projectpath, datapath)

