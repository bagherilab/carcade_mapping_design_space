import ABM
import os
import pickle
import json
import re
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
from colour import Color
from matplotlib import cm as mplcm
from matplotlib import colors as mcolors
from matplotlib.colors import to_rgba
from argparse import ArgumentParser

'''
ABM_STATS takes a directory of (or a single) .pkl simulation files that result from ABM_ANALYZE and
does statistics on the DATA features in the dataframe over time.

Usage:
    python abm_stats.py FILES [--color COLOR] [--marker MARKER] [--saveLoc SAVELOC]

    FILES
        Path to .pkl or directory
    [--norm NORM]
        Normalize final values by initial (INIT) or untreated (UNTREATED)
    [--average]
        Average across seed replicates
    [--saveLoc SAVELOC]
        Location of where to save file, default will save here
'''

# -------------- PARSING AND INPUT/OUTPUT FUNCTIONS -------------

def get_parser():

    # Setup argument parser.
    parser = ArgumentParser(description="Plot ABM data from dataframe.")
    parser.add_argument(dest="files", help="Path to .pkl file or directory.")
    parser.add_argument("--norm", default="INIT", dest="norm",
                        help="Normalize final values by initial (INIT) or untreated (UNTREATED).")
    parser.add_argument("--score", default="SUM", dest="score",
                        help="Score outputs for simulations with cancer and healthy cells using SUM score or RATIO score.")
    parser.add_argument("--average", default=False, dest="average",
                        action='store_true', help="Average across seed replicates.")
    parser.add_argument("--partial", default=False, dest="partial",
                        action='store_true', help="Flag indicating only partial dataset present instead of full combinatoral set.")
    parser.add_argument("--saveLoc", default="", dest="saveLoc",
                        help="Location of where to save file, default will save here.")

    return parser

# Function: get_pkl_files
# Author: Andrea Collins
# If arg is a directory, returns a list of all the files in that
# directory that are pkl files. If arg is a single file, checks
# if that file is a pkl file and if it is, returns that filename.
def get_pkl_files(arg):
    if arg[-1] == "/" or arg[-1] == "\\":
        return [arg + f for f in os.listdir(arg) if ABM.is_pkl(f)]
    else:
        assert ABM.is_pkl(arg)
        return [arg]

# ---------------- GLOBALS ------------------------

FILEID_SPLIT_INDICES = {'EXP': 0, 'PLATE': 1, 'TREAT': 2, 'POPS': 3, 'DIMENSION': 4, 'DOSE': 5, 'TREAT RATIO': 6, 'CAR AFFINITY': 7, 'ANTIGENS CANCER': 8, 'ANTIGENS HEALTHY': 9 }
STRING_FEATURES = ['SEED', 'PLATE', 'TREAT RATIO']
POP_NAMES = ['CANCER', 'CANCER LIVE', 'HEALTHY', 'HEALTHY LIVE', 'T-CELL', 'T-CELL LIVE', 'CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']
POP_FRAC_NAMES = ['CANCER LIVE', 'HEALTHY LIVE']
POP_LYSIS_NAMES = ['TISSUE', 'CANCER', 'HEALTHY']
POP_CODES = {'CANCER': 0, 'HEALTHY': 1, 'CD4 T-CELL': 2, 'CD8 T-CELL': 3}
STATES_CANCER = ['NEUTR', 'APOPT', 'QUIES', 'MIGRA', 'PROLI', 'SENES', 'NECRO']
STATES_TISSUELIVE = ['NEUTR', 'QUIES', 'MIGRA', 'PROLI', 'SENES']
STATES_TCELL = ['NEUTR', 'APOPT', 'MIGRA', 'PROLI', 'SENES', 'CYTOT', 'STIMU', 'EXHAU', 'ANERG', 'STARV', 'PAUSE']
MOL_NAMES = ['GLUCOSE', 'OXYGEN', 'TGFA', 'IL-2']
MOL_CONC_UNITS = {'GLUCOSE': 'fmol/ml', 'OXYGEN': 'mmHg', 'TGFA': 'pg/ml', 'IL-2': 'pg/ml'}

DISH_TIMES = [0, 1, 4, 7]
TISSUE_TIMES = [22, 26, 29, 31]

FIG_SIZE_X = 7
FIG_SIZE_Y = 7

TICKSIZE = 20
FONTSIZE_AXES_VALUES = 20
FONTSIZE_AXES_TITLES = 15
LABELPAD = 7
SCORE_MIN_HEALTHY_THRESHOLD = 0.5

# Set DOSE color scale
DCOLORS = []
dcmap = mplcm.get_cmap('Greens')
dvalues = [0.33, 0.66, 1.0]
for v in dvalues:
    DCOLORS.append(dcmap(v))

# Set AFFINITY color scale
ACOLORS = []
acmap = mplcm.get_cmap('Oranges')
avalues = [0.2, 0.4, 0.6, 0.8, 1.0]
for v in avalues:
    ACOLORS.append(acmap(v))

# Set TREAT RATIO color scale
TRCOLORS = []
trcmap = mplcm.get_cmap('Reds')
trvalues = [0.2, 0.4, 0.6, 0.8, 1.0]
for v in trvalues:
    TRCOLORS.append(trcmap(v))

# Set ANTIGENS CANCER color scale
ACCOLORS = []
accmap = mplcm.get_cmap('Blues')
acvalues = [0.2, 0.4, 0.6, 0.8, 1.0]
for v in acvalues:
    ACCOLORS.append(accmap(v))

# Set ANTIGENS CANCER color scale
AHCOLORS = []
ahcmap = mplcm.get_cmap('Purples')
ahvalues = [0.5, 1.0]
for v in ahvalues:
    AHCOLORS.append(ahcmap(v))

CMAPS = {'DOSE': 'Greens',
         'CAR AFFINITY': 'Oranges',
         'TREAT RATIO': 'Reds',
         'ANTIGENS CANCER': 'Blues',
         'ANTIGENS HEALTHY': 'Purples'}

CMAPS_VALUES = {'DOSE': dvalues,
                'CAR AFFINITY': avalues,
                'TREAT RATIO': trvalues,
                'ANTIGENS CANCER': acvalues,
                'ANTIGENS HEALTHY': ahvalues}

doseColorDict = {
    "0": 'black',
    "250": DCOLORS[0],
    "500": DCOLORS[1],
    "1000": DCOLORS[2]
}

trColorDict = {
    "NA": 'black',
    "0:100": TRCOLORS[0],
    "25:75": TRCOLORS[1],
    "50:50": TRCOLORS[2],
    "75:25": TRCOLORS[3],
    "100:0": TRCOLORS[4],
}

trColorDictHeatMap = {
    "0.0": TRCOLORS[0],
    "0.25": TRCOLORS[1],
    "0.50": TRCOLORS[2],
    "0.75": TRCOLORS[3],
    "1.0": TRCOLORS[4],
}

affinityColorDict = {
    "NA": 'black',
    "0.0": ACOLORS[4],
    "1e-10": ACOLORS[4],
    "1e-09": ACOLORS[3],
    "1e-08": ACOLORS[2],
    "1e-07": ACOLORS[1],
    "1e-06": ACOLORS[0],
}

acColorDict = {
    "0": 'black',
    "100": ACCOLORS[0],
    "500": ACCOLORS[1],
    "1000": ACCOLORS[2],
    "5000": ACCOLORS[3],
    "10000": ACCOLORS[4]
}

ahColorDict = {
    "CONTROL": "black",
    "0": AHCOLORS[0],
    "100": AHCOLORS[1]
}

COLOR_DICT = {
    "DOSE": doseColorDict,
    "TREAT RATIO": trColorDict,
    "CAR AFFINITY": affinityColorDict,
    "ANTIGENS CANCER": acColorDict,
    "ANTIGENS HEALTHY": ahColorDict
}

COLOR_DICT_HEATMAP = {
    "DOSE": doseColorDict,
    "TREAT RATIO": trColorDictHeatMap,
    "CAR AFFINITY": affinityColorDict,
    "ANTIGENS CANCER": acColorDict,
    "ANTIGENS HEALTHY": ahColorDict
}

doseLineDict = {
    "0": 'solid',
    "250": "dotted",
    "500": "solid",
    "1000": "dashed"
}

liveLineDict = {
    "LIVE": 'dashed',
    "TOTAL": 'solid'
}

lysedLineDict = {
    'CANCER': 'solid',
    'HEALTHY': 'dashed'
}

doseMarkerDict = {
    "0": 'X',
    "250": 's',
    "500": '^',
    "1000": 'o'
}

trMarkerDict = {
    "NA": 'X',
    "0:100": 's',
    "25:75": '^',
    "50:50": 'o',
    "75:25": 'v',
    "100:0": 'D',
}

affinityMarkerDict = {
    "NA": 'X',
    "0.0": 'X',
    "1e-06": 's',
    "1e-07": '^',
    "1e-08": 'o',
    "1e-09": 'v',
    "1e-10": 'D'
}

acMarkerDict = {
    "0": 'X',
    "100": 's',
    "500": '^',
    "1000": 'o',
    "5000": 'v',
    "10000": 'D',
}

ahMarkerDict = {
    "0": 'X',
    "100": 'o'
}

noMarkerDict = {
    'o': 'o'
}

MARKER_DICT = {
    'o': noMarkerDict,
    "DOSE": doseMarkerDict,
    "TREAT RATIO": trMarkerDict,
    "CAR AFFINITY": affinityMarkerDict,
    "ANTIGENS CANCER": acMarkerDict,
    "ANTIGENS HEALTHY": ahMarkerDict
}

AXES_SETS = {
    "DOSE": [250, 500, 1000],
    "TREAT RATIO": ["0:100", "25:75", "50:50", "75:25", "100:0"],
    "CAR AFFINITY": [1e-6, 1e-7, 1e-8, 1e-9],
    "ANTIGENS CANCER": [100, 500, 1000, 5000, 10000],
    "ANTIGENS HEALTHY": [0, 100]
    }

TREAT_RATIO_DICT = {'0:100': 0.0, '25:75': 0.25, '50:50': 0.5, '75:25': 0.75, '100:0': 1.0}
TREAT_RATIO_DICT_REVERSE = {'0.0': '0:100', '0.25': '25:75', '0.50': '50:50', '0.75': '75:25', '1.0': '100:0'}

CONSTANTS = {
        "DOSE": 500,
        "TREAT RATIO": 0.5,
        "CAR AFFINITY": 1e-7,
        "ANTIGENS CANCER": 1000,
        "ANTIGENS HEALTHY": [0, 100]
        }

# ------------- COLOR FUNCTIONS -----------------

def NonLinCdict(steps, col_array):
    cdict = {'red': (), 'green': (), 'blue': ()}
    for s, col in zip(steps, col_array):
        rgb =mcolors.to_rgb(col)
        cdict['red'] = cdict['red'] + ((s, rgb[0], rgb[0]),)
        cdict['green'] = cdict['green'] + ((s, rgb[1], rgb[1]),)
        cdict['blue'] = cdict['blue'] + ((s, rgb[2], rgb[2]),)
    return cdict

def maximum_absolute_scaling(df):
    # copy the dataframe
    df_scaled = df.copy()
    # apply maximum absolute scaling
    for column in df_scaled.columns:
        if 'NORM' not in column and 'SCORE' not in column:
            df_scaled[column] = df_scaled[column] / df_scaled[column].abs().max()
    return df_scaled

# --------- STATS/ANALYSIS FUNCTIONS ------------

def average_conditions(simsDF, FILEID):

    simsDFavg = pd.DataFrame({'DOSE': pd.Series([], dtype='float'),
                              'TREAT_RATIO': pd.Series([], dtype='float'),
                              'CAR_AFFINITY': pd.Series([], dtype='float'),
                              'ANTIGENS_CANCER': pd.Series([], dtype='float'),
                              'ANTIGENS_HEALTHY': pd.Series([], dtype='float'),
                              'Y_NORM_CANCER_LIVE': pd.Series([], dtype='float'),
                              'Y_NORM_HEALTHY_LIVE': pd.Series([], dtype='float'),
                              'Y_NORM_TCELL_LIVE': pd.Series([], dtype='float'),
                              'SCORE': pd.Series([], dtype=float)})

    features = ['DOSE', 'TREAT RATIO', 'CAR AFFINITY', 'ANTIGENS CANCER']
    if '_CH_' in FILEID:
        features.append('ANTIGENS HEALTHY')

    cancerDict = {}
    healthyDict = {}
    tcellDict = {}
    scoreDict = {}

    for d in AXES_SETS['DOSE']:
        d = str(d)
        cancerDict[d] = {}
        healthyDict[d] = {}
        tcellDict[d] = {}
        scoreDict[d] = {}
        for tr in AXES_SETS['TREAT RATIO']:
            cancerDict[d][tr] = {}
            healthyDict[d][tr] = {}
            tcellDict[d][tr] = {}
            scoreDict[d][tr] = {}
            for ca in AXES_SETS['CAR AFFINITY']:
                ca = str(ca)
                cancerDict[d][tr][ca] = {}
                healthyDict[d][tr][ca] = {}
                tcellDict[d][tr][ca] = {}
                scoreDict[d][tr][ca] = {}
                for ac in AXES_SETS['ANTIGENS CANCER']:
                    ac = str(ac)
                    if '_CH_' in FILEID:
                        cancerDict[d][tr][ca][ac] = {}
                        healthyDict[d][tr][ca][ac] = {}
                        tcellDict[d][tr][ca][ac] = {}
                        scoreDict[d][tr][ca][ac] = {}

                        for ah in AXES_SETS['ANTIGENS HEALTHY']:
                            ah = str(ah)
                            cancerDict[d][tr][ca][ac][ah] = []
                            healthyDict[d][tr][ca][ac][ah] = []
                            tcellDict[d][tr][ca][ac][ah] = []
                            scoreDict[d][tr][ca][ac][ah] = []

                    else:
                        cancerDict[d][tr][ca][ac] = []
                        healthyDict[d][tr][ca][ac] = []
                        tcellDict[d][tr][ca][ac] = []
                        scoreDict[d][tr][ca][ac] = []

    for i in range(0, len(simsDF)):
        dose = str(int(simsDF.iloc[i]['DOSE']))
        treatratio = str(simsDF.iloc[i]['TREAT_RATIO'])
        caraffinity = str(simsDF.iloc[i]['CAR_AFFINITY'])
        antigenscancer = str(int(simsDF.iloc[i]['ANTIGENS_CANCER']))

        if treatratio == '0': treatratio = '0.0'
        elif treatratio == '0.5': treatratio = '0.50'
        elif treatratio == '1': treatratio = '1.0'

        treatratio = TREAT_RATIO_DICT_REVERSE[treatratio]

        if '_CH_' in FILEID:
            antigenshealthy = str(int(simsDF.iloc[i]['ANTIGENS_HEALTHY']))
            cancerDict[dose][treatratio][caraffinity][antigenscancer][antigenshealthy].append(simsDF.iloc[i]['Y_NORM_CANCER_LIVE'])
            healthyDict[dose][treatratio][caraffinity][antigenscancer][antigenshealthy].append(simsDF.iloc[i]['Y_NORM_HEALTHY_LIVE'])
            tcellDict[dose][treatratio][caraffinity][antigenscancer][antigenshealthy].append(simsDF.iloc[i]['Y_NORM_TCELL_LIVE'])
            scoreDict[dose][treatratio][caraffinity][antigenscancer][antigenshealthy].append(simsDF.iloc[i]['SCORE'])
        else:
            cancerDict[dose][treatratio][caraffinity][antigenscancer].append(simsDF.iloc[i]['Y_NORM_CANCER_LIVE'])
            tcellDict[dose][treatratio][caraffinity][antigenscancer].append(simsDF.iloc[i]['Y_NORM_TCELL_LIVE'])

    for d in AXES_SETS['DOSE']:
        for tr in AXES_SETS['TREAT RATIO']:
            for ca in AXES_SETS['CAR AFFINITY']:
                for ac in AXES_SETS['ANTIGENS CANCER']:
                    if '_CH_' in FILEID:
                        for ah in AXES_SETS['ANTIGENS HEALTHY']:
                            if len(cancerDict[str(d)][tr][str(ca)][str(ac)][str(ah)]) == 0:
                                continue
                            else:
                                simsDictAvg = {}
                                simsDictAvg['DOSE'] = float(d)
                                simsDictAvg['TREAT_RATIO'] = float(TREAT_RATIO_DICT[tr])
                                simsDictAvg['CAR_AFFINITY'] = float(ca)
                                simsDictAvg['ANTIGENS_CANCER'] = float(ac)
                                simsDictAvg['ANTIGENS_HEALTHY'] = float(ah)
                                d = str(d)
                                ca = str(ca)
                                ac = str(ac)
                                ah = str(ah)
                                canceravg = float(sum(cancerDict[d][tr][ca][ac][ah]) / len(cancerDict[d][tr][ca][ac][ah]))
                                healthyavg = float(sum(healthyDict[d][tr][ca][ac][ah]) / len(healthyDict[d][tr][ca][ac][ah]))
                                tcellavg = float(sum(tcellDict[d][tr][ca][ac][ah]) / len(tcellDict[d][tr][ca][ac][ah]))
                                scoreavg = float(sum(scoreDict[d][tr][ca][ac][ah]) / len(scoreDict[d][tr][ca][ac][ah]))
                                simsDictAvg['Y_NORM_CANCER_LIVE'] = float(canceravg)
                                simsDictAvg['Y_NORM_HEALTHY_LIVE'] = float(healthyavg)
                                simsDictAvg['Y_NORM_TCELL_LIVE'] = float(tcellavg)
                                simsDictAvg['SCORE'] = float(scoreavg)
                                simsDFavg = simsDFavg.append(simsDictAvg, ignore_index=True)

                    else:
                        if len(cancerDict[str(d)][tr][str(ca)][str(ac)]) == 0:
                            continue
                        else:
                            simsDictAvg = {}
                            simsDictAvg['DOSE'] = float(d)
                            simsDictAvg['TREAT_RATIO'] = float(TREAT_RATIO_DICT[tr])
                            simsDictAvg['CAR_AFFINITY'] = float(ca)
                            simsDictAvg['ANTIGENS_CANCER'] = float(ac)
                            d = str(d)
                            ca = str(ca)
                            ac = str(ac)
                            canceravg = float(sum(cancerDict[d][tr][ca][ac]) / len(cancerDict[d][tr][ca][ac]))
                            tcellavg = float(sum(tcellDict[d][tr][ca][ac]) / len(tcellDict[d][tr][ca][ac]))
                            simsDictAvg['Y_NORM_CANCER_LIVE'] = float(canceravg)
                            simsDictAvg['Y_NORM_TCELL_LIVE'] = float(tcellavg)
                            simsDFavg = simsDFavg.append(simsDictAvg, ignore_index=True)
    return simsDFavg

def make_empty_features_dict():
    doseDict = {
        "0": 0,
        "250": 0,
        "500": 0,
        "1000": 0
    }

    trDict = {
        "NA": 0,
        "0:100": 0,
        "25:75": 0,
        "50:50": 0,
        "75:25": 0,
        "100:0": 0,
    }

    affinityDict = {
        "NA": 0,
        "0.0": 0,
        "1e-06": 0,
        "1e-07": 0,
        "1e-08": 0,
        "1e-09": 0,
    }

    acDict = {
        "0": 0,
        "100": 0,
        "500": 0,
        "1000": 0,
        "5000": 0,
        "10000": 0
    }

    ahDict = {
        "0": 0,
        "100": 0
    }

    FEATURE_DICT = {
        "DOSE": doseDict,
        "TREAT RATIO": trDict,
        "CAR AFFINITY": affinityDict,
        "ANTIGENS CANCER": acDict,
        "ANTIGENS HEALTHY": ahDict
    }

    return FEATURE_DICT

def anova(simsDF, FILEID, NORM, SCORE, SAVELOC):

    if '_CH_' not in FILEID:
        add = ''

        print('ANOVA FOR CANCER CELLS.')
        # model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + '
                           'C(DOSE):C(TREAT_RATIO) + C(DOSE):C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) +'
                           'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) +'
                           'C(CAR_AFFINITY):C(ANTIGENS_CANCER)', data=simsDF).fit()
        anova_CANCER = sm.stats.anova_lm(model_CANCER, typ=2)
        print(anova_CANCER)

        f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'ANOVA_CANCER.xlsx'
        anova_CANCER.to_excel(f, header=True)

        print('ANOVA FOR T-CELLS.')
        # model_TCELL = ols('Y_NORM_TCELL_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_TCELL = ols('Y_NORM_TCELL_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + '
                          'C(DOSE):C(TREAT_RATIO) + C(DOSE):C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) +'
                          'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) +'
                          'C(CAR_AFFINITY):C(ANTIGENS_CANCER)', data=simsDF).fit()
        anova_TCELL = sm.stats.anova_lm(model_TCELL, typ=2)
        print(anova_TCELL)

        f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'ANOVA_TCELL.xlsx'
        anova_TCELL.to_excel(f, header=True)

    else:
        add = SCORE + '_'

        print('ANOVA FOR CANCER CELLS.')
        #model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + C(ANTIGENS_HEALTHY) + '
                           'C(DOSE):C(TREAT_RATIO) + C(DOSE):C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) + C(DOSE):C(ANTIGENS_HEALTHY) +'
                           'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) + C(TREAT_RATIO):C(ANTIGENS_HEALTHY) +'
                           'C(CAR_AFFINITY):C(ANTIGENS_CANCER) + C(CAR_AFFINITY):C(ANTIGENS_HEALTHY) + '
                           'C(ANTIGENS_CANCER):C(ANTIGENS_HEALTHY)', data=simsDF).fit()
        anova_CANCER = sm.stats.anova_lm(model_CANCER, typ=2)
        print(anova_CANCER)

        f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'ANOVA_CANCER.xlsx'
        anova_CANCER.to_excel(f, header=True)

        print('ANOVA FOR HEALTHY CELLS.')
        # model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_HEALTHY = ols('Y_NORM_HEALTHY_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + C(ANTIGENS_HEALTHY) + '
                            'C(DOSE):C(TREAT_RATIO) + C(DOSE):C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) + C(DOSE):C(ANTIGENS_HEALTHY) +'
                            'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) + C(TREAT_RATIO):C(ANTIGENS_HEALTHY) +'
                            'C(CAR_AFFINITY):C(ANTIGENS_CANCER) + C(CAR_AFFINITY):C(ANTIGENS_HEALTHY) +'
                            'C(ANTIGENS_CANCER):C(ANTIGENS_HEALTHY)', data=simsDF).fit()
        anova_HEALTHY = sm.stats.anova_lm(model_HEALTHY, typ=2)
        print(anova_HEALTHY)

        f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'ANOVA_HEALTHY.xlsx'
        anova_HEALTHY.to_excel(f, header=True)

        print('ANOVA FOR T-CELLS.')
        # model_TCELL = ols('Y_NORM_TCELL_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_TCELL = ols('Y_NORM_TCELL_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + C(ANTIGENS_HEALTHY) + '
                          'C(DOSE):C(TREAT_RATIO) + C(DOSE):C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) + C(DOSE):C(ANTIGENS_HEALTHY) +'
                          'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) + C(TREAT_RATIO):C(ANTIGENS_HEALTHY) +'
                          'C(CAR_AFFINITY):C(ANTIGENS_CANCER) + C(CAR_AFFINITY):C(ANTIGENS_HEALTHY) +'
                          'C(ANTIGENS_CANCER):C(ANTIGENS_HEALTHY)', data=simsDF).fit()
        anova_TCELL = sm.stats.anova_lm(model_TCELL, typ=2)
        print(anova_TCELL)

        f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'ANOVA_TCELL.xlsx'
        anova_TCELL.to_excel(f, header=True)

        print('ANOVA FOR SCORE.')
        # model_CANCER = ols('SCORE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_SCORE = ols('SCORE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + C(ANTIGENS_HEALTHY) + '
                          'C(DOSE):C(TREAT_RATIO) + C(DOSE):C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) + C(DOSE):C(ANTIGENS_HEALTHY) +'
                          'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) + C(TREAT_RATIO):C(ANTIGENS_HEALTHY) +'
                          'C(CAR_AFFINITY):C(ANTIGENS_CANCER) + C(CAR_AFFINITY):C(ANTIGENS_HEALTHY) +'
                          'C(ANTIGENS_CANCER):C(ANTIGENS_HEALTHY)', data=simsDF).fit()
        anova_SCORE = sm.stats.anova_lm(model_SCORE, typ=2)
        print(anova_SCORE)

        f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'ANOVA_SCORE.xlsx'
        anova_SCORE.to_excel(f, header=True)

    return

def feature_analysis(simsDF, RANK, FILEID, NORM, SCORE, SAVELOC):

    if RANK == 'Y_NORM_CANCER_LIVE':
        simsDF = simsDF[simsDF[RANK] < 1]

    elif RANK == 'SCORE':
        simsDF = simsDF[simsDF[RANK] > 0]

    elif RANK == 'Y_NORM_HEALTHY_LIVE':

        simsDF = simsDF[simsDF[RANK] >= SCORE_MIN_HEALTHY_THRESHOLD]

    featuresDict = make_empty_features_dict()

    for i in range(0, len(simsDF)):
        for feature in featuresDict:
            if '_CH_' not in FILEID and feature == 'ANTIGENS HEALTHY':
                continue
            else:
                if feature != 'TREAT RATIO' and feature != 'CAR AFFINITY':
                    key = int(simsDF.iloc[i][feature.replace(' ','_')])
                    key = str(key)
                else:
                    key = str(simsDF.iloc[i][feature.replace(' ','_')])

                if key == '0' and feature == 'TREAT_RATIO': key = '0.0'
                elif key == '0.5': key = '0.50'
                elif key == '1': key = '1.0'

                if feature == 'TREAT RATIO':
                    key = TREAT_RATIO_DICT_REVERSE[str(key)]

                featuresDict[feature][key] += 1

    if '_CH_' in FILEID:
        add = SCORE + '_'
    else:
        add = ''

    with open(SAVELOC + FILEID + '_' + NORM + '_' + add + 'FEATURE_VALUES' + '_' + RANK.replace('_','') + '.json', 'w') as f:
        json_obj = json.dumps(featuresDict, indent=4)
        f.write(json_obj)

    return

def save_sorted_df_xlsx(simsDF, RANK, FILEID, NORM, SCORE, SAVELOC):

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']

    if '_CH_' in FILEID:
        outputs = ['SCORE', 'Y_NORM_CANCER_LIVE']
    else:
        outputs = ['Y_NORM_CANCER_LIVE']

    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')
        outputs.append('Y_NORM_HEALTHY_LIVE')

    outputs.remove(RANK)

    sortby = [RANK] + outputs + columns

    ascending = []
    for c in sortby:
        if c == 'CAR_AFFINITY' or c == 'SCORE' or c == 'Y_NORM_HEALTHY_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    if '_CH_' in FILEID:
        add = SCORE + '_'
    else:
        add = ''
    f = SAVELOC + FILEID + '_BEFORERANK_' + NORM + '_' + add + 'SORTED' + '_' + RANK.replace('_', '') + '.xlsx'
    simsDF.to_excel(f, header=True)

    if RANK == 'Y_NORM_CANCER_LIVE':
        simsDF = simsDF[simsDF[RANK] < 1]

    elif RANK == 'SCORE':
        simsDF = simsDF[simsDF[RANK] > 0]

    elif RANK == 'Y_NORM_HEALTHY_LIVE':

        simsDF = simsDF[simsDF[RANK] >= SCORE_MIN_HEALTHY_THRESHOLD]

    if '_CH_' in FILEID:
        add = SCORE + '_'
    else:
        add = ''

    f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'SORTED' + '_' + RANK.replace('_','') + '.xlsx'
    simsDF.to_excel(f, header=True)

    return

def save_sorted_df_xlsx_all(simsDF, RANK, FILEID, NORM, SCORE, SAVELOC):

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']

    if '_CH_' in FILEID:
        outputs = ['SCORE', 'Y_NORM_CANCER_LIVE']
    else:
        outputs = ['Y_NORM_CANCER_LIVE']

    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')
        outputs.append('Y_NORM_HEALTHY_LIVE')

    outputs.remove(RANK)

    sortby = [RANK] + outputs + columns

    ascending = []
    for c in sortby:
        if c == 'CAR_AFFINITY' or c == 'SCORE' or c == 'Y_NORM_HEALTHY_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    if '_CH_' in FILEID:
        add = SCORE + '_'
    else:
        add = ''

    f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'SORTED_ALL' + '_' + RANK.replace('_','') + '.xlsx'
    simsDF.to_excel(f, header=True)

    return

# ----------- PLOTTING FUNCTIONS ----------------

def plot_two_axis_heatmaps_diverging_all(simsDF, X, Y, COLOR, ANNOTATE, FILEID, NORM, SCORE, SAVELOC):

    X_AXIS = AXES_SETS[X]
    Y_AXIS = AXES_SETS[Y]

    if Y == 'ANTIGENS CANCER' or Y == 'CAR AFFINITY':
        Y_AXIS = list(reversed(Y_AXIS))

    response = 'CANCER'
    cmap = 'coolwarm'
    center = 1.0
    vmax = 3.0
    vmin = 0
    if '_CH_' in FILEID:
        vmax = 5
        vmin = 0
    if COLOR == 'Y_NORM_HEALTHY_LIVE':
        response = 'HEALTHY'
        cmap = 'coolwarm_r'
        vmax = 1.5
    if COLOR == 'Y_NORM_TCELL_LIVE':
        response = 'TCELL'
        vmax = 300
        cmap = 'coolwarm_r'
    if COLOR == 'SCORE':
        response = 'SCORE'
        cmap = 'coolwarm_r'
        vmax = 2
        center = 0

    heatmapDict = {}

    for y in Y_AXIS:
        heatmapDict[y] = {}
        for x in X_AXIS:
            heatmapDict[y][x] = []

    for y in Y_AXIS:
        for x in X_AXIS:
            if Y == 'TREAT RATIO':
                subsetDF = simsDF[simsDF[Y.replace(" ", "_")] == TREAT_RATIO_DICT[y]]
                subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]
            elif X == 'TREAT RATIO':
                subsetDF = simsDF[simsDF[Y.replace(" ", "_")] == y]
                subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == TREAT_RATIO_DICT[x]]
            else:
                subsetDF = simsDF[simsDF[Y.replace(" ", "_")] == y]
                subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]

            for i in range(0,len(subsetDF)):
                heatmapDict[y][x].append(subsetDF.iloc[i][COLOR])

    heatmapDF = pd.DataFrame(index=Y_AXIS, columns=X_AXIS, dtype=float)

    for y in Y_AXIS:
        for x in X_AXIS:
            heatmapDF.at[y, x] = float(sum(heatmapDict[y][x]) / len(heatmapDict[y][x]))

    figHeatmap = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figHeatmap.add_subplot(1, 1, 1)
    ax = sns.heatmap(heatmapDF, cmap=cmap, annot=ANNOTATE, square=True,
                     xticklabels=X_AXIS, yticklabels=Y_AXIS, center=center,
                     cbar=True, cbar_kws={"shrink": 0.82}, vmax=vmax, vmin=vmin)
    ax.set_ylabel(Y, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_xlabel(X, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(COLOR + "\nRESPONSE - ALL", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    if ANNOTATE:
        annotate = '_ANNOTATED'
    else:
        annotate = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_HEATMAP_DIVERGING_ALL_" + response + "_" + X.replace(" ",'') + "_" + Y.replace(' ','') + annotate + '.svg', bbox_inches='tight')

    return

def plot_two_axis_heatmaps_diverging_constant(simsDF, X, Y, COLOR, ANNOTATE, FILEID, NORM, SCORE, SAVELOC):

    X_AXIS = AXES_SETS[X]
    Y_AXIS = AXES_SETS[Y]

    if Y == 'ANTIGENS CANCER' or Y == 'CAR AFFINITY':
        Y_AXIS = list(reversed(Y_AXIS))

    for key in AXES_SETS:
        key_ = key.replace(" ",'_')
        if key != X and key != Y and key != 'ANTIGENS HEALTHY':
            simsDF = simsDF[simsDF[key_] == CONSTANTS[key]]

    DFlist = []

    if X != 'ANTIGENS HEALTHY' and Y != 'ANTIGENS HEALTHY' and '_CH_' in FILEID:
        for i in range(0, len(CONSTANTS['ANTIGENS HEALTHY'])):
            DFlist.append(simsDF[simsDF['ANTIGENS_HEALTHY'] == CONSTANTS['ANTIGENS HEALTHY'][i]])
    else:
        DFlist.append(simsDF)

    response = 'CANCER'
    cmap = 'coolwarm'
    center = 1.0
    vmax = 3.0
    vmin = 0
    if '_CH_' in FILEID:
        vmax = 5
        vmin = 0
    if COLOR == 'Y_NORM_HEALTHY_LIVE':
        response = 'HEALTHY'
        cmap = 'coolwarm_r'
        vmax = 1.5
    if COLOR == 'Y_NORM_TCELL_LIVE':
        response = 'TCELL'
        cmap = 'coolwarm_r'
        vmax = 300
    if COLOR == 'SCORE':
        response = 'SCORE'
        cmap = 'coolwarm_r'
        vmax = 2
        center = 0

    for d in range(0,len(DFlist)):

        DF = DFlist[d]
        heatmapDict = {}

        for y in Y_AXIS:
            heatmapDict[y] = {}
            for x in X_AXIS:
                heatmapDict[y][x] = []

        for y in Y_AXIS:
            for x in X_AXIS:
                if Y == 'TREAT RATIO':
                    subsetDF = DF[DF[Y.replace(" ", "_")] == TREAT_RATIO_DICT[y]]
                    subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]
                elif X == 'TREAT RATIO':
                    subsetDF = DF[DF[Y.replace(" ", "_")] == y]
                    subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == TREAT_RATIO_DICT[x]]
                else:
                    subsetDF = DF[DF[Y.replace(" ", "_")] == y]
                    subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]

                for i in range(0, len(subsetDF)):
                    heatmapDict[y][x].append(subsetDF.iloc[i][COLOR])

        heatmapDF = pd.DataFrame(index=Y_AXIS, columns=X_AXIS, dtype=float)

        for y in Y_AXIS:
            for x in X_AXIS:
                heatmapDF.at[y, x] = float(sum(heatmapDict[y][x]) / len(heatmapDict[y][x]))

        figHeatmap = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
        ax = figHeatmap.add_subplot(1, 1, 1)
        ax = sns.heatmap(heatmapDF, cmap=cmap, annot=ANNOTATE, square=True,
                         xticklabels=X_AXIS, yticklabels=Y_AXIS, center=center,
                         cbar=True, cbar_kws={"shrink": 0.82}, vmax=vmax, vmin=vmin)
        ax.set_ylabel(Y, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_xlabel(X, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_title(COLOR + "\nRESPONSE - CONSTANT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
        add_HA = ''
        if X != 'ANTIGENS HEALTHY' and Y != 'ANTIGENS HEALTHY' and '_CH_' in FILEID:
            add_HA = '_HA_' + str(CONSTANTS['ANTIGENS HEALTHY'][d])
        ax.set_title(COLOR + " RESPONSE - CONSTANT" + add_HA.replace('_',''), fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)

        if '_CH_' in FILEID:
            add_score = "_" + SCORE
        else:
            add_score = ''

        if ANNOTATE:
            annotate = '_ANNOTATED'
        else:
            annotate = ''

        plt.xticks(fontsize=TICKSIZE)
        plt.yticks(fontsize=TICKSIZE)

        plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_DIVERGING_CONSTANT_" + response + "_" + X.replace(" ",'') + "_" + Y.replace(' ','') + add_HA + annotate + '.svg', bbox_inches='tight')

    return

def plot_two_axis_heatmaps_grayscale_all(simsDF, X, Y, COLOR, ANNOTATE, FILEID, NORM, SCORE, SAVELOC):

    X_AXIS = AXES_SETS[X]
    Y_AXIS = AXES_SETS[Y]

    if Y == 'ANTIGENS CANCER' or Y == 'CAR AFFINITY':
        Y_AXIS = list(reversed(Y_AXIS))

    response = 'CANCER'
    cmap = 'Greys'
    vmax = 3.0
    vmin = 0
    if '_CH_' in FILEID:
        vmax = 5
    if COLOR == 'Y_NORM_HEALTHY_LIVE':
        response = 'HEALTHY'
        vmax = 1.5
    if COLOR == 'Y_NORM_TCELL_LIVE':
        response = 'TCELL'
        vmax = 300
    if COLOR == 'SCORE':
        response = 'SCORE'
        vmax = 2

    heatmapDict = {}

    for y in Y_AXIS:
        heatmapDict[y] = {}
        for x in X_AXIS:
            heatmapDict[y][x] = []

    for y in Y_AXIS:
        for x in X_AXIS:
            if Y == 'TREAT RATIO':
                subsetDF = simsDF[simsDF[Y.replace(" ", "_")] == TREAT_RATIO_DICT[y]]
                subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]
            elif X == 'TREAT RATIO':
                subsetDF = simsDF[simsDF[Y.replace(" ", "_")] == y]
                subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == TREAT_RATIO_DICT[x]]
            else:
                subsetDF = simsDF[simsDF[Y.replace(" ", "_")] == y]
                subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]

            for i in range(0,len(subsetDF)):
                heatmapDict[y][x].append(subsetDF.iloc[i][COLOR])

    heatmapDF = pd.DataFrame(index=Y_AXIS, columns=X_AXIS, dtype=float)

    for y in Y_AXIS:
        for x in X_AXIS:
            heatmapDF.at[y, x] = float(sum(heatmapDict[y][x]) / len(heatmapDict[y][x]))

    figHeatmap = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figHeatmap.add_subplot(1, 1, 1)
    ax = sns.heatmap(heatmapDF, cmap=cmap, annot=ANNOTATE, square=True,
                     xticklabels=X_AXIS, yticklabels=Y_AXIS, cbar=True, cbar_kws={"shrink": 0.82}, vmax=vmax, vmin=vmin)
    ax.set_ylabel(Y, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_xlabel(X, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(COLOR + "\nRESPONSE - ALL", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    if ANNOTATE:
        annotate = '_ANNOTATED'
    else:
        annotate = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_HEATMAP_GRAYSCALE_ALL_" + response + "_" + X.replace(" ",'') + "_" + Y.replace(' ','') + annotate + '.svg', bbox_inches='tight')

    return

def plot_two_axis_heatmaps_grayscale_constant(simsDF, X, Y, COLOR, ANNOTATE, FILEID, NORM, SCORE, SAVELOC):

    X_AXIS = AXES_SETS[X]
    Y_AXIS = AXES_SETS[Y]

    if Y == 'ANTIGENS CANCER' or Y == 'CAR AFFINITY':
        Y_AXIS = list(reversed(Y_AXIS))

    for key in AXES_SETS:
        key_ = key.replace(" ",'_')
        if key != X and key != Y and key != 'ANTIGENS HEALTHY':
            simsDF = simsDF[simsDF[key_] == CONSTANTS[key]]

    DFlist = []

    if X != 'ANTIGENS HEALTHY' and Y != 'ANTIGENS HEALTHY' and '_CH_' in FILEID:
        for i in range(0, len(CONSTANTS['ANTIGENS HEALTHY'])):
            DFlist.append(simsDF[simsDF['ANTIGENS_HEALTHY'] == CONSTANTS['ANTIGENS HEALTHY'][i]])
    else:
        DFlist.append(simsDF)

    response = 'CANCER'
    cmap = 'Greys'
    vmax = 3.0
    vmin = 0
    if '_CH_' in FILEID:
        vmax = 5
    if COLOR == 'Y_NORM_HEALTHY_LIVE':
        response = 'HEALTHY'
        vmax = 1.5
    if COLOR == 'Y_NORM_TCELL_LIVE':
        response = 'TCELL'
        vmax = 300
    if COLOR == 'SCORE':
        response = 'SCORE'
        vmax = 2

    for d in range(0,len(DFlist)):

        DF = DFlist[d]
        heatmapDict = {}

        for y in Y_AXIS:
            heatmapDict[y] = {}
            for x in X_AXIS:
                heatmapDict[y][x] = []

        for y in Y_AXIS:
            for x in X_AXIS:
                if Y == 'TREAT RATIO':
                    subsetDF = DF[DF[Y.replace(" ", "_")] == TREAT_RATIO_DICT[y]]
                    subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]
                elif X == 'TREAT RATIO':
                    subsetDF = DF[DF[Y.replace(" ", "_")] == y]
                    subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == TREAT_RATIO_DICT[x]]
                else:
                    subsetDF = DF[DF[Y.replace(" ", "_")] == y]
                    subsetDF = subsetDF[subsetDF[X.replace(" ", "_")] == x]

                for i in range(0, len(subsetDF)):
                    heatmapDict[y][x].append(subsetDF.iloc[i][COLOR])

        heatmapDF = pd.DataFrame(index=Y_AXIS, columns=X_AXIS, dtype=float)

        for y in Y_AXIS:
            for x in X_AXIS:
                heatmapDF.at[y, x] = float(sum(heatmapDict[y][x]) / len(heatmapDict[y][x]))

        figHeatmap = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
        ax = figHeatmap.add_subplot(1, 1, 1)
        ax = sns.heatmap(heatmapDF, cmap=cmap, annot=ANNOTATE, square=True,
                         xticklabels=X_AXIS, yticklabels=Y_AXIS, cbar=True, cbar_kws={"shrink": 0.82}, vmax=vmax, vmin=vmin)
        ax.set_ylabel(Y, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_xlabel(X, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_title(COLOR + "\nRESPONSE - CONSTANT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
        add_HA = ''
        if X != 'ANTIGENS HEALTHY' and Y != 'ANTIGENS HEALTHY' and '_CH_' in FILEID:
            add_HA = '_HA_' + str(CONSTANTS['ANTIGENS HEALTHY'][d])
        ax.set_title(COLOR + " RESPONSE - CONSTANT" + add_HA.replace('_',''), fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)

        if '_CH_' in FILEID:
            add_score = "_" + SCORE
        else:
            add_score = ''

        if ANNOTATE:
            annotate = '_ANNOTATED'
        else:
            annotate = ''

        plt.xticks(fontsize=TICKSIZE)
        plt.yticks(fontsize=TICKSIZE)

        plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_GRAYSCALE_CONSTANT_" + response + "_" + X.replace(" ",'') + "_" + Y.replace(' ','') + add_HA + annotate + '.svg', bbox_inches='tight')

    return

def plot_output_heatmap_diverging(simsDF, OUTPUT, FILEID, NORM, SCORE, SAVELOC):
    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1.0)
            plt.register_cmap(name=cm_name, data=cm)

    simsDF = maximum_absolute_scaling(simsDF)

    ascending = []
    for c in columns:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_CANCER_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=columns, ascending=ascending)

    columns.append(OUTPUT)
    if 'CANCER' in OUTPUT:
        cmap_response = 'coolwarm'
        center = 1.0
        vmax = 3
        if '_CH_' in FILEID:
            vmax = 5
        if 'AVG' in FILEID and 'VIVO' in FILEID:
            vmax = 0.35
    elif 'HEALTHY' in OUTPUT:
        cmap_response = 'coolwarm_r'
        center = 1.0
        vmax = 1.5
        if 'AVG' in FILEID and 'VIVO' in FILEID:
            vmax = 1
    elif 'TCELL' in OUTPUT:
        cmap_response = 'coolwarm_r'
        center = 1.0
        vmax = 300
        if 'AVG' in FILEID and 'VIVO' in FILEID:
            vmax = 100
    else:
        cmap_response = 'coolwarm_r'
        center = 0
        vmax = 2
    cmap.append(cmap_response)
    vmaxes.append(vmax)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0})

    for i, (s, a, c, v) in enumerate(zip(columns, ax, cmap, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=True, center=center, vmin=0, vmax=v)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_HEATMAP_DIVERGING_" + OUTPUT + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_grayscale(simsDF, OUTPUT, FILEID, NORM, SCORE, SAVELOC):
    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1.0)
            plt.register_cmap(name=cm_name, data=cm)

    simsDF = maximum_absolute_scaling(simsDF)

    ascending = []
    for c in columns:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_CANCER_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=columns, ascending=ascending)

    columns.append(OUTPUT)
    cmap_response = 'Greys'
    if 'CANCER' in OUTPUT:
        vmax = 3
        if '_CH_' in FILEID:
            vmax = 5
        if 'AVG' in FILEID and 'VIVO' in FILEID:
            vmax = 0.35
    elif 'HEALTHY' in OUTPUT:
        vmax = 1.5
        if 'AVG' in FILEID and 'VIVO' in FILEID:
            vmax = 1
    elif 'TCELL' in OUTPUT:
        vmax = 300
        if 'AVG' in FILEID and 'VIVO' in FILEID:
            vmax = 100
    else:
        vmax = 2
    cmap.append(cmap_response)
    vmaxes.append(vmax)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0})

    for i, (s, a, c, v) in enumerate(zip(columns, ax, cmap, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=True, vmin=0, vmax=v)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_HEATMAP_GRAYSCALE_" + OUTPUT + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_diverging_multiple_outputs(simsDF, OUTPUTS, SORT_SCORES, FILEID, NORM, SCORE, SAVELOC):
    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    centers = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1)
            plt.register_cmap(name=cm_name, data=cm)
            centers.append(None)

    simsDF = maximum_absolute_scaling(simsDF)

    if SORT_SCORES: sortby = ['SCORE'] + columns
    else: sortby = columns

    ascending = []
    for c in sortby:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_HEALTHY_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)


    for output in OUTPUTS:
        if 'CANCER' in output:
            cmap_response = 'coolwarm'
            center = 1.0
            vmax = 3
            if '_CH_' in FILEID:
                vmax = 5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 0.35
        elif 'HEALTHY' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
            vmax = 1.5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 1
        elif 'TCELL' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
            vmax = 300
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 100
        else:
            cmap_response = 'coolwarm_r'
            center = 0
            vmax = 2
        cmap.append(cmap_response)
        vmaxes.append(vmax)
        columns.append(output)
        centers.append(center)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0}, figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    for i, (s, a, c, m, v) in enumerate(zip(columns, ax, cmap, centers, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, center=m, vmin=0, vmax=v)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])
    add_sort = '_SORTED_SCORES' if SORT_SCORES else ''

    if '_CH_' in FILEID:
        add_score = "_" + SCORE
    else:
        add_score = ''

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_DIVERGING_ALLOUTPUTS" + add_sort + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_grayscale_multiple_outputs(simsDF, OUTPUTS, SORT_SCORES, FILEID, NORM, SCORE, SAVELOC):
    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    centers = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1)
            plt.register_cmap(name=cm_name, data=cm)
            centers.append(None)

    simsDF = maximum_absolute_scaling(simsDF)

    if SORT_SCORES: sortby = ['SCORE'] + columns
    else: sortby = columns

    ascending = []
    for c in sortby:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_HEALTHY_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    cmap_response = 'Greys'
    for output in OUTPUTS:
        if 'CANCER' in output:
            vmax = 3
            if '_CH_' in FILEID:
                vmax = 5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 0.35
        elif 'HEALTHY' in output:
            vmax = 1.5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 1
        elif 'TCELL' in output:
            vmax = 300
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 100
        else:
            vmax = 2
        cmap.append(cmap_response)
        vmaxes.append(vmax)
        columns.append(output)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0}, figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    for i, (s, a, c, v) in enumerate(zip(columns, ax, cmap, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0, vmax=v)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])
    add_sort = '_SORTED_SCORES' if SORT_SCORES else ''

    if '_CH_' in FILEID:
        add_score = "_" + SCORE
    else:
        add_score = ''

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_GRAYSCALE_ALLOUTPUTS" + add_sort + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_diverging_multiple_outputs_choose_sort(simsDF, OUTPUTS, SORT, FILEID, NORM, SCORE, SAVELOC):
    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    centers = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1)
            plt.register_cmap(name=cm_name, data=cm)
            centers.append(None)

    simsDF = maximum_absolute_scaling(simsDF)

    sortby = [SORT] + columns

    ascending = []
    for c in sortby:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_HEALTHY_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)


    for output in OUTPUTS:
        if 'CANCER' in output:
            cmap_response = 'coolwarm'
            center = 1.0
            vmax = 3
            if '_CH_' in FILEID:
                vmax = 5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 0.35
        elif 'HEALTHY' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
            vmax = 1.5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 1
        elif 'TCELL' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
            vmax = 300
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 100
        else:
            cmap_response = 'coolwarm_r'
            center = 0
            vmax = 2
        cmap.append(cmap_response)
        vmaxes.append(vmax)
        columns.append(output)
        centers.append(center)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0}, figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    for i, (s, a, c, m, v) in enumerate(zip(columns, ax, cmap, centers, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, center=m, vmin=0, vmax=v)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])
    add_sort = '_SORTED_'+SORT
    if '_CH_' in FILEID:
        add_score = "_" + SCORE
    else:
        add_score = ''

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_DIVERGING_ALLOUTPUTS" + add_sort + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_grayscale_multiple_outputs_choose_sort(simsDF, OUTPUTS, SORT, FILEID, NORM, SCORE, SAVELOC):
    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    centers = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1)
            plt.register_cmap(name=cm_name, data=cm)
            centers.append(None)

    simsDF = maximum_absolute_scaling(simsDF)

    sortby = [SORT] + columns

    ascending = []
    for c in sortby:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_HEALTHY_LIVE':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    cmap_response = 'Greys'
    for output in OUTPUTS:
        if 'CANCER' in output:
            vmax = 3
            if '_CH_' in FILEID:
                vmax = 5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 0.35
        elif 'HEALTHY' in output:
            vmax = 1.5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 1
        elif 'TCELL' in output:
            vmax = 300
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 100
        else:
            vmax = 2
        cmap.append(cmap_response)
        vmaxes.append(vmax)
        columns.append(output)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0}, figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    for i, (s, a, c, v) in enumerate(zip(columns, ax, cmap, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0, vmax=v)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])
    add_sort = '_SORTED_'+SORT
    if '_CH_' in FILEID:
        add_score = "_" + SCORE
    else:
        add_score = ''

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_GRAYSCALE_ALLOUTPUTS" + add_sort + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_with_subplots_line_multiple_outputs(simsDF, OUTPUTS, SORT_SCORES, FILEID, NORM, SCORE, SAVELOC):
    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    centers = []
    vmaxes = []

    for colordict in COLOR_DICT_HEATMAP:
        if '_CH_' not in FILEID and colordict == 'ANTIGENS HEALTHY':
            continue
        else:
            values = []
            cols = []
            cdict = COLOR_DICT_HEATMAP[colordict]
            for key in cdict:
                try:
                    values.append(float(key))
                    cols.append(cdict[key])
                except ValueError:
                    continue
            values = [val/max(values) for val in values]
            cm_name = 'cmap_' + colordict.replace(' ', '_')
            cmapdict = NonLinCdict(values, cols)
            cm = mcolors.LinearSegmentedColormap(cm_name, cmapdict)
            cmap.append(cm)
            vmaxes.append(1)
            plt.register_cmap(name=cm_name, data=cm)
            centers.append(None)

    simsDF = maximum_absolute_scaling(simsDF)

    if SORT_SCORES:
        if '_CH_' in FILEID:
            sortby = ['SCORE'] + columns
        else:
            sortby = ['Y_NORM_CANCER_LIVE'] + columns
    else: sortby = columns

    ascending = []
    for c in sortby:
        if c == 'SCORE' or c == 'CAR_AFFINITY' or c == 'Y_NORM_HEALTHY_LIVE':
            ascending.append(True)
        else:
            ascending.append(False)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)


    for output in OUTPUTS:
        if 'CANCER' in output:
            cmap_response = 'coolwarm'
            center = 1.0
            vmax = 3
            if '_CH_' in FILEID:
                vmax = 5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 0.35
        elif 'HEALTHY' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
            vmax = 1.5
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 1
        elif 'TCELL' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
            vmax = 300
            if 'AVG' in FILEID and 'VIVO' in FILEID:
                vmax = 100
        else:
            cmap_response = 'coolwarm_r'
            center = 0
            vmax = 2
        cmap.append(cmap_response)
        vmaxes.append(vmax)
        columns.append(output)
        centers.append(center)

    f, ax = plt.subplots(len(columns),1, gridspec_kw={'wspace': 0}, figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    for i, (s, a, c, m, v) in enumerate(zip(columns, ax, cmap, centers, vmaxes)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            x = np.array([simsDF[columns[i]].values])
            sns.heatmap(x, ax=a, xticklabels=[], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        elif 'SCORE' in columns[i]:
            x = [index for index in range(0, len(simsDF))]
            sns.lineplot(x, simsDF[columns[i]].values, ax=ax[i-1], color='black')
            ax[i-1].set_xlim([0,len(simsDF)-1])
            ax[i-1].set_ylim([0,max(simsDF['SCORE'])])
            ax[i-1].set_yticks([0,max(simsDF[columns[i]])])
            ax[i-1].set_yticklabels([0,max(simsDF[columns[i]])])
            a.set_visible(False)
            ax[i-1].set_ylabel('SCORE', rotation='horizontal', labelpad=15)
        elif 'CANCER' in columns[i]:
            x = [index for index in range(0, len(simsDF))]
            sns.lineplot(x, simsDF[columns[i]].values, ax=a, color='darkgray')
            a.set_xlim([0, len(simsDF)-1])
            a.set_ylim([0, max(simsDF[columns[i]])])
            a.set_yticks([0,max(simsDF[columns[i]])])
            a.set_yticklabels([0,max(simsDF[columns[i]])])
            a.set_ylabel('NORMALIZED CANCER AND HEALTHY VALUE',rotation='horizontal', labelpad=125)
            if '_C_' in FILEID:
                y = [1 for i in x]
                sns.lineplot(x, y, ax=a, color='lightgray')
        elif 'HEALTHY' in columns[i]:
            x = [index for index in range(0, len(simsDF))]
            sns.lineplot(x, simsDF[columns[i]].values, ax=ax[i-1], color='lightgray')
            if max(simsDF[columns[i]]) > max(simsDF[columns[i-1]]):
                ax[i-1].set_ylim([0, max(simsDF[columns[i]])])
                ax[i-1].set_yticks([0,max(simsDF[columns[i]])])
                ax[i-1].set_yticklabels([0, max(simsDF[columns[i]])])
        else:
            x = [index for index in range(0, len(simsDF))]
            if '_CH_' in FILEID:
                sns.lineplot(x, simsDF[columns[i]].values, ax=ax[i - 1], color='darkgray')
                ax[i-1].set_xlim([0, len(simsDF)-1])
                ax[i-1].set_ylim([0, max(simsDF[columns[i]])])
                ax[i-1].set_yticks([0, max(simsDF[columns[i]])])
                ax[i-1].set_yticklabels([0, max(simsDF[columns[i]])])
                ax[i-1].set_ylabel(columns[i], rotation='horizontal', labelpad=75)
            else:
                sns.lineplot(x, simsDF[columns[i]].values, ax=a, color='darkgray')
                a.set_xlim([0, len(simsDF)-1])
                a.set_ylim([0, max(simsDF[columns[i]])])
                a.set_yticks([0, max(simsDF[columns[i]])])
                a.set_yticklabels([0, max(simsDF[columns[i]])])
                a.set_ylabel(columns[i], rotation='horizontal', labelpad=75)
        if i < 5 and '_CH_' in FILEID:
            a.set_yticklabels([s], rotation=0)
        elif i < 4 and '_C_' in FILEID:
            a.set_yticklabels([s], rotation=0)
        a.xaxis.set_ticks([])


    if '_CH_' in FILEID:
        add_sort = '_SORTED_SCORES' if SORT_SCORES else ''
    else:
        add_sort = '_SORTED_CANCER' if SORT_SCORES else ''

    if '_CH_' in FILEID:
        add_score = "_" + SCORE
    else:
        add_score = ''
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig(SAVELOC + FILEID + "_" + NORM + add_score + "_HEATMAP_ALLOUTPUTS_SUBPLOTS" + add_sort + '.pdf', bbox_inches='tight')

    return

def plot_output_elbow_sort_output_and_color(simsDF, OUTPUT, COLOR, FILEID, NORM, SCORE, SAVELOC):

    colorDict = COLOR_DICT_HEATMAP[COLOR.replace('_', ' ')]

    sortby = [OUTPUT] + [COLOR]

    ascending = []
    for c in sortby:
        if c == 'CAR_AFFINITY':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    figScore = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figScore.add_subplot(1, 1, 1)

    for i in range(0, len(simsDF)):
        if COLOR != 'TREAT_RATIO' and COLOR != 'CAR_AFFINITY':
            key = int(simsDF.iloc[i][COLOR])
            key = str(key)
        else:
            key = str(simsDF.iloc[i][COLOR])

        if key == '0' and COLOR == 'TREAT_RATIO': key = '0.0'
        elif key == '0.5': key = '0.50'
        elif key == '1': key = '1.0'

        ax.scatter(i, simsDF.iloc[i][OUTPUT], color=colorDict[key])

    ax.set_xlabel("INDEX", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(OUTPUT + " VALUE", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(OUTPUT + " VALUE ELBOW PLOT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([0, len(simsDF)])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_ELBOW_OCSORT_" + OUTPUT.replace('_', '') + '_' + COLOR.replace('_', '') + '.svg', bbox_inches='tight')

    return

def plot_output_elbow_sort_all(simsDF, OUTPUT, COLOR, FILEID, NORM, SCORE, SAVELOC):

    colorDict = COLOR_DICT_HEATMAP[COLOR.replace('_', ' ')]

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    sortby = [OUTPUT] + columns

    ascending = []
    for c in sortby:
        if c == 'CAR_AFFINITY':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    figScore = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figScore.add_subplot(1, 1, 1)

    for i in range(0, len(simsDF)):
        if COLOR != 'TREAT_RATIO' and COLOR != 'CAR_AFFINITY':
            key = int(simsDF.iloc[i][COLOR])
            key = str(key)
        else:
            key = str(simsDF.iloc[i][COLOR])

        if key == '0' and COLOR == 'TREAT_RATIO': key = '0.0'
        elif key == '0.5': key = '0.50'
        elif key == '1': key = '1.0'

        ax.scatter(i, simsDF.iloc[i][OUTPUT], color=colorDict[key])

    ax.set_xlabel("INDEX", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(OUTPUT + "\nVALUE", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(OUTPUT + "\nVALUE ELBOW PLOT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([0, len(simsDF)])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')

    if '_CH_' in FILEID:
        add = '_' + SCORE
    else:
        add = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_ELBOW_ALLSORT_" + OUTPUT.replace('_', '') + '_' + COLOR.replace('_', '') + '.svg', bbox_inches='tight')

def plot_output_elbow_CH_sort_cancer_and_color(simsDF, COLOR, FILEID, NORM, SCORE, SAVELOC):

    colorDict = COLOR_DICT_HEATMAP[COLOR.replace('_', ' ')]

    sortby = ['Y_NORM_CANCER_LIVE'] + [COLOR]

    ascending = []
    for c in sortby:
        if c == 'CAR_AFFINITY':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    figScore = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figScore.add_subplot(1, 1, 1)

    for i in range(0, len(simsDF)):
        if COLOR != 'TREAT_RATIO' and COLOR != 'CAR_AFFINITY':
            key = int(simsDF.iloc[i][COLOR])
            key = str(key)
        else:
            key = str(simsDF.iloc[i][COLOR])

        if key == '0' and COLOR == 'TREAT_RATIO': key = '0.0'
        elif key == '0.5': key = '0.50'
        elif key == '1': key = '1.0'

        ax.scatter(i, simsDF.iloc[i]['Y_NORM_CANCER_LIVE'], color=colorDict[key], marker='o')
        ax.scatter(i, simsDF.iloc[i]['Y_NORM_HEALTHY_LIVE'], color=colorDict[key], marker='s')

    ax.set_xlabel("INDEX", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("NORMALIZED CANCER AND\nHEALTHY LIVE FRACTION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title("NORMALIZED CANCER AND\nHEALTHY LIVE FRACTION\nELBOW PLOT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([0, len(simsDF)])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')

    if '_CH_' in FILEID:
        add = '_' + SCORE
    else:
        add = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_ELBOW_OCSORT_CH_CANCER_" + COLOR.replace('_', '') + '.svg', bbox_inches='tight')

    return

def plot_output_elbow_CH_sort_healthy_and_color(simsDF, COLOR, FILEID, NORM, SCORE, SAVELOC):

    colorDict = COLOR_DICT_HEATMAP[COLOR.replace('_', ' ')]

    sortby = ['Y_NORM_HEALTHY_LIVE'] + [COLOR]

    ascending = []
    for c in sortby:
        if c == 'CAR_AFFINITY':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    figScore = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figScore.add_subplot(1, 1, 1)

    for i in range(0, len(simsDF)):
        if COLOR != 'TREAT_RATIO' and COLOR != 'CAR_AFFINITY':
            key = int(simsDF.iloc[i][COLOR])
            key = str(key)
        else:
            key = str(simsDF.iloc[i][COLOR])

        if key == '0' and COLOR == 'TREAT_RATIO': key = '0.0'
        elif key == '0.5': key = '0.50'
        elif key == '1': key = '1.0'

        ax.scatter(i, simsDF.iloc[i]['Y_NORM_CANCER_LIVE'], color=colorDict[key], marker='o')
        ax.scatter(i, simsDF.iloc[i]['Y_NORM_HEALTHY_LIVE'], color=colorDict[key], marker='s')

    ax.set_xlabel("INDEX", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("NORMALIZED CANCER AND\nHEALTHY LIVE FRACTION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title("NORMALIZED CANCER AND\nHEALTHY LIVE FRACTION\nELBOW PLOT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([0, len(simsDF)])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_ELBOW_OCSORT_CH_HEALTHY_" + COLOR.replace('_', '') + '.svg', bbox_inches='tight')

    return

def plot_output_elbow_CH_sort_all(simsDF, COLOR, FILEID, NORM, SCORE, SAVELOC):

    colorDict = COLOR_DICT_HEATMAP[COLOR.replace('_', ' ')]

    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    sortby = ['Y_NORM_HEALTHY_LIVE', 'Y_NORM_CANCER_LIVE'] + columns

    ascending = []
    for c in sortby:
        if c == 'CAR_AFFINITY':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)

    figScore = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figScore.add_subplot(1, 1, 1)

    for i in range(0, len(simsDF)):
        if COLOR != 'TREAT_RATIO' and COLOR != 'CAR_AFFINITY':
            key = int(simsDF.iloc[i][COLOR])
            key = str(key)
        else:
            key = str(simsDF.iloc[i][COLOR])

        if key == '0' and COLOR == 'TREAT_RATIO': key = '0.0'
        elif key == '0.5': key = '0.50'
        elif key == '1': key = '1.0'

        ax.scatter(i, simsDF.iloc[i]['Y_NORM_CANCER_LIVE'], color=colorDict[key], marker='o')
        ax.scatter(i, simsDF.iloc[i]['Y_NORM_HEALTHY_LIVE'], color=colorDict[key], marker='s')

    ax.set_xlabel("INDEX", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("NORMALIZED CANCER AND\nHEALTHY LIVE FRACTION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title("NORMALIZED CANCER AND\nHEALTHY LIVE FRACTION\nELBOW PLOT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([0, len(simsDF)])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')

    if '_CH_' in FILEID:
        add = "_" + SCORE
    else:
        add = ''

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    plt.savefig(SAVELOC + FILEID + "_" + NORM + add + "_ELBOW_ALLSORT_CH_" + COLOR.replace('_', '') + '.svg', bbox_inches='tight')

# ------------- MAIN STATS FUNCTIONS -----------------

def stats_data(simsDF, NORM, SCORE, AVG, PARTIAL, FILEID, SAVELOC):

    filesplit = FILEID.split('_')
    # Count Xs in filesplit
    X = 0
    for i in range(5, 10):
        if filesplit[i] == 'X': X += 1

    print('\t\t' + 'Cleaning up data...')
    simsDF = simsDF.replace({'TREAT RATIO': TREAT_RATIO_DICT})

    simsDF = simsDF.rename(columns={'TREAT RATIO': 'TREAT_RATIO', 'CAR AFFINITY': 'CAR_AFFINITY',
                                    'ANTIGENS CANCER': 'ANTIGENS_CANCER', 'ANTIGENS HEALTHY': 'ANTIGENS_HEALTHY'})
    simsDF['Y_NORM_CANCER_LIVE'] = None
    simsDF['Y_NORM_HEALTHY_LIVE'] = None
    simsDF['Y_NORM_TCELL_LIVE'] = None

    untreatedDF = simsDF[simsDF['DOSE'] == 0]

    simsDF = simsDF[simsDF['DOSE'] != 0]

    # Normalize data by INITIALIZATION quantity
    if NORM == 'INIT':

        if 'VITRO' in FILEID:
            for i in range(0,len(simsDF)):
                simsDF.at[i, 'Y_NORM_CANCER_LIVE'] = float(simsDF.iloc[i]['CANCER LIVE'][-1]/simsDF.iloc[i]['CANCER LIVE'][0])
                simsDF.at[i, 'Y_NORM_TCELL_LIVE'] = float(simsDF.iloc[i]['T-CELL LIVE'][-1]/simsDF.iloc[i]['DOSE'])
                if '_CH_' in FILEID:
                    simsDF.at[i, 'Y_NORM_HEALTHY_LIVE'] = float(simsDF.iloc[i]['HEALTHY LIVE'][-1]/simsDF.iloc[i]['HEALTHY LIVE'][0])

        else:
            for i in range(0,len(simsDF)):
                simsDF.at[i, 'Y_NORM_CANCER_LIVE'] = float(simsDF.iloc[i]['CANCER LIVE'][-1]/simsDF.iloc[i]['CANCER LIVE'][44])
                simsDF.at[i, 'Y_NORM_TCELL_LIVE'] = float(simsDF.iloc[i]['T-CELL LIVE'][-1]/simsDF.iloc[i]['DOSE'])
                if '_CH_' in FILEID:
                    simsDF.at[i, 'Y_NORM_HEALTHY_LIVE'] = float(simsDF.iloc[i]['HEALTHY LIVE'][-1]/simsDF.iloc[i]['HEALTHY LIVE'][44])

    # Normalize data by UNTREATED quantity
    else:

        if 'VITRO' in FILEID:
            for i in range(0,len(simsDF)):

                for j in range(0,len(untreatedDF)):
                    if untreatedDF.iloc[j]['SEED'] == simsDF.iloc[i]['SEED']:
                        cancerUntreated = untreatedDF.iloc[j]['CANCER LIVE'][-1]
                        if '_CH_' in FILEID:
                            healthyUntreated = untreatedDF.iloc[j]['HEALTHY LIVE'][-1]
                        break

                simsDF.at[i, 'Y_NORM_CANCER_LIVE'] = float(simsDF.iloc[i]['CANCER LIVE'][-1]/cancerUntreated)
                simsDF.at[i, 'Y_NORM_TCELL_LIVE'] = float(simsDF.iloc[i]['T-CELL LIVE'][-1]/simsDF.iloc[i]['DOSE'])
                if '_CH_' in FILEID:
                    simsDF.at[i, 'Y_NORM_HEALTHY_LIVE'] = float(simsDF.iloc[i]['HEALTHY LIVE'][-1]/healthyUntreated)

        else:
            for i in range(0,len(simsDF)):

                for j in range(0,len(untreatedDF)):
                    if untreatedDF.iloc[j]['SEED'] == simsDF.iloc[i]['SEED']:
                        cancerUntreated = untreatedDF.iloc[j]['CANCER LIVE'][-1]
                        if '_CH_' in FILEID:
                            healthyUntreated = untreatedDF.iloc[j]['HEALTHY LIVE'][-1]
                        break

                simsDF.at[i, 'Y_NORM_CANCER_LIVE'] = float(simsDF.iloc[i]['CANCER LIVE'][-1]/cancerUntreated)
                simsDF.at[i, 'Y_NORM_TCELL_LIVE'] = float(simsDF.iloc[i]['T-CELL LIVE'][-1]/simsDF.iloc[i]['DOSE'])
                if '_CH_' in FILEID:
                    simsDF.at[i, 'Y_NORM_HEALTHY_LIVE'] = float(simsDF.iloc[i]['HEALTHY LIVE'][-1]/healthyUntreated)

    # Create anova dataframe for use in rest of code (not just anova functions)
    simsDFanova = pd.DataFrame({'DOSE': pd.Series([],dtype='float'),
                                'TREAT_RATIO': pd.Series([],dtype='float'),
                                'CAR_AFFINITY': pd.Series([],dtype='float'),
                                'ANTIGENS_CANCER': pd.Series([],dtype='float'),
                                'ANTIGENS_HEALTHY': pd.Series([],dtype='float'),
                                'Y_NORM_CANCER_LIVE': pd.Series([],dtype='float'),
                                'Y_NORM_HEALTHY_LIVE': pd.Series([], dtype='float'),
                                'Y_NORM_TCELL_LIVE': pd.Series([],dtype='float') })

    # Fill in simsDFanova
    if '_C_' in FILEID:
        for i in range(0, len(simsDF)):
            simDict = {}
            for column in ['DOSE','TREAT_RATIO','CAR_AFFINITY','ANTIGENS_CANCER','Y_NORM_CANCER_LIVE','Y_NORM_TCELL_LIVE']:
                simDict[column] = float(simsDF.at[i, column])
            simsDFanova = simsDFanova.append(simDict, ignore_index=True)
    else:
        for i in range(0, len(simsDF)):
            simDict = {}
            for column in ['DOSE','TREAT_RATIO','CAR_AFFINITY','ANTIGENS_CANCER','ANTIGENS_HEALTHY','Y_NORM_CANCER_LIVE','Y_NORM_HEALTHY_LIVE','Y_NORM_TCELL_LIVE']:
                simDict[column] = float(simsDF.at[i, column])
            simsDFanova = simsDFanova.append(simDict, ignore_index=True)

    # Create list of combinations of features and list of responses
    if '_C_' in FILEID:
        comb = list(combinations(['DOSE','TREAT RATIO','CAR AFFINITY','ANTIGENS CANCER'], 2))
        response_list = ['Y_NORM_CANCER_LIVE', 'Y_NORM_TCELL_LIVE']
    else:
        comb = list(combinations(['DOSE','TREAT RATIO','CAR AFFINITY','ANTIGENS CANCER','ANTIGENS HEALTHY'], 2))
        response_list = ['Y_NORM_CANCER_LIVE', 'Y_NORM_HEALTHY_LIVE', 'Y_NORM_TCELL_LIVE', 'SCORE']

    # Create SCORE if both cancer and healthy cells present
    if '_CH_' in FILEID:
        scorelist = []
        for i in range(0, len(simsDFanova)):
            cancer = simsDFanova.iloc[i]['Y_NORM_CANCER_LIVE']
            healthy = simsDFanova.iloc[i]['Y_NORM_HEALTHY_LIVE']

            if SCORE == 'SUM':
                # calculate cancer score
                delta = 1
                if cancer > 1 or healthy < SCORE_MIN_HEALTHY_THRESHOLD:
                    delta = 0

                score = (1 - cancer + healthy) * delta
            else:
                if cancer == 0:
                    print(simsDF.iloc[i]['TUMOR ID'] + str(simsDF.iloc[i]['SEED']) + ' CANCER 0')
                if healthy == 0:
                    print(simsDF.iloc[i]['TUMOR ID'] + str(simsDF.iloc[i]['SEED']) + ' HEALTHY 0')
                score = 1

            scorelist.append(score)
        simsDFanova['SCORE'] = scorelist

    # Create color reference lists
    colorlist = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID: colorlist.append('ANTIGENS_HEALTHY')

    # Conduct analyses for data that is not averaged across replicates
    if not AVG:

        # Do ANOVA only if full combinatoral dataset present
        # if (('_CH_' in FILEID and X == 5) or ('_C_' in FILEID and X == 4)) and not PARTIAL:
        #     print('\t\t' + 'Running ANOVA.')
        #     anova(simsDFanova, FILEID, NORM, SCORE, SAVELOC)

        # Do two axis heatmaps only if full combinatoral dataset present
        # if (('_CH_' in FILEID and X == 5) or ('_C_' in FILEID and X == 4)) and not PARTIAL:
        #     print('\t\tPlotting output two-axis heatmaps.')
        #     for response in response_list:
        #         for combo in comb:
        #             # Make diverging scale plots
        #             plot_two_axis_heatmaps_diverging_all(simsDFanova, combo[0], combo[1], response, True, FILEID, NORM, SCORE, SAVELOC)
        #             plot_two_axis_heatmaps_diverging_constant(simsDFanova, combo[0], combo[1], response, True, FILEID, NORM, SCORE, SAVELOC)
        #             plot_two_axis_heatmaps_diverging_all(simsDFanova, combo[0], combo[1], response, False, FILEID, NORM, SCORE, SAVELOC)
        #             plot_two_axis_heatmaps_diverging_constant(simsDFanova, combo[0], combo[1], response, False, FILEID, NORM, SCORE, SAVELOC)
        #             plt.close('all')
        #             # Make grayscale plots
        #             plot_two_axis_heatmaps_grayscale_all(simsDFanova, combo[0], combo[1], response, True, FILEID, NORM, SCORE, SAVELOC)
        #             plot_two_axis_heatmaps_grayscale_constant(simsDFanova, combo[0], combo[1], response, True, FILEID, NORM, SCORE, SAVELOC)
        #             plot_two_axis_heatmaps_grayscale_all(simsDFanova, combo[0], combo[1], response, False, FILEID, NORM, SCORE, SAVELOC)
        #             plot_two_axis_heatmaps_grayscale_constant(simsDFanova, combo[0], combo[1], response, False, FILEID, NORM, SCORE, SAVELOC)

        # Plot output heatmaps (unsorted)
        print('\t\tPlotting output all-axis heatmaps.')
        # for response in response_list:
        #     plot_output_heatmap_diverging(simsDFanova, response, FILEID, NORM, SCORE, SAVELOC)
        #     plot_output_heatmap_grayscale(simsDFanova, response, FILEID, NORM, SCORE, SAVELOC)
        # if '_CH_' in FILEID:
            # plot_output_heatmap_diverging_multiple_outputs(simsDFanova, response_list, False, FILEID, NORM, SCORE, SAVELOC)
            # plot_output_heatmap_diverging_multiple_outputs(simsDFanova, response_list, True, FILEID, NORM, SCORE, SAVELOC)
            # plot_output_heatmap_grayscale_multiple_outputs(simsDFanova, response_list, False, FILEID, NORM, SCORE, SAVELOC)
            # plot_output_heatmap_grayscale_multiple_outputs(simsDFanova, response_list, True, FILEID, NORM, SCORE, SAVELOC)
        plot_output_heatmap_with_subplots_line_multiple_outputs(simsDFanova, response_list, True, FILEID, NORM, SCORE, SAVELOC)
        plt.close('all')

        # Plot output heatmaps (sorted)
        # print('\t\tPlotting output all-axis heatmaps with chosen sort.')
        # plot_output_heatmap_diverging_multiple_outputs_choose_sort(simsDFanova, response_list, 'Y_NORM_CANCER_LIVE', FILEID, NORM, SCORE, SAVELOC)
        # plot_output_heatmap_grayscale_multiple_outputs_choose_sort(simsDFanova, response_list, 'Y_NORM_CANCER_LIVE', FILEID, NORM, SCORE, SAVELOC)
        # if '_CH_' in FILEID:
        #     plot_output_heatmap_diverging_multiple_outputs_choose_sort(simsDFanova, response_list, 'Y_NORM_HEALTHY_LIVE', FILEID, NORM, SCORE, SAVELOC)
        #     plot_output_heatmap_grayscale_multiple_outputs_choose_sort(simsDFanova, response_list, 'Y_NORM_HEALTHY_LIVE', FILEID, NORM, SCORE, SAVELOC)
        '''       
        # Make elbow plots
        print('\t\tPlotting output elbow plots.')
        for response in response_list:
            for color in colorlist:
                plot_output_elbow_sort_output_and_color(simsDFanova, response, color, FILEID, NORM, SCORE, SAVELOC)
                plot_output_elbow_sort_all(simsDFanova, response, color, FILEID, NORM, SCORE, SAVELOC)
                plot_output_elbow_CH_sort_cancer_and_color(simsDFanova, color, FILEID, NORM, SCORE, SAVELOC)
                plot_output_elbow_CH_sort_healthy_and_color(simsDFanova, color, FILEID, NORM, SCORE, SAVELOC)
                plot_output_elbow_CH_sort_all(simsDFanova, color, FILEID, NORM, SCORE, SAVELOC)
                plt.close('all')

        # Create excel files sorted by output
        print('\t\tCreating Excel files for all threshold passing simulations.')
        for response in response_list:
            if response != 'Y_NORM_TCELL_LIVE':
                feature_analysis(simsDFanova, response, FILEID, NORM, SCORE, SAVELOC)
                save_sorted_df_xlsx(simsDFanova, response, FILEID, NORM, SCORE, SAVELOC)
                save_sorted_df_xlsx_all(simsDFanova, response, FILEID, NORM, SCORE, SAVELOC)
    '''
    # Conduct analyses for data that is averaged across replicates
    else:

        # Average data across conditions
        print('\t\tAveraging datta across conditions.')
        simsDFavg = average_conditions(simsDFanova, FILEID)
        FILEIDAVG = FILEID + '_AVG'

        # Plot output heatmaps (unsorted)
        print('\t\t\tPlotting output all-axis heatmaps.')
        # for response in response_list:
        #     plot_output_heatmap_diverging(simsDFavg, response, FILEIDAVG, NORM, SCORE, SAVELOC)
        #     plot_output_heatmap_grayscale(simsDFavg, response, FILEIDAVG, NORM, SCORE, SAVELOC)
        #     plot_output_heatmap_with_subplots_line_multiple_outputs(simsDFavg, response_list, True, FILEIDAVG, NORM, SCORE, SAVELOC)
        # if '_CH_' in FILEID:
        #     plot_output_heatmap_diverging_multiple_outputs(simsDFavg, response_list, False, FILEIDAVG, NORM, SCORE, SAVELOC)
        #     plot_output_heatmap_diverging_multiple_outputs(simsDFavg, response_list, True, FILEIDAVG, NORM, SCORE, SAVELOC)
        #     plot_output_heatmap_grayscale_multiple_outputs(simsDFavg, response_list, False, FILEIDAVG, NORM, SCORE, SAVELOC)
        #     plot_output_heatmap_grayscale_multiple_outputs(simsDFavg, response_list, True, FILEIDAVG, NORM, SCORE, SAVELOC)
        plot_output_heatmap_with_subplots_line_multiple_outputs(simsDFavg, response_list, True, FILEIDAVG, NORM, SCORE, SAVELOC)
        plt.close('all')

        # Plot output heatmaps (sorted)
        # print('\t\t\tPlotting output all-axis heatmaps with chosen sort.')
        # plot_output_heatmap_diverging_multiple_outputs_choose_sort(simsDFavg, response_list, 'Y_NORM_CANCER_LIVE', FILEIDAVG, NORM, SCORE, SAVELOC)
        # plot_output_heatmap_grayscale_multiple_outputs_choose_sort(simsDFavg, response_list, 'Y_NORM_CANCER_LIVE', FILEIDAVG, NORM, SCORE, SAVELOC)
        # if '_CH_' in FILEID:
        #     plot_output_heatmap_diverging_multiple_outputs_choose_sort(simsDFavg, response_list, 'Y_NORM_HEALTHY_LIVE', FILEIDAVG, NORM, SCORE, SAVELOC)
        #     plot_output_heatmap_grayscale_multiple_outputs_choose_sort(simsDFavg, response_list, 'Y_NORM_HEALTHY_LIVE', FILEIDAVG, NORM, SCORE, SAVELOC)

        # # Plot elbow plots
        # print('\t\t\tPlotting output elbow plots.')
        # for response in response_list:
        #     for color in colorlist:
        #         plot_output_elbow_sort_output_and_color(simsDFavg, response, color, FILEIDAVG, NORM, SCORE, SAVELOC)
        #         plot_output_elbow_sort_all(simsDFavg, response, color, FILEIDAVG, NORM, SCORE, SAVELOC)
        #         plot_output_elbow_CH_sort_cancer_and_color(simsDFavg, color, FILEIDAVG, NORM, SCORE, SAVELOC)
        #         plot_output_elbow_CH_sort_healthy_and_color(simsDFavg, color, FILEIDAVG, NORM, SCORE, SAVELOC)
        #         plot_output_elbow_CH_sort_all(simsDFavg, color, FILEIDAVG, NORM, SCORE, SAVELOC)
        #         plt.close('all')
        #
        # # Create excel files sorted by output
        # print('\t\t\tCreating Excel files for all threshold passing simulations.')
        # for response in response_list:
        #     if response != 'Y_NORM_TCELL_LIVE':
        #         feature_analysis(simsDFavg, response, FILEIDAVG, NORM, SCORE, SAVELOC)
        #         save_sorted_df_xlsx(simsDFavg, response, FILEIDAVG, NORM, SCORE, SAVELOC)
        #         save_sorted_df_xlsx_all(simsDFavg, response, FILEIDAVG, NORM, SCORE, SAVELOC)
        #     plt.close('all')

    return

# ---------------- MAIN FUNCTION ---------------------

if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()

    # Get files
    PKLFILES = get_pkl_files(args.files)

    print("Running stats for the following files:")

    for file in PKLFILES:

        if 'VITRO' in file:
            fileName = re.sub('.*VITRO', 'VITRO', file)

        else:
            fileName = re.sub('.*VIVO', 'VIVO', file)

        print('\t' + fileName)

        FILEID = ''

        if 'ANALYZED' in fileName:
            FILEID = fileName.replace('_ANALYZED.pkl', '')

        if 'ENVIRONMENT' in fileName:
            FILEID = fileName.replace('_ENVIRONMENT.pkl', '')

        if 'SPATIAL' in fileName:
            FILEID = fileName.replace('_SPATIAL.pkl', '')

        if 'LYSED' in fileName:
            FILEID = fileName.replace('_LYSED.pkl', '')

        with open(file, 'rb') as f:
            simsDF = pickle.load(f)

        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', -1)
        #print(simsDF)

        if args.norm == 'INIT':
            NORM = 'INIT'
        elif args.norm == 'UNTREATED':
            NORM = 'UNTREATED'
        else:
            NORM = 'INIT'

        if args.score == 'SUM':
            SCORE = 'SUM'
        elif args.score == 'RATIO':
            SCORE = 'RATIO'
        else:
            SCORE = 'SUM'

        stats_data(simsDF, NORM, SCORE, args.average, args.partial, FILEID, args.saveLoc)