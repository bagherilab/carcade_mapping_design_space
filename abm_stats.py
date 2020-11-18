import ABM
import os
import pickle
import re
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
from colour import Color
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
    [--color COLOR]
        Feature by which to color data by. If X given, will color by axis along which data varies..
        If two featuers vary, default will be used. (default: CAR AFFINITY)
    [--marker MARKER]
        Feature by which to change marker of data by (default: TREAT RATIO)
    [--saveLoc SAVELOC]
        Location of where to save file, default will save here
'''

# -------------- PARSING AND INPUT/OUTPUT FUNCTIONS -------------

def get_parser():

    # Setup argument parser.
    parser = ArgumentParser(description="Plot ABM data from dataframe")
    parser.add_argument(dest="files", help="Path to .pkl file or directory")
    parser.add_argument("--color", default="CAR AFFINITY", dest="color",
                        help="Feature by which to color data by (default: CAR AFFINITY)")
    parser.add_argument("--marker", default="TREAT RATIO", dest="marker",
                        help="Feature by which to change marker of data by (default: TREAT RATIO)")
    parser.add_argument("--saveLoc", default="", dest="saveLoc",
                        help="Location of where to save file, default will save here")

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

# Set AFFINITY color scale
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
acolor1 = Color(colors['peachpuff'])
ACOLORS = list(acolor1.range_to((Color(colors['saddlebrown'])), 5))
for i in range(0, len(ACOLORS)):
    ACOLORS[i] = Color(ACOLORS[i]).rgb

# Set TREAT RATIO color scale
trcolor1 = Color(colors['lightcoral'])
TRCOLORS = list(trcolor1.range_to((Color(colors['maroon'])), 5))
for i in range(0, len(TRCOLORS)):
    TRCOLORS[i] = Color(TRCOLORS[i]).rgb

doseColorDict = {
    "0": 'black',
    "250": 'lightgreen',
    "500": "green",
    "1000": "darkgreen"
}

trColorDict = {
    "NA": 'black',
    "0:100": 'pink',
    "25:75": 'lightcoral',
    "50:50": 'crimson',
    "75:25": 'indianred',
    "100:0": 'maroon',
}

trColorDictHeatMap = {
    "0.0": 'pink',
    "0.25": 'lightcoral',
    "0.50": 'crimson',
    "0.75": 'indianred',
    "1.0": 'maroon',
}

affinityColorDict = {
    "NA": 'black',
    "0.0": ACOLORS[3],
    "1e-09": ACOLORS[3],
    "1e-08": ACOLORS[2],
    "1e-07": ACOLORS[1],
    "1e-06": ACOLORS[0],
}

acColorDict = {
    "0": 'black',
    "100": 'skyblue',
    "500": "dodgerblue",
    "1000": "cornflowerblue",
    "5000": "royalblue",
    "10000": "navy"
}

ahColorDict = {
    "CONTROL": "black",
    "0": 'plum',
    "100": 'blueviolet'
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

CONSTANTS = {
        "DOSE": 500,
        "TREAT RATIO": 0.5,
        "CAR AFFINITY": 1e-7,
        "ANTIGENS CANCER": 1000,
        "ANTIGENS HEALTHY": [0, 100]
        }

FIG_SIZE_Y = 7
FIG_SIZE_X = 7

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

# ------------- STATS FUNCTIONS -----------------

def anova(simsDF, FILEID, SAVELOC):

    print('ANOVA FOR CANCER CELLS.')
    #model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
    model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + '
                       'C(DOSE):C(TREAT_RATIO) + C(DOSE): C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) +'
                       'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) +'
                       'C(CAR_AFFINITY):C(ANTIGENS_CANCER)', data=simsDF).fit()
    anova_CANCER = sm.stats.anova_lm(model_CANCER, typ=2)
    print(anova_CANCER)

    if '_CH_' in FILEID:
        print('ANOVA FOR HEALTHY CELLS.')
        # model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_HEALTHY = ols('Y_NORM_HEALTHY_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + '
                           'C(DOSE):C(TREAT_RATIO) + C(DOSE): C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) +'
                           'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) +'
                           'C(CAR_AFFINITY):C(ANTIGENS_CANCER)', data=simsDF).fit()
        anova_HEALTHY = sm.stats.anova_lm(model_HEALTHY, typ=2)
        print(anova_HEALTHY)

    print('ANOVA FOR T-CELLS.')
    # model_TCELL = ols('Y_NORM_TCELL_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
    model_TCELL = ols('Y_NORM_TCELL_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + '
                       'C(DOSE):C(TREAT_RATIO) + C(DOSE): C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) +'
                       'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) +'
                       'C(CAR_AFFINITY):C(ANTIGENS_CANCER)', data=simsDF).fit()
    anova_TCELL = sm.stats.anova_lm(model_TCELL, typ=2)
    print(anova_TCELL)

    return

def plot_two_axis_heatmaps_all(simsDF, X, Y, COLOR, FILEID, SAVELOC):

    X_AXIS = AXES_SETS[X]
    Y_AXIS = AXES_SETS[Y]

    if Y == 'ANTIGENS CANCER' or Y == 'CAR AFFINITY':
        Y_AXIS = list(reversed(Y_AXIS))

    response = 'CANCER'
    cmap = 'coolwarm'
    center = 1.0
    vmax = 3.0
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
    ax = sns.heatmap(heatmapDF, cmap=cmap, annot=True, square=True,
                     xticklabels=X_AXIS, yticklabels=Y_AXIS, center=center,
                     cbar=True, cbar_kws={"shrink": 0.82}, vmax=vmax, vmin=vmin)
    ax.set_ylabel(Y, fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
    ax.set_xlabel(X, fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
    ax.set_title(COLOR + " RESPONSE - ALL", fontname='Arial', fontweight='bold', fontsize=12, pad=5)
    plt.savefig(SAVELOC + FILEID + "_HEATMAP_ALL_" + response + "_" + X.replace(" ",'') + "_" + Y.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_two_axis_heatmaps_constant(simsDF, X, Y, COLOR, FILEID, SAVELOC):

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
        ax = sns.heatmap(heatmapDF, cmap=cmap, annot=True, square=True,
                         xticklabels=X_AXIS, yticklabels=Y_AXIS, center=center,
                         cbar=True, cbar_kws={"shrink": 0.82}, vmax=vmax, vmin=vmin)
        ax.set_ylabel(Y, fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
        ax.set_xlabel(X, fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
        ax.set_title(COLOR + " RESPONSE - CONSTANT", fontname='Arial', fontweight='bold', fontsize=12, pad=5)
        add = ''
        if X != 'ANTIGENS HEALTHY' and Y != 'ANTIGENS HEALTHY' and '_CH_' in FILEID:
            add = '_HA_' + str(CONSTANTS['ANTIGENS HEALTHY'][d])
        ax.set_title(COLOR + " RESPONSE - CONSTANT" + add.replace('_',''), fontname='Arial', fontweight='bold', fontsize=12, pad=5)
        plt.savefig(SAVELOC + FILEID + "_HEATMAP_CONSTANT_" + response + "_" + X.replace(" ",'') + "_" + Y.replace(' ','') + add + '.svg', bbox_inches='tight')

    return

def plot_output_heatmap(simsDF, OUTPUT, FILEID, SAVELOC):
    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []

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
            plt.register_cmap(name=cm_name, data=cm)

    simsDF = maximum_absolute_scaling(simsDF)

    ascending = []
    for c in columns:
        if c == 'SCORE' or c == 'CAR_AFFINITY':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=columns, ascending=ascending)

    columns.append(OUTPUT)

    if 'CANCER' in OUTPUT:
        cmap_response = 'coolwarm'
        center = 1.0
    elif 'HEALTHY' in OUTPUT:
        cmap_response = 'coolwarm_r'
        center = 1.0
    elif 'T_CELL' in OUTPUT:
        cmap_response = 'coolwarm'
        center = 1.0
    else:
        cmap_response = 'coolwarm_r'
        center = 0
    cmap.append(cmap_response)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0})

    for i, (s, a, c) in enumerate(zip(columns, ax, cmap)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, center=center)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])

    plt.savefig(SAVELOC + FILEID + "_HEATMAP_" + OUTPUT + '.pdf', bbox_inches='tight')

    return

def plot_output_heatmap_multiple_outputs(simsDF, OUTPUTS, FILEID, SORT_SCORES, SAVELOC):
    columns = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID:
        columns.append('ANTIGENS_HEALTHY')

    cmap = []
    centers = []

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
            plt.register_cmap(name=cm_name, data=cm)
            centers.append(None)

    simsDF = maximum_absolute_scaling(simsDF)

    if SORT_SCORES: sortby = ['SCORE'] + columns
    else: sortby = columns

    ascending = []
    for c in sortby:
        if c == 'SCORE' or c == 'CAR_AFFINITY':
            ascending.append(False)
        else:
            ascending.append(True)

    simsDF = simsDF.sort_values(by=sortby, ascending=ascending)


    for output in OUTPUTS:
        if 'CANCER' in output:
            cmap_response = 'coolwarm'
            center = 1.0
        elif 'HEALTHY' in output:
            cmap_response = 'coolwarm_r'
            center = 1.0
        elif 'T_CELL' in output:
            cmap_response = 'coolwarm'
            center = 1.0
        else:
            cmap_response = 'coolwarm_r'
            center = 0
        cmap.append(cmap_response)
        columns.append(output)
        centers.append(center)

    f, ax = plt.subplots(1, len(columns), gridspec_kw={'wspace': 0}, figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    for i, (s, a, c, m) in enumerate(zip(columns, ax, cmap, centers)):
        if 'NORM' not in columns[i] and 'SCORE' not in columns[i]:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, vmin=0.0, vmax=1.0)
        else:
            sns.heatmap(np.array([simsDF[columns[i]].values]).T, ax=a, xticklabels=[s], annot=False, cmap=c, cbar=False, center=m)
        ax[i].set_xticklabels([s], rotation=90)
        if i > 0:
            a.yaxis.set_ticks([])
    add = '_SORTED_SCORES' if SORT_SCORES else ''
    plt.savefig(SAVELOC + FILEID + "_HEATMAP_ALLOUTPUTS" + add + '.pdf', bbox_inches='tight')

    return

def plot_output_elbow_sort_output_and_color(simsDF, OUTPUT, COLOR, FILEID, SAVELOC):

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

    ax.set_xlabel("INDEX", fontname='Arial', fontweight='bold', fontsize=14, labelpad=5)
    ax.set_ylabel(OUTPUT + " VALUE", fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
    ax.set_title(OUTPUT + " VALUE ELBOW PLOT", fontname='Arial', fontweight='bold', fontsize=12, pad=5)
    ax.set_xlim([0, len(simsDF)])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')

    plt.savefig(SAVELOC + FILEID + "_ELBOW_OCSORT_" + OUTPUT.replace('_', '') + '_' + COLOR.replace('_', '') + '.svg', bbox_inches='tight')

    return

def plot_output_elbow_sort_all(simsDF, OUTPUT, COLOR, FILEID, SAVELOC):

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

    ax.set_xlabel("INDEX", fontname='Arial', fontweight='bold', fontsize=14, labelpad=5)
    ax.set_ylabel(OUTPUT + " VALUE", fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
    ax.set_title(OUTPUT + " VALUE ELBOW PLOT", fontname='Arial', fontweight='bold', fontsize=12, pad=5)
    ax.set_xlim([0, len(simsDF)])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')

    plt.savefig(SAVELOC + FILEID + "_ELBOW_ALLSORT_" + OUTPUT.replace('_', '') + '_' + COLOR.replace('_', '') + '.svg', bbox_inches='tight')

# ------------- MAIN STATS FUNCTIONS -----------------

def stats_data(simsDF, COLOR, MARKER, FILEID, SAVELOC):

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

    simsDF = simsDF[simsDF['DOSE'] != 0]

    # TO DO: Convert dataframe to save final value/time point of list
    if 'VITRO' in FILEID:
        for i in range(0,len(simsDF)):
            simsDF.at[i, 'Y_NORM_CANCER_LIVE'] = float(simsDF.iloc[i]['CANCER LIVE'][-1]/simsDF.iloc[i]['CANCER LIVE'][0])
            simsDF.at[i, 'Y_NORM_TCELL_LIVE'] = float(simsDF.iloc[i]['T-CELL LIVE'][-1]/simsDF.iloc[i]['DOSE'])
            if '_CH_' in FILEID:
                simsDF.at[i, 'Y_NORM_HEALTHY_LIVE'] = float(simsDF.iloc[i]['HEALTHY LIVE'][-1]/simsDF.iloc[i]['HEALTHY LIVE'][0])

    simsDFanova = pd.DataFrame({'DOSE': pd.Series([],dtype='float'),
                                'TREAT_RATIO': pd.Series([],dtype='float'),
                                'CAR_AFFINITY': pd.Series([],dtype='float'),
                                'ANTIGENS_CANCER': pd.Series([],dtype='float'),
                                'ANTIGENS_HEALTHY': pd.Series([],dtype='float'),
                                'Y_NORM_CANCER_LIVE': pd.Series([],dtype='float'),
                                'Y_NORM_HEALTHY_LIVE': pd.Series([], dtype='float'),
                                'Y_NORM_TCELL_LIVE': pd.Series([],dtype='float') })

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

    #print('\t\t' + 'Running ANOVA.')
    #anova(simsDFanova, FILEID, SAVELOC)

    if '_C_' in FILEID:
        comb = list(combinations(['DOSE','TREAT RATIO','CAR AFFINITY','ANTIGENS CANCER'], 2))
        response_list = ['Y_NORM_CANCER_LIVE', 'Y_NORM_TCELL_LIVE']
    else:
        comb = list(combinations(['DOSE','TREAT RATIO','CAR AFFINITY','ANTIGENS CANCER','ANTIGENS HEALTHY'], 2))
        response_list = ['Y_NORM_CANCER_LIVE', 'Y_NORM_HEALTHY_LIVE', 'Y_NORM_TCELL_LIVE']

    # if ('_CH_' in FILEID and X == 5) or ('_C_' in FILEID and X == 4):
    #     print('\t\tPlotting output two-axis heatmaps.')
    #     for response in response_list:
    #         for combo in comb:
    #             plot_two_axis_heatmaps_all(simsDF, combo[0], combo[1], response, FILEID, SAVELOC)
    #             plot_two_axis_heatmaps_constant(simsDF, combo[0], combo[1], response, FILEID, SAVELOC)
    #             plt.close('all')

    if '_CH_' in FILEID:
        scorelist = []
        for i in range(0, len(simsDFanova)):
            cancer = simsDFanova.iloc[i]['Y_NORM_CANCER_LIVE']
            healthy = simsDFanova.iloc[i]['Y_NORM_HEALTHY_LIVE']

            # calculate cancer score
            delta = 1
            if cancer > 1 or healthy < 0.8:
                delta = 0

            score = (1-cancer + healthy)*delta

            scorelist.append(score)
        simsDFanova['SCORE'] = scorelist
        response_list.append('SCORE')

    print('\t\tPlotting output all-axis heatmaps.')
    for response in response_list:
        plot_output_heatmap(simsDFanova, response, FILEID, SAVELOC)
    # plot_output_heatmap_multiple_outputs(simsDFanova, response_list, FILEID, False, SAVELOC)
    # plot_output_heatmap_multiple_outputs(simsDFanova, response_list, FILEID, True, SAVELOC)
    plt.close('all')

    # colorlist = ['DOSE','TREAT_RATIO','CAR_AFFINITY','ANTIGENS_CANCER']
    # if '_CH_' in FILEID: colorlist.append('ANTIGENS_HEALTHY')
    #
    # print('\t\tPlotting output elbow plots.')
    # for response in response_list:
    #     for color in colorlist:
    #         plot_output_elbow_sort_output_and_color(simsDFanova, response, color, FILEID, SAVELOC)
    #         plot_output_elbow_sort_all(simsDFanova, response, color, FILEID, SAVELOC)
    #         plt.close('all')

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

        stats_data(simsDF, args.color, args.marker, FILEID, args.saveLoc)