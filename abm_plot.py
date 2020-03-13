import ABM
import os
import pickle
import re
import matplotlib.pyplot as plt
import seaborn as sns
from colour import Color
from matplotlib import colors as mcolors
from matplotlib.colors import to_rgba
from argparse import ArgumentParser

'''
ABM_PLOT takes a directory of (or a single) .pkl simulation files that result from ABM_ANALYZE and
plots the DATA features in the dataframe over time.

Usage:
    python abm_plot.py FILES [--color COLOR] [--saveLoc SAVELOC]

    FILES
        Path to .pkl or directory
    [--color COLOR]
        Feature by which to color data by (default: CAR AFFINITY)
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

STRING_FEATURES = ['SEED', 'PLATE', 'TREAT RATIO']
POP_NAMES = ['CANCER', 'CANCER LIVE', 'HEALTHY', 'HEALTHY LIVE', 'T-CELL', 'T-CELL LIVE', 'CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']
POP_FRAC_NAMES = ['CANCER LIVE', 'HEALTHY LIVE']
STATES_CANCER = ['NEUTR', 'APOPT', 'QUIES', 'MIGRA', 'PROLI', 'SENES', 'NECRO']
STATES_TISSUELIVE = ['NEUTR', 'QUIES', 'MIGRA', 'PROLI', 'SENES']
STATES_TCELL = ['NEUTR', 'APOPT', 'MIGRA', 'PROLI', 'SENES', 'CYTOT', 'STIMU', 'EXHAU', 'ANERG', 'STARV', 'PAUSE']

# Set AFFINITY color scale
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
color1 = Color(colors['peachpuff'])
ACOLORS = list(color1.range_to((Color(colors['saddlebrown'])), 5))
for i in range(0, len(ACOLORS)):
    ACOLORS[i] = Color(ACOLORS[i]).rgb

affinityColorDict = {
    "0.0": 'black',
    "1e-06": ACOLORS[0],
    "1e-07": ACOLORS[1],
    "1e-08": ACOLORS[2],
    "1e-09": ACOLORS[3],
    "1e-10": ACOLORS[4]
}

doseColorDict = {
    "0": 'black',
    "250": 'lightgreen',
    "500": "green",
    "1000": "darkgreen"
}

trColorDict = {
    "NA": 'black',
    "0:100": 'lightcoral',
    "50:50": "indianred",
    "100:0": "maroon",
}

caColorDict = {
    "0": 'black',
    "100": 'skyblue',
    "500": "dodgerblue",
    "1000": "cornflowerblue",
    "5000": "royalblue",
    "10000": "navy"
}

COLOR_DICT = {
    "CAR AFFINITY": affinityColorDict,
    "DOSE": doseColorDict,
    "TREAT RATIO": trColorDict,
    "ANTIGENS CANCER": caColorDict
}

doseLineDict = {
    "0": 'solid',
    "250": "dotted",
    "500": "solid",
    "1000": "dashed"
}

# ---------------- LITERATURE DATA ------------------------

# NOTE: This data (Arcangeli 2017) came from co-cultures. Didn't seem to do single culture experiments.
Arcangeli2017 = {
    "CARS": {
        "WT": {
            "Data": {
                "Antigens": [7500, 1600, 1600, 0],
                "Ant Err": [2000, 200, 300, 0],
                "Kill %": [0.725, 0.625, 0.40, 0.15],
                "Kill % Err": [0.075, 0.025, 0.10, 0.05]
            },
            "KD": 1e-9
        },
        "CAMH1": {
            "Data": {
                "Antigens": [7500, 1600, 1600, 0],
                "Ant Err": [2000, 200, 300, 0],
                "Kill %": [0.725, 0.675, 0.35, 0.15],
                "Kill % Err": [0.075, 0.025, 0.05, 0.05]
            },
            "KD": 1e-9
        },
        "CAML": {
            "Data": {
                "Antigens": [7500, 1600, 1600, 0],
                "Ant Err": [2000, 200, 300, 0],
                "Kill %": [0.70, 0.60, 0.35, 0.15],
                "Kill % Err": [0.10, 0.05, 0.05, 0.05]
            },
            "KD": 9e-7
        },
    },
    "CITATION": "Arcangeli 2017 (CD123-CD28-OX40-CD3\u03B6, in vitro, E:T Ratio: 5:1)",
    "MARKER": "o"

}
Wantabe2014 = {
    "CARS": {
        "CD20": {
            "Data": {
                "Antigens": [240, 5320, 26990, 142722],
                "Ant Err": [0, 0, 0, 0],
                "Kill %": [0.18, 0.3, 0.32, 0.37],
                "Kill % Err": [0.04, 0.1, 0.1, 0.06]
            },
            "KD": 1e-10
        }
    },
    "CITATION": "Wantabe 2014 (CD20-CD28-CD3\u03B6, in vitro, E:T Ratio: 10:1,3:1,1:1)",
    "MARKER": "^"
}
Ghorashian2019 = {
    "CARS": {
        "CAT": {
            "Data": {
                "Antigens": [27, 2285, 22307],
                "Ant Err": [0, 0, 0],
                "Kill %": [0.05, 0.21, 0.45],
                "Kill % Err": [0.01, 0.05, 0.03]
            },
            "KD": 1.4e-8
        },
        "FMC": {
            "Data": {
                "Antigens": [27, 2285, 22307],
                "Ant Err": [0, 0, 0],
                "Kill %": [0.05, 0.2857, 0.5],
                "Kill % Err": [0.01, 0.066, 0.03]
            },
            "KD": 3.23e-10
        }
    },
    "CITATION": "Ghorashian 2019 (CD19-41BB-CD3\u03B6, in vitro, E:T Ratio: 6.4:1,1:1)",
    "MARKER": "s"
}
Caruso2015 = {
    "CARS": {
        "Cetux-CAR": {
            "Data": {
                "Antigens": [1e6],
                "Ant Err": [0],
                "Kill %": [0.6],
                "Kill % Err": [0.1]
            },
            "KD": 8.5e-9
        },
        "Nimo-CAR": {
            "Data": {
                "Antigens": [1e6],
                "Ant Err": [0],
                "Kill %": [0.5],
                "Kill % Err": [0.1]
            },
            "KD": 4.53e-8
        }
    },
    "CITATION": "Caruso (EGFR-CD28-CD3\u03B6, in vitro, E:T Ratio: 5:1)",
    "MARKER": "v"
}

litDict = [Arcangeli2017, Wantabe2014, Ghorashian2019, Caruso2015]

# ---------------- PLOTTING FUNCTIONS ---------------------

def plot_counts(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure()
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        ax.plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i][POP_NAME],
                linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=14, labelpad=5)
    ax.set_ylabel(POP_NAME + " CELL COUNTS (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
    ax.set_title(POP_NAME + " COUNT OVER TIME", fontname='Arial', fontweight='bold', fontsize=12, pad=5)
    ax.set_xlim([min(simsDF.iloc[0]['TIME']), max(simsDF.iloc[0]['TIME'])])
    ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_count_fracs(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCountsFrac = plt.figure()
    ax = figCountsFrac.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        ax.plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i][POP_NAME],
                linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=14, labelpad=5)
    ax.set_ylabel(POP_NAME + " CELL FRACTION", fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
    ax.set_title(POP_NAME + " OVER TIME", fontname='Arial', fontweight='bold', fontsize=12, pad=5)
    ax.set_xlim([min(simsDF.iloc[0]['TIME']), max(simsDF.iloc[0]['TIME'])])
    ax.set_ylim([0, 1])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSFRAC_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

# STUB
def plot_states(simsDF, COLOR, FILEID, SAVELOC):

    return

def plot_state_fracs(simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]
    figStateFracs, axs = plt.subplots(12, 7, figsize=(10, 15))

    for i in range(0, len(simsDF)):
        # Plot CANCER cell state fracs
        axs[0, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['QUIES ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['APOPT ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['NECRO ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[6, 0].set_visible(False)
        axs[7, 0].set_visible(False)
        axs[8, 0].set_visible(False)
        axs[9, 0].set_visible(False)
        axs[10, 0].set_visible(False)
        axs[11, 0].set_visible(False)

        # Plot CANCER LIVE cell state fracs
        axs[0, 1].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[1] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 1].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[1] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 1].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['QUIES ' + POP_NAMES[1] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 1].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[1] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 1].set_visible(False)
        axs[5, 1].set_visible(False)
        axs[6, 1].set_visible(False)
        axs[7, 1].set_visible(False)
        axs[8, 1].set_visible(False)
        axs[9, 1].set_visible(False)
        axs[10, 1].set_visible(False)
        axs[11, 1].set_visible(False)

        # Plot HEALTHY cell state fracs
        axs[0, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[2] + ' %'],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[2] + ' %'],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['QUIES ' + POP_NAMES[2] + ' %'],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[2] + ' %'],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['APOPT ' + POP_NAMES[2] + ' %'],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['NECRO ' + POP_NAMES[2] + ' %'],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[6, 2].set_visible(False)
        axs[7, 2].set_visible(False)
        axs[8, 2].set_visible(False)
        axs[9, 2].set_visible(False)
        axs[10, 2].set_visible(False)
        axs[11, 2].set_visible(False)

        # Plot HEALTHY LIVE cell state fracs
        axs[0, 3].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[3] + ' %'],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 3].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[3] + ' %'],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 3].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['QUIES ' + POP_NAMES[3] + ' %'],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 3].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[3] + ' %'],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 3].set_visible(False)
        axs[5, 3].set_visible(False)
        axs[6, 3].set_visible(False)
        axs[7, 3].set_visible(False)
        axs[8, 3].set_visible(False)
        axs[9, 3].set_visible(False)
        axs[10, 3].set_visible(False)
        axs[11, 3].set_visible(False)

        # Plot T-CELL state fracs
        axs[0, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 4].set_visible(False)
        axs[3, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['APOPT ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 4].set_visible(False)
        axs[6, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['CYTOT ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[7, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['STIMU ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[8, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['EXHAU ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[9, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['ANERG ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[10, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['STARV ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[11, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PAUSE ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])

        # Plot CD4 T-CELL state fracs
        axs[0, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 5].set_visible(False)
        axs[3, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['APOPT ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 5].set_visible(False)
        axs[6, 5].set_visible(False)
        axs[7, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['STIMU ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[8, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['EXHAU ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[9, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['ANERG ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[10, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['STARV ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[11, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PAUSE ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])

        # Plot CD8 T-CELL state fracs
        axs[0, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 6].set_visible(False)
        axs[3, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['APOPT ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 6].set_visible(False)
        axs[6, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['CYTOT ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[7, 6].set_visible(False)
        axs[8, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['EXHAU ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[9, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['ANERG ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[10, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['STARV ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[11, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PAUSE ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])

    for a in range(0, 10):
        for b in range(0, 7):
            axs[a, b].set_xlim([0, max(simsDF.iloc[0]['TIME'])])
            axs[a, b].set_ylim([0, 1])

    YLABELS = ['MIGRAGORY', 'PROLIFERATIVE', 'QUIESCENT', 'SENESCENT', 'APOPTOTIC',
               'NECROTIC', 'CYTOTOXIC', 'STIMULATORY', 'EXHAUSTED', 'ANERGIC', 'STARVED', 'PAUSED']
    TITLES = ['CANCER CELLS', 'ALIVE CANCER CELLS', 'HEALTHY CELLS', 'ALIVE HEALTHY CELLS',
              'CAR T-CELLS', 'CD4 CAR T-CELLS', 'CD8 CAR T-CELLS']

    listYAxes = [0, 7, 14, 21, 28, 30, 32, 35, 37, 46, 48, 53, 60, 67, 74, 81]
    listXAxes = [11, 12, 13, 22, 24, 32, 33, 34, 35, 37, 48, 81, 82, 83]
    listYlabels = [0, 7, 14, 21, 28, 35, 46, 53, 60, 67, 74, 81]
    j = 0
    k = 0
    for i, ax in enumerate(axs.flat):
        if i in listYAxes and i in listYlabels:
            ax.set_ylabel(YLABELS[j], fontname='Arial', fontweight='bold', fontsize=7, labelpad=10)
            ax.get_yaxis().set_visible(True)
            j = j + 1
        elif i in listYAxes and i not in listYlabels:
            ax.get_yaxis().set_visible(True)
        else:
            ax.get_yaxis().set_visible(False)
        if i < 7:
            ax.set_title(TITLES[k], fontname='Arial', fontweight='bold', fontsize=7, pad=10)
            k = k + 1
        if i in listXAxes:
            ax.get_xaxis().set_visible(True)
        else:
            ax.get_xaxis().set_visible(False)

    figStateFracs.text(0.5, 0.08, "TIME", ha='center', fontname='Arial', fontweight='bold', fontsize=13)
    figStateFracs.text(0.01, 0.5, "STATE (FRACTION)", va='center', rotation='vertical', fontname='Arial',
                     fontweight='bold', fontsize=13)

    for color in colorDict:
        axs[0,0].plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        axs[0,0].plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    axs[0, 0].legend(bbox_to_anchor=(9.5, 1), frameon=False)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_STATESFRAC.svg', bbox_inches='tight')

    return

def plot_state_fracs_neutral(simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]
    figStateFracs, axs = plt.subplots(13, 7, figsize=(10, 15))

    for i in range(0, len(simsDF)):
        # Plot CANCER cell state fracs
        axs[0, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['NEUTR ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['QUIES ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['APOPT ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[6, 0].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['NECRO ' + POP_NAMES[0] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[7, 0].set_visible(False)
        axs[8, 0].set_visible(False)
        axs[9, 0].set_visible(False)
        axs[10, 0].set_visible(False)
        axs[11, 0].set_visible(False)
        axs[12, 0].set_visible(False)

        # Plot CANCER LIVE cell state fracs
        axs[0, 1].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['NEUTR ' + POP_NAMES[1] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 1].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[1] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 1].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[1] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 1].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['QUIES ' + POP_NAMES[1] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 1].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[1] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 1].set_visible(False)
        axs[6, 1].set_visible(False)
        axs[7, 1].set_visible(False)
        axs[8, 1].set_visible(False)
        axs[9, 1].set_visible(False)
        axs[10, 1].set_visible(False)
        axs[11, 1].set_visible(False)
        axs[12, 1].set_visible(False)

        # Plot HEALTHY cell state fracs
        axs[0, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['NEUTR ' + POP_NAMES[2] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[2] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[2] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['QUIES ' + POP_NAMES[2] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[2] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['APOPT ' + POP_NAMES[2] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[6, 2].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['NECRO ' + POP_NAMES[2] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[7, 2].set_visible(False)
        axs[8, 2].set_visible(False)
        axs[9, 2].set_visible(False)
        axs[10, 2].set_visible(False)
        axs[11, 2].set_visible(False)
        axs[12, 2].set_visible(False)

        # Plot HEALTHY LIVE cell state fracs
        axs[0, 3].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['NEUTR ' + POP_NAMES[3] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 3].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[3] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 3].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[3] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 3].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['QUIES ' + POP_NAMES[3] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 3].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[3] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 3].set_visible(False)
        axs[6, 3].set_visible(False)
        axs[7, 3].set_visible(False)
        axs[8, 3].set_visible(False)
        axs[9, 3].set_visible(False)
        axs[10, 3].set_visible(False)
        axs[11, 3].set_visible(False)
        axs[12, 3].set_visible(False)

        # Plot T-CELL state fracs
        axs[0, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['NEUTR ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 4].set_visible(False)
        axs[4, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['APOPT ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[6, 4].set_visible(False)
        axs[7, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['CYTOT ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[8, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['STIMU ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[8, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['EXHAU ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[10, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['ANERG ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[11, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['STARV ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[12, 4].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PAUSE ' + POP_NAMES[4] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])

        # Plot CD4 T-CELL state fracs
        axs[0, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['NEUTR ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 5].set_visible(False)
        axs[4, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['APOPT ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[6, 5].set_visible(False)
        axs[7, 5].set_visible(False)
        axs[8, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['STIMU ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[8, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['EXHAU ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[10, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['ANERG ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[11, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['STARV ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[12, 5].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PAUSE ' + POP_NAMES[6] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])

        # Plot CD8 T-CELL state fracs
        axs[0, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['NEUTR ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['MIGRA ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PROLI ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 6].set_visible(False)
        axs[4, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['SENES ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['APOPT ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[6, 6].set_visible(False)
        axs[7, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['CYTOT ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[8, 6].set_visible(False)
        axs[9, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['EXHAU ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[10, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['ANERG ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[11, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['STARV ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[12, 6].plot(simsDF.iloc[i]['TIME'], simsDF.iloc[i]['PAUSE ' + POP_NAMES[8] + ' %'], color=colorDict[str(simsDF.iloc[i][COLOR])])

    for a in range(0, 11):
        for b in range(0, 7):
            axs[a, b].set_xlim([0, max(simsDF.iloc[0]['TIME'])])
            axs[a, b].set_ylim([0, 1])

    YLABELS = ['NEUTRAL', 'MIGRAGORY', 'PROLIFERATIVE', 'QUIESCENT', 'SENESCENT', 'APOPTOTIC',
               'NECROTIC', 'CYTOTOXIC', 'STIMULATORY', 'EXHAUSTED', 'ANERGIC', 'STARVED', 'PAUSED']
    TITLES = ['CANCER CELLS', 'ALIVE CANCER CELLS', 'HEALTHY CELLS', 'ALIVE HEALTHY CELLS',
              'CAR T-CELLS', 'CD4 CAR T-CELLS', 'CD8 CAR T-CELLS']

    listYAxes = [0, 7, 14, 21, 28, 35, 37, 39, 42, 44, 53, 55, 60, 67, 74, 81, 88]
    listXAxes = [18, 19, 20, 29, 31, 39, 40, 41, 42, 44, 55, 88, 89, 90]
    listYlabels = [0, 7, 14, 21, 28, 35, 42, 53, 60, 67, 74, 81, 88]
    j = 0
    k = 0
    for i, ax in enumerate(axs.flat):
        if i in listYAxes and i in listYlabels:
            ax.set_ylabel(YLABELS[j], fontname='Arial', fontweight='bold', fontsize=7, labelpad=10)
            ax.get_yaxis().set_visible(True)
            j = j + 1
        elif i in listYAxes and i not in listYlabels:
            ax.get_yaxis().set_visible(True)
        else:
            ax.get_yaxis().set_visible(False)
        if i < 7:
            ax.set_title(TITLES[k], fontname='Arial', fontweight='bold', fontsize=7, pad=10)
            k = k + 1
        if i in listXAxes:
            ax.get_xaxis().set_visible(True)
        else:
            ax.get_xaxis().set_visible(False)

    figStateFracs.text(0.5, 0.08, "TIME", ha='center', fontname='Arial', fontweight='bold', fontsize=13)
    figStateFracs.text(0.01, 0.5, "STATE (FRACTION)", va='center', rotation='vertical', fontname='Arial',
                     fontweight='bold', fontsize=13)

    for color in colorDict:
        axs[0,0].plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        axs[0,0].plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    axs[0, 0].legend(bbox_to_anchor=(9.5, 1), frameon=False)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_STATESFRACNEUTR.svg', bbox_inches='tight')

    return

def plot_kill_curve_sim(simsDF, FILEID, SAVELOC, TIME):

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0,len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figKC = plt.figure()
    ax = figKC.add_subplot(1, 1, 1)

    KD6D250 = []
    KD7D250 = []
    KD8D250 = []
    KD9D250 = []

    KD6D500 = []
    KD7D500 = []
    KD8D500 = []
    KD9D500 = []

    KD6D1000 = []
    KD7D1000 = []
    KD8D1000 = []
    KD9D1000 = []


    # Make ANTIGEN and % KILLING CANCER CELLS array sorted by dose
    for i in range(0, len(simsDF)):
        antigen = simsDF.iloc[i]['ANTIGENS CANCER']
        dose = simsDF.iloc[i]['DOSE']
        affinity = simsDF.iloc[i]['CAR AFFINITY']
        killed = 1-simsDF.iloc[i]['CANCER LIVE %'][index]

        if dose == 250:
            if affinity >= 1e-6:
                KD6D250.append([antigen, killed])
            elif 1e-6 > affinity >= 1e-7:
                KD7D250.append([antigen, killed])
            elif 1e-7 > affinity >= 1e-8:
                KD8D250.append([antigen, killed])
            elif 1e-8 > affinity >= 1e-9:
                KD9D250.append([antigen, killed])
            elif 1e-9 > affinity >= 1e-10:
                continue
            else:
                continue

        elif dose == 500:
            if affinity >= 1e-6:
                KD6D500.append([antigen, killed])
            elif 1e-6 > affinity >= 1e-7:
                KD7D500.append([antigen, killed])
            elif 1e-7 > affinity >= 1e-8:
                KD8D500.append([antigen, killed])
            elif 1e-8 > affinity >= 1e-9:
                KD9D500.append([antigen, killed])
            elif 1e-9 > affinity >= 1e-10:
                continue
            else:
                continue

        elif dose == 1000:
            if affinity >= 1e-6:
                KD6D1000.append([antigen, killed])
            elif 1e-6 > affinity >= 1e-7:
                KD7D1000.append([antigen, killed])
            elif 1e-7 > affinity >= 1e-8:
                KD8D1000.append([antigen, killed])
            elif 1e-8 > affinity >= 1e-9:
                KD9D1000.append([antigen, killed])
            elif 1e-9 > affinity >= 1e-10:
                continue
            else:
                continue

    dataLists = [KD6D250, KD7D250, KD8D250, KD9D250,
                 KD6D500, KD7D500, KD8D500, KD9D500,
                 KD6D1000, KD7D1000, KD8D1000, KD9D1000]

    KD6D250AVG = []
    KD7D250AVG = []
    KD8D250AVG = []
    KD9D250AVG = []

    KD6D500AVG = []
    KD7D500AVG = []
    KD8D500AVG = []
    KD9D500AVG = []

    KD6D1000AVG = []
    KD7D1000AVG = []
    KD8D1000AVG = []
    KD9D1000AVG = []

    dataListsSorted = []

    avgLists = [KD6D250AVG, KD7D250AVG, KD8D250AVG, KD9D250AVG,
                KD6D500AVG, KD7D500AVG, KD8D500AVG, KD9D500AVG,
                KD6D1000AVG, KD7D1000AVG, KD8D1000AVG, KD9D1000AVG]

    for dList in dataLists:
        dList.sort(key=lambda x: x[0])
        dataListsSorted.append(dList)

    x = [0,5,10,15,20,25]
    for a in range(0,len(avgLists)):
        if dataListsSorted[a] != []:
            for i in range(0, 5):
                antSum = 0
                killSum = 0
                for j in range(x[i], x[i+1]):
                    antSum += dataListsSorted[a][j][0]
                    killSum += dataListsSorted[a][j][1]
                antAvg = antSum/len(dataListsSorted[a][x[i]:x[i+1]])
                killAvg = killSum / len(dataListsSorted[a][x[i]:x[i + 1]])
                avgLists[a].append([antAvg, killAvg])

    ax.plot([a[0] for a in KD6D250AVG], [a[1] for a in KD6D250AVG], linestyle=doseLineDict['250'],
            color=COLOR_DICT["CAR AFFINITY"]["1e-06"])
    ax.plot([a[0] for a in KD7D250AVG], [a[1] for a in KD7D250AVG], linestyle=doseLineDict['250'],
            color=COLOR_DICT["CAR AFFINITY"]["1e-07"])
    ax.plot([a[0] for a in KD8D250AVG], [a[1] for a in KD8D250AVG], linestyle=doseLineDict['250'],
            color=COLOR_DICT["CAR AFFINITY"]["1e-08"])
    ax.plot([a[0] for a in KD9D250AVG], [a[1] for a in KD9D250AVG], linestyle=doseLineDict['250'],
            color=COLOR_DICT["CAR AFFINITY"]["1e-09"])

    ax.plot([a[0] for a in KD6D500AVG], [a[1] for a in KD6D500AVG], linestyle=doseLineDict['500'],
            color=COLOR_DICT["CAR AFFINITY"]["1e-06"])
    ax.plot([a[0] for a in KD7D500AVG], [a[1] for a in KD7D500AVG], linestyle=doseLineDict['500'],
            color=COLOR_DICT["CAR AFFINITY"]["1e-07"])
    ax.plot([a[0] for a in KD8D500AVG], [a[1] for a in KD8D500AVG], linestyle=doseLineDict['500'],
            color=COLOR_DICT["CAR AFFINITY"]["1e-08"])
    ax.plot([a[0] for a in KD9D500AVG], [a[1] for a in KD9D500AVG], linestyle=doseLineDict['500'],
            color=COLOR_DICT["CAR AFFINITY"]["1e-09"])

    ax.plot([a[0] for a in KD6D1000AVG], [a[1] for a in KD6D1000AVG], linestyle=doseLineDict['1000'],
            color=COLOR_DICT["CAR AFFINITY"]["1e-06"])
    ax.plot([a[0] for a in KD7D1000AVG], [a[1] for a in KD7D1000AVG], linestyle=doseLineDict['1000'],
            color=COLOR_DICT["CAR AFFINITY"]["1e-07"])
    ax.plot([a[0] for a in KD8D1000AVG], [a[1] for a in KD8D1000AVG], linestyle=doseLineDict['1000'],
            color=COLOR_DICT["CAR AFFINITY"]["1e-08"])
    ax.plot([a[0] for a in KD9D1000AVG], [a[1] for a in KD9D1000AVG], linestyle=doseLineDict['1000'],
            color=COLOR_DICT["CAR AFFINITY"]["1e-09"])

    for color in affinityColorDict:
        ax.plot([0],[0], color=affinityColorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    ax.set_title("KILL CURVE (SIMULATED DATA)", fontname='Arial',
                 fontweight='bold', fontsize=14, pad=5)
    ax.set_xlabel("ANTIGENS", fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
    ax.set_ylabel("% LYSIS", fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
    ax.set_ylim([0, 1])
    ax.set_xlim([0, 10000])
    ax.legend(bbox_to_anchor=(1, 1), frameon=False)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_KILLCURVESIM_' + str(TIME) + '.svg', bbox_inches='tight')

    return

def plot_kill_curve_exp(SAVELOC):

    figKC = plt.figure()
    ax = figKC.add_subplot(1, 1, 1)

    # Plot literature values
    for d in litDict:
        for CAR in d["CARS"]:
            if d["CARS"][CAR]["KD"] > 1e-6:
                d["CARS"][CAR].update({"Color": COLOR_DICT["CAR AFFINITY"]["1e-06"]})
            elif 1e-6 > d["CARS"][CAR]["KD"] >= 1e-7:
                d["CARS"][CAR].update({"Color": COLOR_DICT["CAR AFFINITY"]["1e-07"]})
            elif 1e-7 > d["CARS"][CAR]["KD"] >= 1e-8:
                d["CARS"][CAR].update({"Color": COLOR_DICT["CAR AFFINITY"]["1e-08"]})
            elif 1e-8 > d["CARS"][CAR]["KD"] >= 1e-9:
                d["CARS"][CAR].update({"Color": COLOR_DICT["CAR AFFINITY"]["1e-09"]})
            elif 1e-9 > d["CARS"][CAR]["KD"] >= 1e-10:
                d["CARS"][CAR].update({"Color": COLOR_DICT["CAR AFFINITY"]["1e-10"]})
            else:
                d["CARS"][CAR].update({"Color": "black"})

            ax.errorbar(d["CARS"][CAR]["Data"]["Antigens"], d["CARS"][CAR]["Data"]["Kill %"],
                        xerr=d["CARS"][CAR]["Data"]["Ant Err"], fmt='None', yerr=d["CARS"][CAR]["Data"]["Kill % Err"],
                        ecolor='lightgray', zorder=0)
            # ax.scatter(d["CARS"][CAR]["Data"]["Antigens"], d["CARS"][CAR]["Data"]["Kill %"],
            #            marker=d["MARKER"], color=d["CARS"][CAR]["Color"], zorder=1)
            ax.plot(d["CARS"][CAR]["Data"]["Antigens"], d["CARS"][CAR]["Data"]["Kill %"],
                        marker=d["MARKER"], color=d["CARS"][CAR]["Color"], zorder=1)

        ax.scatter(None, None, marker=d["MARKER"], label=d["CITATION"], color='black')

    for color in affinityColorDict:
        ax.plot([0],[0], color=affinityColorDict[color], label=color, linestyle='solid')

    ax.set_title("KILL CURVE (EXPERIMENTAL DATA)", fontname='Arial',
                 fontweight='bold', fontsize=14, pad=5)
    ax.set_xlabel("ANTIGENS", fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
    ax.set_ylabel("% LYSIS", fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
    ax.set_ylim([0, 1])
    ax.set_xlim([0, 10000])
    ax.legend(bbox_to_anchor=(1, 0.40), frameon=False)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'KILLCURVEEXP.svg', bbox_inches='tight')

    return

def plot_ET_ratio(simsDF, FILEID, SAVELOC):
    figKC = plt.figure()
    ax = figKC.add_subplot(1, 1, 1)

    CONTROL = []

    KD6A100 = []
    KD7A100 = []
    KD8A100 = []
    KD9A100 = []

    KD6A500 = []
    KD7A500 = []
    KD8A500 = []
    KD9A500 = []

    KD6A1000 = []
    KD7A1000 = []
    KD8A1000 = []
    KD9A1000 = []

    KD6A5000 = []
    KD7A5000 = []
    KD8A5000 = []
    KD9A5000 = []

    KD6A10000 = []
    KD7A10000 = []
    KD8A10000 = []
    KD9A10000 = []

    # Make ANTIGEN and % KILLING CANCER CELLS array sorted by dose
    for i in range(0, len(simsDF)):
        antigen = simsDF.iloc[i]['ANTIGENS CANCER']
        dose = simsDF.iloc[i]['DOSE']
        affinity = simsDF.iloc[i]['CAR AFFINITY']
        killed = 1 - simsDF.iloc[i]['CANCER LIVE %'][-1]
        init_cancer = simsDF.iloc[i]['CANCER'][0]
        ETratio = dose/init_cancer

        if dose == 0:
            CONTROL.append([ETratio,killed])
        else:
            if antigen == 100:
                if affinity >= 1e-6:
                    KD6A100.append([ETratio, killed])
                elif 1e-6 > affinity >= 1e-7:
                    KD7A100.append([ETratio, killed])
                elif 1e-7 > affinity >= 1e-8:
                    KD8A100.append([ETratio, killed])
                elif 1e-8 > affinity >= 1e-9:
                    KD9A100.append([ETratio, killed])
                elif 1e-9 > affinity >= 1e-10:
                    continue
                else:
                    continue

            elif antigen == 500:
                if affinity >= 1e-6:
                    KD6A500.append([ETratio, killed])
                elif 1e-6 > affinity >= 1e-7:
                    KD7A500.append([ETratio, killed])
                elif 1e-7 > affinity >= 1e-8:
                    KD8A500.append([ETratio, killed])
                elif 1e-8 > affinity >= 1e-9:
                    KD9A500.append([ETratio, killed])
                elif 1e-9 > affinity >= 1e-10:
                    continue
                else:
                    continue

            elif antigen == 1000:
                if affinity >= 1e-6:
                    KD6A1000.append([ETratio, killed])
                elif 1e-6 > affinity >= 1e-7:
                    KD7A1000.append([ETratio, killed])
                elif 1e-7 > affinity >= 1e-8:
                    KD8A1000.append([ETratio, killed])
                elif 1e-8 > affinity >= 1e-9:
                    KD9A1000.append([ETratio, killed])
                elif 1e-9 > affinity >= 1e-10:
                    continue
                else:
                    continue

            elif antigen == 5000:
                if affinity >= 1e-6:
                    KD6A5000.append([ETratio, killed])
                elif 1e-6 > affinity >= 1e-7:
                    KD7A5000.append([ETratio, killed])
                elif 1e-7 > affinity >= 1e-8:
                    KD8A5000.append([ETratio, killed])
                elif 1e-8 > affinity >= 1e-9:
                    KD9A5000.append([ETratio, killed])
                elif 1e-9 > affinity >= 1e-10:
                    continue
                else:
                    continue

            elif antigen == 10000:
                if affinity >= 1e-6:
                    KD6A10000.append([ETratio, killed])
                elif 1e-6 > affinity >= 1e-7:
                    KD7A10000.append([ETratio, killed])
                elif 1e-7 > affinity >= 1e-8:
                    KD8A10000.append([ETratio, killed])
                elif 1e-8 > affinity >= 1e-9:
                    KD9A10000.append([ETratio, killed])
                elif 1e-9 > affinity >= 1e-10:
                    continue
                else:
                    continue

    dataLists = [KD6A100, KD7A100, KD8A100, KD9A100,
                 KD6A500, KD7A500, KD8A500, KD9A500,
                 KD6A1000, KD7A1000, KD8A1000, KD9A1000,
                 KD6A5000, KD7A5000, KD8A5000, KD9A5000,
                 KD6A10000, KD7A10000, KD8A10000, KD9A10000]

    KD6A100AVG = []
    KD7A100AVG = []
    KD8A100AVG = []
    KD9A100AVG = []

    KD6A500AVG = []
    KD7A500AVG = []
    KD8A500AVG = []
    KD9A500AVG = []

    KD6A1000AVG = []
    KD7A1000AVG = []
    KD8A1000AVG = []
    KD9A1000AVG = []

    KD6A5000AVG = []
    KD7A5000AVG = []
    KD8A5000AVG = []
    KD9A5000AVG = []

    KD6A10000AVG = []
    KD7A10000AVG = []
    KD8A10000AVG = []
    KD9A10000AVG = []

    dataListsSorted = []

    avgLists = [KD6A100AVG, KD7A100AVG, KD8A100AVG, KD9A100AVG,
                KD6A500AVG, KD7A500AVG, KD8A500AVG, KD9A500AVG,
                KD6A1000AVG, KD7A1000AVG, KD8A1000AVG, KD9A1000AVG,
                KD6A5000AVG, KD7A5000AVG, KD8A5000AVG, KD9A5000AVG,
                KD6A10000AVG, KD7A10000AVG, KD8A10000AVG, KD9A10000AVG]

    for dList in dataLists:
        if dList != []:
            dList = CONTROL + dList
        dList.sort(key=lambda x: x[0])
        dataListsSorted.append(dList)

    x = [0, 5, 10, 15, 20]
    for a in range(0, len(avgLists)):
        if dataListsSorted[a] != []:
            for i in range(0, 4):
                ETSum = 0
                killSum = 0
                for j in range(x[i], x[i + 1]):
                    ETSum += dataListsSorted[a][j][0]
                    killSum += dataListsSorted[a][j][1]
                ETAvg = ETSum / len(dataListsSorted[a][x[i]:x[i + 1]])
                killAvg = killSum / len(dataListsSorted[a][x[i]:x[i + 1]])
                avgLists[a].append([ETAvg, killAvg])

    ax.plot([a[0] for a in KD6A100AVG], [a[1] for a in KD6A100AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-06"])
    ax.plot([a[0] for a in KD7A100AVG], [a[1] for a in KD7A100AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-07"])
    ax.plot([a[0] for a in KD8A100AVG], [a[1] for a in KD8A100AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-08"])
    ax.plot([a[0] for a in KD9A100AVG], [a[1] for a in KD9A100AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-09"])

    ax.plot([a[0] for a in KD6A500AVG], [a[1] for a in KD6A500AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-06"])
    ax.plot([a[0] for a in KD7A500AVG], [a[1] for a in KD7A500AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-07"])
    ax.plot([a[0] for a in KD8A500AVG], [a[1] for a in KD8A500AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-08"])
    ax.plot([a[0] for a in KD9A500AVG], [a[1] for a in KD9A500AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-09"])

    ax.plot([a[0] for a in KD6A1000AVG], [a[1] for a in KD6A1000AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-06"])
    ax.plot([a[0] for a in KD7A1000AVG], [a[1] for a in KD7A1000AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-07"])
    ax.plot([a[0] for a in KD8A1000AVG], [a[1] for a in KD8A1000AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-08"])
    ax.plot([a[0] for a in KD9A1000AVG], [a[1] for a in KD9A1000AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-09"])

    ax.plot([a[0] for a in KD6A5000AVG], [a[1] for a in KD6A5000AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-06"])
    ax.plot([a[0] for a in KD7A5000AVG], [a[1] for a in KD7A5000AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-07"])
    ax.plot([a[0] for a in KD8A5000AVG], [a[1] for a in KD8A5000AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-08"])
    ax.plot([a[0] for a in KD9A5000AVG], [a[1] for a in KD9A5000AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-09"])

    ax.plot([a[0] for a in KD6A10000AVG], [a[1] for a in KD6A10000AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-06"])
    ax.plot([a[0] for a in KD7A10000AVG], [a[1] for a in KD7A10000AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-07"])
    ax.plot([a[0] for a in KD8A10000AVG], [a[1] for a in KD8A10000AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-08"])
    ax.plot([a[0] for a in KD9A10000AVG], [a[1] for a in KD9A10000AVG],
            color=COLOR_DICT["CAR AFFINITY"]["1e-09"])

    for color in affinityColorDict:
        ax.plot([0],[0], color=affinityColorDict[color], label=color, linestyle='solid')
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    ax.set_title("Killing Ability as a Function of E:T Ratio", fontname='Arial',
                 fontweight='bold', fontsize=14, pad=5)
    ax.set_xlabel("E:T Ratio", fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
    ax.set_ylabel("% LYSIS", fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
    ax.set_ylim([0, 1])
    ax.set_xlim(left=0)
    ax.legend(bbox_to_anchor=(1, 1), frameon=False)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_ETRATIO.svg', bbox_inches='tight')

    return

def plot_volumes(simsDF, COLOR, FILEID, SAVELOC, TIME):

    colorDict = COLOR_DICT[COLOR]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    for p in range(0,len(POP_NAMES)):
        if 'LIVE' not in POP_NAMES[p]:
            key = 'CELL VOLUMES ' + POP_NAMES[p]
            Y = []
            X = []
            order = []

            for i in range(0, len(simsDF)):
                y = simsDF.iloc[i][key][index]
                x = [str(simsDF.iloc[i][COLOR])]*len(y)
                Y = Y + y
                X = X + x

            for key in colorDict:
                if key in X:
                    order.append(key)

            if X != [] and Y != []:
                figV = plt.figure()
                ax = figV.add_subplot(1, 1, 1)
                ax = sns.violinplot(x=X, y=Y, order=order, palette=colorDict)
                ax.set_xlabel(POP_NAMES[p] + " POPULATION", fontname='Arial', fontweight='bold', fontsize=14,
                              labelpad=5)
                ax.set_ylabel("VOLUME (um^3)", fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
                ax.set_title(POP_NAMES[p] + " VOLUME DISTRIBUTIONS AT TIME " + str(TIME), fontname='Arial',
                             fontweight='bold', fontsize=12, pad=5)
                ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)

                if SAVELOC == '':
                    plt.show()
                else:
                    plt.savefig(SAVELOC + FILEID + '_VOLUMES_' + POP_NAMES[p].replace(' ', '') + '_DAY_' + str(TIME) + '.svg', bbox_inches='tight')

    return


def plot_cycles(simsDF, COLOR, FILEID, SAVELOC, TIME):
    colorDict = COLOR_DICT[COLOR]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    for p in range(0, len(POP_NAMES)):
        if 'LIVE' not in POP_NAMES[p]:
            key = 'AVG CELL CYCLES ' + POP_NAMES[p]
            Y = []
            X = []
            order = []

            for i in range(0, len(simsDF)):
                y = simsDF.iloc[i][key][index]
                x = [str(simsDF.iloc[i][COLOR])] * len(y)
                Y = Y + y
                X = X + x

            for key in colorDict:
                if key in X:
                    order.append(key)

            if X != [] and Y != []:
                figV = plt.figure()
                ax = figV.add_subplot(1, 1, 1)
                ax = sns.violinplot(x=X, y=Y, order=order, palette=colorDict)
                ax.set_xlabel(POP_NAMES[p] + " POPULATION", fontname='Arial', fontweight='bold', fontsize=14, labelpad=5)
                ax.set_ylabel("CYCLE LENGTH (min)", fontname='Arial', fontweight='bold', fontsize=12, labelpad=5)
                ax.set_title(POP_NAMES[p] + " CELL CYCLE DISTRIBUTIONS AT TIME " + str(TIME), fontname='Arial',
                             fontweight='bold', fontsize=12, pad=5)
                ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)

                if SAVELOC == '':
                    plt.show()
                else:
                    plt.savefig(SAVELOC + FILEID + '_CYCLES_' + POP_NAMES[p].replace(' ', '') + '_DAY_' + str(TIME) + '.svg', bbox_inches='tight')

    return

def plot_data(simsDF, COLOR, FILEID, SAVELOC):

    # # Plot counts data
    # for p in range(0, len(POP_NAMES)):
    #     plot_counts(POP_NAMES[p], simsDF, COLOR, FILEID, SAVELOC)
    #
    #
    # # Plot count fracs data
    # for p in range(0,len(POP_FRAC_NAMES)):
    #     POP_NAME = POP_FRAC_NAMES[p] + ' %'
    #     plot_count_fracs(POP_NAME, simsDF, COLOR, FILEID, SAVELOC)

    # Plot types


    # Plot type fracs
    # plot_state_fracs(simsDF, COLOR, FILEID, SAVELOC)
    # plot_state_fracs_neutral(simsDF, COLOR, FILEID, SAVELOC)

    # Plot kill curve
    # plot_kill_curve_sim(simsDF, FILEID, SAVELOC, 4)
    # plot_kill_curve_sim(simsDF, FILEID, SAVELOC, 7)
    # plot_kill_curve_exp(SAVELOC)

    # Plot ET Ratio
    # plot_ET_ratio(simsDF, FILEID, SAVELOC)

    # Plot volume distributions
    plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 4)
    plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 7)
    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 4)
    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 7)

    return

# ---------------- MAIN FUNCTION ---------------------

if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()

    # Get files
    PKLFILES = get_pkl_files(args.files)

    for file in PKLFILES:

        if 'VITRO' in file:
            fileName = re.sub('.*VITRO', 'VITRO', file)

        else:
            fileName = re.sub('.*VIVO', 'VIVO', file)

        FILEID = fileName.replace('_ANALYZED.pkl', '')

        simsDF = pickle.load(open(file, 'rb'))
        plot_data(simsDF, args.color, FILEID, args.saveLoc)