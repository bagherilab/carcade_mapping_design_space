import ABM
import os
import pickle
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from colour import Color
from matplotlib import cm as mplcm
from matplotlib import colors as mcolors
from matplotlib.colors import to_rgba
from argparse import ArgumentParser
from math import exp
import numpy as np

'''
ABM_PLOT takes a directory of (or a single) .pkl simulation files that result from ABM_ANALYZE and
plots the DATA features in the dataframe over time.

Usage:
    python abm_plot.py FILES [--color COLOR] [--marker MARKER] [--saveLoc SAVELOC]

    FILES
        Path to .pkl or directory
    [--color COLOR]
        Feature by which to color data by. If X given, will color by axis along which data varies.
        If two featuers vary, default will be used. (default: CAR AFFINITY)
    [--marker MARKER]
        Feature by which to change marker of data by (default: TREAT RATIO)
    [--partial]
        Flag indicating only partial dataset present instead of full combinatoral set.
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
    parser.add_argument("--partial", default=False, dest="partial",
                        action='store_true', help="Flag indicating only partial dataset present instead of full combinatoral set.")
    parser.add_argument("--rank", default=False, dest="rank", action='store_true',
                        help="Flag indicating to make rank parody plot using rank file.")
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

TICKSIZE = 30
FONTSIZE_AXES_VALUES = 20
FONTSIZE_AXES_TITLES = 25
LABELPAD = 10

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
trvalues = [0.2, 0.3, 0.4, 0.6, 0.8, 0.9, 1.0] # CD4%: 0, 10, 25, 50, 75, 90, 100
#trvalues = np.linspace(0.2, 1.0, num=7)
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

doseColorDict = {
    "0": 'black',
    "250": DCOLORS[0],
    "500": DCOLORS[1],
    "1000": DCOLORS[2]
}

# trColorDict = {
#     "NA": 'black',
#     "0:100": TRCOLORS[0],
#     "25:75": TRCOLORS[1],
#     "50:50": TRCOLORS[2],
#     "75:25": TRCOLORS[3],
#     "100:0": TRCOLORS[4],
# }

trColorDict = {
    "NA": 'black',
    "0:100": TRCOLORS[0],
    "10:90": TRCOLORS[1],
    "25:75": TRCOLORS[2],
    "50:50": TRCOLORS[3],
    "75:25": TRCOLORS[4],
    "90:10": TRCOLORS[5],
    "100:0": TRCOLORS[6],
}

affinityColorDict = {
    "NA": 'black',
    "0.0": 'black',
    "1e-06": ACOLORS[0],
    "1e-07": ACOLORS[1],
    "1e-08": ACOLORS[2],
    "1e-09": ACOLORS[3],
    "1e-10": ACOLORS[4]
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

TREAT_RATIO_DICT = {'0:100': 0.0, '10:90': 0.1, '25:75': 0.25, '50:50': 0.5, '75:25': 0.75, '90:10': 0.9, '100:0': 1.0}
TREAT_RATIO_DICT_REVERSE = {'0.0': '0:100', '0.1': '10:90', '0.25': '25:75', '0.5': '50:50', '0.75': '75:25', '0.9': '90:10', '1.0': '100:0'}
# ---------------- LITERATURE DATA ------------------------

# NOTE: This data (Arcangeli 2017) came from co-cultures with a low antigen expresser.
#       Didn't seem to do single culture experiments.

# NOTE: Data from Fig 4 Short-Term Assay (4 hrs)
Arcangeli2017 = {
    "CARS": {
        "WT": {
            "Data": {
                "Antigens": [7500, 1600, 1600, 0],  # THP-1, U937, TIME, MHH-CALL-4
                "Ant Err": [2000, 200, 300, 0],
                "Kill %": [0.78, 0.68, 0.33, 0.10],
                "Kill % Err": [0.02, 0.05, 0.10, 0.02]
            },
            "Normalized Data": {
                "Antigens": [1, 0.213, 0.213, 0],
                "Ant Err": [0.377, 0.0628, 0.0695, 0],
                "Kill %": [0.918, 0.800, 0.388, 0.118],
                "Kill % Err": [0.0319, 0.0618, 0.118, 0.0237]
            },
            "KD": 1e-9
        },
        "CAMH1": {
            "Data": {
                "Antigens": [7500, 1600, 1600, 0],
                "Ant Err": [2000, 200, 300, 0],
                "Kill %": [0.85, 0.72, 0.37, 0.15],
                "Kill % Err": [0.02, 0.05, 0.05, 0.02]
            },
            "Normalized Data": {
                "Antigens": [1, 0.213, 0.213, 0],
                "Ant Err": [0.377, 0.0628, 0.0695, 0],
                "Kill %": [1, 0.847, 0.435, 0.176],
                "Kill % Err": [0.0332, 0.0621, 0.0597, 0.0240]
            },
            "KD": 1e-9
        },
        "CAML": {
            "Data": {
                "Antigens": [7500, 1600, 1600, 0],
                "Ant Err": [2000, 200, 300, 0],
                "Kill %": [0.8, 0.57, 0.33, 0.15],
                "Kill % Err": [0.02, 0.05, 0.05, 0.02]
            },
            "Normalized Data": {
                "Antigens": [1, 0.213, 0.213, 0],
                "Ant Err": [0.377, 0.0628, 0.0695, 0],
                "Kill %": [0.941, 0.671, 0.389, 0.176],
                "Kill % Err": [0.0323, 0.0609, 0.0595, 0.0239]
            },
            "KD": 9e-7
        },
    },
    "CITATION": "Arcangeli 2017 (CD123-CD28-OX40-CD3\u03B6, in vitro, co-culture, E:T Ratio: 5:1)",
    "SAVE": "Arcangeli2017",
    "MARKER": "o"

}

# NOTE: Data from Figure 2C, used 1:1 data
Wantanabe2015 = {
    "CARS": {
        "CD20": {
            "Data": {
                "Antigens": [142722, 26990, 5320, 240],
                "Ant Err": [0, 0, 0, 0],
                "Kill %": [0.37, 0.32, 0.3, 0.18],
                "Kill % Err": [0.06, 0.1, 0.1, 0.04]
            },
            "Normalized Data": {
                "Antigens": [1, 0.189, 0.0373, 0.00168],
                "Ant Err": [0, 0, 0, 0],
                "Kill %": [1, 0.865, 0.811, 0.486],
                "Kill % Err": [0.229, 0.304, 0.301, 0.132]
            },
            "KD": 1e-10
        }
    },
    "CITATION": "Wantabe 2014 (CD20-CD28-CD3\u03B6, in vitro, E:T Ratio: 10:1,3:1,1:1)",
    "SAVE": "Watanabe2015",
    "MARKER": "^"
}

# NOTE: Lower affinity CAR showed better killing than high affinity CAR
# NOTE: Data from Fig 1
Ghorashian2019 = {
    "CARS": {
        "CAT": {
            "Data": {
                "Antigens": [22307, 27],
                "Ant Err": [0, 0],
                "Kill %": [0.5, 0.05],
                "Kill % Err": [0.03, 0.01]
            },
            "Normalized Data": {
                "Antigens": [1, 0.00121],
                "Ant Err": [0, 0],
                "Kill %": [1, 0.1],
                "Kill % Err": [0.0849, 0.0209]
            },
            "KD": 1.4e-8
        },
        "FMC": {
            "Data": {
                "Antigens": [22307, 27],
                "Ant Err": [0, 0],
                "Kill %": [0.45, 0.05],
                "Kill % Err": [0.03, 0.01]
            },
            "Normalized Data": {
                "Antigens": [1, 0.00121],
                "Ant Err": [0, 0],
                "Kill %": [0.9, 0.1],
                "Kill % Err": [0.0807, 0.0209]
            },
            "KD": 3.23e-10
        }
    },
    "CITATION": "Ghorashian 2019 (CD19-41BB-CD3\u03B6, in vitro, E:T Ratio: 6.4:1)",
    "SAVE": "Ghorashian2019",
    "MARKER": "s"
}

# NOTE: Data for 1e6 and 0 line come from Fig 1, Other data from Fig 4
Caruso2015 = {
    "CARS": {
        "Cetux-CAR": {
            "Data": {
                "Antigens": [1e6, 628265, 340181, 134791, 30899, 0],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.6, 0.45, 0.35, 0.5, 0.55, 0.2],
                "Kill % Err": [0.1, 0.05, 0.2, 0.15, 0.05, 0]
            },
            "Normalized Data": {
                "Antigens": [1, 0.628, 0.340, 0.138, 0.0309, 0],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [1, 0.75, 0.583, 0.833, 0.917, 0.333],
                "Kill % Err": [0.236, 0.150, 0.347, 0.286, 0.174, 0.0556]
            },
            "KD": 8.5e-9
        },
        "Nimo-CAR": {
            "Data": {
                "Antigens": [1e6, 628265, 340181, 134791, 30899, 0],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.5, 0.42, 0.3, 0.3, 0.3, 0.2],
                "Kill % Err": [0.1, 0.05, 0.1, 0.15, 0.05, 0]
            },
            "Normalized Data": {
                "Antigens": [1, 0.628, 0.340, 0.138, 0.0309, 0],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.833, 0.7, 0.5, 0.5, 0.5, 0.333],
                "Kill % Err": [0.217, 0.143, 0.186, 0.264, 0.118, 0.0556]
            },
            "KD": 4.53e-8
        }
    },
    "CITATION": "Caruso 2015 (EGFR-CD28-CD3\u03B6, in vitro, E:T Ratio: 5:1)",
    "SAVE": "Caruso2015",
    "MARKER": "v"
}

# NOTE: Data from Fig 4. Used viability data -- kill % = 1 - viability %
# NOTE: Antigen level MFI from Fig 1.
Chmielewski2004 = {
    "CARS": {
        "C6-B1D2": {
            "Data": {
                "Antigens": [200, 25, 17, 9, 8, 7.5],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.825, 0.5, 0.3, 0, 0.2, 0],
                "Kill % Err": [0.05, 0.05, 0.05, 0.1, 0.05, 0.5]
            },
            "Normalized Data": {
                "Antigens": [1, 0.125, 0.085, 0.045, 0.04, 0.0375],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [1, 0.606, 0.364, 0.0, 0.242, 0],
                "Kill % Err": [0.0857, 0.0709, 0.0645, 0, 0.0624, 0]
            },
            "KD": 1.5e-11
        },
        "CGMH3-B1": {
            "Data": {
                "Antigens": [200, 25, 17, 9, 8, 7.5],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.825, 0.6, 0.3, 0.1, 0.2, 0.1],
                "Kill % Err": [0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
            },
            "Normalized Data": {
                "Antigens": [1, 0.125, 0.085, 0.045, 0.04, 0.0375],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [1, 0.727, 0.364, 0.121, 0.242, 0.121],
                "Kill % Err": [0.0857, 0.0749, 0.0645, 0.0611, 0.0624, 0.0610]
            },
            "KD": 1.2e-10
        },
        "C6ML3-9": {
            "Data": {
                "Antigens": [200, 25, 17, 9, 8, 7.5],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.825, 0.55, 0.3, 0.1, 0.2, 0.1],
                "Kill % Err": [0.05, 0.05, 0.05, 0.05, 0.05, 0.05]
            },
            "Normalized Data": {
                "Antigens": [1, 0.125, 0.085, 0.045, 0.04, 0.0375],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [1, 0.667, 0.364, 0.121, 0.242, 0.121],
                "Kill % Err": [0.0857, 0.0728, 0.0645, 0.0610, 0.0624, 0.0610]
            },
            "KD": 1e-9
        },
        "C6.5": {
            "Data": {
                "Antigens": [200, 25, 17, 9, 8, 7.5],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.825, 0.45, 0.2, 0.1, 0.2, 0],
                "Kill % Err": [0.05, 0.1, 0.05, 0.05, 0.05, 0.1]
            },
            "Normalized Data": {
                "Antigens": [1, 0.125, 0.085, 0.045, 0.04, 0.0375],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [1, 0.545, 0.242, 0.121, 0.242, 0],
                "Kill % Err": [0.0857, 0.126, 0.0624, 0.0610, 0.0624, 0]
            },
            "KD": 1.6e-8
        },
        "C6.5G98A": {
            "Data": {
                "Antigens": [200, 25, 17, 9, 8, 7.5],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.7, 0.05, 0, 0.1, 0.2, 0.1],
                "Kill % Err": [0.05, 0.1, 0.05, 0.05, 0.05, 0.05]
            },
            "Normalized Data": {
                "Antigens": [1, 0.125, 0.085, 0.045, 0.04, 0.0375],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.848, 0.0606, 0, 0.121, 0.242, 0.121],
                "Kill % Err": [0.0795, 0.121, 0, 0.0610, 0.0624, 0.0610]
            },
            "KD": 3.2e-7
        },
        "PBL": {
            "Data": {
                "Antigens": [200, 25, 17, 9, 8, 7.5],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.1, 0, 0, 0.1, 0, 0.1],
                "Kill % Err": [0.05, 0, 0.1, 0.05, 0, 0.05]
            },
            "Normalized Data": {
                "Antigens": [1, 0.125, 0.085, 0.045, 0.04, 0.0375],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.121, 0, 0, 0.121, 0, 0.121],
                "Kill % Err": [0.0610, 0, 0, 0.0610, 0, 0.0610]
            },
            "KD": 0
        }
    },
    "CITATION": "Chmielewski 2004 (ErbB2-Fc-CD3\u03B6, in vitro, E:T Ratio: 1:1)",
    "SAVE": "Chmielewski2004",
    "MARKER": "D"
}

# NOTE: Data from Fig 2D, chose to use 1:1 data from each plot
# NOTE: Antigens in ug ErbB2 RNA
Liu2015 = {
    "CARS": {
        "4D5.BBZ": {
            "Data": {
                "Antigens": [10, 1, 0.1],
                "Ant Err": [0, 0, 0],
                "Kill %": [0.45, 0.5, 0.175],
                "Kill % Err": [0.05, 0.05, 0.1]
            },
            "Normalized Data": {
                "Antigens": [1, 0.1, 0.01],
                "Ant Err": [0, 0, 0],
                "Kill %": [0.75, 0.833, 0.292],
                "Kill % Err": [0.0914, 0.0932, 0.0167]
            },
            "KD": 5.0e-10
        },
        "4D5-7.BBZ": {
            "Data": {
                "Antigens": [10, 1, 0.1],
                "Ant Err": [0, 0, 0],
                "Kill %": [0.5, 0.6, 0.175],
                "Kill % Err": [0, 0.03, 0.05]
            },
            "Normalized Data": {
                "Antigens": [1, 0.1, 0.01],
                "Ant Err": [0, 0, 0],
                "Kill %": [0.83, 1, 0.292],
                "Kill % Err": [0.0417, 0.0707, 0.0846]
            },
            "KD": 3.2e-9
        },
        "4D5-5.BBZ": {
            "Data": {
                "Antigens": [10, 1, 0.1],
                "Ant Err": [0, 0, 0],
                "Kill %": [0.35, 0.175, 0.03],
                "Kill % Err": [0, 0.05, 0]
            },
            "Normalized Data": {
                "Antigens": [1, 0.1, 0.01],
                "Ant Err": [0, 0, 0],
                "Kill %": [0.583, 0.292, 0.05],
                "Kill % Err": [0.0292, 0.0846, 0.0025]
            },
            "KD": 1.12e-6
        },
        "4D5-3.BBZ": {
            "Data": {
                "Antigens": [10, 1, 0.1],
                "Ant Err": [0, 0, 0],
                "Kill %": [0.3, 0.1, 0.03],
                "Kill % Err": [0.02, 0.05, 0]
            },
            "Normalized Data": {
                "Antigens": [1, 0.1, 0.01],
                "Ant Err": [0, 0, 0],
                "Kill %": [0.5, 0.167, 0.05],
                "Kill % Err": [0.0417, 0.0837, 0.0025]
            },
            "KD": 3.91e-6
        },
        "Non-EP": {
            "Data": {
                "Antigens": [10, 1, 0.1],
                "Ant Err": [0, 0, 0],
                "Kill %": [0.07, 0.02, 0.03],
                "Kill % Err": [0, 0, 0]
            },
            "Normalized Data": {
                "Antigens": [1, 0.1, 0.01],
                "Ant Err": [0, 0, 0],
                "Kill %": [0.117, 0.0333, 0.05],
                "Kill % Err": [0.00583, 0.00167, 0.0025]
            },
            "KD": 0
        },
    },
    "CITATION": "Liu 2015 (ErbB2-41BB-CD3\u03B6, in vitro, E:T Ratio: 1:1)",
    "SAVE": "Liu2015",
    "MARKER": "X"
}

# NOTE: Data from Fig 2A upper pannel
HernandezLopez2021 = {
    "CARS": {
        "Low Affinity": {
            "Data": {
                "Antigens": [1000, 50118.72, 158489.3, 501187.2, 1584893.2, 7943282.3],  # THP-1, U937, TIME, MHH-CALL-4
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0, 0.7, 0.65, 0.67, 0.78, 0.85],
                "Kill % Err": [0.15, 0, 0.05, 0.02, 0, 0]
            },
            "Normalized Data": {
                "Antigens": [0.00013, 0.0063, 0.02, 0.063, 0.2, 1],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0, 0.71, 0.66, 0.68, 0.79, 0.86],
                "Kill % Err": [0, 0, 0.05, 0.02, 0, 0]
            },
            "KD": 2.1e-7
        },
        "High Affinity": {
            "Data": {
                "Antigens": [1000, 50118.72, 158489.3, 501187.2, 1584893.2, 7943282.3],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.43, 0.98, 0.85, 0.99, 0.985, 0.9],
                "Kill % Err": [0.05, 0, 0.1, 0, 0, 0]
            },
            "Normalized Data": {
                "Antigens": [0.00013, 0.0063, 0.02, 0.063, 0.2, 1],
                "Ant Err": [0, 0, 0, 0, 0, 0],
                "Kill %": [0.43, 0.99, 0.86, 1, 0.995, 0.91],
                "Kill % Err": [0.05, 0, 0.1, 0, 0 , 0]
            },
            "KD": 1.76e-8
        },
    },
    "CITATION": "Hernandez-Lopez 2021 (HER2-41BB-CD3\u03B6, in vitro)",
    "SAVE": "HernandezLopez2021",
    "MARKER": "p"

}

litDict = [Arcangeli2017, Caruso2015, Chmielewski2004, Ghorashian2019, HernandezLopez2021, Liu2015, Wantanabe2015]

def get_max(LIST, INDEX):

    values = [i[INDEX] for i in LIST]
    maxval = max(values)

    return maxval

# ---------------- DICT AND DF FUNCTIONS ---------------------

def make_env_dict():

    envDict = {'TIME': None,
               'AXIS': None,
               'CONC': None}

    return envDict

def make_env_df():
    columns=['TIME',
             'AXIS',
             'CONC']
    envDF = pd.DataFrame(columns=columns)

    return envDF

# ---------------- PLOTTING FUNCTIONS ---------------------

# ------------- ANALYZE PLOTTING FUNCTIONS -----------------

def plot_counts(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            counts = simsDF.iloc[i][POP_NAME]
        else:
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][2:]]
            counts = simsDF.iloc[i][POP_NAME][2:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, counts,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, counts,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, counts,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS\n(NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 100000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 80000
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1400
                else:
                    ymax = 6000
            else:
                ymax = 7000
    if 'VIVO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 60000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 50000
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 4000
        else:
            ymax = 2000
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_dose(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            counts = simsDF.iloc[i][POP_NAME]
        else:
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][2:]]
            counts = simsDF.iloc[i][POP_NAME][2:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, counts,
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, counts,
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, counts,
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS\n(NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])
    ax.set_xticks(ticks=plot_time[::2])

    if 'VITRO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 100000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 80000
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1400
                else:
                    ymax = 6000
            else:
                ymax = 7000
    if 'VIVO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 60000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 50000
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 4000
        else:
            ymax = 2000
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_' + POP_NAME.replace(' ','') + '_DOSE.svg', bbox_inches='tight')

    return

def plot_counts_merge(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            counts = simsDF.iloc[i][POP_NAME]
            counts_live = simsDF.iloc[i][POP_NAME + ' LIVE']
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
            counts = simsDF.iloc[i][POP_NAME][2:]
            counts_live = simsDF.iloc[i][POP_NAME + ' LIVE'][2:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, counts,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(0)])
            ax.plot(plot_time, counts_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, counts,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict["CONTROL"])
            ax.plot(plot_time, counts_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, counts,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
            ax.plot(plot_time, counts_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS\n(NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 100000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 100000
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1400
                else:
                    ymax = 6000
            else:
                ymax = 7000
    if 'VIVO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 60000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 50000
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 4000
        else:
            ymax = 2000
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for type in liveLineDict:
        ax.plot([0],[0], color='black', label=type, linestyle=liveLineDict[type])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_MERGE_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_treat(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            counts = simsDF.iloc[i][POP_NAME]
        else:
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][44:]]
            counts = simsDF.iloc[i][POP_NAME][44:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, counts,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, counts,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, counts,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS\n(NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 100000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 80000
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1400
                else:
                    ymax = 6000
            else:
                ymax = 7000
    if 'VIVO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 60000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 50000
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 4000
        else:
            ymax = 2000
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSTREAT_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_treat_dose(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            counts = simsDF.iloc[i][POP_NAME]
        else:
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][44:]]
            counts = simsDF.iloc[i][POP_NAME][44:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, counts,
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, counts,
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, counts,
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS\n(NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 100000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 80000
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1400
                else:
                    ymax = 6000
            else:
                ymax = 7000
    if 'VIVO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 60000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 50000
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 4000
        else:
            ymax = 2000
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSTREAT_' + POP_NAME.replace(' ','') + '_DOSE.svg', bbox_inches='tight')

    return

def plot_counts_treat_merge(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            counts = simsDF.iloc[i][POP_NAME]
            counts_live = simsDF.iloc[i][POP_NAME + ' LIVE']
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][44:]]
            counts = simsDF.iloc[i][POP_NAME][44:]
            counts_live = simsDF.iloc[i][POP_NAME + ' LIVE'][44:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, counts,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(0)])
            ax.plot(plot_time, counts_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, counts,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict["CONTROL"])
            ax.plot(plot_time, counts_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, counts,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
            ax.plot(plot_time, counts_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS\n(NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 100000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 80000
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1400
                else:
                    ymax = 6000
            else:
                ymax = 7000
    if 'VIVO' in FILEID:
        if POP_NAME in ['T-CELL', 'T-CELL LIVE']:
            ymax = 60000
        elif POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
            ymax = 50000
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 4000
        else:
            ymax = 2000
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for type in liveLineDict:
        ax.plot([0],[0], color='black', label=type, linestyle=liveLineDict[type])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSTREAT_MERGE_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_norm(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE'] and simsDF.iloc[i]['DOSE'] == 0:
            continue
        elif POP_NAME in ['CD4', 'CD4 LIVE'] and TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']] == 0.0:
            continue
        elif POP_NAME in ['CD8', 'CD8 LIVE'] and TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']] == 1.0:
            continue
        else:
            if 'VITRO' in FILEID:
                plot_time = simsDF.iloc[i]['TIME']
                if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
                    frac = TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']]
                    if POP_NAME in ['CD8', 'CD8 LIVE']:
                        frac = 1 - frac
                    counts = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME]]
                elif POP_NAME in ['T-CELL', 'T-CELL LIVE']:
                    counts = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME]]
                else:
                    counts = [n/simsDF.iloc[i][POP_NAME][0] for n in simsDF.iloc[i][POP_NAME]]
            else:
                plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][2:]]
                if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
                    frac = TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']]
                    if POP_NAME in ['CD8', 'CD8 LIVE']:
                        frac = 1-frac
                    counts = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME][2:]]
                elif POP_NAME in ['T-CELL', 'T-CELL LIVE']:
                    counts = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME][2:]]
                else:
                    counts = [n/simsDF.iloc[i][POP_NAME][2] for n in simsDF.iloc[i][POP_NAME][2:]]

            if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
                ax.plot(plot_time, counts,
                        linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                        color=colorDict[str(0)])
            elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
                ax.plot(plot_time, counts,
                        linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                        color=colorDict["CONTROL"])
            else:
                ax.plot(plot_time, counts,
                        linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                        color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL\nCOUNTS NORAMLIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME NORAMLIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
            ymax = 250
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1.4
                else:
                    ymax = 6
            else:
                ymax = 3.5
    if 'VIVO' in FILEID:
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
            ymax = 300
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 1.2
        else:
            ymax = 150
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSNORM_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_norm_dose(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]
    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE'] and simsDF.iloc[i]['DOSE'] == 0:
            continue
        elif POP_NAME in ['CD4', 'CD4 LIVE'] and TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']] == 0.0:
            continue
        elif POP_NAME in ['CD8', 'CD8 LIVE'] and TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']] == 1.0:
            continue
        else:

            if 'VITRO' in FILEID:
                plot_time = simsDF.iloc[i]['TIME']
                if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
                    frac = TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']]
                    if POP_NAME in ['CD8', 'CD8 LIVE']:
                        frac = 1.0-frac
                    counts = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME]]
                elif POP_NAME in ['T-CELL', 'T-CELL LIVE']:
                    counts = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME]]
                else:
                    counts = [n/simsDF.iloc[i][POP_NAME][0] for n in simsDF.iloc[i][POP_NAME]]
            else:
                plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][2:]]
                if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
                    frac = TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']]
                    if POP_NAME in ['CD8', 'CD8 LIVE']:
                        frac = 1-frac
                    counts = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME][2:]]
                elif POP_NAME in ['T-CELL', 'T-CELL LIVE']:
                    counts = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME][2:]]
                else:
                    counts = [n/simsDF.iloc[i][POP_NAME][2] for n in simsDF.iloc[i][POP_NAME][2:]]

            if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
                ax.plot(plot_time, counts,
                        color=colorDict[str(0)])
            elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
                ax.plot(plot_time, counts,
                        color=colorDict["CONTROL"])
            else:
                ax.plot(plot_time, counts,
                        color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL\nCOUNTS NORAMLIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME NORAMLIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
            ymax = 250
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1.4
                else:
                    ymax = 6
            else:
                ymax = 3.5
    if 'VIVO' in FILEID:
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
            ymax = 300
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 1.2
        else:
            ymax = 150
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSNORM_' + POP_NAME.replace(' ','') + '_DOSE.svg', bbox_inches='tight')

    return

def plot_counts_norm_merge(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE'] and simsDF.iloc[i]['DOSE'] == 0:
            continue
        elif POP_NAME in ['CD4', 'CD4 LIVE'] and TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']] == 0.0:
            continue
        elif POP_NAME in ['CD8', 'CD8 LIVE'] and TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']] == 1.0:
            continue
        else:
            if 'VITRO' in FILEID:
                plot_time = simsDF.iloc[i]['TIME']

                if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']:
                    frac = TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']]
                    if POP_NAME in ['CD8']:
                        frac = 1-frac
                    counts = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME]]
                    counts_live = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME + ' LIVE']]
                elif POP_NAME in ['T-CELL']:
                    counts = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME]]
                    counts_live = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME + ' LIVE']]
                else:
                    counts = [n/simsDF.iloc[i][POP_NAME][0] for n in simsDF.iloc[i][POP_NAME]]
                    counts_live = [n/simsDF.iloc[i][POP_NAME + ' LIVE'][0] for n in simsDF.iloc[i][POP_NAME + ' LIVE']]
            else:
                plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
                if POP_NAME in ['CD4', 'CD8']:
                    frac = TREAT_RATIO_DICT[simsDF.iloc[i]['TREAT RATIO']]
                    if POP_NAME in ['CD8']:
                        frac = 1-frac
                    counts = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME][2:]]
                    counts_live = [n/(simsDF.iloc[i]['DOSE']*frac) for n in simsDF.iloc[i][POP_NAME + ' LIVE'][2:]]
                elif POP_NAME in ['T-CELL']:
                    counts = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME][2:]]
                    counts_live = [n/simsDF.iloc[i]['DOSE'] for n in simsDF.iloc[i][POP_NAME + ' LIVE'][2:]]
                else:
                    counts = [n/simsDF.iloc[i][POP_NAME][2] for n in simsDF.iloc[i][POP_NAME][2:]]
                    counts_live = [n/simsDF.iloc[i][POP_NAME + ' LIVE'][2] for n in simsDF.iloc[i][POP_NAME + ' LIVE'][2:]]

            if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
                ax.plot(plot_time, counts,
                        linestyle=liveLineDict['TOTAL'],
                        color=colorDict[str(0)])
                ax.plot(plot_time, counts_live,
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict[str(0)])
            elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
                ax.plot(plot_time, counts,
                        linestyle=liveLineDict['TOTAL'],
                        color=colorDict["CONTROL"])
                ax.plot(plot_time, counts_live,
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict["CONTROL"])
            else:
                ax.plot(plot_time, counts,
                        linestyle=liveLineDict['TOTAL'],
                        color=colorDict[str(simsDF.iloc[i][COLOR])])
                ax.plot(plot_time, counts_live,
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL\nCOUNTS NORAMLIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nOVER TIME NORAMLIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
            ymax = 250
        else:
            if '_CH_' in FILEID:
                if 'HEALTHY' in POP_NAME:
                    ymax = 1.4
                else:
                    ymax = 6
            else:
                ymax = 3.5
    if 'VIVO' in FILEID:
        if POP_NAME in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
            ymax = 300
        elif POP_NAME in ['HEALTHY', 'HEALTHY LIVE']:
            ymax = 1.2
        else:
            ymax = 150
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for type in liveLineDict:
        ax.plot([0],[0], color='black', label=type, linestyle=liveLineDict[type])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSNORM_MERGE_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_count_fracs(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCountsFrac = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCountsFrac.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            count_frac = simsDF.iloc[i][POP_NAME]
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
            count_frac = simsDF.iloc[i][POP_NAME][2:]

        if int(simsDF.iloc[i]['DOSE']) == 0:
            if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
                ax.plot(plot_time, count_frac,
                        linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                        color=colorDict[str(0)])
            elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
                ax.plot(plot_time, count_frac,
                        linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                        color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, count_frac,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL FRACTION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])
    ax.set_ylim([0, 1])

    if 'VITRO' in FILEID:
        ax.set_xticks(simsDF.iloc[0]['TIME'])
    else:
        ax.set_xticks(plot_time[::10])

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSFRAC_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_frac_remaining(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if 'DISH' in FILEID:
            killed = [simsDF.iloc[i][POP_NAME][a]/simsDF.iloc[i][POP_NAME][0] for a in range(0, len(simsDF.iloc[i][POP_NAME]))]
            plot_time = simsDF.iloc[i]['TIME']
        else:
            killed = [simsDF.iloc[i][POP_NAME][a]/simsDF.iloc[i][POP_NAME][44] for a in range(44, len(simsDF.iloc[i][POP_NAME]))]
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][44:]]
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, killed,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, killed,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, killed,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL\nFRACTION REMAINING", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " FRACTION\n REMAINING OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        ax.set_xticks(simsDF.iloc[0]['TIME'])
        if '_CH_' in FILEID:
            if 'HEALTHY' in POP_NAME:
                ymax = 1.4
            else:
                ymax = 6
        else:
            ymax = 4
    else:
        ax.set_xticks([int(i) for i in plot_time[::2]])
        if 'CANCER' in POP_NAME:
            ymax = 3
        else:
            ymax = 1.2
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSREMAININGNORM_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_counts_frac_remaining_dose(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if 'DISH' in FILEID:
            killed = [simsDF.iloc[i][POP_NAME][a]/simsDF.iloc[i][POP_NAME][0] for a in range(0, len(simsDF.iloc[i][POP_NAME]))]
            plot_time = simsDF.iloc[i]['TIME']
        else:
            killed = [simsDF.iloc[i][POP_NAME][a]/simsDF.iloc[i][POP_NAME][44] for a in range(44, len(simsDF.iloc[i][POP_NAME]))]
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][44:]]
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, killed,
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, killed,
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, killed,
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL\nFRACTION REMAINING", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " FRACTION\n REMAINING OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        ax.set_xticks(simsDF.iloc[0]['TIME'])
        if '_CH_' in FILEID:
            if 'HEALTHY' in POP_NAME:
                ymax = 1.4
            else:
                ymax = 6
        else:
            ymax = 4
    else:
        ax.set_xticks([int(i) for i in plot_time[::2]])
        if 'CANCER' in POP_NAME:
            ymax = 3
        else:
            ymax = 1.2
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSREMAININGNORM_' + POP_NAME.replace(' ','') + '_DOSE.svg', bbox_inches='tight')

    return

def plot_counts_frac_remaining_merge(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if 'DISH' in FILEID:
            killed = [simsDF.iloc[i][POP_NAME][a]/simsDF.iloc[i][POP_NAME][0] for a in range(0, len(simsDF.iloc[i][POP_NAME]))]
            killed_live = [simsDF.iloc[i][POP_NAME + ' LIVE'][a] / simsDF.iloc[i][POP_NAME + ' LIVE'][0] for a in
                      range(0, len(simsDF.iloc[i][POP_NAME + ' LIVE']))]
            plot_time = simsDF.iloc[i]['TIME']
        else:
            killed = [simsDF.iloc[i][POP_NAME][a]/simsDF.iloc[i][POP_NAME][44] for a in range(44, len(simsDF.iloc[i][POP_NAME]))]
            killed_live = [simsDF.iloc[i][POP_NAME + ' LIVE'][a] / simsDF.iloc[i][POP_NAME + ' LIVE'][44] for a in
                      range(44, len(simsDF.iloc[i][POP_NAME + ' LIVE']))]
            plot_time = [t-1 for t in simsDF.iloc[i]['TIME'][44:]]
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, killed,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(0)])
            ax.plot(plot_time, killed_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, killed,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict["CONTROL"])
            ax.plot(plot_time, killed_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, killed,
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
            ax.plot(plot_time, killed_live,
                    linestyle=liveLineDict['LIVE'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL\nFRACTION REMAINING", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " FRACTION\n REMAINING OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        ax.set_xticks(simsDF.iloc[0]['TIME'])
        if '_CH_' in FILEID:
            if 'HEALTHY' in POP_NAME:
                ymax = 1.4
            else:
                ymax = 6
        else:
            ymax = 4
    else:
        ax.set_xticks([int(i) for i in plot_time[::2]])
        if 'CANCER' in POP_NAME:
            ymax = 3
        else:
            ymax = 1.2
    ax.set_ylim(bottom=0, top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for type in liveLineDict:
        ax.plot([0],[0], color='black', label=type, linestyle=liveLineDict[type])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSREMAININGNORM_MERGE_' + POP_NAME.replace(' ','') + '.svg', bbox_inches='tight')

    return

# STUB
def plot_states(simsDF, COLOR, FILEID, SAVELOC):

    return

def plot_state_fracs(simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figStateFracs, axs = plt.subplots(12, 7, figsize=(10, 15))

    for i in range(0, len(simsDF)):

        if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
            color = 0
        elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
            color = "CONTROL"
        else:
            color = simsDF.iloc[i][COLOR]

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            start = 0
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
            start = 2

        # Plot CANCER cell state fracs
        axs[0, 0].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 0].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 0].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 0].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 0].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 0].plot(plot_time, simsDF.iloc[i]['NECRO ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 0].set_visible(False)
        axs[7, 0].set_visible(False)
        axs[8, 0].set_visible(False)
        axs[9, 0].set_visible(False)
        axs[10, 0].set_visible(False)
        axs[11, 0].set_visible(False)

        # Plot CANCER LIVE cell state fracs
        axs[0, 1].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 1].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 1].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 1].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 1].set_visible(False)
        axs[5, 1].set_visible(False)
        axs[6, 1].set_visible(False)
        axs[7, 1].set_visible(False)
        axs[8, 1].set_visible(False)
        axs[9, 1].set_visible(False)
        axs[10, 1].set_visible(False)
        axs[11, 1].set_visible(False)

        # Plot HEALTHY cell state fracs
        axs[0, 2].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 2].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 2].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 2].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[4, 2].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[5, 2].plot(plot_time, simsDF.iloc[i]['NECRO ' + POP_NAMES[2] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[6, 2].set_visible(False)
        axs[7, 2].set_visible(False)
        axs[8, 2].set_visible(False)
        axs[9, 2].set_visible(False)
        axs[10, 2].set_visible(False)
        axs[11, 2].set_visible(False)

        # Plot HEALTHY LIVE cell state fracs
        axs[0, 3].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[3] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[1, 3].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[3] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[2, 3].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[3] + ' %'][start:],
                       color=colorDict[str(simsDF.iloc[i][COLOR])])
        axs[3, 3].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[3] + ' %'][start:],
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
        axs[0, 4].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 4].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 4].set_visible(False)
        axs[3, 4].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 4].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 4].set_visible(False)
        axs[6, 4].plot(plot_time, simsDF.iloc[i]['CYTOT ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[7, 4].plot(plot_time, simsDF.iloc[i]['STIMU ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 4].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[9, 4].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 4].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 4].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])

        # Plot CD4 T-CELL state fracs
        axs[0, 5].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 5].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 5].set_visible(False)
        axs[3, 5].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 5].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 5].set_visible(False)
        axs[6, 5].set_visible(False)
        axs[7, 5].plot(plot_time, simsDF.iloc[i]['STIMU ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 5].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[9, 5].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 5].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 5].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])

        # Plot CD8 T-CELL state fracs
        axs[0, 6].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 6].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 6].set_visible(False)
        axs[3, 6].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 6].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 6].set_visible(False)
        axs[6, 6].plot(plot_time, simsDF.iloc[i]['CYTOT ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[7, 6].set_visible(False)
        axs[8, 6].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[9, 6].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 6].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 6].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])

    for a in range(0, 12):
        for b in range(0, 7):
            axs[a, b].set_xlim([0, max(plot_time)])
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

        if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
            color = 0
        elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
            color = "CONTROL"
        else:
            color = simsDF.iloc[i][COLOR]

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            start = 0
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
            start = 2


        # Plot CANCER cell state fracs
        axs[0, 0].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 0].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 0].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 0].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 0].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 0].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 0].plot(plot_time, simsDF.iloc[i]['NECRO ' + POP_NAMES[0] + ' %'][start:], color=colorDict[str(color)])
        axs[7, 0].set_visible(False)
        axs[8, 0].set_visible(False)
        axs[9, 0].set_visible(False)
        axs[10, 0].set_visible(False)
        axs[11, 0].set_visible(False)
        axs[12, 0].set_visible(False)

        # Plot CANCER LIVE cell state fracs
        axs[0, 1].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 1].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 1].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 1].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 1].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[1] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 1].set_visible(False)
        axs[6, 1].set_visible(False)
        axs[7, 1].set_visible(False)
        axs[8, 1].set_visible(False)
        axs[9, 1].set_visible(False)
        axs[10, 1].set_visible(False)
        axs[11, 1].set_visible(False)
        axs[12, 1].set_visible(False)

        # Plot HEALTHY cell state fracs
        axs[0, 2].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 2].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 2].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 2].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 2].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 2].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 2].plot(plot_time, simsDF.iloc[i]['NECRO ' + POP_NAMES[2] + ' %'][start:], color=colorDict[str(color)])
        axs[7, 2].set_visible(False)
        axs[8, 2].set_visible(False)
        axs[9, 2].set_visible(False)
        axs[10, 2].set_visible(False)
        axs[11, 2].set_visible(False)
        axs[12, 2].set_visible(False)

        # Plot HEALTHY LIVE cell state fracs
        axs[0, 3].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[3] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 3].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[3] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 3].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[3] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 3].plot(plot_time, simsDF.iloc[i]['QUIES ' + POP_NAMES[3] + ' %'][start:], color=colorDict[str(color)])
        axs[4, 3].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[3] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 3].set_visible(False)
        axs[6, 3].set_visible(False)
        axs[7, 3].set_visible(False)
        axs[8, 3].set_visible(False)
        axs[9, 3].set_visible(False)
        axs[10, 3].set_visible(False)
        axs[11, 3].set_visible(False)
        axs[12, 3].set_visible(False)

        # Plot T-CELL state fracs
        axs[0, 4].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 4].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 4].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 4].set_visible(False)
        axs[4, 4].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 4].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 4].set_visible(False)
        axs[7, 4].plot(plot_time, simsDF.iloc[i]['CYTOT ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 4].plot(plot_time, simsDF.iloc[i]['STIMU ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 4].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 4].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 4].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])
        axs[12, 4].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[4] + ' %'][start:], color=colorDict[str(color)])

        # Plot CD4 T-CELL state fracs
        axs[0, 5].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 5].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 5].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 5].set_visible(False)
        axs[4, 5].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 5].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 5].set_visible(False)
        axs[7, 5].set_visible(False)
        axs[8, 5].plot(plot_time, simsDF.iloc[i]['STIMU ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 5].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 5].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 5].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])
        axs[12, 5].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[6] + ' %'][start:], color=colorDict[str(color)])

        # Plot CD8 T-CELL state fracs
        axs[0, 6].plot(plot_time, simsDF.iloc[i]['NEUTR ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[1, 6].plot(plot_time, simsDF.iloc[i]['MIGRA ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[2, 6].plot(plot_time, simsDF.iloc[i]['PROLI ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[3, 6].set_visible(False)
        axs[4, 6].plot(plot_time, simsDF.iloc[i]['SENES ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[5, 6].plot(plot_time, simsDF.iloc[i]['APOPT ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[6, 6].set_visible(False)
        axs[7, 6].plot(plot_time, simsDF.iloc[i]['CYTOT ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[8, 6].set_visible(False)
        axs[9, 6].plot(plot_time, simsDF.iloc[i]['EXHAU ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[10, 6].plot(plot_time, simsDF.iloc[i]['ANERG ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[11, 6].plot(plot_time, simsDF.iloc[i]['STARV ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])
        axs[12, 6].plot(plot_time, simsDF.iloc[i]['PAUSE ' + POP_NAMES[8] + ' %'][start:], color=colorDict[str(color)])

    for a in range(0, 13):
        for b in range(0, 7):
            axs[a, b].set_xlim([0, max(plot_time)])
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

def plot_volumes(simsDF, COLOR, FILEID, SAVELOC, TIME):

    colorDict = COLOR_DICT[COLOR]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    filesplit = FILEID.split('_')

    for p in range(0,len(POP_NAMES)):
        if filesplit[8] == 'NA' and 'HEALTHY' in POP_NAMES[p]:
            continue
        else:
            if 'LIVE' not in POP_NAMES[p]:
                key = 'CELL VOLUMES ' + POP_NAMES[p]
                Y = []
                X = []
                order = []

                for i in range(0, len(simsDF)):
                    y = simsDF.iloc[i][key][index]

                    if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
                        x = [str(0)] * len(y)
                    elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
                        x = ["CONTROL"] * len(y)
                    else:
                        x = [str(simsDF.iloc[i][COLOR])] * len(y)

                    Y = Y + y
                    X = X + x

                for key in colorDict:
                    if key in X:
                        order.append(key)

                if X != [] and Y != []:
                    figV = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
                    ax = figV.add_subplot(1, 1, 1)
                    ax = sns.violinplot(x=X, y=Y, order=order, palette=colorDict)
                    ax.set_xlabel(POP_NAMES[p] + " POPULATION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES,
                                  labelpad=LABELPAD)
                    ax.set_ylabel("VOLUME (um^3)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
                    ax.set_title(POP_NAMES[p] + " VOLUME\nDISTRIBUTIONS AT TIME " + str(TIME), fontname='Arial',
                                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
                    # ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
                    ymax = 7000
                    if POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                        ymax = 700
                    ax.set_ylim(bottom=0, top=ymax)

                    plt.xticks(fontsize=TICKSIZE, rotation=45)
                    plt.yticks(fontsize=TICKSIZE)

                    if SAVELOC == '':
                        plt.show()
                    else:
                        plt.savefig(SAVELOC + FILEID.replace('_VOLUMES','') + '_VOLUMES_' + POP_NAMES[p].replace(' ', '') + '_DAY_' + str(TIME) + '.svg', bbox_inches='tight')

    return

def plot_cycles(simsDF, COLOR, FILEID, SAVELOC, TIME, HOURS):
    colorDict = COLOR_DICT[COLOR]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    filesplit = FILEID.split('_')

    for p in range(0, len(POP_NAMES)):
        if filesplit[8] == 'NA' and 'HEALTHY' in POP_NAMES[p]:
            continue
        else:
            if 'LIVE' not in POP_NAMES[p]:
                key = 'AVG CELL CYCLES ' + POP_NAMES[p]
                Y = []
                X = []
                order = []

                for i in range(0, len(simsDF)):
                    y = simsDF.iloc[i][key][index]

                    if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
                        x = [str(0)] * len(y)
                    elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
                        x = ["CONTROL"] * len(y)
                    else:
                        x = [str(simsDF.iloc[i][COLOR])] * len(y)

                    Y = Y + y
                    X = X + x

                for key in colorDict:
                    if key in X:
                        order.append(key)

                if X != [] and Y != []:

                    if HOURS:
                        Y = [float(c)/float(60) for c in Y]

                    figC = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
                    ax = figC.add_subplot(1, 1, 1)
                    ax = sns.violinplot(x=X, y=Y, order=order, palette=colorDict)
                    unit = "(hrs)" if HOURS else "(min)"
                    ax.set_xlabel(POP_NAMES[p] + " POPULATION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
                    ax.set_ylabel("CYCLE LENGTH " + unit, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
                    ax.set_title(POP_NAMES[p] + " CELL CYCLE\nDISTRIBUTIONS AT TIME " + str(TIME), fontname='Arial',
                                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
                    # ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
                    if not HOURS:
                        if 'VITRO' in FILEID:
                            if POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                                ymax = 1600
                            else:
                                ymax = 2000
                        else:
                            ymax = 2000
                    else:
                        if 'VITRO' in FILEID:
                            if POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                                ymax = 25
                            else:
                                ymax = 35
                        else:
                            ymax = 35
                    ax.set_ylim(bottom=0, top=ymax)
                    plt.xticks(fontsize=TICKSIZE, rotation=45)
                    plt.yticks(fontsize=TICKSIZE)
                    if SAVELOC == '':
                        plt.show()
                    else:
                        if HOURS:
                            units = '_HOURS'
                        else:
                            units = '_MINUTES'
                        plt.savefig(SAVELOC + FILEID.replace('_CYCLES','') + '_CYCLES_' + POP_NAMES[p].replace(' ', '') + '_DAY_' + str(TIME) + units + '.svg', bbox_inches='tight')

    return

def plot_volumes_split(simsDF, COLOR, FILEID, SAVELOC, TIME):

    colorDict = COLOR_DICT[COLOR]

    filesplit = FILEID.split('_')

    for p in range(0,len(POP_NAMES)):
        if filesplit[8] == 'NA' and 'HEALTHY' in POP_NAMES[p]:
            continue
        else:
            if 'LIVE' not in POP_NAMES[p]:

                TIMES_SIM = simsDF.iloc[0]['TIME']
                index = -1

                key = 'CELL VOLUMES ' + POP_NAMES[p]
                Y = []
                X = []
                order = []
                T = []

                for time in TIME:
                    for t in range(0, len(TIMES_SIM)):
                        if float(TIMES_SIM[t]) == float(time):
                            index = t
                            for i in range(0, len(simsDF)):
                                y = simsDF.iloc[i][key][index]
                                times = [time] * len(y)
                                if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
                                    x = [str(0)] * len(y)
                                elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
                                    x = ["CONTROL"] * len(y)
                                else:
                                    x = [str(simsDF.iloc[i][COLOR])] * len(y)

                                Y = Y + y
                                X = X + x
                                T = T + times

                dictVols = {'TIME': T, 'FEATURE': X, 'VOLUME': Y}
                dfVols = pd.DataFrame(dictVols, columns={'TIME': pd.Series([], dtype='float'),
                                                         'FEATURE': pd.Series([], dtype='str'),
                                                         'VOLUME': pd.Series([], dtype='float')})
                for key in colorDict:
                    if key in X:
                        order.append(key)

                if X != [] and Y != []:
                    figV = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
                    ax = figV.add_subplot(1, 1, 1)
                    sns.violinplot(ax=ax, data=dfVols, x="FEATURE", y="VOLUME", order=order, hue="TIME", split=True)#, palette=colorDict)
                    ax.set_xlabel(POP_NAMES[p] + " POPULATION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES,
                                  labelpad=LABELPAD)
                    ax.set_ylabel("VOLUME (um^3)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
                    ax.set_title(POP_NAMES[p] + " VOLUME\nDISTRIBUTIONS AT TIME " + str(TIME[0]) + " AND " + str(TIME[1]), fontname='Arial',
                                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
                    ax.legend(bbox_to_anchor=(1.15, 1.0), frameon=False)
                    ymax = 7000
                    if POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                        ymax = 700
                    ax.set_ylim(bottom=0, top=ymax)

                    plt.xticks(fontsize=TICKSIZE, rotation=45)
                    plt.yticks(fontsize=TICKSIZE)

                    if SAVELOC == '':
                        plt.show()
                    else:
                        plt.savefig(SAVELOC + FILEID.replace('_VOLUMES','') + '_VOLUMES_' + POP_NAMES[p].replace(' ', '') + '_DAY_' + str(TIME[0]) + str(TIME[1]) + '.svg', bbox_inches='tight')

    return

def plot_cycles_split(simsDF, COLOR, FILEID, SAVELOC, TIME, HOURS):

    colorDict = COLOR_DICT[COLOR]

    filesplit = FILEID.split('_')

    for p in range(0,len(POP_NAMES)):
        if filesplit[8] == 'NA' and 'HEALTHY' in POP_NAMES[p]:
            continue
        else:
            if 'LIVE' not in POP_NAMES[p]:

                TIMES_SIM = simsDF.iloc[0]['TIME']
                index = -1

                key = 'AVG CELL CYCLES ' + POP_NAMES[p]
                Y = []
                X = []
                order = []
                T = []

                for time in TIME:
                    for t in range(0, len(TIMES_SIM)):
                        if float(TIMES_SIM[t]) == float(time):
                            index = t
                            for i in range(0, len(simsDF)):
                                y = simsDF.iloc[i][key][index]
                                times = [time] * len(y)
                                if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
                                    x = [str(0)] * len(y)
                                elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
                                    x = ["CONTROL"] * len(y)
                                else:
                                    x = [str(simsDF.iloc[i][COLOR])] * len(y)

                                Y = Y + y
                                X = X + x
                                T = T + times
                if HOURS:
                    Y = [float(c) / float(60) for c in Y]

                dictVols = {'TIME': T, 'FEATURE': X, 'CYCLE': Y}
                dfVols = pd.DataFrame(dictVols, columns={'TIME': pd.Series([], dtype='float'),
                                                         'FEATURE': pd.Series([], dtype='str'),
                                                         'CYCLE': pd.Series([], dtype='float')})
                for key in colorDict:
                    if key in X:
                        order.append(key)

                if X != [] and Y != []:
                    figV = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
                    ax = figV.add_subplot(1, 1, 1)
                    sns.violinplot(ax=ax, data=dfVols, x="FEATURE", y="CYCLE", order=order, hue="TIME", split=True)#, palette=colorDict)
                    ax.set_xlabel(POP_NAMES[p] + " POPULATION", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES,
                                  labelpad=LABELPAD)
                    unit = "(hrs)" if HOURS else "(min)"
                    ax.set_ylabel("CYCLE LENGTH " + unit, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
                    ax.set_title(POP_NAMES[p] + " CELL CYCLE\nDISTRIBUTIONS AT TIME " + str(TIME[0]) + " AND " + str(TIME[1]), fontname='Arial',
                                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
                    ax.legend(bbox_to_anchor=(1.15, 1.0), frameon=False)

                    if not HOURS:
                        if 'VITRO' in FILEID:
                            if POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                                ymax = 1600
                            else:
                                ymax = 2000
                        else:
                            ymax = 2000
                    else:
                        if 'VITRO' in FILEID:
                            if POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                                ymax = 25
                            else:
                                ymax = 35
                        else:
                            ymax = 35
                    ax.set_ylim(bottom=0, top=ymax)

                    plt.xticks(fontsize=TICKSIZE, rotation=45)
                    plt.yticks(fontsize=TICKSIZE)

                    if SAVELOC == '':
                        plt.show()
                    else:
                        if HOURS:
                            units = '_HOURS'
                        else:
                            units = '_MINUTES'
                        plt.savefig(SAVELOC + FILEID.replace('_CYCLES','') + '_CYCLES_' + POP_NAMES[p].replace(' ', '') + '_DAY_' + str(TIME[0]) + str(TIME[1]) + units + '.svg', bbox_inches='tight')

    return

def plot_kill_curve_sim(simsDF, FILEID, SAVELOC, TIME):

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0,len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figKC = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
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

    x = [0,10,20,30,40,50]
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
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("ANTIGENS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("% LYSIS\n(1 - CANCER LIVE %)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylim([0, 1])
    ax.set_xlim([0, 10000])
    ax.legend(bbox_to_anchor=(1, 1), frameon=False)

    plt.xticks(fontsize=TICKSIZE, rotation=45)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_KILLCURVESIM_' + str(TIME) + '.svg', bbox_inches='tight')

    return

def plot_kill_curve_relative_sim(simsDF, FILEID, SAVELOC, TIME):

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0,len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figKC = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
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
        if 'DISH' in FILEID:
            killed = 1-(simsDF.iloc[i]['CANCER LIVE'][index]/simsDF.iloc[i]['CANCER LIVE'][0])
        else:
            killed = 1-(simsDF.iloc[i]['CANCER LIVE'][index]/simsDF.iloc[i]['CANCER LIVE'][44])

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

    x = [0,10,20,30,40,50]
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
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("ANTIGENS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("% LYSIS\n[1 - (CANCER LIVE END/CANCER LIVE SEEDED)]", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylim([0, 1])
    ax.set_xlim([0, 10000])
    ax.legend(bbox_to_anchor=(1, 1), frameon=False)

    plt.xticks(fontsize=TICKSIZE, rotation=45)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_KILLCURVEREALATIVESIM_' + str(TIME) + '.svg', bbox_inches='tight')

    return

def plot_kill_curve_normalized_sim(simsDF, FILEID, SAVELOC, TIME):

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0,len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figKC = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
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
        killed = 1-(simsDF.iloc[i]['CANCER LIVE'][index]/simsDF.iloc[i]['CANCER LIVE'][0])

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

    x = [0,10,20,30,40,50]
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

    ax.plot([a[0]/get_max(KD6D250AVG, 0) for a in KD6D250AVG], [a[1]/get_max(KD6D250AVG, 1) for a in KD6D250AVG],
            linestyle=doseLineDict['250'], color=COLOR_DICT["CAR AFFINITY"]["1e-06"], linewidth=4)
    ax.plot([a[0]/get_max(KD7D250AVG, 0) for a in KD7D250AVG], [a[1]/get_max(KD7D250AVG, 1) for a in KD7D250AVG],
            linestyle=doseLineDict['250'],color=COLOR_DICT["CAR AFFINITY"]["1e-07"], linewidth=4)
    ax.plot([a[0]/get_max(KD8D250AVG, 0) for a in KD8D250AVG], [a[1]/get_max(KD8D250AVG, 1) for a in KD8D250AVG],
            linestyle=doseLineDict['250'], color=COLOR_DICT["CAR AFFINITY"]["1e-08"], linewidth=4)
    ax.plot([a[0]/get_max(KD9D250AVG, 0) for a in KD9D250AVG], [a[1]/get_max(KD9D250AVG, 1) for a in KD9D250AVG],
            linestyle=doseLineDict['250'], color=COLOR_DICT["CAR AFFINITY"]["1e-09"], linewidth=4)

    ax.plot([a[0]/get_max(KD6D500AVG, 0) for a in KD6D500AVG], [a[1]/get_max(KD6D500AVG, 1) for a in KD6D500AVG],
            linestyle=doseLineDict['500'], color=COLOR_DICT["CAR AFFINITY"]["1e-06"], linewidth=4)
    ax.plot([a[0]/get_max(KD7D500AVG, 0) for a in KD7D500AVG], [a[1]/get_max(KD7D500AVG, 1) for a in KD7D500AVG],
            linestyle=doseLineDict['500'], color=COLOR_DICT["CAR AFFINITY"]["1e-07"], linewidth=4)
    ax.plot([a[0]/get_max(KD8D500AVG, 0) for a in KD8D500AVG], [a[1]/get_max(KD8D500AVG, 1) for a in KD8D500AVG],
            linestyle=doseLineDict['500'], color=COLOR_DICT["CAR AFFINITY"]["1e-08"], linewidth=4)
    ax.plot([a[0]/get_max(KD9D500AVG, 0) for a in KD9D500AVG], [a[1]/get_max(KD9D500AVG, 1) for a in KD9D500AVG],
            linestyle=doseLineDict['500'], color=COLOR_DICT["CAR AFFINITY"]["1e-09"], linewidth=4)

    ax.plot([a[0]/get_max(KD6D1000AVG, 0) for a in KD6D1000AVG], [a[1]/get_max(KD6D1000AVG, 1) for a in KD6D1000AVG],
            linestyle=doseLineDict['1000'], color=COLOR_DICT["CAR AFFINITY"]["1e-06"], linewidth=4)
    ax.plot([a[0]/get_max(KD7D1000AVG, 0) for a in KD7D1000AVG], [a[1]/get_max(KD7D1000AVG, 1) for a in KD7D1000AVG],
            linestyle=doseLineDict['1000'], color=COLOR_DICT["CAR AFFINITY"]["1e-07"], linewidth=4)
    ax.plot([a[0]/get_max(KD8D1000AVG, 0) for a in KD8D1000AVG], [a[1]/get_max(KD8D1000AVG, 1) for a in KD8D1000AVG],
            linestyle=doseLineDict['1000'], color=COLOR_DICT["CAR AFFINITY"]["1e-08"], linewidth=4)
    ax.plot([a[0]/get_max(KD9D1000AVG, 0) for a in KD9D1000AVG], [a[1]/get_max(KD9D1000AVG, 1) for a in KD9D1000AVG],
            linestyle=doseLineDict['1000'], color=COLOR_DICT["CAR AFFINITY"]["1e-09"], linewidth=4)

    for color in affinityColorDict:
        ax.plot([0],[0], color=affinityColorDict[color], label=color, linestyle='solid', linewidth=4)
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose], linewidth=4)
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    ax.set_title("KILL CURVE (SIMULATED DATA)", fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("ANTIGENS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("% LYSIS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylim([-2.5, 1])
    ax.set_xlim([0, 1])
    ax.legend(bbox_to_anchor=(1, 1), frameon=False)

    plt.xticks([0.0, 0.5, 1.0], fontsize=TICKSIZE, rotation=45)
    plt.yticks([0.0, 0.5, 1.0], fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_KILLCURVESIMNORM_' + str(TIME) + '.svg', bbox_inches='tight')

    return

def plot_kill_curve_exp(SAVELOC):

    figKC = plt.figure()
    ax = figKC.add_subplot(1, 1, 1)

    # Plot literature values
    for d in litDict:
        for CAR in d["CARS"]:
            if 0 >= d["CARS"][CAR]["KD"] >= 1e-5:
                d["CARS"][CAR].update({"Color": "grey"})
            elif 1e-5 > d["CARS"][CAR]["KD"] >= 1e-6:
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

    ax.set_title("KILL CURVE\nEXPERIMENTAL DATA)", fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("ANTIGENS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("% LYSIS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylim([0, 1])
    ax.set_xlim([0, 10000])
    ax.legend(bbox_to_anchor=(1, 0.40), frameon=False)

    plt.xticks(fontsize=TICKSIZE, rotation=45)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'KILLCURVEEXP.svg', bbox_inches='tight')

    return

def plot_kill_curve_exp_separated(SAVELOC):


    # Plot literature values
    for d in litDict:

        figKC = plt.figure()
        ax = figKC.add_subplot(1, 1, 1)

        for CAR in d["CARS"]:
            if 0 >= d["CARS"][CAR]["KD"] >= 1e-5:
                d["CARS"][CAR].update({"Color": "grey"})
            elif 1e-5 > d["CARS"][CAR]["KD"] >= 1e-6:
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

        ax.set_title("KILL CURVE\n(EXPERIMENTAL DATA)\n" + d["SAVE"], fontname='Arial',
                     fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
        ax.set_xlabel("ANTIGENS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_ylabel("% LYSIS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_ylim([0, 1])
        ax.set_xlim([0, 10000])
        ax.legend(bbox_to_anchor=(1, 0.40), frameon=False)

        plt.xticks(fontsize=TICKSIZE, rotation=45)
        plt.yticks(fontsize=TICKSIZE)

        if SAVELOC == '':
            plt.show()
        else:
            plt.savefig(SAVELOC + 'KILLCURVEEXP' + d["SAVE"] + '.svg', bbox_inches='tight')

    return

def plot_kill_curve_normalized_exp(SAVELOC):

    figKC = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figKC.add_subplot(1, 1, 1)

    # Plot literature values
    for d in litDict:
        for CAR in d["CARS"]:
            if 0 >= d["CARS"][CAR]["KD"] >= 1e-5:
                d["CARS"][CAR].update({"Color": "grey"})
            elif 1e-5 > d["CARS"][CAR]["KD"] >= 1e-6:
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

            ax.errorbar(d["CARS"][CAR]["Normalized Data"]["Antigens"], d["CARS"][CAR]["Normalized Data"]["Kill %"],
                        xerr=d["CARS"][CAR]["Normalized Data"]["Ant Err"], fmt='None', yerr=d["CARS"][CAR]["Normalized Data"]["Kill % Err"],
                        ecolor='lightgray', zorder=0)
            # ax.scatter(d["CARS"][CAR]["Data"]["Antigens"], d["CARS"][CAR]["Data"]["Kill %"],
            #            marker=d["MARKER"], color=d["CARS"][CAR]["Color"], zorder=1)
            ax.plot(d["CARS"][CAR]["Normalized Data"]["Antigens"], d["CARS"][CAR]["Normalized Data"]["Kill %"],
                        marker=d["MARKER"], color=d["CARS"][CAR]["Color"], zorder=1)

        ax.scatter(None, None, marker=d["MARKER"], label=d["CITATION"], color='black')

    for color in affinityColorDict:
        ax.plot([0],[0], color=affinityColorDict[color], label=color, linestyle='solid')

    ax.set_title("NORMALIZED KILL CURVE\n(EXPERIMENTAL DATA)", fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("NORMALIZED ANTIGENS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("NORMALIZED % LYSIS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylim([0, 1])
    ax.set_xlim([0, 1])
    ax.legend(bbox_to_anchor=(1, 0.40), frameon=False)

    plt.xticks(fontsize=TICKSIZE, rotation=45)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'KILLCURVENORMEXP.svg', bbox_inches='tight')

    return

def plot_kill_curve_normalized_exp_separated(SAVELOC):

    # Plot literature values
    for d in litDict:

        figKC = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
        ax = figKC.add_subplot(1, 1, 1)

        for CAR in d["CARS"]:
            if 0 >= d["CARS"][CAR]["KD"] >= 1e-5:
                d["CARS"][CAR].update({"Color": "grey"})
            elif 1e-5 > d["CARS"][CAR]["KD"] >= 1e-6:
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

            ax.errorbar(d["CARS"][CAR]["Normalized Data"]["Antigens"], d["CARS"][CAR]["Normalized Data"]["Kill %"],
                        xerr=d["CARS"][CAR]["Normalized Data"]["Ant Err"], fmt='None', yerr=d["CARS"][CAR]["Normalized Data"]["Kill % Err"],
                        ecolor='lightgray', zorder=0, linewidth=4)
            # ax.scatter(d["CARS"][CAR]["Data"]["Antigens"], d["CARS"][CAR]["Data"]["Kill %"],
            #            marker=d["MARKER"], color=d["CARS"][CAR]["Color"], zorder=1)
            ax.plot(d["CARS"][CAR]["Normalized Data"]["Antigens"], d["CARS"][CAR]["Normalized Data"]["Kill %"],
                        marker=d["MARKER"], color=d["CARS"][CAR]["Color"], zorder=1, linewidth=4, markersize=15)

        ax.scatter(None, None, marker=d["MARKER"], label=d["CITATION"], color='black', linewidth=4)

        for color in affinityColorDict:
            ax.plot([0],[0], color=affinityColorDict[color], label=color, linestyle='solid', linewidth=4)

        ax.set_title("NORMALIZED KILL CURVE\n(EXPERIMENTAL DATA)\n" + d["SAVE"], fontname='Arial',
                     fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
        ax.set_xlabel("NORMALIZED ANTIGENS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_ylabel("NORMALIZED % LYSIS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_ylim([0, 1])
        ax.set_xlim([0, 1])
        ax.legend(bbox_to_anchor=(1, 0.40), frameon=False)

        plt.xticks([0.0, 0.5, 1.0], fontsize=TICKSIZE, rotation=45)
        plt.yticks([0.0, 0.5, 1.0], fontsize=TICKSIZE)

        if SAVELOC == '':
            plt.show()
        else:
            plt.savefig(SAVELOC + 'KILLCURVENORMEXP' + d["SAVE"] + '.svg', bbox_inches='tight')

    return

# FIX? MAYBE USE LYSIS FILE?
def plot_ET_ratio(simsDF, FILEID, SAVELOC):
    figKC = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
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
        if 'DISH' in FILEID:
            init_cancer = simsDF.iloc[i]['CANCER'][0]
        else:
            init_cancer = simsDF.iloc[i]['CANCER'][44]
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

#FIX
def plot_ET_ratio_relative(simsDF, FILEID, SAVELOC):
    figKC = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
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
        if 'DISH' in FILEID:
            killed = 1 - (simsDF.iloc[i]['CANCER LIVE'][-1]/simsDF.iloc[i]['CANCER LIVE'][0])
            init_cancer = simsDF.iloc[i]['CANCER'][0]
        else:
            killed = 1 - (simsDF.iloc[i]['CANCER LIVE'][-1] / simsDF.iloc[i]['CANCER LIVE'][44])
            init_cancer = simsDF.iloc[i]['CANCER'][44]
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
        plt.savefig(SAVELOC + FILEID + '_ETRATIORELATIVE.svg', bbox_inches='tight')

    return

def plot_CH_scatter(simsDF, COLOR, MARKER, FILEID, SAVELOC, TIME):
    colorDict = COLOR_DICT[COLOR]
    markerDict = MARKER_DICT[MARKER]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figCH = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCH.add_subplot(1, 1, 1)

    for i in range(0, len(simsDF)):
        cancer = simsDF.iloc[i]['CANCER LIVE'][index]
        healthy = simsDF.iloc[i]['HEALTHY LIVE'][index]
        if int(simsDF.iloc[i]['DOSE']) == 0:
            color = 'black'
        else:
            color = colorDict[str(simsDF.iloc[i][COLOR])]
        if MARKER != 'o':
            marker = markerDict[simsDF.iloc[i][MARKER]]
        else:
            marker = markerDict[MARKER]

        ax.scatter(cancer, healthy, color=color, marker=marker, s=150) # facecolor='none', edgecolor=color, marker=marker, s=200, linewidth=2)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    if MARKER != 'o':
        for marker in markerDict:
            ax.plot([0],[0], color='black', label=marker, linestyle=markerDict[marker])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)
    ax.set_title("LIVING CANCER VS HEALTHY CELLS", fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("LIVING CANCER CELL COUNT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("LIVING HEALTHY CELL COUNT", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_SCATTER_CH_' + COLOR.replace(' ','') + '_' + MARKER.replace(' ','') + '_' + str(TIME) + '.svg', bbox_inches='tight')

    return

def plot_CH_scatter_normalized(simsDF, COLOR, MARKER, FILEID, SAVELOC, TIME):
    colorDict = COLOR_DICT[COLOR]
    markerDict = MARKER_DICT[MARKER]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figCH = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCH.add_subplot(1, 1, 1)

    if 'VITRO' in FILEID:
        xmax = 6
        if COLOR == 'CAR AFFINITY' and '_CH_' in FILEID:
            ymax = 1.4
        elif COLOR == 'ANTIGENS CANCER' and '_CH_' in FILEID:
            ymax = 1.4
        else:
            ymax = 1.4
    else:
        xmax = 2.5
        ymax = 1


    for i in range(0, len(simsDF)):
        if 'DISH' in FILEID:
            cancer = simsDF.iloc[i]['CANCER LIVE'][index]/simsDF.iloc[i]['CANCER LIVE'][0]
            healthy = simsDF.iloc[i]['HEALTHY LIVE'][index]/simsDF.iloc[i]['HEALTHY LIVE'][0]
        else:
            cancer = simsDF.iloc[i]['CANCER LIVE'][index]/simsDF.iloc[i]['CANCER LIVE'][44]
            healthy = simsDF.iloc[i]['HEALTHY LIVE'][index]/simsDF.iloc[i]['HEALTHY LIVE'][44]
        if int(simsDF.iloc[i]['DOSE']) == 0:
            color = 'black'
        else:
            color = colorDict[str(simsDF.iloc[i][COLOR])]
        if MARKER != 'o':
            marker = markerDict[simsDF.iloc[i][MARKER]]
        else:
            marker = markerDict[MARKER]

        ax.scatter(cancer, healthy,  color=color, marker=marker, s=150) # facecolor='none', edgecolor=color, marker=marker, s=200, linewidth=2)

    ax.plot([1, 1], [0, ymax], color='lightgray', zorder=-1, linewidth=5)
    ax.plot([0, xmax], [1, 1], color='lightgray', zorder=-1, linewidth=5)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    if MARKER != 'o':
        for marker in markerDict:
            ax.plot([0],[0], color='black', label=marker, linestyle=markerDict[marker])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)
    ax.set_title("LIVING CANCER VS HEALTHY\nCELLS NORMALIZED", fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("LIVING CANCER CELL\nCOUNT NORMALIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("LIVING HEALTHY CELL\nCOUNT NORMALIZED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    ax.set_xlim([0,xmax])
    ax.set_ylim([0,ymax])

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_SCATTERNORM_CH_' + COLOR.replace(' ','') + '_' + MARKER.replace(' ','') + '_' + str(TIME) + '.svg', bbox_inches='tight')

    return

# ------------- ENVIRONMENT PLOTTING FUNCTIONS -----------------

def plot_env_conc_times_bar(MOL_NAME, simsDF, COLOR, FILEID, SAVELOC, TIMES):

    colorDict = COLOR_DICT[COLOR]
    HUES = {}

    SIM_TIMES = simsDF.iloc[0]['TIME']
    indicies = []
    TITLE_TIMES = ''
    for t in range(0, len(TIMES)):
        if 'VITRO' in FILEID:
            TITLE_TIMES += str(TIMES[t]) + ', '
        else:
            TITLE_TIMES += str(TIMES[t]-1) + ', '
        for s in range(0, len(SIM_TIMES)):
            if float(TIMES[t]) == float(SIM_TIMES[s]):
                indicies.append(s)
                if 'VITRO' in FILEID:
                    HUES[s] = TIMES[t]
                else:
                    HUES[s] = TIMES[t]-1

    TITLE_TIMES = TITLE_TIMES[0:-2]

    figEnv = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figEnv.add_subplot(1, 1, 1)
    YLABEL = MOL_NAME + ' CONC (' + MOL_CONC_UNITS[MOL_NAME] + ')'

    key = MOL_NAME + ' TOTAL CONC'
    envDF = make_env_df()
    order = []

    for i in range(0, len(simsDF)):
        for index in indicies:
            envDict = make_env_dict()
            y = simsDF.iloc[i][key][index]
            h = HUES[index]

            if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
                x = str(0)
            elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
                x = "CONTROL"
            else:
                x = str(simsDF.iloc[i][COLOR])

            envDict['TIME'] = h
            envDict['AXIS'] = x
            envDict['CONC'] = y

            envDF = envDF.append(envDict, ignore_index=True)

    if 'CONTROL' in envDF['AXIS'].to_list() and COLOR != 'ANTIGENS HEALTHY':
        order.append('CONTROL')

    for key in colorDict:
        if key in envDF['AXIS'].to_list():
            order.append(key)

    if 'VITRO' in FILEID:
        times = TIMES
    else:
        times = [t - 1 for t in TIMES]

    if MOL_NAME == 'GLUCOSE':
        ymax = 5e9
    elif MOL_NAME == 'IL-2':
        if 'VITRO' in FILEID:
            ymax = 4e7
        else:
            ymax = 1.4e7
    else:
        ymax = 160

    sns.barplot(ax=ax, x='TIME', y='CONC', hue='AXIS', palette=colorDict, data=envDF, order=times, hue_order=order, capsize=0.1, ci="sd")
    sns.swarmplot(ax=ax, x='TIME', y='CONC', hue='AXIS', color='black', dodge=True, order=times, hue_order=order, data=envDF)
    ax.set_xlabel(COLOR, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, labelpad=LABELPAD)
    ax.set_ylabel(YLABEL, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(MOL_NAME + "  CONCENTRATION\nOVER TIME " + TITLE_TIMES, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, pad=LABELPAD+7)
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
    ax.set_ylim(bottom=0, top=ymax)

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_CONC_' + MOL_NAME + '_XTIME_DAYS_' + TITLE_TIMES.replace(',','').replace(' ','')
                    + '.svg', bbox_inches='tight')
    return

def plot_evn_conc_times_line(MOL_NAME, simsDF, COLOR, FILEID, SAVELOC):

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            conc = simsDF.iloc[i][MOL_NAME + ' TOTAL CONC']
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
            conc = simsDF.iloc[i][MOL_NAME + ' TOTAL CONC'][2:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, conc,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, conc,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, conc,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(MOL_NAME + " CONCENTRATION\n(" + MOL_CONC_UNITS[MOL_NAME] + ")", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(MOL_NAME + " CONCENTRATION\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        ax.set_xticks(plot_time)
    else:
        ax.set_xticks(plot_time[::10])

    if MOL_NAME == 'GLUCOSE':
        ymax = 5e9
    elif MOL_NAME == 'IL-2':
        if 'VITRO' in FILEID:
            ymax = 4e7
        else:
            ymax = 1.4e7
    else:
        ymax = 160

    ax.set_ylim(bottom=0,top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_CONC_' + MOL_NAME + '.svg', bbox_inches='tight')

    return

def plot_evn_concs_ideal_realistic_parity(MOL_NAME, simsDF, XAXIS, COLOR, TIMES, FILEID, SAVELOC):

    simsDFtreated = simsDF[simsDF['DOSE'] != 0]
    simsDFideal = simsDFtreated[simsDFtreated['ANTIGENS HEALTHY'] == 0]
    simsDFrealistic = simsDFtreated[simsDFtreated['ANTIGENS HEALTHY'] == 100]

    colorDict = COLOR_DICT[COLOR]

    figParity = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figParity.add_subplot(1, 1, 1)

    if MOL_NAME == 'GLUCOSE':
        ymax = 5e9
    elif MOL_NAME == 'IL-2':
        if 'VITRO' in FILEID:
            ymax = 7e7
        else:
            ymax = 1.4e7
    else:
        ymax = 160

    timeIndicies = []
    SIM_TIMES = simsDF.iloc[0]['TIME']
    for t in TIMES:
        for s in range(0,len(SIM_TIMES)):
            if float(t) == float(SIM_TIMES[s]):
                timeIndicies.append(s)

    timeMarkers = {'1': 's',
                   '4': 'd',
                   '7': 'o'}

    for time in timeIndicies:
        ideal = []
        realistic = []
        color = []
        marker = []
        for i in range(0,len(simsDFideal)):
            ideal.append(simsDFideal.iloc[i][MOL_NAME + ' TOTAL CONC'][time])
            simRealisticMatch = simsDFrealistic[simsDFrealistic['DOSE'] == simsDFideal.iloc[i]['DOSE']]
            simRealisticMatch = simRealisticMatch[simRealisticMatch['TREAT RATIO'] == simsDFideal.iloc[i]['TREAT RATIO']]
            simRealisticMatch = simRealisticMatch[simRealisticMatch['CAR AFFINITY'] == simsDFideal.iloc[i]['CAR AFFINITY']]
            simRealisticMatch = simRealisticMatch[simRealisticMatch['ANTIGENS CANCER'] == simsDFideal.iloc[i]['ANTIGENS CANCER']]
            simRealisticMatch = simRealisticMatch[simRealisticMatch['SEED'] == simsDFideal.iloc[i]['SEED']]
            realistic.append(simRealisticMatch.iloc[0][MOL_NAME + ' TOTAL CONC'][time])
            color.append(colorDict[str(simsDFideal.iloc[i][COLOR])])
            marker.append(timeMarkers[str(int(time/2))])

        if XAXIS == 'REALISTIC':
            YAXIS = 'IDEAL'
            ax.plot([0, ymax], [0, ymax], color='lightgray', zorder=0)
            ax.set_xlim([0, ymax])
            ax.set_ylim([0, ymax])
            for i in range(0, len(ideal)):
                ax.scatter(realistic[i], ideal[i], s=100, zorder=1, color=color[i], marker=marker[i])
        else:
            YAXIS = 'REALISTIC'
            ax.plot([0, ymax], [0, ymax], color='lightgray', zorder=0)
            ax.set_xlim([0, ymax])
            ax.set_ylim([0, ymax])
            for i in range(0, len(ideal)):
                ax.scatter(ideal[i], realistic[i], s=100, zorder=1, color=color[i], marker=marker[i])

    ax.set_title(MOL_NAME + " CONCENTRATION (" + MOL_CONC_UNITS[MOL_NAME] + ")\n" + "IN " + XAXIS + " vs " + YAXIS + " CO-CULTURE", fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel(MOL_NAME + " CONCENTRATION IN " + XAXIS, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(MOL_NAME + " CONCENTRATION IN " + YAXIS, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_CONC_' + MOL_NAME + '_COCULTURE_PARITY_' + XAXIS + '_' + YAXIS + '_' + COLOR.replace(' ', '') + '.svg', bbox_inches='tight')

    return

def plot_env_conc_axis_bar(MOL_NAME, simsDF, COLOR, FILEID, SAVELOC, TIMES):

    colorDict = COLOR_DICT[COLOR]
    HUES = {}

    SIM_TIMES = simsDF.iloc[0]['TIME']
    indicies = []
    TITLE_TIMES = ''
    for t in range(0, len(TIMES)):
        if 'VITRO' in FILEID:
            TITLE_TIMES += str(TIMES[t]) + ', '
        else:
            TITLE_TIMES += str(TIMES[t]-1) + ', '
        for s in range(0, len(SIM_TIMES)):
            if float(TIMES[t]) == float(SIM_TIMES[s]):
                indicies.append(s)
                if 'VITRO' in FILEID:
                    HUES[s] = TIMES[t]
                else:
                    HUES[s] = TIMES[t]-1

    TITLE_TIMES = TITLE_TIMES[0:-2]

    figEnv = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figEnv.add_subplot(1, 1, 1)
    YLABEL = MOL_NAME + ' CONC (' + MOL_CONC_UNITS[MOL_NAME] + ')'

    key = MOL_NAME + ' TOTAL CONC'
    envDF = make_env_df()
    order = []

    for i in range(0, len(simsDF)):
        for index in indicies:
            envDict = make_env_dict()
            y = simsDF.iloc[i][key][index]
            h = HUES[index]

            if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
                x = str(0)
            elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
                x = "CONTROL"
            else:
                x = str(simsDF.iloc[i][COLOR])

            envDict['TIME'] = h
            envDict['AXIS'] = x
            envDict['CONC'] = y

            envDF = envDF.append(envDict, ignore_index=True)

    if 'CONTROL' in envDF['AXIS'].to_list() and COLOR != 'ANTIGENS HEALTHY':
        order.append('CONTROL')

    for key in colorDict:
        if key in envDF['AXIS'].to_list():
            order.append(key)

    if 'VITRO' in FILEID:
        times = TIMES
    else:
        times = [t - 1 for t in TIMES]

    if MOL_NAME == 'GLUCOSE':
        ymax = 5e9
    elif MOL_NAME == 'IL-2':
        if 'VITRO' in FILEID:
            ymax = 4e7
        else:
            ymax = 1.4e7
    else:
        ymax = 160

    sns.barplot(ax=ax, x='AXIS', y='CONC', hue='TIME', palette='Greys', data=envDF, order=order, hue_order=times, capsize=0.1, ci="sd")
    sns.swarmplot(ax=ax, x='AXIS', y='CONC', hue='TIME', color='black', dodge=True, order=order, hue_order=times, data=envDF)
    ax.set_xlabel(COLOR, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(YLABEL, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(MOL_NAME + "  CONCENTRATION\nOVER TIME " + TITLE_TIMES, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD+7)
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
    ax.set_ylim(bottom=0, top=ymax)

    plt.xticks(fontsize=TICKSIZE,  rotation = 45)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_CONC_' + MOL_NAME + '_XX_DAYS_' + TITLE_TIMES.replace(',','').replace(' ','')
                    + '.svg', bbox_inches='tight')
    return

# def plot_evn_conc_times_trend(MOL_NAME, simsDF, COLOR, FILEID, SAVELOC, TIMES):
#
#     colorDict = COLOR_DICT[COLOR]
#
#     figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
#     ax = figCounts.add_subplot(1, 1, 1)
#     for i in range(0, len(simsDF)):
#
#         if 'VITRO' in FILEID:
#             plot_time = simsDF.iloc[i]['TIME']
#             conc = simsDF.iloc[i][MOL_NAME + ' TOTAL CONC']
#         else:
#             plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
#             conc = simsDF.iloc[i][MOL_NAME + ' TOTAL CONC'][2:]
#
#         if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
#             ax.plot(plot_time, conc,
#                     linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
#                     color=colorDict[str(0)])
#         elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
#             ax.plot(plot_time, conc,
#                     linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
#                     color=colorDict["CONTROL"])
#         else:
#             ax.plot(plot_time, conc,
#                     linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
#                     color=colorDict[str(simsDF.iloc[i][COLOR])])
#     ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
#     ax.set_ylabel(MOL_NAME + " CONCENTRATION (" + MOL_CONC_UNITS[MOL_NAME] + ")", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
#     ax.set_title(MOL_NAME + " CONCENTRATION OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
#     ax.set_xlim([min(plot_time), max(plot_time)])
#
#     if 'VITRO' in FILEID:
#         ax.set_xticks(plot_time)
#     else:
#         ax.set_xticks(plot_time[::10])
#
#     if MOL_NAME == 'GLUCOSE':
#         ymax = 5e9
#     elif MOL_NAME == 'IL-2':
#         if 'VITRO' in FILEID:
#             ymax = 4e7
#         else:
#             ymax = 1.4e7
#     else:
#         ymax = 160
#
#     ax.set_ylim(bottom=0,top=ymax)
#
#     for color in colorDict:
#         ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
#     for dose in doseLineDict:
#         ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
#     ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)
#
#     plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
#     plt.yticks(fontsize=TICKSIZE)
#
#     if SAVELOC == '':
#         plt.show()
#     else:
#         plt.savefig(SAVELOC + FILEID + '_CONC_' + MOL_NAME + '.svg', bbox_inches='tight')
#
#     return

# ------------- SPATIAL PLOTTING FUNCTIONS -----------------

def plot_counts_radius(POP_NAME, simsDF, COLOR, FILEID, SAVELOC, TIME):

    colorDict = COLOR_DICT[COLOR]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figCountsRad = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCountsRad.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("RADIUS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + "\nCELL COUNTS (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nACROSS RADIUS AT TIME " + str(TIME), fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, pad=LABELPAD)
    ax.set_xlim([min(simsDF.iloc[0]['RADIUS']), max(simsDF.iloc[0]['RADIUS'])])

    if 'NORMALIZED' in POP_NAME:
        if 'CANCER' in POP_NAME or 'HEALTHY' in POP_NAME:
            ymax = 3
        else:
            ymax = 30
        ax.set_ylim(bottom=0, top=ymax)
    else:
        ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_' + POP_NAME.replace(' ','').replace('NORMALIZED','_NORM') + '_DAY_' + str(TIME) +  '.svg', bbox_inches='tight')

    return

def plot_counts_radius_dose(POP_NAME, simsDF, COLOR, FILEID, SAVELOC, TIME):

    colorDict = COLOR_DICT[COLOR]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figCountsRad = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCountsRad.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("RADIUS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + "\nCELL COUNTS (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nACROSS RADIUS AT TIME " + str(TIME), fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, pad=LABELPAD)
    ax.set_xlim([min(simsDF.iloc[0]['RADIUS']), max(simsDF.iloc[0]['RADIUS'])])

    if 'NORMALIZED' in POP_NAME:
        if 'CANCER' in POP_NAME or 'HEALTHY' in POP_NAME:
            ymax = 3
        else:
            ymax = 30
        ax.set_ylim(bottom=0, top=ymax)
    else:
        ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_' + POP_NAME.replace(' ','').replace('NORMALIZED','_NORM') + '_DAY_' + str(TIME) +  '_DOSE.svg', bbox_inches='tight')

    return

def plot_counts_radius_merge(POP_NAME, simsDF, COLOR, FILEID, SAVELOC, TIME):

    colorDict = COLOR_DICT[COLOR]

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0, len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figCountsRad = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCountsRad.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(0)])
            if 'NORMALIZED' in POP_NAME:
                pop = POP_NAME.replace(' NORMALIZED','')
                ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][pop + ' LIVE NORMALIZED'][index],
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict[str(0)])
            else:
                ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME + ' LIVE'][index],
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict["CONTROL"])
            if 'NORMALIZED' in POP_NAME:
                pop = POP_NAME.replace(' NORMALIZED','')
                ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][pop + ' LIVE NORMALIZED'][index],
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict["CONTROL"])
            else:
                ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME + ' LIVE'][index],
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict["CONTROL"])
        else:
            ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME][index],
                    linestyle=liveLineDict['TOTAL'],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
            if 'NORMALIZED' in POP_NAME:
                pop = POP_NAME.replace(' NORMALIZED','')
                ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][pop + ' LIVE NORMALIZED'][index],
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict[str(simsDF.iloc[i][COLOR])])
            else:
                ax.plot(simsDF.iloc[i]['RADIUS'], simsDF.iloc[i][POP_NAME + ' LIVE'][index],
                        linestyle=liveLineDict['LIVE'],
                        color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("RADIUS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + "\nCELL COUNTS (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT\nACROSS RADIUS AT TIME " + str(TIME), fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(simsDF.iloc[0]['RADIUS']), max(simsDF.iloc[0]['RADIUS'])])

    if 'NORMALIZED' in POP_NAME:
        if 'CANCER' in POP_NAME or 'HEALTHY' in POP_NAME:
            ymax = 3
        else:
            ymax = 30
        ax.set_ylim(bottom=0, top=ymax)
    else:
        ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for type in liveLineDict:
        ax.plot([0],[0], color='black', label=type, linestyle=liveLineDict[type])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTS_MERGE_' + POP_NAME.replace(' ','').replace('NORMALIZED','_NORM') + '_DAY_' + str(TIME) + '.svg', bbox_inches='tight')

    return

# ------------- LYSIS PLOTTING FUNCTIONS -----------------

def plot_counts_lysed_time(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            lysed = simsDF.iloc[i][POP_NAME + ' LYSED TOTAL']
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
            lysed = simsDF.iloc[i][POP_NAME + ' LYSED TOTAL'][2:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, lysed,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, lysed,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, lysed,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    if 'HEALTHY' in POP_NAME or 'DISH' in FILEID:
        xInit = [simsDF.iloc[0]['TIME'][0], simsDF.iloc[0]['TIME'][-1]]
        yInit = [simsDF.iloc[i][POP_NAME + ' SEEDED'], simsDF.iloc[i][POP_NAME + ' SEEDED']]
        ax.plot(xInit, yInit, linestyle='solid', label='SEEDED', color='gray')
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS LYSED (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])
    ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0], [0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0], [0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSLYSED_' + POP_NAME.replace(' ', '') + '.svg', bbox_inches='tight')

    return

def plot_counts_lysed_time_merge(POP_NAMES, simsDF, COLOR, FILEID, SAVELOC):
    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        for POP_NAME in POP_NAMES:
            if POP_NAME != 'TISSUE':
                if 'VITRO' in FILEID:
                    plot_time = simsDF.iloc[i]['TIME']
                    lysed = simsDF.iloc[i][POP_NAME + ' LYSED TOTAL']
                else:
                    plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
                    lysed = simsDF.iloc[i][POP_NAME + ' LYSED TOTAL'][2:]
                if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
                    ax.plot(plot_time, lysed,
                            linestyle=lysedLineDict[POP_NAME],
                            color=colorDict[str(0)])
                elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
                    ax.plot(plot_time, lysed,
                            linestyle=lysedLineDict[POP_NAME],
                            color=colorDict["CONTROL"])
                else:
                    ax.plot(plot_time, lysed,
                            linestyle=lysedLineDict[POP_NAME],
                            color=colorDict[str(simsDF.iloc[i][COLOR])])
    if 'HEALTHY' in POP_NAME or 'DISH' in FILEID:
        xInit = [plot_time[0], plot_time[-1]]
        yInit = [simsDF.iloc[i][POP_NAME + ' SEEDED'], simsDF.iloc[i][POP_NAME + ' SEEDED']]
        ax.plot(xInit, yInit, linestyle='solid', label='SEEDED', color='gray')
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("TISSUE CELL COUNTS LYSED (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title("CELL COUNT OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])
    ax.set_ylim(bottom=0)

    if 'VITRO' in FILEID:
        ax.set_xticks(plot_time)
    else:
        ax.set_xticks(plot_time[::10])

    for color in colorDict:
        ax.plot([0], [0], color=colorDict[color], label=color, linestyle='solid')
    for type in lysedLineDict:
        ax.plot([0], [0], color='black', label=type, linestyle=lysedLineDict[type])
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)

    if 'VITRO' in FILEID:
        plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    else:
        plt.xticks(plot_time[::10], [int(i) for i in plot_time[::10]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSLYSED_MERGE' + '.svg', bbox_inches='tight')

    return

def plot_counts_lysed_time_exact(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(simsDF.iloc[i]['TIME EXACT'], simsDF.iloc[i][POP_NAME + ' LYSED CUMULATIVE EXACT'],
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(simsDF.iloc[i]['TIME EXACT'], simsDF.iloc[i][POP_NAME + ' LYSED CUMULATIVE EXACT'],
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(simsDF.iloc[i]['TIME EXACT'], simsDF.iloc[i][POP_NAME + ' LYSED CUMULATIVE EXACT'],
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    xInit = [simsDF.iloc[0]['TIME'][0], simsDF.iloc[0]['TIME'][-1]*1440]
    yInit = [simsDF.iloc[i][POP_NAME + ' SEEDED'], simsDF.iloc[i][POP_NAME + ' SEEDED']]
    ax.plot(xInit, yInit, linestyle='solid', label='SEEDED', color='gray')
    ax.set_xlabel("TIME (MINUTES)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " CELL COUNTS LYSED (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " COUNT OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(simsDF.iloc[0]['TIME']), max(simsDF.iloc[0]['TIME'])*1440])
    ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0], [0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0], [0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE, rotation=90)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSLYSEDEXACT_' + POP_NAME.replace(' ', '') + '.svg', bbox_inches='tight')

    return

def plot_counts_lysed_time_exact_merge(POP_NAMES, simsDF, COLOR, FILEID, SAVELOC):
    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):
        for POP_NAME in POP_NAMES:
            if POP_NAME != 'TISSUE':
                    if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
                        ax.plot(simsDF.iloc[i]['TIME EXACT'], simsDF.iloc[i][POP_NAME + ' LYSED CUMULATIVE EXACT'],
                                linestyle=lysedLineDict[POP_NAME],
                                color=colorDict[str(0)])
                    elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
                        ax.plot(simsDF.iloc[i]['TIME EXACT'], simsDF.iloc[i][POP_NAME + ' LYSED CUMULATIVE EXACT'],
                                linestyle=lysedLineDict[POP_NAME],
                                color=colorDict["CONTROL"])
                    else:
                        ax.plot(simsDF.iloc[i]['TIME EXACT'], simsDF.iloc[i][POP_NAME + ' LYSED CUMULATIVE EXACT'],
                                linestyle=lysedLineDict[POP_NAME],
                                color=colorDict[str(simsDF.iloc[i][COLOR])])
    xInit = [simsDF.iloc[0]['TIME'][0], simsDF.iloc[0]['TIME'][-1]*1440]
    yInit = [simsDF.iloc[i][POP_NAME + ' SEEDED'], simsDF.iloc[i][POP_NAME + ' SEEDED']]
    ax.plot(xInit, yInit, linestyle='solid', label='SEEDED', color='gray')
    ax.set_xlabel("TIME (MINUTES)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("TISSUE CELL COUNTS LYSED (NUMBERS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title("CELL COUNT OVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(simsDF.iloc[0]['TIME']), max(simsDF.iloc[0]['TIME'])*1440])
    ax.set_ylim(bottom=0)

    for color in colorDict:
        ax.plot([0], [0], color=colorDict[color], label=color, linestyle='solid')
    for type in lysedLineDict:
        ax.plot([0], [0], color='black', label=type, linestyle=lysedLineDict[type])
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE, rotation=90)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSLYSEDEXACT_MERGE' + '.svg', bbox_inches='tight')

    return

def plot_counts_lysed_radius_exact(POP_NAME, simsDF, COLOR, FILEID, SAVELOC):
    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        radius = []
        time = []

        for c in range(0, len(simsDF.iloc[i]['TISSUE LYSED EXACT'])):
            if POP_NAME != 'TISSUE':
                if POP_CODES[POP_NAME] == simsDF.iloc[i]['TISSUE LYSED EXACT'][c]:
                    radius.append(simsDF.iloc[i]['RADIUS LYSED EXACT'][c])
                    time.append(simsDF.iloc[i]['TIME EXACT'][c])
            else:
                radius.append(simsDF.iloc[i]['RADIUS LYSED EXACT'][c])
                time.append(simsDF.iloc[i]['TIME EXACT'][c])
        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(time, radius,
                    # marker=doseMarkerDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(time, radius,
                    # marker=doseMarkerDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(time, radius,
                    # marker=doseMarkerDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])

    xInit = [simsDF.iloc[0]['TIME'][0], simsDF.iloc[0]['TIME'][-1]*1440]
    yInit = [simsDF.iloc[i][POP_NAME + ' SEEDED'], simsDF.iloc[i][POP_NAME + ' SEEDED']]
    ax.plot(xInit, yInit, linestyle='solid', label='SEEDED', color='gray')
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(POP_NAME + " RADIUS LYSED", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(POP_NAME + " KILLED OVER TIME AND RADIUS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(simsDF.iloc[0]['TIME']), max(simsDF.iloc[0]['TIME'])*1440])
    ax.set_ylim(bottom=0, top=max(simsDF.iloc[0]['RADIUS']))

    for color in colorDict:
        ax.plot([0], [0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0], [0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_COUNTSLYSEDRADIUS_' + POP_NAME.replace(' ', '') + '.svg', bbox_inches='tight')

    return

# ------------- RANK PLOTTING FUNCTION -----------------

def plot_rank_parity(rankDF, COLOR, XRANK, SAVELOC):

    figRank = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figRank.add_subplot(1, 1, 1)
    ax.plot([0, len(rankDF)+1], [0, len(rankDF)+1], color='lightgray', zorder=0)

    colorDict = COLOR_DICT[COLOR]

    if XRANK == 'VITRO':
        YRANK = 'VIVO'
        rankDFsorted = rankDF.sort_values(by='RANK DISH', ascending=True)
        x = rankDFsorted['RANK DISH']
        y = rankDFsorted['RANK TISSUE']
        for i in range(0,len(rankDFsorted)):
            ax.scatter(x.iloc[i], y.iloc[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])
    else:
        YRANK = 'VITRO'
        rankDFsorted = rankDF.sort_values(by='RANK TISSUE', ascending=True)
        x = rankDFsorted['RANK TISSUE']
        y = rankDFsorted['RANK DISH']
        for i in range(0, len(rankDFsorted)):
            ax.scatter(x.iloc[i], y.iloc[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])


    ax.set_title("RANK IN " + XRANK + " vs RANK IN " + YRANK, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("RANK IN " + XRANK, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("RANK IN" + YRANK, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    ax.set_xlim([0,len(rankDF)+1])
    ax.set_ylim([0,len(rankDF)+1])

    plt.xticks([i for i in range(1,15)], fontsize=TICKSIZE-10, rotation=45)
    plt.yticks([i for i in range(1,15)], fontsize=TICKSIZE-10)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'RANK_' + XRANK + '_' + YRANK + '_' + COLOR.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_rank_ladder(rankDF, COLOR, XRANK, SAVELOC):

    figRank = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figRank.add_subplot(1, 1, 1)

    colorDict = COLOR_DICT[COLOR]

    if XRANK == 'VITRO':
        YRANK = 'VIVO'
        rankDFsorted = rankDF.sort_values(by='RANK DISH', ascending=False)
        x = rankDFsorted['RANK DISH']
        y = rankDFsorted['RANK TISSUE']
        for i in range(0,len(rankDFsorted)):
            ax.plot([1, 2], [x.iloc[i], y.iloc[i]], marker='o', zorder=1, markerfacecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])], markeredgecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])], color=colorDict[str(rankDFsorted.iloc[i][COLOR])])
    else:
        YRANK = 'VITRO'
        rankDFsorted = rankDF.sort_values(by='RANK TISSUE', ascending=False)
        x = rankDFsorted['RANK TISSUE']
        y = rankDFsorted['RANK DISH']
        for i in range(0, len(rankDFsorted)):
            ax.plot([1, 2], [x.iloc[i], y.iloc[i]], marker='o', zorder=1, markerfacecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])], markeredgecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])], color=colorDict[str(rankDFsorted.iloc[i][COLOR])])


    ax.set_title("RANK IN " + XRANK + " vs RANK IN " + YRANK, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("CONTEXT", fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("RANK", fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    ax.set_xlim([0.7,2.3])
    ax.set_ylim([0,len(rankDF)+1])

    plt.xticks([1, 2], labels=[XRANK, YRANK], fontsize=TICKSIZE-10)
    #plt.xlabel([XRANK, YRANK], fontsize=TICKSIZE-10, rotation=45)
    plt.yticks([i for i in range(1,15)], fontsize=TICKSIZE-10)
    plt.gca().invert_yaxis()

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'RANK_LADDER_' + XRANK + '_' + YRANK + '_' + COLOR.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_score_parity(rankDF, COLOR, XRANK, SAVELOC):

    figRank = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figRank.add_subplot(1, 1, 1)

    colorDict = COLOR_DICT[COLOR]

    if XRANK == 'VITRO':
        YRANK = 'VIVO'
        rankDFsorted = rankDF.sort_values(by='SCORE DISH', ascending=True)
        x = rankDFsorted['SCORE DISH']
        y = rankDFsorted['SCORE TISSUE']
        ax.plot([-0.3, 0.8], [-0.3, 0.8], color='lightgray', zorder=0)
        ax.set_xlim([-0.3, 0.8])
        ax.set_ylim([-0.3, 0.8])
        for i in range(0,len(rankDFsorted)):
            ax.scatter(x.iloc[i], y.iloc[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])
    else:
        YRANK = 'VITRO'
        rankDFsorted = rankDF.sort_values(by='SCORE TISSUE', ascending=True)
        x = rankDFsorted['SCORE TISSUE']
        y = rankDFsorted['SCORE DISH']
        ax.plot([-0.3, 0.8], [-0.3, 0.8], color='lightgray', zorder=0)
        ax.set_xlim([-0.3, 0.8])
        ax.set_ylim([-0.3, 0.8])
        for i in range(0, len(rankDFsorted)):
            ax.scatter(x.iloc[i], y.iloc[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])


    ax.set_title("SCORE IN " + XRANK + " vs SCORE IN " + YRANK, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("SCORE IN " + XRANK, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("SCORE IN" + YRANK, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    plt.xticks(fontsize=TICKSIZE-10)
    plt.yticks(fontsize=TICKSIZE-10)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'SCORE_' + XRANK + '_' + YRANK + '_' + COLOR.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_score_norm_parity(rankDF, COLOR, XRANK, SAVELOC):

    figRank = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figRank.add_subplot(1, 1, 1)
    ax.plot([0, 1], [0, 1], color='lightgray', zorder=0)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])

    colorDict = COLOR_DICT[COLOR]

    if XRANK == 'VITRO':
        YRANK = 'VIVO'
        rankDFsorted = rankDF.sort_values(by='SCORE DISH', ascending=True)
        a = rankDFsorted['SCORE DISH']
        b = rankDFsorted['SCORE TISSUE']
        x = [(i - min(a)) / (max(a) - min(a)) for i in a]
        y = [(i - min(b)) / (max(b) - min(b)) for i in b]
        for i in range(0,len(x)):
            ax.scatter(x[i], y[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])
    else:
        YRANK = 'VITRO'
        rankDFsorted = rankDF.sort_values(by='SCORE TISSUE', ascending=True)
        a = rankDFsorted['SCORE TISSUE']
        b = rankDFsorted['SCORE DISH']

        x = [(i - min(a)) / (max(a) - min(a)) for i in a]
        y = [(i - min(b)) / (max(b) - min(b)) for i in b]
        for i in range(0, len(x)):
            ax.scatter(x[i], y[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])


    ax.set_title("SCORE IN " + XRANK + " vs SCORE IN " + YRANK, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("SCORE IN " + XRANK, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("SCORE IN" + YRANK, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    plt.xticks(fontsize=TICKSIZE-10, rotation=45)
    plt.yticks(fontsize=TICKSIZE-10)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'SCORENORM_' + XRANK + '_' + YRANK + '_' + COLOR.replace(' ','') + '.svg', bbox_inches='tight')

    return

def plot_score_ladder(rankDF, COLOR, XRANK, SAVELOC):

    figRank = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figRank.add_subplot(1, 1, 1)

    colorDict = COLOR_DICT[COLOR]

    if XRANK == 'VITRO':
        YRANK = 'VIVO'
        rankDFsorted = rankDF.sort_values(by='SCORE DISH', ascending=True)
        x = rankDFsorted['SCORE DISH']
        y = rankDFsorted['SCORE TISSUE']
        for i in range(0,len(rankDFsorted)):
            ax.plot([1, 2], [x.iloc[i], y.iloc[i]], marker='o', zorder=1,
                    markerfacecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])],
                    markeredgecolor=colorDict[str(rankDFsorted.iloc[i][COLOR])],
                    color=colorDict[str(rankDFsorted.iloc[i][COLOR])])
    else:
        YRANK = 'VITRO'
        rankDFsorted = rankDF.sort_values(by='SCORE TISSUE', ascending=True)
        x = rankDFsorted['SCORE TISSUE']
        y = rankDFsorted['SCORE DISH']
        for i in range(0, len(rankDFsorted)):
            ax.scatter(x.iloc[i], y.iloc[i], s=100, zorder=1, color=colorDict[str(rankDFsorted.iloc[i][COLOR])])


    ax.set_title("SCORE IN " + XRANK + " vs SCORE IN " + YRANK, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("CONTEXT", fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("SCORE", fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    ax.set_xlim([0.7, 2.3])
    ax.set_ylim([-0.3, 0.8])
    plt.yticks([-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    plt.xticks([1, 2], labels=[XRANK, YRANK], fontsize=TICKSIZE - 10)
    plt.xticks(fontsize=TICKSIZE-10)
    plt.yticks(fontsize=TICKSIZE-10)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'SCORE_LADDER_' + XRANK + '_' + YRANK + '_' + COLOR.replace(' ','') + '.svg', bbox_inches='tight')

    return

# -------- BINDING HEURISTIC PLOTTING FUNCTION ---------

def plot_binding_heuristic_CAR(SAVELOC):

    Vloc = 6.7e-12
    Na = 6.022e23
    CARS = 50000
    CARSavg = 50000

    alpha = 3
    beta = 0.01
    gamma = 0.2

    KDs = [1e-6, 1e-7, 1e-8, 1e-9]

    figBinding = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figBinding.add_subplot(1, 1, 1)

    ligands = [i for i in range(0,10000)]

    for KD in KDs:
        y = []
        for L in ligands:
            h = ((gamma*L)/((beta*KD*Vloc*Na) + (gamma*L)))*(CARS/CARSavg)*alpha
            p = (2*(1/(1 + exp(-1*h)))) - 1
            y.append(p)
        ax.plot(ligands, y, color=affinityColorDict[str(KD)], linewidth=4, label=str(KD))

    ax.set_xlabel("ANTIGENS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("PROBABILITY OF\nBINDING AND KILLING", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title("BINDING HEURISTIC\nFOR CAR-ANTIGEN BINDING", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(ligands), max(ligands)])
    ax.set_ylim([0,1])
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
    plt.xticks(fontsize=TICKSIZE, rotation=45)
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'BINDING_HEURISTIC_CAR' + '.svg', bbox_inches='tight')

    return

def plot_binding_heuristic_self(SAVELOC):

    Vloc = 6.7e-12
    Na = 6.022e23
    KDself = 7.8e-6
    PD1_START = 150

    alpha = 3
    beta = 0.02
    gamma = 0.2

    PD1s = [150, 1000, 3000, 9000]
    PD1_colors = { "150": 'lightgray',
                     "1000": 'darkgray',
                     "3000": 'gray',
                     "9000": 'black'
    }

    figBinding = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figBinding.add_subplot(1, 1, 1)

    ligands = [i for i in range(0, 50000)]

    for PD1 in PD1s:
        y = []
        for L in ligands:
            h = ((gamma * L) / ((beta * KDself * Vloc * Na) + (gamma * L))) * (PD1 / PD1_START) * alpha
            p = (2 * (1 / (1 + exp(-1 * h)))) - 1
            y.append(p)
        ax.plot(ligands, y, color=PD1_colors[str(PD1)], linewidth=4, label='PD1 = ' + str(PD1))

    ax.set_xlabel("SELF LIGANDS (PDL1s)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("PROBABILITY OF\nBINDING", fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title("BINDING HEURISTIC\nFOR PD1-PDL1 BINDING", fontname='Arial', fontweight='bold',
                 fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(ligands), max(ligands)])
    ax.set_ylim([0, 1])
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
    plt.xticks(fontsize=TICKSIZE, rotation=45)
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + 'BINDING_HEURISTIC_SELF' + '.svg', bbox_inches='tight')

    return

# ------------- MAIN PLOTTING FUNCTION -----------------

def plot_data(simsDF, ANALYSIS, COLOR, MARKER, PARTIAL, FILEID, SAVELOC):

    filesplit = FILEID.split('_')

    if ANALYSIS == 'ANALYZED':

        # Count Xs in filesplit
        X = 0
        for i in range(5, 10):
            if filesplit[i] == 'X': X += 1

        if X == 1 or PARTIAL:
            if COLOR == 'X':
                if X == 1:
                    if filesplit[FILEID_SPLIT_INDICES['DOSE']] == 'X': COLOR = 'DOSE'
                    if filesplit[FILEID_SPLIT_INDICES['TREAT RATIO']] == 'X': COLOR = 'TREAT RATIO'
                    if filesplit[FILEID_SPLIT_INDICES['CAR AFFINITY']] == 'X': COLOR = 'CAR AFFINITY'
                    if filesplit[FILEID_SPLIT_INDICES['ANTIGENS CANCER']] == 'X': COLOR = 'ANTIGENS CANCER'
                    if filesplit[FILEID_SPLIT_INDICES['ANTIGENS HEALTHY']] == 'X': COLOR = 'ANTIGENS HEALTHY'

                if PARTIAL:
                    COLOR = 'CAR AFFINITY'

            fileid = FILEID.replace('_STATES','')

            # Plot counts data
            if 'CYCLES' not in FILEID and 'VOLUMES' not in FILEID:

                print('\t\t' + 'Plotting count data over time')
                for p in range(0, len(POP_NAMES)):
                    if filesplit[FILEID_SPLIT_INDICES['ANTIGENS HEALTHY']] == 'NA' and 'HEALTHY' in POP_NAMES[p]:
                        continue
                    else:
                        plot_counts(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)
                        plot_counts_norm(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)

                        if 'VIVO' in FILEID and POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                            plot_counts_treat(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)

                        if COLOR == 'DOSE' or PARTIAL:
                            plot_counts_dose(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)
                            plot_counts_norm_dose(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)

                            if 'VIVO' in FILEID and POP_NAMES[p] in ['CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE', 'T-CELL', 'T-CELL LIVE']:
                                plot_counts_treat_dose(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)

                        if 'LIVE' not in POP_NAMES[p]:
                            plot_counts_merge(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)
                            plot_counts_norm_merge(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)

                            if 'VIVO' in FILEID and POP_NAMES[p] in ['CD4', 'CD8', 'T-CELL']:
                                plot_counts_treat_merge(POP_NAMES[p], simsDF, COLOR, fileid, SAVELOC)

                        plt.close("all")
                
                # Plot count fracs data
                for p in range(0, len(POP_FRAC_NAMES)):
                    if filesplit[FILEID_SPLIT_INDICES['ANTIGENS HEALTHY']] == 'NA' and 'HEALTHY' in POP_FRAC_NAMES[p]:
                        continue
                    else:
                        print('\t\t' + 'Plotting count killed fraction data')
                        plot_counts_frac_remaining(POP_FRAC_NAMES[p].replace(' LIVE', ''), simsDF, COLOR, fileid, SAVELOC)
                        plot_counts_frac_remaining(POP_FRAC_NAMES[p], simsDF, COLOR, fileid, SAVELOC)
                        if COLOR == 'DOSE' or PARTIAL:
                            plot_counts_frac_remaining_dose(POP_FRAC_NAMES[p], simsDF, COLOR, fileid, SAVELOC)
                            plot_counts_frac_remaining_dose(POP_FRAC_NAMES[p].replace(' LIVE', ''), simsDF, COLOR, fileid, SAVELOC)
                        plot_counts_frac_remaining_merge(POP_FRAC_NAMES[p].replace(' LIVE', ''), simsDF, COLOR, fileid, SAVELOC)
                        POP_NAME = POP_FRAC_NAMES[p] + ' %'
                        print('\t\t' + 'Plotting count fraction data')
                        plot_count_fracs(POP_NAME, simsDF, COLOR, fileid, SAVELOC)

                plt.close("all")

                # Plot scatter
                print('\t\t' + 'Plotting scatter data')
                if filesplit[FILEID_SPLIT_INDICES['POPS']] == 'CH':
                    MARKER='o'
                    if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'DISH':
                        plot_CH_scatter(simsDF, COLOR, MARKER, fileid, SAVELOC, 7)
                        plot_CH_scatter_normalized(simsDF, COLOR, MARKER, fileid, SAVELOC, 7)
                    if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'TISSUE':
                        plot_CH_scatter(simsDF, COLOR, MARKER, fileid, SAVELOC, 29)
                        plot_CH_scatter_normalized(simsDF, COLOR, MARKER, fileid, SAVELOC, 29)
                        plot_CH_scatter(simsDF, COLOR, MARKER, fileid, SAVELOC, 31)
                        plot_CH_scatter_normalized(simsDF, COLOR, MARKER, fileid, SAVELOC, 31)
                    plt.close("all")

                # Plot type fracs
                print('\t\t' + 'Plotting type fraction data')
                plot_state_fracs(simsDF, COLOR, fileid, SAVELOC)
                plot_state_fracs_neutral(simsDF, COLOR, fileid, SAVELOC)
                plt.close("all")

            # Plot volume distributions
            if 'STATES' not in FILEID:
                if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'DISH' and 'CYCLES' not in FILEID:
                    print('\t\t' + 'Plotting volume distribution data')
                    plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 1)
                    plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 4)
                    plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 7)
                    plot_volumes_split(simsDF, COLOR, FILEID, SAVELOC, [4,7])
                    plt.close("all")
                if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'TISSUE' and 'CYCLES' not in FILEID:
                    print('\t\t' + 'Plotting volume distribution data')
                    plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 1)
                    plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 5)
                    plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 22)
                    plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 26)
                    plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 29)
                    plot_volumes(simsDF, COLOR, FILEID, SAVELOC, 31)
                    plt.close("all")
                    
            # Plot cycle distributions
            if 'STATES' not in FILEID:
                if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'DISH' and 'VOLUMES' not in FILEID:
                    print('\t\t' + 'Plotting cycle distribution data')
                    # Plot in units of minutes
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 1, False)
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 4, False)
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 7, False)
                    plt.close("all")
                    # Plot in units of hours
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 1, True)
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 4, True)
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 7, True)
                    plot_cycles_split(simsDF, COLOR, FILEID, SAVELOC, [4,7], True)
                    plt.close("all")
                if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'TISSUE' and 'VOLUMES' not in FILEID:
                    print('\t\t' + 'Plotting cycle distribution data')
                    # Plot in units of minutes
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 22, False)
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 26, False)
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 29, False)
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 31, False)
                    plt.close("all")
                    # Plot in units of hours
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 22, True)
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 26, True)
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 29, True)
                    plot_cycles(simsDF, COLOR, FILEID, SAVELOC, 31, True)
                    plt.close("all")
        
        if 'VOLUMES' not in FILEID and 'CYCLES' not in FILEID and not PARTIAL:
            fileid = FILEID.replace('_STATES', '')
            if filesplit[FILEID_SPLIT_INDICES['CAR AFFINITY']] == 'X' and filesplit[FILEID_SPLIT_INDICES['ANTIGENS CANCER']] == 'X':
                if COLOR == 'X': COLOR = 'CAR AFFINITY'

                # Plot kill curve
                print('\t\t' + 'Plotting kill curve data')
                if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'DISH':
                    plot_kill_curve_sim(simsDF, fileid, SAVELOC, 4)
                    plot_kill_curve_sim(simsDF, fileid, SAVELOC, 7)
                    plot_kill_curve_relative_sim(simsDF, fileid, SAVELOC, 4)
                    plot_kill_curve_relative_sim(simsDF, fileid, SAVELOC, 7)
                    plot_kill_curve_normalized_sim(simsDF, FILEID, SAVELOC, 7)
                if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'TISSUE':
                    plot_kill_curve_sim(simsDF, fileid, SAVELOC, 26)
                    plot_kill_curve_sim(simsDF, fileid, SAVELOC, 29)
                    plot_kill_curve_sim(simsDF, fileid, SAVELOC, 31)
                    plot_kill_curve_relative_sim(simsDF, fileid, SAVELOC, 26)
                    plot_kill_curve_relative_sim(simsDF, fileid, SAVELOC, 29)
                    plot_kill_curve_relative_sim(simsDF, fileid, SAVELOC, 31)
                plot_kill_curve_exp(SAVELOC)
                plot_kill_curve_normalized_exp(SAVELOC)
                plot_kill_curve_exp_separated(SAVELOC)
                plot_kill_curve_normalized_exp_separated(SAVELOC)
                plt.close("all")
            
            if filesplit[FILEID_SPLIT_INDICES['DOSE']] == 'X' and filesplit[FILEID_SPLIT_INDICES['CAR AFFINITY']] == 'X':

                # Plot ET Ratio
                print('\t\t' + 'Plotting E:T ratio data')
                plot_ET_ratio(simsDF, fileid, SAVELOC)
                plot_ET_ratio_relative(simsDF, fileid, SAVELOC)
                plt.close("all")

            if 'CH' in filesplit[FILEID_SPLIT_INDICES['DOSE']]:
                if X >= 2:
                    if COLOR == 'X': COLOR = 'CAR AFFINITY'
                    if MARKER == 'X': MARKER = 'TREAT RATIO'

                    # Plot counts data
                    print('\t\t' + 'Plotting cancer-healthy data')
                    if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'DISH':
                        plot_CH_scatter(simsDF, COLOR, MARKER, fileid, SAVELOC, 7)
                    if filesplit[FILEID_SPLIT_INDICES['PLATE']] == 'TISSUE':
                        plot_CH_scatter(simsDF, COLOR, MARKER, fileid, SAVELOC, 29)
                        plot_CH_scatter(simsDF, COLOR, MARKER, fileid, SAVELOC, 31)
                    plt.close("all")

    if ANALYSIS == 'ENVIRONMENT':

        # Count Xs in filesplit
        X = 0
        for i in range(5, 10):
            if filesplit[i] == 'X': X += 1

        if X == 1 or PARTIAL:
            if COLOR == 'X':
                if X == 1:
                    if filesplit[FILEID_SPLIT_INDICES['DOSE']] == 'X': COLOR = 'DOSE'
                    if filesplit[FILEID_SPLIT_INDICES['TREAT RATIO']] == 'X': COLOR = 'TREAT RATIO'
                    if filesplit[FILEID_SPLIT_INDICES['CAR AFFINITY']] == 'X': COLOR = 'CAR AFFINITY'
                    if filesplit[FILEID_SPLIT_INDICES['ANTIGENS CANCER']] == 'X': COLOR = 'ANTIGENS CANCER'
                    if filesplit[FILEID_SPLIT_INDICES['ANTIGENS HEALTHY']] == 'X': COLOR = 'ANTIGENS HEALTHY'

                if PARTIAL:
                    COLOR = 'CAR AFFINITY'

            if filesplit[FILEID_SPLIT_INDICES['EXP']] == 'VITRO':
                TIMES = [i for i in DISH_TIMES if i != 0]
            else:
                TIMES = TISSUE_TIMES

            for m in range(0, len(MOL_NAMES)):
                if MOL_NAMES[m] == 'OXYGEN':
                    continue
                else:
                    plot_env_conc_times_bar(MOL_NAMES[m], simsDF, COLOR, FILEID, SAVELOC, TIMES)
                    plot_evn_conc_times_line(MOL_NAMES[m], simsDF, COLOR, FILEID, SAVELOC)
                    plot_env_conc_axis_bar(MOL_NAMES[m], simsDF, COLOR, FILEID, SAVELOC, TIMES)
                    plt.close("all")

        if X == 5 and 'VITRO' in FILEID:
            for m in range(0, len(MOL_NAMES)):
                if MOL_NAMES[m] == 'OXYGEN':
                    continue
                else:
                    for XAXIS in ['IDEAL', 'REALISTIC']:
                        if COLOR == 'X':
                            for color in ['DOSE', 'TREAT RATIO', 'CAR AFFINITY', 'ANTIGENS CANCER']:
                                plot_evn_concs_ideal_realistic_parity(MOL_NAMES[m], simsDF, XAXIS, color, [7], FILEID, SAVELOC)
                                plt.close("all")
                        else:
                            plot_evn_concs_ideal_realistic_parity(MOL_NAMES[m], simsDF, XAXIS, COLOR, [7], FILEID, SAVELOC)
                            plt.close("all")

    if ANALYSIS == 'SPATIAL':

        if filesplit[FILEID_SPLIT_INDICES['EXP']] == 'VITRO':
            TIMES = DISH_TIMES
        else:
            TIMES = TISSUE_TIMES

        # Count Xs in filesplit
        X = 0
        for i in range(5, 10):
            if filesplit[i] == 'X': X += 1

        if X == 1 or PARTIAL:
            if COLOR == 'X':
                if X == 1:
                    if filesplit[FILEID_SPLIT_INDICES['DOSE']] == 'X': COLOR = 'DOSE'
                    if filesplit[FILEID_SPLIT_INDICES['TREAT RATIO']] == 'X': COLOR = 'TREAT RATIO'
                    if filesplit[FILEID_SPLIT_INDICES['CAR AFFINITY']] == 'X': COLOR = 'CAR AFFINITY'
                    if filesplit[FILEID_SPLIT_INDICES['ANTIGENS CANCER']] == 'X': COLOR = 'ANTIGENS CANCER'
                    if filesplit[FILEID_SPLIT_INDICES['ANTIGENS HEALTHY']] == 'X': COLOR = 'ANTIGENS HEALTHY'

                if PARTIAL:
                    COLOR = 'CAR AFFINITY'

            # Plot counts data
            print('\t\t' + 'Plotting count data across radius')
            for p in range(0, len(POP_NAMES)):
                for TIME in TIMES:
                    if filesplit[FILEID_SPLIT_INDICES['ANTIGENS HEALTHY']] == 'NA' and 'HEALTHY' in POP_NAMES[p]:
                        continue
                    else:
                        plot_counts_radius(POP_NAMES[p], simsDF, COLOR, FILEID, SAVELOC, TIME)
                        plot_counts_radius(POP_NAMES[p] + ' NORMALIZED', simsDF, COLOR, FILEID, SAVELOC, TIME)
                        if COLOR == 'DOSE' or PARTIAL:
                            plot_counts_radius_dose(POP_NAMES[p], simsDF, COLOR, FILEID, SAVELOC, TIME)
                            plot_counts_radius_dose(POP_NAMES[p] + ' NORMALIZED', simsDF, COLOR, FILEID, SAVELOC, TIME)

                        plt.close("all")

                    if 'LIVE' not in POP_NAMES[p]:
                        plot_counts_radius_merge(POP_NAMES[p], simsDF, COLOR, FILEID, SAVELOC, TIME)
                        plot_counts_radius_merge(POP_NAMES[p] + ' NORMALIZED', simsDF, COLOR, FILEID, SAVELOC, TIME)

                        plt.close("all")

    if ANALYSIS == 'LYSED':

        # Count Xs in filesplit
        X = 0
        for i in range(5, 10):
            if filesplit[i] == 'X': X += 1

        if X == 1 or PARTIAL:
            if COLOR == 'X':
                if X == 1:
                    if filesplit[FILEID_SPLIT_INDICES['DOSE']] == 'X': COLOR = 'DOSE'
                    if filesplit[FILEID_SPLIT_INDICES['TREAT RATIO']] == 'X': COLOR = 'TREAT RATIO'
                    if filesplit[FILEID_SPLIT_INDICES['CAR AFFINITY']] == 'X': COLOR = 'CAR AFFINITY'
                    if filesplit[FILEID_SPLIT_INDICES['ANTIGENS CANCER']] == 'X': COLOR = 'ANTIGENS CANCER'
                    if filesplit[FILEID_SPLIT_INDICES['ANTIGENS HEALTHY']] == 'X': COLOR = 'ANTIGENS HEALTHY'

                if PARTIAL:
                    COLOR = 'CAR AFFINITY'

        # Plot counts data
        print('\t\t' + 'Plotting lysis count data over time')
        for p in range(0, len(POP_LYSIS_NAMES)):
            if filesplit[FILEID_SPLIT_INDICES['ANTIGENS HEALTHY']] == 'NA' and 'HEALTHY' in POP_LYSIS_NAMES[p]:
                continue
            else:
                plot_counts_lysed_time(POP_LYSIS_NAMES[p], simsDF, COLOR, FILEID, SAVELOC)
                plot_counts_lysed_time_exact(POP_LYSIS_NAMES[p], simsDF, COLOR, FILEID, SAVELOC)
                # plot_counts_lysed_radius_exact(POP_LYSIS_NAMES[p], simsDF, COLOR, FILEID, SAVELOC)
            plt.close("all")

        if filesplit[FILEID_SPLIT_INDICES['POPS']] == 'CH':
            plot_counts_lysed_time_merge(POP_LYSIS_NAMES, simsDF, COLOR, FILEID, SAVELOC)
            plot_counts_lysed_time_exact_merge(POP_LYSIS_NAMES, simsDF, COLOR, FILEID, SAVELOC)
            plt.close("all")

    return

def color_test():
    x = [0, 10]
    y1 = [0, 10]
    y2 = [0, 9]
    y3 = [0, 8]
    y4 = [0, 7]
    y5 = [0, 6]

    y = [y1, y2, y3, y4, y5]

    TESTCOLORS = []
    cmap = mplcm.get_cmap('Purples')
    values = [0.5, 1.0]
    for v in values:
        TESTCOLORS.append(cmap(v))

    figTest = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figTest.add_subplot(1, 1, 1)
    for i in range(0,2):
        ax.plot(x, y[i], color=TESTCOLORS[i], label=i, linewidth=5)

    plt.show()
    return

# ---------------- MAIN FUNCTION ---------------------

if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()
    analysis = ''

    print("Making figures for the following files:")

    if args.rank:

        fileName = re.sub('.*RANK', 'RANK', args.files)
        print('\t' + fileName)

        with open(args.files, 'rb') as f:
            rankDF = pd.read_excel(f, header=0)

        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', -1)

        rankDF['TREAT RATIO'] = rankDF['TREAT RATIO'].astype(str)
        rankDF['CAR AFFINITY'] = rankDF['CAR AFFINITY'].astype(str)
        for i in range(0,len(rankDF)):
            rankDF['TREAT RATIO'].iloc[i] = TREAT_RATIO_DICT_REVERSE[rankDF['TREAT RATIO'].iloc[i]]
        #
        # print(rankDF)
        # print(rankDF.dtypes)

        if args.color == 'X':
            for COLOR in ['DOSE', 'TREAT RATIO', 'CAR AFFINITY', 'ANTIGENS CANCER']:
                plot_rank_parity(rankDF, COLOR, 'VITRO', args.saveLoc)
                plot_rank_parity(rankDF, COLOR, 'VIVO', args.saveLoc)
                plot_rank_ladder(rankDF, COLOR, 'VITRO', args.saveLoc)
                plot_rank_ladder(rankDF, COLOR, 'VIVO', args.saveLoc)
                plot_score_parity(rankDF, COLOR, 'VITRO', args.saveLoc)
                plot_score_parity(rankDF, COLOR, 'VIVO', args.saveLoc)
                plot_score_norm_parity(rankDF, COLOR, 'VITRO', args.saveLoc)
                plot_score_norm_parity(rankDF, COLOR, 'VIVO', args.saveLoc)
                plot_score_ladder(rankDF, COLOR, 'VITRO', args.saveLoc)
                plot_score_ladder(rankDF, COLOR, 'VIVO', args.saveLoc)
                plt.close("all")

    else:

        # Get files
        PKLFILES = get_pkl_files(args.files)

        for file in PKLFILES:

            if 'VITRO' in file:
                fileName = re.sub('.*VITRO', 'VITRO', file)

            else:
                fileName = re.sub('.*VIVO', 'VIVO', file)

            print('\t' + fileName)

            FILEID = ''

            if 'ANALYZED' in fileName:
                FILEID = fileName.replace('_ANALYZED.pkl', '')
                analysis = 'ANALYZED'

            if 'ENVIRONMENT' in fileName:
                FILEID = fileName.replace('_ENVIRONMENT.pkl', '')
                analysis = 'ENVIRONMENT'

            if 'SPATIAL' in fileName:
                FILEID = fileName.replace('_SPATIAL.pkl', '')
                analysis = 'SPATIAL'

            if 'LYSED' in fileName:
                FILEID = fileName.replace('_LYSED.pkl', '')
                analysis = 'LYSED'

            with open(file, 'rb') as f:
                simsDF = pickle.load(f)

            pd.set_option('display.max_rows', None)
            pd.set_option('display.max_columns', None)
            pd.set_option('display.width', None)
            pd.set_option('display.max_colwidth', -1)

            #print(simsDF.iloc[0])

            #plot_binding_heuristic_CAR(args.saveLoc)
            #plot_binding_heuristic_self(args.saveLoc)

            plot_data(simsDF, analysis, args.color, args.marker, args.partial, FILEID, args.saveLoc)

