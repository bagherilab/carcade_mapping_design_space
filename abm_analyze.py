import ABM
from abm_parse import load as ABM_load
from abm_parse import get_hex_rings
from abm_parse import get_radius
import os
import re
import pickle
import pandas as pd
import numpy as np
from argparse import ArgumentParser

__author__ = "Alexis N. Prybutok"
__email__ = "aprybutok@u.northwestern.edu"

'''
ABM_ANALYZE takes a directory of (or a single) .pkl simulation files
and extracts the data into a dataframe. If the analyze data is desired, 
the resulting file will contain a data frame in the form:

    TUMOR ID | SEED | PLATE | DAMAGE | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

and saves it to a .pkl where the DATA are the following list of information:
    
    TIME
    CANCER
    CANCER LIVE
    HEALTHY
    HEALTHY LIVE
    T-CELL
    T-CELL LIVE
    CD4
    CD4 LIVE
    CD8
    CD8 LIVE
    CANCER LIVE %
    HEALTHY LIVE %
    NEUTR CANCER
    NEUTR CANCER LIVE
    NEUTR HEALTHY
    NEUTR HEALTHY LIVE
    NEUTR T-CELL
    NEUTR CD4
    NEUTR CD8
    NEUTR CANCER %
    NEUTR CANCER LIVE %
    NEUTR HEALTHY %
    NEUTR HEALTHY LIVE %
    NEUTR T-CELL %
    NEUTR CD4 %
    NEUTR CD8 %
    APOPT CANCER
    APOPT HEALTHY
    APOPT T-CELL
    APOPT CD4
    APOPT CD8
    APOPT CANCER %
    APOPT HEALTHY %
    APOPT T-CELL %
    APOPT CD4 %
    APOPT CD8 %
    QUIES CANCER
    QUIES CANCER LIVE
    QUIES HEALTHY
    QUIES HEALTHY LIVE
    QUIES CANCER %
    QUIES CANCER LIVE %
    QUIES HEALTHY %
    QUIES HEALTHY LIVE %
    MIGRA CANCER
    MIGRA CANCER LIVE
    MIGRA HEALTHY
    MIGRA HEALTHY LIVE
    MIGRA T-CELL
    MIGRA CD4
    MIGRA CD8
    MIGRA CANCER %
    MIGRA CANCER LIVE %
    MIGRA HEALTHY %
    MIGRA HEALTHY LIVE %
    MIGRA T-CELL %
    MIGRA CD4 %
    MIGRA CD8 %
    PROLI CANCER
    PROLI CANCER LIVE
    PROLI HEALTHY
    PROLI HEALTHY LIVE
    PROLI T-CELL
    PROLI CD4
    PROLI CD8
    PROLI CANCER %
    PROLI CANCER LIVE %
    PROLI HEALTHY LIVE %
    PROLI HEALTHY %
    PROLI T-CELL %
    PROLI CD4 %
    PROLI CD8 %
    SENES CANCER
    SENES CANCER LIVE
    SENES HEALTHY
    SENES HEALTHY LIVE
    SENES T-CELL
    SENES CD4
    SENES CD8
    SENES CANCER %
    SENES CANCER LIVE %
    SENES HEALTHY %
    SENES HEALTHY LIVE %
    SENES T-CELL %
    SENES CD4 %
    SENES CD8 %
    NECRO CANCER
    NECRO HEALTHY
    NECRO CANCER %
    NECRO HEALTHY %
    CYTOT T-CELL
    CYTOT CD8
    CYTOT T-CELL %
    CYTOT CD8 %
    STIMU T-CELL
    STIMU CD4
    STIMU T-CELL %
    STIMU CD4 %
    EXHAU T-CELL
    EXHAU CD4
    EXHAU CD8
    EXHAU T-CELL %
    EXHAU CD4 %
    EXHAU CD8 %
    ANERG T-CELL
    ANERG CD4
    ANERG CD8
    ANERG T-CELL %
    ANERG CD4 %
    ANERG CD8 %
    STARV T-CELL %
    STARV CD4 %
    STARV CD8 %
    PAUSE T-CELL %
    PAUSE CD4 %
    PUASE CD8%
    AVG CELL CYCLES CANCER
    AVG CELL CYCLES HEALTHY
    AVG CELL CYCLES T-CELL
    AVG CELL CYCLES CD4
    AVG CELL CYCLES CD8
    CELL VOLUMES CANCER
    CELL VOLUMES HEALTHY
    CELL VOLUMES T-CELL
    CELL VOLUMES CD4
    CELL VOLUMES CD8

each in the format of a list of the value of the specified information at each point in time.

If environment data is desired, and environment flag used, will produce a corresponding environment file that contains
extracts the data into a dataframe in the form:

    TUMOR ID | SEED | PLATE | DAMAGE | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

and saves it to a .pkl where the DATA are the following list of information:

    TIME
    RADIUS
    SEED
    PLATE
    NUTRIENTS
    DOSE
    TREAT RATIO
    CAR AFFINITY
    ANTIGENS CANCER
    ANTIGENS HEALTHY
    GLUCOSE
    OXYGEN
    TGFA
    IL-2
    GLUCOSE TOTAL
    OXYGEN TOTAL
    TGFA TOTAL
    IL-2 TOTAL
    GLUCOSE TOTAL CONC
    TGFA TOTAL CONC
    IL-2 TOTAL CONC
    
where GLUCOSE, OXYGEN, TGFA, and IL-2 are in the format of a list (per time point) of a list of concentrations at each radius.
TOTAL columns are the total (absolute) amount of molecule in total area at each time point.
TOTAL CONC columns are the total concentration amount of molecule across the entire simulation at each time point.

If spatial data is desired, and spatial flag used, it will produce a corresponding spatial file that contains
extracts the data into a dataframe in the form:

    TUMOR ID | SEED | PLATE | DAMAGE | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

and saves it to a .pkl where the DATA are the following list of information:

    TIME
    RADIUS
    SEED
    PLATE
    NUTRIENTS
    DOSE
    TREAT RATIO
    CAR AFFINITY
    ANTIGENS CANCER
    ANTIGENS HEALTHY
    CANCER
    CANCER LIVE
    CANCER NORMALIZED
    CANCER LIVE NORMALIZED
    HEALTHY
    HEALTHY LIVE
    HEALTHY NORMALIZED
    HEALTHY LIVE NORMALIZED
    T-CELL
    T-CELL LIVE
    T-CELL NORMALIZED
    T-CELL LIVE NORMALIZED
    CD4
    CD4 LIVE
    CD4 NORMALIZED
    CD4 LIVE NORMALIZED
    CD8
    CD8 LIVE
    CD8 NORMALIZED
    CD8 LIVE NORMALIZED
    
where each cell population is in the format of a list of a list of counts at each radius.

Usage:
    python abm_analyze.py FILES [--saveLoc SAVELOC] [--analyze] [--environment] [--spatial]

    FILES
        Path to .pkl or directory
    [--analyze]
        Flag indicating cell counts and state data will be analyzed (_ANALYZE)
    [--environment]
        Flag indicating environment data will be saved to another file (_ENVIRONMENT)
    [--spatial]
        Flag indicating spatial data will be saved to another file (_SPATIAL)
    [--saveLoc SAVELOC]
        Location of where to save file, default won't save files
'''

# -------------- PARSING AND INPUT/OUTPUT FUNCTIONS -------------

def get_parser():

    # Setup argument parser.
    parser = ArgumentParser(description="Record all ABM CART simulation data into a dataframe")
    parser.add_argument(dest="files", help="Path to .pkl file or directory")
    parser.add_argument("--analyze", default=False, dest="analyze",
                        action='store_true', help="Save counts/states data to separate file (_ANALYZED)")
    parser.add_argument("--environment", default=False, dest="environment",
                        action='store_true', help="Save environment data to separate file (_ENVIRONMENT")
    parser.add_argument("--spatial", default=False, dest="spatial",
                        action='store_true', help="Save spatial data to separate file (_SPATIAL)")
    parser.add_argument("--saveLoc", default='', dest="saveLoc",
                        help="Location of where to save file, default won't save file")

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

def get_json_files(arg):
    if arg[-1] == "/" or arg[-1] == "\\":
        return [arg + f for f in os.listdir(arg) if ABM.is_json(f)]
    else:
        assert ABM.is_json(arg)
        return [arg]



# ---------------- MAKE EMPTY CONTAINERS ---------------------

def make_analyze_df():
    columns = ['TUMOR ID',
               'SEED',
               'PLATE',
               'NUTRIENTS',
               'DOSE',
               'TREAT RATIO',
               'CAR AFFINITY',
               'ANTIGENS CANCER',
               'ANTIGENS HEALTHY',
               'TIME',
               'CANCER',
               'CANCER LIVE',
               'HEALTHY',
               'HEALTHY LIVE',
               'T-CELL',
               'T-CELL LIVE',
               'CD4',
               'CD4 LIVE',
               'CD8',
               'CD8 LIVE',
               'CANCER LIVE %',
               'HEALTHY LIVE %',
               'NEUTR CANCER',
               'NEUTR CANCER LIVE',
               'NEUTR HEALTHY',
               'NEUTR HEALTHY LIVE',
               'NEUTR T-CELL',
               'NEUTR CD4',
               'NEUTR CD8',
               'NEUTR CANCER %',
               'NEUTR CANCER LIVE %',
               'NEUTR HEALTHY %',
               'NEUTR HEALTHY LIVE %',
               'NEUTR T-CELL %',
               'NEUTR CD4 %',
               'NEUTR CD8 %',
               'APOPT CANCER',
               'APOPT HEALTHY',
               'APOPT T-CELL',
               'APOPT CD4',
               'APOPT CD8',
               'APOPT CANCER %',
               'APOPT HEALTHY %',
               'APOPT T-CELL %',
               'APOPT CD4 %',
               'APOPT CD8 %',
               'QUIES CANCER',
               'QUIES CANCER LIVE',
               'QUIES HEALTHY',
               'QUIES HEALTHY LIVE',
               'QUIES CANCER %',
               'QUIES CANCER LIVE %',
               'QUIES HEALTHY %',
               'QUIES HEALTHY LIVE %',
               'MIGRA CANCER',
               'MIGRA CANCER LIVE',
               'MIGRA HEALTHY',
               'MIGRA HEALTHY LIVE',
               'MIGRA T-CELL',
               'MIGRA CD4',
               'MIGRA CD8',
               'MIGRA CANCER %',
               'MIGRA CANCER LIVE %',
               'MIGRA HEALTHY %',
               'MIGRA HEALTHY LIVE %',
               'MIGRA T-CELL %',
               'MIGRA CD4 %',
               'MIGRA CD8 %',
               'PROLI CANCER',
               'PROLI CANCER LIVE',
               'PROLI HEALTHY',
               'PROLI HEALTHY LIVE',
               'PROLI T-CELL',
               'PROLI CD4',
               'PROLI CD8',
               'PROLI CANCER %',
               'PROLI CANCER LIVE %',
               'PROLI HEALTHY LIVE %',
               'PROLI HEALTHY %',
               'PROLI T-CELL %',
               'PROLI CD4 %',
               'PROLI CD8 %',
               'SENES CANCER',
               'SENES CANCER LIVE',
               'SENES HEALTHY',
               'SENES HEALTHY LIVE',
               'SENES T-CELL',
               'SENES CD4',
               'SENES CD8',
               'SENES CANCER %',
               'SENES CANCER LIVE %',
               'SENES HEALTHY %',
               'SENES HEALTHY LIVE %',
               'SENES T-CELL %',
               'SENES CD4 %',
               'SENES CD8 %',
               'NECRO CANCER',
               'NECRO HEALTHY',
               'NECRO CANCER %',
               'NECRO HEALTHY %',
               'CYTOT T-CELL',
               'CYTOT CD8',
               'CYTOT T-CELL %',
               'CYTOT CD8 %',
               'STIMU T-CELL',
               'STIMU CD4',
               'STIMU T-CELL %',
               'STIMU CD4 %',
               'EXHAU T-CELL',
               'EXHAU CD4',
               'EXHAU CD8',
               'EXHAU T-CELL %',
               'EXHAU CD4 %',
               'EXHAU CD8 %',
               'ANERG T-CELL',
               'ANERG CD4',
               'ANERG CD8',
               'ANERG T-CELL %',
               'ANERG CD4 %',
               'ANERG CD8 %',
               'STARV T-CELL',
               'STARV CD4',
               'STARV CD8',
               'STARV T-CELL %',
               'STARV CD4 %',
               'STARV CD8 %',
               'PAUSE T-CELL',
               'PAUSE CD4',
               'PAUSE CD8',
               'PAUSE T-CELL %',
               'PAUSE CD4 %',
               'PAUSE CD8 %',
               'AVG CELL CYCLES CANCER',
               'AVG CELL CYCLES HEALTHY',
               'AVG CELL CYCLES T-CELL',
               'AVG CELL CYCLES CD4',
               'AVG CELL CYCLES CD8',
               'CELL VOLUMES CANCER',
               'CELL VOLUMES HEALTHY',
               'CELL VOLUMES T-CELL',
               'CELL VOLUMES CD4',
               'CELL VOLUMES CD8',
    ]
    simsDF = pd.DataFrame(columns=columns)

    return simsDF

def make_env_df():
    columns = ['TUMOR ID',
               'SEED',
               'PLATE',
               'NUTRIENTS',
               'DOSE',
               'TREAT RATIO',
               'CAR AFFINITY',
               'ANTIGENS CANCER',
               'ANTIGENS HEALTHY',
               'TIME',
               'RADIUS',
               'GLUCOSE',
               'OXYGEN',
               'TGFA',
               'IL-2',
               'GLUCOSE TOTAL',
               'OXYGEN TOTAL',
               'TGFA TOTAL',
               'IL-2 TOTAL',
               'GLUCOSE TOTAL CONC',
               'TGFA TOTAL CONC',
               'IL-2 TOTAL CONC'
    ]
    envDF = pd.DataFrame(columns=columns)

    return envDF

def make_spatial_df():
    columns = ['TUMOR ID',
               'SEED',
               'PLATE',
               'NUTRIENTS',
               'DOSE',
               'TREAT RATIO',
               'CAR AFFINITY',
               'ANTIGENS CANCER',
               'ANTIGENS HEALTHY',
               'TIME',
               'RADIUS',
               'CANCER',
               'CANCER LIVE',
               'CANCER NORMALIZED',
               'CANCER LIVE NORMALIZED',
               'HEALTHY',
               'HEALTHY LIVE',
               'HEALTHY NORMALIZED',
               'HEALTHY LIVE NORMALIZED',
               'T-CELL',
               'T-CELL LIVE',
               'T-CELL NORMALIZED',
               'T-CELL LIVE NORMALIZED',
               'CD4',
               'CD4 LIVE',
               'CD4 NORMALIZED',
               'CD4 LIVE NORMALIZED',
               'CD8',
               'CD8 LIVE',
               'CD8 NORMALIZED',
               'CD8 LIVE NORMALIZED'
    ]
    spatialsDF = pd.DataFrame(columns=columns)

    return spatialsDF

def make_analyze_dict():
    simsDict = {    'TUMOR ID': None,
                    'SEED': None,
                    'PLATE': None,
                    'NUTRIENTS': None,
                    'DOSE': None,
                    'TREAT RATIO': None,
                    'CAR AFFINITY': None,
                    'ANTIGENS CANCER': None,
                    'ANTIGENS HEALTHY': None,
                    'TIME': [],
                    'CANCER': [],
                    'CANCER LIVE': [],
                    'HEALTHY': [],
                    'HEALTHY LIVE': [],
                    'T-CELL': [],
                    'T-CELL LIVE': [],
                    'CD4': [],
                    'CD4 LIVE': [],
                    'CD8': [],
                    'CD8 LIVE': [],
                    'CANCER LIVE %': [],
                    'HEALTHY LIVE %': [],
                    'NEUTR CANCER': [],
                    'NEUTR CANCER LIVE': [],
                    'NEUTR HEALTHY': [],
                    'NEUTR HEALTHY LIVE': [],
                    'NEUTR T-CELL': [],
                    'NEUTR CD4': [],
                    'NEUTR CD8': [],
                    'NEUTR CANCER %': [],
                    'NEUTR CANCER LIVE %': [],
                    'NEUTR HEALTHY %': [],
                    'NEUTR HEALTHY LIVE %': [],
                    'NEUTR T-CELL %': [],
                    'NEUTR CD4 %': [],
                    'NEUTR CD8 %': [],
                    'APOPT CANCER': [],
                    'APOPT HEALTHY': [],
                    'APOPT T-CELL': [],
                    'APOPT CD4': [],
                    'APOPT CD8': [],
                    'APOPT CANCER %': [],
                    'APOPT HEALTHY %': [],
                    'APOPT T-CELL %': [],
                    'APOPT CD4 %': [],
                    'APOPT CD8 %': [],
                    'QUIES CANCER': [],
                    'QUIES CANCER LIVE': [],
                    'QUIES HEALTHY': [],
                    'QUIES HEALTHY LIVE': [],
                    'QUIES CANCER %': [],
                    'QUIES CANCER LIVE %': [],
                    'QUIES HEALTHY %': [],
                    'QUIES HEALTHY LIVE %': [],
                    'MIGRA CANCER': [],
                    'MIGRA CANCER LIVE': [],
                    'MIGRA HEALTHY': [],
                    'MIGRA HEALTHY LIVE': [],
                    'MIGRA T-CELL': [],
                    'MIGRA CD4': [],
                    'MIGRA CD8': [],
                    'MIGRA CANCER %': [],
                    'MIGRA CANCER LIVE %': [],
                    'MIGRA HEALTHY %': [],
                    'MIGRA HEALTHY LIVE %': [],
                    'MIGRA T-CELL %': [],
                    'MIGRA CD4 %': [],
                    'MIGRA CD8 %': [],
                    'PROLI CANCER': [],
                    'PROLI CANCER LIVE': [],
                    'PROLI HEALTHY': [],
                    'PROLI HEALTHY LIVE': [],
                    'PROLI T-CELL': [],
                    'PROLI CD4': [],
                    'PROLI CD8': [],
                    'PROLI CANCER %': [],
                    'PROLI CANCER LIVE %': [],
                    'PROLI HEALTHY %': [],
                    'PROLI HEALTHY LIVE %': [],
                    'PROLI T-CELL %': [],
                    'PROLI CD4 %': [],
                    'PROLI CD8 %': [],
                    'SENES CANCER': [],
                    'SENES CANCER LIVE': [],
                    'SENES HEALTHY': [],
                    'SENES HEALTHY LIVE': [],
                    'SENES T-CELL': [],
                    'SENES CD4': [],
                    'SENES CD8': [],
                    'SENES CANCER %': [],
                    'SENES CANCER LIVE %': [],
                    'SENES HEALTHY %': [],
                    'SENES HEALTHY LIVE %': [],
                    'SENES T-CELL %': [],
                    'SENES CD4 %': [],
                    'SENES CD8 %': [],
                    'NECRO CANCER': [],
                    'NECRO HEALTHY': [],
                    'NECRO CANCER %': [],
                    'NECRO HEALTHY %': [],
                    'CYTOT T-CELL': [],
                    'CYTOT CD8': [],
                    'CYTOT T-CELL %': [],
                    'CYTOT CD8 %': [],
                    'STIMU T-CELL': [],
                    'STIMU CD4': [],
                    'STIMU T-CELL %': [],
                    'STIMU CD4 %': [],
                    'EXHAU T-CELL': [],
                    'EXHAU CD4': [],
                    'EXHAU CD8': [],
                    'EXHAU T-CELL %': [],
                    'EXHAU CD4 %': [],
                    'EXHAU CD8 %': [],
                    'ANERG T-CELL': [],
                    'ANERG CD4': [],
                    'ANERG CD8': [],
                    'ANERG T-CELL %': [],
                    'ANERG CD4 %': [],
                    'ANERG CD8 %': [],
                    'STARV T-CELL': [],
                    'STARV CD4': [],
                    'STARV CD8': [],
                    'STARV T-CELL %': [],
                    'STARV CD4 %': [],
                    'STARV CD8 %': [],
                    'PAUSE T-CELL': [],
                    'PAUSE CD4': [],
                    'PAUSE CD8': [],
                    'PAUSE T-CELL %': [],
                    'PAUSE CD4 %': [],
                    'PAUSE CD8 %': [],
                    'AVG CELL CYCLES CANCER': [],
                    'AVG CELL CYCLES HEALTHY': [],
                    'AVG CELL CYCLES T-CELL': [],
                    'AVG CELL CYCLES CD4': [],
                    'AVG CELL CYCLES CD8': [],
                    'CELL VOLUMES CANCER': [],
                    'CELL VOLUMES HEALTHY': [],
                    'CELL VOLUMES T-CELL': [],
                    'CELL VOLUMES CD4': [],
                    'CELL VOLUMES CD8': [],
    }

    return simsDict

def make_env_dict():
    envDict = { 'TUMOR ID': None,
                'SEED': None,
                'PLATE': None,
                'NUTRIENTS': None,
                'DOSE': None,
                'TREAT RATIO': None,
                'CAR AFFINITY': None,
                'ANTIGENS CANCER': None,
                'ANTIGENS HEALTHY': None,
                'TIME': [],
                'RADIUS': [],
                'GLUCOSE': [],
                'OXYGEN': [],
                'TGFA': [],
                'IL-2': [],
                'GLUCOSE TOTAL': [],
                'OXYGEN TOTAL': [],
                'TGFA TOTAL': [],
                'IL-2 TOTAL': [],
                'GLUCOSE TOTAL CONC': [],
                'TGFA TOTAL CONC': [],
                'IL-2 TOTAL CONC': []
    }

    return envDict

def make_spatial_dict():
    spatialDict = { 'TUMOR ID': None,
                    'SEED': None,
                    'PLATE': None,
                    'NUTRIENTS': None,
                    'DOSE': None,
                    'TREAT RATIO': None,
                    'CAR AFFINITY': None,
                    'ANTIGENS CANCER': None,
                    'ANTIGENS HEALTHY': None,
                    'TIME': [],
                    'RADIUS': [],
                    'CANCER': [],
                    'CANCER LIVE': [],
                    'CANCER NORMALIZED': [],
                    'CANCER LIVE NORMALIZED': [],
                    'HEALTHY': [],
                    'HEALTHY LIVE': [],
                    'HEALTHY NORMALIZED': [],
                    'HEALTHY LIVE NORMALIZED': [],
                    'T-CELL': [],
                    'T-CELL LIVE': [],
                    'T-CELL NORMALIZED': [],
                    'T-CELL LIVE NORMALIZED': [],
                    'CD4': [],
                    'CD4 LIVE': [],
                    'CD4 NORMALIZED': [],
                    'CD4 LIVE NORMALIZED': [],
                    'CD8': [],
                    'CD8 LIVE': [],
                    'CD8 NORMALIZED': [],
                    'CD8 LIVE NORMALIZED': []
    }

    return spatialDict

# ---------------------- ANALYZE TUMOR ---------------------------

def analyze_sim(simsDF, agents, T, C, TUMORID, SEED):

    # Make simulation dict
    simDict = make_analyze_dict()

    # Add simulation information to dict
    simDict['TUMOR ID'] = TUMORID
    simDict['SEED'] = int(SEED)

    tsplit = TUMORID.split('_')

    simDict['PLATE'] = tsplit[1]
    simDict['NUTRIENTS'] = 'CONSTANT' if tsplit[0] == 'VITRO' else 'GRAPH VASCULATURE'
    simDict['DOSE'] = int(tsplit[4])
    simDict['TREAT RATIO'] = tsplit[5].replace('-',':') if tsplit[5] != "NA" else tsplit[5]
    simDict['CAR AFFINITY'] = float(tsplit[6]) if tsplit[6] != "NA" else tsplit[6]
    simDict['ANTIGENS CANCER'] = int(tsplit[7])
    simDict['ANTIGENS HEALTHY'] = int(tsplit[8]) if tsplit[8] != "NA" else tsplit[8]

    # CELL POP NUMBERS
    POP_CANCER = 0
    POP_HEALTHY = 1
    POP_CD4 = 2
    POP_CD8 = 3

    # These don't change even if pop number does
    INDEX_CANCER = 0
    INDEX_CANCERLIVE = 1
    INDEX_HEALTHY = 2
    INDEX_HEALTHYLIVE = 3
    INDEX_TCELL = 4
    INDEX_TCELLLIVE = 5
    INDEX_CD4 = 6
    INDEX_CD4LIVE = 7
    INDEX_CD8 = 8
    INDEX_CD8LIVE = 9

    POP_NAMES = ['CANCER', 'CANCER LIVE', 'HEALTHY', 'HEALTHY LIVE', 'T-CELL', 'T-CELL LIVE', 'CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']
    POP_INDICES = [INDEX_CANCER, INDEX_CANCERLIVE, INDEX_HEALTHY, INDEX_HEALTHYLIVE, INDEX_TCELL, INDEX_TCELLLIVE, INDEX_CD4, INDEX_CD4LIVE, INDEX_CD8, INDEX_CD8LIVE]
    STATES_CANCER = ['NEUTR', 'APOPT', 'QUIES', 'MIGRA', 'PROLI', 'SENES', 'NECRO']
    STATES_CANCER_INDICES = [0, 1, 2, 3, 4, 5, 6]
    STATES_TISSUELIVE = ['NEUTR', 'QUIES', 'MIGRA', 'PROLI', 'SENES']
    STATES_TISSUELIVE_INDICES = [0, 2, 3, 4, 5]
    STATES_TCELL = ['NEUTR', 'APOPT', 'MIGRA', 'PROLI', 'SENES', 'CYTOT', 'STIMU', 'EXHAU', 'ANERG', 'STARV', 'PAUSE']
    STATES_TCELL_INDICES = [0, 1, 3, 4, 5, 7, 8, 9, 10, 11, 12]


    # Set height
    H = 0

    for time in range(0, len(T)):

        # Reset counts and types
        counts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # [CANCER, CANCER LIVE, HEALTHY, HEALTHY LIVE, T-CELL, T-CELL LIVE, CD4, CD4 LIVE, CD8, CD8 LIVE]

        types = [   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # CANCER
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # CANCER LIVE
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # HEALTHY
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # HEALTHY LIVE
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # T-CELL
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # T-CELL LIVE
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # CD4
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # CD4 LIVE
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],  # CD8
                    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]   # CD8 LIVE
        ]

        typesFrac = []

        cycles = [[], [], [], [], [], [], [], [], [], []]   # [CANCER, CANCER LIVE, HEALTHY, HEALTHY LIVE, T-CELL, T-CELL LIVE, CD4, CD4 LIVE, CD8, CD8 LIVE]
        volumes = [[], [], [], [], [], [], [], [], [], []]  # [CANCER, CANCER LIVE, HEALTHY, HEALTHY LIVE, T-CELL, T-CELL LIVE, CD4, CD4 LIVE, CD8, CD8 LIVE]
        #3367
        for loc in range(0, len(C)):
            for pos in range(0, 54):
                if agents[time][H][loc]['pop'][pos] != -1:

                    # Find CANCER cells
                    if agents[time][H][loc]['pop'][pos] == POP_CANCER:
                        counts[INDEX_CANCER] += 1
                        types[INDEX_CANCER][agents[time][H][loc]['type'][pos]] += 1
                        volumes[INDEX_CANCER].append(agents[time][H][loc]['volume'][pos])
                        if agents[time][H][loc]['cycle'][pos] != -1:
                            cycles[INDEX_CANCER].append(agents[time][H][loc]['cycle'][pos])

                        # Find CANCER LIVE cells
                        if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                            counts[INDEX_CANCERLIVE] += 1
                            types[INDEX_CANCERLIVE][agents[time][H][loc]['type'][pos]] += 1
                            volumes[INDEX_CANCERLIVE].append(agents[time][H][loc]['volume'][pos])
                            if agents[time][H][loc]['cycle'][pos] != -1:
                                cycles[INDEX_CANCERLIVE].append(agents[time][H][loc]['cycle'][pos])

                    # Find HEALTHY cells
                    if agents[time][H][loc]['pop'][pos] == POP_HEALTHY:
                        counts[INDEX_HEALTHY] += 1
                        types[INDEX_HEALTHY][agents[time][H][loc]['type'][pos]] += 1
                        volumes[INDEX_HEALTHY].append(agents[time][H][loc]['volume'][pos])
                        if agents[time][H][loc]['cycle'][pos] != -1:
                            cycles[INDEX_HEALTHY].append(agents[time][H][loc]['cycle'][pos])

                        # Find HEALTHY LIVE cells
                        if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                            counts[INDEX_HEALTHYLIVE] += 1
                            types[INDEX_HEALTHYLIVE][agents[time][H][loc]['type'][pos]] += 1
                            volumes[INDEX_HEALTHYLIVE].append(agents[time][H][loc]['volume'][pos])
                            if agents[time][H][loc]['cycle'][pos] != -1:
                                cycles[INDEX_HEALTHYLIVE].append(agents[time][H][loc]['cycle'][pos])

                    # Find CD4 cells
                    if agents[time][H][loc]['pop'][pos] == POP_CD4:
                        counts[INDEX_CD4] += 1
                        counts[INDEX_TCELL] += 1
                        types[INDEX_CD4][agents[time][H][loc]['type'][pos]] += 1
                        types[INDEX_TCELL][agents[time][H][loc]['type'][pos]] += 1
                        volumes[INDEX_CD4].append(agents[time][H][loc]['volume'][pos])
                        volumes[INDEX_TCELL].append(agents[time][H][loc]['volume'][pos])
                        if agents[time][H][loc]['cycle'][pos] != -1:
                            cycles[INDEX_CD4].append(agents[time][H][loc]['cycle'][pos])
                            cycles[INDEX_TCELL].append(agents[time][H][loc]['cycle'][pos])

                        # Find CD4 LIVE cells
                        if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                            counts[INDEX_CD4LIVE] += 1
                            counts[INDEX_TCELLLIVE] += 1
                            types[INDEX_CD4LIVE][agents[time][H][loc]['type'][pos]] += 1
                            types[INDEX_TCELLLIVE][agents[time][H][loc]['type'][pos]] += 1
                            volumes[INDEX_CD4LIVE].append(agents[time][H][loc]['volume'][pos])
                            volumes[INDEX_TCELLLIVE].append(agents[time][H][loc]['volume'][pos])
                            if agents[time][H][loc]['cycle'][pos] != -1:
                                cycles[INDEX_CD4LIVE].append(agents[time][H][loc]['cycle'][pos])
                                cycles[INDEX_TCELLLIVE].append(agents[time][H][loc]['cycle'][pos])

                    # Find CD8 cells
                    if agents[time][H][loc]['pop'][pos] == POP_CD8:
                        counts[INDEX_CD8] += 1
                        counts[INDEX_TCELL] += 1
                        types[INDEX_CD8][agents[time][H][loc]['type'][pos]] += 1
                        types[INDEX_TCELL][agents[time][H][loc]['type'][pos]] += 1
                        volumes[INDEX_CD8].append(agents[time][H][loc]['volume'][pos])
                        volumes[INDEX_TCELL].append(agents[time][H][loc]['volume'][pos])
                        if agents[time][H][loc]['cycle'][pos] != -1:
                            cycles[INDEX_CD8].append(agents[time][H][loc]['cycle'][pos])
                            cycles[INDEX_TCELL].append(agents[time][H][loc]['cycle'][pos])

                        # Find CD8 LIVE cells
                        if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                            counts[INDEX_CD8LIVE] += 1
                            counts[INDEX_TCELLLIVE] += 1
                            types[INDEX_CD8LIVE][agents[time][H][loc]['type'][pos]] += 1
                            types[INDEX_TCELLLIVE][agents[time][H][loc]['type'][pos]] += 1
                            volumes[INDEX_CD8LIVE].append(agents[time][H][loc]['volume'][pos])
                            volumes[INDEX_TCELLLIVE].append(agents[time][H][loc]['volume'][pos])
                            if agents[time][H][loc]['cycle'][pos] != -1:
                                cycles[INDEX_CD8LIVE].append(agents[time][H][loc]['cycle'][pos])
                                cycles[INDEX_TCELLLIVE].append(agents[time][H][loc]['cycle'][pos])

        # Calculate state fractions for each pop
        for c in range(0, len(counts)):
            if counts[c] == 0:
                typesFrac.append([0 for x in types[c]])
            else:
                typesFrac.append([x / counts[c] for x in types[c]])

        # Add information to dictionary
        simDict['TIME'].append(T[time])

        for p in range(0, len(POP_NAMES)):
            simDict[POP_NAMES[p]].append(counts[POP_INDICES[p]])

            if 'LIVE' not in POP_NAMES[p]:

                simDict['AVG CELL CYCLES ' + POP_NAMES[p]].append(cycles[POP_INDICES[p]])
                simDict['CELL VOLUMES ' + POP_NAMES[p]].append(volumes[POP_INDICES[p]])

                if POP_NAMES[p] == 'CANCER' or POP_NAMES[p] == 'HEALTHY':
                    for s in range(0, len(STATES_CANCER)):
                        key = STATES_CANCER[s] + ' ' + POP_NAMES[p]
                        simDict[key].append(types[POP_INDICES[p]][STATES_CANCER_INDICES[s]])
                        simDict[key + ' %'].append(typesFrac[POP_INDICES[p]][STATES_CANCER_INDICES[s]])

                if POP_NAMES[p] == 'T-CELL' or POP_NAMES[p] == 'CD4' or POP_NAMES[p] == 'CD8':
                    for s in range(0, len(STATES_TCELL)):
                        if (POP_NAMES[p] == 'CD4' and STATES_TCELL[s] == 'CYTOT'):
                            continue
                        elif (POP_NAMES[p] == 'CD8' and STATES_TCELL[s] == 'STIMU'):
                            continue
                        else:
                            key = STATES_TCELL[s] + ' ' + POP_NAMES[p]
                            simDict[key].append(types[POP_INDICES[p]][STATES_TCELL_INDICES[s]])
                            simDict[key + ' %'].append(typesFrac[POP_INDICES[p]][STATES_TCELL_INDICES[s]])

            else:
                if POP_NAMES[p] == 'CANCER LIVE' or POP_NAMES[p] == 'HEALTHY LIVE':
                    if counts[POP_INDICES[p-1]] == 0:
                        simDict[POP_NAMES[p] + ' %'].append(0)
                    else:
                        simDict[POP_NAMES[p] + ' %'].append(counts[POP_INDICES[p]] / counts[POP_INDICES[p-1]])

                    for s in range(0, len(STATES_TISSUELIVE)):
                        key = STATES_TISSUELIVE[s] + ' ' + POP_NAMES[p]
                        simDict[key].append(types[POP_INDICES[p]][STATES_TISSUELIVE_INDICES[s]])
                        simDict[key + ' %'].append(typesFrac[POP_INDICES[p]][STATES_TISSUELIVE_INDICES[s]])

    # Add tumor information to full simulation dataframe
    simsDF = simsDF.append(simDict, ignore_index=True)

    return simsDF

def env_sim(envDF, environments, T, R, TUMORID, SEEDS):
    HEIGHT = 0
    HEX_VOL_UM = 6780.97    # um^3

    for SEED in range(0, SEEDS):
        print("\t\t  >" + TUMORID + "_0" + str(SEED) + " ENVIRONMENT")

        # Make simulation dict
        envDict = make_env_dict()

        # Add simulation information to dict
        envDict['TUMOR ID'] = TUMORID

        tsplit = TUMORID.split('_')

        envDict['PLATE'] = tsplit[1]
        envDict['NUTRIENTS'] = 'CONSTANT' if tsplit[0] == 'VITRO' else 'GRAPH VASCULATURE'
        envDict['DOSE'] = int(tsplit[4])
        envDict['TREAT RATIO'] = tsplit[5].replace('-', ':') if tsplit[5] != "NA" else tsplit[5]
        envDict['CAR AFFINITY'] = float(tsplit[6]) if tsplit[6] != "NA" else tsplit[6]
        envDict['ANTIGENS CANCER'] = int(tsplit[7])
        envDict['ANTIGENS HEALTHY'] = int(tsplit[8]) if tsplit[8] != "NA" else tsplit[8]
        envDict['SEED'] = int(SEED)

        hexRings = get_hex_rings(R)
        radii = [i for i in range(0, len(hexRings))]
        envDict['RADIUS'] = [r+1 for r in radii]

        glucose = environments['glucose']   # fmol/um^3
        oxygen = environments['oxygen']     # mmHg
        tgfa = environments['tgfa']         # pg/cm^3
        IL2 = environments['IL-2']          # molecules IL-2/cm^3

        for time in range(0, len(T)):

            # Add information to dictionary
            envDict['TIME'].append(T[time])
            envDict['GLUCOSE'].append(glucose[SEED][time][HEIGHT])
            envDict['OXYGEN'].append(oxygen[SEED][time][HEIGHT])
            envDict['TGFA'].append(tgfa[SEED][time][HEIGHT])
            envDict['IL-2'].append(IL2[SEED][time][HEIGHT])

            # Initialize counters
            glucTotal = 0   # fmol
            oxyTotal = 0    # mmHg
            tgfaTotal = 0   # pg
            IL2Total = 0    # molecules

            volTotal = 0    # um^3

            for n in range(0, len(hexRings)):
                numHexes = hexRings[n]
                r = radii[n]
                volTotal += numHexes*HEX_VOL_UM
                glucTotal += numHexes*HEX_VOL_UM*glucose[SEED][time][HEIGHT][r]
                oxyTotal += numHexes*HEX_VOL_UM*oxygen[SEED][time][HEIGHT][r]
                tgfaTotal += numHexes*HEX_VOL_UM*1E-12*tgfa[SEED][time][HEIGHT][r]
                IL2Total += numHexes*HEX_VOL_UM*1E-12*IL2[SEED][time][HEIGHT][r]

            oxyTotal = oxyTotal/volTotal

            envDict['GLUCOSE TOTAL'].append(glucTotal)  # fmol
            envDict['OXYGEN TOTAL'].append(oxyTotal)    # average mmHg
            envDict['TGFA TOTAL'].append(tgfaTotal)     # pg
            envDict['IL-2 TOTAL'].append(IL2Total)      # molecules

            IL2pg = (1E12*15500*IL2Total)/(6.022E23)

            envDict['GLUCOSE TOTAL CONC'].append(glucTotal/volTotal/1E-12)  # fmol/ml
            envDict['TGFA TOTAL CONC'].append(tgfaTotal/volTotal/1E-12)     # pg/ml
            envDict['IL-2 TOTAL CONC'].append(IL2pg/volTotal/1E-12)         # pg/ml

        # Add tumor information to full simulation dataframe
        envDF = envDF.append(envDict, ignore_index=True)

    return envDF

def spatial_sim(spatialsDF, agents, T, R, C, TUMORID, SEED):

    # Make simulation dict
    spatialDict = make_spatial_dict()

    # Add simulation information to dict
    spatialDict['TUMOR ID'] = TUMORID
    spatialDict['SEED'] = int(SEED)

    tsplit = TUMORID.split('_')

    spatialDict['PLATE'] = tsplit[1]
    spatialDict['NUTRIENTS'] = 'CONSTANT' if tsplit[0] == 'VITRO' else 'GRAPH VASCULATURE'
    spatialDict['DOSE'] = int(tsplit[4])
    spatialDict['TREAT RATIO'] = tsplit[5].replace('-',':') if tsplit[5] != "NA" else tsplit[5]
    spatialDict['CAR AFFINITY'] = float(tsplit[6]) if tsplit[6] != "NA" else tsplit[6]
    spatialDict['ANTIGENS CANCER'] = int(tsplit[7])
    spatialDict['ANTIGENS HEALTHY'] = int(tsplit[8]) if tsplit[8] != "NA" else tsplit[8]
    spatialDict['RADIUS'] = [r+1 for r in range(0, R)]

    # CELL POP NUMBERS
    POP_CANCER = 0
    POP_HEALTHY = 1
    POP_CD4 = 2
    POP_CD8 = 3

    # These don't change even if pop number does
    INDEX_CANCER = 0
    INDEX_CANCERLIVE = 1
    INDEX_HEALTHY = 2
    INDEX_HEALTHYLIVE = 3
    INDEX_TCELL = 4
    INDEX_TCELLLIVE = 5
    INDEX_CD4 = 6
    INDEX_CD4LIVE = 7
    INDEX_CD8 = 8
    INDEX_CD8LIVE = 9

    POP_NAMES = ['CANCER', 'CANCER LIVE', 'HEALTHY', 'HEALTHY LIVE', 'T-CELL', 'T-CELL LIVE', 'CD4', 'CD4 LIVE', 'CD8', 'CD8 LIVE']
    POP_INDICES = [INDEX_CANCER, INDEX_CANCERLIVE, INDEX_HEALTHY, INDEX_HEALTHYLIVE, INDEX_TCELL, INDEX_TCELLLIVE, INDEX_CD4, INDEX_CD4LIVE, INDEX_CD8, INDEX_CD8LIVE]

    # Set height
    H = 0

    # Get hex rings
    hexRings = get_hex_rings(R)

    for time in range(0, len(T)):

        # Reset list of counts at each radii for each species
        counts = [[0 for j in range(0, R)] for i in range(0, 10)]  # [CANCER, CANCER LIVE, HEALTHY, HEALTHY LIVE, T-CELL, T-CELL LIVE, CD4, CD4 LIVE, CD8, CD8 LIVE]
        countsNorm = [[0 for j in range(0, R)] for i in range(0, 10)]  # [CANCER NORMALIZED, CANCER LIVE NORMALIZED, HEALTHY NORMALIZED, HEALTHY LIVE NORMALIZED, T-CELL NORMALIZED,
                                                                       # T-CELL LIVE NORMALIZED, CD4 NORMALIZED, CD4 LIVE NORMALIZED, CD8 NORMALIZED, CD8 LIVE NORMALIZED]

        for loc in range(0, len(C)):

            location = C[loc]
            radius = get_radius(location)
            radLocs = hexRings[radius]

            for pos in range(0, 54):
                if agents[time][H][loc]['pop'][pos] != -1:

                    # Find CANCER cells
                    if agents[time][H][loc]['pop'][pos] == POP_CANCER:
                        counts[INDEX_CANCER][radius] += 1

                        # Find CANCER LIVE cells
                        if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                            counts[INDEX_CANCERLIVE][radius] += 1

                    # Find HEALTHY cells
                    if agents[time][H][loc]['pop'][pos] == POP_HEALTHY:
                        counts[INDEX_HEALTHY][radius] += 1

                        # Find HEALTHY LIVE cells
                        if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                            counts[INDEX_HEALTHYLIVE][radius] += 1

                    # Find CD4 cells
                    if agents[time][H][loc]['pop'][pos] == POP_CD4:
                        counts[INDEX_CD4][radius] += 1
                        counts[INDEX_TCELL][radius] += 1

                        # Find CD4 LIVE cells
                        if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                            counts[INDEX_CD4LIVE][radius] += 1
                            counts[INDEX_TCELLLIVE][radius] += 1

                    # Find CD8 cells
                    if agents[time][H][loc]['pop'][pos] == POP_CD8:
                        counts[INDEX_CD8][radius] += 1
                        counts[INDEX_TCELL][radius] += 1

                        # Find CD8 LIVE cells
                        if agents[time][H][loc]['type'][pos] != 1 and agents[time][H][loc]['type'][pos] != 6:
                            counts[INDEX_CD8LIVE][radius] += 1
                            counts[INDEX_TCELLLIVE][radius] += 1

            for pop in POP_INDICES:
                countsNorm[pop][radius] = counts[pop][radius]/radLocs

        # Add information to dictionary
        spatialDict['TIME'].append(T[time])

        for p in range(0, len(POP_NAMES)):
            spatialDict[POP_NAMES[p]].append(counts[POP_INDICES[p]])
            spatialDict[POP_NAMES[p] + ' NORMALIZED'].append(countsNorm[POP_INDICES[p]])

    # Add tumor information to full simulation dataframe
    spatialsDF = spatialsDF.append(spatialDict, ignore_index=True)

    return spatialsDF

def analyze_sims(file):
    # Make simulation dataframe
    simsDF = make_analyze_df()

    if 'VITRO' in file:
        fileName = re.sub('.*VITRO', 'VITRO', file)

    else:
        fileName = re.sub('.*VIVO', 'VIVO', file)

    TUMORID = fileName.replace('.pkl', '')

    # Load tumor
    D, d, R, H, T, N, C, POPS, TYPES = ABM_load(file)

    for s in range(N):
        print("\t >" + TUMORID + "_0" + str(s) + "_ANALYZED")
        agents = D['agents'][s]
        simsDF = analyze_sim(simsDF, agents, T, C, TUMORID, s)

    return simsDF

def env_sims(file):
    # Make simulation environment dataframe
    envDF = make_env_df()

    if 'VITRO' in file:
        fileName = re.sub('.*VITRO', 'VITRO', file)

    else:
        fileName = re.sub('.*VIVO', 'VIVO', file)

    TUMORID = fileName.replace('.pkl', '')

    # Load tumor
    D, d, R, H, T, N, C, POPS, TYPES = ABM_load(file)

    # Load full pickle file
    parsedFile = pickle.load(open(file, "rb"))
    environments = parsedFile['environments']

    envDF = env_sim(envDF, environments, T, R, TUMORID, N)

    return envDF

def spatial_sims(file):
    # Make simulation environment dataframe
    spatialsDF = make_spatial_df()

    if 'VITRO' in file:
        fileName = re.sub('.*VITRO', 'VITRO', file)

    else:
        fileName = re.sub('.*VIVO', 'VIVO', file)

    TUMORID = fileName.replace('.pkl', '')

    # Load tumor
    D, d, R, H, T, N, C, POPS, TYPES = ABM_load(file)

    for s in range(N):
        print("\t >" + TUMORID + "_0" + str(s) + "_SPATIAL")
        agents = D['agents'][s]
        spatialsDF = spatial_sim(spatialsDF, agents, T, R, C, TUMORID, s)

    return spatialsDF

# ---------------------- MAIN ---------------------------------

if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()

    # Get files
    PKLFILES = get_pkl_files(args.files)

    # pd.set_option('display.max_rows', None)
    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.width', None)
    # pd.set_option('display.max_colwidth', -1)

    for file in PKLFILES:

        if 'VITRO' in file:
            fileName = re.sub('.*VITRO', 'VITRO', file)

        else:
            fileName = re.sub('.*VIVO', 'VIVO', file)

        TUMORID = fileName.replace('.pkl', '')

        print(TUMORID)

        if args.analyze:
            simsDF = analyze_sims(file)
            # print(simsDF)

            if args.saveLoc != '':
                with open(args.saveLoc + TUMORID + '_ANALYZED.pkl', 'wb') as f:
                    pickle.dump(simsDF, f)

        if args.environment:
            envsDF = env_sims(file)
            # print(envsDF)

            if args.saveLoc != '':
                with open(args.saveLoc + TUMORID + '_ENVIRONMENT.pkl', 'wb') as f:
                    pickle.dump(envsDF, f)

        if args.spatial:
            spatialsDF = spatial_sims(file)
            # print(spatialsDF)

            if args.saveLoc != '':
                with open(args.saveLoc + TUMORID + '_SPATIAL.pkl', 'wb') as f:
                    pickle.dump(spatialsDF, f)