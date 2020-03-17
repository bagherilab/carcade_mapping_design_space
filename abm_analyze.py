import ABM
from abm_parse import load as ABM_load
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
and extracts the data into a dataframe in the form:

    TUMOR ID | SEED | PLATE | DAMAGE | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

and saves it to a .pkl where the DATA are the following list of inormation:
    
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

each in the format of a list of the value of the specified informat at each point in time.

Usage:
    python abm_analyze.py FILES [--saveLoc SAVELOC]

    FILES
        Path to .pkl or directory
    [--saveLoc]
        Location of where to save file, default will save here
'''

# -------------- PARSING AND INPUT/OUTPUT FUNCTIONS -------------

def get_parser():

    # Setup argument parser.
    parser = ArgumentParser(description="Receord all ABM CART simulation data into a dataframe")
    parser.add_argument(dest="files", help="Path to .pkl file or directory")
    parser.add_argument("--saveLoc", default='', dest="saveLoc",
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

def get_json_files(arg):
    if arg[-1] == "/" or arg[-1] == "\\":
        return [arg + f for f in os.listdir(arg) if ABM.is_json(f)]
    else:
        assert ABM.is_json(arg)
        return [arg]



# ---------------- MAKE EMPTY CONTAINERS ---------------------

def make_df():
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

def make_dict():
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

# ---------------------- ANALYZE TUMOR ---------------------------

def analyze_sim(simsDF, agents, T, TUMORID, SEED):

    # Make simulation dict
    simsDict = make_dict()

    # Add simulation information to dict
    simsDict['TUMOR ID'] = TUMORID
    simsDict['SEED'] = int(SEED)

    tsplit = TUMORID.split('_')

    simsDict['PLATE'] = tsplit[1]
    simsDict['NUTRIENTS'] = 'CONSTANT' if tsplit[1] == 'VITRO' else 'GRAPH VASCULATURE'
    simsDict['DOSE'] = int(tsplit[4])
    simsDict['TREAT RATIO'] = tsplit[5].replace('-',':') if tsplit[5] != "NA" else tsplit[5]
    simsDict['CAR AFFINITY'] = float(tsplit[6]) if tsplit[6] != "NA" else tsplit[6]
    simsDict['ANTIGENS CANCER'] = int(tsplit[7])
    simsDict['ANTIGENS HEALTHY'] = int(tsplit[8]) if tsplit[8] != "NA" else tsplit[8]

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

        for loc in range(0, 3367):
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
        simsDict['TIME'].append(T[time])

        for p in range(0, len(POP_NAMES)):
            simsDict[POP_NAMES[p]].append(counts[POP_INDICES[p]])

            if 'LIVE' not in POP_NAMES[p]:

                simsDict['AVG CELL CYCLES ' + POP_NAMES[p]].append(cycles[POP_INDICES[p]])
                simsDict['CELL VOLUMES ' + POP_NAMES[p]].append(volumes[POP_INDICES[p]])

                if POP_NAMES[p] == 'CANCER' or POP_NAMES[p] == 'HEALTHY':
                    for s in range(0, len(STATES_CANCER)):
                        key = STATES_CANCER[s] + ' ' + POP_NAMES[p]
                        simsDict[key].append(types[POP_INDICES[p]][STATES_CANCER_INDICES[s]])
                        simsDict[key + ' %'].append(typesFrac[POP_INDICES[p]][STATES_CANCER_INDICES[s]])

                if POP_NAMES[p] == 'T-CELL' or POP_NAMES[p] == 'CD4' or POP_NAMES[p] == 'CD8':
                    for s in range(0, len(STATES_TCELL)):
                        if (POP_NAMES[p] == 'CD4' and STATES_TCELL[s] == 'CYTOT'):
                            continue
                        elif (POP_NAMES[p] == 'CD8' and STATES_TCELL[s] == 'STIMU'):
                            continue
                        else:
                            key = STATES_TCELL[s] + ' ' + POP_NAMES[p]
                            simsDict[key].append(types[POP_INDICES[p]][STATES_TCELL_INDICES[s]])
                            simsDict[key + ' %'].append(typesFrac[POP_INDICES[p]][STATES_TCELL_INDICES[s]])


            else:
                if POP_NAMES[p] == 'CANCER LIVE' or POP_NAMES[p] == 'HEALTHY LIVE':
                    if counts[POP_INDICES[p-1]] == 0:
                        simsDict[POP_NAMES[p] + ' %'].append(0)
                    else:
                        simsDict[POP_NAMES[p] + ' %'].append(counts[POP_INDICES[p]] / counts[POP_INDICES[p-1]])

                    for s in range(0, len(STATES_TISSUELIVE)):
                        key = STATES_TISSUELIVE[s] + ' ' + POP_NAMES[p]
                        simsDict[key].append(types[POP_INDICES[p]][STATES_TISSUELIVE_INDICES[s]])
                        simsDict[key + ' %'].append(typesFrac[POP_INDICES[p]][STATES_TISSUELIVE_INDICES[s]])

    # Add tumor information to full simulation dataframe
    simsDF = simsDF.append(simsDict, ignore_index=True)

    return simsDF

def analyze_sims(file):
    # Make simulation dataframe
    simsDF = make_df()

    if 'VITRO' in file:
        fileName = re.sub('.*VITRO', 'VITRO', file)

    else:
        fileName = re.sub('.*VIVO', 'VIVO', file)

    TUMORID = fileName.replace('.pkl', '')

    # Load tumor
    D, d, R, H, T, N, C, POPS, TYPES = ABM_load(file)

    for s in range(N):
        print("\t >" + TUMORID + "_0" + str(s))
        agents = D['agents'][s]
        simsDF = analyze_sim(simsDF, agents, T, TUMORID, s)

    return simsDF


# ---------------------- MAIN ---------------------------------

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

        TUMORID = fileName.replace('.pkl', '')

        print(TUMORID)

        simsDF = analyze_sims(file)

        if args.saveLoc != '':
            pickle.dump(simsDF, open(args.saveLoc + TUMORID + '_ANALYZED.pkl', "wb"))