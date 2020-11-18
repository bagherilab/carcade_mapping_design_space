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
ABM_ANALYZE takes a directory of (or a single) .LYSIS.json simulation files
and extracts the data into a dataframe. The resulting file will contain a data frame in the form:

    TUMOR ID | SEED | PLATE | DAMAGE | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

and saves it to a .pkl where the DATA are the following list of information:
    
    TIME
    RADIUS
    TISSUE SEEDED
    CANCER SEEDED
    HEALTHY SEEDED
    T-CELL SEEDED
    CD4 SEEDED
    CD8 SEEDED
    TIME EXACT                              list of exact times of cell death
    RADIUS LYSED EXACT                      list of exact radius cell killed over time 
    TISSUE LYSED EXACT                      list of exact cells killed over time
    TISSUE LYSED CUMULATIVE EXACT           list of cumulative tissue cells killed at exact times
    CANCER LYSED CUMULATIVE EXACT           list of cumulative cancer cells killed at exact times
    HEALTHY LYSED CUMULATIVE EXACT          list of cumulative healthy cells killed at exact times
    TISSUE LYSED TOTAL                      list of total tissue cells killed at collection time points
    CANCER LYSED TOTAL                      list of total cancer cells killed at collection time points
    HEALTHY LYSED TOTAL                     list of total healthy cells killed at collection time points
    
each in the format of a list of the value of the specified information at each point in time.

Usage:
    python abm_lysis.py FILES [--saveLoc SAVELOC]

    FILES
        Path to .LYSIS.json or directory
    [--saveLoc SAVELOC]
        Location of where to save file, default won't save files
'''

# -------------- PARSING AND INPUT/OUTPUT FUNCTIONS -------------

def get_parser():

    # Setup argument parser.
    parser = ArgumentParser(description="Record all ABM CART simulation data into a dataframe")
    parser.add_argument(dest="files", help="Path to .pkl file or directory")
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

def make_lysis_df():
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
               'TISSUE SEEDED',
               'CANCER SEEDED',
               'HEALTHY SEEDED',
               'T-CELL SEEDED',
               'CD4 SEEDED',
               'CD8 SEEDED',
               'TIME EXACT',
               'RADIUS LYSED EXACT',
               'TISSUE LYSED EXACT',
               'TISSUE LYSED CUMULATIVE EXACT',
               'CANCER LYSED CUMULATIVE EXACT',
               'HEALTHY LYSED CUMULATIVE EXACT',
               'TISSUE LYSED TOTAL',
               'CANCER LYSED TOTAL',
               'HEALTHY LYSED TOTAL'
    ]
    lysisDF = pd.DataFrame(columns=columns)

    return lysisDF

def make_lysis_dict():
    lysisDict = {   'TUMOR ID': None,
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
                    'TISSUE SEEDED': [],
                    'CANCER SEEDED': [],
                    'HEALTHY SEEDED': [],
                    'T-CELL SEEDED': [],
                    'CD4 SEEDED': [],
                    'CD8 SEEDED': [],
                    'TIME EXACT': [],
                    'RADIUS LYSED EXACT': [],
                    'TISSUE LYSED EXACT': [],
                    'TISSUE LYSED CUMULATIVE EXACT': [],
                    'CANCER LYSED CUMULATIVE EXACT': [],
                    'HEALTHY LYSED CUMULATIVE EXACT': [],
                    'TISSUE LYSED TOTAL': [],
                    'CANCER LYSED TOTAL': [],
                    'HEALTHY LYSED TOTAL': []
    }

    return lysisDict

# ---------------------- ANALYZE TUMOR ---------------------------

def lysis_sim(lysisDF, LYSIS_JSON, TUMORID):

    # Make simulation dict
    lysisDict = make_lysis_dict()

    # CELL POP NUMBERS
    POP_CANCER = 0
    POP_HEALTHY = 1
    POP_CD4 = 2
    POP_CD8 = 3

    # Add simulation information to dict
    lysisDict['TUMOR ID'] = TUMORID
    lysisDict['SEED'] = LYSIS_JSON['seed']

    tsplit = TUMORID.split('_')

    lysisDict['PLATE'] = tsplit[1]
    lysisDict['NUTRIENTS'] = 'CONSTANT' if tsplit[0] == 'VITRO' else 'GRAPH VASCULATURE'
    lysisDict['DOSE'] = int(tsplit[4])
    lysisDict['TREAT RATIO'] = tsplit[5].replace('-',':') if tsplit[5] != "NA" else tsplit[5]
    lysisDict['CAR AFFINITY'] = float(tsplit[6]) if tsplit[6] != "NA" else tsplit[6]
    lysisDict['ANTIGENS CANCER'] = int(tsplit[7])
    lysisDict['ANTIGENS HEALTHY'] = int(tsplit[8]) if tsplit[8] != "NA" else tsplit[8]
    lysisDict['RADIUS'] = [r for r in range(1, int(LYSIS_JSON['config']['size']['radius']) + 1)]
    lysisDict['CANCER SEEDED'] = int(LYSIS_JSON['config']['pops'][POP_CANCER][-1])
    lysisDict['HEALTHY SEEDED'] = int(LYSIS_JSON['config']['pops'][POP_HEALTHY][-1])
    lysisDict['TISSUE SEEDED'] = lysisDict['CANCER SEEDED'] + lysisDict['HEALTHY SEEDED']

    for helper in LYSIS_JSON['helpers']:
        if helper['type'] == 'TREAT':
            for pop in helper['pops']:
                if pop[0] == POP_CD4:
                    lysisDict['CD4 SEEDED'] = int(pop[1])
                if pop[0] == POP_CD8:
                    lysisDict['CD8 SEEDED'] = int(pop[1])
            lysisDict['T-CELL SEEDED'] = lysisDict['CD4 SEEDED'] + lysisDict['CD8 SEEDED']

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

    for timepoint in LYSIS_JSON['timepoints']:

        tissueTotal = 0
        cancerTotal = 0
        healthyTotal = 0

        if float(timepoint['time']) == float(LYSIS_JSON['config']['days']):

            timeExact = []
            radiusExact = []
            tissueExact = []
            tissueTotalExact = []
            cancerTotalExact = []
            healthyTotalExact = []

            for cell in timepoint['cells']:
                timeExact.append(int(cell[0]))
                coord = cell[1]
                coord = coord[0:-2]
                radius = get_radius(coord)
                radiusExact.append(radius)
                tissueExact.append(cell[2][1])

                if cell[2][1] == POP_CANCER:
                    cancerTotal += 1
                    tissueTotal += 1
                    cancerTotalExact.append(cancerTotal)
                    healthyTotalExact.append(healthyTotal)
                    tissueTotalExact.append(tissueTotal)

                if cell[2][1] == POP_HEALTHY:
                    healthyTotal += 1
                    tissueTotal += 1
                    cancerTotalExact.append(cancerTotal)
                    healthyTotalExact.append(healthyTotal)
                    tissueTotalExact.append(tissueTotal)

            lysisDict['TIME EXACT'] = timeExact
            lysisDict['RADIUS LYSED EXACT'] = radiusExact
            lysisDict['TISSUE LYSED EXACT'] = tissueExact
            lysisDict['TISSUE LYSED CUMULATIVE EXACT'] = tissueTotalExact
            lysisDict['CANCER LYSED CUMULATIVE EXACT'] = cancerTotalExact
            lysisDict['HEALTHY LYSED CUMULATIVE EXACT'] = healthyTotalExact

        else:
            for cell in timepoint['cells']:
                if cell[2][1] == POP_CANCER:
                    cancerTotal += 1
                    tissueTotal += 1
                if cell[2][1] == POP_HEALTHY:
                    healthyTotal += 1
                    tissueTotal += 1

        # Add information to dictionary
        lysisDict['TIME'].append(timepoint['time'])
        lysisDict['TISSUE LYSED TOTAL'].append(tissueTotal)
        lysisDict['CANCER LYSED TOTAL'].append(cancerTotal)
        lysisDict['HEALTHY LYSED TOTAL'].append(healthyTotal)

    # Add tumor information to full simulation dataframe
    lysisDF = lysisDF.append(lysisDict, ignore_index=True)

    return lysisDF

def lysis_sims(file):
    # Make simulation environment dataframe
    lysisDF = make_lysis_df()

    if 'VITRO' in file:
        fileName = re.sub('.*VITRO', 'VITRO', file)

    else:
        fileName = re.sub('.*VIVO', 'VIVO', file)

    TUMORID = fileName.replace('.LYSIS.json', '')

    # Load tumor
    LYSIS_JSON = ABM.load_json(file)

    lysisDF = lysis_sim(lysisDF, LYSIS_JSON, TUMORID)

    return lysisDF

# ---------------------- MAIN ---------------------------------

if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()

    # Get files
    LYSISFILES = get_json_files(args.files)

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', -1)

    for file in LYSISFILES:

        if 'VITRO' in file:
            fileName = re.sub('.*VITRO', 'VITRO', file)

        else:
            fileName = re.sub('.*VIVO', 'VIVO', file)

        TUMORID = fileName.replace('.LYSIS.json', '')

        print(TUMORID)

        lysisDF = lysis_sims(file)

        if args.saveLoc != '':
            pickle.dump(lysisDF, open(args.saveLoc + TUMORID + '_LYSED.pkl', "wb"))