import ABM
from abm_analyze import make_analyze_df
from abm_analyze import make_env_df
from abm_analyze import make_spatial_df
from abm_lysis import make_lysis_df
import os
import pickle
import pandas as pd
import bz2
from argparse import ArgumentParser

__author__ = "Alexis N. Prybutok"
__email__ = "aprybutok@u.northwestern.edu"

'''
ABM_SUBSET takes a directory of (or a single) .pkl analyzed files
and extracts the specified subset of the data into a dataframe. If the ANALYZE type is specified,
data is given in the form:

    TUMOR ID | SEED | PLATE | NUTRIENTS | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

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
    STARV T-CELL
    STARV CD4
    STARV CD8
    STARV T-CELL %
    STARV CD4 %
    STARV CD8 %
    PAUSE T-CELL
    PAUSE CD4
    PAUSE CD8
    PAUSE T-CELL %
    PAUSE CD4 %
    PAUSE CD8 %
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

If ENVIRONMENT type is specified, will produce a corresponding environment file that contains
extracts the data into a dataframe in the form:

    TUMOR ID | SEED | PLATE | DAMAGE | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

and saves it to a .pkl where the DATA are the following list of inormation:

    TIME
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
    
where GLUCOSE, OXYGEN, TGFA, and IL-2 are in the format of a list of a list of concentrations at each radius.
TOTAL columns are the total (absolute) amount of molecule in total area.

If SPATIAL type is specified, it will produce a corresponding spatial file that contains
extracts the data into a dataframe in the form:

    TUMOR ID | SEED | PLATE | DAMAGE | DOSE | TREAT RATIO | CAR AFFINITY | ANTIGENS CANCER | ANTIGENS HEALTHY | DATA

and saves it to a .pkl where the DATA are the following list of information:

    TIME
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
    HEALTHY
    HEALTHY LIVE
    T-CELL
    T-CELL LIVE
    CD4
    CD4 LIVE
    CD8
    CD8 LIVE
    
where each cell population is in the format of a list of a list of counts at each radius.

Usage:
    python abm_subset.py FILES XMLNAME TYPE SAVELOC [--subset SUBSET]

    FILES
        Path to .pkl or directory
    XMLNAME
        Name of XML file used to generate dataset for purpose of saving
    TYPE
        Type of analysis file given (analyzed, environment, spatial, or lysis)
    SAVELOC
        Location of where to save file, default will save here
    [--subset SUBSET]
        Feature(s) by which to subselect data and plot only subset of data
        with given specified value (ex: [DOSE,500] will only plot all data
        with a dose of 500 cells) given as semicolon separated list of sets of 
        feature name and value of feature by which to plot separated by a colon.
        More than one feature can be chosen per plot, as indicated by commas 
        between sets, but each new list will make a different set of plots. (ex: 
        [(DOSE:500),(TREAT RATIO:50-50)];[(DOSE:500),(CAR ANTIGEN:1000)] will 
        make two sets of plot, where the first is only data with both DOSE 
        500 and TREAT RATIO 50:50, and the second is only data with DOSE 500
        and ANTIGEN 1000) (default: none, all data plotted on one graph)
    [--states]
        Flag indicating only to save state data (not list column data).
'''

# -------------- PARSING AND INPUT/OUTPUT FUNCTIONS -------------

def get_parser():

    # Setup argument parser.
    parser = ArgumentParser(description="Plot ABM data from dataframe")
    parser.add_argument(dest="files", help="Path to .pkl file or directory")
    parser.add_argument(dest="xmlName", help="Name of XML file used to generate dataset for purpose of saving")
    parser.add_argument(dest="type", help="Type of analyzed file given (analyzed, environment, or subset")
    parser.add_argument(dest="saveLoc", help="Location of where to save file, default will save here")
    parser.add_argument("--subset", default='', dest="subset",
                        help="Feature(s) by which to subselect data and plot only subset of data "
                             "with given specified value (ex: [DOSE,500] will only plot all data"
                             "with a dose of 500 cells) given as semicolon separated list of sets of "
                             "feature name and value of feature by which to plot separated by a colon. "
                             "More than one feature can be chosen per plot, as indicated by commas "
                             "between sets, but each new list will make a different set of plots. (ex: "
                             "[(DOSE:500),(TREAT RATIO:50-50)];[(DOSE:500),(CAR ANTIGENS:1000)] will "
                             "make two sets of plot, where the first is only data with both DOSE "
                             "500 and TREAT RATIO 50:50, and the second is only data with DOSE 500"
                             "and CAR ANTIGEN 1000) (default: none, all data plotted on one graph)")
    parser.add_argument("--states", default=False, dest="states", action='store_true',
                        help="Flag indicating only to save state data (not list column data).")

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

# ------------------VARIABLES--------------------------------

LIST_COLUMNS = ["AVG CELL CYCLES CANCER", "AVG CELL CYCLES HEALTHY", "AVG CELL CYCLES T-CELL", "AVG CELL CYCLES CD4", "AVG CELL CYCLES CD8",
                "CELL VOLUMES CANCER","CELL VOLUMES HEALTHY", "CELL VOLUMES T-CELL", "CELL VOLUMES CD4", "CELL VOLUMES CD8"]

# ---------------------- MAIN ---------------------------------

if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()

    # Parse feature plotting
    SUBSETS = []
    if args.subset != '':
        ssplit = str(args.subset).split(";")
        for s in ssplit:
            subset = []
            psplit = s.split(",")
            for p in psplit:
                vsplit = p.split(":")
                label = vsplit[0].strip('[]()')
                value = vsplit[1].strip('[]()')
                if label == 'ANTIGENS CANCER' or label == 'DOSE':
                    value = int(value)
                elif label == 'CAR AFFINITY':
                    value = float(value)
                elif label == 'TREAT RATIO':
                    value = value.replace('-', ':')
                subset.append([label, value])
            SUBSETS.append(subset)

    # Get files
    PKLFILES = get_pkl_files(args.files)

    # Set type
    TYPE = args.type.upper()
    print(TYPE)

    # If no subsets specified, put all data in one file
    # Warning: Highly likely you will run into memory error
    if args.subset == '':

        # Make dataframe based on type given
        if TYPE == 'ENVIRONMENT':
            simsDF = make_env_df()
            untreatedDF = make_env_df()
        elif TYPE == 'SPATIAL':
            simsDF = make_spatial_df()
            untreatedDF = make_spatial_df()
        elif TYPE == 'LYSED':
            simsDF = make_lysis_df()
            untreatedDF = make_lysis_df()
        else:
            if TYPE != 'ANALYZED':
                print('No valid type specified. Assuming ANALYZE.')
                TYPE == 'ANALYZED'
            simsDF = make_analyze_df()
            untreatedDF = make_analyze_df()

        if TYPE == 'ANALYZED' and args.states:
            simsDF = simsDF.drop(LIST_COLUMNS, axis=1)
            untreatedDF = untreatedDF.drop(LIST_COLUMNS, axis=1)

        untreatedName = ''

        # Grab all files
        for file in PKLFILES:
            if '0_NA_NA' in file:
                with open(file, 'rb') as f:
                    simDF = pickle.load(f)
                if args.states:
                    simDF = simDF.loc[:,'TUMOR ID':'PAUSE CD8 %']
                untreatedDF = untreatedDF.append(simDF, ignore_index=True)
                untreatedName = file.replace(args.files, '')
            else:
                with open(file, 'rb') as f:
                    simDF = pickle.load(f)
                if args.states:
                    simDF = simDF.loc[:, 'TUMOR ID':'PAUSE CD8 %']
                simsDF = simsDF.append(simDF, ignore_index=True)

        # Set up options dictionary for naming save file
        optionsDict = {
            'DOSE': 'X',
            'TREAT RATIO': 'X',
            'CAR AFFINITY': 'X',
            'ANTIGENS CANCER': 'X',
            'ANTIGENS HEALTHY': 'X'
        }

        # Fill out optionsDict with NA if true in untreated file
        untreatedSplit = untreatedName.split('_')
        if untreatedSplit[-2] == 'NA' and TYPE != 'LYSED':
            optionsDict['ANTIGENS HEALTHY'] = 'NA'

        if untreatedSplit[-3] == 'NA' and TYPE == 'LYSED':
            optionsDict['ANTIGENS HEALTHY'] = 'NA'

        # Add untreated control to set
        simsDF = simsDF.append(untreatedDF, ignore_index=True)

        # Construct save file name
        name = optionsDict['DOSE'] + '_' + optionsDict['TREAT RATIO'] + '_' + \
               optionsDict['CAR AFFINITY'] + '_' + optionsDict['ANTIGENS CANCER'] + '_' + optionsDict['ANTIGENS HEALTHY'] + '_'

        # Save file
        try:
            if args.states:
                with open(args.saveLoc + args.xmlName + '_' + name + 'STATES_' + TYPE + '.pkl', 'wb') as f:
                    pickle.dump(simsDF, f)
            else:
                with open(args.saveLoc + args.xmlName + '_' + name + TYPE + '.pkl', 'wb') as f:
                    pickle.dump(simsDF, f)
            print('Data saved!')
        except MemoryError:
            # Remove empty pickle file
            os.remove(args.saveLoc + args.xmlName + '_' + name + TYPE + '.pkl')
            print('MemoryError when attempting to save full pkl.')
            if TYPE == 'ANALYZED' and not args.states:
                print('Saving all non-list columns.')
                with open(args.saveLoc + args.xmlName + '_' + name + 'STATES_ANALYZED.pkl', 'wb') as f:
                    simsDFstates = simsDF.loc[:, 'TUMOR ID':'PAUSE CD8 %']
                    pickle.dump(simsDFstates, f)
                    #pickle.dump(simsDF[:, 0:-10], f)
                print('Saving all list columns in separate files.')
                simsDFcycle = simsDF[["TUMOR ID", "SEED", "PLATE", "NUTRIENTS", "DOSE", "TREAT RATIO", "CAR AFFINITY",
                     "ANTIGENS CANCER", "ANTIGENS HEALTHY", "TIME", "AVG CELL CYCLES CANCER", "AVG CELL CYCLES HEALTHY",
                     "AVG CELL CYCLES T-CELL", "AVG CELL CYCLES CD4", "AVG CELL CYCLES CD8"]]
                simsDFvol = simsDF[["TUMOR ID", "SEED", "PLATE", "NUTRIENTS", "DOSE", "TREAT RATIO", "CAR AFFINITY",
                                    "ANTIGENS CANCER", "ANTIGENS HEALTHY", "TIME", "CELL VOLUMES CANCER","CELL VOLUMES HEALTHY",
                                    "CELL VOLUMES T-CELL", "CELL VOLUMES CD4", "CELL VOLUMES CD8"]]
                with open(args.saveLoc + args.xmlName + '_' + name + 'CYCLES_ANALYZED.pkl', 'wb') as f:
                    pickle.dump(simsDFcycle, f)
                with open(args.saveLoc + args.xmlName + '_' + name + 'VOLUMES_ANALYZED.pkl', 'wb') as f:
                    pickle.dump(simsDFvol, f)
                print('Data saved!')

    else:

        for subset in SUBSETS:

            fileCount = 0

            # Make dataframe based on type given
            if args.type.upper() == 'ENVIRONMENT':
                simsDF = make_env_df()
                untreatedDF = make_env_df()
            elif args.type.upper() == 'SPATIAL':
                simsDF = make_spatial_df()
                untreatedDF = make_spatial_df()
            elif args.type.upper() == 'LYSED':
                simsDF = make_lysis_df()
                untreatedDF = make_lysis_df()
            else:
                if args.type.upper() != 'ANALYZED':
                    print('No valid type specified. Assuming ANALYZED.')
                simsDF = make_analyze_df()
                untreatedDF = make_analyze_df()

            if TYPE == 'ANALYZED' and args.states:
                simsDF = simsDF.drop(LIST_COLUMNS, axis=1)
                untreatedDF = untreatedDF.drop(LIST_COLUMNS, axis=1)

            # Set up options dictionary for sorting files and naming save file
            optionsDict = {
                'DOSE': 'X',
                'TREAT RATIO': 'X',
                'CAR AFFINITY': 'X',
                'ANTIGENS CANCER': 'X',
                'ANTIGENS HEALTHY': 'X'
            }

            # Fill out options dictionary with set specifications
            for s in subset:
                optionsDict[s[0]] = str(s[1]).replace(':','-')

            untreatedName = ''

            print('Adding the following files to set:')

            # Find only files with specificed set information
            for file in PKLFILES:
                if '0_NA_NA' in file:
                    print('\t' + file.replace(args.files, ''))
                    with open(file, 'rb') as f:
                        fileCount += 1
                        simDF = pickle.load(f)
                        if args.states:
                            simDF = simDF.loc[:, 'TUMOR ID':'PAUSE CD8 %']
                    untreatedDF = untreatedDF.append(simDF, ignore_index=True)
                    untreatedName = file.replace(args.files, '')
                else:

                    filename = file.replace(args.files, '')
                    fileSplit = filename.split("_")

                    # Check file name for set specifications and open file and add to set only if match
                    if fileSplit[4] == optionsDict['DOSE'] or optionsDict['DOSE'] == 'X':
                        if fileSplit[5] == optionsDict['TREAT RATIO'] or optionsDict['TREAT RATIO'] == 'X':
                            try:
                                optionsAffinity = float(optionsDict['CAR AFFINITY'])
                            except ValueError:
                                optionsAffinity = 'X'
                            if float(fileSplit[6]) == optionsAffinity or optionsDict['CAR AFFINITY'] == 'X':
                                if fileSplit[7] == optionsDict['ANTIGENS CANCER'] or optionsDict['ANTIGENS CANCER'] == 'X':
                                    if fileSplit[8] == optionsDict['ANTIGENS HEALTHY'] or optionsDict['ANTIGENS HEALTHY'] == 'X' or fileSplit[8] == 'NA':

                                        print('\t' + filename)

                                        # Open file and save to set
                                        with open(file, 'rb') as f:
                                            simDF = pickle.load(f)
                                            fileCount += 1
                                            if args.states:
                                                simDF = simDF.loc[:, 'TUMOR ID':'PAUSE CD8 %']
                                        simsDF = simsDF.append(simDF, ignore_index=True)

            # Fill out optionsDict with NA if true in untreated file
            untreatedSplit = untreatedName.split('_')
            if untreatedSplit[-2] == 'NA' and TYPE != 'LYSED':
                optionsDict['ANTIGENS HEALTHY'] = 'NA'

            if untreatedSplit[-3] == 'NA' and TYPE == 'LYSED':
                optionsDict['ANTIGENS HEALTHY'] = 'NA'


            # Add untreated control to set
            simsDF = simsDF.append(untreatedDF, ignore_index=True)

            # Construct save file name
            name = optionsDict['DOSE'] + '_' + optionsDict['TREAT RATIO'] + '_' + \
                   optionsDict['CAR AFFINITY'] + '_' + optionsDict['ANTIGENS CANCER'] + '_' + \
                   optionsDict['ANTIGENS HEALTHY'] + '_'

            # Save set
            print('Saving set...')
            if (fileCount <= 6 and TYPE == 'ANALYZED') or (TYPE != 'ANALYZED'):
                try:
                    if args.states:
                        with open(args.saveLoc + args.xmlName + '_' + name + 'STATES_' + TYPE + '.pkl', 'wb') as f:
                            pickle.dump(simsDF, f)
                    else:
                        with open(args.saveLoc + args.xmlName + '_' + name + TYPE + '.pkl', 'wb') as f:
                            pickle.dump(simsDF, f)
                    print('Set saved!')
                except MemoryError:
                    # Remove empty pickle file
                    os.remove(args.saveLoc + args.xmlName + '_' + name + TYPE + '.pkl')
                    print('MemoryError when attempting to save full pkl.')
                    if TYPE == 'ANALYZED':
                        print('Saving all non-list columns.')
                        with open(args.saveLoc + args.xmlName + '_' + name + 'STATES_ANALYZED.pkl', 'wb') as f:
                            pickle.dump(simsDF.iloc[:,0:-10], f)
                        print('Saving all list columns in separate files.')
                        simsDFcycle = simsDF[
                            ["TUMOR ID", "SEED", "PLATE", "NUTRIENTS", "DOSE", "TREAT RATIO", "CAR AFFINITY",
                             "ANTIGENS CANCER", "ANTIGENS HEALTHY", "TIME", "AVG CELL CYCLES CANCER", "AVG CELL CYCLES HEALTHY",
                             "AVG CELL CYCLES T-CELL", "AVG CELL CYCLES CD4", "AVG CELL CYCLES CD8"]]
                        simsDFvol = simsDF[
                            ["TUMOR ID", "SEED", "PLATE", "NUTRIENTS", "DOSE", "TREAT RATIO", "CAR AFFINITY",
                             "ANTIGENS CANCER", "ANTIGENS HEALTHY", "TIME", "CELL VOLUMES CANCER", "CELL VOLUMES HEALTHY",
                             "CELL VOLUMES T-CELL", "CELL VOLUMES CD4", "CELL VOLUMES CD8"]]
                        with open(args.saveLoc + args.xmlName + '_' + name + 'CYCLES_ANALYZED.pkl', 'wb') as f:
                            pickle.dump(simsDFcycle, f)
                        with open(args.saveLoc + args.xmlName + '_' + name + 'VOLUMES_ANALYZED.pkl', 'wb') as f:
                            pickle.dump(simsDFvol, f)
                        print('Set saved!')
            else:
                print('Large number of files.')
                if TYPE == 'ANALYZED':
                    print('Saving all non-list columns.')
                    with open(args.saveLoc + args.xmlName + '_' + name + 'STATES_ANALYZED.pkl', 'wb') as f:
                        pickle.dump(simsDF.iloc[:, 0:-10], f)
                    print('Set saved!')