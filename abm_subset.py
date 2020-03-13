import ABM
from abm_analyze import make_df
import os
import pickle
from argparse import ArgumentParser

__author__ = "Alexis N. Prybutok"
__email__ = "aprybutok@u.northwestern.edu"

'''
ABM_SUBSET takes a directory of (or a single) .pkl analyzed files
and extracts the specified subset of the data into a dataframe in the form:

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
    python abm_subset.py FILES XMLNAME SAVELOC [--subset SUBSET]

    FILES
        Path to .pkl or directory
    XMLNAME
        Name of XML file used to generate dataset for purpose of saving
    SAVELOC
        Location of where to save file, default will save here
    [--subset SUBSET]
        Feature(s) by which to subselect data and plot only subset of data
        with given specified value (ex: [DOSE,500] will only plot all data
        with a dose of 500 cells) given as semicolon separated list of sets of 
        feature name and value of feature by which to plot separated by a colon.
        More than one feature can be chosen per plot, as indicated by commas 
        between sets, but each new list will make a different set of plots. (ex: 
        [(DOSE:500),(TREAT RATIO:50-50)];[(DOSE:500),(ANTIGEN:1000)] will 
        make two sets of plot, where the first is only data with both DOSE 
        500 and TREAT RATIO 50:50, and the second is only data with DOSE 500
        and ANTIGEN 1000) (default: none, all data plotted on one graph)
'''

# -------------- PARSING AND INPUT/OUTPUT FUNCTIONS -------------

def get_parser():

    # Setup argument parser.
    parser = ArgumentParser(description="Plot ABM data from dataframe")
    parser.add_argument(dest="files", help="Path to .pkl file or directory")
    parser.add_argument(dest="xmlName", help="Name of XML file used to generate dataset for purpose of saving")
    parser.add_argument(dest="saveLoc", help="Location of where to save file, default will save here")
    parser.add_argument("--subset", default='', dest="subset",
                        help="Feature(s) by which to subselect data and plot only subset of data "
                             "with given specified value (ex: [DOSE,500] will only plot all data"
                             "with a dose of 500 cells) given as semicolon separated list of sets of "
                             "feature name and value of feature by which to plot separated by a colon. "
                             "More than one feature can be chosen per plot, as indicated by commas "
                             "between sets, but each new list will make a different set of plots. (ex: "
                             "[(DOSE:500),(TREAT RATIO:50-50)];[(DOSE:500),(ANTIGEN:1000)] will "
                             "make two sets of plot, where the first is only data with both DOSE "
                             "500 and TREAT RATIO 50:50, and the second is only data with DOSE 500"
                             "and ANTIGEN 1000) (default: none, all data plotted on one graph)")

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


# ---------------------- MAIN ---------------------------------

if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()

    # Parse feature plotting
    SUBSETS = []
    if args.subset != '':
        ssplit = str(args.subset).split(";")
        for s in ssplit:
            set = []
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
                set.append([label, value])
            SUBSETS.append(set)

    # Get files
    PKLFILES = get_pkl_files(args.files)

    # Make dataframe
    simsDF = make_df()
    untreatedDF = make_df()


    for file in PKLFILES:
        if '0_NA_NA' in file:
            simDF = pickle.load(open(file, 'rb'))
            untreatedDF = untreatedDF.append(simDF, ignore_index=True)
        else:
            simDF = pickle.load(open(file, 'rb'))
            simsDF = simsDF.append(simDF, ignore_index=True)

    if args.subset == '':
        optionsDict = {
            'DOSE': 'X',
            'TREAT RATIO': 'X',
            'CAR AFFINITY': 'X',
            'ANTIGENS CANCER': 'X'
        }
        simsDF = simsDF.append(untreatedDF, ignore_index=True)
        name = optionsDict['DOSE'] + '_' + optionsDict['TREAT RATIO'] + '_' + \
               optionsDict['CAR AFFINITY'] + '_' + optionsDict['ANTIGENS CANCER'] + '_'
        pickle.dump(simsDF, open(args.saveLoc + args.xmlName + name + 'ANALYZED.pkl', "wb"))
    else:
        for set in SUBSETS:

            optionsDict = {
                'DOSE': 'X',
                'TREAT RATIO': 'X',
                'CAR AFFINITY': 'X',
                'ANTIGENS CANCER': 'X'
            }

            setDF = make_df()
            setDF = setDF.append(simsDF, ignore_index=True)
            for s in set:
                optionsDict[s[0]] = str(s[1]).replace(':','-')
                setDF = setDF.loc[setDF[s[0]] == s[1]]

            setDF = setDF.append(untreatedDF, ignore_index=True)
            name = optionsDict['DOSE'] + '_' + optionsDict['TREAT RATIO'] + '_' + \
                   optionsDict['CAR AFFINITY'] + '_' + optionsDict['ANTIGENS CANCER'] + '_'
            pickle.dump(setDF, open(args.saveLoc + args.xmlName + '_' + name + 'ANALYZED.pkl', "wb"))