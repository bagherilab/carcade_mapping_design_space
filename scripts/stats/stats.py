import scripts.analyze.analyze_utilities
import scripts.plot.plot_utilities
import scripts.stats.stats_utilities
import scripts.stats.stats_heatmaps
import scripts.stats.stats_features
import pickle
import pandas as pd
from itertools import combinations
import matplotlib.pyplot as plt

def get_file_id(fileName):
    """Get file ID based on file name."""

    FILEID = ''

    if 'ANALYZED' in fileName:
        FILEID = fileName.replace('_ANALYZED.pkl', '')

    if 'ENVIRONMENT' in fileName:
        FILEID = fileName.replace('_ENVIRONMENT.pkl', '')

    if 'SPATIAL' in fileName:
        FILEID = fileName.replace('_SPATIAL.pkl', '')

    if 'LYSED' in fileName:
        FILEID = fileName.replace('_LYSED.pkl', '')

    return FILEID

def check_norm_arg(norm):
    """Check normalize argument for how to normalize data."""

    if norm == 'INIT':
        NORM = 'INIT'
    elif norm == 'UNTREATED':
        NORM = 'UNTREATED'
    else:
        NORM = 'INIT'

    return NORM

def check_score_arg(score):
    """Check score argument for how to calculate score."""

    if score == 'SUM':
        SCORE = 'SUM'
    elif score == 'RATIO':
        SCORE = 'RATIO'
    else:
        SCORE = 'SUM'

    return SCORE

def clean_data(simsDF):
    """Clean up dataframe, add rows for nomalized cell populations, and separate out untreated simulations."""

    print('\t\t' + 'Cleaning up data...')
    TREAT_RATIO_DICT = scripts.plot.plot_utilities.make_treat_ratio_key_dict()
    simsDF = simsDF.replace({'TREAT RATIO': TREAT_RATIO_DICT})

    simsDF = simsDF.rename(columns={'TREAT RATIO': 'TREAT_RATIO', 'CAR AFFINITY': 'CAR_AFFINITY',
                                    'ANTIGENS CANCER': 'ANTIGENS_CANCER', 'ANTIGENS HEALTHY': 'ANTIGENS_HEALTHY'})
    simsDF['Y_NORM_CANCER_LIVE'] = None
    simsDF['Y_NORM_HEALTHY_LIVE'] = None
    simsDF['Y_NORM_TCELL_LIVE'] = None

    untreatedDF = simsDF[simsDF['DOSE'] == 0]

    simsDF = simsDF[simsDF['DOSE'] != 0]

    return simsDF, untreatedDF

def normalize_data_by_initial(simsDF, FILEID):
    """Normalize data by initial cell count at start of treatment."""

    if 'VITRO' in FILEID:
        for i in range(0, len(simsDF)):
            simsDF.at[i, 'Y_NORM_CANCER_LIVE'] = float(
                simsDF.iloc[i]['CANCER LIVE'][-1] / simsDF.iloc[i]['CANCER LIVE'][0])
            simsDF.at[i, 'Y_NORM_TCELL_LIVE'] = float(simsDF.iloc[i]['T-CELL LIVE'][-1] / simsDF.iloc[i]['DOSE'])
            if '_CH_' in FILEID:
                simsDF.at[i, 'Y_NORM_HEALTHY_LIVE'] = float(
                    simsDF.iloc[i]['HEALTHY LIVE'][-1] / simsDF.iloc[i]['HEALTHY LIVE'][0])

    else:
        for i in range(0, len(simsDF)):
            simsDF.at[i, 'Y_NORM_CANCER_LIVE'] = float(
                simsDF.iloc[i]['CANCER LIVE'][-1] / simsDF.iloc[i]['CANCER LIVE'][44])
            simsDF.at[i, 'Y_NORM_TCELL_LIVE'] = float(simsDF.iloc[i]['T-CELL LIVE'][-1] / simsDF.iloc[i]['DOSE'])
            if '_CH_' in FILEID:
                simsDF.at[i, 'Y_NORM_HEALTHY_LIVE'] = float(
                    simsDF.iloc[i]['HEALTHY LIVE'][-1] / simsDF.iloc[i]['HEALTHY LIVE'][44])

    return simsDF

def normalize_data_by_untreated(simsDF, untreatedDF, FILEID):
    """Normalize data by number of cells in corresponding untreated seed."""

    if 'VITRO' in FILEID:
        for i in range(0, len(simsDF)):

            for j in range(0, len(untreatedDF)):
                if untreatedDF.iloc[j]['SEED'] == simsDF.iloc[i]['SEED']:
                    cancerUntreated = untreatedDF.iloc[j]['CANCER LIVE'][-1]
                    if '_CH_' in FILEID:
                        healthyUntreated = untreatedDF.iloc[j]['HEALTHY LIVE'][-1]
                    break

            simsDF.at[i, 'Y_NORM_CANCER_LIVE'] = float(simsDF.iloc[i]['CANCER LIVE'][-1] / cancerUntreated)
            simsDF.at[i, 'Y_NORM_TCELL_LIVE'] = float(simsDF.iloc[i]['T-CELL LIVE'][-1] / simsDF.iloc[i]['DOSE'])
            if '_CH_' in FILEID:
                simsDF.at[i, 'Y_NORM_HEALTHY_LIVE'] = float(simsDF.iloc[i]['HEALTHY LIVE'][-1] / healthyUntreated)

    else:
        for i in range(0, len(simsDF)):

            for j in range(0, len(untreatedDF)):
                if untreatedDF.iloc[j]['SEED'] == simsDF.iloc[i]['SEED']:
                    cancerUntreated = untreatedDF.iloc[j]['CANCER LIVE'][-1]
                    if '_CH_' in FILEID:
                        healthyUntreated = untreatedDF.iloc[j]['HEALTHY LIVE'][-1]
                    break

            simsDF.at[i, 'Y_NORM_CANCER_LIVE'] = float(simsDF.iloc[i]['CANCER LIVE'][-1] / cancerUntreated)
            simsDF.at[i, 'Y_NORM_TCELL_LIVE'] = float(simsDF.iloc[i]['T-CELL LIVE'][-1] / simsDF.iloc[i]['DOSE'])
            if '_CH_' in FILEID:
                simsDF.at[i, 'Y_NORM_HEALTHY_LIVE'] = float(simsDF.iloc[i]['HEALTHY LIVE'][-1] / healthyUntreated)

    return simsDF

def normalize_data(simsDF, untreatedDF, NORM, FILEID):
    """Call approrpiate normalize data function based on type of normalization requested."""

    # Normalize data by INITIALIZATION quantity
    if NORM == 'INIT':

        simsDF = normalize_data_by_initial(simsDF, FILEID)

    # Normalize data by UNTREATED quantity
    else:

        simsDF = normalize_data_by_untreated(simsDF, untreatedDF, FILEID)

    return simsDF

def create_feature_and_response_combo_lists(FILEID):
    """Create lists with combinations of all features and list of possible responses based on simualtion type."""

    # Create list of combinations of features and list of responses
    if '_C_' in FILEID:
        comb = list(combinations(['DOSE', 'TREAT RATIO', 'CAR AFFINITY', 'ANTIGENS CANCER'], 2))
        response_list = ['Y_NORM_CANCER_LIVE', 'Y_NORM_TCELL_LIVE']
    else:
        comb = list(combinations(['DOSE', 'TREAT RATIO', 'CAR AFFINITY', 'ANTIGENS CANCER', 'ANTIGENS HEALTHY'], 2))
        response_list = ['Y_NORM_CANCER_LIVE', 'Y_NORM_HEALTHY_LIVE', 'Y_NORM_TCELL_LIVE', 'SCORE']

    return comb, response_list

def calculate_score_sum(simsDF, i, cancer, healthy, FILEID):
    """Calculate score as subtraction of normalized healthy and cancer counts."""

    if 'VITRO' in FILEID:
        norm = float(simsDF.iloc[i]['CANCER LIVE'][0] / simsDF.iloc[i]['HEALTHY LIVE'][0])

    else:
        norm = float(simsDF.iloc[i]['CANCER LIVE'][44] / simsDF.iloc[i]['HEALTHY LIVE'][44])

    score = (healthy * norm) - cancer

    return score

def calculate_score(simsDF, simsDFanova, SCORE, FILEID):
    """Calculate score based on score type selected."""

    scorelist = []
    for i in range(0, len(simsDFanova)):
        cancer = simsDFanova.iloc[i]['Y_NORM_CANCER_LIVE']
        healthy = simsDFanova.iloc[i]['Y_NORM_HEALTHY_LIVE']

        if SCORE == 'SUM':
            score = calculate_score_sum(simsDF, i, cancer, healthy, FILEID)

        else:
            score = None

        scorelist.append(score)
    simsDFanova['SCORE'] = scorelist

    return simsDFanova

def plot_all_output_heatmaps_lineplot(simsDFanova, NORM, SCORE, FILEID, SAVELOC):
    """Call plotters for heatmaps with all outputs showing as line plots and sorted by score."""

    comb, response_list = create_feature_and_response_combo_lists(FILEID)

    scripts.stats.stats_heatmaps.plot_output_heatmap_with_subplots_line_multiple_outputs(simsDFanova, response_list, True, FILEID, NORM, SCORE, SAVELOC)

    plt.close('all')

    return

def create_csv_files_with_simulations_sorted_by_response(simsDFanova, NORM, SCORE, FILEID, SAVELOC):
    """Call printers to make excel files with data analyzed by feature and score values."""

    comb, response_list = create_feature_and_response_combo_lists(FILEID)

    # Create excel files sorted by output
    print('\t\tCreating Excel files for all threshold passing simulations.')
    for response in response_list:
        if response != 'Y_NORM_TCELL_LIVE':
            scripts.stats.stats_features.feature_analysis(simsDFanova, response, FILEID, NORM, SCORE, SAVELOC)
            scripts.stats.stats_features.save_sorted_df_csv(simsDFanova, response, FILEID, NORM, SCORE, SAVELOC)

    return

def plot_heatmaps_and_make_csv_files(simsDFanova, NORM, SCORE, FILEID, SAVELOC):
    """Call functions that will make heatmaps, elbow plots, and make excel files."""

    # Plot output heatmaps
    plot_all_output_heatmaps_lineplot(simsDFanova, NORM, SCORE, FILEID, SAVELOC)

    # Create excel files sorted by output
    create_csv_files_with_simulations_sorted_by_response(simsDFanova, NORM, SCORE, FILEID, SAVELOC)

    return

def average_conditions(simsDF, FILEID):
    """Average conditions across replicates (seeds)."""

    AXES_SETS = scripts.stats.stats_utilities.define_axes_sets_dict()
    TREAT_RATIO_DICT = scripts.plot.plot_utilities.make_treat_ratio_key_dict()
    TREAT_RATIO_DICT_REVERSE = scripts.stats.stats_utilities.make_treat_ratio_reverse_key_dict_stats()

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
        elif treatratio == '0.1': treatratio = '0.10'
        elif treatratio == '0.9': treatratio = '0.90'

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

def conduct_stats_analyses_all_data(simsDFanova, NORM, SCORE, FILEID, SAVELOC):
    """Conduct stats analysis on all data (not averaged)."""

    plot_heatmaps_and_make_csv_files(simsDFanova, NORM, SCORE, FILEID, SAVELOC)

    return

def conduct_stats_analysis_averaged_data(simsDFanova, NORM, SCORE, FILEID, SAVELOC):
    """Conduct analysis on averaged data."""

    # Average data across conditions
    print('\t\tAveraging datta across conditions.')
    simsDFavg = average_conditions(simsDFanova, FILEID)
    FILEIDAVG = FILEID + '_AVG'

    plot_heatmaps_and_make_csv_files(simsDFavg, NORM, SCORE, FILEIDAVG, SAVELOC)

    return

def stats_data(simsDF, NORM, SCORE, AVG, FILEID, SAVELOC):
    """Run stats analysis on given file."""

    simsDF, untreatedDF = clean_data(simsDF)

    # Create normalized data based on NORM
    simsDF = normalize_data(simsDF, untreatedDF, NORM, FILEID)

    # Create anova dataframe for use in rest of code (not just anova functions)
    simsDFanova = scripts.stats.stats_features.make_anova_df()
    simsDFanova = scripts.stats.stats_features.populate_anova_df(simsDF, simsDFanova, FILEID)

    # Create SCORE if both cancer and healthy cells present
    if '_CH_' in FILEID:
        simsDFanova = calculate_score(simsDF, simsDFanova, SCORE, FILEID)

    # Conduct analyses for data that is not averaged across replicates
    if not AVG:

        conduct_stats_analyses_all_data(simsDFanova, NORM, SCORE, FILEID, SAVELOC)

    # Conduct analyses for data that is averaged across replicates
    else:
        conduct_stats_analysis_averaged_data(simsDFanova, NORM, SCORE, FILEID, SAVELOC)

    return

def stats(files, saveLoc, norm='INIT', score='SUM', average=False):
    """Run stats analysis on all given files.

    stats.py takes a directory of (or a single) .pkl simulation files that result from analyze_cells.py and
    does statistics on the DATA features in the dataframe over time.

    Usage:
        stats(files, saveLoc, norm='INIT', score='SUM', average=False)

        files
            Path to .pkl or directory.
        savLoc
            Location of where to save file, default will save here.
        [norm]
            Normalize final values by initial (INIT) or untreated (UNTREATED) (default: INIT).
        [score]
            Dictate which score type to use (default and only current option: SUM).
        [average]
            Average across seed replicates (default: False).
    """

    # Get files
    PKLFILES = scripts.analyze.analyze_utilities.get_pkl_files(files)

    print("Running stats for the following files:")

    for file in PKLFILES:

        fileName = scripts.plot.plot_utilities.get_file_name(file)

        print('\t' + fileName)

        FILEID = get_file_id(fileName)

        with open(file, 'rb') as f:
            simsDF = pickle.load(f)

        NORM = check_norm_arg(norm)
        SCORE = check_score_arg(score)

        stats_data(simsDF, NORM, SCORE, average, FILEID, saveLoc)

    return