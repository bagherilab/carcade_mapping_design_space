import scripts.stats.stats_utilities
import json
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols

def make_empty_features_dict():
    """Initialize empty dictionary of feature counts."""

    doseDict = {
        "0": 0,
        "250": 0,
        "500": 0,
        "1000": 0,
        "5000": 0,
        "10000": 0
    }

    trDict = {
        "NA": 0,
        "0:100": 0,
        "10:90": 0,
        "25:75": 0,
        "50:50": 0,
        "75:25": 0,
        "90:10": 0,
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

def make_anova_df():
    """Initialize dataframe for ANOVA and score analysis."""

    simsDFanova = pd.DataFrame({'DOSE': pd.Series([], dtype='float'),
                                'TREAT_RATIO': pd.Series([], dtype='float'),
                                'CAR_AFFINITY': pd.Series([], dtype='float'),
                                'ANTIGENS_CANCER': pd.Series([], dtype='float'),
                                'ANTIGENS_HEALTHY': pd.Series([], dtype='float'),
                                'Y_NORM_CANCER_LIVE': pd.Series([], dtype='float'),
                                'Y_NORM_HEALTHY_LIVE': pd.Series([], dtype='float'),
                                'Y_NORM_TCELL_LIVE': pd.Series([], dtype='float')})

    return simsDFanova

def populate_anova_df(simsDF, simsDFanova, FILEID):
    """Populate ANOVA dataframe with analyzed data."""

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

    return simsDFanova

def anova(simsDF, FILEID, NORM, SCORE, SAVELOC):
    """Conduct ANOVA statistical analysis."""

    if '_CH_' not in FILEID:
        add = ''

        #print('ANOVA FOR CANCER CELLS.')
        # model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + '
                           'C(DOSE):C(TREAT_RATIO) + C(DOSE):C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) +'
                           'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) +'
                           'C(CAR_AFFINITY):C(ANTIGENS_CANCER)', data=simsDF).fit()
        anova_CANCER = sm.stats.anova_lm(model_CANCER, typ=2)
        #print(anova_CANCER)

        f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'ANOVA_CANCER.xlsx'
        anova_CANCER.to_excel(f, header=True)

        #print('ANOVA FOR T-CELLS.')
        # model_TCELL = ols('Y_NORM_TCELL_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_TCELL = ols('Y_NORM_TCELL_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + '
                          'C(DOSE):C(TREAT_RATIO) + C(DOSE):C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) +'
                          'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) +'
                          'C(CAR_AFFINITY):C(ANTIGENS_CANCER)', data=simsDF).fit()
        anova_TCELL = sm.stats.anova_lm(model_TCELL, typ=2)
        #print(anova_TCELL)

        f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'ANOVA_TCELL.xlsx'
        anova_TCELL.to_excel(f, header=True)

    else:
        add = SCORE + '_'

        #print('ANOVA FOR CANCER CELLS.')
        #model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + C(ANTIGENS_HEALTHY) + '
                           'C(DOSE):C(TREAT_RATIO) + C(DOSE):C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) + C(DOSE):C(ANTIGENS_HEALTHY) +'
                           'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) + C(TREAT_RATIO):C(ANTIGENS_HEALTHY) +'
                           'C(CAR_AFFINITY):C(ANTIGENS_CANCER) + C(CAR_AFFINITY):C(ANTIGENS_HEALTHY) + '
                           'C(ANTIGENS_CANCER):C(ANTIGENS_HEALTHY)', data=simsDF).fit()
        anova_CANCER = sm.stats.anova_lm(model_CANCER, typ=2)
        #print(anova_CANCER)

        f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'ANOVA_CANCER.xlsx'
        anova_CANCER.to_excel(f, header=True)

        #print('ANOVA FOR HEALTHY CELLS.')
        # model_CANCER = ols('Y_NORM_CANCER_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_HEALTHY = ols('Y_NORM_HEALTHY_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + C(ANTIGENS_HEALTHY) + '
                            'C(DOSE):C(TREAT_RATIO) + C(DOSE):C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) + C(DOSE):C(ANTIGENS_HEALTHY) +'
                            'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) + C(TREAT_RATIO):C(ANTIGENS_HEALTHY) +'
                            'C(CAR_AFFINITY):C(ANTIGENS_CANCER) + C(CAR_AFFINITY):C(ANTIGENS_HEALTHY) +'
                            'C(ANTIGENS_CANCER):C(ANTIGENS_HEALTHY)', data=simsDF).fit()
        anova_HEALTHY = sm.stats.anova_lm(model_HEALTHY, typ=2)
        #print(anova_HEALTHY)

        f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'ANOVA_HEALTHY.xlsx'
        anova_HEALTHY.to_excel(f, header=True)

        #print('ANOVA FOR T-CELLS.')
        # model_TCELL = ols('Y_NORM_TCELL_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_TCELL = ols('Y_NORM_TCELL_LIVE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + C(ANTIGENS_HEALTHY) + '
                          'C(DOSE):C(TREAT_RATIO) + C(DOSE):C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) + C(DOSE):C(ANTIGENS_HEALTHY) +'
                          'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) + C(TREAT_RATIO):C(ANTIGENS_HEALTHY) +'
                          'C(CAR_AFFINITY):C(ANTIGENS_CANCER) + C(CAR_AFFINITY):C(ANTIGENS_HEALTHY) +'
                          'C(ANTIGENS_CANCER):C(ANTIGENS_HEALTHY)', data=simsDF).fit()
        anova_TCELL = sm.stats.anova_lm(model_TCELL, typ=2)
        #print(anova_TCELL)

        f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'ANOVA_TCELL.xlsx'
        anova_TCELL.to_excel(f, header=True)

        #print('ANOVA FOR SCORE.')
        # model_CANCER = ols('SCORE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER)', data=simsDF).fit()
        model_SCORE = ols('SCORE ~ C(DOSE) + C(TREAT_RATIO) + C(CAR_AFFINITY) + C(ANTIGENS_CANCER) + C(ANTIGENS_HEALTHY) + '
                          'C(DOSE):C(TREAT_RATIO) + C(DOSE):C(CAR_AFFINITY) + C(DOSE):C(ANTIGENS_CANCER) + C(DOSE):C(ANTIGENS_HEALTHY) +'
                          'C(TREAT_RATIO):C(CAR_AFFINITY) + C(TREAT_RATIO):C(ANTIGENS_CANCER) + C(TREAT_RATIO):C(ANTIGENS_HEALTHY) +'
                          'C(CAR_AFFINITY):C(ANTIGENS_CANCER) + C(CAR_AFFINITY):C(ANTIGENS_HEALTHY) +'
                          'C(ANTIGENS_CANCER):C(ANTIGENS_HEALTHY)', data=simsDF).fit()
        anova_SCORE = sm.stats.anova_lm(model_SCORE, typ=2)
        #print(anova_SCORE)

        f = SAVELOC + FILEID + '_' + NORM + '_' + add + 'ANOVA_SCORE.xlsx'
        anova_SCORE.to_excel(f, header=True)

    return

def feature_analysis(simsDF, RANK, FILEID, NORM, SCORE, SAVELOC):
    """Conduct analysis to count number of each features value across simulations that meet a specified threshold."""

    SCORE_MIN_HEALTHY_THRESHOLD = scripts.stats.stats_utilities.define_score_min_healthy_threshold()
    TREAT_RATIO_DICT_REVERSE = scripts.stats.stats_utilities.make_treat_ratio_reverse_key_dict_stats()

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
                elif key == '0.1': key = '0.10'
                elif key == '0.9': key = '0.90'

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
    """Save sorted dataframe with analyzed simulations that meet a desired threshold sorted by desired output as excel file."""

    SCORE_MIN_HEALTHY_THRESHOLD = scripts.stats.stats_utilities.define_score_min_healthy_threshold()

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
        simsDF = simsDF[simsDF['Y_NORM_HEALTHY_LIVE'] >= SCORE_MIN_HEALTHY_THRESHOLD]
        simsDF = simsDF[simsDF['Y_NORM_CANCER_LIVE'] < 1]

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
    """Save sorted dataframe with all analyzed simulations sorted by desired output regardless of desired threshold as excel file."""

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