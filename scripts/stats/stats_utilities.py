import scripts.plot.plot_utilities
from matplotlib import colors as mcolors

def define_axes_sets_dict():
    """Define axes sets dictionary."""

    AXES_SETS = {
        "DOSE": [250, 500, 1000, 5000, 10000],
        "TREAT RATIO": ["0:100", "10:90", "25:75", "50:50", "75:25", "90:10", "100:0"],
        "CAR AFFINITY": [1e-6, 1e-7, 1e-8, 1e-9],
        "ANTIGENS CANCER": [100, 500, 1000, 5000, 10000],
        "ANTIGENS HEALTHY": [0, 100]
    }

    return AXES_SETS

def define_constants_dict():
    """Define dictionary of intermediate values of each feature held constant in some analyses."""

    CONSTANTS = {
        "DOSE": 500,
        "TREAT RATIO": 0.5,
        "CAR AFFINITY": 1e-7,
        "ANTIGENS CANCER": 1000,
        "ANTIGENS HEALTHY": [0, 100]
    }

    return CONSTANTS

def define_score_min_healthy_threshold():
    """Define minimum normalized healthy threshold that make treatment viable option."""

    SCORE_MIN_HEALTHY_THRESHOLD = 0.5

    return SCORE_MIN_HEALTHY_THRESHOLD

def define_color_reference_list(FILEID):
    """Define list of features to reference color maps for based on simuatlion type."""

    # Create color reference lists
    colorlist = ['DOSE', 'TREAT_RATIO', 'CAR_AFFINITY', 'ANTIGENS_CANCER']
    if '_CH_' in FILEID: colorlist.append('ANTIGENS_HEALTHY')

    return colorlist

def NonLinCdict(steps, col_array):
    """Create non-linear color dictionary."""

    cdict = {'red': (), 'green': (), 'blue': ()}

    for s, col in zip(steps, col_array):
        rgb = mcolors.to_rgb(col)
        cdict['red'] = cdict['red'] + ((s, rgb[0], rgb[0]),)
        cdict['green'] = cdict['green'] + ((s, rgb[1], rgb[1]),)
        cdict['blue'] = cdict['blue'] + ((s, rgb[2], rgb[2]),)

    return cdict

def maximum_absolute_scaling(df):
    """Calcualte minimum absolute scaling of given dataframe."""

    # Copy the dataframe
    df_scaled = df.copy()

    # Apply maximum absolute scaling
    for column in df_scaled.columns:
        if 'NORM' not in column and 'SCORE' not in column:
            if df_scaled[column].abs().max() != 0.0:
                df_scaled[column] = df_scaled[column] / df_scaled[column].abs().max()

    return df_scaled

def make_dose_color_dict():

    DCOLORS = scripts.plot.plot_utilities.make_dose_color_scale()

    doseColorDict = {
        "0": 'black',
        "250": DCOLORS[0],
        "500": DCOLORS[1],
        "1000": DCOLORS[2]
    }
    return doseColorDict

def make_treat_ratio_heatmap_color_dict():
    """Make dictionary of feature value to color value for treat ratio feature for use in heatmap creation."""

    TRCOLORS = scripts.plot.plot_utilities.make_treat_ratio_color_scale()

    trColorDictHeatMap = {
        "0.0": TRCOLORS[0],
        "0.10": TRCOLORS[1],
        "0.25": TRCOLORS[2],
        "0.50": TRCOLORS[3],
        "0.75": TRCOLORS[4],
        "0.90": TRCOLORS[5],
        "1.0": TRCOLORS[6],
    }

    return trColorDictHeatMap

def make_affinity_heatmap_color_dcit():

    ACOLORS = scripts.plot.plot_utilities.make_car_affinity_color_scale()

    affinityColorDictHeatMap = {
        "NA": 'black',
        "0.0": ACOLORS[4],
        "1e-10": ACOLORS[4],
        "1e-09": ACOLORS[3],
        "1e-08": ACOLORS[2],
        "1e-07": ACOLORS[1],
        "1e-06": ACOLORS[0],
    }
    return affinityColorDictHeatMap

def make_features_heatmap_color_dict():
    """Make dictionary of features to feature-color dictionaries with treat ratio heatmap dictionary for use in heatmap creation."""

    doseColorDict = make_dose_color_dict()
    trColorDictHeatMap = make_treat_ratio_heatmap_color_dict()
    affinityColorDictHeatMap = make_affinity_heatmap_color_dcit()
    acColorDict = scripts.plot.plot_utilities.make_antigens_cancer_color_dict()
    ahColorDict = scripts.plot.plot_utilities.make_antigens_healthy_color_dict()

    COLOR_DICT_HEATMAP = {
        "DOSE": doseColorDict,
        "TREAT RATIO": trColorDictHeatMap,
        "CAR AFFINITY": affinityColorDictHeatMap,
        "ANTIGENS CANCER": acColorDict,
        "ANTIGENS HEALTHY": ahColorDict
    }

    return COLOR_DICT_HEATMAP

def make_features_heatmap_color_dict_extended_dose():
    """Make dictionary of features to feature-color dictionaries with treat ratio heatmap dictionary and extended dose values for use in heatmap creation."""

    doseColorDict = scripts.plot.plot_utilities.make_dose_color_dict()
    trColorDictHeatMap = make_treat_ratio_heatmap_color_dict()
    affinityColorDictHeatMap = make_affinity_heatmap_color_dcit()
    acColorDict = scripts.plot.plot_utilities.make_antigens_cancer_color_dict()
    ahColorDict = scripts.plot.plot_utilities.make_antigens_healthy_color_dict()

    COLOR_DICT_HEATMAP = {
        "DOSE": doseColorDict,
        "TREAT RATIO": trColorDictHeatMap,
        "CAR AFFINITY": affinityColorDictHeatMap,
        "ANTIGENS CANCER": acColorDict,
        "ANTIGENS HEALTHY": ahColorDict
    }

    return COLOR_DICT_HEATMAP

def make_treat_ratio_reverse_key_dict_stats():
    """Make dictionary of CD4+ fraction values to treat ratio value."""

    TREAT_RATIO_DICT_REVERSE = {'0.0': '0:100', '0.1': '10:90', '0.25': '25:75', '0.50': '50:50', '0.75': '75:25',
                                '0.9': '90:10', '1.0': '100:0'}

    return TREAT_RATIO_DICT_REVERSE