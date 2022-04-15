import re
from matplotlib import cm as mplcm
from colour import Color

def define_fileid_split_indices_dict():
    """Make dictionary of file split information contained at each index after splitting based on _."""

    FILEID_SPLIT_INDICES = {'EXP': 0, 'PLATE': 1, 'TREAT': 2, 'POPS': 3, 'DIMENSION': 4, 'DOSE': 5, 'TREAT RATIO': 6,
                            'CAR AFFINITY': 7, 'ANTIGENS CANCER': 8, 'ANTIGENS HEALTHY': 9}

    return FILEID_SPLIT_INDICES

def define_pop_frac_names_list():
    """Define list of cell populations with fractions calculated over time."""

    POP_FRAC_NAMES = ['CANCER LIVE', 'HEALTHY LIVE']

    return POP_FRAC_NAMES

def define_pop_lysis_names_list():
    """Define list of cell populations lysed over time."""

    POP_LYSIS_NAMES = ['TISSUE', 'CANCER', 'HEALTHY']

    return POP_LYSIS_NAMES

def define_pop_codes_dict():
    """Define population codes dictionary based on cell codes for each population in json."""

    POP_CODES = {'CANCER': 0, 'HEALTHY': 1, 'CD4 T-CELL': 2, 'CD8 T-CELL': 3}

    return POP_CODES

def define_mol_names_list():
    """Define list of environment species names."""

    MOL_NAMES = ['GLUCOSE', 'OXYGEN', 'TGFA', 'IL-2']

    return MOL_NAMES

def define_mol_concentrations_units_dict():
    """Define dictionary of environment species concentration/value units."""

    MOL_CONC_UNITS = {'GLUCOSE': 'fmol/ml', 'OXYGEN': 'mmHg', 'TGFA': 'pg/ml', 'IL-2': 'pg/ml'}

    return MOL_CONC_UNITS

def define_dish_times_list():
    """Define list of times to make some plots at for dish simulations."""

    DISH_TIMES = [0, 1, 4, 7]

    return DISH_TIMES

def define_tissue_times_list():
    """Define list of times to make some plots at for tissue simulations."""

    TISSUE_TIMES = [22, 26, 29, 31]

    return TISSUE_TIMES

def make_dose_color_scale():
    """Make color scale for dose feature."""

    # Set DOSE color scale
    DCOLORS = []
    dcmap = mplcm.get_cmap('Greens')
    dvalues = [0.33, 0.66, 1.0]
    for v in dvalues:
        DCOLORS.append(dcmap(v))

    return DCOLORS

def make_treat_ratio_color_scale():
    """Make color scale for treat ratio feature."""

    TRCOLORS = []
    trcmap = mplcm.get_cmap('Reds')
    trvalues = [0.2, 0.3, 0.4, 0.6, 0.8, 0.9, 1.0]  # CD4%: 0, 10, 25, 50, 75, 90, 100
    for v in trvalues:
        TRCOLORS.append(trcmap(v))

    return TRCOLORS

def make_car_affinity_color_scale():
    """Make color scale for CAR affinity feature."""

    # Set AFFINITY color scale
    ACOLORS = []
    acmap = mplcm.get_cmap('Oranges')
    avalues = [0.2, 0.4, 0.6, 0.8, 1.0]
    for v in avalues:
        ACOLORS.append(acmap(v))

    return ACOLORS

def make_antigens_cancer_color_scale():
    """Make color scale for antigens per cancer cell feature."""

    # Set ANTIGENS CANCER color scale
    ACCOLORS = []
    accmap = mplcm.get_cmap('Blues')
    acvalues = [0.2, 0.4, 0.6, 0.8, 1.0]
    for v in acvalues:
        ACCOLORS.append(accmap(v))

    return ACCOLORS

def make_antigens_healthy_color_scale():
    """Make color scale for antigens per healthy cell feature."""

    # Set ANTIGENS HEALTHY color scale
    AHCOLORS = []
    ahcmap = mplcm.get_cmap('Purples')
    ahvalues = [0.5, 1.0]
    for v in ahvalues:
        AHCOLORS.append(ahcmap(v))

    return AHCOLORS

def make_dose_color_dict():
    """Make dictionary of feature value to color value for dose feature."""

    DCOLORS = make_dose_color_scale()

    doseColorDict = {
        "0": 'black',
        "250": DCOLORS[0],
        "500": DCOLORS[1],
        "1000": DCOLORS[2],
        "5000": Color('#142516').get_hex(),
        "10000": Color('#0b150c').get_hex(),
    }

    return doseColorDict

def make_treat_ratio_color_dict():
    """Make dictionary of feature value to color value for treat ratio feature."""

    TRCOLORS = make_treat_ratio_color_scale()

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
    return trColorDict

def make_car_affinity_color_dict():
    """Make dictionary of feature value to color value for CAR affinity feature."""

    ACOLORS = make_car_affinity_color_scale()

    affinityColorDict = {
        "NA": 'black',
        "0.0": 'black',
        "1e-06": ACOLORS[0],
        "1e-07": ACOLORS[1],
        "1e-08": ACOLORS[2],
        "1e-09": ACOLORS[3],
        "1e-10": ACOLORS[4]
    }
    return affinityColorDict

def make_antigens_cancer_color_dict():
    """Make dictionary of feature value to color value for antigens per cancer cell feature."""

    ACCOLORS = make_antigens_cancer_color_scale()

    acColorDict = {
        "0": 'black',
        "100": ACCOLORS[0],
        "500": ACCOLORS[1],
        "1000": ACCOLORS[2],
        "5000": ACCOLORS[3],
        "10000": ACCOLORS[4]
    }

    return acColorDict

def make_antigens_healthy_color_dict():
    """Make dictionary of feature value to color value for antigens per healthy cell feature."""

    AHCOLORS = make_antigens_healthy_color_scale()

    ahColorDict = {
        "CONTROL": "black",
        "0": AHCOLORS[0],
        "100": AHCOLORS[1]
    }

    return ahColorDict

def make_features_color_dict():
    """Make dictionary of features to feature-color dictionaries."""

    doseColorDict = make_dose_color_dict()
    trColorDict = make_treat_ratio_color_dict()
    affinityColorDict = make_car_affinity_color_dict()
    acColorDict = make_antigens_cancer_color_dict()
    ahColorDict = make_antigens_healthy_color_dict()

    COLOR_DICT = {
        "DOSE": doseColorDict,
        "TREAT RATIO": trColorDict,
        "CAR AFFINITY": affinityColorDict,
        "ANTIGENS CANCER": acColorDict,
        "ANTIGENS HEALTHY": ahColorDict
    }

    return COLOR_DICT

def make_dose_line_dict():
    """Make dictionary of feature value to linestyle for dose feature."""

    doseLineDict = {
        "0": 'solid',
        "250": "dotted",
        "500": "solid",
        "1000": "dashed",
        "5000": 'dashdot',
        "10000": 'dashdot'
    }

    return doseLineDict

def make_live_line_dict():
    """Make dictionary of cell states included to linestyle."""

    liveLineDict = {
        "LIVE": 'dashed',
        "TOTAL": 'solid'
    }

    return liveLineDict

def make_lysed_line_dict():
    """Make dictionary of lysed cell population to linestyle."""

    lysedLineDict = {
        'CANCER': 'solid',
        'HEALTHY': 'dashed'
    }

    return lysedLineDict

def make_dose_marker_dict():
    """Make dictionary of feature value to marker style for dose feature."""

    doseMarkerDict = {
        "0": 'X',
        "250": 's',
        "500": '^',
        "1000": 'o'
    }

    return doseMarkerDict

def make_treat_ratio_marker_dict():
    """Make dictionary of feature value to marker style for treat ratio feature."""

    trMarkerDict = {
        "NA": 'X',
        "0:100": 's',
        "25:75": '^',
        "50:50": 'o',
        "75:25": 'v',
        "100:0": 'D',
    }

    return trMarkerDict

def make_car_affinity_marker_dict():
    """Make dictionary of feature value to marker style for CAR affinity feature."""

    affinityMarkerDict = {
        "NA": 'X',
        "0.0": 'X',
        "1e-06": 's',
        "1e-07": '^',
        "1e-08": 'o',
        "1e-09": 'v',
        "1e-10": 'D'
    }

    return affinityMarkerDict

def make_antigens_cancer_maker_dict():
    """Make dictionary of feature value to marker style for antigens per cancer cell feature."""

    acMarkerDict = {
        "0": 'X',
        "100": 's',
        "500": '^',
        "1000": 'o',
        "5000": 'v',
        "10000": 'D',
    }

    return acMarkerDict

def make_antigens_healthy_marker_dict():
    """Make dictionary of feature value to marker style for antigens per healthy cell feature."""

    ahMarkerDict = {
        "0": 'X',
        "100": 'o'
    }

    return ahMarkerDict

def make_no_marker_dict():
    """Make dictionary of marker style for no marker selected."""

    noMarkerDict = {
        'o': 'o'
    }
    return noMarkerDict

def make_features_marker_dict():
    """Make dictionary of features to feature-markerstyle dictionaries."""

    noMarkerDict = make_no_marker_dict()
    doseMarkerDict = make_dose_marker_dict()
    trMarkerDict = make_treat_ratio_marker_dict()
    affinityMarkerDict = make_car_affinity_marker_dict()
    acMarkerDict = make_antigens_cancer_maker_dict()
    ahMarkerDict = make_antigens_healthy_marker_dict()

    MARKER_DICT = {
        'o': noMarkerDict,
        "DOSE": doseMarkerDict,
        "TREAT RATIO": trMarkerDict,
        "CAR AFFINITY": affinityMarkerDict,
        "ANTIGENS CANCER": acMarkerDict,
        "ANTIGENS HEALTHY": ahMarkerDict
    }

    return MARKER_DICT

def make_treat_ratio_key_dict():
    """Make dictionary of treat ratio value to CD4+ fraction values."""

    TREAT_RATIO_DICT = {'0:100': 0.0, '10:90': 0.1, '25:75': 0.25, '50:50': 0.5, '75:25': 0.75, '90:10': 0.9,
                        '100:0': 1.0}

    return TREAT_RATIO_DICT

def make_treat_ratio_reverse_key_dict():
    """Make dictionary of CD4+ fraction values to treat ratio value."""

    TREAT_RATIO_DICT_REVERSE = {'0.0': '0:100', '0.1': '10:90', '0.25': '25:75', '0.5': '50:50', '0.75': '75:25',
                                '0.9': '90:10', '1.0': '100:0'}

    return TREAT_RATIO_DICT_REVERSE

def make_literature_data_list():
    """Make dictionaries containing CAR and cell lysis information for all experimental literature publications."""

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
                    "Antigens": [1000, 50118.72, 158489.3, 501187.2, 1584893.2, 7943282.3],
                # THP-1, U937, TIME, MHH-CALL-4
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
                    "Kill % Err": [0.05, 0, 0.1, 0, 0, 0]
                },
                "KD": 1.76e-8
            },
        },
        "CITATION": "Hernandez-Lopez 2021 (HER2-41BB-CD3\u03B6, in vitro)",
        "SAVE": "HernandezLopez2021",
        "MARKER": "p"

    }

    litDict = [Arcangeli2017, Caruso2015, Chmielewski2004, Ghorashian2019, HernandezLopez2021, Liu2015, Wantanabe2015]

    return litDict

def define_plotting_globals():
    """Define plotting globals for plot stylistic features."""

    FIG_SIZE_X = 7
    FIG_SIZE_Y = 7

    TICKSIZE = 30
    FONTSIZE_AXES_VALUES = 20
    FONTSIZE_AXES_TITLES = 25
    LABELPAD = 10

    return FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD

def get_file_name(file):
    """Get plot file name based on name of file being plotted."""

    if 'VITRO' in file:
        fileName = re.sub('.*VITRO', 'VITRO', file)

    else:
        fileName = re.sub('.*VIVO', 'VIVO', file)

    return fileName

def count_x_in_file_name(filesplit):
    """Count the number of X's (or features) where all values included in subset based on file name."""

    # Count Xs in filesplit
    X = 0
    for i in range(5, 10):
        if filesplit[i] == 'X': X += 1

    return X