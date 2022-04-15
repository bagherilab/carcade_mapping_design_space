import scripts.plot.plot_utilities
from math import exp
import matplotlib.pyplot as plt

def plot_binding_heuristic_CAR(SAVELOC):
    """Plot CAR-antigen binding heuristic."""

    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    affinityColorDict = scripts.plot.plot_utilities.make_car_affinity_color_dict()

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
    """Plot PD1-PDL1 binding heuristic."""

    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

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