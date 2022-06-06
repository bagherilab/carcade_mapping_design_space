import scripts.plot.plot_utilities
import matplotlib.pyplot as plt

def get_max(LIST, INDEX):
    """Find max value in list of values."""

    values = [i[INDEX] for i in LIST]
    maxval = max(values)

    return maxval

def plot_kill_curve_normalized_sim(simsDF, FILEID, SAVELOC, TIME):
    """Plot simulated kill curve data where cells killed is calculated as 1-(living/living at treatment start) and values of cancer cells killed and antigen expression level are normalized to max value."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    affinityColorDict = scripts.plot.plot_utilities.make_car_affinity_color_dict()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    TIMES = simsDF.iloc[0]['TIME']
    index = -1
    for t in range(0,len(TIMES)):
        if float(TIMES[t]) == float(TIME):
            index = t

    figKC = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figKC.add_subplot(1, 1, 1)

    KD6D250 = []
    KD7D250 = []
    KD8D250 = []
    KD9D250 = []

    KD6D500 = []
    KD7D500 = []
    KD8D500 = []
    KD9D500 = []

    KD6D1000 = []
    KD7D1000 = []
    KD8D1000 = []
    KD9D1000 = []


    # Make ANTIGEN and % KILLING CANCER CELLS array sorted by dose
    for i in range(0, len(simsDF)):
        antigen = simsDF.iloc[i]['ANTIGENS CANCER']
        dose = simsDF.iloc[i]['DOSE']
        affinity = simsDF.iloc[i]['CAR AFFINITY']
        killed = 1-(simsDF.iloc[i]['CANCER LIVE'][index]/simsDF.iloc[i]['CANCER LIVE'][0])

        if dose == 250:
            if affinity >= 1e-6:
                KD6D250.append([antigen, killed])
            elif 1e-6 > affinity >= 1e-7:
                KD7D250.append([antigen, killed])
            elif 1e-7 > affinity >= 1e-8:
                KD8D250.append([antigen, killed])
            elif 1e-8 > affinity >= 1e-9:
                KD9D250.append([antigen, killed])
            elif 1e-9 > affinity >= 1e-10:
                continue
            else:
                continue

        elif dose == 500:
            if affinity >= 1e-6:
                KD6D500.append([antigen, killed])
            elif 1e-6 > affinity >= 1e-7:
                KD7D500.append([antigen, killed])
            elif 1e-7 > affinity >= 1e-8:
                KD8D500.append([antigen, killed])
            elif 1e-8 > affinity >= 1e-9:
                KD9D500.append([antigen, killed])
            elif 1e-9 > affinity >= 1e-10:
                continue
            else:
                continue

        elif dose == 1000:
            if affinity >= 1e-6:
                KD6D1000.append([antigen, killed])
            elif 1e-6 > affinity >= 1e-7:
                KD7D1000.append([antigen, killed])
            elif 1e-7 > affinity >= 1e-8:
                KD8D1000.append([antigen, killed])
            elif 1e-8 > affinity >= 1e-9:
                KD9D1000.append([antigen, killed])
            elif 1e-9 > affinity >= 1e-10:
                continue
            else:
                continue

    dataLists = [KD6D250, KD7D250, KD8D250, KD9D250,
                 KD6D500, KD7D500, KD8D500, KD9D500,
                 KD6D1000, KD7D1000, KD8D1000, KD9D1000]

    KD6D250AVG = []
    KD7D250AVG = []
    KD8D250AVG = []
    KD9D250AVG = []

    KD6D500AVG = []
    KD7D500AVG = []
    KD8D500AVG = []
    KD9D500AVG = []

    KD6D1000AVG = []
    KD7D1000AVG = []
    KD8D1000AVG = []
    KD9D1000AVG = []

    dataListsSorted = []

    avgLists = [KD6D250AVG, KD7D250AVG, KD8D250AVG, KD9D250AVG,
                KD6D500AVG, KD7D500AVG, KD8D500AVG, KD9D500AVG,
                KD6D1000AVG, KD7D1000AVG, KD8D1000AVG, KD9D1000AVG]

    for dList in dataLists:
        dList.sort(key=lambda x: x[0])
        dataListsSorted.append(dList)

    x = [0,10,20,30,40,50]
    for a in range(0,len(avgLists)):
        if dataListsSorted[a] != []:
            for i in range(0, 5):
                antSum = 0
                killSum = 0
                for j in range(x[i], x[i+1]):
                    antSum += dataListsSorted[a][j][0]
                    killSum += dataListsSorted[a][j][1]
                antAvg = antSum/len(dataListsSorted[a][x[i]:x[i+1]])
                killAvg = killSum / len(dataListsSorted[a][x[i]:x[i + 1]])
                avgLists[a].append([antAvg, killAvg])

    ax.plot([a[0]/get_max(KD6D250AVG, 0) for a in KD6D250AVG], [a[1]/get_max(KD6D250AVG, 1) for a in KD6D250AVG],
            linestyle=doseLineDict['250'], color=COLOR_DICT["CAR AFFINITY"]["1e-06"], linewidth=4)
    ax.plot([a[0]/get_max(KD7D250AVG, 0) for a in KD7D250AVG], [a[1]/get_max(KD7D250AVG, 1) for a in KD7D250AVG],
            linestyle=doseLineDict['250'],color=COLOR_DICT["CAR AFFINITY"]["1e-07"], linewidth=4)
    ax.plot([a[0]/get_max(KD8D250AVG, 0) for a in KD8D250AVG], [a[1]/get_max(KD8D250AVG, 1) for a in KD8D250AVG],
            linestyle=doseLineDict['250'], color=COLOR_DICT["CAR AFFINITY"]["1e-08"], linewidth=4)
    ax.plot([a[0]/get_max(KD9D250AVG, 0) for a in KD9D250AVG], [a[1]/get_max(KD9D250AVG, 1) for a in KD9D250AVG],
            linestyle=doseLineDict['250'], color=COLOR_DICT["CAR AFFINITY"]["1e-09"], linewidth=4)

    ax.plot([a[0]/get_max(KD6D500AVG, 0) for a in KD6D500AVG], [a[1]/get_max(KD6D500AVG, 1) for a in KD6D500AVG],
            linestyle=doseLineDict['500'], color=COLOR_DICT["CAR AFFINITY"]["1e-06"], linewidth=4)
    ax.plot([a[0]/get_max(KD7D500AVG, 0) for a in KD7D500AVG], [a[1]/get_max(KD7D500AVG, 1) for a in KD7D500AVG],
            linestyle=doseLineDict['500'], color=COLOR_DICT["CAR AFFINITY"]["1e-07"], linewidth=4)
    ax.plot([a[0]/get_max(KD8D500AVG, 0) for a in KD8D500AVG], [a[1]/get_max(KD8D500AVG, 1) for a in KD8D500AVG],
            linestyle=doseLineDict['500'], color=COLOR_DICT["CAR AFFINITY"]["1e-08"], linewidth=4)
    ax.plot([a[0]/get_max(KD9D500AVG, 0) for a in KD9D500AVG], [a[1]/get_max(KD9D500AVG, 1) for a in KD9D500AVG],
            linestyle=doseLineDict['500'], color=COLOR_DICT["CAR AFFINITY"]["1e-09"], linewidth=4)

    ax.plot([a[0]/get_max(KD6D1000AVG, 0) for a in KD6D1000AVG], [a[1]/get_max(KD6D1000AVG, 1) for a in KD6D1000AVG],
            linestyle=doseLineDict['1000'], color=COLOR_DICT["CAR AFFINITY"]["1e-06"], linewidth=4)
    ax.plot([a[0]/get_max(KD7D1000AVG, 0) for a in KD7D1000AVG], [a[1]/get_max(KD7D1000AVG, 1) for a in KD7D1000AVG],
            linestyle=doseLineDict['1000'], color=COLOR_DICT["CAR AFFINITY"]["1e-07"], linewidth=4)
    ax.plot([a[0]/get_max(KD8D1000AVG, 0) for a in KD8D1000AVG], [a[1]/get_max(KD8D1000AVG, 1) for a in KD8D1000AVG],
            linestyle=doseLineDict['1000'], color=COLOR_DICT["CAR AFFINITY"]["1e-08"], linewidth=4)
    ax.plot([a[0]/get_max(KD9D1000AVG, 0) for a in KD9D1000AVG], [a[1]/get_max(KD9D1000AVG, 1) for a in KD9D1000AVG],
            linestyle=doseLineDict['1000'], color=COLOR_DICT["CAR AFFINITY"]["1e-09"], linewidth=4)

    for color in affinityColorDict:
        ax.plot([0],[0], color=affinityColorDict[color], label=color, linestyle='solid', linewidth=4)
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose], linewidth=4)
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    ax.set_title("KILL CURVE (SIMULATED DATA)", fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel("ANTIGENS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel("% LYSIS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylim([-2.5, 1])
    ax.set_xlim([0, 1])
    ax.legend(bbox_to_anchor=(1, 1), frameon=False)

    plt.xticks([0.0, 0.5, 1.0], fontsize=TICKSIZE, rotation=45)
    # plt.yticks([0.0, 0.5, 1.0], fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_KILLCURVESIMNORM_' + str(TIME) + '.svg', bbox_inches='tight')

    return

def plot_kill_curve_normalized_exp_separated(SAVELOC):
    """Plot kill curve of experimental literature data where lysed % and antigen value are normalized to max value and all publications are on a separate plot."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    affinityColorDict = scripts.plot.plot_utilities.make_car_affinity_color_dict()
    litDict = scripts.plot.plot_utilities.make_literature_data_list()

    # Plot literature values
    for d in litDict:

        figKC = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
        ax = figKC.add_subplot(1, 1, 1)

        for CAR in d["CARS"]:
            if 0 >= d["CARS"][CAR]["KD"] >= 1e-5:
                d["CARS"][CAR].update({"Color": "grey"})
            elif 1e-5 > d["CARS"][CAR]["KD"] >= 1e-6:
                d["CARS"][CAR].update({"Color": COLOR_DICT["CAR AFFINITY"]["1e-06"]})
            elif 1e-6 > d["CARS"][CAR]["KD"] >= 1e-7:
                d["CARS"][CAR].update({"Color": COLOR_DICT["CAR AFFINITY"]["1e-07"]})
            elif 1e-7 > d["CARS"][CAR]["KD"] >= 1e-8:
                d["CARS"][CAR].update({"Color": COLOR_DICT["CAR AFFINITY"]["1e-08"]})
            elif 1e-8 > d["CARS"][CAR]["KD"] >= 1e-9:
                d["CARS"][CAR].update({"Color": COLOR_DICT["CAR AFFINITY"]["1e-09"]})
            elif 1e-9 > d["CARS"][CAR]["KD"] >= 1e-10:
                d["CARS"][CAR].update({"Color": COLOR_DICT["CAR AFFINITY"]["1e-10"]})
            else:
                d["CARS"][CAR].update({"Color": "black"})

            ax.errorbar(d["CARS"][CAR]["Normalized Data"]["Antigens"], d["CARS"][CAR]["Normalized Data"]["Kill %"],
                        xerr=d["CARS"][CAR]["Normalized Data"]["Ant Err"], fmt='None', yerr=d["CARS"][CAR]["Normalized Data"]["Kill % Err"],
                        ecolor='lightgray', zorder=0, linewidth=4)
            # ax.scatter(d["CARS"][CAR]["Data"]["Antigens"], d["CARS"][CAR]["Data"]["Kill %"],
            #            marker=d["MARKER"], color=d["CARS"][CAR]["Color"], zorder=1)
            ax.plot(d["CARS"][CAR]["Normalized Data"]["Antigens"], d["CARS"][CAR]["Normalized Data"]["Kill %"],
                        marker=d["MARKER"], color=d["CARS"][CAR]["Color"], zorder=1, linewidth=4, markersize=15)

        ax.scatter(None, None, marker=d["MARKER"], label=d["CITATION"], color='black', linewidth=4)

        for color in affinityColorDict:
            ax.plot([0],[0], color=affinityColorDict[color], label=color, linestyle='solid', linewidth=4)

        ax.set_title("NORMALIZED KILL CURVE\n(EXPERIMENTAL DATA)\n" + d["SAVE"], fontname='Arial',
                     fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
        ax.set_xlabel("NORMALIZED ANTIGENS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_ylabel("NORMALIZED % LYSIS", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
        ax.set_ylim([0, 1])
        ax.set_xlim([0, 1])
        ax.legend(bbox_to_anchor=(1, 0.40), frameon=False)

        plt.xticks([0.0, 0.5, 1.0], fontsize=TICKSIZE, rotation=45)
        plt.yticks([0.0, 0.5, 1.0], fontsize=TICKSIZE)

        if SAVELOC == '':
            plt.show()
        else:
            plt.savefig(SAVELOC + 'KILLCURVENORMEXP' + d["SAVE"] + '.svg', bbox_inches='tight')

    return