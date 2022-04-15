import scripts.plot.plot_utilities
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def make_env_dict():
    """Initialize empty enviornment informaton dictionary for plotting bar plots."""

    envDict = {'TIME': None,
               'AXIS': None,
               'CONC': None}

    return envDict

def make_env_df():
    """Initialize empty enviornment informaton dataframe for plotting bar plots."""

    columns=['TIME',
             'AXIS',
             'CONC']
    envDF = pd.DataFrame(columns=columns)

    return envDF

def plot_env_conc_times_bar(MOL_NAME, simsDF, COLOR, FILEID, SAVELOC, TIMES):
    """Plot bar plot of environment concentrations where data is grouped by time and colored by feature value."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    MOL_CONC_UNITS = scripts.plot.plot_utilities.define_mol_concentrations_units_dict()

    colorDict = COLOR_DICT[COLOR]
    HUES = {}

    SIM_TIMES = simsDF.iloc[0]['TIME']
    indicies = []
    TITLE_TIMES = ''
    for t in range(0, len(TIMES)):
        if 'VITRO' in FILEID:
            TITLE_TIMES += str(TIMES[t]) + ', '
        else:
            TITLE_TIMES += str(TIMES[t]-1) + ', '
        for s in range(0, len(SIM_TIMES)):
            if float(TIMES[t]) == float(SIM_TIMES[s]):
                indicies.append(s)
                if 'VITRO' in FILEID:
                    HUES[s] = TIMES[t]
                else:
                    HUES[s] = TIMES[t]-1

    TITLE_TIMES = TITLE_TIMES[0:-2]

    figEnv = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figEnv.add_subplot(1, 1, 1)
    YLABEL = MOL_NAME + ' CONC (' + MOL_CONC_UNITS[MOL_NAME] + ')'

    key = MOL_NAME + ' TOTAL CONC'
    envDF = make_env_df()
    order = []

    for i in range(0, len(simsDF)):
        for index in indicies:
            envDict = make_env_dict()
            y = simsDF.iloc[i][key][index]
            h = HUES[index]

            if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
                x = str(0)
            elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
                x = "CONTROL"
            else:
                x = str(simsDF.iloc[i][COLOR])

            envDict['TIME'] = h
            envDict['AXIS'] = x
            envDict['CONC'] = y

            envDF = envDF.append(envDict, ignore_index=True)

    if 'CONTROL' in envDF['AXIS'].to_list() and COLOR != 'ANTIGENS HEALTHY':
        order.append('CONTROL')

    for key in colorDict:
        if key in envDF['AXIS'].to_list():
            order.append(key)

    if 'VITRO' in FILEID:
        times = TIMES
    else:
        times = [t - 1 for t in TIMES]

    if MOL_NAME == 'GLUCOSE':
        ymax = 5e9
    elif MOL_NAME == 'IL-2':
        if 'VITRO' in FILEID:
            ymax = 4e7
        else:
            ymax = 1.4e7
    else:
        ymax = 160

    sns.barplot(ax=ax, x='TIME', y='CONC', hue='AXIS', palette=colorDict, data=envDF, order=times, hue_order=order, capsize=0.1, ci="sd")
    sns.swarmplot(ax=ax, x='TIME', y='CONC', hue='AXIS', color='black', dodge=True, order=times, hue_order=order, data=envDF)
    ax.set_xlabel(COLOR, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, labelpad=LABELPAD)
    ax.set_ylabel(YLABEL, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(MOL_NAME + "  CONCENTRATION\nOVER TIME " + TITLE_TIMES, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, pad=LABELPAD+7)
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
    ax.set_ylim(bottom=0, top=ymax)

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_CONC_' + MOL_NAME + '_XTIME_DAYS_' + TITLE_TIMES.replace(',','').replace(' ','')
                    + '.svg', bbox_inches='tight')
    return

def plot_evn_conc_times_line(MOL_NAME, simsDF, COLOR, FILEID, SAVELOC):
    """Plot given species concentrations over time."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    MOL_CONC_UNITS = scripts.plot.plot_utilities.define_mol_concentrations_units_dict()

    colorDict = COLOR_DICT[COLOR]

    figCounts = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figCounts.add_subplot(1, 1, 1)
    for i in range(0, len(simsDF)):

        if 'VITRO' in FILEID:
            plot_time = simsDF.iloc[i]['TIME']
            conc = simsDF.iloc[i][MOL_NAME + ' TOTAL CONC']
        else:
            plot_time = [t - 1 for t in simsDF.iloc[i]['TIME'][2:]]
            conc = simsDF.iloc[i][MOL_NAME + ' TOTAL CONC'][2:]

        if int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS CANCER':
            ax.plot(plot_time, conc,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(0)])
        elif int(simsDF.iloc[i]['DOSE']) == 0 and COLOR == 'ANTIGENS HEALTHY':
            ax.plot(plot_time, conc,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict["CONTROL"])
        else:
            ax.plot(plot_time, conc,
                    linestyle=doseLineDict[str(simsDF.iloc[i]['DOSE'])],
                    color=colorDict[str(simsDF.iloc[i][COLOR])])
    ax.set_xlabel("TIME (DAYS)", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(MOL_NAME + " CONCENTRATION\n(" + MOL_CONC_UNITS[MOL_NAME] + ")", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(MOL_NAME + " CONCENTRATION\nOVER TIME", fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlim([min(plot_time), max(plot_time)])

    if 'VITRO' in FILEID:
        ax.set_xticks(plot_time)
    else:
        ax.set_xticks(plot_time[::10])

    if MOL_NAME == 'GLUCOSE':
        ymax = 5e9
    elif MOL_NAME == 'IL-2':
        if 'VITRO' in FILEID:
            ymax = 4e7
        else:
            ymax = 1.4e7
    else:
        ymax = 160

    ax.set_ylim(bottom=0,top=ymax)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(plot_time[::2], [int(i) for i in plot_time[::2]], fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_CONC_' + MOL_NAME + '.svg', bbox_inches='tight')

    return

def plot_evn_concs_ideal_realistic_parity(MOL_NAME, simsDF, XAXIS, COLOR, TIMES, FILEID, SAVELOC):
    """Make parity plot of given species concentration in realistic vs ideal co-culture simulations for each simulation."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()
    doseLineDict = scripts.plot.plot_utilities.make_dose_line_dict()

    MOL_CONC_UNITS = scripts.plot.plot_utilities.define_mol_concentrations_units_dict()

    simsDFtreated = simsDF[simsDF['DOSE'] != 0]
    simsDFideal = simsDFtreated[simsDFtreated['ANTIGENS HEALTHY'] == 0]
    simsDFrealistic = simsDFtreated[simsDFtreated['ANTIGENS HEALTHY'] == 100]

    colorDict = COLOR_DICT[COLOR]

    figParity = plt.figure(figsize=(FIG_SIZE_X, FIG_SIZE_Y))
    ax = figParity.add_subplot(1, 1, 1)

    if MOL_NAME == 'GLUCOSE':
        ymax = 5e9
    elif MOL_NAME == 'IL-2':
        if 'VITRO' in FILEID:
            ymax = 7e7
        else:
            ymax = 1.4e7
    else:
        ymax = 160

    timeIndicies = []
    SIM_TIMES = simsDF.iloc[0]['TIME']
    for t in TIMES:
        for s in range(0,len(SIM_TIMES)):
            if float(t) == float(SIM_TIMES[s]):
                timeIndicies.append(s)

    timeMarkers = {'1': 's',
                   '4': 'd',
                   '7': 'o'}

    for time in timeIndicies:
        ideal = []
        realistic = []
        color = []
        marker = []
        for i in range(0,len(simsDFideal)):
            ideal.append(simsDFideal.iloc[i][MOL_NAME + ' TOTAL CONC'][time])
            simRealisticMatch = simsDFrealistic[simsDFrealistic['DOSE'] == simsDFideal.iloc[i]['DOSE']]
            simRealisticMatch = simRealisticMatch[simRealisticMatch['TREAT RATIO'] == simsDFideal.iloc[i]['TREAT RATIO']]
            simRealisticMatch = simRealisticMatch[simRealisticMatch['CAR AFFINITY'] == simsDFideal.iloc[i]['CAR AFFINITY']]
            simRealisticMatch = simRealisticMatch[simRealisticMatch['ANTIGENS CANCER'] == simsDFideal.iloc[i]['ANTIGENS CANCER']]
            simRealisticMatch = simRealisticMatch[simRealisticMatch['SEED'] == simsDFideal.iloc[i]['SEED']]
            realistic.append(simRealisticMatch.iloc[0][MOL_NAME + ' TOTAL CONC'][time])
            color.append(colorDict[str(simsDFideal.iloc[i][COLOR])])
            marker.append(timeMarkers[str(int(time/2))])

        if XAXIS == 'REALISTIC':
            YAXIS = 'IDEAL'
            ax.plot([0, ymax], [0, ymax], color='lightgray', zorder=0)
            ax.set_xlim([0, ymax])
            ax.set_ylim([0, ymax])
            for i in range(0, len(ideal)):
                ax.scatter(realistic[i], ideal[i], s=100, zorder=1, color=color[i], marker=marker[i])
        else:
            YAXIS = 'REALISTIC'
            ax.plot([0, ymax], [0, ymax], color='lightgray', zorder=0)
            ax.set_xlim([0, ymax])
            ax.set_ylim([0, ymax])
            for i in range(0, len(ideal)):
                ax.scatter(ideal[i], realistic[i], s=100, zorder=1, color=color[i], marker=marker[i])

    ax.set_title(MOL_NAME + " CONCENTRATION (" + MOL_CONC_UNITS[MOL_NAME] + ")\n" + "IN " + XAXIS + " vs " + YAXIS + " CO-CULTURE", fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD)
    ax.set_xlabel(MOL_NAME + " CONCENTRATION IN " + XAXIS, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(MOL_NAME + " CONCENTRATION IN " + YAXIS, fontname='Arial', fontweight='bold',
                  fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)

    for color in colorDict:
        ax.plot([0],[0], color=colorDict[color], label=color, linestyle='solid')
    for dose in doseLineDict:
        ax.plot([0],[0], color='black', label=dose, linestyle=doseLineDict[dose])
    ax.legend(bbox_to_anchor=(1.0,1.0), frameon=False)

    plt.xticks(fontsize=TICKSIZE)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_CONC_' + MOL_NAME + '_COCULTURE_PARITY_' + XAXIS + '_' + YAXIS + '_' + COLOR.replace(' ', '') + '.svg', bbox_inches='tight')

    return

def plot_env_conc_axis_bar(MOL_NAME, simsDF, COLOR, FILEID, SAVELOC, TIMES):
    """Plot bar plot of environment concentrations where data is grouped by feature value and colored by time."""

    COLOR_DICT = scripts.plot.plot_utilities.make_features_color_dict()
    FIG_SIZE_X, FIG_SIZE_Y, TICKSIZE, FONTSIZE_AXES_VALUES, FONTSIZE_AXES_TITLES, LABELPAD = scripts.plot.plot_utilities.define_plotting_globals()

    MOL_CONC_UNITS = scripts.plot.plot_utilities.define_mol_concentrations_units_dict()

    colorDict = COLOR_DICT[COLOR]
    HUES = {}

    SIM_TIMES = simsDF.iloc[0]['TIME']
    indicies = []
    TITLE_TIMES = ''
    for t in range(0, len(TIMES)):
        if 'VITRO' in FILEID:
            TITLE_TIMES += str(TIMES[t]) + ', '
        else:
            TITLE_TIMES += str(TIMES[t]-1) + ', '
        for s in range(0, len(SIM_TIMES)):
            if float(TIMES[t]) == float(SIM_TIMES[s]):
                indicies.append(s)
                if 'VITRO' in FILEID:
                    HUES[s] = TIMES[t]
                else:
                    HUES[s] = TIMES[t]-1

    TITLE_TIMES = TITLE_TIMES[0:-2]

    figEnv = plt.figure(figsize=(FIG_SIZE_X,FIG_SIZE_Y))
    ax = figEnv.add_subplot(1, 1, 1)
    YLABEL = MOL_NAME + ' CONC (' + MOL_CONC_UNITS[MOL_NAME] + ')'

    key = MOL_NAME + ' TOTAL CONC'
    envDF = make_env_df()
    order = []

    for i in range(0, len(simsDF)):
        for index in indicies:
            envDict = make_env_dict()
            y = simsDF.iloc[i][key][index]
            h = HUES[index]

            if COLOR == 'ANTIGENS CANCER' and int(simsDF.iloc[i]['DOSE']) == 0:
                x = str(0)
            elif COLOR == 'ANTIGENS HEALTHY' and int(simsDF.iloc[i]['DOSE']) == 0:
                x = "CONTROL"
            else:
                x = str(simsDF.iloc[i][COLOR])

            envDict['TIME'] = h
            envDict['AXIS'] = x
            envDict['CONC'] = y

            envDF = envDF.append(envDict, ignore_index=True)

    if 'CONTROL' in envDF['AXIS'].to_list() and COLOR != 'ANTIGENS HEALTHY':
        order.append('CONTROL')

    for key in colorDict:
        if key in envDF['AXIS'].to_list():
            order.append(key)

    if 'VITRO' in FILEID:
        times = TIMES
    else:
        times = [t - 1 for t in TIMES]

    if MOL_NAME == 'GLUCOSE':
        ymax = 5e9
    elif MOL_NAME == 'IL-2':
        if 'VITRO' in FILEID:
            ymax = 4e7
        else:
            ymax = 1.4e7
    else:
        ymax = 160

    sns.barplot(ax=ax, x='AXIS', y='CONC', hue='TIME', palette='Greys', data=envDF, order=order, hue_order=times, capsize=0.1, ci="sd")
    sns.swarmplot(ax=ax, x='AXIS', y='CONC', hue='TIME', color='black', dodge=True, order=order, hue_order=times, data=envDF)
    ax.set_xlabel(COLOR, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_ylabel(YLABEL, fontname='Arial', fontweight='bold', fontsize=FONTSIZE_AXES_VALUES, labelpad=LABELPAD)
    ax.set_title(MOL_NAME + "  CONCENTRATION\nOVER TIME " + TITLE_TIMES, fontname='Arial',
                 fontweight='bold', fontsize=FONTSIZE_AXES_TITLES, pad=LABELPAD+7)
    ax.legend(bbox_to_anchor=(1.0, 1.0), frameon=False)
    ax.set_ylim(bottom=0, top=ymax)

    plt.xticks(fontsize=TICKSIZE,  rotation = 45)
    plt.yticks(fontsize=TICKSIZE)

    if SAVELOC == '':
        plt.show()
    else:
        plt.savefig(SAVELOC + FILEID + '_CONC_' + MOL_NAME + '_XX_DAYS_' + TITLE_TIMES.replace(',','').replace(' ','')
                    + '.svg', bbox_inches='tight')
    return

