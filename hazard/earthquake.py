#############
import logging
logger = logging.getLogger("afad_backend.%s" % __name__)
logger.setLevel(logging.DEBUG) ## CHANGE THIS FOR LOG LEVEL
#logger.debug("Loading module...")
######################


import numpy as np
import pandas as pd
import scipy.optimize as opt
import matplotlib.pyplot as plt
from data.constants import *
from config import config

inf = float('inf')


seismic_data = pd.read_csv(config.get("PATH","earthquake_DB")) # DATA / "Earthquake/earthquake_DB.csv")
# Si vede che con l'aumentare dell'indice aumenta la LONGITUDINE !!! 
# ma va scelta anche la latitudine perche è da sud a nord!!
Seismic_list = seismic_data.values.tolist()

Long = []
Lat = []
for r in range(len(Seismic_list)):
    Long.append(Seismic_list[r][0])
    Lat.append(Seismic_list[r][1])


def parametriccurve(datiPGA, alpha, MAF_asy, PGA_asy):
    MAF_PGA = MAF_asy * np.exp(alpha * ((np.log(datiPGA / PGA_asy)) ** (-1)))
    return MAF_PGA


def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


def getIndexPositions(listOfElements, element):
    ''' Returns the indexes of all occurrences of give element in
    the list- listOfElements '''
    indexPosList = []
    indexPos = 0
    while True:
        try:
            # Search for item in list from indexPos to the end of list
            indexPos = listOfElements.index(element, indexPos)
            # Add the index position in list
            indexPosList.append(indexPos)
            indexPos += 1
        except ValueError as e:
            break
    return indexPosList


def find_pgas(coordinates):

    Closest_Long = closest(Long, coordinates[0])
    elements = getIndexPositions(Long, Closest_Long)
    Lat_Closest = []
    for j in elements:
        Lat_Closest.append(Lat[j])

    Closest_Lat = closest(Lat_Closest, coordinates[1])
    ind = Lat_Closest.index(Closest_Lat)
    Index = elements[0] + ind
    # the value of the relative Lat value
    PGAs = []
    for f in range(len(TR)):
        PGAs.append(round(seismic_data.iloc[(Index)][TR[f]], 5))
    return PGAs


def exponential(x, a, k, b):
    return a * np.exp(-x * k)


def Seismic_hazard(Main_Folder, po, pga, tn, pga_list, category, assetID,
                   draw=False):


    popt, pcov = opt.curve_fit(exponential, pga, tn, p0=[1, 1, 1])
    a = (popt[0])
    k = (popt[1])

    y_model = []
    for i in range(len(pga_list)):
        y_model.append(a * np.exp(-pga_list[i] * k))
    y_model[0] = 0.1  # To be decided if applied or not???

    if draw:
        plt.rcParams["figure.figsize"] = [8, 5]
        plt.title('Earthquake HAZARD CURVE: ' + category + ' - ' + str(assetID))
        plt.xlabel('MEAN IM')
        plt.ylabel('MAF [1/TR]')
        plt.ylim(0, 0.11)

        plt.plot(pga_list, y_model, 'b',
                 label='Hazard curve fitted from Turksih Seismic DB',
                 linestyle='-.')
        plt.legend()
        plt.savefig(
            "Results/Hazard/" + "Earthquake_Hazard_Curve_for_" + str(
                assetID) + ".png")
        plt.show()
    #    if SHOW_PLOT=="NO":
    #       plt.close()
    return y_model, [a, k]


def Plot_hazard_curve(a, k, Range, max_value, hazard, assetID, curve_type, draw=False):
    logger.debug("Plot_hazard_curve called")#(%s,%s,%s,%s,%s,%s,%s)" % (a, k, Range, max_value, hazard, assetID, curve_type))
    y_value = []

    if curve_type == "Exponential":
        for i in range(len(Range)):
            if a * np.exp(-Range[i] * k) < 10 ** -100:
                y_value.append(0)
            else:
                y_value.append(a * np.exp(-Range[i] * k))
    #        y_value[0]=0.1
    else:
        for i in range(len(Range)):
            y_value.append(a * Range[i] ** -k)

    if draw:
        plt.rcParams["figure.figsize"] = [8, 5]
        plt.title(hazard + ' HAZARD CURVE: ' + assetID)
        plt.xlabel('MEAN IM')
        plt.ylabel('MAF [1/TR]')
        #        plt.ylim(0,0.11)
        plt.xlim(0, max_value)

        plt.plot(Range, y_value, 'b', label='Hazard curve', linestyle='-.')
        plt.legend()
        plt.savefig(
            "Results/Hazard/" + hazard + " Hazard Curve for " + str(
                assetID) + ".png")
        plt.show()
    return y_value

#logger.debug("Module loaded.")

# if __name__ == '__main__':
#     lon,lat = 37,37
#     df = Hazard_avalanche_MAFPowerLaw(lon,lat)
#     print("a")