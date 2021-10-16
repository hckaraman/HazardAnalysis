# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 12:25:12 2020

@author: RVA04
"""

#############
import logging

logger = logging.getLogger("afad_backend.%s" % __name__)
logger.setLevel(logging.DEBUG)  ## CHANGE THIS FOR LOG LEVEL
# logger.debug("Loading module...")
######################


import time
import geopandas as gpd
import pandas as pd
import rasterio
import numpy as np
# from rasterio.features import shapes
# import descartes
from bisect import bisect_left
import math
import scipy.optimize
import matplotlib.pyplot as plt
import os
from shapely.geometry import Point

from config import config

# LOAD RASTER IMAGE
# def load_raster(path):
#     logger.debug(f"loading {path}")
#     with rasterio.Env():
#         with rasterio.open(path) as src:  
#             image = src.read(1)  ## first and only band  
#     #src.crs
#     #src.transform
#     #src.bounds
#     lat_image = [src.xy(i,0)[1] for i in range(src.height)]
#     lon_image = [src.xy(0,i)[0] for i in range(src.width)]

#     lat_image.reverse() # ordered list from max to min

#     return image, lat_image, lon_image
from load_raster import _load_raster


def load_raster(path):
    logger.debug(f"loading {path}")
    image, lat, lon = _load_raster(path, True)
    return image, lat, lon


def great_circle_vec(lat1, lng1, lat2, lng2, earth_radius=6371009):
    """
    Vectorized function to calculate the great-circle distance between two
    points or between vectors of points, using haversine.

    Parameters
    ----------
    lat1 : float or array of float
    lng1 : float or array of float
    lat2 : float or array of float
    lng2 : float or array of float
    earth_radius : numeric
        radius of earth in units in which distance will be returned (default is
        meters)

    Returns
    -------
    distance : float or vector of floats
        distance or vector of distances from (lat1, lng1) to (lat2, lng2) in
        units of earth_radius
    """

    phi1 = np.deg2rad(lat1)
    phi2 = np.deg2rad(lat2)
    d_phi = phi2 - phi1

    theta1 = np.deg2rad(lng1)
    theta2 = np.deg2rad(lng2)
    d_theta = theta2 - theta1

    h = np.sin(d_phi / 2) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(d_theta / 2) ** 2
    h = np.minimum(1.0, h)  # protect against floating point errors

    arc = 2 * np.arcsin(np.sqrt(h))

    # return distance in units of earth_radius
    distance = arc * earth_radius
    return distance


def take_closest(myList, myNumber):  # take closest point
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    # print("Before-After: {}, {}".format(before, after))
    if after - myNumber < myNumber - before:
        return after
    else:
        return before


def getRasterValueAtCoo(raster, x, y, x_ref, y_ref, crs='epsg:4326', out_index=True):
    """
    # raster image as numpy array
    # y_ref: latitude of the reference point
    # x_ref: longitude of the reference point
    # y: list of latitudes corresponding to raster grid y-axis (ASSUMED ORDERED LIST max2min--> iy will be then inverted in the output to follow image order)
    # x: list of longitudes corresponding to raster grid x-axis
    #crs = {'init': 'epsg:4326'}
    #closest_x = min(x, key=lambda list_value : abs(list_value - x_ref))
    #closest_y = min(y, key=lambda list_value : abs(list_value - y_ref))
    """
    closest_x = take_closest(x, x_ref)
    closest_y = take_closest(y, y_ref)
    #    print("Reference point: {}, {}".format(x_ref,y_ref))
    #    print("Closest point: {}, {}".format(closest_x, closest_y))
    ix = x.index(closest_x)
    iy = y.index(closest_y)
    if out_index: return raster[-iy][ix], ix, iy
    return raster[-iy][ix]


def getContextMorphology(raster, x, y, x_ref, y_ref, crs='epsg:4326', n=10, l_min=250):  # Altitude and Slope evaluation
    """
    Get sections at reference point
    n: number of points in each direction
    l_min: minimum section length
    """
    altitude, ix, iy = getRasterValueAtCoo(raster, x, y, x_ref, y_ref, crs=crs, out_index=True)
    dy = great_circle_vec(x[ix], y[iy], x[ix], y[iy + 1])  # distance in m along y
    dx = great_circle_vec(x[ix], y[iy], x[ix + 1], y[iy])  # distance in m along x
    #    print("dx = {} m, dy = {} m".format(dx,dy))

    AltNS = []
    for i in range(n * 2 + 1):
        if -iy + i > int(len(y)):
            AltNS.append(raster[-iy][ix])
        else:
            AltNS.append(raster[-iy + i][ix])

    SlopeNS = []
    for i in range(n * 2):
        if -iy + i > int(len(y)):
            SlopeNS.append(0)
        else:
            #            SlopeNS.append(((-raster[-iy+i][ix]+raster[-iy+i+1][ix])/dy)/math.pi*180)
            SlopeNS.append(np.degrees(
                np.arctan(
                    (-raster[-iy + i][ix] + raster[-iy + i + 1][ix]) / dy)
            )
            )

    SlopeVarNS = []
    for i in range(n * 2 - 1):
        if -iy + i > int(len(y)):
            SlopeVarNS.append(0)
        else:
            SlopeVarNS.append(SlopeNS[i] - SlopeNS[i + 1])

    #    print('SlopeVariationNS = '+str(SlopeVarNS)+' ;SlopeNS = '+str(SlopeNS)+' ;AltitudeNS = '+str(AltNS))

    Col1 = AltNS

    Col2_prov = []
    for i in range(len(SlopeNS) - 1):
        Col2_prov.append(max(SlopeNS[i + 1], SlopeNS[i]))
    Col2 = [Col2_prov[0]]
    for i in range(len(Col2_prov)):
        Col2.append(Col2_prov[i])
    Col2.append(Col2_prov[len(Col2_prov) - 1])

    Col3 = [SlopeVarNS[0]]
    for i in range(len(SlopeVarNS)):
        Col3.append(SlopeVarNS[i])
    Col3.append(SlopeVarNS[len(SlopeVarNS) - 1])

    Morph = [Col1, Col2, Col3]

    return Morph


def getSlopeOrientation(SNsec, EWsec):
    pass


# TEST - check altitude at a known location
def SlopeDataFrame(p):
    image, lat_image, lon_image = load_raster('./data/Flood/SRTM_TR_30_shaped.tif')
    Morph = getContextMorphology(image, lon_image, lat_image, float(p["Longitude"]), float(p["Latitude"]), n=10,
                                 l_min=250)
    return Morph


'''************************************************************'''
'''                 FRAGILITY CURVE OF SLOPE                   '''
'''************************************************************'''


def Fragility_Avalanche(p):
    Morph = SlopeDataFrame(p)

    pr1 = 0
    for i in range(len(Morph[0])):
        if Morph[0][i] > 2000 and Morph[1][i] > 30 and Morph[1][i] < 50 and Morph[2][
            i] > 10:  # 2700 montagne del caucaso
            pr1 = 1
            break
    return pr1


'''************************************************************'''
'''                    SLOPE HAZARD CURVE EQ                      '''
'''************************************************************'''


def HazardSlopeEQ(p):
    df_provinceID = pd.read_excel("./data/Avalanche/MAF_Provinces.xlsx")  # directory / 'MAF_Provinces.xlsx')
    provinceID = df_provinceID.values.tolist()
    df_Turkey_province = gpd.read_file("./data/Avalanche/gadm36_TUR_1.shp")  # directory / 'gadm36_TUR_1.shx')
    ProvinceShape = df_Turkey_province.values.tolist()
    asset = Point(float(p["Longitude"]), float(p["Latitude"]))
    pr2 = []
    maf = "maf"
    Province = "Prov"
    for i in range(len(df_Turkey_province)):
        if asset.within(ProvinceShape[i][0]) == True:
            ProvinceNum = provinceID[i][0]
            Province = provinceID[i][1]
            events = provinceID[i][2]
            maf = provinceID[i][3]

            if events == 0:
                #                print("no eventi")
                pr2.append(0)
            else:
                pr2.append(1)
    #            print("trovato")
    if maf == "maf":
        maf = 0
        Province = "Not Found"
        pr2 = [0]

    return maf / 100, Province, sum(pr2)


'''************************************************************'''
'''                  LANDSLIDE HAZARD CURVE EQ                    '''
'''************************************************************'''


def HazardCurve(p):
    Maf, Province, pr2 = HazardSlopeEQ(p)
    pr1 = Fragility_Avalanche(p)
    if pr1 * pr2 == 1:
        vel_max = 350  # km/h
        rho_max = 500  # kg/m3
        Pressure = (vel_max / 3.6) ** 2 * rho_max / 2 / 1000  # kPa (dynamic pressure formula)

        # Plot the hazard curve
        x_Pr = np.linspace(0.01, Pressure, 100)
        y_Maf = []
        for i in range(len(x_Pr)):
            y_Maf.append(Maf * 0.1345 * x_Pr[i] ** -0.627)
        '''
        plt.title("Avalanceh Hazard - Asset [Long:"+str(round(float(p["Longitude"]),2))+";Lat:"+str(round(float(p["Latitude"]),2))+"]")
        plt.loglog(x_Pr,y_Maf,'r')
        plt.xlabel('Pressure [kPa]')
        plt.ylabel('MAF [1/y]')
        plt.ylim(0,round(max(y_Maf)*1.1,1))
        plt.xlim(0,round(Pressure*1.1,-2))'''
        Maf = Maf * 0.1345
    else:
        #        print('NO HAZARD')

        Pressure = "NO HAZARD"
    return Pressure, Maf, 0.627


''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


def Hazard_avalanche_MAFPowerLaw(Longitude, Latitude):
    start_time = time.time()

    p = pd.DataFrame({'Longitude': [Longitude],
                      'Latitude': [Latitude]})
    Pressure, a, k = HazardCurve(p)
    return [a, k, Pressure]


''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# logger.debug("Module loaded.")

if __name__ == '__main__':
    lon, lat = 37, 37
    df = Hazard_avalanche_MAFPowerLaw(lon, lat)
    print("a")
