# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 09:50:03 2020

@author: RVA04
"""
#############
import logging
logger = logging.getLogger("afad_backend.%s" % __name__)
logger.setLevel(logging.DEBUG) ## CHANGE THIS FOR LOG LEVEL
#logger.debug("Loading module...")
######################

from data.constants import *
import geopandas as gpd
import rasterio
import numpy as np
from bisect import bisect_left
import pandas as pd
import scipy.optimize
import os
import matplotlib.pyplot as plt
from config import config

from load_raster import _load_raster
def load_raster(path):
    logger.debug(f"loading {path}")
    image, lat, lon = _load_raster(path,True)
    return image, lat, lon
    
# TAKING THE NEAREST EXISTING POINT IN THE LOADED VECTOR
def take_closest(myList, myNumber): #take closest point
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    #print("Before-After: {}, {}".format(before, after))
    if after - myNumber < myNumber - before:
       return after
    else:
       return before

#WIND INFLUENCE
def wind_influence(a_wind, k_wind):
    logger.debug("wind_influence called")
    #udner 5m/s 1st value, under 10m/s 2nd value, over 10 m/s the last
    m = [0.2936, 0.1937, 0.0346]
    q = [1.9767, 2.4704, 4.0965]
    
    MAF_1 = []
    FII = np.linspace(3,6,100)
    for i in range(len(FII)):
        if FII[i] < 3.4447:
            MAF_1.append(a_wind*np.exp(q[0]*k_wind/m[0]) * np.exp(-k_wind/m[0]*FII[i]))
        elif FII[i] < 4.4074:
            MAF_1.append(a_wind*np.exp(q[1]*k_wind/m[1]) * np.exp(-k_wind/m[1]*FII[i]))
        else:
            MAF_1.append(a_wind*np.exp(q[2]*k_wind/m[2]) * np.exp(-k_wind/m[2]*FII[i]))

    return MAF_1, FII

#LANDCOVER INFLUENCE
def landcover(Longitude,Latitude):
    logger.debug("landcover called")
    image, lat_image, lon_image = load_raster(config.get("PATH","LandCoverMap"))
    valx = take_closest(lon_image, Longitude)
    valy = take_closest(lat_image, Latitude)
    idx = lon_image.index(valx)
    idy = lat_image.index(valy)

    step = 15  
    land = []
    for i in range(idy-step,idy+step):
        for j in range(idx-step,idx+step):
            iim = image[-i][j]
            if iim == 244 or iim > 310 and iim < 314:
                land.append(3)
            else:
                land.append(0)
    tmp = max(land)
    return tmp

def landcover_influence(Longitude, Latitude):
    cover = landcover(Longitude,Latitude)
    logger.debug("landcover_influence called")
    if cover == 3:
        moisture = pd.read_excel(config.get("PATH","Moisture_Turkey"))
        valx = take_closest(list(moisture["longitude"]), Longitude)
        valy = take_closest(list(moisture["latitude"]), Latitude)
        tmp = moisture[moisture["latitude"] == valy]
        tmp = tmp[moisture["longitude"] == valx]
        moist = float(tmp['moisture'])  # in percentage
        del tmp

        relation = pd.read_excel(config.get("PATH","Moisture_CLC_Fire"))
        if cover == 1:
            tmp = relation['Cropland']
        elif cover == 2:
            tmp = relation['Pasture']
        elif cover == 3:
            tmp = relation['Forest']
        elif cover == 4:
            tmp = relation['Grass']

        if moist < 0.475:
            diff = []
            for i in range(len(tmp)):
                dif = moist - relation['Moisture'].loc[i]
                if dif > 0:
                    diff.append(dif)
                else:
                    diff.append(1)
            row = diff.index(min(diff))
            MAF_2 = tmp[row] + (moist - relation['Moisture'][row]) / (
                        relation['Moisture'][row + 1] - relation['Moisture'][row]) * (tmp[row + 1] - tmp[row])
        else:
            MAF_2 = 0.01

    elif cover == 0:
        MAF_2 = "NO HAZARD"
    return MAF_2


#HISTORICAL PROBABILITY
def historical_prob():
    A_forest = 210000
    year = 50
    A_burn = 6084
    MAF_3 = (A_burn / A_forest) / year

    return MAF_3

#EXPONENTIAL PARAMETERS
def WildfireHazardCurve(a_wind, k_wind,Longitude, Latitude):
    logger.debug("WildfireHazardCurve called")
    MAF_1, FII = wind_influence(a_wind, k_wind)
    MAF_2 = landcover_influence(Longitude, Latitude)
    MAF_3 = historical_prob()

    if MAF_2 == "NO HAZARD":
        a = "NO HAZARD"
        k = "NO HAZARD"
    else:
        MAF = []
        for i in range(len(MAF_1)):
            MAF.append(MAF_1[i] * MAF_2 * MAF_3)

        x = FII
        y = MAF

        def exponential(x,a,k,b):
            return a*np.exp(-x*k)
        popt, pcov = scipy.optimize.curve_fit(exponential, x, y, p0=[100000,10, 1])

        a=(popt[0])
        k=(popt[1])
    
    return a,k

def hazard_wildfire(a_wind, k_wind, Longitude, Latitude):
    HazardBase, HazardExponent = WildfireHazardCurve(a_wind, k_wind, Longitude, Latitude)
    return  HazardBase, HazardExponent

#logger.debug("Module loaded.")



