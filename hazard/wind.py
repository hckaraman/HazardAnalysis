# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 16:50:13 2020

@author: RVA04
"""
#############
import logging
logger = logging.getLogger("afad_backend.%s" % __name__)
logger.setLevel(logging.DEBUG) ## CHANGE THIS FOR LOG LEVEL
#logger.debug("Loading module...")
######################

import time
import rasterio
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
from bisect import bisect_left
import numpy as np
import scipy.optimize 
import math
import os

from config import config

from load_raster import _load_raster
def load_raster(path):
    logger.debug(f"loading {path}")
    image, lat, lon = _load_raster(path,True)
    return image, lat, lon

# TAKING THE NEAREST EXISTING POINT IN THE LOADED VECTOR
def take_closest(myList, myNumber): #take closest point
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
    #print("Before-After: {}, {}".format(before, after))
    if after - myNumber < myNumber - before:
       return after
    else:
       return before

#DEFINING THE ASSET POSITION
def Chose_Asset(p):
    image, lat, lon = load_raster(config.get('PATH','WindMap'))
    
    closest_x = take_closest(lon, float(p["Longitude"]))
    closest_y = take_closest(lat,  float(p["Latitude"]))
    ix=lon.index(closest_x)
    iy=lat.index(closest_y)
    wp=image[-iy][ix]/0.9
    return wp         

#4 POINT HAZARD CURVE
def FourPointWind(p):
    wp=Chose_Asset(p)
    
    Percent=[0.68,0.5,0.10,0.02]
    tr=[]
    maf=[]
    wind_peak=[]
    for i in range(len(Percent)):
        Tr=-50/(np.log(1-Percent[i]))
        MAF=-(np.log(1-Percent[i]))/50
        tr.append(Tr)
        maf.append(MAF)
        if -50/(np.log(1-Percent[i]))<50:
            cr=0.75*(1-0.2*np.log(-np.log(1-1/Tr)))**(0.5)
        else:
            cr=0.65*(1-0.138*np.log(-np.log(1-1/Tr)))
        wind_peak.append(cr*wp*1.5)
    return wind_peak,tr,maf

# PLOTTING 
def WindHazardCurve(p):
    wind_peak,tr,maf=FourPointWind(p)
    
    #interpolating the hazard curve
    x = np.array([round(i,6) for i in wind_peak])

    y = np.array([round(i,6) for i in maf])

    
#    print(x,y)
    def exponential(x,a,k,b):
        return a*np.exp(-x*k)
    popt, pcov = scipy.optimize.curve_fit(exponential, x, y, p0=[100,1, 1])
    
    a=(popt[0])
    k=(popt[1])
    
    return a,k


    
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def hazard_wind(Longitude,Latitude):
    p=pd.DataFrame({'Longitude':[Longitude],
                'Latitude':[Latitude]})
    HazardBase, HazardExponent = WindHazardCurve(p)  

    return  [HazardBase, HazardExponent]


#logger.debug("Module loaded.")