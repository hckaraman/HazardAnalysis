# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:42:22 2020

@author: RVA04
"""
#############
import logging
logger = logging.getLogger("afad_backend.%s" % __name__)
logger.setLevel(logging.DEBUG) ## CHANGE THIS FOR LOG LEVEL
#logger.debug("Loading module...")
######################

# SRTM (Risoluzione 90 m x 80 m)
# From RASTER (e.g. TIFF) to grid
import rasterio
import numpy as np
from rasterio.features import shapes
import pandas as pd
import geopandas as gpd
from bisect import bisect_left
import math
import scipy.optimize
import matplotlib.pyplot as plt
import os

from config import config


#%%
'''************************************************************'''
'''                   MORPHOLOGY DEFINITION                    '''
'''************************************************************'''

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
    image, lat, lon = _load_raster(path,True)
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


def euclidean_dist_vec(y1, x1, y2, x2):
    """
    Vectorized function to calculate the euclidean distance between two points
    or between vectors of points.

    Parameters
    ----------
    y1 : float or array of float
    x1 : float or array of float
    y2 : float or array of float
    x2 : float or array of float

    Returns
    -------
    distance : float or array of float
        distance or vector of distances from (x1, y1) to (x2, y2) in graph units
    """

    # euclid's formula
    distance = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
    return distance

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
    '''print("Reference point: {}, {}".format(x_ref,y_ref))
    print("Closest point: {}, {}".format(closest_x, closest_y))'''
    ix=x.index(closest_x)
    iy=y.index(closest_y)
    if out_index: return raster[-iy][ix], ix, iy
    return raster[-iy][ix]

def getSlopeMorphology(a, b, th1=0.05, th2=0.02): #morphology of the slope

    if a>th1:
        if b>th1: return 1, "++.++"
        elif b>th2: return 1, "++.+"
        elif b<-th1: return 3, "++.--"
        elif b<-th2: return 3, "++.-"
        else: return 2, "++.0"
    elif a>th2:
        if b>th1: return 1, "+.++"
        elif b>th2: return 1, "+.+"
        elif b<-th1: return 3, "+.--"
        elif b<-th2: return 3, "+.-"
        else: return 2, "+.0"
    elif a<-th1:
        if b>th1: return 5, "--.++"
        elif b>th2: return 5, "--.+"
        elif b<-th1: return 7, "--.--"
        elif b<-th2: return 7, "--.-"
        else: return 8, "--.0"
    elif a<-th2:
        if b>th1: return 5, "-.++"
        elif b>th2: return 5, "-.+"
        elif b<-th1: return 7, "-.--"
        elif b<-th2: return 7, "-.-"
        else: return 8, "-.0"
    else:
        if b>th1: return 9, "0.++"
        elif b>th2: return 9, "0.+"
        elif b<-th1: return 6, "0.--"
        elif b<-th2: return 6, "0.-"
        else: return 4, "0.0"

def getContextMorphology(raster, x, y, x_ref, y_ref, crs='epsg:4326', n=3, l_min = 500): # Altitude and Slope evaluation
    """
    Get sections at reference point
    n: number of points in each direction
    l_min: minimum section length
    """
    altitude, ix, iy = getRasterValueAtCoo(raster, x, y, x_ref, y_ref, crs=crs, out_index=True)

    dy=great_circle_vec(x[ix],y[iy],x[ix],y[iy+1]) # distance in m along y
    dx=great_circle_vec(x[ix],y[iy],x[ix+1],y[iy]) # distance in m along x
    '''print("dx = {} m, dy = {} m".format(dx,dy))'''

    n1 = math.ceil(max(l_min/(dx*2),l_min/(dy*2)))
    if n1>n: n=n1

    SN_section = [raster[i][ix] for i in range(-iy-n,-iy+n+1)]
    WE_section = [raster[-iy][i] for i in range(ix-n,ix+n+1)]
    SN_slope = [(SN_section[i]-SN_section[i+1])/dy for i in range(len(SN_section)-1)]
    WE_slope = [(WE_section[i]-WE_section[i+1])/dx for i in range(len(WE_section)-1)]

    Altitude_SN=[]
    for i in range(ix-n,ix+n+1):
        Altitude_SN.append([])
        ii=len(Altitude_SN)-1
        for j in range(-iy-n,-iy+n+1):
            Altitude_SN[ii].append(raster[j][i])

    SE_slope=[]
    NE_slope=[]
    for i in range(len(Altitude_SN)-1):
        SE_slope.append(-(Altitude_SN[i+1][i+1]-Altitude_SN[i][i])/(dx**2+dy**2)**(1/2))
        NE_slope.append(-(Altitude_SN[i+1][len(Altitude_SN)-2-i]-Altitude_SN[i][len(Altitude_SN)-1-i])/(dx**2+dy**2)**(1/2))

    for i in range(len(NE_slope)):
        NE_slope[i]=np.degrees(np.arctan(NE_slope[i]))
        SE_slope[i]=np.degrees(np.arctan(SE_slope[i]))
        SN_slope[i]=np.degrees(np.arctan(SN_slope[i]))
        WE_slope[i]=np.degrees(np.arctan(WE_slope[i]))

    for i in range(int(len(NE_slope)/2)):
        ii=i+int(len(NE_slope)/2)
        NE_slope[ii]=-NE_slope[ii]
        SE_slope[ii]=-SE_slope[ii]
        SN_slope[ii]=-SN_slope[ii]
        WE_slope[ii]=-WE_slope[ii]

    d_slope={'NE':NE_slope,
       'SE':SE_slope,
       'SN':SN_slope,
       'WE':WE_slope}

    Slope=pd.DataFrame(data=d_slope)

    d_area={'NE':[2*dx*dy,1.5*dx*dy,dx*dy,dx*dy,1.5*dx*dy,2*dx*dy],
            'SE':[2*dx*dy,1.5*dx*dy,dx*dy,dx*dy,1.5*dx*dy,2*dx*dy],
            'SN':[1.5*dx*dy,dx*dy,dx*dy,dx*dy,dx*dy,1.5*dx*dy],
            'WE':[1.5*dx*dy,dx*dy,dx*dy,dx*dy,dx*dy,1.5*dx*dy]}

    Area=pd.DataFrame(data=d_area)

    """
    case 1   case 2    case 3    case 4   case 5   case 6   case 7   case 8   case 9
    s\.      s\._      s\./n     s_._n      .       _./n      ./n      ._n    s_.
      \n         n                        s/ \n    s        s/       s/          \n

    ++.++    ++.0+     ++.--     0+.0-    --.++    0-.--    --.--    --.0-    0+.++
    ++.+     ++.0-     ++.-      0-.0+    --.+     0+.--    --.-     --.0+    0-.++
     +.++               +.--     0-.0-     -.++              -.--     -.0      0.+
     +.+
    thresholds:
        0+ = 0% to 2%
        +  = 2% to 5%
        ++ = more than 5%
    """
    SNn, SNmorph = getSlopeMorphology(SN_slope[n-1], SN_slope[n])
    WEn, WEmorph = getSlopeMorphology(WE_slope[n-1], WE_slope[n])
    '''print("SNn = {}, EWn = {}".format(SNn,WEn))'''
    SlopeOrientation = ""
    if SNn in [1,2,9]: SlopeOrientation+="N"
    elif SNn in [6,7,8]: SlopeOrientation+="S"
    if WEn in [1,2,9]: SlopeOrientation+="E"
    elif WEn in [6,7,8]: SlopeOrientation+="W"

    return altitude,Altitude_SN, Slope,SlopeOrientation,dx,dy
    #return altitude, ix, iy, SN_section, WE_section, SN_slope, WE_slope, SNmorph, WEmorph

def getSlopeOrientation(SNsec, EWsec):
    pass

'''************************************************************'''
'''************************************************************'''
'''************************************************************'''

#%% HWSD  (Risoluzione 900 m x 800 m)

'''************************************************************'''
'''              SOIL CHATATERISTICS SELECTION                 '''
'''************************************************************'''

def getSoil(df_HWSD, val, full_HWSD_row=False):

    """
    ASSUMPTIONS:We consider just the main soil unit (i.e. SEQ=1)

    if REF_DEPTH == 10 cm  =>  Lithosols of FAO-74 and Lithic Leptosols of FAO-90       =>  ROCK
    if REF_DEPTH == 30 cm  =>  Rendzinas and Rankers of FAO-74 and Leptosols of FAO-90  =>  ROCK???
    if REF_DEPTH == 100 cm =>  get S_USDA_TEX_CLASS

    Return a Code Texture:
        0 ROCK
        1 clay (heavy)
        2 silty clay
        3 clay
        4 silty clay loam
        5 clay loam
        6 silt
        7 silt loam
        8 sandy clay
        9 loam
        10 sandy clay loam
        11 sandy loam
        12 loamy sand
        13 sand

    """
    row=df_HWSD[df_HWSD["MU_GLOBAL"]==val]
    #seq=row["SEQ"].tolist()
    #i=seq.pos(1)

    depth_row=row[row["SEQ"]==1]
    depth=depth_row["REF_DEPTH"]
    depth=depth.tolist()[0]
    if depth == 10: soil=0
    elif depth == 30: soil=0
    else: soil=depth_row["S_USDA_TEX_CLASS"].tolist()[0]

    if full_HWSD_row: return soil, row
    return soil

'''************************************************************'''
'''************************************************************'''
'''************************************************************'''


'''************************************************************'''
'''                 FRAGILITY CURVE OF SLOPE EQ                  '''
'''************************************************************'''

def Fragility_PvsEQ(Slope,Soil,soil):
    logger.debug('Fragility_PvsEQ called')
    IMc=[]
    FS_eq=[]
    Hl=[]
    for i in range(4):
        IMc.append([])
        FS_eq.append([])
        Hl.append([])
        for j in range(len(Slope)):
            if i==0 :
                k="NE"
            else:
                if i==1:
                    k="SE"
                else:
                    if i==2:
                        k="SN"
                    else:
                        k="WE"

            if soil!=0:
                #Soil paramters
                gam=Soil[3]     #kN/m3
                phi=Soil[4]     #degree
                c=Soil[5]       #kPa
                gam_w=9.81      #kN/m3
                tb15=Soil[10]
                tb30=Soil[11]
                tb45=Soil[12]
                tb60=Soil[13]
                lb15=Soil[14]
                lb30=Soil[15]
                lb45=Soil[16]
                lb60=Soil[17]
                #Slope paramters
                if Slope[k][j]<15:
                    dip=0.0001
                    phi=90          #discarded from the computation
                else:
                    dip=Slope[k][j]        #degree

                L=1             #m
    #                    Linc=L/math.cos(math.radians(dip))
                if dip<15:
                    H=tb15
                    L=lb15/2
                else:
                    if dip<30:
                        H=tb15-(dip-15)/15*(tb15-tb30)
                        L=(lb15-(dip-15)/15*(lb15-lb30))/2
    #                    H=tb30
    #                    L=lb30/2
                    else:
                        if dip<45:
                            H=tb30-(dip-30)/15*(tb30-tb45)
                            L=(lb30-(dip-30)/15*(lb30-lb45))/2
                        else:
                            if dip<60:
                                H=tb45-(dip-45)/15*(tb45-tb60)
                                L=(lb45-(dip-45)/15*(lb45-lb60))/2
                            else:
    #                            if (dip==60 and dip>60):
                                H=tb60
                                L=lb60/2

                Ip=(H**2+L**2)**(1/2)

                depthW=H       #m

                ''' ********************************* '''
                '''FRAGILITY CURVE OF SLOPE'''

                if depthW>H:
                    hw=0
                else:
                    hw=(H-depthW)/H

                alfa_dw=math.degrees(math.atan(H/L))
                psi_dw=dip-alfa_dw

                psi_up=dip+alfa_dw
                if psi_up>90:
                    psi_up=90

                #FS evelautino
                C=c*2*Ip
                Wdw=H*L/2*gam
                Wup=H*L/2*gam
                Ndw=Wdw*math.cos(math.radians(psi_dw))
                Nup=Wup*math.cos(math.radians(psi_up))
                FST=(Ndw*math.tan(math.radians(phi))+Nup*math.tan(math.radians(phi))+C)/((Wdw*math.sin(math.radians(psi_dw)))+(Wup*math.sin(math.radians(psi_up))))
                FS_eq[i].append(FST)

                #Intensity measure evaluation
                IMc[i].append((FST-1)*math.tan(math.radians(dip))/(math.tan(math.radians(phi))*math.tan(math.radians(dip))+1)) #the formula has *g, but the IM is in g and not m/s2
                Hl[i].append(H/3)
            else:
                IMc[i].append(0)
                Hl[i].append(0)
                FS_eq[i].append(0)

    df_H=pd.DataFrame(data={'NE':Hl[0],'SE':Hl[1],'SN':Hl[2],'WE':Hl[3]})
    df_FS=pd.DataFrame(data={'NE':FS_eq[0],'SE':FS_eq[1],'SN':FS_eq[2],'WE':FS_eq[3]})

    return   IMc,Hl

#%%

'''************************************************************'''
'''                    SLOPE HAZARD CURVE EQ                      '''
'''************************************************************'''

def HazardSlopeEQ(seismic_data,Seismic_list,p):
    '''List of PGAs forthe TR'''

    '''SELECTION OF THE ASSET HAZARD CURVE'''
    Percent=[0.02,0.10,0.50,0.68]
    MAF=[]
    for i in range(len(Percent)):
        MAF.append(-np.log(1-Percent[i])/50)
    Long=[]
    Lat=[]
    for r in range(len(Seismic_list)):
        Long.append(Seismic_list[r][0])
        Lat.append(Seismic_list[r][1])

    SeismicLon=list(seismic_data['Longitude'])
    SeismicLat=list(seismic_data['Latitude'])
    Asset_Lon=take_closest(SeismicLon,float(p["Longitude"]))
    Asset_Lat=take_closest(SeismicLat,float(p["Latitude"]))
    AssetRow=seismic_data.loc[seismic_data['Longitude'] == Asset_Lon]
    k=AssetRow.loc[AssetRow['Latitude'] == Asset_Lat]
    AssetPGA=k
    IM=[list(k['2%'])[0],list(k['10%'])[0],list(k['50%'])[0],list(k['68%'])[0]]
    ''' ********************************* '''

    #interpolating the hazard curve
    x = np.array(IM)
    y = np.array(MAF)
    def exponential(x,a,k,b):
        return a*np.exp(-x*k)
    popt, pcov = scipy.optimize.curve_fit(exponential, x, y, p0=[0.5,5, 1])
    a=(popt[0])
    k=(popt[1])

    '''
    #Plot the 4point curve
    plt.figure()
    x=np.linspace(x[0],x[len(x)-1],20)
    plt.title("Earthquake Hazard - Asset [Long:"+str(round(float(p["Longitude"]),2))+";Lat:"+str(round(float(p["Latitude"]),2))+"]")
    plt.plot(IM,MAF,'bo',x,a*np.exp(-x*k),'r')
    plt.xlabel('Acceleration [g]')
    plt.ylabel('MAF [1/y]')
   '''

    return AssetPGA,a,k

'''************************************************************'''
'''                  LANDSLIDE HAZARD CURVE EQ                    '''
'''************************************************************'''

def LanslideHazardCurveEQ(a,k,Slope,IMc,Hl,p,d_soilDB,dx,dy,altitude,Altitude_SN):

    #Escluding values not realistic and Putting in one single list
    haz=[]
    for i in range(4):
        for j in range(len(Slope)):
            IM=IMc[i][j]
            h=a*np.exp(-IM*k)
            if h<1*10**(-10):
                haz.append(0)
            else:
                haz.append(h)

    #Putting in one single list
    IM=[]
    for i in range(4):
        for j in range(len(Slope)):
            IM.append(IMc[i][j])
    Height=[]
    for i in range(4):
        for j in range(len(Slope)):
            Height.append(Hl[i][j])

    #Ordering the values from highest to lowest hazard
    hazard_long=pd.DataFrame(data={'Hazard':haz,'IntensityMeasure':IM,'Height':Height})
    hazard=hazard_long.sort_values('Hazard',ascending=False)

    #Defining new name for some of the output
    Resolution= pd.DataFrame({'dx':[dx],'dy':[dy]})
    Asset_Altitude=altitude
    Asset_Location=p
    SoilProperties_DB=d_soilDB
    AltitudeSN=pd.DataFrame(data={'0':Altitude_SN[0],'1':Altitude_SN[1],'2':Altitude_SN[2],'3':Altitude_SN[3],'4':Altitude_SN[4],'5':Altitude_SN[5],'6':Altitude_SN[6]})

    #creating list from DB
    h=list(hazard['Hazard'])
    d=list(hazard['Height'])
    im=list(hazard['IntensityMeasure'])

    #discarding all nil values
    hazard=[]
    depth=[]
    intensity=[]
    for i in range(len(h)):
        if h[i]==0:
            i
        else:
            hazard.append(h[i])
            depth.append(d[i])
            intensity.append(im[i])

    h_eq=hazard
    d_eq=depth

    if len(h_eq)>2:
        x=depth
        y=hazard
        def exponential(x,a,k,b):
            return a*np.exp(-x*k)
        popt, pcov = scipy.optimize.curve_fit(exponential, x, y, p0=[100,1, 1])
        IndBaseEQHazard=(popt[0])
        IndExponentEQHazard=(popt[1])
    else:
        IndBaseEQHazard=(None)
        IndExponentEQHazard=(None)
    return AltitudeSN,Resolution,Asset_Altitude,Asset_Location,SoilProperties_DB,hazard,h_eq,d_eq,IndBaseEQHazard,IndExponentEQHazard

#%%

'''************************************************************'''
'''                 FRAGILITY CURVE OF SLOPE rain               '''
'''************************************************************'''

def Fragility_PvsRain(Slope,Soil,soil):

    IMc=[]
    Hl=[]
    FS_rain=[]
    for i in range(4):
        IMc.append([])
        Hl.append([])
        FS_rain.append([])
        for j in range(len(Slope)):
            IMc[i].append([])
            Hl[i].append([])
            if i==0 :
                k="NE"
            else:
                if i==1:
                    k="SE"
                else:
                    if i==2:
                        k="SN"
                    else:
                        k="WE"
            if soil!=0:
                #Soil paramters
                Sr0=0.3
                n=Soil[2]/100
                gam=Soil[3]     #kN/m3
                phi=Soil[4]     #degree
                c=Soil[5]       #kPa
                kt=Soil[6]       #m/s
                f0=Soil[7]      #mm/h
                fc=Soil[8]      #mm/h
                kf=Soil[9]      #1/h
                gam_w=9.81      #kN/m3
                tb15=Soil[10]
                tb30=Soil[11]
                tb45=Soil[12]
                tb60=Soil[13]
                lb15=Soil[14]
                lb30=Soil[15]
                lb45=Soil[16]
                lb60=Soil[17]
                #Slope paramters
                if Slope[k][j]<15:
                    dip=0.0001
                    phi=90          #discarded from the computation
                else:
                    dip=Slope[k][j]        #degree

                L=1             #m
    #                    Linc=L/math.cos(math.radians(dip))
                if dip<15:
                    H=tb15
                    L=lb15/2
                else:
                    if dip<30:
                        H=tb15-(dip-15)/15*(tb15-tb30)
                        L=(lb15-(dip-15)/15*(lb15-lb30))/2
    #                    H=tb30
    #                    L=lb30/2
                    else:
                        if dip<45:
                            H=tb30-(dip-30)/15*(tb30-tb45)
                            L=(lb30-(dip-30)/15*(lb30-lb45))/2
                        else:
                            if dip<60:
                                H=tb45-(dip-45)/15*(tb45-tb60)
                                L=(lb45-(dip-45)/15*(lb45-lb60))/2
                            else:
    #                            if (dip==60 and dip>60):
                                H=tb60
                                L=lb60/2

                Ip=(H**2+L**2)**(1/2)

                depthW=H       #m

                ''' ********************************* '''
                '''FRAGILITY CURVE OF SLOPE'''

                if depthW>H:
                    hw=0
                else:
                    hw=(H-depthW)/H

                alfa_dw=math.degrees(math.atan(H/L))
                psi_dw=dip-alfa_dw

                psi_up=dip+alfa_dw
                if psi_up>90:
                    psi_up=90

                #FS evelautino
                C=c*2*Ip
                Wdw=H*L/2*gam
                Wup=H*L/2*gam
                Ndw=Wdw*math.cos(math.radians(psi_dw))
                Nup=Wup*math.cos(math.radians(psi_up))
                FST=(Ndw*math.tan(math.radians(phi))+Nup*math.tan(math.radians(phi))+C)/((Wdw*math.sin(math.radians(psi_dw)))+(Wup*math.sin(math.radians(psi_up))))
                FS_rain[i].append(FST)

                minutes=[60,120,240,720]
                for z in range(len(minutes)):
                    def Funct(x):
                        F=gam_w*(math.sin(math.radians(dip)))*(math.cos(math.radians(dip)))*x*H*L
                        return (Ndw*math.tan(math.radians(phi))+Nup*math.tan(math.radians(phi))+C)/((Wdw*math.sin(math.radians(psi_dw)))+(Wup*math.sin(math.radians(psi_up)))+F) - [1]
                    x = scipy.optimize.broyden1(Funct, [1], f_tol=1e-14)
                    if x<0:
                        h_rain=0
                    else:
                        if x>1:
                            h_rain=0
                        else:
                            h_rain=(x-hw)*(n*H*(1-Sr0))/(math.exp(-kt*(((math.sin(math.radians(dip)))/(n*L*(1-Sr0)))*(minutes[z]*60))))
                            h_rain=h_rain/(minutes[z]/60)*1000
                            h_rain=h_rain[0]
    ##                        Infiltration=fc+(f0-fc)*math.exp(-kf*minutes/60)
    #                        Infiltration=fc*minutes[z]/60+(f0-fc)/kf*(1-math.exp(-kf*minutes[z]/60))
    #                        if h_rain*minutes[z]/60<Infiltration:
    #                            h_rain=h_rain
    #                        else:
    #                            h_rain=0
                    IMc[i][j].append(h_rain)
                    Hl[i][j].append(H/3)
            else:
                    c=0
                    H=0
                    gam=0
                    minutes=[60,120,240,720]
                    IMc[i][j].append(0)
                    Hl[i][j].append(0)
                    FS_rain[i][j].append(0)

    return   FS_rain,IMc,H,c,gam,Hl,minutes

#%%

'''************************************************************'''
'''                    SLOPE HAZARD CURVE RAIN                     '''
'''************************************************************'''

def HazardSlopeRain(seismic_data,Seismic_list,p,minutes,df_RainCoeff,df_Latitude,df_Longitude) :

    '''List of coefficents for rain'''
    RainCoeff=df_RainCoeff.values.tolist()
    RainLat=list(df_Latitude[0])
    RainLon=list(df_Longitude[0])

    closest_x = take_closest(RainLon, float(p["Longitude"]))
    closest_y = take_closest(RainLat, float(p["Latitude"]))
    ix=RainLon.index(closest_x)
    iy=RainLat.index(closest_y)
    CoefRain=RainCoeff[iy][ix]
    '''
    plt.figure()
    '''
    Tr=[10,20,50,100,200]
    MAF=[]
    IM=[]

    base=[]
    expon=[]
    for i in range(len(Tr)):
        MAF.append(1/Tr[i])

    for z in range(len(minutes)):
        IM.append([])
        for ii in range(len(Tr)):
            IM[z].append(CoefRain*34.5*((Tr[ii])**0.19)/(minutes[z]/60)**0.56)

        x = np.array(IM[z])
        y = np.array(MAF)

        def exponential(x,a,k,b):
            return a*np.exp(-x*k)

        popt, pcov = scipy.optimize.curve_fit(exponential, x, y, p0=[5,0.05, 1])
        a=(popt[0])
        k=(popt[1])
        base.append(a)
        expon.append(k)

        '''
        #Plot the 4point curve
        plt.subplot(int(str("22"+str(z+1))))
        x=np.linspace(x[0],x[len(x)-1],20)
        plt.plot(IM[z],MAF,'bo',x,a*np.exp(-x*k),'r')
        plt.ylabel('MAF [1/y]')
        plt.xlabel('Rainfall [mm/h]')
    plt.suptitle("Rainfall Hazard - Asset [Long:"+str(round(float(p["Longitude"]),2))+";Lat:"+str(round(float(p["Latitude"]),2))+"]")
    '''
    return base,expon

'''************************************************************'''
'''                  LANDSLIDE HAZARD CURVE RAIN                   '''
'''************************************************************'''

def LanslideHazardCurveRain(base,expon,Slope,IMc,c,gam,H,p,d_soilDB,dx,dy,altitude,Altitude_SN,Hl,minutes):
    '''
    plt.figure()
    '''
    HAZ=[]
    h_rain=[]
    d_rain=[]
    for zz in range(len(base)):
        d_rain.append([])
        h_rain.append([])
        HAZ.append([])                                  #Escluding values not realistic and Putting in one single list
        for i in range(4):
            for j in range(len(Slope)):
                IM=IMc[i][j][zz]
                h=base[zz]*np.exp(-IM*expon[zz])
                if h<1*10**(-20):
                    HAZ[zz].append(0)
                else:
                    if IM==0:
                        HAZ[zz].append(0)
                    else:
                        HAZ[zz].append(h)

        Height=[]
        for i in range(4):
            for j in range(len(Slope)):
                Height.append(Hl[i][j][zz])
        IM=[]                                   #Putting in one single list
        for i in range(4):
            for j in range(len(Slope)):
                IM.append(IMc[i][j][zz])

        hazard_long=pd.DataFrame(data={'Hazard':HAZ[zz],'Height':Height,'IntensityMeasure':IM})
        hazard=hazard_long.sort_values('Hazard',ascending=False)
        Resolution= pd.DataFrame({'dx':[dx],'dy':[dy]})
        Asset_Altitude=altitude
        Asset_Location=p
        SoilProperties_DB=d_soilDB

        AltitudeSN=pd.DataFrame(data={'0':Altitude_SN[0],'1':Altitude_SN[1],'2':Altitude_SN[2],'3':Altitude_SN[3],'4':Altitude_SN[4],'5':Altitude_SN[5],'6':Altitude_SN[6]})

        #creating list from DB
        h=list(hazard['Hazard'])
        d=list(hazard['Height'])
        im=list(hazard['IntensityMeasure'])

        #discarding all nil values
        hazard=[]
        depth=[]
        intensity=[]
        for i in range(len(h)):
            if h[i]==0:
                i
            else:
                hazard.append(h[i])
                depth.append(d[i])
                intensity.append(im[i])

                h_rain[zz].append(h[i])
                d_rain[zz].append(d[i])

        IndBaseRainHazard=[]
        IndExponentRainHazard=[]
        IndThirdCoeffRainHazard=[]
        for i in range(len(h_rain)):
            IndBaseRainHazard.append(None)
            IndExponentRainHazard.append(None)
            IndThirdCoeffRainHazard.append(None)

    return AltitudeSN,Resolution,Asset_Altitude,Asset_Location,SoilProperties_DB,h_rain,d_rain,IndBaseRainHazard,IndExponentRainHazard,IndThirdCoeffRainHazard

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def hazard_landslide(Longitude, Latitude):


    #TEST - check altitude at a known location
    image, lat_image, lon_image = load_raster('./data/Flood/SRTM_TR_30_shaped.tif')
    # Acwa Power Kırıkkale Doğalgaz Santrali = item 222  --> Altitude from Google is 682 m @ Lat 39.884866 Lon 33.40826
    #p=df_power_plant[["Latitude","Longitude"]].loc[[222]]
    p=pd.DataFrame({'Longitude':[Longitude],
                    'Latitude':[Latitude]})
    #lat_image1.reverse()
    #val, ix, iy =getRasterValueAtCoo(image,lon_image,lat_image,float(p["Longitude"]),float(p["Latitude"]))
    altitude,Altitude_SN, Slope,SlopeOrientation,dx,dy = getContextMorphology(image, lon_image,lat_image,float(p["Longitude"]),float(p["Latitude"]), n=3, l_min = 500)
    #altitude, ix, iy, SN_section, WE_section, SN_slope, WE_slope, SNmorph, WEmorph, SlopeOrientation = getContextMorphology(image, lon_image,lat_image,float(p["Longitude"]),float(p["Latitude"]), n=3, l_min = 250)
    # Read HWSD data.
    df_HWSD=pd.read_excel('./data/Landslide/HWSD_DATA.xlsx')#directory / 'HWSD_DATA.xlsx')

    df_RainCoeff=pd.read_excel('./data/Landslide/wc2.1_5m_prec/Rainfall_Data.xlsx')#directory / 'wc2.1_5m_prec/Rainfall_Data.xlsx')
    df_Latitude=pd.read_excel('./data/Landslide/wc2.1_5m_prec/Latitude_Rainfall_Data.xlsx')#directory / 'wc2.1_5m_prec/Latitude_Rainfall_Data.xlsx')
    df_Longitude=pd.read_excel('./data/Landslide/wc2.1_5m_prec/Longitude_Rainfall_Data.xlsx')#directory / 'wc2.1_5m_prec/Longitude_Rainfall_Data.xlsx')

    imageHWSD, lat_imageHWSD, lon_imageHWSD = load_raster('./data/Flood/HWSD_TR_shaped.tif')
    val, ix, iy =getRasterValueAtCoo(imageHWSD,lon_imageHWSD,lat_imageHWSD,float(p["Longitude"]),float(p["Latitude"]))

    ''' IF THE ASSET IS TO CLOSE TO THE SEE THIS LOOP TAKE THE FIRST POINT IN LAND 
    (AND NOT THE NEAREST, THAT CAN LEAD TO ERROR BACAUSE OF THE NO DATA IN THE SEE'''
    ccc=np.linspace(0.001,0.1,10)
    val_new=[0,0,0]
    for j in range(len(ccc)):
        if sum(val_new)==0:
            cc=ccc[j]
            change=[[-cc,-cc],[+cc,+cc],[-cc,+cc],[+cc,-cc]]
            val_new=[]
            if val==0:
                for i in range(3):
                    p_new=pd.DataFrame({'Longitude':[Longitude+change[i][0]],
                            'Latitude':[Latitude+change[i][1]]})
                    val, ix, iy =getRasterValueAtCoo(imageHWSD,lon_imageHWSD,lat_imageHWSD,float(p_new["Longitude"]),float(p_new["Latitude"]))
                    val_new.append(val)
    for i in range(len(val_new)):
        if val_new[i]!=0:
            val=val_new[i]
    '''                                                                         '''

    soil = getSoil(df_HWSD, val, full_HWSD_row=False)

    d_soilDB=pd.read_excel('./data/Landslide/Soil_Rock_Charateristics_DB.xlsx')#directory / 'Soil_Rock_Charateristics_DB.xlsx')
    Soil=d_soilDB.iloc[int(soil)]

    if soil!=0:
        IMc,Hl=Fragility_PvsEQ(Slope,Soil,soil)
        '''SEISMIC DATA (PGA) FOR 50 Years (TURKEY)'''
        seismic_data=pd.read_csv('./data/Landslide/EQ_Data_50years_PGA.csv')#directory / "EQ_Data_50years_PGA.csv")
        # Si vede che con l'aumentare dell'indice aumenta la LONGITUDINE !!! 
        # ma va scelta anche la latitudine perche è da sud a nord!!
        Seismic_list=seismic_data.values.tolist()
        AssetPGA,BaseEQHazard,ExponentEQHazard=HazardSlopeEQ(seismic_data,Seismic_list,p)
        AltitudeSN,Resolution,Asset_Altitude,Asset_Location,SoilProperties_DB,hazard,h_eq,d_eq,IndBaseEQHazard,IndExponentEQHazard=LanslideHazardCurveEQ(BaseEQHazard,ExponentEQHazard,Slope,IMc,Hl,p,d_soilDB,dx,dy,altitude,Altitude_SN)

        FS_rain,IMc,H,c,gam,Hl,minutes=Fragility_PvsRain(Slope,Soil,soil)
        BaseRainHazard,ExponentRainHazard=HazardSlopeRain(seismic_data,Seismic_list,p,minutes,df_RainCoeff,df_Latitude,df_Longitude)
        AltitudeSN,Resolution,Asset_Altitude,Asset_Location,SoilProperties_DB,h_rain,d_rain,IndBaseRainHazard,IndExponentRainHazard,IndThirdCoeffRainHazard=LanslideHazardCurveRain(BaseRainHazard,ExponentRainHazard,Slope,IMc,c,gam,H,p,d_soilDB,dx,dy,altitude,Altitude_SN,Hl,minutes)

        '''************************************************************'''
        '''                    SLOPE HAZARD CURVE EQ+RAIN                    '''
        '''************************************************************'''



        def truncate(n, decimals=0):
            multiplier = 10 ** decimals
            return int(n * multiplier) / multiplier

        Landslide_Hazard=pd.DataFrame(data={'Earthquake':[h_eq,d_eq],
                                            'rain60':[h_rain[0],d_rain[0]],
                                            'rain120':[h_rain[1],d_rain[1]],
                                            'rain240':[h_rain[2],d_rain[2]],
                                            'rain720':[h_rain[3],d_rain[3]]
                                            })
        if h_eq == [] and h_rain[3] == []:
            IndBaseEQHazard = None
            IndExponentEQHazard =None
        else:
            maxlen = max(len(h_eq),len(h_rain[3]))
            if maxlen < 3:
                IndBaseEQHazard = None
                IndExponentEQHazard =None
            else:
                for i in range(len(h_rain)):
                    h_rain[i] = h_rain[i] + [0] * int(maxlen-len(h_rain[i]))
                h_eq = h_eq + [0] * int(maxlen-len(h_eq))

                haz_val = []
                for i in range(maxlen):
                    haz_val.append(max(h_eq[i],h_rain[0][i],h_rain[1][i],h_rain[2][i],h_rain[3][i]))

                if len(d_eq) > len(d_rain[3]):
                    x = d_eq
                else: x = d_rain[3]

                y = haz_val


                def exponential(x,a,k,b):
                    return a*np.exp(-x*k)
                popt, pcov = scipy.optimize.curve_fit(exponential, x, y, p0=[100,1, 1])
                IndBaseEQHazard=(popt[0])
                IndExponentEQHazard=(popt[1])
    else:
        IndBaseEQHazard = None
        IndExponentEQHazard =None
        d_eq = None
    return [IndBaseEQHazard,IndExponentEQHazard],d_eq#,[IndBaseRainHazard,IndExponentRainHazard,IndThirdCoeffRainHazard],DoubleCurve,FunctionEquation,Landslide_Hazard




#logger.debug("Module loaded.")

if __name__ == '__main__':
    lon,lat = 37,37
    df = hazard_landslide(lon,lat)
    print("a")