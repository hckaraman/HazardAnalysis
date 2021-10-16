#############
import logging
logger = logging.getLogger("afad_backend.%s" % __name__)
logger.setLevel(logging.DEBUG) ## CHANGE THIS FOR LOG LEVEL
#logger.debug("Loading module...")
######################


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
from data.constants import *

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
    #    d_area={'NE':[dx*dy,dx*dy,dx*dy,dx*dy,dx*dy,dx*dy], 
#            'SE':[dx*dy,dx*dy,dx*dy,dx*dy,dx*dy,dx*dy],
#            'SN':[dx*dy,dx*dy,dx*dy,dx*dy,dx*dy,dx*dy],
#            'WE':[dx*dy,dx*dy,dx*dy,dx*dy,dx*dy,dx*dy]} 
    
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

#%%

'''************************************************************'''
'''                 FRAGILITY CURVE OF SLOPE EQ                  '''
'''************************************************************'''

def Fragility_PvsEQ(Slope,soil,d_soilDB):
    logger.debug('Fragility_PvsEQ called')
    Slope_list=Slope.values.tolist()
    slopes=[]
    for i in range(len(Slope_list)):
        for j in range(len(Slope_list[0])):
            slopes.append(Slope_list[i][j])
        
    if soil==str(''):
        k=0
    else:
        k=soil
    
    if k!=0:
        measure=0                
    else:
        Soil=d_soilDB.iloc[int(0)]
        if max(slopes)<45:
            measure=0         
        else:
            measure=Soil[3]            
    return measure
        
#%%

'''************************************************************'''
'''                    SLOPE HAZARD CURVE EQ                      '''
'''************************************************************'''

def HazardSlopeEQ(seismic_data,Seismic_list,p,measure):
        
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
    
    if measure!=0:
        #interpolating the hazard curve
        x = np.array(IM)
        y = np.array(MAF)
        def exponential(x,a,k,b):
            return a*np.exp(-x*k)
        popt, pcov = scipy.optimize.curve_fit(exponential, x, y, p0=[0.5,5, 1])
        a_eq=(popt[0])
        k_eq=(popt[1])
        
        Intensity = a_eq*np.exp(-measure*k_eq)

        
        a = Intensity*0.00844228532124152
        k = 0.469
        
        
    else:
        a="NO HAZARD"
        k="NO HAZARD"
    return a,k

def hazard_rockfall(Longitude,Latitude):

    image, lat_image, lon_image = load_raster('./data/Flood/SRTM_TR_30_shaped.tif')

    p=pd.DataFrame({'Longitude':[Longitude],
                    'Latitude':[Latitude]})
   
    altitude,Altitude_SN, Slope,SlopeOrientation,dx,dy = getContextMorphology(image, lon_image,lat_image,float(p["Longitude"]),float(p["Latitude"]), n=3, l_min = 500)
    df_HWSD=pd.read_excel('./data/Rockfall/HWSD_DATA.xlsx')#DATA / 'Rockfall/HWSD_DATA.xlsx')
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
    logger.debug(val)
    soil = getSoil(df_HWSD, val)
    d_soilDB=pd.read_excel('./data/Rockfall/Soil_Rock_Charateristics_DB.xlsx')#DATA / 'Rockfall/Soil_Rock_Charateristics_DB.xlsx')
    seismic_data=pd.read_csv('./data/Rockfall/EQ_Data_50years_PGA.csv')#DATA / "Rockfall/EQ_Data_50years_PGA.csv")
    Seismic_list=seismic_data.values.tolist()    
    
    measure = Fragility_PvsEQ(Slope,soil,d_soilDB)
    
    BaseHazard,ExponentHazard=HazardSlopeEQ(seismic_data,Seismic_list,p,measure)
    
    return [BaseHazard,ExponentHazard]

#BaseHazard,ExponentHazard = hazard_rockfall(29.76,40.75)
if __name__ == '__main__':
    df = hazard_rockfall(41.5257998166669,40.7867457566693)
    print("a")

# logger.debug("Module loaded.")