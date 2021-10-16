# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 14:21:49 2020

@author: RVA04
"""
#############
import logging
logger = logging.getLogger("afad_backend.%s" % __name__)
logger.setLevel(logging.DEBUG) ## CHANGE THIS FOR LOG LEVEL
#logger.debug("Loading module...")
######################

import geopandas as gpd
import pandas as pd
import math
from bisect import bisect_left
from shapely.geometry import Point
import matplotlib.pyplot as plt
import rasterio
import numpy as np
import scipy.optimize
from scipy.stats import norm
import os
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
    ix=x.index(closest_x)
    iy=y.index(closest_y)
    if out_index: return raster[-iy][ix], ix, iy
    return raster[-iy][ix] 

def getContextMorphology(raster, x, y, x_ref, y_ref, crs='epsg:4326', n=20, l_min = 250): # Altitude and Slope evaluation
    """
    Get sections at reference point
    n: number of points in each direction
    l_min: minimum section length
    """
    altitude, ix, iy = getRasterValueAtCoo(raster, x, y, x_ref, y_ref, crs=crs, out_index=True)
 
    dy=great_circle_vec(x[ix],y[iy],x[ix],y[iy+1]) # distance in m along y
    dx=great_circle_vec(x[ix],y[iy],x[ix+1],y[iy]) # distance in m along x
    
    n1 = math.ceil(max(l_min/(dx*2),l_min/(dy*2)))
    if n1>n: n=n1
    
    Altitude_SN=[]
    for i in range(ix-n,ix+n+1):
        Altitude_SN.append([])
        ii=len(Altitude_SN)-1
        for j in range(-iy-n,-iy+n+1):
            Altitude_SN[ii].append(raster[j][i])
            
    df_MeasureSN=pd.DataFrame(data=Altitude_SN)
    
    return df_MeasureSN,dx,dy
                
def twoDF(p):
    logger.debug('twoDF called')
    # Qavg_image, lat_image, lon_image = load_raster(config.get('PATH','AvgDischargeMap'))
    Qavg_image, lat_image, lon_image = load_raster('./data/Flood/logQavg.tif')
    # Qvar_image, lat_image, lon_image = load_raster(config.get('PATH','VarDischargeMap'))
    Qvar_image, lat_image, lon_image = load_raster('./data/Flood/logQvar.tif')

    Altitude_image, lat_alt, lon_alt = load_raster('./data/Flood/SRTM_TR_30_shaped.tif')
    
    qavg,ix,iy=getRasterValueAtCoo(Qavg_image, lon_image, lat_image,float(p["Longitude"]),float(p["Latitude"]), crs='epsg:4326', out_index=True)
    lon_river=lon_image[ix]
    lat_river=lat_image[iy]
    qvar,ix,iy=getRasterValueAtCoo(Qvar_image, lon_image, lat_image,float(p["Longitude"]),float(p["Latitude"]), crs='epsg:4326', out_index=True)
    
    AltitudeAsset,ix,iy=getRasterValueAtCoo(Altitude_image, lon_alt, lat_alt,float(p["Longitude"]),float(p["Latitude"]), crs='epsg:4326', out_index=True)
    lon_asset=lon_alt[ix]
    lat_asset=lat_alt[iy]
    
    distance = great_circle_vec(lat_asset, lon_asset, lat_river, lon_river, earth_radius=6371009)     
    
    df_AltitudeSN,dx,dy = getContextMorphology(Altitude_image, lon_alt, lat_alt,float(p["Longitude"]),float(p["Latitude"]), crs='epsg:4326', n=20, l_min = 250)
    df_DischargeSN,dxr,dyr = getContextMorphology(Qavg_image, lon_image, lat_image,float(p["Longitude"]),float(p["Latitude"]), crs='epsg:4326', n=20, l_min = 250)
    df_DischargeVarSN,dxrv,dyrv = getContextMorphology(Qvar_image, lon_image, lat_image,float(p["Longitude"]),float(p["Latitude"]), crs='epsg:4326', n=20, l_min = 250)
    
    DischargeSN=df_DischargeSN.values.tolist()
    DischargeVarSN=df_DischargeVarSN.values.tolist()
    q=[]
    qv=[]
    for i in range(len(DischargeSN)):
        q.append([])
        qv.append([])
        for j in range(len(DischargeSN[0])):
            if DischargeSN[i][j]==0:
                q[i].append(0)
                qv[i].append(0)
            else:                
                q[i].append(10**DischargeSN[i][j])
                qv[i].append(10**DischargeVarSN[i][j])
    
    df_DischargeSN=pd.DataFrame(q)
    df_DischargeVarSN=pd.DataFrame(qv) 
    return df_AltitudeSN,df_DischargeSN,df_DischargeVarSN,dx,dy,dxr,dyr,distance,AltitudeAsset

def riverdirection(df_DischargeSN):
    s=abs(sum(df_DischargeSN[0].values.tolist()))
    n=abs(sum(df_DischargeSN[int(len(df_DischargeSN)-1)].values.tolist()))
    alt=df_DischargeSN.values.tolist()
    w=abs(sum(alt[0]))
    e=abs(sum(alt[int(len(df_DischargeSN)-1)]))
    if w+e-s-n>0:
        direction='WE'
        if e*w==0:
            direction='unk'
    else:
        if w+e-s-n<0:
            direction='SN'
        if s*n==0:
            direction='unk'
        else:
            if w+e-s-n==0:
                direction='unk'
    return direction

def PlotWEdirection(p,twoDF_result):
    logger.debug('PlotWEdirection called')
    df_AltitudeSN,df_DischargeSN,df_DischargeVarSN,dx,dy,dxr,dyr,distance,AltitudeAsset = twoDF_result
    #altimetria
    a=df_AltitudeSN[int(len(df_AltitudeSN)/2)].values.tolist()
    al=int(len(df_AltitudeSN))*dx
    xa=np.linspace(-al/2,al/2,len(df_AltitudeSN[int(len(df_AltitudeSN)/2)]))
    #river position
    q=df_DischargeSN[int(len(df_AltitudeSN)/2)].values.tolist()
    qv=df_DischargeVarSN[int(len(df_AltitudeSN)/2)].values.tolist()
    aq=int(len(df_DischargeSN))*dxr
    xq=np.linspace(-aq/2,aq/2,len(df_DischargeSN[int(len(df_DischargeSN)/2)]))
    q_4plot=[]
    for i in range(len(q)):
        if q[i]==0:
            q[i]=0
            qv[i]=0
            q_4plot.append(None)
            xq[i]=None
        else:
            q_4plot.append(a[i])
            xq[i]=xq[i]#-distance
            Riverposition=xq[i]
    Riverposition = q.index(max(q))           
    OneRiver=q.index(max(q))  
    for i in range(len(q)):
        if i==OneRiver:
           q[i]=q[i]
           q_4plot[i]=q_4plot[i]
        else:
           q[i]=None
           q_4plot[i]=None  
      
    #asset position
    x=0
    y=a[int(len(df_AltitudeSN)/2)]        
        
    qa=q[OneRiver]
    qv=qv[q.index(qa)]
      
    x=[]
    for i in range(len(xa)):
        x.append(xa[i])
            
    if a.index(min(a))<10:
        arrai=np.linspace(0,2,9)
    else:
        if a.index(min(a))>(len(a)-1)-10:
            arrai=np.linspace(len(a)-9,len(a)-1,9)
        else: 
            arrai=np.linspace(a.index(min(a))-5,a.index(min(a))+5,9) 
    minimum=[]
    for i in range(len(arrai)):
        minimum.append(a[int(arrai[i])])
    AltitudeRiver=min(minimum)

    IndexRiver=a.index(AltitudeRiver)
    IndexAsset=int((len(xa)-1)/2)
    
    return qa,qv,x,a,AltitudeRiver,IndexRiver,IndexAsset
    
def PlotSNdirection(p,twoDF_result):
    logger.debug('PlotSNdirection called')
    df_AltitudeSN,df_DischargeSN,df_DischargeVarSN,dx,dy,dxr,dyr,distance,AltitudeAsset = twoDF_result
    #altimetria
    a=df_AltitudeSN.values.tolist()
    a=a[int(len(df_AltitudeSN)/2)]
    al=int(len(df_AltitudeSN))*dy
    xa=np.linspace(-al/2,al/2,len(df_AltitudeSN[int(len(df_AltitudeSN)/2)]))
    #river position
    q=df_DischargeSN.values.tolist()
    qv=df_DischargeVarSN.values.tolist()
    q=q[int(len(df_AltitudeSN)/2)]
    qv=qv[int(len(df_AltitudeSN)/2)]
    aq=int(len(df_DischargeSN))*dyr
    xq=np.linspace(-aq/2,aq/2,len(df_DischargeSN[int(len(df_DischargeSN)/2)]))
    q_4plot=[]
    for i in range(len(q)):
        if q[i]==0:
            q[i]=0
            qv[i]=0
            q_4plot.append(None)
            xq[i]=None
        else:
            q_4plot.append(a[i])
            xq[i]=xq[i]#-distance
            Riverposition=xq[i]
    Riverposition = q.index(max(q))
    OneRiver=q.index(max(q))  
    for i in range(len(q)):
        if i==OneRiver:
           q[i]=q[i]
           q_4plot[i]=q_4plot[i]
        else:
           q[i]=None
           q_4plot[i]=None       
    #asset position
    x=0
    y=a[int(len(df_AltitudeSN)/2)]        
    
    qa=q[OneRiver]
    qv=qv[q.index(qa)]
    
    x=[]
    for i in range(len(xa)):
        x.append(xa[i])
            
    if a.index(min(a))<10:
        arrai=np.linspace(0,2,9)
    else:
        if a.index(min(a))>(len(a)-1)-10:
            arrai=np.linspace(len(a)-9,len(a)-1,9)
        else: 
            arrai=np.linspace(a.index(min(a))-5,a.index(min(a))+5,9) 
    minimum=[]
    for i in range(len(arrai)):
        minimum.append(a[int(arrai[i])])
    AltitudeRiver=min(minimum)

    IndexRiver=a.index(AltitudeRiver)
    IndexAsset=int((len(xa)-1)/2)
    
    return qa,qv,x,a,AltitudeRiver,IndexRiver,IndexAsset

def FloodingHazard(p):
    logger.debug('FloodingHazard called')
    twoDF_result = twoDF(p)
    df_AltitudeSN,df_DischargeSN,df_DischargeVarSN,dx,dy,dxr,dyr,distance,AltitudeAsset = twoDF_result
    RiverOrientation=riverdirection(df_DischargeSN)
    
    ###CHECK IF THE STREAM EXIST, IF NOT THE COMPUTATION IS NOT PERFORMED
    DischargeSN=df_DischargeSN.values.tolist()
    partSN=abs(sum(DischargeSN[int(len(df_AltitudeSN)/2)]))             #hor
    partWE=abs(sum(df_DischargeSN[int(len(df_AltitudeSN)/2)]))          #vert
                
    if partSN==0 and partWE==0:
        Qaverage, Qvariance, x_profile, y_profile, AltitudeRiver,IndexRiver,IndexAsset = PlotSNdirection(p,twoDF_result) 
        Qaverage, Qvariance, x_profile, y_profile, AltitudeRiver,IndexRiver,IndexAsset = PlotWEdirection(p,twoDF_result)   
        logger.debug('NO STREAM - HAZARD NOT COMPUTED')
        
        
    else:
        if RiverOrientation=='WE':
            Qaverage, Qvariance, x_profile, y_profile, AltitudeRiver,IndexRiver,IndexAsset = PlotSNdirection(p,twoDF_result)   
        else:
            if RiverOrientation=='SN':
                Qaverage, Qvariance, x_profile, y_profile, AltitudeRiver,IndexRiver,IndexAsset = PlotWEdirection(p,twoDF_result)
            else:
                if RiverOrientation=='unk':
                    if partSN>partWE:
                        Qaverage, Qvariance, x_profile, y_profile, AltitudeRiver,IndexRiver,IndexAsset = PlotSNdirection(p,twoDF_result)
                    else:
                        Qaverage, Qvariance, x_profile, y_profile, AltitudeRiver,IndexRiver,IndexAsset = PlotWEdirection(p,twoDF_result)
                        
        def DischargeFrequecies(Qvariance,Qaverage):
            def Funct(sigma):
                return abs((0.5*(1+math.erf((Qvariance*Qaverage-Qaverage)/(sigma*2**0.5)))) - 0.9)
            logger.debug('DischargeFrequecies optimization Funct')
            sigma = scipy.optimize.minimize(Funct, 60, method='nelder-mead',options={'xatol': 1e-14, 'disp': False})
            sigma=float(sigma.x)
            return sigma
            
        sigma = DischargeFrequecies(Qvariance,Qaverage)    
            
        def MAFvsQ(sigma,Qaverage):
            Prob=[0.0060,0.0311,0.1877,0.2895]
            Tr=[]
            MAF=[]
            Q_haz=[]
            for i in range(len(Prob)):
                Tr.append(-15/math.log(1-Prob[i]))
                MAF.append(-math.log(1-Prob[i])/15)
                
                def Funct(x):
                    return abs((0.5*(1+math.erf((x-Qaverage)/(sigma*2**0.5)))) - (1-(Prob[i])))

                logger.debug('MAFvsQ optimization Funct %d of %d' % (i,len(Prob)))
                x = scipy.optimize.minimize(Funct, Qaverage, method='nelder-mead',options={'xatol': 1e-14, 'disp': False})
                x=float(x.x)
                Q_haz.append(x)
                
            #interpolating the hazard curve
            xx = np.array(Q_haz)
            yy = np.array(MAF)
            basefirst=0.3
            expfirst=0.03 #0.0067*math.log(Qaverage)-0.0258
            def exponential(xx,a,k,b):
                return a*np.exp(-xx*k)
            logger.debug('MAFvsQ optimization exponential')
            popt, pcov = scipy.optimize.curve_fit(exponential, xx, yy, p0=[basefirst,expfirst, 1])
            a=(popt[0])
            k=(popt[1]) 
            return Q_haz,a,k
        
        
        Q_haz,HazardBase,HazardExponent = MAFvsQ(sigma,Qaverage)
        
        x_profile=x_profile-x_profile[0] 
        
        def InducedHazard(Qaverage,AltitudeAsset,AltitudeRiver,HazardBase,HazardExponent,x_profile,y_profile,IndexRiver,IndexAsset):
            vel=Qaverage**0.1

            ia=IndexAsset
            ir=IndexRiver
            new_y=y_profile
            if ia<ir:
                y_profile_shrink=(y_profile[ia:ir])
                peak=max(y_profile_shrink)
                ip=y_profile_shrink.index(max(y_profile_shrink))
                for i in range(len(new_y)):
                    if peak>AltitudeAsset:
                        if new_y[i]>=peak:
                            new_y[i]=9999
                        else:
                            if x_profile[i]<x_profile[ip+ia]:
                                new_y[i]=9999
            else: 
                if ia>ir:
                    y_profile_shrink=(y_profile[ir:ia])
                    peak=max(y_profile_shrink)
                    ip=y_profile_shrink.index(max(y_profile_shrink))
                    # print(peak,AltitudeAsset)
                    for i in range(len(new_y)):
                        if peak>AltitudeAsset:
                            if new_y[i]>=peak:
                                new_y[i]=9999
                            else:
                                if x_profile[i]>x_profile[ip+ir]:
                                    new_y[i]=9999
                else:
                    peak=y_profile[ir]

                    logger.warn('ASSET LOCATED IN THE RIVER')

            if peak>AltitudeAsset:
                Hvector=np.linspace(AltitudeRiver+0.01,AltitudeAsset+5.01,100)
                maf=[]
                q_4df=[]
                for j in range(len(Hvector)):
                    trap=[]
                    for i in range(len(x_profile)-1):
                        if Hvector[j]<peak:
                            if Hvector[j]>new_y[i+1]-0.001 and Hvector[j]>new_y[i]-0.001:
                                trap.append((x_profile[i+1]-x_profile[i])*((Hvector[j]-new_y[i+1])+(Hvector[j]-new_y[i]))/2)
                            else:trap.append(0)
                        else: 
                            if Hvector[j]>y_profile[i+1]-0.001 and Hvector[j]>y_profile[i]-0.001:
                                trap.append((x_profile[i+1]-x_profile[i])*((Hvector[j]-y_profile[i+1])+(Hvector[j]-y_profile[i]))/2)
                            else:trap.append(0)
                    discharge=sum(trap)*vel
                    q_4df.append(discharge)
                    maf.append(HazardBase*np.exp(-discharge*HazardExponent))                
            else:
                Hvector=np.linspace(AltitudeRiver+0.01,AltitudeAsset+5.01,100)
                maf=[]
                q_4df=[]
                for j in range(len(Hvector)):
                    trap=[]
                    for i in range(len(x_profile)-1):
                        if Hvector[j]>y_profile[i+1]-0.001 and Hvector[j]>y_profile[i]-0.001:
                            trap.append((x_profile[i+1]-x_profile[i])*((Hvector[j]-y_profile[i+1])+(Hvector[j]-y_profile[i]))/2)
                        else:trap.append(0)
                    discharge=sum(trap)*vel
                    q_4df.append(discharge)
                    maf.append(HazardBase*np.exp(-discharge*HazardExponent))
            
            MAF=[]
            height=[]
            masl=[]
            for i in range(len(maf)):
                # if maf[i]>10**-15:
                if maf[i] != 0:
                    masl.append(Hvector[i])
                    MAF.append(maf[i])
                    if Hvector[i]>AltitudeAsset:
                        if peak>Hvector[i]:
                            height.append(None)
                        else:
                            height.append(Hvector[i]-AltitudeAsset)
                    else:
                        height.append(None) 
                else:
                    masl.append(Hvector[i])
                    MAF.append(None)
                    if Hvector[i]>AltitudeAsset:
                        if peak>Hvector[i]:
                            height.append(None)
                        else:
                            height.append(Hvector[i]-AltitudeAsset)
                    else:
                        height.append(None)  
            
            
            kk=[]
            for i in range(len(height)):
                if height[i]==None:
                    height[i]=0
                if MAF[i]==None:
                    MAF[i]=0
                kk.append(height[i]*MAF[i])
            k=sum(kk)
        
            if k!=0:
                MAFint=[]
                heigthint=[]
                for i in range(len(height)):
                    if MAF[i]!=0:
                        if height[i]!=0:
                            MAFint.append(MAF[i])
                            heigthint.append(height[i])

                    x = np.array(heigthint)
                    y = np.array(MAFint)

                    if len(x) > 2:
                        def exponential(x,a,k,b):
                            return a*np.exp(-x*k)
                        #logger.debug("optimizing exponential function")
                        popt, pcov = scipy.optimize.curve_fit(exponential, x, y, p0=[1,1, 1])

                        IndHazardBase=(popt[0])
                        IndHazardExponent=(popt[1])

                    elif len(x) == 2:
                        IndHazardBase = max(y)
                        IndHazardExponent = 0
                    elif len(x) == 1:
                        IndHazardBase = y[0]
                        IndHazardExponent = 0
                
            if k==0:
                IndHazardBase=str('NaN')
                IndHazardExponent=str('NaN')
               
            df_hazard=pd.DataFrame(data={'Meter above sea level [masl]':masl,'Height [m]':height,'MAF [1/y]':MAF,'Q [m3/s]':q_4df})
                 
            return df_hazard,IndHazardBase,IndHazardExponent
        
        df_Hazard,IndHazardBase,IndHazardExponent = InducedHazard(Qaverage,AltitudeAsset,AltitudeRiver,HazardBase,HazardExponent,x_profile,y_profile,IndexRiver,IndexAsset)
        
    #    AssetPosition=p
    if partSN==0 and partWE==0:
        return "NO_HAZARD","NO_HAZARD","NO_HAZARD"
    else:    
        return IndHazardBase,IndHazardExponent,df_Hazard
    

    
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def hazard_flood(Longitude,Latitude):
#    directory=os.getcwd()+"/"
    # print(directory)
    p=pd.DataFrame({'Longitude':[Longitude],
                'Latitude':[Latitude]})  
    IndHazardBase,IndHazardExponent,df_hazard = FloodingHazard(p)
    return  [IndHazardBase,IndHazardExponent,df_hazard]

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''



#logger.debug("Module loaded.")

if __name__ == '__main__':
    lon,lat = 37,37
    df = hazard_flood(lon,lat)
    print("a")