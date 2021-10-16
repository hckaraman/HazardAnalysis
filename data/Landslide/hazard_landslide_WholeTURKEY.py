# -*- coding: utf-8 -*-
"""
Created on Wed May 27 13:38:47 2020

@author: RVA04
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 16:03:12 2020

@author: RVA04
"""

# SRTM (Risoluzione 90 m x 80 m)
# From RASTER (e.g. TIFF) to grid
import rasterio
import numpy as np
#from rasterio.features import shapes
#import descartes
from bisect import bisect_left
import math
import scipy.optimize 
import matplotlib.pyplot as plt
import time
import geopandas as gpd
import pandas as pd
import os

'''************************************************************'''
'''                   MORPHOLOGY DEFINITION                    '''
'''************************************************************'''

# LOAD RASTER IMAGE
def load_raster(path):
    with rasterio.Env():
        with rasterio.open(path) as src:  
            image = src.read(1)  ## first and only band  
    #src.crs
    #src.transform
    #src.bounds
    lat_image = [src.xy(i,0)[1] for i in range(src.height)]
    lon_image = [src.xy(0,i)[0] for i in range(src.width)]
    
    lat_image.reverse() # ordered list from max to min
    
    return image, lat_image, lon_image
    
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
#    print("Reference point: {}, {}".format(x_ref,y_ref))
#    print("Closest point: {}, {}".format(closest_x, closest_y))
    ix=x.index(closest_x)
    iy=y.index(closest_y)
    if out_index: return raster[-iy][ix], ix, iy
    return raster[-iy][ix] 

def getContextMorphology(raster, x, y, x_ref, y_ref, crs='epsg:4326', n=1, l_min = 250): # Altitude and Slope evaluation
    """
    Get sections at reference point
    n: number of points in each direction
    l_min: minimum section length
    """
    altitude, ix, iy = getRasterValueAtCoo(raster, x, y, x_ref, y_ref, crs=crs, out_index=True)
    dy=great_circle_vec(x[ix],y[iy],x[ix],y[iy+1]) # distance in m along y
    dx=great_circle_vec(x[ix],y[iy],x[ix+1],y[iy]) # distance in m along x
#    print("dx = {} m, dy = {} m".format(dx,dy))
    
    Altitude_SN=[]
    for i in range(ix-n,ix+n+1):
        Altitude_SN.append([])
        ii=len(Altitude_SN)-1
        for j in range(-iy-n,-iy+n+1):
            Altitude_SN[ii].append(raster[j][i])
    
    midpoint=int((len(Altitude_SN)-1)/2)
    Slope=[]
    Dist=[]
    for i in range(len(Altitude_SN)):
        ii=-n+i
        Dist.append([])
        for j in range(len(Altitude_SN[0])):
            jj=-n+j
            if ii==0:
                dist=dy
            else:
                if jj==0:
                    dist=dx
                else:
                    dist=(dx**2+dy**2)**(1/2)
            Slope.append(abs((Altitude_SN[midpoint][midpoint]-Altitude_SN[midpoint+ii][midpoint+jj])/dist)/math.pi*180)
            Dist[i].append(dist)
        Dip=max(Slope)
    
    return Dip
                 
def getSlopeOrientation(SNsec, EWsec): 
    pass

#TEST - check altitude at a known location
   
def SlopeDataFrame(lat_vect,lon_vect,image,lon_image,lat_image,directory):

    Point_slope=[]
    for i in range(len(lon_vect)):
        Point_slope.append([])
        for j in range(len(lat_vect)):
            jj=len(lat_vect)-j-1
            p_lon=lon_vect[i]
            p_lat=lat_vect[jj]
       
            Dip = getContextMorphology(image, lon_image,lat_image,p_lon,p_lat, n=1, l_min = 250)
            if Dip>60:
                Point_slope[i].append(60)
            else:
                Point_slope[i].append(Dip)
        
        if i==0:
            df_Slope=pd.DataFrame(data=Point_slope[i])
        else:
            df_Slope.insert(int(i),str(i),Point_slope[i])
    return df_Slope,Point_slope

'''************************************************************'''
'''************************************************************'''
'''************************************************************'''

#HWSD  (Risoluzione 900 m x 800 m)

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
    if not list(depth):
        soil=0 #this shoul be water
    else:
        depth=depth.tolist()[0]
        if depth == 10: soil=0
        elif depth == 30: soil=0
        else: soil=depth_row["S_USDA_TEX_CLASS"].tolist()[0]
            
        if full_HWSD_row: return soil, row
    return soil


def SoilDataFrame(lat_vect,lon_vect,lon_imageHWSD,lat_imageHWSD,imageHWSD,df_HWSD):

    Point_soil=[]
    for i in range(len(lon_vect)):
        Point_soil.append([])
        for j in range(len(lat_vect)):
            jj=len(lat_vect)-j-1
            p_lon= take_closest(lon_imageHWSD,lon_vect[i])
            p_lat= take_closest(lat_imageHWSD,lat_vect[jj])
    
            val, ix, iy =getRasterValueAtCoo(imageHWSD,lon_imageHWSD,lat_imageHWSD,p_lon,p_lat)
            soil = getSoil(df_HWSD, val)
            Point_soil[i].append(soil)
        
        if i==0:
            df_Soil=pd.DataFrame(data=Point_soil[i])
        else:
            df_Soil.insert(int(i),str(i),Point_soil[i])
    
    df_Soil.fillna('',inplace=True)
    Point_soil=df_Soil.values.tolist()
    
    return df_Soil,Point_soil



'''************************************************************'''
'''************************************************************'''
'''************************************************************'''


'''************************************************************'''
'''                 FRAGILITY CURVE OF SLOPE EQ                  '''
'''************************************************************'''


def Fragility_PvsEQ(Point_slope,Point_soil,d_soilDB):
    logger.debug('Fragility_PvsEQ called')
    IMc=[]
    FS=[]
    for i in range(len(Point_soil)):
        IMc.append([])
        FS.append([])
        for j in range(len(Point_soil[0])):
            if Point_soil[i][j]==str(''):
                k=0
            else:
                k=Point_soil[i][j]
            Soil=d_soilDB.iloc[int(k)]
            if Soil[0]==0:
                FST=0
                IMT=0
            else:
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
                if Point_slope[j][i]<15:
                    FST=0
                    IMT=0         
                else:
                    dip=Point_slope[j][i]        #degree
                    #Model dimension
                    if dip<15:
                        H=tb15
                        L=lb15/2
                    else:
                        if dip<30:
                            H=tb15-(dip-15)/15*(tb15-tb30)
                            L=lb15-(dip-15)/15*(lb15-lb30)
#                            H=tb30
#                            L=lb30/2
                        else:
                            if dip<45:
                                H=tb30-(dip-30)/15*(tb30-tb45)
                                L=(lb30-(dip-30)/15*(lb30-lb45))/2
                            else:
                                if dip<60:
                                    H=tb45-(dip-45)/15*(tb45-tb60)
                                    L=(lb45-(dip-45)/15*(lb45-lb60))/2
                                else:
#                                    if (dip==60 and dip>60):
                                    H=tb60
                                    L=lb60/2
                    
                    Ip=(H**2+L**2)**(1/2)
                    
                    depthW=H       #m
                    if depthW>H:
                        hw=0
                    else:
                        hw=(H-depthW)/H
                    
                    alfa_dw=math.degrees(math.atan(H/L))
                    psi_dw=dip-alfa_dw
                    
                    psi_up=dip+alfa_dw
                    if psi_up>90:
                        psi_up=90
                    
                    #METODO2
                    C=c*2*Ip
                    Wdw=H*L/2*gam
                    Wup=H*L/2*gam
                    Ndw=Wdw*math.cos(math.radians(psi_dw))
                    Nup=Wup*math.cos(math.radians(psi_up))
                    FST=(Ndw*math.tan(math.radians(phi))+Nup*math.tan(math.radians(phi))+C)/((Wdw*math.sin(math.radians(psi_dw)))+(Wup*math.sin(math.radians(psi_up))))                    
                    IMT=(FST-1)*math.tan(math.radians(dip))/(math.tan(math.radians(phi))*math.tan(math.radians(dip))+1)
            FS[i].append(FST)
            IMc[i].append(IMT)            
#    df_FS=pd.DataFrame(FS)
#    df_IMc=pd.DataFrame(IMc)
    return IMc


'''************************************************************'''
'''                    SLOPE HAZARD CURVE                      '''
'''************************************************************'''


def HazardSlopeEQ(directory,lon_vect,lat_vect,Point_slope,Point_soil,d_soilDB):
    
    '''SEISMIC DATA (PGA) FOR 50 Years (TURKEY)'''
    seismic_data=pd.read_csv(directory+"EQ_Data_50years_PGA.csv") 
    # Si vede che con l'aumentare dell'indice aumenta la LONGITUDINE !!! 
    # ma va scelta anche la latitudine perche è da sud a nord!!
    Seismic_list=seismic_data.values.tolist() 
    
    IMc = Fragility_PvsEQ(Point_slope,Point_soil,d_soilDB)

    '''List of PGAs forthe TR'''
#    PGAs_2_list = seismic_data['2%'].tolist()
#    PGAs_10_list = seismic_data['10%'].tolist()
#    PGAs_50_list = seismic_data['50%'].tolist()
#    PGAs_68_list = seismic_data['68%'].tolist()
#    PGAs_LISTS=[PGAs_2_list,PGAs_10_list,PGAs_50_list,PGAs_68_list]
    
    Percent=[0.02,0.10,0.50,0.68]
    MAF=[]
    for i in range(len(Percent)):
        MAF.append(-np.log(1-Percent[i])/50)
    
    Long=[]
    Lat=[]
    for r in range(len(Seismic_list)):
        Long.append(Seismic_list[r][0])
        Lat.append(Seismic_list[r][1])
    
    Point_haz=[]
    for i in range(len(lon_vect)):
        Point_haz.append([])
        for j in range(len(lat_vect)):
            p_lon= take_closest(Long,lon_vect[i])
            
            AssetRow=seismic_data.loc[seismic_data['Longitude'] == p_lon]
            p_lat= take_closest(AssetRow['Latitude'].values.tolist(),lat_vect[j])
            k=AssetRow.loc[AssetRow['Latitude'] == p_lat]
            
            IM=[list(k['2%'])[0],list(k['10%'])[0],list(k['50%'])[0],list(k['68%'])[0]]
            x = np.array(IM)
            y = np.array(MAF)
            
            def exponential(x,a,k,b):
                return a*np.exp(-x*k)
            
            popt, pcov = scipy.optimize.curve_fit(exponential, x, y, p0=[0.5,5, 1])
            
            if IMc[j][i]==0:
                haz=0
            else:
                haz=popt[0]*np.exp(-IMc[j][i]*popt[1])
            
            if haz<1*10**(-5):
                haz=0
            
            Point_haz[i].append(haz)
    
        if i==0:
            df_Haz_eq=pd.DataFrame(data=Point_haz[i])
        else:
            df_Haz_eq.insert(int(i),str(i),Point_haz[i])
    
    yy=[]
    for i in range(len(lon_vect)):
        for j in range(len(lat_vect)):
            jj=len(lat_vect)-j-1
            yy.append(lat_vect[jj])    
    xx=[]
    for i in range(len(lon_vect)):
        for j in range(len(lat_vect)):
            xx.append(lon_vect[i])   
    
    z=df_Haz_eq.values.tolist()
    zz=[]
    for i in range(len(lon_vect)):
        for j in range(len(lat_vect)):
            zz.append(z[j][i])
    
    xxx=[]
    yyy=[]
    zzz=[]
    for i in range(len(zz)):
        if zz[i]==0:
            0
        else:
            xxx.append(xx[i])
            yyy.append(yy[i])
            zzz.append(math.log10(zz[i])+5)
            
    return df_Haz_eq


'''************************************************************'''
'''                 FRAGILITY CURVE OF SLOPE rain1                   '''
'''************************************************************'''


def Fragility_PvsRain1(Point_slope,Point_soil,d_soilDB):
    
    IMc=[]
    FS=[]
    for i in range(len(Point_soil)):
        IMc.append([])
        FS.append([])
        for j in range(len(Point_soil[0])):
            if Point_soil[i][j]==str(''):
                k=0
            else:
                k=Point_soil[i][j]
            Soil=d_soilDB.iloc[int(k)]            
            k=int(k)
            if k==12 and Point_slope[j][i]>68:
                FST=0
                h_rain=0
            else:
                if k==13 and Point_slope[j][i]>49:
                    FST=0
                    h_rain=0 
                else:
                    if Soil[0]==0:
                        FST=0
                        h_rain=0
                    else:
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
                        H=Soil[10]
                        tb15=Soil[10]
                        tb30=Soil[11]
                        tb45=Soil[12]
                        tb60=Soil[13]
                        lb15=Soil[14]
                        lb30=Soil[15]
                        lb45=Soil[16]
                        lb60=Soil[17]
                        #Slope paramters
                        if Point_slope[j][i]<15:
                            FST=0
                            h_rain=0         
                        else:
                            dip=Point_slope[j][i]        #degree
                            if dip==90:
                                dip=89
                            #Model dimension
                            if dip<15:
                                H=tb15
                                L=lb15/2
                            else:
                                if dip<30:
                                    H=tb15-(dip-15)/15*(tb15-tb30)
                                    L=lb15-(dip-15)/15*(lb15-lb30)
                                    if k>8:
                                        H=tb30
                                        L=lb30/2
                                else:
                                    if dip<45:
                                        H=tb30-(dip-30)/15*(tb30-tb45)
                                        L=(lb30-(dip-30)/15*(lb30-lb45))/2
                                    else:
                                        if dip<60:
                                            H=tb45-(dip-45)/15*(tb45-tb60)
                                            L=(lb45-(dip-45)/15*(lb45-lb60))/2
                                        else:
        #                                    if (dip==60 and dip>60):
                                            H=tb60
                                            L=lb60/2
                                                        
                            Ip=(H**2+L**2)**(1/2)
                            
                            depthW=H       #m
                            if depthW>H:
                                hw=0
                            else:
                                hw=(H-depthW)/H
                    
                            alfa_dw=math.degrees(math.atan(H/L))
                            psi_dw=dip-alfa_dw
                            
                            psi_up=dip+alfa_dw
                            if psi_up>90:
                                psi_up=90
                                
                            #METODO2
                            C=c*2*Ip
                            Wdw=H*L/2*gam
                            Wup=H*L/2*gam
                            Ndw=Wdw*math.cos(math.radians(psi_dw))
                            Nup=Wup*math.cos(math.radians(psi_up))
                            FST=(Ndw*math.tan(math.radians(phi))+Nup*math.tan(math.radians(phi))+C)/((Wdw*math.sin(math.radians(psi_dw)))+(Wup*math.sin(math.radians(psi_up))))  
                            
                            #minutes=[15,30,60,120]
                            minutes=720            #Rain
#                            N=(math.cos(math.radians(dip))**2)*H*gam_w*((hw*(n-1))+(gam/gam_w*(1-n))+(n*Sr0*(1-hw)))
#                            W=((math.cos(math.radians(dip)))*H*gam_w*((hw*(n-1))+(gam/gam_w*(1-n))+(n*Sr0*(1-hw))))
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
                                    h_rain=(x-hw)*(n*H*(1-Sr0))/(math.exp(-kt*(((math.sin(math.radians(dip)))/(n*L*(1-Sr0)))*(minutes*60))))
                                    h_rain=h_rain/(minutes/60)*1000
                                    h_rain=h_rain[0]
                                         
#                                    Infiltration=fc+(f0-fc)*math.exp(-kf*minutes/60)  
                                    Infiltration=fc*minutes/60+(f0-fc)/kf*(1-math.exp(-kf*minutes/60))
                                    if h_rain*minutes/60<Infiltration:
                                        h_rain=h_rain
                                    else:
                                        h_rain=0
            FS[i].append(FST)                
            IMc[i].append(h_rain)         
#    df_FS=pd.DataFrame(FS)
#    df_IMc=pd.DataFrame(IMc)
    return IMc


'''************************************************************'''
'''                    SLOPE HAZARD CURVE                      '''
'''************************************************************'''

def HazardSlopeRain1(lon_vect,lat_vect,Point_slope,Point_soil,d_soilDB,df_RainCoeff,df_Latitude,df_Longitude):
    IMc = Fragility_PvsRain1(Point_slope,Point_soil,d_soilDB)
    
    '''List of coefficents for rain'''
    RainCoeff=df_RainCoeff.values.tolist()
    RainLat=list(df_Latitude[0])
    RainLon=list(df_Longitude[0])
    
    minutes=720
    Tr=[2,5,10,20,50]
    MAF=[]
    IM=[]
    for ii in range(len(Tr)):
        IM.append(34.5*((Tr[ii])**0.19)/(minutes/60)**0.56)
        
    for i in range(len(Tr)):
        MAF.append(1/Tr[i])
    
    Point_haz=[]
    for i in range(len(lon_vect)):
        Point_haz.append([])
        for j in range(len(lat_vect)):
            
            closest_x = take_closest(RainLon, lon_vect[i])
            closest_y = take_closest(RainLat, lat_vect[j])
            ix=RainLon.index(closest_x)
            iy=RainLat.index(closest_y)
            CoefRain=RainCoeff[iy][ix]
        
            x = np.array(IM)*CoefRain
            y = np.array(MAF)
            
            def exponential(x,a,k,b):
                return a*np.exp(-x*k)
            
            popt, pcov = scipy.optimize.curve_fit(exponential, x, y, p0=[5,0.05, 1])
    
            if IMc[j][i]==0:
                haz=0
            else:
                haz=popt[0]*np.exp(-IMc[j][i]*popt[1])
            if haz<1*10**(-5):
                haz=0
            Point_haz[i].append(haz)
        if i==0:
            df_Haz_r1=pd.DataFrame(data=Point_haz[i])
        else:
            df_Haz_r1.insert(int(i),str(i),Point_haz[i])
    
    yy=[]
    for i in range(len(lon_vect)):
        for j in range(len(lat_vect)):
            jj=len(lat_vect)-j-1
            yy.append(lat_vect[jj])    
    xx=[]
    for i in range(len(lon_vect)):
        for j in range(len(lat_vect)):
            xx.append(lon_vect[i])   
    z=df_Haz_r1.values.tolist()
    zz=[]
    for i in range(len(lon_vect)):
        for j in range(len(lat_vect)):
            zz.append(z[j][i]) 
    xxx=[]
    yyy=[]
    zzz=[]
    for i in range(len(zz)):
        if zz[i]==0:
            0
        else:
            xxx.append(xx[i])
            yyy.append(yy[i])
            zzz.append(math.log10(zz[i])+5)
    
    return df_Haz_r1


'''************************************************************'''
'''                 FRAGILITY CURVE OF SLOPE rain2                   '''
'''************************************************************'''


def Fragility_PvsRain2(Point_slope,Point_soil,d_soilDB):
    
    IMc=[]
    FS=[]
    for i in range(len(Point_soil)):
        IMc.append([])
        FS.append([])
        for j in range(len(Point_soil[0])):
            if Point_soil[i][j]==str(''):
                k=0
            else:
                k=Point_soil[i][j]
            Soil=d_soilDB.iloc[int(k)]            
            k=int(k)
            if k==12 and Point_slope[j][i]>68:
                FST=0
                h_rain=0
            else:
                if k==13 and Point_slope[j][i]>49:
                    FST=0
                    h_rain=0
                else:
                    if Soil[0]==0:
                        FST=0
                        h_rain=0
                    else:
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
                        H=Soil[10]
                        tb15=Soil[10]
                        tb30=Soil[11]
                        tb45=Soil[12]
                        tb60=Soil[13]
                        lb15=Soil[14]
                        lb30=Soil[15]
                        lb45=Soil[16]
                        lb60=Soil[17]
                        #Slope paramters
                        if Point_slope[j][i]<15:
                            FST=0
                            h_rain=0         
                        else:
                            dip=Point_slope[j][i]        #degree
                            if dip==90:
                                dip=89
                            #Model dimension
                            if dip<15:
                                H=tb15
                                L=lb15/2
                            else:
                                if dip<30:
                                    H=tb15-(dip-15)/15*(tb15-tb30)
                                    L=lb15-(dip-15)/15*(lb15-lb30)
                                    if k>8:
                                        H=tb30
                                        L=lb30/2
                                else:
                                    if dip<45:
                                        H=tb30-(dip-30)/15*(tb30-tb45)
                                        L=(lb30-(dip-30)/15*(lb30-lb45))/2
                                    else:
                                        if dip<60:
                                            H=tb45-(dip-45)/15*(tb45-tb60)
                                            L=(lb45-(dip-45)/15*(lb45-lb60))/2
                                        else:
        #                                    if (dip==60 and dip>60):
                                            H=tb60
                                            L=lb60/2
                                                        
                            Ip=(H**2+L**2)**(1/2)
                            
                            depthW=H       #m
                            if depthW>H:
                                hw=0
                            else:
                                hw=(H-depthW)/H
                    
                            alfa_dw=math.degrees(math.atan(H/L))
                            psi_dw=dip-alfa_dw
                            
                            psi_up=dip+alfa_dw
                            if psi_up>90:
                                psi_up=90
                                
                            #METODO2
                            C=c*2*Ip
                            Wdw=H*L/2*gam
                            Wup=H*L/2*gam
                            Ndw=Wdw*math.cos(math.radians(psi_dw))
                            Nup=Wup*math.cos(math.radians(psi_up))
                            FST=(Ndw*math.tan(math.radians(phi))+Nup*math.tan(math.radians(phi))+C)/((Wdw*math.sin(math.radians(psi_dw)))+(Wup*math.sin(math.radians(psi_up))))  
                            
                            #minutes=[15,30,60,120]
                            minutes=1440           #Rain
#                            N=(math.cos(math.radians(dip))**2)*H*gam_w*((hw*(n-1))+(gam/gam_w*(1-n))+(n*Sr0*(1-hw)))
#                            W=((math.cos(math.radians(dip)))*H*gam_w*((hw*(n-1))+(gam/gam_w*(1-n))+(n*Sr0*(1-hw))))
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
                                    h_rain=(x-hw)*(n*H*(1-Sr0))/(math.exp(-kt*(((math.sin(math.radians(dip)))/(n*L*(1-Sr0)))*(minutes*60))))
                                    h_rain=h_rain/(minutes/60)*1000
                                    h_rain=h_rain[0]
                                         
#                                    Infiltration=fc+(f0-fc)*math.exp(-kf*minutes/60)  
                                    Infiltration=fc*minutes/60+(f0-fc)/kf*(1-math.exp(-kf*minutes/60))
                                    if h_rain*minutes/60<Infiltration:
                                        h_rain=h_rain
                                    else:
                                        h_rain=0
            FS[i].append(FST)                
            IMc[i].append(h_rain)         
#    df_FS=pd.DataFrame(FS)
#    df_IMc=pd.DataFrame(IMc)
    return IMc


'''************************************************************'''
'''                    SLOPE HAZARD CURVE                      '''
'''************************************************************'''

def HazardSlopeRain2(lon_vect,lat_vect,Point_slope,Point_soil,d_soilDB,df_RainCoeff,df_Latitude,df_Longitude):
    IMc = Fragility_PvsRain2(Point_slope,Point_soil,d_soilDB)

    '''List of coefficents for rain'''
    RainCoeff=df_RainCoeff.values.tolist()
    RainLat=list(df_Latitude[0])
    RainLon=list(df_Longitude[0])
    
    minutes=1440
    Tr=[2,5,10,20,50]
    MAF=[]
    IM=[]
    for ii in range(len(Tr)):
        IM.append(34.5*((Tr[ii])**0.19)/(minutes/60)**0.56)
        
    for i in range(len(Tr)):
        MAF.append(1/Tr[i])
    
    Point_haz=[]
    for i in range(len(lon_vect)):
        Point_haz.append([])
        for j in range(len(lat_vect)):

            closest_x = take_closest(RainLon, lon_vect[i])
            closest_y = take_closest(RainLat, lat_vect[j])
            ix=RainLon.index(closest_x)
            iy=RainLat.index(closest_y)
            CoefRain=RainCoeff[iy][ix]            
            
            x = np.array(IM)*CoefRain
            y = np.array(MAF)
            
            def exponential(x,a,k,b):
                return a*np.exp(-x*k)
            
            popt, pcov = scipy.optimize.curve_fit(exponential, x, y, p0=[5,0.05, 1])

            if IMc[j][i]==0:
                haz=0
            else:
                haz=popt[0]*np.exp(-IMc[j][i]*popt[1])
            if haz<1*10**(-5):
                haz=0     
            Point_haz[i].append(haz)
        if i==0:
            df_Haz_r2=pd.DataFrame(data=Point_haz[i])
        else:
            df_Haz_r2.insert(int(i),str(i),Point_haz[i])
    
    return df_Haz_r2


'''************************************************************'''
'''                        MAP OF TURKEY                       '''
'''************************************************************'''



'''************************************************************'''
'''                   TOTAL HAZARD                      '''
'''************************************************************'''

def MaxHaz(directory,lon_vect,lat_vect,start_time,Point_slope,Point_soil,d_soilDB,df_RainCoeff,df_Latitude,df_Longitude):
    
    df_Haz_eq = HazardSlopeEQ(directory,lon_vect,lat_vect,Point_slope,Point_soil,d_soilDB)
    df_Haz_r1 = HazardSlopeRain1(lon_vect,lat_vect,Point_slope,Point_soil,d_soilDB,df_RainCoeff,df_Latitude,df_Longitude)
    df_Haz_r2 = HazardSlopeRain2(lon_vect,lat_vect,Point_slope,Point_soil,d_soilDB,df_RainCoeff,df_Latitude,df_Longitude)
    
    eq=df_Haz_eq.values.tolist()
    r1=df_Haz_r1.values.tolist()
    r2=df_Haz_r2.values.tolist()
    
    hazard=[]
    for i in range(len(eq)):
        hazard.append([])
        for j in range(len(eq[0])):
            hazard[i].append(max(eq[i][j],r1[i][j],r2[i][j]))
    
    df_hazard=pd.DataFrame(hazard)
    
    yy=[]
    for i in range(len(lon_vect)):
        for j in range(len(lat_vect)):
            jj=len(lat_vect)-j-1
            yy.append(lat_vect[jj])    
    xx=[]
    for i in range(len(lon_vect)):
        for j in range(len(lat_vect)):
            xx.append(lon_vect[i])   
    z=df_hazard.values.tolist()
    zz=[]
    exc=[]
    for i in range(len(lon_vect)):
        exc.append([])
        for j in range(len(lat_vect)):
            if z[j][i]==0:
                zz.append(None)
            else:
                zz.append(z[j][i])
       
    xxx=[]
    yyy=[]
    zzz=[]
    for i in range(len(zz)):
        if zz[i]==0:
            zzz.append(None)
        else:
            xxx.append(xx[i])
            yyy.append(yy[i])
            if zz[i]==None:
                zzz.append(None) 
            else:
                zzz.append(math.log10(zz[i])+5) 
    
    Turkey=gpd.read_file(directory+'gadm36_TUR_0.shx')
    ax = Turkey.plot(color='white', edgecolor='black',figsize=(24,16),linewidth=0.5)
 
    df_latlon=pd.DataFrame(data={'Longitude':xx,'Latitude':yy})

    def Plot_Turkey_Sismicity(Long,Lat,Intensity):      
        points = ax.scatter(Long,Lat, c=Intensity, cmap=plt.cm.Reds, linewidth=0)
        # All the Long points and then Lat points
        plt.colorbar(points, label='Landslide_Hazard')
        plt.title("Landslide Hazard Map (total) - Turkey")
        plt.tight_layout()  # Makes the figure bigger 
        plt.savefig(directory+"Landslide Map of Turkey\\Landslide Hazard Map (total) - "+str(round(700/len(lat_vect)*1650/len(lon_vect)))+" km2 - Turkey.png")
        plt.show()
        
    Plot_Turkey_Sismicity(xxx,yyy,zzz)
    
        #size>4 per poter salvare l'excel    
#    df_hazard.to_excel(directory+"Landslide_Hazard_50years.xlsx")
#    df_latlon.to_excel(directory+"Landslide_longitude&latitude.xlsx")

    print(pd.DataFrame(data={'time [sec]': [(time.time() - start_time)], 
                   '#elements [-]':[len(lat_vect)*len(lon_vect)], 
                   'element size [km2]':[700/len(lat_vect)*1650/len(lon_vect)]}))  

    return df_hazard

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def hazard_landslide_wholeTURKEY(Discretization):
    # discretization = numero adimensionale che scala la grandezza dei vettori latitude and longitude (maggiore il numero e più è grezza la mappa)
    start_time = time.time()
    directory=os.getcwd()+"\\"
    
    image, lat_image, lon_image = load_raster(directory+'SRTM_TR_30_shaped.tif') #ElevationTR_all_small_shaped.tif , SlopeTR1.tif
    #Cover, lat_cover, lon_cover = load_raster(directory+'CLC_LANDCOVER_TR.tif')
    lat_vect=np.linspace(lat_image[10],lat_image[int(len(lat_image)-10)],int(len(lat_image)/Discretization-1))
    lon_vect=np.linspace(lon_image[10],lon_image[int(len(lon_image)-10)],int(len(lon_image)/Discretization-1))  
    
    df_Slope,Point_slope=SlopeDataFrame(lat_vect,lon_vect,image,lon_image,lat_image,directory)
    # Read HWSD data.
    df_HWSD=pd.read_excel(directory+'HWSD_DATA.xlsx')
    
    df_RainCoeff=pd.read_excel(directory+'wc2.1_5m_prec\\Rainfall_Data.xlsx')
    df_Latitude=pd.read_excel(directory+'wc2.1_5m_prec\\Latitude_Rainfall_Data.xlsx')
    df_Longitude=pd.read_excel(directory+'wc2.1_5m_prec\\Longitude_Rainfall_Data.xlsx')
    
    imageHWSD, lat_imageHWSD, lon_imageHWSD = load_raster(directory+'HWSD_TR_shaped.tif')
    df_Soil,Point_soil=SoilDataFrame(lat_vect,lon_vect,lon_imageHWSD,lat_imageHWSD,imageHWSD,df_HWSD)
    d_soilDB=pd.read_excel(directory+'Soil_Rock_Charateristics_DB.xlsx')
    
    df_hazard=MaxHaz(directory,lon_vect,lat_vect,start_time,Point_slope,Point_soil,d_soilDB,df_RainCoeff,df_Latitude,df_Longitude)
   
    return directory,lat_vect,lon_vect,df_hazard

if __name__ == '__main__':
	hazard_landslide_wholeTURKEY(50)
    

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''    