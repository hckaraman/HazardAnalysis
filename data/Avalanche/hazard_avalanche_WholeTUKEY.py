# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 13:55:12 2020

@author: RVA04
"""

import time
import geopandas as gpd
import pandas as pd
import rasterio
import numpy as np
#from rasterio.features import shapes
#import descartes
from bisect import bisect_left
import math
import scipy.optimize 
import matplotlib.pyplot as plt
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
       
    AltNS=[]
    for i in range(3):
        if -iy+i>int(len(y)):
            AltNS.append(raster[-iy][ix])
        else:
            AltNS.append(raster[-iy+i][ix])
    
        
    SlopeNS=[]
    for i in range(2):
        if -iy+i>int(len(y)):
            SlopeNS.append(0)
        else:
            SlopeNS.append(((-raster[-iy+i][ix]+raster[-iy+i+1][ix])/dy)/math.pi*180)
    
    if -iy+i>int(len(y)):
        SlopeVarNS=0
    else:  
        SlopeVarNS=(SlopeNS[0]-SlopeNS[1])
    
    SlopeVarSN=SlopeVarNS
    SlopeSN=max(SlopeNS)
    AltitudeSN=(AltNS[0]+AltNS[1]+AltNS[2])/3
    
    return SlopeSN,SlopeVarSN,AltitudeSN
                
def getSlopeOrientation(SNsec, EWsec): 
    pass

#TEST - check altitude at a known location    
def SlopeDataFrame(lat_vect,lon_vect,image,lon_image,lat_image):
   
    Point_slope=[]
    Point_slopevar=[]
    Point_alt=[]
    for i in range(len(lon_vect)):
        Point_slope.append([])
        Point_slopevar.append([])
        Point_alt.append([])
        for j in range(len(lat_vect)):
            jj=len(lat_vect)-j-1
            p_lon=lon_vect[i]
            p_lat=lat_vect[jj]
       
            SlopeSN,SlopeVarSN,AltitudeSN = getContextMorphology(image, lon_image,lat_image,p_lon,p_lat, n=1, l_min = 250)

            Point_slope[i].append(SlopeSN)
            Point_slopevar[i].append(SlopeVarSN)
            Point_alt[i].append(AltitudeSN)
            #se negativo guarda a sud
            
        if i==0:
            df_Slope=pd.DataFrame(data=Point_slope[i])
            df_SlopeVar=pd.DataFrame(data=Point_slopevar[i])
            df_Alt=pd.DataFrame(data=Point_alt[i])
        else:
            df_Slope.insert(int(i),str(i),Point_slope[i])
            df_SlopeVar.insert(int(i),str(i),Point_slopevar[i])
            df_Alt.insert(int(i),str(i),Point_alt[i])
            
    return df_Slope,df_SlopeVar,df_Alt

'''************************************************************'''
'''                 FRAGILITY CURVE OF SLOPE                   '''
'''************************************************************'''



def Fragility_Avalanche(directory,lat_vect,lon_vect,image,lon_image,lat_image):
    df_Slope,df_SlopeVar,df_Alt=SlopeDataFrame(lat_vect,lon_vect,image,lon_image,lat_image)
#    d_soilDB=pd.read_excel(directory+'Soil_Rock_Charateristics_DB.xlsx')
    Point_slope=df_Slope.values.tolist()
    Point_slopeVar=df_SlopeVar.values.tolist()
    Point_alt=df_Alt.values.tolist()
    
#    IMc=[]
    FS=[]
    for i in range(len(Point_slope)):
#        IMc.append([])
        FS.append([])
        for j in range(len(Point_slope[0])):
            #Soil paramters
            gam=1.96    #kN/m3
            phi=22.5    #degree
            c=0.9       #kPa
            tb30=3.8
            tb45=2.2
            tb60=1.9
            lb30=46
            lb45=10
            lb60=5
            #Slope paramters
            if Point_alt[i][j]<2000: #2700 montagne del caucaso
                FST=0
            else:
                if Point_slope[i][j]<30: 
                    FST=0
                else:
                    if Point_slope[i][j]>50:
                        FST=0
                    else:
                        if Point_slopeVar[i][j]<10:
                            FST=0
#                IMT=0         
                        else:
                            dip=Point_slope[i][j]        #degree
                            #Model dimension 
                            if dip<45:
                                H=tb30-(dip-30)/15*(tb30-tb45)                  #depth of failure surface
                                L=(lb30-(dip-30)/15*(lb30-lb45))/2              #L is half of the lenght of the failure on the terrain level (it is l/2 in excel)
                            else:
                                H=tb45-(dip-45)/15*(tb45-tb60)
                                L=(lb45-(dip-45)/15*(lb45-lb60))/2              #metà lunghezza perchè il pendio è diviso in due parti per il calcolo
                                
                            Ip=(H**2+L**2)**(1/2)
                            
                            alfa_dw=math.degrees(math.atan(H/L))
                            psi_dw=dip-alfa_dw
                            if psi_dw<0:
                                psi_dw=0
                            
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
            #               IMT=(FST-1)*math.tan(math.radians(dip))/(math.tan(math.radians(phi))*math.tan(math.radians(dip))+1) #the formula has *g, but the IM is in g and not m/s2      
            FS[i].append(FST)
#            IMc[i].append(IMT)        
    df_FS=pd.DataFrame(FS)
#    df_IMc=pd.DataFrame(IMc)

    return df_FS



'''************************************************************'''
'''************************************************************'''
'''************************************************************'''


def HazardSlope(directory,start_time,lat_vect,lon_vect,image,lon_image,lat_image):
    df_FS=Fragility_Avalanche(directory,lat_vect,lon_vect,image,lon_image,lat_image)
    FS_list=df_FS.values.tolist()
    
    yy=[]
    for i in range(len(lon_vect)):
        for j in range(len(lat_vect)):
            jj=len(lat_vect)-j-1
            yy.append(lat_vect[jj])    
    xx=[]
    for i in range(len(lon_vect)):
        for j in range(len(lat_vect)):
            xx.append(lon_vect[i])
            
    z=FS_list
    zz=[]
    for i in range(len(lon_vect)):
        for j in range(len(lat_vect)):
            if z[j][i]==0:
                zz.append(None)
            else:
                zz.append(1/z[j][i])
                print(lon_vect[i],lat_vect[j])
    
    df_hazard=pd.DataFrame(FS_list)
    df_latlon=pd.DataFrame(data={'Longitude':xx,'Latitude':yy})
    
    Turkey=gpd.read_file(directory+'gadm36_TUR_0.shx')
    ax = Turkey.plot(color='white', edgecolor='black',figsize=(24,16),linewidth=0.5)
    
    def Plot_Turkey_Sismicity(Long,Lat,Intensity):
        points = ax.scatter(Long,Lat, c=Intensity, cmap=plt.cm.Reds, linewidth=0)
        # All the Long points and then Lat points
        plt.colorbar(points, label='Avalanches_Hazard')
        plt.title("Snow FS - Avalanches Hazard Map - Turkey")
        plt.tight_layout()  # Makes the figure bigge
        plt.savefig(directory+"Avalanches Hazard Map - "+str(round(700/len(lat_vect)*1650/len(lon_vect)))+" km2 - Turkey.png")
        plt.show()       
    Plot_Turkey_Sismicity(xx,yy,zz)
    df_AvalancheHazard=pd.DataFrame({'Longitude':xx,'Latitude':yy,'Wind':zz})    
    return df_AvalancheHazard

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def hazard_avalanche_wholeTURKEY(Discretization):
    start_time = time.time()
    directory=os.getcwd()+"\\"
    
    image, lat_image, lon_image = load_raster(directory+'SRTM_TR_30_shaped.tif') #ElevationTR_all_small_shaped.tif , SlopeTR1.tif
    #Cover, lat_cover, lon_cover = load_raster(directory+'CLC_LANDCOVER_TR.tif')
    lat_vect=np.linspace(lat_image[5],lat_image[int(len(lat_image)-5)],int(len(lat_image)/Discretization-1))
    lon_vect=np.linspace(lon_image[5],lon_image[int(len(lon_image)-5)],int(len(lon_image)/Discretization-1))
    
    df_AvalancheHazard = HazardSlope(directory,start_time,lat_vect,lon_vect,image,lon_image,lat_image)
    return directory, df_AvalancheHazard,lat_vect,lon_vect

if __name__ == '__main__':
	hazard_avalanche_wholeTURKEY(100)    

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''