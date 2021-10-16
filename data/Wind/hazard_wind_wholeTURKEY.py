# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 16:50:13 2020

@author: RVA04
"""

import time
import rasterio
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import os

# LOAD RASTER IMAGE
def load_raster(path):
    with rasterio.Env():
        with rasterio.open(path) as src:  
            image = src.read(1)  ## first and only band  
    #src.crs
    #src.transform
    #src.bounds
    lat = [src.xy(i,0)[1] for i in range(src.height)]
    lon = [src.xy(0,i)[0] for i in range(src.width)]    
    return image, lat, lon

# COMPUTING PEAK VELOCITY
def peakvel(directory):
    image, lat, lon = load_raster(directory+'TUR_wind-speed_10m.tif') #https://globalwindatlas.info/
    
    size=10
    wind=[]
    x=[]
    y=[]
    exc=[]
    for i in range(0,len(lat),size):
        exc.append([])
        for j in range(0,len(lon),size):
            if image[i][j]>0 and image[i][j]<10**(3):
                w=image[i][j]/0.9*1.5
                wind.append(w)
                #diviso per 0.9 perchè i dati sono basati su Tr=10 e per passare ai valori di Tr=50 vanno aumentati da 0.9 a 1
                #moltiplicato per 1.5 per ottenere il tempo di picco
                x.append(lon[j])
                y.append(lat[i])
            else:
                w=None
            exc[int(i/size)].append(w)
            
    df_excel=pd.DataFrame(data=exc)
    df_latlon=pd.DataFrame(data={'Longitude':x,'Latitude':y})
    MaxWind=wind.index(max(wind))
    LatMax=y[MaxWind]
    LonMax=x[MaxWind]
    WindMax=wind[MaxWind]
    return x,y,wind,df_excel,WindMax,LonMax,LatMax,df_latlon

# PLOTTING          
def Plot_Turkey_wind(directory,start_time):
    x,y,wind,df_excel,WindMax,LonMax,LatMax,df_latlon = peakvel(directory)   
    
    Turkey=gpd.read_file(directory+'gadm36_TUR_0.shx')
    ax = Turkey.plot(color='white', edgecolor='black',figsize=(24,16),linewidth=0.5)
    
    points = ax.scatter(x,y,c = wind, cmap=plt.cm.Reds, linewidth=0)
    # All the Long points and then Lat points
    plt.colorbar(points, label='Wind_Map')
    plt.title("Wind Map - Turkey")
    plt.tight_layout()  # Makes the figure bigger 
    plt.savefig(directory+"Wind Map - Turkey.png")
    plt.show()
    print(pd.DataFrame(data={'time [sec]': [(time.time() - start_time)], 
                       '#elements [-]':[len(wind)], 
                       'element size [km2]':[700*1650/len(wind)],
                       'MaxSpeed [m/s]':[WindMax]})) 
    
    df_WindHazard=pd.DataFrame({'Longitude':x,'Latitude':y,'Wind':wind})
    return df_WindHazard
    
#size>4 per poter salvare l'excel    
#df_excel.to_excel(directory+"Wind_Peak_Velocity_50years.xlsx")
#df_latlon.to_excel(directory+"Wind_Peak_longitude&latitude.xlsx") 

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

def hazard_wind_wholeTURKEY():
    start_time = time.time()
    directory=os.getcwd()+"\\"
    df_WindHazard = Plot_Turkey_wind(directory,start_time)
    return directory, df_WindHazard

if __name__ == '__main__':
	hazard_wind_wholeTURKEY()    

''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''


