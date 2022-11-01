# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 16:44:16 2022

@author: SKirillov 
"""

def gradapproach(X, Y, xb=0, yb=0, max_approaches=4):
    
    import numpy as np
    
    nc = X.shape[1]-1
    nr = X.shape[0]-1
    
    def inner(x, y):
        nonlocal xb, yb
        status = 0
        count = 0
        
        for run in range(1, max_approaches+1):
            
            Y0 = Y[yb, xb]
            X0 = X[yb, xb]
            Xy = X[yb+1, xb] - X0
            Yy = Y[yb+1, xb] - Y0
            Xx = X[yb, xb+1] - X0
            Yx = Y[yb, xb+1] - Y0
            
            if Yx == 0:
                dy = (y-Y0) / Yy
                dx = (x-X0) / Xx
            else:
                dy = (x - X0 - (y-Y0)*Xx/Yx) / (Xy - Yy/Yx*Xx)
                dx = (y-Y0 - dy*Yy) / Yx
    
            xb_n = np.floor(xb+dx).astype(int)
            yb_n = np.floor(yb+dy).astype(int)
            
            # print(run, dx, dy, xb, yb, xb_n, yb_n)
    
            if run == max_approaches:
                if dx > 1:  xb = xb;      dx = 2-dx
                if dy > 1:  yb = yb;      dy = 2-dy
                if dx < 0:  xb = xb - 1;  dx = dx + 1
                if dy < 0:  yb = yb - 1;  dy = dy + 1
    
                if (yb < 0) | (yb >= nr):
                    status = 0
                    break
                if (xb < 0) | (xb >= nc):
                    status = 0
                    break
                status = 1
                break
    
            if (xb_n == xb) & (yb_n == yb):
                status = 1
                break
    
            if run < max_approaches:
                if xb_n >= nc:  xb_n = nc-1;  count += 1
                if xb_n < 0:    xb_n = 0;     count += 1
                if yb_n >= nr:  yb_n = nr-1;  count += 1
                if yb_n < 0:    yb_n = 0;     count += 1
    
            if count > 2:
                status = 0
                break
    
            xb = xb_n
            yb = yb_n
    
        return status, xb, yb, dx, dy, run
    return inner

def interpolate(xb, yb, dx, dy, values, min_points=1):
    
    import numpy as np
    
    d = np.full([2,2], np.nan)
    d[0, 0] = np.sqrt(dx**2 + dy**2) + 1e-24
    d[1, 0] = np.sqrt(dx**2 + (1-dy)**2) + 1e-24
    d[1, 1] = np.sqrt((1-dx)**2 + (1-dy)**2) + 1e-24
    d[0, 1] = np.sqrt((1-dx)**2 + dy**2) + 1e-24
    
    z = []
    dd = 0
    
    for (j, i) in ([(0, 0), (1, 0), (1, 1), (0, 1)]):
        if np.isfinite(values[yb+j, xb+i]):
            z.append(values[yb+j, xb+i]/d[j, i])
            dd += 1/d[j, i]
    
    if len(z) < min_points:
        return np.nan
    else:
        return sum(z)/dd


def interpolate_lonlat(lat, xb, yb, dx, dy, values, min_points=1):
    
    import numpy as np
    
    scale = np.cos(np.deg2rad(lat))
        
    d = np.full([2,2], np.nan)
    d[0, 0] = np.sqrt((dx*scale)**2 + dy**2) + 1e-6
    d[1, 0] = np.sqrt((dx*scale)**2 + (1-dy)**2) + 1e-6
    d[1, 1] = np.sqrt(((1-dx)*scale)**2 + (1-dy)**2) + 1e-6
    d[0, 1] = np.sqrt(((1-dx)*scale)**2 + dy**2) + 1e-6
    
    z = []
    dd = 0
    
    for (j, i) in ([(0, 0), (1, 0), (1, 1), (0, 1)]):
        if np.isfinite(values[yb+j, xb+i]):
            z.append(values[yb+j, xb+i]/d[j, i])
            dd += 1/d[j, i]
    
    if len(z) < min_points:
        return np.nan
    else:
        return sum(z)/dd




def inside_outside(x, y, x_point, y_point, tolerance=1e-6):
    """
    The function that returns the following statused of 
    (x_point, y_point) point within (x, y) polygon
    
    The function returns:    
        0 - corner
        +1/-1 - inner  -1 = CCW, +1 = CW
        2 - edge
        3 - outter
        4 - unknown (somethin is wrong)
    
    tolerance in degrees
    """
    
    import numpy as np
    
    if (x[0] != x[-1]) | (y[0] != y[-1]):
        x = np.append(x, x[0])
        y = np.append(y, y[0])
        
    dx = x - x_point
    dy = y - y_point
    directions = np.arctan2(dy, dx)*180/np.pi
    angles = np.full([len(x)-1], np.nan)
    
    for a in range(0, len(x)-1):
        if (dx[a]==0) & (dy[a]==0):
            status = 0
            return status
        if a < len(x):
            angles[a] = directions[a+1]-directions[a]
            if abs(abs(angles[a]) - 180) <= tolerance:
                status = 2
                return status
            if angles[a] > 180:
                angles[a] -= 360
            if angles[a] < -180:
                angles[a] += 360
    
    d = np.sum(angles)
    
    if abs(d) <= tolerance:
        status = 3
    elif abs(d+360) <= tolerance:
        status = 1
    elif abs(d-360) <= tolerance:
        status = -1
    else:
        status = 4
    
    return status





def read_coastlines(min_points=1000, roughness=1, limits=[-180, 180, -90, 90]):
    
    import pickle
    import numpy as np
    import os

    path = 'd:/work/GeoData/OSM coastlines/'
    lat = []
    lon = []
    
    with open(path+'Coastlines_lon.pkl', 'rb') as file:
        x_res = pickle.load(file)
    with open(path+'Coastlines_lat.pkl', 'rb') as file:
        y_res = pickle.load(file)      
    volume = os.path.getsize(path + 'Coastlines_lon.pkl') + os.path.getsize(path + 'Coastlines_lat.pkl')
      
    for x, y in zip(x_res, y_res):
        if (len(x) >= min_points):
            positions = np.where((x >= limits[0]) & (x <= limits[1]) & \
                                 (y >= limits[2]) & (y <= limits[3]))[0]
            if len(positions) > 1:
                lon.append(x[positions][::roughness])
                lat.append(y[positions][::roughness])

    n_in = np.array([len(i) for i in y_res])
    n_out = np.array([len(i) for i in lon])
    compression = np.sum(n_out)/np.sum(n_in)
    
    print('read_coastlines():')
    print('          With given parameters and within given limits')
    print(f'          {compression*100:.2f}% ({compression*volume*1e-6:.1f}Mb of {volume*1e-6:.0f}Mb) of the raw data is used')
    print('')
        
    if compression>0.05:
        print('          With the given parameters the coastline subset may be too large for plotting')
        print('          Think of attenuating the parameters to reduce the subset size')
                      
        
    return lon, lat


def read_coastlines2(min_length=50, averaging=2, limits=[-180, 180, -90, 90]):
    """
    averaging - is a distance in km
    """

    import pickle
    import numpy as np
    import os
    import geopy.distance as dist
    from utils import printProgressBar
    from pathlib import Path
    
    def smooth_coast(lon_subset, lat_subset, step, averaging):
        
        lat = []
        lon = []
        
        for count, (x, y) in enumerate(zip(lon_subset, lat_subset)):
            positions = np.where((x >= limits[0]) & (x <= limits[1]) & \
                                 (y >= limits[2]) & (y <= limits[3]))[0]
            
            printProgressBar(count+1, len(lon_subset), 
                             prefix = f'coastline element {count+1} out of {len(lon_subset)}:', suffix = 'Complete', length = 50) 
            
            x = x[positions]
            y = y[positions]
                
            start = 0    
            end = 1
            lon_av = []
            lat_av = []
                
            while end < len(x) - (step+1):
                   
                while (dist.distance((y[start], x[start]),(y[end], x[end])).km <= averaging) & (end < (len(x) - (step+1))):
                    end += step
                    
                lon_av.append(np.mean(x[start : end]))        
                lat_av.append(np.mean(y[start : end]))   
                start = end
                                    
            if len(lon_av) >= 3:

                lon_av = np.array(lon_av)
                lat_av = np.array(lat_av)
                
                diff = (lon_av[2:len(lon_av)] + lon_av[0:len(lon_av)-2])/2 - lon_av[1:len(lon_av)-1]
                index = np.where(abs(diff) > 1)[0]
                lon_av = np.delete(lon_av, 2+index)
                lat_av = np.delete(lat_av, 2+index)
                
                lon_av = np.append(lon_av, lon_av[0])
                lat_av = np.append(lat_av, lat_av[0])
                
                lon.append(lon_av)
                lat.append(lat_av)
                
            else:
                
                lon_av = []
                lat_av = []
                
                for i in (0,1,2):
                    delta = np.floor(len(x)/3).astype(int)
                    lon_av.append(np.mean(x[delta*i : delta*(i+1)]))
                    lat_av.append(np.mean(y[delta*i : delta*(i+1)]))
                                
                lon_av.append(lon_av[0])
                lat_av.append(lat_av[0])

                lon.append(np.array(lon_av))
                lat.append(np.array(lat_av))
                
        print('')    
            
        return lon, lat


    path = './coastlines/'

    if Path(path + 'min_length=' + str(min_length) + 'km_averaging=' + str(averaging)+ 'km.pkl').exists():
        with open(path + 'min_length=' + str(min_length) + 'km_averaging=' + str(averaging)+ 'km.pkl', 'rb') as file:
            lon, lat = pickle.load(file)
    else:

        with open(path+'Coastlines.pkl', 'rb') as file:
            x_res, y_res, length_res = pickle.load(file)
        
        lat_subset = [y for (i, y) in enumerate(y_res) if length_res[i] >= min_length]
        lon_subset = [x for (i, x) in enumerate(x_res) if length_res[i] >= min_length]
          
        lon, lat = smooth_coast(lon_subset, lat_subset, np.round(averaging/0.05/10).astype(int), averaging)
                          
        with open(path + 'min_length=' + str(min_length) + 'km_averaging=' + str(averaging)+ 'km.pkl', 'wb') as file:
            pickle.dump([lon, lat], file)   
        
    return lon, lat



def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd, flush=True)
    # Print New Line on Complete
    if iteration == total: 
        print()
      
