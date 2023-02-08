# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 13:00:19 2022

@author: SKirillov
"""

import os
import pickle
import jigsawpy
import yaml
import multiprocessing

import cartopy.crs as ccrs
import cartopy.feature as cf
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.cm as cmp
import geopy.distance as dist
import netCDF4 as nc
import numpy as np
import geopy.distance as dist

from fastkml import geometry, kml
from pathlib import Path
from utils import (gradapproach, inside_outside, interpolate_lonlat,
                   printProgressBar, read_coastlines2)


def plotting(lon, lat, triangles, depths, coastnode, settings):
    
    domain=[  
        settings["plot_region"]["min_lon"], settings["plot_region"]["max_lon"], 
        settings["plot_region"]["min_lat"], settings["plot_region"]["max_lat"] ]
    
    proj_work = ccrs.__dict__[ settings["plot_region"]["projection"] ] \
                 (central_longitude= settings["plot_region"]["central_longitude"])
    
    
    active_nodes = np.where((domain[0] < lon) & (lon < domain[1]) & \
                            (domain[2] < lat) & (lat < domain[3]) & \
                            (depths > 0)     )[0]
    
    temp_set = set(active_nodes)
       
    active_elements = [n for n in range(0, len(triangles))
                           if ((triangles[n][0] in temp_set)
                             & (triangles[n][1] in temp_set)
                             & (triangles[n][2] in temp_set))
                            ]
    
    # PLOT MESH
    
    fig = plt.figure(figsize=(15, 15))
    ax = plt.axes(projection=proj_work)
    ax.set_extent(domain)
    ax.gridlines(draw_labels=True, xlocs=None, ylocs=None)
    ax.hold_limits()
    
    for triangle in triangles[active_elements]:
        plt.fill( [ lon[triangle[0]], lon[triangle[1]], lon[triangle[2]], lon[triangle[0]] ],
                  [ lat[triangle[0]], lat[triangle[1]], lat[triangle[2]], lat[triangle[0]] ], 
                    linewidth=0.1, edgecolor='g', facecolor='g', alpha=0.2, transform=ccrs.Geodetic())   
    
    lon_p = lon[active_nodes]
    lat_p = lat[active_nodes]
    coastnode_p = coastnode[active_nodes]
    
    plt.plot(lon_p, lat_p, 'o', markerfacecolor='g', markeredgewidth=0, markersize=0.75, transform=ccrs.Geodetic())
    index = np.where(coastnode_p == 1)[0]
    plt.plot(lon_p[index], lat_p[index], 'o', markerfacecolor='b', markeredgewidth=0, markersize=1.00, transform=ccrs.Geodetic())
        
    plt.savefig('_mesh.jpg', dpi=600, format='jpg')
    
    
    # PLOT DEPTHS
    
    fig = plt.figure(figsize=(15, 15))
    ax = plt.axes(projection=proj_work)
    ax.set_extent(domain)
    ax.gridlines(draw_labels=True, xlocs=None, ylocs=None)
    ax.hold_limits()
    
    ind = {n:index for index, n in enumerate(active_nodes)} 
    
    triangles_subset = triangles[active_elements]
    for index, j in enumerate(triangles_subset):
        for i in (0,1,2):
            triangles_subset[index, i] = ind[j[i]]    
    
    triangles_subset = [list(tria) for tria in triangles_subset]
    
    triang = mtri.Triangulation(lon[active_nodes], lat[active_nodes], triangles_subset)
    
    contourf = ax.tricontourf(triang, depths[active_nodes], transform=ccrs.Geodetic(), 
                                cmap = cmp.jet)
    
    cax = fig.add_axes([ax.get_position().x1+0.01, ax.get_position().y0, 0.02,ax.get_position().height])
    plt.colorbar(contourf, cax=cax, orientation='vertical', pad=0.03,)
    
    plt.savefig('_depths.jpg', dpi=600, format='jpg')
    
    
def find_closest_point_of_polygon(x, y, x_poly, y_poly):    
    
    distance_to_poly = 1e10
        
    for index in range(0, len(x_poly) - 1):
            
        distance_1 = np.sqrt((x_poly[index] - x) ** 2 + (y_poly[index] - y) ** 2)
        if distance_1 < distance_to_poly:
           (nearest_xp, nearest_yp) = (x_poly[index], y_poly[index])
           distance_to_poly = distance_1
            
        
        if x_poly[index+1] - x_poly[index] != 0:
            a = (y_poly[index+1] - y_poly[index])/(x_poly[index+1] - x_poly[index])
            b = y_poly[index] - a * x_poly[index]
        else:
            a = np.nan
            
        
        if not np.isnan(a):
            if a != 0:
                xp = (y + x / a - b) / (a + 1/a)
                
            else:
                xp = x
                
            yp = a * xp + b
        
        else:
            xp = x_poly[index]
            yp = y
            
            
        distance_2 = np.sqrt((x - xp) ** 2 + (y - yp) ** 2)
        if (distance_2 < distance_to_poly) and (0 < (xp - x_poly[index])/(x_poly[index+1] - x_poly[index]) < 1):
           (nearest_xp, nearest_yp) = (xp, yp)
           distance_to_poly = distance_2
               
    return nearest_xp, nearest_yp, distance_to_poly
            

def create_tria_in_node_dictionary(triangles):
    
    """
    Creating a dictionary of nodes with indices of linked elements 
    i.e. { 0: [123, 0, 45, 345], ..., 234: [23, 24, 265, 4350, 4352], etc. }
    """
    
    tria_in_node = {}
    for n, triangle in enumerate(triangles):
        if triangle[0] > -1:
            for i in (0, 1, 2):
                if triangle[i] not in tria_in_node:
                    tria_in_node[triangle[i]] = [n]
                else:
                    tria_in_node[triangle[i]].append(n)

    return tria_in_node


# CREATETHE FIND_TOPO CLASS WITH AN ENCLOSED GRADAPPROACH METHOD --------------


class Find_topo:
    def __init__(self, settings):
        topo_path = settings[
            "topo_path"
        ]  # "./topo/RTopo-2.0.1_30sec_bos_fix_lowres_D3.nc"
        self.__ds = nc.Dataset(topo_path)
        self.__lat, self.__lon = (
            self.__ds.variables["lat"][:, 1],
            self.__ds.variables["lon"][:, 0],
        )
        self.topo = self.__ds.variables["topo"][:, :]

        self.topo[21600, :] = self.topo[21599, :]

        self.lon, self.lat = np.meshgrid(self.__lon, self.__lat)
        self.ga = gradapproach(self.lon, self.lat)


# READ KML FILE AND SET RESOLUTION FOR THE INTERNAL AREA ----------------------
#   presision sets the amount of nodes for each segment of kml contour for
#   the followed linear interpolation between internal and external contours


def from_kml(file, order=True):
    """
    Opens kml files that define regions with different resolution from the base one.

    Parameters
    ----------
    file: str
        path to the .kml file
    order: bool
        If True, fist polygone is internal and second is external. I False, then other way round.

    """

    with open(file) as kml_file:
        doc = kml_file.read().encode("utf-8")
        k = kml.KML()
        k.from_string(doc)
        for feature_0 in k.features():
            for step, feature_1 in enumerate(feature_0.features()):
                for step, feature_2 in enumerate(feature_1.features()):
                    lon, lat = np.squeeze(
                        np.vsplit(np.array(feature_2.geometry.coords.xy), 2)
                    )
                    lon = np.append(lon, lon[0])
                    lat = np.append(lat, lat[0])
                    if step == 0:
                        internal = [lon, lat]
                    else:
                        outer = [lon, lat]

        proj_work = ccrs.Stereographic(central_longitude=internal[0][0], 
                                       central_latitude=internal[1][0])
        x, y = proj_work.transform_point(internal[0][0], internal[1][0], ccrs.Geodetic())
        (x_out, y_out, _) = np.hsplit(
            proj_work.transform_points(ccrs.Geodetic(), outer[0], outer[1]), 3
        )
        if inside_outside(x_out, y_out, x, y) == 3:
            
            internal, outer = outer, internal

    return internal, outer


# THE FUNCTION FOR ADJUSTING RESOLUTIONS TO THE BATHYMETRY
#    Need to be called if needed!


def bathymetry_adjustment(settings, latitudes, longitudes, result):

    print("")
    print("adjusting resolutions over bathymetry")
    
    topo = Find_topo(settings)

    for j, lat in enumerate(latitudes):
        printProgressBar(
            j + 1, len(latitudes), prefix="Progress:", suffix="Complete", length=50
        )

        for i, lon in enumerate(longitudes):

            [status, xb, yb, dx, dy, run] = topo.ga(lon, lat)

            depth = interpolate_lonlat(lat, xb, yb, dx, dy, topo.topo, min_points=1)
            if depth > 0:
                pass
            elif depth >= -settings["refine_over_bathymetry"]["min_depth"]:
                
                result[j, i] = settings["refine_over_bathymetry"]["min_resolution"]
                
            elif depth >= -settings["refine_over_bathymetry"]["max_depth"]:
                
                result[j, i] = settings["refine_over_bathymetry"]["min_resolution"] +       \
                    (result[j, i] - settings["refine_over_bathymetry"]["min_resolution"]) * \
                    (abs(depth) - settings["refine_over_bathymetry"]["min_depth"])/         \
                    (settings["refine_over_bathymetry"]["max_depth"] - settings["refine_over_bathymetry"]["min_depth"])     
                
    return result


def estimate_number_of_nodes(result, longitudes, latitudes):
    
    total = 0
    
    dlat = abs(latitudes[1] - latitudes[0])
    dlon = abs(longitudes[1] - longitudes[0])
    
    for j, lat in enumerate(latitudes):
        for i, lon in enumerate(longitudes):
            area = dlat * 111 * dlon * 111 * np.cos(np.deg2rad(lat))
            element_area = 0.7 * (result[j, i] ** 2) / 2
            total += area/element_area
    
    print("")
    print(f'The estimated # of elements to build with jigsaw (including land) is { np.round(total).astype(int)}')
    
    return


def init_multiprocessing(args):
    global counter
    counter = args
    

def inner(parameters):
    
    with counter.get_lock():
        counter.value += 1
    
    result_at_lon = parameters[0]    
    longitudes = parameters[1]
    latitudes = parameters[2]
    i = parameters[3]
    grid_lon = parameters[4]
    min_resolution = parameters[5]
    max_distance = parameters[6]
    lon_coast = parameters[7]
    lat_coast = parameters[8]
    
    output = result_at_lon
    
    for j, grid_lat in enumerate(latitudes):
        
        if (grid_lat >- 85) & (grid_lat < 85):
        
            min_distance = 1e+10
            d1 = max_distance/111
            d2 = max_distance/111/np.cos(np.deg2rad(grid_lat))
            variant = 1
            if grid_lon - d2 < -180: variant = 2
            if grid_lon + d2 > 180: variant = 3
                        
            for polygon in range(0, len(lon_coast)):
                
                if variant == 1:
                    cut = np.where((lat_coast[polygon] >= grid_lat - d1) & \
                                   (lat_coast[polygon] <= grid_lat + d1) & \
                                   (lon_coast[polygon] >= grid_lon - d2) & \
                                   (lon_coast[polygon] <= grid_lon + d2) )[0]
                elif variant == 2:
                    cut = np.where((lat_coast[polygon] >= grid_lat - d1) & \
                                   (lat_coast[polygon] <= grid_lat + d1) & \
                                   ((lon_coast[polygon] >= 360 + grid_lon - d2) | \
                                    (lon_coast[polygon] <= grid_lon + d2)) )[0]
                elif variant == 3:
                    cut = np.where((lat_coast[polygon] >= grid_lat - d1) & \
                                   (lat_coast[polygon] <= grid_lat + d1) & \
                                   ((lon_coast[polygon] >= grid_lon - d2) | \
                                    (lon_coast[polygon] <= grid_lon + d2 - 360)) )[0]
                    
                               
                for y, x in zip(lat_coast[polygon][cut], lon_coast[polygon][cut]):
                    distance = dist.distance((y, x),(grid_lat, grid_lon)).km 
                    if distance < min_distance:
                        min_distance = distance
            
            if min_distance < max_distance:        
                output[j] = min_resolution + (result_at_lon[j] - min_resolution) * \
                               min_distance / max_distance
    
    printProgressBar(counter.value, len(longitudes), 
                     prefix = 'Progress:', suffix = 'Complete', length = 50) 
    
    return(i, output)


# ADD LINEARLY INTERPOLATED RESOLUTIONS (BETWEEN "min_resolution" in km AT THE COAST
#      AND BACKGROUND RESOLUTION AT THE "max_distance" in km FROM THE COAST)
#      "min_length" is the minimum length of coatline in km used for calculations
#      "averaging" is the smoothing parameter in km

def refine_along_coastlines(
    longitudes, latitudes, result, min_resolution, max_distance, min_length, averaging
):
    
    counter = multiprocessing.Value('i', 0)
    proc_num = multiprocessing.cpu_count()
  
    print('')    
    print('preparing smoothed coastlines')
    
    lon_coast, lat_coast = read_coastlines2(min_length=min_length, averaging=averaging)
    
    print('')    
    print(f'refining along coastlines on multiple ({proc_num}) CPUs')

    parameters = tuple([
        (np.squeeze(result[:, i]), longitudes, latitudes, i, grid_lon,
         min_resolution, max_distance, lon_coast, lat_coast)
        for i, grid_lon in enumerate(longitudes)
        ])
    
    pool = multiprocessing.Pool(processes=proc_num, initializer=init_multiprocessing, 
                                initargs = (counter, ))
    temp_results = pool.map(inner, parameters)
    pool.close()

    for temp in temp_results:
        result[:, temp[0]] = temp[1]   

    return result


# DEFINE RESOLUTIONS WITHIN POLYGONS FROM *.KML FILES


def refine(region, longitudes, latitudes, result):

    (poly_in, poly_out, resolution) = (
        region["Polygon inside"],
        region["Polygon outside"],
        region["resolution"]  )

    min_lon = np.min(poly_out[0])
    max_lon = np.max(poly_out[0])
    min_lat = np.min(poly_out[1])
    max_lat = np.max(poly_out[1])

    # checking if the current polygon encircles poles for adjusting min_lat and max_lat
    proj_work = ccrs.Stereographic(central_longitude=0, central_latitude=-90)
    x, y = proj_work.transform_point(0, -90, ccrs.Geodetic())
    (x_out, y_out, _) = np.hsplit(
        proj_work.transform_points(ccrs.Geodetic(), poly_out[0], poly_out[1]), 3
    )
    if abs(inside_outside(x_out, y_out, x, y)) == 1:
        min_lat = -90
        min_lon = -180
        max_lon = 180
    proj_work = ccrs.Stereographic(central_longitude=0, central_latitude=90)
    x, y = proj_work.transform_point(0, 90, ccrs.Geodetic())
    (x_out, y_out, _) = np.hsplit(
        proj_work.transform_points(ccrs.Geodetic(), poly_out[0], poly_out[1]), 3
    )
    if abs(inside_outside(x_out, y_out, x, y)) == 1:
        max_lat = 90
        min_lon = -180
        max_lon = 180

    variant = 1

    if max_lon > 180:
        variant = 2
        index = np.where(poly_out[0] > 180)[0]
        poly_out[0][index] = poly_out[0][index] - 360
        max_lon = np.min(poly_out[0][poly_out[0] > 0])
        min_lon = np.max(poly_out[0][poly_out[0] < 0])

    print("")
    print(f"{region['name']}: resolution {region['resolution']} km")
    print(
        f"{min_lon:.1f}   ", f"{max_lon:.1f}   ", f"{min_lat:.1f}   ", f"{max_lat:.1f}"
    )

    for i, grid_lon in enumerate(longitudes):

        printProgressBar(
            i + 1, len(longitudes), prefix="Progress:", suffix="Complete", length=50
        )

        for j, grid_lat in enumerate(latitudes):

            if ((variant == 1) & (min_lon <= grid_lon <= max_lon) & 
                (min_lat <= grid_lat <= max_lat)) | \
                ((variant == 2) & (max_lon <= grid_lon) & (min_lat <= grid_lat <= max_lat)) | \
                ((variant == 2) & (min_lon >= grid_lon) & (min_lat <= grid_lat <= max_lat)):

                proj_work = ccrs.Stereographic(
                    central_longitude=grid_lon, central_latitude=grid_lat
                )
                (x_out, y_out, _) = np.hsplit(
                    proj_work.transform_points(
                        ccrs.Geodetic(), poly_out[0], poly_out[1]
                    ),
                    3,
                )
                (x_in, y_in, _) = np.hsplit(
                    proj_work.transform_points(ccrs.Geodetic(), poly_in[0], poly_in[1]),
                    3,
                )

                x, y = proj_work.transform_point(grid_lon, grid_lat, ccrs.Geodetic())

                if abs(inside_outside(x_in, y_in, x, y)) < 3:

                    result[j, i] = resolution

                elif abs(inside_outside(x_out, y_out, x, y)) < 3:
                    
                    # DISTANCE TO THE INTERNAL CONTOUR
                    _, _, distance_to_in = find_closest_point_of_polygon(x, y, x_in, y_in)
                        
                    # DISTANCE TO THE EXTERNAL CONTOUR
                    nearest_xp, nearest_yp, distance_to_out = find_closest_point_of_polygon(x, y, x_out, y_out)     
                        
                    llon, llat = ccrs.Geodetic().transform_point(
                        nearest_xp, nearest_yp, proj_work
                    )
                    ii = np.where(
                        abs(llon - longitudes) == np.min(abs(llon - longitudes))
                    )[0]
                    jj = np.where(
                        abs(llat - latitudes) == np.min(abs(llat - latitudes))
                    )[0]
                    base_resolution = result[jj[0], ii[0]]

                    result[j, i] = resolution + (base_resolution - resolution) * (
                        distance_to_in / (distance_to_in + distance_to_out)
                    )

    return result


# THE BLOCK WHERE THE BACKGROUND RESOLUTIONS ARE SET, REFINED ALONG COASTLINES,
# REFINED WITHINN KML POLYGONS AND __CAN_BE__ ADJUSTED TO BATHYMETRY IF NEEDED


def define_resolutions(settings):

    latitudes = np.linspace(-90, 90, settings["n_latitudes"] + 1)  # 180*16
    longitudes = np.linspace(
        -180, 180, np.round(settings["n_longitudes"]).astype(int) + 1
    )  # np.round(360*16/5.75).astype(int)

    # base_resolution = 111
    base_resolution = settings["base_resolution"]

    result = np.full([len(latitudes), len(longitudes)], base_resolution).astype(float)

    # resolution in the entire Arctic above 60N
    # arctic_resolution = 25
    # arctic_low = 60
    # arctic_high = 65
    # for i in range(0, len(longitudes)):
    #    for j in range(0, len(latitudes)):
    #        if (latitudes[j] > arctic_low) & (latitudes[j] < arctic_high):
    #            result[j, i] = base_resolution - (base_resolution - arctic_resolution)*(latitudes[j] - arctic_low)/5
    #        elif latitudes[j] >= arctic_high:
    #            result[j, i] = arctic_resolution

    # resolution in the entire Arctic above 60N
    # for j in range(0, len(latitudes)):
    #     result[j, :] = base_resolution * np.cos(np.deg2rad(latitudes[j]))
    #     if abs(latitudes[j]) >= 77:
    #         result[j, :] = 25

    if settings["mercator_resolution"]["do_mercator_refinement"]:
        for j in range(0, len(latitudes)):
            result[j, :] = base_resolution * np.cos(np.deg2rad(latitudes[j]))
            if latitudes[j] >= settings["mercator_resolution"]["norhtern_boundary"]:
                result[j, :] = settings["mercator_resolution"][
                    "norhtern_lowest_resolution"
                ]
            elif latitudes[j] <= settings["mercator_resolution"]["southern_boundary"]:
                result[j, :] = settings["mercator_resolution"][
                    "southern_lowest_resolution"
                ]

    # refine resolution along coastlines
    
    if settings["refine_along_coastlines"]["do_refine_along_coastlines"]:

        min_resolution = settings["refine_along_coastlines"]["min_resolution"]
        max_distance = settings["refine_along_coastlines"]["max_distance"]
        min_length = settings["refine_along_coastlines"]["min_length"]
        averaging = settings["refine_along_coastlines"]["averaging"]

        result = refine_along_coastlines(
            longitudes,
            latitudes,
            result,
            min_resolution=min_resolution,
            max_distance=max_distance,
            min_length=min_length,
            averaging=averaging,
        )

    
    regions = []
    
    if settings["do_refine_within_regions"]:
        
        # create polygons from .kml files
    
        for reg in settings["regions"]:
            internal, outer = from_kml( reg["path"] )
            regions.append({})
            regions[-1]["name"] = reg["name"]
            regions[-1]["Polygon inside"] = internal[0], internal[1]
            regions[-1]["Polygon outside"] = outer[0], outer[1]
            regions[-1]["resolution"] = reg["resolution"]
            # regions[-1]["precision"] = reg["precision"]
    
        # refining within the regions
        
        for region in regions:
            result = refine(region, longitudes, latitudes, result)

    # bathymetry_adjustment()
    
    if settings["refine_over_bathymetry"]["do_refine_over_bathymetry"]:

        result = bathymetry_adjustment(settings, latitudes, longitudes, result)
        

    # saving
    with open("_result_temp.pkl", "wb") as file:
        pickle.dump([regions, result, latitudes, longitudes], file)

    return result, regions, longitudes, latitudes


# JIGSAW TRIANGULATION IS CALLED BY THIS FUNCTION


def triangulation(src_path, dst_path, longitudes, latitudes, result):

    opts = jigsawpy.jigsaw_jig_t()

    geom = jigsawpy.jigsaw_msh_t()
    hmat = jigsawpy.jigsaw_msh_t()
    mesh = jigsawpy.jigsaw_msh_t()

    # ------------------------------------ setup files for JIGSAW

    opts.geom_file = os.path.join(src_path, "geom.msh")
    opts.jcfg_file = os.path.join(dst_path, "log.jig")
    opts.mesh_file = os.path.join(dst_path, "mesh.msh")
    opts.hfun_file = os.path.join(dst_path, "hfun.msh")

    # ------------------------------------ define JIGSAW geometry

    geom.mshID = "ellipsoid-mesh"
    geom.radii = 6371 * np.ones(3)
    geom.ndims = +2

    # ------------------------------------ add coastline edges to geometry

    # inp1, inp2 = add_coastline()
    # geom.vert2 = np.array(inp1, dtype = geom.VERT2_t)
    # geom.edge2 = np.array(inp2, dtype = geom.EDGE2_t)

    jigsawpy.savemsh(opts.geom_file, geom)

    # ------------------------------------ compute HFUN over GEOM

    hmat.mshID = "ellipsoid-grid"
    hmat.ndims = +2
    hmat.radii = geom.radii
    hmat.xgrid = np.array(longitudes * np.pi / 180, dtype=hmat.REALS_t)
    hmat.ygrid = np.array(latitudes * np.pi / 180, dtype=hmat.REALS_t)

    hmat.value = np.array(np.array(result), dtype=hmat.REALS_t)

    jigsawpy.savemsh(opts.hfun_file, hmat)

    # ------------------------------------ build mesh via JIGSAW!

    print("")
    print("Call libJIGSAW")

    opts.hfun_scal = "absolute"

    opts.mesh_dims = +2  # 2-dim. simplexes
    # opts.mesh_kern = "delfront"         # DELFRONT kernel
    opts.optm_qlim = +0.95
    opts.optm_iter = +32

    opts.mesh_top1 = True  # for sharp feat's
    opts.geom_feat = True

    # jigsawpy.lib.jigsaw(opts, geom, mesh, None, hmat)

    jigsawpy.cmd.jigsaw(opts, mesh)

    #print("")
    #print("Saving .vtk file")
    #jigsawpy.savevtk(os.path.join(dst_path, "_result.vtk"), mesh)

    # saving mesh file
    with open("_mesh_temp.pkl", "wb") as file:
        pickle.dump(mesh, file)

    return mesh


# CUTTING OFF THE LAND (POSITIVE TOPOGRAPHY) FROM THE RESULTING JIGSAW MESH


def cut_land(settings, mesh):
    
    # TRANSFORMING CARTESIAN MESH TO LONGITUDES AND LATITUDES -----------------

    print("")
    print("transform jigsaw cartesian mesh to Lon-Lat")

    radii = 6371 * np.ones(3)

    ans = jigsawpy.tools.mathutils.R3toS2(radii, mesh.vert3["coord"])
    lon = np.rad2deg(ans[:, 0])
    lat = np.rad2deg(ans[:, 1])
    coastnode = np.full(len(lon), 0)
    depths = np.full(len(lon), -1)
    triangles = np.copy(mesh.tria3["index"])
    

    # CREATING DICTIONARY OF POINTS WITH INDICES OF LINKED TRIANGLES ----------

    tria_in_node = create_tria_in_node_dictionary(triangles)
    
    # DELETING VERTICES WITH POSITIVE ELEVATIONS AND LINKED TRIANGLES ---------

    print("")
    print("reading global topography")

    topo = Find_topo(settings)

    print("")
    print("deleting triangles over land: step 1 of 2")

    for index in range(0, len(lon)):

        if (index % 1000 == 0) | (index == len(lon) - 1):
            printProgressBar(
                index + 1, len(lon), prefix=f"Progress:", suffix="Complete", length=50
            )

        [status, xb, yb, dx, dy, run] = topo.ga(lon[index], lat[index])

        depths[index] = interpolate_lonlat(
            lat[index], xb, yb, dx, dy, topo.topo, min_points=3
        )

    print("")
    print("deleting triangles over land: step 2 of 2")

    for node in range(0, len(lon)):

        if (node % 1000 == 0) | (node == len(lon) - 1):
            printProgressBar(
                node + 1, len(lon), prefix=f"Progress:", suffix="Complete", length=50
            )

        if depths[node] >= 0:                    # land
            depths[node] = -1
            for j in tria_in_node[node]:

                for i in (0, 1, 2):
                    if triangles[j, i] > -1:
                        coastnode[triangles[j, i]] = 1
                
                triangles[j, :] = -1
                
        elif depths[node] >= -settings["minimum_depth"]:        # shallow
            depths[node] = settings["minimum_depth"]
            
        else:                                     # ocean
            depths[node] = abs(depths[node])

            

    # DELETING "LOOSE" TRIANGLES ----------------------------------------------

    print("")
    print('deleting "loose" elements')

    
    loose = True
    max_triangles = 6
    

    while loose:

        tria_in_node = create_tria_in_node_dictionary(triangles)
        total = 0

        loose_nodes = [n for n in tria_in_node if len(tria_in_node[n]) <= max_triangles]

        if max_triangles > 2:
            max_triangles -= 1

        if len(loose_nodes) > 0:
            for node in loose_nodes:

                # find all other nodes
                nodes = []
                for j in tria_in_node[node]:
                    for i in (0, 1, 2):
                        nodes.append(triangles[j, i])

                nodes = len(set(nodes))

                if (
                    ((len(tria_in_node[node]) == 1) & (nodes == 3))
                    | ((len(tria_in_node[node]) == 2) & (nodes == 5))
                    | (
                        (len(tria_in_node[node]) >= 3)
                        & (nodes >= len(tria_in_node[node]) + 3)
                    )
                ):

                    total += 1
                    depths[node] = -1

                    for j in tria_in_node[node]:
                        for i in (0, 1, 2):
                            if triangles[j, i] >= 0:
                                coastnode[triangles[j, i]] = 1

                        triangles[j, :] = -1

            if total == 0:
                loose = False

            print(f'{total} "loose" elements were deleted at this step')

        else:
            loose = False

# =============================================================================
#     DELETE NARROWS: ELEMENTS WITH ALL NODES == COASTNODES
#
#     pass 
# =============================================================================


    for node in range(len(lon)):
        n = 0
        if node in tria_in_node:
            for tria in tria_in_node[node]:
                if triangles[tria, 0] > -1: 
                    n += 1
        if n == 0:
            coastnode[node] = 0
            

    # FIND A STARTING NODE NEAR (OE, 0N) in Atlantic Ocean --------------------

    index = np.where((abs(lon) < 3) & (abs(lat) < 3))[0][0]

    # FINDING ALL NODES HAVING CONNECTION TO THE STARTING NODE ----------------
    # i.e. ELIMINATING LAKES --------------------------------------------------

    print("")
    print("eliminating inland lakes")

    active_triangles = []
    for i in tria_in_node[index]:
        active_triangles.append(i)

    dictionary = {}
    number_of_oceanic_nodes = len(depths[depths>-1])

    while active_triangles:

        tmp = []
        for node in triangles[active_triangles[-1]]:
            if (node > -1) & (node not in dictionary):
                dictionary[node] = 1
                for j in tria_in_node[node]:
                    tmp.append(j)

        active_triangles.pop()
        for j in tmp:
            active_triangles.append(j)
        
        if (len(dictionary) % 500 == 0):
            printProgressBar(
                len(dictionary), number_of_oceanic_nodes, 
                prefix="Progress:", suffix="Complete", length=50
            )
    
    printProgressBar( 10, 10, prefix="Progress:", suffix="Complete", length=50)
                    
    temp_set = set(dictionary)

    indices = []    
    for node in range(len(depths)):
        if node not in temp_set:
            depths[node] = -1
        else:
            indices.append(node)
    
    if settings["plot_region"]["do_plotting"]:
        print("")
        print("preparing figures")
        plotting(lon, lat, triangles, depths, coastnode, settings)
    
    
    # PREPARING FINAL ARRAYS --------------------------------------------------

    print("")
    print("re-sorting elements and nodes")
    
    lon_new = lon[indices]
    lat_new = lat[indices]
    depths_new = depths[indices]
    coastnode_new = coastnode[indices]

    shift = {}
    p = 0
    for index in np.arange(0, len(lon)):
        if index not in temp_set:
            p += 1
        else:
            shift[index] = p

    triangles_new = []
    for n, triangle in enumerate(triangles):
        if (
            (triangle[0] in temp_set)
            & (triangle[1] in temp_set)
            & (triangle[2] in temp_set)
        ):
            triangles_new.append(
                (
                    triangle[0] - shift[triangle[0]],
                    triangle[1] - shift[triangle[1]],
                    triangle[2] - shift[triangle[2]],
                )
            )

    # SAVING CONFIG FILES FOR FESOM -------------------------------------------

    print("")
    print("saving FESOM2 files")

    with open("elem2d.out", "w") as file:
        file.write(f"{len(triangles_new)}\n")
        for index in range(0, len(triangles_new)):
            file.write(
                f"{triangles_new[index][0]+1} {triangles_new[index][1]+1} {triangles_new[index][2]+1}\n"
            )

    with open("nod2d.out", "w") as file:
        file.write(f"{len(lon_new)}\n")
        for index in range(0, len(lon_new)):
            file.write(
                f"{index+1} {lon_new[index]} {lat_new[index]} {coastnode_new[index]}\n"
            )

    with open("aux3d.out", "w") as file:
        file.write(f'{int(len(settings["levels"]))}\n')
        for level in settings["levels"]:
            file.write(f"{level:.2f}\n")
        for index in range(0, len(lon_new)):
            file.write(f"{depths_new[index]:.2f}\n")

    # SAVING .VKT FILE FOR PARAVIEW -------------------------------------------

    print("")
    print("saving .vkt file for ParaView app")

    radii = 6371 * np.ones(3)
    ans = np.squeeze(np.dstack((np.deg2rad(lon_new), np.deg2rad(lat_new))))

    E3 = jigsawpy.tools.mathutils.S2toR3(radii, ans)

    with open("_result.vtk", "w") as file:
        file.write("# vtk DataFile Version 3.0\n")
        file.write("_result.vtk\n")
        file.write("ASCII\n")
        file.write("DATASET UNSTRUCTURED_GRID \n")
        file.write(f"POINTS {len(lon_new)} double \n")
        for index in range(0, len(lon_new)):
            file.write(f"{E3[index, 0]} {E3[index, 1]} {E3[index, 2]}\n")
        file.write(f"CELLS {len(triangles_new)} {4*len(triangles_new)}\n")
        for index in range(0, len(triangles_new)):
            file.write(
                f"3 {triangles_new[index][0]} {triangles_new[index][1]} {triangles_new[index][2]}\n"
            )
        file.write(f"CELL_TYPES {len(triangles_new)}\n")
        for index in range(0, len(triangles_new)):
            file.write("5\n")


def main():
    with open("./configure_example.yaml") as file:
        settings = yaml.load(file, Loader=yaml.FullLoader)

    print(settings["levels"])

    if (settings["use_existed_refinements"]) & (Path("_result_temp.pkl").exists()):
        with open("_result_temp.pkl", "rb") as file:
            regions, result, latitudes, longitudes = pickle.load(file)
    else:
        result, regions, longitudes, latitudes = define_resolutions(settings)

    
    estimate_number_of_nodes(result, longitudes, latitudes)
    
    if (settings["use_existed_jigsaw_mesh"]) & (Path("_mesh_temp.pkl").exists()):
        with open("_mesh_temp.pkl", "rb") as file:
            mesh = pickle.load(file)
    else:
        mesh = triangulation("jigsaw/", "jigsaw/", longitudes, latitudes, result)

    cut_land(settings, mesh)


if __name__ == "__main__":
    main()
