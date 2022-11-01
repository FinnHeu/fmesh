# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 13:00:19 2022

@author: SKirillov
"""

import os
import pickle
import sys
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cf
import geopy.distance as dist
import jigsawpy
import netCDF4 as nc
import numpy as np
import yaml
import multiprocessing
from fastkml import geometry, kml

from utils import (gradapproach, inside_outside, interpolate_lonlat,
                   printProgressBar, read_coastlines2)

# CREATING THE DICTIONARY OF NODES WITH INDICES OF LINKED ELEMENTS ------------
# i.e. { 234: [23, 24, 265, 4350, 4352], 235: [23, 24, 123, 990], etc. }


def create_tria_in_node_dictionary(triangles):

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

    if not order:
        internal, outer = outer, internal

    return internal, outer


# THE FUNCTION FOR ADJUSTING RESOLUTIONS TO THE BATHYMETRY
#    Need to be called if needed!


def bathymetry_adjustment(settings, latitudes, longitudes, result):

    print("read topography")
    topo = Find_topo(settings)

    minimum = np.min(topo.topo)
    maximum = np.max(topo.topo)

    for j, lat in enumerate(latitudes):
        printProgressBar(
            j + 1, len(latitudes), prefix=f"Progress:", suffix="Complete", length=50
        )

        for i, lon in enumerate(longitudes):

            [status, xb, yb, dx, dy, run] = topo.ga(lon, lat)

            depth = interpolate_lonlat(lat, xb, yb, dx, dy, topo.topo, min_points=1)
            if depth > 0:
                result[j, i] = 150 - 149 * (depth - 0) / (maximum - 0)
    return result


def init_multiprocessing(args):
    global counter
    counter = args
    

def inner(parameters):
    
    import geopy.distance as dist
    
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
                     prefix = f'Progress:', suffix = 'Complete', length = 50) 
    
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
    print(f'preparing smoothed coastlines')
    
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

    (poly_in, poly_out, resolution, precision) = (
        region["Polygon inside"],
        region["Polygon outside"],
        region["resolution"],
        region["precision"],
    )

    min_lon = np.min(poly_out[0])
    max_lon = np.max(poly_out[0])
    min_lat = np.min(poly_out[1])
    max_lat = np.max(poly_out[1])

    # checking if the current polygon encircles poles for adjusting min_lat and max_lat
    proj_work = ccrs.Orthographic(central_longitude=0, central_latitude=-90)
    x, y = proj_work.transform_point(0, -90, ccrs.Geodetic())
    (x_out, y_out, _) = np.hsplit(
        proj_work.transform_points(ccrs.Geodetic(), poly_out[0], poly_out[1]), 3
    )
    if abs(inside_outside(x_out, y_out, x, y)) == 1:
        min_lat = -90
        min_lon = -180
        max_lon = 180
    proj_work = ccrs.Orthographic(central_longitude=0, central_latitude=90)
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
    print(region["name"])
    print(
        f"{min_lon:.1f}   ", f"{max_lon:.1f}   ", f"{min_lat:.1f}   ", f"{max_lat:.1f}"
    )

    for i, grid_lon in enumerate(longitudes):

        printProgressBar(
            i + 1, len(longitudes), prefix=f"Progress:", suffix="Complete", length=50
        )

        for j, grid_lat in enumerate(latitudes):

            if ((variant == 1) & (min_lon <= grid_lon <= max_lon) & 
                (min_lat <= grid_lat <= max_lat)) | \
                ((variant == 2) & (max_lon <= grid_lon) & (min_lat <= grid_lat <= max_lat)) | \
                ((variant == 2) & (min_lon >= grid_lon) & (min_lat <= grid_lat <= max_lat)):

                proj_work = ccrs.Orthographic(
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
                    distance_to_in = 1e10

                    for index in range(0, len(x_in) - 1):
                        dist = np.sqrt((x_in[index] - x) ** 2 + (y_in[index] - y) ** 2)
                        for (xp, yp) in zip(
                            np.linspace(x_in[index], x_in[index + 1], precision),
                            np.linspace(y_in[index], y_in[index + 1], precision),
                        ):
                            dist = np.sqrt((xp - x) ** 2 + (yp - y) ** 2)
                            if dist < distance_to_in:
                                distance_to_in = dist

                    # DISTANCE TO THE EXTERNAL CONTOUR
                    distance_to_out = 1e10

                    #   dist = np.sqrt((x_out[:-1]-x)**2 + (y_out[:-1]-y)**2)
                    #   nearest_corner = np.where(dist == np.min(dist))[0][0]
                    #   if nearest_corner == 0: nc2 = len(x_out)-2
                    #   else: nc2 = nearest_corner - 1

                    for index in range(0, len(x_out) - 1):
                        for (xp, yp) in zip(
                            np.linspace(x_out[index], x_out[index + 1], precision),
                            np.linspace(y_out[index], y_out[index + 1], precision),
                        ):
                            dist = np.sqrt((xp - x) ** 2 + (yp - y) ** 2)
                            if dist < distance_to_out:
                                distance_to_out = dist
                                (nearest_xp, nearest_yp) = (xp, yp)

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

    min_resolution = settings["refine_along_coastlines"]["min_resolution"]
    max_distance = settings["refine_along_coastlines"]["max_distance"]
    min_length = settings["refine_along_coastlines"]["min_length"]
    averaging = settings["refine_along_coastlines"]["averaging"]

    if settings["refine_along_coastlines"]["do_refine_along_coastlines"]:
        result = refine_along_coastlines(
            longitudes,
            latitudes,
            result,
            min_resolution=min_resolution,
            max_distance=max_distance,
            min_length=min_length,
            averaging=averaging,
        )

    # create polygons from .kml files
    regions = []
    # from_kml('_kml/CAA.kml','CAA', res=5, precision=10, order=True)
    # from_kml('./kml/Baffin_Bay.kml','Baffin', res=5, precision=20, order=True)
    # from_kml('./kml/fram.kml','fram', res=5, precision=20, order=True)
    # from_kml('./kml/Nares.kml','Nares', res=1, precision=20, order=True)
    # from_kml('./kml/Peabody.kml','Peabody', res=0.18, precision=20, order=True)
    # from_kml('_kml/Cardigan.kml','Cardigan', res=0.2, precision=10, order=True)
    # from_kml('_kml/Fury.kml','Fury', res=0.2, precision=10, order=True)
    # from_kml('./kml/Denmark.kml','Denmark', res=4, precision=10, order=True)
    # from_kml('./kml/Gibraltar.kml','Gibraltar', res=4, precision=10, order=True)

    for reg in settings["regions"]:
        internal, outer = from_kml(
            reg["path"],
            #    reg['name'],
            #    res=reg['resolution'],
            #    precision=reg['precision'],
            order=reg["order"],
        )
        regions.append({})
        regions[-1]["name"] = reg["name"]
        regions[-1]["Polygon inside"] = internal[0], internal[1]
        regions[-1]["Polygon outside"] = outer[0], outer[1]
        regions[-1]["resolution"] = reg["resolution"]
        regions[-1]["precision"] = reg["precision"]

    # refining
    for region in regions:
        result = refine(region, longitudes, latitudes, result)

    # bathymetry_adjustment()
    pass

    # saving
    with open("_result_temp.pkl", "wb") as file:
        pickle.dump([regions, result, latitudes, longitudes], file)

    with open("_config.txt", "w") as file:
        file.write(f"amount of meridional nodes = {len(latitudes)}\n")
        file.write(f"amount of zonal nodes = {len(longitudes)}\n")
        file.write(f"base resolution = {base_resolution} km\n")
        file.write("\n")

        file.write("Coastlines:\n")
        file.write(f"min_resolution = {min_resolution} km\n")
        file.write(f"max_distance = {max_distance} km\n")
        file.write(f"min_length = {min_length} km\n")
        file.write(f"averaging = {averaging} km\n")
        file.write("\n")

        for region in regions:
            file.write(f'Region: {region["name"]}\n')
            file.write(f'   resolution: {region["resolution"]} km\n')
            file.write(f'   precision: {region["precision"]}\n')

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

    print("")
    print("Saving .vtk file")
    jigsawpy.savevtk(os.path.join(dst_path, "_result.vtk"), mesh)

    # saving mesh file
    with open("_mesh_temp.pkl", "wb") as file:
        pickle.dump(mesh, file)

    return mesh


# CUTTING OFF THE LAND (POSITIVE TOPOGRAPHY) FROM THE RESULTING JIGSAW MESH


def cut_land(settings, mesh, depth_limit=-20):

    # TRANSFORMING CARTESIAN MESH TO LONGITUDES AND LATITUDES -----------------

    print("")
    print("transform jigsaw cartesian mesh to Lon-Lat")

    radii = 6371 * np.ones(3)

    ans = jigsawpy.tools.mathutils.R3toS2(radii, mesh.vert3["coord"])
    lon = np.rad2deg(ans[:, 0])
    lat = np.rad2deg(ans[:, 1])
    coastnode = np.full(len(lon), 0)
    depths = np.full(len(lon), -1)
    depth_codes = np.full(len(lon), 0)
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

    for index in range(0, len(lon)):

        if (index % 1000 == 0) | (index == len(lon) - 1):
            printProgressBar(
                index + 1, len(lon), prefix=f"Progress:", suffix="Complete", length=50
            )

        if depths[index] >= 0:

            n = 0
            for j in tria_in_node[index]:
                for i in (0, 1, 2):
                    if (depths[triangles[j, i]] >= 0) & (triangles[j, i] != index):
                        n += 1

            if n >= 2:
                depth_codes[index] = 1  # land
                for j in tria_in_node[index]:

                    for i in (0, 1, 2):
                        if triangles[j, i] > 0:
                            coastnode[triangles[j, i]] = 1

                    # coastnode[ triangles[j, :] ] = 1
                    triangles[j, :] = -1
            else:
                depth_codes[index] = 2  # shallow

        elif depths[index] >= depth_limit:
            depth_codes[index] = 2  # shallow

        else:
            depth_codes[index] = 3  # ocean

    for index in range(0, len(lon)):
        if depth_codes[index] == 1:
            depths[index] = -1
        elif depth_codes[index] == 2:
            depths[index] = abs(depth_limit)
        else:
            depths[index] = abs(depths[index])

    # DELETING "LOOSE" TRIANGLES ----------------------------------------------

    print("")
    print('deliting "loose" triangles')

    loose = True
    max_triangles = 6

    while loose:

        tria_in_node = create_tria_in_node_dictionary(triangles)
        total = 0

        loose_nodes = [n for n in tria_in_node if len(tria_in_node[n]) <= max_triangles]

        if max_triangles > 2:
            max_triangles -= 1

        if len(loose_nodes) > 0:
            for index in loose_nodes:

                # find all other nodes
                nodes = []
                for j in tria_in_node[index]:
                    for i in (0, 1, 2):
                        nodes.append(triangles[j, i])

                nodes = len(set(nodes))

                if (
                    ((len(tria_in_node[index]) == 1) & (nodes == 3))
                    | ((len(tria_in_node[index]) == 2) & (nodes == 5))
                    | (
                        (len(tria_in_node[index]) >= 3)
                        & (nodes >= len(tria_in_node[index]) + 3)
                    )
                ):

                    total += 1
                    depths[index] = -1

                    for j in tria_in_node[index]:
                        for i in (0, 1, 2):
                            if triangles[j, i] > 0:
                                coastnode[triangles[j, i]] = 1

                        # coastnode[ triangles[j, :] ] = 1
                        triangles[j, :] = -1

            if total == 0:
                loose = False

            print(f'{total} "loose" triangles have been deleted at this step')

        else:
            loose = False

    # FIND A STARTING NODE NEAR (OE, 0N) in Atlantic Ocean --------------------

    index = np.where((abs(lon) < 3) & (abs(lat) < 3))[0][0]

    # FINDING ALL NODES HAVING CONNECTION TO THE STARTING NODE ----------------
    # i.e. ELIMINATING LAKES --------------------------------------------------

    print("")
    print("eliminating lakes")

    active_triangles = []
    for i in tria_in_node[index]:
        active_triangles.append(i)

    dictionary = {}

    while active_triangles:

        tmp = []
        for i in triangles[active_triangles[-1]]:
            if (i > -1) & (i not in dictionary):
                dictionary[i] = 1
                for j in tria_in_node[i]:
                    tmp.append(j)

        active_triangles.pop()
        for j in tmp:
            active_triangles.append(j)

    # PREPARING FINAL ARRAYS --------------------------------------------------

    print("")
    print("re-sorting vertices and triangles")

    temp_set = set(dictionary)

    lon_new = lon[list(temp_set)]
    lat_new = lat[list(temp_set)]
    depths_new = depths[list(temp_set)]
    coastnode_new = coastnode[list(temp_set)]

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
    print("saving .vkt file")

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
    with open("./configure.yaml") as file:
        settings = yaml.load(file, Loader=yaml.FullLoader)

    print(settings["levels"])

    if Path("_result_temp.pkl").exists():
        with open("_result_temp.pkl", "rb") as file:
            regions, result, latitudes, longitudes = pickle.load(file)
    else:
        result, regions, longitudes, latitudes = define_resolutions(settings)

    if Path("_mesh_temp.pkl").exists():
        with open("_mesh_temp.pkl", "rb") as file:
            mesh = pickle.load(file)
    else:
        mesh = triangulation("jigsaw/", "jigsaw/", longitudes, latitudes, result)

    cut_land(settings, mesh, depth_limit=-30)


if __name__ == "__main__":
    main()
