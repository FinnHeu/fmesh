# `fmesh`: A `FESOM` Mesh Generator

This is a simple mesh generation for regional `FESOM2` setups. Originally created by Sergey Kirillov.

We have a [video](https://youtu.be/8sxx6JwoLsg) demonstrating the basic process of defining regions and looking at the results. 
[<img src="https://user-images.githubusercontent.com/3407313/196788021-4a86fc91-0883-4e8d-8e4c-21c0b6bcdb65.png" width="50%">](https://youtu.be/8sxx6JwoLsg "Creating regional FESOM2 mesh with fmesh.")

## Traditional Installation

Create new `conda` environment from `requirements.yml` by executing:

```bash
conda env create -f requirements.yml
```
and activate it
```
conda activate fmesh
```

## Container Installation

You can pull the container and run with your favourite container software. This
has been tested with both `Docker` for local workstation use and `Singulairty`
for use on HPC systems. 

#### Downloading necessary data

You have to download bathymetry and coastlines. This can be done by executing:

```bash
./download_data.sh
```
## Usage

You have to configure your future mesh in the `configure.yaml` and then run:

```bash
python fmesh.py
```
in the `fmesh` directory.

## Container Usage

For `singularity`, you can run the image directly:

```
$ singularity run oras://ghcr.io/FESOM/fmesh.sif:latest
```

A similar approach can be done for `docker`

```
$ docker run docker://ghcr.io/FESOM/fmesh:latest
```

## Configuration

The process of mesh creation is split into several steps, that you can configure. In general you define the base resolution (ether constant, or increasing towards poles) and perform refinement on some regions of the basis mesh.

### Define base resolution

The main step is defining your base resolution at a principal grid. The first step is to set the number of latitudinal and longitudinal lines in the grid. A larger number of lines is required for the FESOM2 mesh with high resolution in certain regions (i.e. along coasts or within narrow straits). The general rule is to reach such a density of grid lines that quite a few grid nodes will happen to be within those regions.
Below is the example of refinement in the equatorial Atlantic with two different principal grids – 0.5 degree (360x720 lines) and 6 degrees (30x60 lines):

<figure to upload>
  
For defining base resolution in the nodes of principal grid there are two options:
- constant resolution (`do_mercator_refinement: false`)
- resolution changing with latitude (`do_mercator_refinement: true`)

For the resolution changing with latitude you can define the latitude where resolution become constant ( `norhtern_boundary`, `southern_boundary`) and the value of constant resolution (`norhtern_lowest_resolution`, `southern_lowest_resolution`).

### Refine resolution along coastlines

This option refines resolutions globally along coastlines of all continents and islands. To turn it on you have to set `do_refine_along_coastlines: true`. You can control:
- `min_resolution` resolution at the coast.
- `max_distance` distance from the coast, where resolution is set back to base resolution. So the resolution will change from `min_resolution` at the coast to base mesh resolution at `max_distance`.
- `min_length` minimum length of the coastline to be considered. For example, if it is set to 200 km the island with 200 km coastline will not be considered.
- `averaging` the smoothing length of the coastline. It’s not recommended to set this parameter too small. The 50 km or so is a good choice. A considerably smaller averaging parameter would help resolving small embayments in your global mesh, but may increase calculation time dramatically. For reaching finer resolutions in the specific geographical areas, it’s better to use “select regions with increased resolution” option rather than to refine coastal resolutions globally.
  
Note that the principal grid (see above) should be dense enough for this part producing a good consistent refinement. For example, if the mesh size of the principal grid is 100 km and the `max_distance` is also 100 km, the considered coastal zone will simply be at subgrid scale and the refinement won’t work properly.


### Select regions with increased resolution

Regions are defined in `.kml` files, that contain two paths for inner region (where the resolution will be set to constant `resolution` value) and outer region, where resolution will be set to the base resolution. In the area between the paths resolution will be linearly changing from `resolution` value to base resolution.

The easiest way to create paths is to use Google Earth, and path tool. The paths should not be closed, this will be done automatically by `fmesh`.
<img width="1607" alt="Screenshot 2022-10-19 at 21 29 04" src="https://user-images.githubusercontent.com/3407313/196785898-251c66ee-5822-4042-b157-f66318d09962.png">

The program applies all regions successively. Therefore, it is important to control an order of regions in your configuration file. The dependence of the final resolution from the order is demonstrated below.

<figure to upload>
<figure to upload>

There is an additional parameters which may be set for each region:
- `precision` each section of polygon will by split in `precision` number of
  times, to make the linear interpolation more precise. Using 10 is a good default value. See the difference between the final resolution maps for precision=10 and precision=3. 

<figure to upload>

### Select a region for a control plotting

If the parameter `do_plotting` is set to true, the program also creates two figures showing the resulting mesh (“_mesh.jpg”) and depths (“_depths.jpg”) within a region of interests limited by `min_lat`, `min_lon`, `max_lat` and `max_lon`. Be aware of using an appropriate cartopy projection name set by projection parameter and try to avoid depicting areas with very large amount of nodes. Otherwise, it may take long time to plot these figures.

### Using pre-calculated meshes 

There are two options to fasten calculations, especially during testing any changes in the program code.
By setting `use_existed_refinements` true, you may skip the refinements along coastlines and within regions with refined resolutions. A “_result_temp.pkl” file created before has to exist. Otherwise, the program will ignore this flag and refine the resolutions on the globe from a scratch.
By setting `use_existed_jigsaw_mesh` true, you may skip calling jigsaw mesh generator. A “_mesh_temp.pkl” file created before has to exist. Otherwise, the program will ignore this flag and the jigsaw mesh generator will be called.
  
### Set vertical levels

For convinience you can set vertical levels in your resulting `aux3d.out` file with `levels` parameter.

### First look at the mesh

The `fmesh` also saves a `.vtk file`, with the resulting mesh. This file can be opened in [ParaView](https://www.paraview.org/):
<img width="1920" alt="Screenshot 2022-10-19 at 21 40 35" src="https://user-images.githubusercontent.com/3407313/196788021-4a86fc91-0883-4e8d-8e4c-21c0b6bcdb65.png">

 
