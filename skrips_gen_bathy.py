# This script generates an MITgcm bathymetry file from ETOPO2 data.
# Update the case_name and file paths below if necessary.
# Author: Yuan-Jen Lin

####################################

save_bathy      = True
plot_bathy      = False

case_name       = "SKRIPS_20230301_0p01deg"
dir_model_input = f"/scratch/alpine/yuli3660/{case_name}/model_input/"
path_bathy      = f"{dir_model_input}bathymetry.bin"

path_geo_em     = f"/scratch/alpine/yuli3660/{case_name}/WPS/geo_em.d01.nc"

path_etopo      = "/projects/yuli3660/model_input/etopo2.nc"
dir_figs        = "/projects/yuli3660/figs/"

####################################


import numpy as np
import xarray as xr
import xesmf as xe
import matplotlib.pyplot as plt
import sys


###### 1. Load WRF info ######
# Rectangular grid (`lat-lon`)
# We can get `lon` by selecting any lat (e.g., 0), and `lat` by selecting any lon.

with xr.open_dataset(f'{path_geo_em}') as ds:
    wrf_lon = ds['XLONG_M'].isel({'south_north': 0}).squeeze()
    wrf_lat = ds['XLAT_M'].isel({'west_east': 0}).squeeze()
    wrf_landmask = ds['LANDMASK'].squeeze()
    wrf_lu_index = ds['LU_INDEX'].squeeze()


###### 2. Load ETOPO2 data ######
# Load data around the target region to speed up regridding and save memory.
# Add a buffer to lon/lat bounds to avoid NA values at the edges during regridding.
# Note: if the target region crosses the dateline, it selects all longitudes (x-direction).
# Note: data do not contain NA values (ds['z'].isnull().sum()=0)

buffer_lon, buffer_lat = 2, 2

with xr.open_dataset(f'{path_etopo}') as ds:
    etopo_z_full = ds['z']
    etopo_z = ds['z'].sel({'x': slice(wrf_lon.min()-buffer_lon, wrf_lon.max()+buffer_lon),
                           'y': slice(wrf_lat.min()-buffer_lat, wrf_lat.max()+buffer_lat)})
    etopo_x = etopo_z['x']
    etopo_y = etopo_z['y']


###### 3. ESMF regrid ######
# Interpolate ETOPO data from its original lon/lat to WRF lon/lat

grid_in   = {'lon': etopo_x, 'lat': etopo_y}
grid_out  = {'lon': wrf_lon, 'lat': wrf_lat}

regridder = xe.Regridder(grid_in, grid_out, 'bilinear',
                         periodic=False,
                         ignore_degenerate=True)
etopo_z_re= regridder(etopo_z, skipna=True, na_thres=1.)


###### 4. Mismatch of coastline between WRF and ETOPO data ######
# Positive/zero `etopo` OR land in WRF OR lake in WRF -> set to zero
# Positive/zero `etopo` AND ocean in WRF (lake excluded) -> set to min depth

min_depth = -1.
condition_zero_depth = (etopo_z_re>=0) | (wrf_landmask==1) | (wrf_lu_index==21)
condition_min_depth  = (etopo_z_re>=0) & (wrf_landmask==0) & (wrf_lu_index!=21)

etopo_z_cond1 = xr.where(condition_zero_depth, 0.,        etopo_z_re)
etopo_z_cond2 = xr.where(condition_min_depth,  min_depth, etopo_z_cond1)


###### 5. Save (optional) ######
# >: big-endian
# f4: 32-bit float

if save_bathy:
    f     = open(f'{path_bathy}', 'wb')
    bathy = etopo_z_cond2.values
    bathy.astype('>f4').tofile(f)
    f.close()


###### 6. Plot (optional) ######
# Sanity check

if plot_bathy:
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(20, 15))
    plt.subplots_adjust(wspace=0.1, hspace=0.2, left=0.15, right=0.9, bottom=0.05, top=0.95)

    i=0,0
    wrf_landmask.plot(ax=axes[i], extend='both', cmap='Blues_r')
    # wrf_landmask.plot.contour(ax=axes[i], levels=[0], colors='k')
    axes[i].set_title('a. WRF Landmask')

    if np.isin([-180., 180.], etopo_x).all():
        # Target region crosses the dateline. Plot global map.
        etopo_z_region = etopo_z_full.copy()
    else:
        etopo_z_region = etopo_z.sel({'x': slice(wrf_lon.min(), wrf_lon.max()),
                                      'y': slice(wrf_lat.min(), wrf_lat.max())})

    levels = np.arange(-8000, 4001, 100)
    
    i=1,0
    etopo_z_region.plot(ax=axes[i], levels=levels, extend='both', cmap='jet')
    etopo_z_region.plot.contour(ax=axes[i], levels=[0], colors='k')
    axes[i].set_title('b. ETOPO2 (original grid)')   
    i=2,0
    etopo_z_re.plot(ax=axes[i], levels=levels, extend='both', cmap='jet')
    # etopo_z_re.plot.contour(ax=axes[i], levels=[0], colors='k')
    axes[i].set_title('c. ETOPO2 (WRF grid)')
    i=0,1
    (etopo_z_cond1 - etopo_z_re).plot(ax=axes[i], extend='both', cmap='seismic')
    axes[i].set_title('d. Condition 1 impact')
    i=1,1
    (etopo_z_cond2 - etopo_z_cond1).plot(ax=axes[i], extend='both', cmap='Reds_r')
    axes[i].set_title('e. Condition 2 impact')
    i=2,1
    etopo_z_cond2.plot(ax=axes[i], levels=levels, extend='both', cmap='jet')
    # etopo_z_cond2.plot.contour(ax=axes[i], levels=[-1e-4], colors='k', linestyles='solid')
    axes[i].set_title('f. ETOPO2 (WRF grid; mismatch coastline solved; final bathy)')
    plt.savefig(f'{dir_figs}bathymetry.png', dpi=300)
    plt.close()