# This script generates IC (initial condition) and OBCS (open boundary conditions)
# of variables U, V, T, S for MITgcm from the GLORYS12V1 product.
# Update the case_name and file paths below if necessary.
# Author: Yuan-Jen Lin

####################################

save_ic_obcs    = True
plot_ic_obcs    = False

case_name       = "SKRIPS_20230301_0p01deg"
dir_model_input = f"/scratch/alpine/yuli3660/{case_name}/model_input/"
dir_figs        = "/projects/yuli3660/figs/"

path_geo_em     = f"/scratch/alpine/yuli3660/{case_name}/WPS/geo_em.d01.nc"
path_namelist   = f"/scratch/alpine/yuli3660/{case_name}/WPS/namelist.wps"

path_dz         = "/projects/yuli3660/model_input/delZ_60v.txt"

# Check the end date of `GLOBAL_MULTIYEAR_PHY_001_030`:
# https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/services
# -> Update `my_end_date` accordingly.

my_end_date     = '2025-12-22'

####################################


import numpy as np
import xarray as xr
import xesmf as xe
import matplotlib.pyplot as plt
import sys
import copernicusmarine
import f90nml
from scipy.io import loadmat
from datetime import datetime, timedelta
import pyinterp.fill
import pyinterp.backends.xarray


def get_z(path_dz):
    
    with open(path_dz, 'r') as f:
        dz = np.array([float(line.strip().rstrip(',')) for line in f])
    zf = np.cumsum(dz)
    zc = 0.5 * (zf + np.concatenate([[0.], zf])[:-1])
    
    return zc


###### 1. Load WRF info ######
# Rectangular grid (`lat-lon`)
# We can get `lon` by selecting any lat (e.g., 0), and `lat` by selecting any lon.

with xr.open_dataset(f'{path_geo_em}') as ds:
    wrf_lon = ds['XLONG_M'].isel({'south_north': 0}).squeeze()
    wrf_lat = ds['XLAT_M'].isel({'west_east': 0}).squeeze()
    wrf_landmask = ds['LANDMASK'].squeeze()
    wrf_lu_index = ds['LU_INDEX'].squeeze()

wrf_nml        = f90nml.read(path_namelist)
wrf_start_date = wrf_nml['share']['start_date'][:10] # yyyy-mm-dd_hh:mm:ss
wrf_end_date   = wrf_nml['share']['end_date'][:10]   # yyyy-mm-dd_hh:mm:ss
wrf_ref_lon    = wrf_nml['geogrid']['ref_lon']
wrf_ref_lat    = wrf_nml['geogrid']['ref_lat']
wrf_e_we       = wrf_nml['geogrid']['e_we']
wrf_e_sn       = wrf_nml['geogrid']['e_sn']
wrf_dx         = wrf_nml['geogrid']['dx']
wrf_dy         = wrf_nml['geogrid']['dy']

# Extend WRF end date by one day for MITgcm OBCS boundary condition
wrf_end_date   = (datetime.strptime(wrf_end_date, "%Y-%m-%d") + timedelta(days=1)).strftime("%Y-%m-%d")

print(f"wrf_start_date: {wrf_start_date}")
print(f"wrf_end_date: {wrf_end_date}")

# OVERWRITE IF NEEDED
# wrf_start_date = '2023-03-01'
# wrf_end_date = '2023-03-01'


###### 2. Z-axis ######

zc = get_z(path_dz)


###### 3. Load T,S,U,V data ######
# Sources:
# a) Operational Mercator global ocean analysis and forecast system
# URL: https://data.marine.copernicus.eu/product/GLOBAL_ANALYSISFORECAST_PHY_001_024/description
# b) GLORYS12V1 product; the CMEMS global ocean eddy-resolving (1/12° horizontal resolution, 50 vertical levels) reanalysis covering the altimetry (1993 onward).
# URL: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description

if not copernicusmarine.login(check_credentials_valid=True):
    copernicusmarine.login(force_overwrite=True)

iv_in_ls  = ['thetao', 'so', 'uo', 'vo']
iv_out_ls = ['T', 'S', 'U', 'V']

iv_in_ls  = ['vo',]
iv_out_ls = ['V',]

wrf_end_date_compare    = datetime.strptime(wrf_end_date, "%Y-%m-%d")
my_end_date_compare     = datetime.strptime(my_end_date, "%Y-%m-%d")

ds_id_ls = {}

if wrf_end_date_compare <= my_end_date_compare:
    print(f'Use multi-year (`my`) data')
    for iv in iv_in_ls:
        ds_id_ls[iv] = "cmems_mod_glo_phy_my_0.083deg_P1D-m"
else:
    print(f'Use anfc data')
    ds_id_ls['thetao'] = "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m"
    ds_id_ls['so']     = "cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m"
    ds_id_ls['uo']     = "cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m"
    ds_id_ls['vo']     = "cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m"

buffer_lon, buffer_lat = 2, 2

ds = {}

for ivN, iv in enumerate(iv_in_ls):

    ds_request = copernicusmarine.open_dataset(
        dataset_id        = ds_id_ls[iv],
        minimum_longitude = wrf_ref_lon - 0.5*wrf_e_we*wrf_dx - buffer_lon,
        maximum_longitude = wrf_ref_lon + 0.5*wrf_e_we*wrf_dx + buffer_lon,
        minimum_latitude  = wrf_ref_lat - 0.5*wrf_e_sn*wrf_dy - buffer_lat,
        maximum_latitude  = wrf_ref_lat + 0.5*wrf_e_sn*wrf_dy + buffer_lat,
        start_datetime    = wrf_start_date,
        end_datetime      = wrf_end_date,
        variables         = [iv]
    )
    
    grid_in     = {'lon': ds_request[iv]['longitude'], 'lat': ds_request[iv]['latitude']}
    grid_out    = {'lon': wrf_lon, 'lat': wrf_lat}

    regridder   = xe.Regridder(grid_in, grid_out, 'bilinear',
                               periodic=False,
                               ignore_degenerate=True)
    ds_re = regridder(ds_request[iv], skipna=True, na_thres=1.)

    # Z interpolation
    ds[iv] = ds_re.interp({'depth': -1.*zc})

    # Fill NaN values
    # Rename axes and units to match pyinterp's requirements
    ds[iv]['XLAT_M'].attrs['units']='degrees_north'
    ds[iv]['XLONG_M'].attrs['units']='degrees_east'
    ds[iv] = ds[iv].rename({'south_north': 'XLAT_M', 'west_east': 'XLONG_M'})
    
    ds[f'{iv}_filled'] = ds[iv].copy()
    for itime in range(len(ds[iv]['time'])):
        grid                  = pyinterp.backends.xarray.Grid3D(ds[iv][itime,:,:,:])
        has_converged, filled = pyinterp.fill.gauss_seidel(grid)

        ds[f'{iv}_filled'][itime,:,:,:] = filled.transpose(2, 1, 0)

        if not has_converged:
            raise RuntimeError("Gauss-Seidel solver did not converge")
        # print(ds[iv])
        # print(ds[f'{iv}_filled'])

    # If zc (destination z) is deeper than the "second-last" ds_request['depth'][-2] (source z),
    # the corresponding levels will be all zeros.
    # The `min_level_zero` is the shallowest level with all-zero values.
    source_z          = ds_request['depth'].values
    destin_z          = ds[f'{iv}_filled']['depth'].values
    destin_z_too_deep = destin_z >= source_z[-2]
    min_level_zero    = -int(np.sum(destin_z_too_deep))
    if min_level_zero < 0:
        print(f"min_level_zero: {min_level_zero}")
        # The below one-line code is extremely slow for large arrays
        # ds[f'{iv}_filled'].data[:, min_level_zero:, :, :] = ds[f'{iv}_filled'].data[:, min_level_zero-1:min_level_zero, :, :]

        # To avoid heavy memory use, I used the following code instead: 
        # (filling the deeper levels by repeating the last non-zero level)
        data_dims   = ds[f'{iv}_filled'].dims
        data_bottom = ds[f'{iv}_filled'].isel({'depth': min_level_zero-1}).expand_dims({'depth': destin_z[min_level_zero:]}).transpose(*data_dims)
        data_top    = ds[f'{iv}_filled'].isel({'depth': slice(0, min_level_zero)})
        ds[f'{iv}_filled'] = xr.concat([data_top, data_bottom], dim='depth')
    else:
        print(f"There are no all-zero levels.")


    ###### 4. Save (optional) ######
    # >: big-endian
    # f4: 32-bit float

    iv_out  = iv_out_ls[ivN]
    lat_str = 'XLAT_M'
    lon_str = 'XLONG_M'
    
    if save_ic_obcs:
        # ic
        f      = open(f'{dir_model_input}{iv_out}_{wrf_start_date}_filled.bin', 'wb')
        out_iv = ds[f'{iv}_filled'].sel({'time': wrf_start_date}).squeeze().values
        out_iv.astype('>f4').tofile(f)
        f.close()

        # obcs
        f      = open(f'{dir_model_input}{iv_out}_{wrf_start_date}_{wrf_end_date}_obcs_S_filled.bin', 'wb')
        out_iv = ds[f'{iv}_filled'].isel({lat_str: 0}).squeeze().values
        out_iv.astype('>f4').tofile(f)
        f.close()

        f      = open(f'{dir_model_input}{iv_out}_{wrf_start_date}_{wrf_end_date}_obcs_N_filled.bin', 'wb')
        out_iv = ds[f'{iv}_filled'].isel({lat_str: -1}).squeeze().values
        out_iv.astype('>f4').tofile(f)
        f.close()

        f      = open(f'{dir_model_input}{iv_out}_{wrf_start_date}_{wrf_end_date}_obcs_W_filled.bin', 'wb')
        out_iv = ds[f'{iv}_filled'].isel({lon_str: 0}).squeeze().values
        out_iv.astype('>f4').tofile(f)
        f.close()

        f      = open(f'{dir_model_input}{iv_out}_{wrf_start_date}_{wrf_end_date}_obcs_E_filled.bin', 'wb')
        out_iv = ds[f'{iv}_filled'].isel({lon_str: -1}).squeeze().values
        out_iv.astype('>f4').tofile(f)
        f.close()


    ###### 5. Plot (optional) ######
    # Sanity check
    
    if plot_ic_obcs:
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(20, 15))
        plt.subplots_adjust(wspace=0.1, hspace=0.2, left=0.15, right=0.9, bottom=0.05, top=0.95)

        # ds[f'{iv}_filled']['XLONG_M'] = np.where(ds[f'{iv}_filled']['XLONG_M']<0, ds[f'{iv}_filled']['XLONG_M']+360., ds[f'{iv}_filled']['XLONG_M'])
        
        i=0,0
        plot_in = ds[f'{iv}_filled'].sel({'time': wrf_start_date}).squeeze()[0, :, :]
        plot_in.plot(ax=axes[i], levels=100, extend='both', cmap='jet')
        axes[i].set_title(f'a. {iv_out}: IC ({wrf_start_date})')

        axes[0,1].axis('off')
        
        i=1,0
        plot_in = ds[f'{iv}_filled'].isel({lat_str: 0}).squeeze()[:, 0, :]
        plot_in.plot(ax=axes[i], levels=100, extend='both', cmap='jet')
        axes[i].set_title(f'b. {iv_out}: S open boundary conditions (obcs)')
        i=1,1
        plot_in = ds[f'{iv}_filled'].isel({lat_str: -1}).squeeze()[:, 0, :]
        plot_in.plot(ax=axes[i], levels=100, extend='both', cmap='jet')
        axes[i].set_title(f'c. {iv_out}: N open boundary conditions (obcs)')
        i=2,0
        plot_in = ds[f'{iv}_filled'].isel({lon_str: 0}).squeeze()[:, 0, :]
        plot_in.plot(ax=axes[i], levels=100, extend='both', cmap='jet')
        axes[i].set_title(f'd. {iv_out}: W open boundary conditions (obcs)')
        i=2,1
        plot_in = ds[f'{iv}_filled'].isel({lon_str: -1}).squeeze()[:, 0, :]
        plot_in.plot(ax=axes[i], levels=100, extend='both', cmap='jet')
        axes[i].set_title(f'e. {iv_out}: E open boundary conditions (obcs)')
        plt.savefig(f'{dir_figs}{iv_out}.png', dpi=300)
        plt.close()
