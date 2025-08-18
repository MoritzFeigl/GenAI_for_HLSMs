#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : plot_parameter_field
Project Name: 2020_GenAI_para_mHM
Description : insert your description here, if applicable
Author      : ottor
Created     : 08.02.21 10:35
"""

# IMPORTS
from typing import Union, Iterable, Optional

# all external modules installed with conda
import cartopy.crs as ccrs  # used version 0.17.0
import matplotlib.pyplot as plt  # 3.2.1
import numpy as np  # 1.18.1
import xarray as xr  # 0.15.0
from rasterio.warp import transform  # 1.1.3
import rioxarray
import os

# GLOBAL VARIABLES
RENAME_DIMS = {
    'lon_out': 'x',
    'lat_out': 'y',
}


# FUNCTIONS
def _all_dims_contained(dims: Iterable[str], obj: Union[xr.DataArray, xr.Dataset]):
    """check if all dimensions of a list are containing in an xarray object"""
    return all([dim in obj.dims for dim in dims])


def read(filename, select_variable: Optional[Union[str, Iterable[str]]] = None, select_slice: Optional[dict] = None):
    """open a netcdf file from disk, optionally selecting variables or slices of dimensions"""
    with xr.open_dataset(filename) as ds:
        # select one or more variables
        if select_variable is not None:
            ds = ds[select_variable]
        # select a part of the coordinates
        if select_slice is not None and _all_dims_contained(select_slice.keys(), ds):
            ds = ds.sel(**select_slice)
        # rename the coordinates contained in RENAME_DIMS
        if _all_dims_contained(RENAME_DIMS.keys(), ds):
            ds = ds.rename(RENAME_DIMS)
        # load only the part needed and close the file handle
        ds.load()
    return ds


def select(ds: Union[xr.DataArray, xr.Dataset], isel: Optional[dict] = None, sel: Optional[dict] = None):
    """again, perform a selection for coordinates based on labels (sel) or indices (isel)"""
    selected = ds
    if isel:
        selected = selected.isel(**isel)
    if sel:
        selected = selected.sel(**sel)
    return selected


def prepare_tiff_out(ds: Union[xr.DataArray, xr.Dataset]):
    """Prepare for tiff generation, assumes 1-D input coordinates named x,y and creates lon, lat"""
    # Compute the lon/lat coordinates with rasterio.warp.transform
    ny, nx = len(ds['y']), len(ds['x'])
    x, y = np.meshgrid(ds['x'], ds['y'])
    # Rasterio works with 1D arrays
    lon, lat = x.flatten(), y.flatten()
    lon = np.asarray(lon).reshape((ny, nx))
    lat = np.asarray(lat).reshape((ny, nx))
    ds.coords['lon'] = (('y', 'x'), lon)
    ds.coords['lat'] = (('y', 'x'), lat)
    ds['lat'].attrs = {'standard_name': 'latitude', 'long_name': 'Latitude',
                       'units': 'degrees_north', '_CoordinateAxisType': 'Lat', 'axis': 'y'}
    ds['lon'].attrs = {'standard_name': 'longitude', 'long_name': 'Longitude',
                       'units': 'degrees_east', '_CoordinateAxisType': 'Lon', 'axis': 'x'}
    return ds


def print_simple(filename: str, ds: Union[xr.DataArray, xr.Dataset], x_coord: str = 'lon', y_coord: str = 'lat'):
    """print using the native xarray plotting, provide a 2D xarray object as ds"""
    ds.plot(x=x_coord, y=y_coord)
    plt.savefig(filename)
    plt.close()


def print_map(filename: str, ds: Union[xr.DataArray, xr.Dataset], x_coord: str = 'lon', y_coord: str = 'lat'):
    """print using the cartopy library, provide a 2D xarray object as ds"""
    # set the map projection
    ax = plt.axes(projection=ccrs.PlateCarree())
    # add coastlines to map, add other items if you like
    ax.coastlines()
    ds.plot.pcolormesh(ax=ax, x=x_coord, y=y_coord,)
    plt.savefig(filename)
    plt.close()


def export_to_geotiff(filename: str, ds: Union[xr.DataArray, xr.Dataset]):
    ds.rio.to_raster(filename)

# CLASSES


# load the dataset, optionally select variables or slices of dimensions
parameters = ["KSat_till", "FieldCap_till", "L1_FieldCap", "fRoots_temp", "L1_fRoots",
                        "ThetaS_till", "L1_SatSoilMoisture", "L1_Max_Canopy_Intercept"]
output_path = f"/gpfs/data/fs71468/GenAI_para_runs/analysis_results/validation_results/parameter_tiffs/"
os.makedirs(output_path, exist_ok=True)
plot_path = f"/gpfs/data/fs71468/GenAI_para_runs/analysis_results/validation_results/parameter_tiffs/plots/"
os.makedirs(plot_path, exist_ok=True)
hor = [2]
for horizon in hor:
    for version in ["exp0-run3", "exp2-run1"]:
        for parameter in [parameters[i] for i in [2, 3, 4, 6]]:
            val_path = f"/gpfs/data/fs71468/GenAI_para_runs/results/Validation/{version}/"
            train_path = f"/gpfs/data/fs71468/GenAI_para_runs/results/Training/{version}/"
            # loop over basins and train/val data
            dataArray_list = []
            min_lat = np.inf
            max_lat = -np.inf
            min_lon = np.inf
            max_lon = -np.inf
            for para_path in [val_path, train_path]:
                basins = [file.split("_")[0] for file in os.listdir(para_path) if "_parameters" in file]
                for basin in basins:
                    print(basin)
                    # load netcdf and parameter
                    filename = para_path + basin + '_parameters.nc'
                    full_array = read(filename, select_variable=parameter)
                    if parameter in ["KSat_till", "FieldCap_till", "ThetaS_till"]:
                        data_array = select(full_array, isel=dict(horizon_till=horizon, land_cover_period=2))
                    elif parameter in ["L1_FieldCap", "fRoots_temp", "L1_fRoots", "L1_SatSoilMoisture"]:
                        data_array = select(full_array, isel=dict(horizon_out=horizon, land_cover_period_out=2))
                        data_array = data_array.rename({'y': 'lat'})
                        data_array = data_array.rename({'x': 'lon'})
                    elif parameter == "L1_Max_Canopy_Intercept":
                        data_array = select(full_array, isel=dict(month_of_year=7))
                        data_array = data_array.rename({'y': 'lat'})
                        data_array = data_array.rename({'x': 'lon'})
                    # get extent
                    min_lat = min(min_lat, np.min(data_array.lat.values))
                    max_lat = max(max_lat, np.max(data_array.lat.values))
                    min_lon = min(min_lon, np.min(data_array.lon.values))
                    max_lon = max(max_lon, np.max(data_array.lon.values))
                    data_array = data_array.rename({'lat': 'y'})
                    data_array = data_array.rename({'lon': 'x'})
                    dataArray_list.append(prepare_tiff_out(data_array))
                    print(len(dataArray_list))
            # calculate resolution
            lat_resolution = abs(data_array.lat.values[0, 0] - data_array.lat.values[1, 0])
            lon_resolution = abs(data_array.lon.values[0, 0] - data_array.lon.values[0, 1])
            # Calculate extent for lon and lat
            lon_extent = max_lon - min_lon
            lat_extent = max_lat - min_lat
            # Calculate num for lon and lat using the resolution
            num_lon = int(lon_extent / lon_resolution) + 1
            num_lat = int(lat_extent / lat_resolution) + 1
            lat = np.linspace(min_lat, max_lat, num=num_lat)
            lon = np.linspace(min_lon, max_lon, num=num_lon)
            # create empty array extending all basins
            if parameter in ["KSat_till", "FieldCap_till", "ThetaS_till"]:
                land_cover_period = data_array.coords['land_cover_period'].values.tolist()
                horizon_till = data_array.coords['horizon_till'].values.tolist()
                full_data_array = xr.DataArray(np.nan, coords=[('lat', lat), ('lon', lon)],
                                              dims=['lat', 'lon'])
                full_data_array = full_data_array.assign_coords(land_cover_period=land_cover_period,
                                                              horizon_till=horizon_till)
            elif parameter in ["L1_FieldCap", "fRoots_temp", "L1_fRoots", "L1_SatSoilMoisture"]:
                land_cover_period_out = data_array.coords['land_cover_period_out'].values.tolist()
                horizon_out = data_array.coords['horizon_out'].values.tolist()
                full_data_array = xr.DataArray(np.nan, coords=[('lat', lat), ('lon', lon)],
                                              dims=['lat', 'lon'])
                full_data_array = full_data_array.assign_coords(land_cover_period_out=land_cover_period_out,
                                                              horizon_out=horizon_out)
            elif parameter == "L1_Max_Canopy_Intercept":
                month_of_year = data_array.coords['month_of_year'].values.tolist()
                full_data_array = xr.DataArray(np.nan, coords=[('lat', lat), ('lon', lon)],
                                              dims=['lat', 'lon'])
                full_data_array = full_data_array.assign_coords(month_of_year=month_of_year)
            full_data_array.attrs['long_name'] = parameter
            for data_array in dataArray_list:
                # Create a temporary DataArray with the same shape as full_data_array but containing the new values
                temp_array = xr.full_like(full_data_array, np.nan)
                temp_array.loc[dict(lat=data_array.lat, lon=data_array.lon)] = data_array.values
                # Use where to mask NaN values in temp_array, then combine_first to fill in non-NaN values
                # This only updates values in full_data_array where temp_array is not NaN
                full_data_array = temp_array.where(~np.isnan(temp_array), full_data_array)
            # prepare and save full data array
            full_data_array = full_data_array.rename({'lat': 'y'})
            full_data_array = full_data_array.rename({'lon': 'x'})
            if parameter == "L1_FieldCap":
                fild_cap_array = full_data_array
            if parameter == "L1_Max_Canopy_Intercept":
                  # Set full_data_array to NaN on all grid points where fild_cap_array is NaN
                full_data_array = full_data_array.where(~np.isnan(fild_cap_array), np.nan)
                # set full_data_array to NA on all grid points where fild_cap_array is NaN
            if version == "exp0-run3":
                file_name = f"mHM_{parameter}"
            elif version == "exp2-run1":
                file_name = f"GenAI_para_{parameter}"
            print_simple(f'{output_path}/plots/{file_name}_horizon{horizon}.png', prepare_tiff_out(full_data_array))
            export_to_geotiff(f'{output_path}/{file_name}.tiff', prepare_tiff_out(full_data_array))









