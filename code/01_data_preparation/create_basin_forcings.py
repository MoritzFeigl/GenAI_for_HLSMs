"""
Script to generate mHM forcings from data_download provided by UFZ that match the formatting of the original
mHM forcings

"""
from pathlib import Path

import numpy as np
import pandas as pd
import xarray
from scipy import spatial


def get_basin_values(variable: xarray.DataArray, forcing_temp: xarray.DataArray) -> xarray.DataArray:
    """ Extract basin data from whole Germany forcing files

    Parameters
    ----------
    variable: xarray.DataArray,
        The data for the new forcings
    forcing_temp: xarray.DataArray,
        The corresponding old mHM forcings that are used as a template


    Returns
    -------
    xarray.DataArray with data cropped to the extent of the provided mHM forcings

    """
    coords = {'time': variable.coords['time'],
              'lon': forcing_temp.coords['lon'],
              'lat': forcing_temp.coords['lat']
              }

    lats = variable.lat.data.flatten()
    lons = variable.lon.data.flatten()
    lat_lons = np.vstack((lats, lons)).T

    if forcing_temp.lon.max() > np.max(lons):
        return None

    else:
        left_top_idx = spatial.KDTree(lat_lons).query([forcing_temp.lat[0, 0], forcing_temp.lon[0, 0]])[1]
        left_top = np.unravel_index(left_top_idx, (variable.shape[1], variable.shape[2]))

        left_bottom_idx = spatial.KDTree(lat_lons).query([forcing_temp.lat[-1, 0], forcing_temp.lon[-1, 0]])[1]
        left_bottom = np.unravel_index(left_bottom_idx, (variable.shape[1], variable.shape[2]))

        right_top_idx = spatial.KDTree(lat_lons).query([forcing_temp.lat[0, -1], forcing_temp.lon[0, -1]])[1]
        right_top = np.unravel_index(right_top_idx, (variable.shape[1], variable.shape[2]))

        left_min = min(left_top[0], left_bottom[0])
        left_max = max(left_top[0], left_bottom[0])

        basin_data = variable.data[:, left_min: left_max + 1, left_top[1]:right_top[1] + 1]
        basin_xr = xarray.DataArray(basin_data[:, ::-1, :], coords=coords, dims=forcing_temp.dims)

    return basin_xr


def save_netdcf(data: xarray.DataArray, forcings_dir: Path, basin_name: str, variable_name: str):
    """ Writes data as xarray Dataset into a netcdf file

    Parameters
    ----------
    data: xarray.DataArray
        data_ that should be saved
    forcings_dir: Path
        path for saving the file
    basin_name: str
        name of the basin to save (used to generate file name)
    variable_name: str
        name of the variable to save (used to generate file name)

    """
    data = xarray.Dataset({variable_name: data}, coords=data.coords)
    data['time'].encoding['units'] = 'days since 1949-12-31'
    outpath = forcings_dir / f'{basin_name}'
    if not outpath.exists():
        outpath.mkdir(parents=True)
    data.to_netcdf(outpath / f'{variable_name}.nc')


if __name__ == '__main__':

    new_forcings_dir = Path('data_dir/required_data/forcings')
    old_forcings_dir = Path('data_dir/required_data/forcings_UFZ')
    output_dir = Path('data_dir/forcings')

    if not new_forcings_dir.exists():
        output_dir.mkdir(parents=True)

    tavg = xarray.open_dataarray(new_forcings_dir / 'tavg.nc')
    pet = xarray.open_dataarray(new_forcings_dir / 'pet.nc')
    pre = xarray.open_dataarray(new_forcings_dir / 'pre.nc')

    pre = pre.assign_coords(lon=tavg.lon)
    pre = pre.assign_coords(lat=tavg.lat)

    basins_small = pd.read_csv('data_dir/required_data/basins/small_basins.csv')['Stat_ID'].tolist()
    basins_large = pd.read_csv('data_dir/required_data/basins/large_basins.csv')['Stat_ID'].tolist()
    basins = [f'sub_{x}' for x in basins_large + basins_small]

    bad_basins = []
    for basin in basins:
        print(f'processing: {basin}')
        old_forcing = xarray.open_dataarray(old_forcings_dir / f'{basin}' / 'tavg.nc')

        print(f'processing tavg')
        tavg_basin = get_basin_values(variable=tavg, forcing_temp=old_forcing)
        if tavg_basin is not None:
            save_netdcf(data=tavg_basin, forcings_dir=output_dir, basin_name=basin, variable_name='tavg')
        else:
            bad_basins.append(basin)
            continue

        print(f'processing pet')
        pet_basin = get_basin_values(variable=pet, forcing_temp=old_forcing)
        save_netdcf(data=pet_basin, forcings_dir=output_dir, basin_name=basin, variable_name='pet')

        print(f'processing pre')
        pre_basin = get_basin_values(variable=pre, forcing_temp=old_forcing)
        save_netdcf(data=pre_basin, forcings_dir=output_dir, basin_name=basin, variable_name='pre')

