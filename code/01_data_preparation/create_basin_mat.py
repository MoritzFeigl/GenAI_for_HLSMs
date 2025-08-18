"""
Script to generate for each basin separately so that the format matches the formatting of the original
mHM static features (i.e. slope).
"""
from pathlib import Path

import pandas as pd
import xarray

if __name__ == '__main__':

    basins_small = pd.read_csv('data_dir/required_data/basins/small_basins.csv')['Stat_ID'].tolist()
    basins_large = pd.read_csv('data_dir/required_data/basins/large_basins.csv')['Stat_ID'].tolist()
    basins = [f'sub_{x}' for x in basins_large + basins_small]

    for split in ['training', 'validation']:

        data_dir = Path('data_dir/MAT_MAT_RANGE') / split
        mat = xarray.open_rasterio(data_dir / 'temp_mean_31468.tif').sel(band=1)
        mat = mat.rename({'x': 'lon', 'y': 'lat'})
        mat = mat.drop('band')

        mat_range = xarray.open_rasterio(data_dir / 'temp_range_31468.tif').sel(band=1)
        mat_range = mat_range.rename({'x': 'lon', 'y': 'lat'})
        mat_range = mat_range.drop('band')

        maprec = xarray.open_rasterio(data_dir / 'precip_sum_31468.tif').sel(band=1)
        maprec = maprec.rename({'x': 'lon', 'y': 'lat'})
        maprec = maprec.drop('band')

        for basin in basins:

            print(f'processing: {basin}')

            slope = xarray.open_dataset(Path('data_dir/required_data/static') / f'{basin}' / 'mpr' / 'slope.nc')
            static = slope['slope']

            basin_mat = mat.sel(lat=slice(static.lat.max(), static.lat.min()),
                                lon=slice(static.lon.min(), static.lon.max()))

            basin_mat_range = mat_range.sel(lat=slice(static.lat.max(), static.lat.min()),
                                            lon=slice(static.lon.min(), static.lon.max()))

            basin_map = maprec.sel(lat=slice(static.lat.max(), static.lat.min()),
                                   lon=slice(static.lon.min(), static.lon.max()))

            dummy = slope.copy()
            dummy['slope'].data = basin_mat.data[::-1]
            dummy = dummy.rename({'slope': 'mat'})
            dummy.to_netcdf(data_dir / f'mat_{basin}.nc',
                            encoding={'mat': {'_FillValue': -9999., 'dtype': 'float64'}})

            dummy = slope.copy()
            dummy['slope'].data = basin_mat_range.data[::-1]
            dummy = dummy.rename({'slope': 'mat_range'})
            dummy.to_netcdf(data_dir / f'mat_range_{basin}.nc',
                            encoding={'mat_range': {'_FillValue': -9999., 'dtype': 'float64'}})

            dummy = slope.copy()
            dummy['slope'].data = basin_map.data[::-1]
            dummy = dummy.rename({'slope': 'map'})
            dummy.to_netcdf(data_dir / f'map_{basin}.nc',
                            encoding={'map': {'_FillValue': -9999., 'dtype': 'float64'}})


