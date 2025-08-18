"""
Script to generate LAI fields for each basin separately so that the format matches the formatting of the original
lai_class of mHM. The LAI is based on MODIS data(MOD15A3H) and we compute seasonal LAI values.
"""
from pathlib import Path

import pandas as pd
import xarray

if __name__ == '__main__':

    basins_small = pd.read_csv('data_dir/required_data/basins/small_basins.csv')['Stat_ID'].tolist()
    basins_large = pd.read_csv('data_dir/required_data/basins/large_basins.csv')['Stat_ID'].tolist()
    basins = [f'sub_{x}' for x in basins_large + basins_small]

    for split in ['training', 'validation']:
        summer_lai_path = Path('data_dir/required_data/MOD15A3H_LAI_31468') / split / 'lai.nc'

        spring_lai = xarray.open_dataset(summer_lai_path)['LAI_spring']
        summer_lai = xarray.open_dataset(summer_lai_path)['LAI_summer']
        autumn_lai = xarray.open_dataset(summer_lai_path)['LAI_autumn']
        winter_lai = xarray.open_dataset(summer_lai_path)['LAI_winter']

        for basin in basins:
            print(f'processing: {basin}')

            outfile_dir = Path(f'data_dir/LAI') / split

            if not outfile_dir.exists():
                outfile_dir.mkdir(parents=True)

            lai_class = xarray.open_dataset(
                Path('data_dir/required_data/static') / f'{basin}' / 'mpr' / 'lai_class.nc')
            static = lai_class['lai_class']

            lai_spring = spring_lai.sel(lat=slice(static.lat.min(), static.lat.max()),
                                        lon=slice(static.lon.min(), static.lon.max()))

            lai_summer = summer_lai.sel(lat=slice(static.lat.min(), static.lat.max()),
                                        lon=slice(static.lon.min(), static.lon.max()))

            lai_autumn = autumn_lai.sel(lat=slice(static.lat.min(), static.lat.max()),
                                        lon=slice(static.lon.min(), static.lon.max()))

            lai_winter = winter_lai.sel(lat=slice(static.lat.min(), static.lat.max()),
                                        lon=slice(static.lon.min(), static.lon.max()))

            if lai_summer.shape == static.sel(month_of_year=1).shape:
                dummy = lai_class.copy()
                for month in [3, 4, 5]:
                    dummy['lai_class'].loc[{'month_of_year': month}] = lai_spring
                for month in [6, 7, 8]:
                    dummy['lai_class'].loc[{'month_of_year': month}] = lai_summer
                for month in [9, 10, 11]:
                    dummy['lai_class'].loc[{'month_of_year': month}] = lai_autumn
                for month in [12, 1, 2]:
                    dummy['lai_class'].loc[{'month_of_year': month}] = lai_winter
                dummy.to_netcdf(outfile_dir / f'lai_{basin}.nc')
