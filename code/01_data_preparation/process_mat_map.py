"""
Script to generate fields of mean annual temperature (mat), mean annual temperature range (mat_range) and mean
annual precipitation (map) from data downloaded from the DWD CDC database. The link to download the data is:
https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/. The data can be used free of charge with
correctly cited (see https://www.dwd.de/DE/service/copyright/copyright_artikel.html?nn=16102&lsbId=625220).
"""
import datetime
from pathlib import Path
from typing import List

import gdal
import numpy as np
import xarray

from rasterio.crs import CRS

from utils.utils import write_numpy2raster


def get_statistics(data_dir: Path, years: List) -> (List, List, List):
    """ Calculates the mean, range and sum of either temperature or precipitation per year (first combining the single
    month and calculating the statistics per year)

    Parameters
    ----------
    data_dir: Path
        Path to the data
    years: List
        List of years to include in the analysis

    Returns
    -------
    Lists with mean, range and sum per year for the given variable
    """
    mean_list = []
    range_list = []
    sum_list = []

    months = list(range(1, 13))

    for year in years:
        month_list = []

        for month in months:
            dummy_date = datetime.datetime(2000, month, 1)
            file = list((data_dir / f'{str(month).zfill(2)}_{dummy_date.strftime("%b")}').glob(
                f'*_{year}{str(month).zfill(2)}.asc'))[0]
            f = open(file, 'r')
            data = np.genfromtxt(f, skip_header=6)
            data[data == nan_value] = np.nan

            if 'temperature' in str(data_dir):
                data = data / 10

            month_list.append(data)

        data = np.stack(month_list)
        data_mean = data.mean(axis=0)
        data_min = data.min(axis=0)
        data_max = data.max(axis=0)
        data_range = data_max - data_min
        data_sum = data.sum(axis=0)

        mean_list.append(data_mean)
        range_list.append(data_range)
        sum_list.append(data_sum)

    return mean_list, range_list, sum_list


if __name__ == '__main__':

    temp_dir = Path('data_dir/required_data/DWD_CDC_data/air_temperature_mean')
    precip_dir = Path('data_dir/required_data/DWD_CDC_data/precipitation')
    output_cellsize = 100

    for split in ['training', 'validation']:
        output_dir = Path('data_dir/MAT_MAT_RANGE') / split

        if not output_dir.exists():
            output_dir.mkdir(parents=True)

        f = open(temp_dir / '01_Jan' / f'grids_germany_monthly_air_temp_mean_201401.asc', 'r')
        lines = f.readlines()
        cellsize = int(lines[4][9:])
        xllcorner = float(lines[2][10:]) + cellsize / 2
        yllcorner = float(lines[3][10:]) + cellsize / 2
        nan_value = float(lines[5][13:])

        if split == 'training':
            years_list = list(range(2014, 2020))
        elif split == 'validation':
            years_list = list(range(2000, 2014))

        temp_mean, temp_range, temp_sum = get_statistics(data_dir=temp_dir, years=years_list)
        precip_mean, precip_range, precip_sum = get_statistics(data_dir=precip_dir, years=years_list)

        temp_mean = np.stack(temp_mean).mean(axis=0)
        temp_range = np.stack(temp_range).mean(axis=0)
        precip_sum = np.stack(precip_sum).mean(axis=0)

        lat = np.arange(yllcorner, yllcorner + cellsize * temp_mean.shape[0], cellsize)
        lon = np.arange(xllcorner, xllcorner + cellsize * temp_mean.shape[1], cellsize)

        crs = CRS.from_epsg('31467')
        write_numpy2raster(temp_mean, output_dir / 'temp_mean.tif', lat, lon, crs)
        write_numpy2raster(temp_range, output_dir / 'temp_range.tif', lat, lon, crs)
        write_numpy2raster(precip_sum, output_dir / 'precip_sum.tif', lat, lon, crs)

        # Create dummy to get extent right
        input_raster = gdal.Open(str(output_dir / 'temp_mean.tif'))
        output_raster = str(output_dir / 'extent_dummy.tif')
        warp = gdal.Warp(output_raster,
                         input_raster,
                         dstSRS='EPSG:31468',
                         format='GTiff',
                         dstNodata=np.nan,
                         xRes=output_cellsize,
                         yRes=output_cellsize)

        warp = None

        extent = xarray.open_rasterio(output_dir / 'extent_dummy.tif').sel(band=1)
        multiplier = output_cellsize
        xmin = round(float(extent.x.min()) / multiplier) * multiplier - multiplier
        xmax = round(float(extent.x.max()) / multiplier) * multiplier + multiplier
        ymin = round(float(extent.y.min()) / multiplier) * multiplier - multiplier
        ymax = round(float(extent.y.max()) / multiplier) * multiplier + multiplier

        # Create mean annual temperature map
        input_raster = gdal.Open(str(output_dir / 'temp_mean.tif'))
        output_raster = str(output_dir / 'temp_mean_31468.tif')
        warp = gdal.Warp(output_raster,
                         input_raster,
                         dstSRS='EPSG:31468',
                         format='GTiff',
                         dstNodata=np.nan,
                         xRes=output_cellsize,
                         yRes=output_cellsize,
                         outputBounds=(xmin, ymin, xmax, ymax))

        warp = None

        # Create mean annual temperature range map
        input_raster = gdal.Open(str(output_dir / 'temp_range.tif'))
        output_raster = str(output_dir / 'temp_range_31468.tif')
        warp = gdal.Warp(output_raster,
                         input_raster,
                         dstSRS='EPSG:31468',
                         format='GTiff',
                         dstNodata=np.nan,
                         xRes=output_cellsize,
                         yRes=output_cellsize,
                         outputBounds=(xmin, ymin, xmax, ymax))
        warp = None

        # Create mean annual precipitation map
        input_raster = gdal.Open(str(output_dir / 'precip_sum.tif'))
        output_raster = str(output_dir / 'precip_sum_31468.tif')
        warp = gdal.Warp(output_raster,
                         input_raster,
                         dstSRS='EPSG:31468',
                         format='GTiff',
                         dstNodata=np.nan,
                         xRes=output_cellsize,
                         yRes=output_cellsize,
                         outputBounds=(xmin, ymin, xmax, ymax))
        warp = None
