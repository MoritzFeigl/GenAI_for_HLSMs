import warnings
from pathlib import Path
from typing import Union, List

import numpy as np
import pandas as pd
import xarray
from numba import jit


class SpatialPredictors(object):
    """ Spatial predictors class to load and process spatial predictors for mHM.

    Parameters
    ----------
    data_dir : Path
        Main directory which contains the basin subfolders (each basin folder contains the folders 'mpr' and 'routing'
    predictor : str
        Name of the spatial predictor, which has to be the same as the corresponding netcdf name, e.g. sand for sand.nc
    basin_names : str, List[str],
        Names of the subbasins for which the predictor should be extracted
    horizon : float, List[float] (optional)
        Horizon of the predictor if available (only soil related predictors have horizons)
    """

    def __init__(self,
                 data_dir: Path,
                 predictor: str,
                 basin_names: Union[List[str], str],
                 horizon: Union[List[float], float] = None):
        self.data_dir = data_dir
        self.predictor_name = predictor
        self.horizon = horizon

        if isinstance(basin_names, List):
            self.basin_names = basin_names
        else:
            self.basin_names = [basin_names]

        self.agg_levels = {}
        self.stats = {}
        self.aggregation_function = None

        temp = xarray.open_dataset(data_dir / self.basin_names[0] / 'mpr' / f'{predictor}.nc')
        self.native_resolution_km = temp['lat'].diff(dim='lat').values[0] / 1000

        self.basins = {}
        for basin_name in self.basin_names:
            print(f'Loading basin {basin_name}')
            self.basins[basin_name] = self._load_basin(basin_name=basin_name)

    def aggregate_basins(self, resolution_km: Union[List[Union[int, float]], int, float], aggregation_fun: str):
        """ Resample the original data to a coarser resolution given a defined aggregation function

        Parameters
        ----------
        resolution_km: float,
            pixel size of the resampled data array in km
        aggregation_fun: str,
            either 'mean' or 'median' can be computed

        """

        self.aggregation_function = aggregation_fun

        if not isinstance(resolution_km, List):
            resolution_km = [resolution_km]

        for basin_name, basin in self.basins.items():
            print(f'Aggregating basin {basin_name}')
            self.agg_levels[basin_name] = {}

            for resolution in resolution_km:
                scale_ratio = int(resolution / self.native_resolution_km)

                lats = basin['lat'].data
                lons = basin['lon'].data
                data = basin.data[::-1, :]

                resampled = _resample_loop(data, scale_ratio, lats, lons, aggregation_fun)

                df = self._calculate_stats(resampled)
                df.name = f'{resolution}km'

                self.agg_levels[basin_name][resolution] = resampled
                self.stats[basin_name] = pd.concat([self.stats[basin_name], df.to_frame()], axis=1)

    def _load_basin(self, basin_name: str) -> xarray.Dataset:
        """ Loads original data and calculates statistics.

        Parameters
        ----------
        basin_name: str
            Name of the subbasin

        Returns
        -------
        xarray Dataset with the original data

        """

        xr = xarray.open_dataset(self.data_dir / basin_name / 'mpr' / f'{self.predictor_name}.nc')[self.predictor_name]

        if 'horizon' in xr.dims:
            horizons = [x for x in self.horizon if x in xr['horizon']]
            xr = xr.sel(horizon=horizons).mean(dim='horizon')

        df = self._calculate_stats(xr.data)
        df.name = f'{self.native_resolution_km}km'

        self.stats[basin_name] = df.to_frame()

        return xr

    @staticmethod
    def _calculate_stats(arr: np.ndarray) -> pd.DataFrame:
        """ Calculates a list of statistics for a given array and returns it as a pd.DataFrame

        Parameters
        ----------
        arr

        """

        df = pd.Series(dtype='float64')
        df['mean'] = np.nanmean(arr)
        df['median'] = np.nanmedian(arr)
        df['std'] = np.nanstd(arr)
        df['min'] = np.nanmin(arr)
        df['q1'] = np.nanpercentile(arr, 1)
        df['q5'] = np.nanpercentile(arr, 5)
        df['q95'] = np.nanpercentile(arr, 95)
        df['q99'] = np.nanpercentile(arr, 99)
        df['max'] = np.nanmax(arr)
        df['range'] = df['max'] - df['min']

        return df

    def get_basin_stats(self, basin_name: str = None) -> pd.DataFrame:
        """ Returns the statistics for a given basin

        Parameters
        ----------
        basin_name: str
            Name of the subbasin

        """
        if (basin_name is None) and len(self.basin_names) == 1:
            return self.stats[self.basin_names[0]].transpose()
        elif (basin_name is not None) and (basin_name in self.basin_names):
            return self.stats[basin_name].transpose()
        else:
            warnings.warn(f'Select either of the following basins as basin_name: {",".join(self.basin_names)}',
                          FutureWarning)

    def get_original_data(self, basin_name: str = None) -> np.ndarray:
        """ Returns the original data for a given basin

        Parameters
        ----------
        basin_name: str
            Name of the subbasin

        """
        if (basin_name is None) and len(self.basin_names) == 1:
            return self.basins[self.basin_names[0]].data[::-1, :]
        elif (basin_name is not None) and (basin_name in self.basin_names):
            return self.basins[basin_name].data[::-1, :]
        else:
            warnings.warn(f'Select either of the following basins as basin_name: {",".join(self.basin_names)}',
                          FutureWarning)

    def get_resampled_data(self, resolution_km: Union[int, float], basin_name: str = None) -> np.ndarray:
        """ Returns the resampled data for a given spatial resolution and basin

        Parameters
        ----------
        resolution_km: float,
            pixel size of the resampled data array in km
        basin_name: str
            Name of the subbasin

        """

        if (basin_name is None) and len(self.basin_names) == 1:
            basin_name = self.basin_names[0]

        if (basin_name is not None) and (basin_name in self.basin_names):
            if self.agg_levels.get(basin_name, {}).get(resolution_km, {}) is not None:
                return self.agg_levels[basin_name][resolution_km]
            else:
                warnings.warn(f'Data has not been resampled to {resolution_km}km yet.', FutureWarning)

        else:
            warnings.warn(f'Select either of the following basins as basin_name: {",".join(self.basin_names)}',
                          FutureWarning)


@jit(nopython=True)
def _resample_loop(data, scale_ratio, lats, lons, aggregation_fun) -> np.ndarray:
    """ Runs the loop for aggregating the original data using numba

    Parameters
    ----------
    data: np.ndarray,
        data of the xarray dataarray to be resampled
    scale_ratio: float
        ratio between the original resolution and the resampled resolution
    lats: List
        List of latitudes from the xarray dataarray
    lons: List
        List of longitudes from the xarray dataarray
    aggregation_fun: str,
                either 'mean' or 'median' can be computed

    Returns
    -------
    np.ndarray of the same shape as the original data but with resampled values

    """

    lat_bins = list(range(0, len(lats), scale_ratio))
    lon_bins = list(range(0, len(lons), scale_ratio))

    resampled = np.zeros(data.shape)
    mask = data * resampled + 1

    for idx_lat, lat in enumerate(lat_bins):
        for idx_lon, lon in enumerate(lon_bins):
            lat_end = min(lat + scale_ratio, len(lats))
            lon_end = min(lon + scale_ratio, len(lons))

            if aggregation_fun == 'mean':
                resampled[lat:lat_end, lon:lon_end] = np.nanmean(data[lat:lat_end, lon:lon_end])
            elif aggregation_fun == 'median':
                resampled[lat:lat_end, lon:lon_end] = np.nanmedian(data[lat:lat_end, lon:lon_end])

    return resampled * mask
