from pathlib import Path

import numpy as np
import pandas as pd
import pyproj
import xarray


def transform_forcing(variable, old_forcings_dir):
    old_forcing = xarray.open_dataarray(old_forcings_dir / 'sub_6340200_original' / f'{variable}.nc')

    wgs84 = pyproj.CRS("EPSG:4326")
    gk = pyproj.CRS("EPSG:31468")

    east = old_forcing.easting.data
    north = old_forcing.northing.data

    mesh = np.meshgrid(north, east, indexing='ij')

    north = mesh[0].flatten()
    east = mesh[1].flatten()

    lat, lon = pyproj.transform(gk, wgs84, north, east)
    lat = lat.reshape(old_forcing.lat.shape)
    lon = lon.reshape(old_forcing.lon.shape)

    old_forcing.lat.data = lat
    old_forcing.lon.data = lon

    return old_forcing


forcings_dir = Path('/data/CFG/UFZ_forcings')

tavg = transform_forcing('tavg', forcings_dir)
tavg.to_netcdf(forcings_dir / 'sub_6340200' / 'tavg.nc')

tavg = transform_forcing('pet', forcings_dir)
tavg.to_netcdf(forcings_dir / 'sub_6340200' / 'pet.nc')

tavg = transform_forcing('pre', forcings_dir)
tavg.to_netcdf(forcings_dir / 'sub_6340200' / 'pre.nc')
