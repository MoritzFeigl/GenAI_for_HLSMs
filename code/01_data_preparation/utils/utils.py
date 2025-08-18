from pathlib import Path

import numpy as np
import rasterio
import rasterio.features
import rasterio.io
import rasterio.mask
from rasterio.transform import Affine
from typing import Union

def write_numpy2raster(arr: np.ndarray, out_filepath: Union[str, Path], latitude: np.ndarray, longitude: np.ndarray,
                       crs: rasterio.crs.CRS):
    """ Write xarray to raster file.

    Parameters
    ----------
    arr: numpy array
        Array to write
    out_filepath: Union[str, Path)
    longitude
    latitude
    crs: rasterio.crs.CRS
        coordinate system of the DataArray

    """
    ps = longitude[1] - longitude[0]
    geotransform = (longitude.min(), ps, 0, latitude.max(), 0, -ps)
    transform = Affine.from_gdal(*geotransform)

    if len(arr.shape) > 2:
        count = len(arr)
        height = arr.shape[1]
        width = arr.shape[2]
    else:
        count = 1
        height = arr.shape[0]
        width = arr.shape[1]

    profile = {
        'driver': 'GTiff',
        'height': height,
        'width': width,
        'count': count,
        'dtype': arr.dtype,
        'crs': crs,
        'transform': transform
    }

    with rasterio.open(out_filepath, 'w', **profile) as out_image:
        out_image.write(arr, count)

