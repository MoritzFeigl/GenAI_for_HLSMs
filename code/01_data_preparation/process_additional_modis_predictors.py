"""
Script to process MODIS LAI data (downloaded from GEE) and write it into a netcdf file. The script for downloading the
data from GEE can be found under /data_download/GEE_lai_download.txt of this repository. The data is processed
separately for a training and validation period.
"""
from pathlib import Path

import numpy as np
import xarray


def get_modis_features(feature_path: Path, mask: xarray.DataArray, clip_zero: bool = True, crop: bool = False) -> \
        xarray.DataArray:
    """

    Parameters
    ----------
    feature_path: Path,
        Path to the MODIS image
    mask: DataArray,
        Mask for cropping and adapting the resolution, template file of an original static mHM feature
    clip_zero: (optional) bool
        If True, the data_download will be clipped to zero (no negative values), default is True
    crop: optional bool
        If True, the data_download will be clipped to the given extent, default is False

    Returns
    -------
    DataArray of the MODIS image with adapted resolution and extent

    """
    feature = xarray.open_rasterio(feature_path)
    feature = feature.isel(band=0).drop("band")

    res_original = (mask['lon'][1] - mask['lon'][0])
    res_feature = (feature['x'][1] - feature['x'][0])
    res_ratio = np.abs(res_feature / res_original)

    feature_data = feature.data[::-1, :]
    if clip_zero:
        feature_data[feature_data < 0] = 0

    high_res_feature = feature_data.repeat(res_ratio, axis=0).repeat(res_ratio, axis=1)
    lon = np.arange(feature['x'].data[0] - res_feature/2 + res_original/2, feature['x'].data[-1] + res_feature/2 +
                    res_original/2, res_original)
    lat = np.arange(feature['y'].data[-1] - res_feature/2 + res_original/2, feature['y'].data[0] + res_feature/2 +
                    res_original/2, res_original)

    data = xarray.DataArray(high_res_feature, dims=("lat", "lon"), coords={"lat": lat, "lon": lon})

    if crop:
        data = data.sel(lat=slice(mask['lat'].data.min(), mask['lat'].data.max()),
                        lon=slice(mask['lon'].data.min(), mask['lon'].data.max())).copy()

    return data


if __name__ == '__main__':
    for split in ['training', 'validation']:

        temp_dir = Path(r'data_dir/required_data/static')
        data_dir = Path('data_dir/required_data/MOD15A3H_LAI_31468') / split

        temp_dir = list(temp_dir.glob('*'))[0]
        original_xr = xarray.open_dataset(temp_dir / 'mpr' / 'bd.nc')['bd']

        xr_mask = xarray.where(original_xr.notnull(), 1, np.nan)
        xr_mask = xr_mask.isel(horizon=0).drop("horizon")

        lai_features = [
            data_dir / 'LAI_spring.tif',
            data_dir / 'LAI_summer.tif',
            data_dir / 'LAI_autumn.tif',
            data_dir / 'LAI_winter.tif'
        ]

        lai_dict = {}

        for lai_feature in lai_features:
            print(f'Creating {lai_feature.stem} DataArray')
            lai_dict[lai_feature.stem] = get_modis_features(lai_feature, xr_mask)

        additional_xr = xarray.Dataset(lai_dict)
        additional_xr.attrs['epsg'] = 31468
        additional_xr.attrs['decription'] = f'Spring (March - May, summer (June - August), autumn (September to' \
                                            f'November and winter LAI (Dec, Jan, Feb) MOD15A3H over the {split} period'
        additional_xr.attrs['original_dataset'] = 'MOD15A3H'
        additional_xr.attrs['original_spatial_resolution'] = '500m'
        additional_xr.to_netcdf(data_dir / 'lai.nc')

