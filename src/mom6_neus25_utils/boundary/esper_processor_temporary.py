
# esper specific
#   open esper file (willbe dims: time, z, lat, lon) 
#   rename coordinates and variables
#   adjust concentrations (factor multiplication)
#   set up time properties
#   regrid tracer onto segment
#   correct negative values
#   add coordinates
#   save

import xarray as xr
import glob

from PyESPER.lir import lir
from PyESPER.nn import nn
from PyESPER.mixed import mixed

import matplotlib.pyplot as plt
import numpy as np


def get_data(year=1995, clim_type='season'):
    # paths are hardcoded
    # read T and S from mom6 boundary conditions
    files = sorted(glob.glob(f'data/esper/so*{year}*nc'))
    ds = []
    for f in files: 
        print(f)
        aux = xr.open_dataset(f)
        aux = aux.set_coords(['lon_segment_002', 'lat_segment_002']) 
        ds.append(aux)
    
    ds_salt = xr.concat(ds, dim='time')

    files = sorted(glob.glob(f'data/esper/thetao*{year}*nc'))
    ds = []
    for f in files:
        print(f)
        aux = xr.open_dataset(f)
        aux = aux.set_coords(['lon_segment_002', 'lat_segment_002']) 
        ds.append(aux)
    
    ds_temp = xr.concat(ds, dim='time')

    ds = xr.merge([ds_salt, ds_temp])
    
    dsmonthly = ds.groupby(f'time.{clim_type}').mean(dim='time')
    dsmonthly = dsmonthly.rename(season='time')

    # preparing output dataset based on salinity file
    dsout = dsmonthly[['so_segment_002',
                       'lon_segment_002',
                       'lat_segment_002',
                       'ny_segment_002',
                       'dz_so_segment_002']].copy(deep=True)
    
    dsout = dsout.rename({'so_segment_002': 'no3_segment_002'})
    dsout.no3_segment_002.values[:] = 0

    # reading topography
    dstopo = xr.open_dataset('data/esper/ocean_topog.nc')

    # reading woa
    dswoa = xr.open_dataset('data/esper/bgc_woa_002.nc',decode_times=False)
    z = ds.dz_thetao_segment_002.cumsum(dim='nz_segment_002')

    return dsmonthly, dsout, dstopo, dswoa, z