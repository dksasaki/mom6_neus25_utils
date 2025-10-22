#!/usr/bin/env python3
"""
How to use
./write_nudging_data.py --config_file config.yaml
"""

import numpy as np
import os
import pandas as pd
import xarray as xr
import xesmf
import glob
import argparse
from pathlib import Path
from yaml import safe_load


# VARIABLES = ['thetao', 'so']

VARIABLES = {
    'tracer': {'varbs': ['so', 'thetao'], # tracer (yh,xh)
               'grid':  {'geolat': 'lat', 'geolon': 'lon'},
               'lon': 'xh',
               'lat': 'yh'}, 
    'uo'    : {'varbs': ['uo'], # u      (yh,xq)
               'grid': {'geolat': 'lat', 'geolon': 'lon'},
               'lon': 'xh',
               'lat': 'yh'},
    'vo'    : {'varbs': ['vo'], # v      (yq,xh)
               'grid': {'geolat': 'lat', 'geolon': 'lon'},
               'lon': 'xh',
               'lat': 'yh'},
}


def add_bounds(ds, lonstr, latstr):
    # Add data points at end of month, since time_bnds aren't used
    # All points extend to 23:59:59 at end of month, except
    # for the end of the year which is padded to 00:00:00 the next Jan 1.
    # normalize=True rolls down to midnight
    mstart = [d - pd.offsets.MonthBegin(normalize=True) if d.day > 1 else d for d in ds['time'].to_pandas()]
    mend = [d + pd.DateOffset(months=1) if d.month == 12 else d + pd.DateOffset(months=1) - pd.Timedelta(seconds=1) for d in mstart]
    starts = ds.copy()
    starts['time'] = mstart
    ends = ds.copy()
    ends['time'] = mend
    bounded = xr.concat((starts, ends), dim='time').sortby('time')
    # Ensure that order is correct so that time can be unlimited dim
    bounded = bounded.transpose('time', 'depth', latstr, lonstr)
    return bounded


def reuse_regrid(*args, **kwargs):
    filename = kwargs.pop('filename', None)
    reuse_weights = kwargs.pop('reuse_weights', False)

    if reuse_weights:
        if filename.is_file():
            return xesmf.Regridder(*args, reuse_weights=True, filename=filename, **kwargs)
        else:
            regrid = xesmf.Regridder(*args, **kwargs)
            regrid.to_netcdf(filename)
            return regrid
    else:
        regrid = xesmf.Regridder(*args, **kwargs)
        return regrid




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config')
    args = parser.parse_args()
    
    with open(args.config, 'r') as file: 
        config = safe_load(file)
    nudging_data = config['filesystem']['monthly_data_nudging']
    # regridder = None  ## TODO: think how to save regrid file when using more than one gtype
    output_dir = Path(config['filesystem']['output_dir']) / 'nudging'
    output_dir.mkdir(exist_ok=True)
    tmpdir = Path(config['filesystem']['tmp'])
    
    for year in range(config['forecasts']['first_year'], config['forecasts']['last_year']+1):
        print(f'{year}')
        print(nudging_data.format(year=year))
        glorys = (
            xr.open_mfdataset(glob.glob(nudging_data.format(year=year)))
            .rename({'latitude': 'lat', 'longitude': 'lon'})
        )
        print('  Filling')
        filled = glorys.ffill('depth', limit=None)
        filled = filled.compute()

        static = xr.open_dataset(config['filesystem']['ocean_static'])

        dataset_list = []
        for v in VARIABLES:

            print(f'  Working on {v}')
            gtype = VARIABLES[v]['grid']   # determine if this is a tracer, or velocity grid
            varbs = VARIABLES[v]['varbs']  # select variables associated with gtype
            lonstr = VARIABLES[v]['lon']   # longitude string
            latstr = VARIABLES[v]['lat']   # latitude string
            
            target_grid = static[gtype.keys()].rename(gtype)

            # if regridder is None:
            print('  Setting up regridder')
            regridder = reuse_regrid(
                glorys[varbs], target_grid, 
                filename= tmpdir / 'regrid_nudging.nc', 
                method='nearest_s2d', 
                reuse_weights=False,
                periodic=False
            )

            print('  Interpolating')

            interped = (
                regridder(filled[varbs])
                .drop(['lon', 'lat'], errors='ignore')
                .compute()
            ) 
            print('  Setting time bounds and coordinates')

            bounded = add_bounds(interped, lonstr, latstr)
            # Add coordinate information
            bounded[lonstr] = ((lonstr, ), target_grid[lonstr].data)
            bounded[latstr] = ((latstr, ), target_grid[latstr].data)
            all_vars = list(bounded.data_vars.keys()) + list(bounded.coords.keys())
            encodings = {v: {'_FillValue': None} for v in all_vars}
            encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian', 'units': 'days since 1993-01-01'})
            bounded['depth'].attrs = {
                'cartesian_axis': 'Z',
                'positive': 'down'
            }
            bounded['time'].attrs['cartesian_axis'] = 'T'
            bounded[lonstr].attrs = {'cartesian_axis': 'X'}
            bounded[latstr].attrs = {'cartesian_axis': 'Y'}
            bounded.compute()
            dataset_list.append(bounded.copy())
            del(bounded)
        print('  Writing')
        dsout = xr.merge(dataset_list)

        for v in ['thetao', 'so', 'uo', 'vo']:
            dsout[v] = dsout[v].bfill(dim='xh')
            dsout[v] = dsout[v].bfill(dim='yh')
        
        dsout.to_netcdf(
            output_dir / f'nudging_monthly_{year}.nc',
            format='NETCDF3_64BIT',
            engine='netcdf4',
            encoding=encodings,
            unlimited_dims='time'
        )
        glorys.close() 

        print(f'{output_dir}/nudging_monthly_{year}.nc saved.')

if __name__ == '__main__':
    main()