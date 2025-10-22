#!/usr/bin/env python3

import numpy as np
from os import path
import xarray

def uvt_hgrid(hgrid):
    u = (
        hgrid
        [['x', 'y']]
        .isel(nxp=slice(0, None, 2), nyp=slice(1, None, 2))
        .rename({'y': 'lat', 'x': 'lon', 'nxp': 'xq', 'nyp': 'yh'})
    )

    v = (
        hgrid
        [['x', 'y']]
        .isel(nxp=slice(1, None, 2), nyp=slice(0, None, 2))
        .rename({'y': 'lat', 'x': 'lon', 'nxp': 'xh', 'nyp': 'yq'})
    )

    t = (
        hgrid
        [['x', 'y']]
        .isel(nxp=slice(1, None, 2), nyp=slice(1, None, 2))
        .rename({'y': 'lat', 'x': 'lon', 'nxp': 'xh', 'nyp': 'yh'})
    )

    return u, v, t

def mult(width, total_width=100):
    i = np.arange(1, total_width+1)
    return 1 - np.tanh((2/np.e)*(i-1)/(width-1))


def create_damping(shape, nsponge, s_width, e_width, n_width, rate):
    assert rate < 1
    mult_south = np.zeros(shape)
    mult_east = np.zeros(shape)
    mult_north = np.zeros(shape)
    mult_south[0:nsponge, :] = mult(s_width, total_width=nsponge)[0:nsponge, np.newaxis]
    mult_east[:, -nsponge:] = np.fliplr(mult(e_width, total_width=nsponge)[np.newaxis, 0:nsponge])
    mult_north[-nsponge:, :] = np.flipud(mult(n_width, total_width=nsponge)[0:nsponge, np.newaxis])
    combined = np.maximum(mult_south, mult_east)
    combined = np.maximum(combined, mult_north)
    tanh = combined * rate
    return tanh


def write_damping(hgrid, output_dir, nsponge, width, rate, suffix=None):
    target_u, target_v, target_t = uvt_hgrid(hgrid)
    
    s_dy = hgrid['dy'].isel(ny=0).mean()
    s_width = int(np.round(width * 2 / s_dy))

    e_dx = hgrid['dx'].isel(nx=-1).mean()
    e_width = int(np.round(width * 4 / e_dx))
    
    # For NWA12: limit to part over ocean
    n_dy = hgrid['dy'].isel(ny=-1, nxp=slice(None, None)).mean()
    n_width = int(np.round(width * 2 / n_dy /2))
    
    uv_ds = xarray.Dataset(
        data_vars=dict(
            Idamp_u=(['yh', 'xq'], create_damping(target_u.lon.shape, nsponge, s_width, e_width, n_width, rate)),
            Idamp_v=(['yq', 'xh'], create_damping(target_v.lon.shape, nsponge, s_width, e_width, n_width, rate))
        ),
        coords=dict(
            xh=target_v.xh,
            xq=target_u.xq,
            yh=target_u.yh,
            yq=target_v.yq
        )
    )
    for v in ['Idamp_u', 'Idamp_v']:
        uv_ds[v].attrs['units'] = 's-1'
        uv_ds[v].attrs['cell_methods'] = 'time: point'
    encodings = {v: {'dtype': np.int32} for v in ['yh', 'yq', 'xh', 'xq']}
    encodings.update({v: {'_FillValue': None} for v in ['Idamp_u', 'Idamp_v']})
    fname = 'damping_tgb_uv_b.nc' if suffix is None else f'damping_tanh_uv_{suffix}.nc'

    

    uv_ds.xh.attrs['_FillValue'] = np.nan
    uv_ds.xh.attrs['units'] = "degrees_east"
    uv_ds.xh.attrs['long_name'] = "h point nominal longitude"
    uv_ds.xh.attrs['axis'] = "X"

    uv_ds.yh.attrs['_FillValue'] = np.nan
    uv_ds.yh.attrs['units'] = "degrees_north"
    uv_ds.yh.attrs['long_name'] = "h point nominal longitude"
    uv_ds.yh.attrs['axis'] = "Y"

    uv_ds.xq.attrs['_FillValue'] = np.nan
    uv_ds.xq.attrs['units'] = "degrees_east"
    uv_ds.xq.attrs['long_name'] = "q point nominal longitude"
    uv_ds.xq.attrs['axis'] = "X"

    uv_ds.yq.attrs['_FillValue'] = np.nan
    uv_ds.yq.attrs['units'] = "degrees_north"
    uv_ds.yq.attrs['long_name'] = "q point nominal longitude"
    uv_ds.yq.attrs['axis'] = "Y"

    uv_ds.Idamp_u.attrs['long_name'] ="0 if land, 1 if ocean at tracer points"
    uv_ds.Idamp_u.attrs['cell_methods'] = "time: point"
    uv_ds.Idamp_u.attrs['cell_measures'] = "area: areacello"
    uv_ds.Idamp_u.attrs['units'] = "s-1"

    uv_ds.Idamp_v.attrs['long_name'] ="0 if land, 1 if ocean at tracer points"
    uv_ds.Idamp_v.attrs['cell_methods'] = "time: point"
    uv_ds.Idamp_v.attrs['cell_measures'] = "area: areacello"
    uv_ds.Idamp_v.attrs['units'] = "s-1"

    
    uv_ds.to_netcdf(
        path.join(output_dir, fname), 
        format='NETCDF3_64BIT',
        engine='netcdf4',
        encoding=encodings
    )
    
    t_ds = xarray.Dataset(
        data_vars=dict(
            Idamp=(['yh', 'xh'], create_damping(target_t.lon.shape, nsponge, s_width, e_width, n_width, rate))
        ),
        coords=dict(
            xh=target_t.xh,
            yh=target_t.yh
        )
    )
    t_ds['Idamp'].attrs['units'] = 's-1'
    t_ds['Idamp'].attrs['cell_methods'] = 'time: point'
    encodings = {v: {'dtype': np.int32} for v in ['yh', 'xh']}
    encodings.update({'Idamp': {'_FillValue': None}})
    fname = 'damping_tgb_t_b.nc' if suffix is None else f'damping_tanh_t_{suffix}.nc'

    t_ds.xh.attrs['units'] = "degrees_east"
    t_ds.xh.attrs['long_name'] = "h point nominal longitude"
    t_ds.xh.attrs['axis'] = "X"

    t_ds.yh.attrs['units'] = "degrees_north"
    t_ds.yh.attrs['long_name'] = "h point nominal longitude"
    t_ds.yh.attrs['axis'] = "Y"


    t_ds.Idamp.attrs['long_name'] ="0 if land, 1 if ocean at tracer points"
    t_ds.Idamp.attrs['cell_methods'] = "time: point"
    t_ds.Idamp.attrs['cell_measures'] = "area: areacello"
    t_ds.Idamp.attrs['units'] = "s-1"


    t_ds.to_netcdf(
        path.join(output_dir, fname),
        format='NETCDF3_64BIT',
        engine='netcdf4',
        encoding=encodings,
        unlimited_dims='time'
    )


if __name__ == '__main__':
    import argparse
    from yaml import safe_load
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config')
    args = parser.parse_args()
    with open(args.config, 'r') as file: 
        config = safe_load(file)
    hgrid = xarray.open_dataset(config['filesystem']['ocean_hgrid'])
    output_dir = config['filesystem']['output_dir']
    write_damping(hgrid, output_dir, 250, 20e3, 1 / (7 * 24 * 3600))
