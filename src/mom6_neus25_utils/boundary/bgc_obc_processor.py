import numpy as np
from os import path
#import warnings
import xarray as xr
import xesmf
import boundary as bnd
import yaml


# esper specific
#   open esper file (willbe dims: time, z, lat, lon) 
#   rename coordinates and variables
#   adjust concentrations (factor multiplication)
#   set up time properties
#   regrid tracer onto segment
#   correct negative values
#   add coordinates
#   save

# cobalt specific
#   read cobalt file
#   rename specific variables
#   select variables
#   set up time properties
#   flood cobalt
#   create medium phytoplankton if needed
#   adjust phosphate values for phytoplankton
#   save

# woa specific
#   open woa
#   regrid tracer onto segment
#   set up time properties
#   clean up negative values
#   save

# argument parser


def read_config(config_file):
    with open(config_file, 'r') as stream:
        config = yaml.safe_load(stream)
    return config
import boundary as bnd


def load_cobalt(fpath_cobalt, grid_file, cobalt_rename, flood_missing_rename):

    # renamed values required by regrid_tracer ('z') and assign_coords/xesmf ('lat', 'lon')
    assert 'z' in cobalt_rename.values(), \
        "cobalt_rename must map some key to 'z': required by regrid_tracer"
    assert 'lat' in cobalt_rename.values() and 'lon' in cobalt_rename.values(), \
        "cobalt_rename must map some keys to 'lat' and 'lon': required by assign_coords/xesmf"

    # xdim/ydim/zdim are passed to flood_kara to locate dimensions by name
    assert all(k in flood_missing_rename for k in ['xdim', 'ydim', 'zdim']), \
        "flood_missing_rename must contain keys: 'xdim', 'ydim', 'zdim'"
    
    ds = xr.open_dataset(fpath_cobalt)
    # rename only geolat_t/geolon_t -> lat/lon for assign_coords/xesmf
    ds = ds.rename(**cobalt_rename)[cobalt_vars]

    # flood land points; xdim/ydim/zdim match native cobalt dim names
    # xdim/ydim/zdim are keyword arguments in a function called within
    # flood_missing
    cobalt_flooded = xr.merge((
        bnd.flood_missing(ds[v], **flood_missing_rename) for v in ds.data_vars))
    
    # cobalt_flooded = cobalt_flooded.squeeze()
    hgrid = xr.open_dataset(grid_file)

    # Need to load or else xesmf will fail when trying to recognize coordinates.
    cobalt_flooded = cobalt_flooded.load()    
    cobalt_flooded = cobalt_flooded.assign_coords(lat=ds['lat'], lon=ds['lon'])
    return cobalt_flooded, hgrid


def cobaltv2_to_v3(cobalt_flooded):
    # For 4P, create medium properties from large
    for v in ['si', 'fe', 'n']:
        cobalt_flooded[f'{v}md'] = cobalt_flooded[f'{v}lg']

    # For variable n:p, create p from n
    cobalt_flooded['psm'] = cobalt_flooded['nsm'] / 24.0
    cobalt_flooded['pmd'] = cobalt_flooded['nmd'] / 20.0
    cobalt_flooded['plg'] = cobalt_flooded['nlg'] / 14.0
    cobalt_flooded['pdi'] = cobalt_flooded['ndi'] / 40.0
    return cobalt_flooded


def export_segments(config, cobalt_flooded):
    """Regrid cobalt tracers onto MOM6 boundary segments and write to file.

    Args:
        config (dict): Must contain:
            - 'segments': list of dicts with 'id' (int) and 'border' (str: 'north', 'south', 'east', or 'west')
            - 'cache': directory for xesmf weight files
            - 'output_dir': directory for output netCDF files
        cobalt_flooded (xarray.Dataset): Flooded cobalt tracer dataset.
    """
    common_kws = dict(write=False)

    def _set_segments(config):
        segments = []
        for seg_config in config.get('segments', []):
            segment = bnd.Segment(seg_config['id'],
                                seg_config['border'],
                                hgrid,
                                regrid_dir=config['cache'],
                                output_dir=config['output_dir'])
            segments.append(segment)
        return segments
    
    def _save_segments(segments, cobalt_flooded):
        for seg in segments:
            cobalt_seg = xr.merge(
                (seg.regrid_tracer(cobalt_flooded[v],
                                regrid_suffix='cobalt',
                                flood=False,
                                periodic=False,
                                **common_kws) for v in cobalt_flooded)
            )
            # Make sure no negative values were produced, just in case.
            for v in cobalt_seg.data_vars:
                cobalt_seg[v] = np.clip(cobalt_seg[v], 0.0, None)
            cobalt_seg = seg.add_coords(cobalt_seg)
            seg.to_netcdf(cobalt_seg, 'bgc_cobalt')
    
    segments = _set_segments(config)
    _save_segments(segments, cobalt_flooded)

# --------------------- cobalt ---------------#
cobalt_vars = [
    # 'alk',
    'cadet_arag',
    'cadet_calc',
    # 'dic',
    # 'dic14',
    # 'do14',
    # 'do14c', 
    # 'di14c',
    'fed', 
    'fedi',
    'felg',
    'fedet',
    'fesm',
    'ldon',
    'ldop', 
    'lith', 
    'lithdet', 
    'nbact', 
    'ndet', 
    'ndi', 
    'nlg', 
    'nsm', 
    # 'nh3', 
    'nh4', 
    # 'no3', 
    # 'o2', 
    'pdet', 
    # 'po4', 
    'srdon', 
    'srdop', 
    'sldon', 
    'sldop', 
    'sidet', 
    'silg', 
    # 'sio4', 
    'nsmz', 
    'nmdz', 
    'nlgz'
]

config = read_config('config.yaml')

cobalt_rename = {'geolat_t': 'lat', 'geolon_t': 'lon', 'st_ocean': 'z'}
flood_missing_rename = dict(xdim='xt_ocean', ydim='yt_ocean', zdim='z')
fpath_cobalt = '/projects/schultz/data/cobalt_global/ocean_cobalt_tracers.1988-2007.ann.nc'
grid_file    = '/projects/schultz/d.sasaki/experiments/v1.1_simulation/tools_and_data/data/source/ocean_hgrid.nc'


cobalt_flooded, hgrid = load_cobalt(fpath_cobalt, grid_file, cobalt_rename, flood_missing_rename)
cobalt_flooded = cobaltv2_to_v3(cobalt_flooded)
export_segments(config, cobalt_flooded)



