import numpy as np
from os import path
#import warnings
import xarray as xr
import xesmf

# try:
#     import .boundary as bnd
# except ImportError:
import boundary as bnd


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

cobalt_rename = {'st_ocean': 'z', 'geolat_t': 'lat', 'geolon_t': 'lon'}
# cobalt


fpath_cobalt = '/projects/schultz/data/cobalt_global/ocean_cobalt_tracers.1988-2007.ann.nc'
grid_file    = '/projects/schultz/d.sasaki/experiments/v1.1_simulation/tools_and_data/data/source/ocean_hgrid.nc'
ds = xr.open_dataset(fpath_cobalt)
ds = ds.rename(**cobalt_rename)[cobalt_vars]
cobalt_flooded = xr.merge((
    bnd.flood_missing(ds[v], xdim='xt_ocean', ydim='yt_ocean', zdim='z') for v in ds.data_vars
))

hgrid = xr.open_dataset(grid_file)


# Need to load or else xesmf will fail when trying to recognize coordinates.
cobalt_flooded = cobalt_flooded.load()    
cobalt_flooded = cobalt_flooded.assign_coords(lat=ds['lat'], lon=ds['lon'])



# For 4P, create medium properties from large
for v in ['si', 'fe', 'n']:
    cobalt_flooded[f'{v}md'] = cobalt_flooded[f'{v}lg']

# For variable n:p, create p from n
cobalt_flooded['psm'] = cobalt_flooded['nsm'] / 24.0
cobalt_flooded['pmd'] = cobalt_flooded['nmd'] / 20.0
cobalt_flooded['plg'] = cobalt_flooded['nlg'] / 14.0
cobalt_flooded['pdi'] = cobalt_flooded['ndi'] / 40.0

common_kws = dict(write=False)



segments = dict(id = 1, border = 'south')

# Load segments
segments = []
# for seg_config in config.get('segments', []):

    # segment = Segment(seg_config['id'], seg_config['border'], hgrid, output_dir=output_dir)
    # segment = bnd.Segment(1,''birtg, hgrid, output_dir='.')
    # segments.append(segment)

segment = bnd.Segment(1,'north', hgrid, output_dir='.')
segments.append(segment)



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

# --------------------- cobalt ---------------#