import numpy as np
from os import path
#import warnings
import xarray as xr
import xesmf
import boundary as bnd
import yaml
import datetime as dtt
import argparse
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


def read_config(config_file):
    with open(config_file, 'r') as stream:
        config = yaml.safe_load(stream)
    return config

class CobaltBoundary:
    """Load, transform, and export COBALT tracers onto MOM6 boundary segments.

    Args:
        fpath_cobalt (str): Path to cobalt netCDF file.
        grid_file (str): Path to ocean_hgrid.nc.
        output_dir (str): Directory for output netCDF files.
        cache_dir (str): Directory for xesmf weight files.
        segments (list): List of dicts with 'id' (int) and 'border' (str: 'north', 'south', 'east', or 'west').
        cobalt_rename (dict): Mapping of native cobalt dim/coord names to required names.
            Must map some key to 'z', 'lat', and 'lon'.
        flood_missing_rename (dict): Must contain keys 'xdim', 'ydim', 'zdim' pointing
            to native cobalt dimension names.
    """

    vars = [
        'cadet_arag', 'cadet_calc',
        'fed', 'fedi', 'felg', 'fedet', 'fesm',
        'ldon', 'ldop', 'lith', 'lithdet',
        'nbact', 'ndet', 'ndi', 'nlg', 'nsm', 'nh4',
        'pdet',
        'srdon', 'srdop', 'sldon', 'sldop',
        'sidet', 'silg',
        'nsmz', 'nmdz', 'nlgz'
    ]

    def __init__(self, fpath_cobalt, grid_file, output_dir, cache_dir,
                 segments, cobalt_rename, flood_missing_rename, vars=None,time0=None):
        self.fpath_cobalt = fpath_cobalt
        self.grid_file = grid_file
        self.output_dir = output_dir
        self.cache_dir = cache_dir
        self.segments = segments
        self.cobalt_rename = cobalt_rename
        self.flood_missing_rename = flood_missing_rename
        self.ds = None
        self.hgrid = None
        self.vars = vars if vars is not None else CobaltBoundary.vars
        self.time0= time0

    def _validate(self):
        assert 'z' in self.cobalt_rename.values(), \
            "cobalt_rename must map some key to 'z': required by regrid_tracer"
        assert 'lat' in self.cobalt_rename.values() and 'lon' in self.cobalt_rename.values(), \
            "cobalt_rename must map some keys to 'lat' and 'lon': required by assign_coords/xesmf"
        assert all(k in self.flood_missing_rename for k in ['xdim', 'ydim', 'zdim']), \
            "flood_missing_rename must contain keys: 'xdim', 'ydim', 'zdim'"

    def load(self):
        self._validate()
        ds = xr.open_dataset(self.fpath_cobalt)
        ds = ds.rename(**self.cobalt_rename)[self.vars]

        if self.time0 is not None:
            ds['time'] = [self.time0]

        # flood land points; xdim/ydim/zdim match native cobalt dim names
        self.ds = xr.merge((
            bnd.flood_missing(ds[v], **self.flood_missing_rename) for v in ds.data_vars))

        # load required before xesmf can recognize coordinates
        self.ds = self.ds.load()
        self.ds = self.ds.assign_coords(lat=ds['lat'], lon=ds['lon'])
        self.hgrid = xr.open_dataset(self.grid_file)
        return self

    def cobaltv2_to_v3(self):
        """Convert COBALTv2 variables to v3."""
        for v in ['si', 'fe', 'n']:
            self.ds[f'{v}md'] = self.ds[f'{v}lg']

        self.ds['psm'] = self.ds['nsm'] / 24.0
        self.ds['pmd'] = self.ds['nmd'] / 20.0
        self.ds['plg'] = self.ds['nlg'] / 14.0
        self.ds['pdi'] = self.ds['ndi'] / 40.0
        return self

    def export(self):
        segments = [
            bnd.Segment(s['id'], s['border'], self.hgrid,
                        regrid_dir=self.cache_dir,
                        output_dir=self.output_dir)
            for s in self.segments
        ]

        for seg in segments:
            cobalt_seg = xr.merge(
                (seg.regrid_tracer(self.ds[v],
                                   regrid_suffix='cobalt',
                                   flood=False,
                                   periodic=False,
                                   write=False) for v in self.ds.data_vars)
            )
            for v in cobalt_seg.data_vars:
                cobalt_seg[v] = np.clip(cobalt_seg[v], 0.0, None)
            cobalt_seg = seg.add_coords(cobalt_seg)
            seg.to_netcdf(cobalt_seg, 'bgc_cobalt')
        return self


class WOABoundary:
    """Load and export WOA climatology tracers onto MOM6 boundary segments.

    Args:
        fpath_woa (str): Path to WOA netCDF file.
        grid_file (str): Path to ocean_hgrid.nc.
        output_dir (str): Directory for output netCDF files.
        cache_dir (str): Directory for xesmf weight files.
        segments (list): List of dicts with 'id' (int) and 'border' (str: 'north', 'south', 'east', or 'west').
        vars (list, optional): List of variables to load. Defaults to all vars in dataset.
    """

    # dim names match WOA natively; no renaming needed
    flood_kws = dict(xdim='lon', ydim='lat', zdim='z')

    # time attributes required for MOM6 climatological forcing
    time_attrs = {
        'units': 'days since 0001-01-01',
        'calendar': 'noleap',
        'modulo': ' ',
        'cartesian_axis': 'T'
    }

    def __init__(self, fpath_woa, grid_file,
                 output_dir, cache_dir, segments,
                 vars=None, flood_kws=None, time_attrs=None):
        self.fpath_woa = fpath_woa
        self.grid_file = grid_file
        self.output_dir = output_dir
        self.cache_dir = cache_dir
        self.segments = segments
        self.vars = vars
        self.ds = None
        self.hgrid = None
        self.flood_kws = flood_kws if flood_kws is not None \
                         else WOABoundary.flood_kws
        self.time_attrs = time_attrs if time_attrs is not None \
                         else WOABoundary.time_attrs
        print(self.flood_kws)

    def load(self):
        self.ds = xr.open_dataset(self.fpath_woa)
        if self.vars is not None:
            self.ds = self.ds[self.vars]
        # depth -> z required by regrid_tracer
        self.ds = self.ds.rename({'depth': 'z'})

        self.hgrid = xr.open_dataset(self.grid_file)
        return self

    def export(self):
        segments = [
            bnd.Segment(s['id'], s['border'], self.hgrid,
                        regrid_dir=self.cache_dir,
                        output_dir=self.output_dir)
            for s in self.segments
        ]

        for seg in segments:
            woa_seg = xr.merge(
                (seg.regrid_tracer(self.ds[v],
                                   regrid_suffix='woa_bgc',
                                   flood=True,
                                   periodic=False,
                                   write=False,
                                   **self.flood_kws) for v in self.ds.data_vars)
            )
            for v in woa_seg.data_vars:
                woa_seg[v] = np.clip(woa_seg[v], 0.0, None)
            woa_seg = seg.add_coords(woa_seg)
            # woa_seg['time'].attrs.update(self.time_attrs)
            seg.to_netcdf(woa_seg, 'bgc_woa')
        return self


def main():
    import argparse
    import yaml



    parser = argparse.ArgumentParser(description='Process WOA/COBALT boundary data')
    parser.add_argument('--config', type=str, default='config.yaml')
    parser.add_argument('--year', type=int,
                        help='Single year to process')

    args = parser.parse_args()


    with open(args.config) as f:
        config = yaml.safe_load(f)
    config = config['boundary']


    cobalt_rename = {'geolat_t': 'lat', 'geolon_t': 'lon', 'st_ocean': 'z'}
    flood_missing_rename = dict(xdim='xt_ocean', ydim='yt_ocean', zdim='z')
    time0 = dtt.datetime.strptime(str(config['time0']), '%Y-%m-%d')

    (CobaltBoundary(fpath_cobalt=config['cobalt_file'],
                    grid_file=config['grid_file'],
                    output_dir=config['output_dir'],
                    cache_dir=config['cache'],
                    segments=config['segments'],
                    cobalt_rename=cobalt_rename,
                    flood_missing_rename=flood_missing_rename,
                    time0=time0)
                        .load()
                        .cobaltv2_to_v3()
                        .export())
    
    (WOABoundary(
        fpath_woa=config['woa_file'],
        grid_file=config['grid_file'],
        output_dir=config['output_dir'],
        cache_dir=config['cache'],
        segments=config['segments'])
            .load().export())

if __name__ == '__main__':
    main()

    # config = read_config('config.yaml')
    # config = config['boundary']

    # cobalt_rename = {'geolat_t': 'lat', 'geolon_t': 'lon', 'st_ocean': 'z'}
    # flood_missing_rename = dict(xdim='xt_ocean', ydim='yt_ocean', zdim='z')
    # time0 = dtt.datetime.strptime(str(config['time0']), '%Y-%m-%d')

    # (CobaltBoundary(fpath_cobalt=config['cobalt_file'],
    #                 grid_file=config['grid_file'],
    #                 output_dir=config['output_dir'],
    #                 cache_dir=config['cache'],
    #                 segments=config['segments'],
    #                 cobalt_rename=cobalt_rename,
    #                 flood_missing_rename=flood_missing_rename,
    #                 time0=time0)
    #                     .load()
    #                     .cobaltv2_to_v3()
    #                     .export())
    
    # (WOABoundary(
    #     fpath_woa=config['woa_file'],
    #     grid_file=config['grid_file'],
    #     output_dir=config['output_dir'],
    #     cache_dir=config['cache'],
    #     segments=config['segments'])
    #         .load().export())

    # cobalt_flooded, hgrid = load_cobalt(fpath_cobalt, grid_file, cobalt_rename, flood_missing_rename)
    # cobalt_flooded = cobaltv2_to_v3(cobalt_flooded)
    # export_segments(config, cobalt_flooded)



