import numpy as np
from os import path
from glob import glob
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


# the variables CobaltBoundary.cobaltv2_to_v3 reads from its own dataset
V2_TO_V3_PARENTS = ['nlg', 'silg', 'felg', 'nsm', 'ndi']


def require_v2_to_v3_parents(annual_vars):
    """Fail early if the v2 -> v3 parents were split across the two sources.

    CobaltBoundary.cobaltv2_to_v3 derives from all five parents at once, so moving
    any of them to the monthly source leaves the annual run unable to convert. This
    says so before load(), rather than surfacing as a bare KeyError once the
    flooding has already been paid for.

    Only relevant when the conversion is actually run; comment out the call in
    main() alongside the cobaltv2_to_v3 line it guards.
    """
    moved = [p for p in V2_TO_V3_PARENTS if p not in annual_vars]
    if moved:
        raise ValueError(
            f'{moved} appear in cobalt_monthly_vars, but '
            'CobaltBoundary.cobaltv2_to_v3 needs all of '
            f'{V2_TO_V3_PARENTS} in the annual source. Either leave them annual, or '
            'skip the conversion by commenting out the cobaltv2_to_v3 line in main().')

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
        cobalt_renamed_dims (dict): Must contain keys 'xdim', 'ydim', 'zdim' pointing
            to native cobalt dimension names.
        chunks (dict, optional): Passed to xr.open_dataset, making the source dask
            backed. Keyed by *native* dimension names, since chunking happens before
            cobalt_rename is applied. Defaults to None, which reads without dask.
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
                 segments, cobalt_rename, cobalt_renamed_dims, vars=None, time0=None,
                 chunks=None):
        self.chunks = chunks
        self.fpath_cobalt = fpath_cobalt
        self.grid_file = grid_file
        self.output_dir = output_dir
        self.cache_dir = cache_dir
        self.segments = segments
        self.cobalt_rename = cobalt_rename
        self.cobalt_renamed_dims = cobalt_renamed_dims
        self.ds = None
        self.hgrid = None
        self.vars = vars if vars is not None else CobaltBoundary.vars
        self.time0= time0

    def _validate(self):
        assert 'z' in self.cobalt_rename.values(), \
            "cobalt_rename must map some key to 'z': required by regrid_tracer"
        assert 'lat' in self.cobalt_rename.values() and 'lon' in self.cobalt_rename.values(), \
            "cobalt_rename must map some keys to 'lat' and 'lon': required by assign_coords/xesmf"
        assert all(k in self.cobalt_renamed_dims for k in ['xdim', 'ydim', 'zdim']), \
            "cobalt_renamed_dims must contain keys: 'xdim', 'ydim', 'zdim'"

    def load(self):
        self._validate()
        ds = xr.open_dataset(self.fpath_cobalt, chunks=self.chunks)
        ds = ds.rename(**self.cobalt_rename)[self.vars]

        if self.time0 is not None:
            ds['time'] = [self.time0]

        # flood land points; xdim/ydim/zdim match native cobalt dim names
        self.ds = xr.merge((
            bnd.flood_missing(ds[v], **self.cobalt_renamed_dims) for v in ds.data_vars))

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


class CobaltBoundaryMonthly:
    """Load, transform, and export a monthly COBALT climatology onto MOM6 segments.

    Companion to CobaltBoundary, for the tracers that are available as a 12-month
    climatology. It handles that source and nothing else: the annual tracers stay
    with CobaltBoundary, and the two write separate files. Keeping them apart means
    neither source constrains the other's grid, so each keeps its own horizontal
    resolution and its own vertical levels.

    The months are put on a 12-step MOM6 modulo time axis, so the forcing cycles
    without reference to a calendar year.

    Note that CobaltBoundary.cobaltv2_to_v3 needs nlg, silg, felg, nsm and ndi in its
    own dataset. If any of those is moved to the monthly source, the annual run can
    no longer do the conversion.

    Args:
        fpath_cobalt_monthly (str or list): Monthly climatology source. Either a single
            path, a glob pattern, or a list of paths (one per variable or one per
            month; combined on their coordinates).
        grid_file (str): Path to ocean_hgrid.nc.
        output_dir (str): Directory for output netCDF files.
        cache_dir (str): Directory for xesmf weight files.
        segments (list): List of dicts with 'id' (int) and 'border' (str: 'north', 'south', 'east', or 'west').
        vars (list): The tracers held in the monthly source. Required; no default.
        cobalt_rename (dict): Mapping of native cobalt dim/coord names to required names.
            Must map some key to 'z', 'lat', and 'lon', and the source must end up with
            a 'time' dimension carrying a coordinate.
        cobalt_renamed_dims (dict): Must contain keys 'xdim', 'ydim', 'zdim' pointing
            to native cobalt dimension names.
        time_attrs (dict, optional): Time attributes written to the output.
            Defaults to CobaltBoundaryMonthly.time_attrs.
        stream (bool, optional): If True, defer flooding to export() and flood one
            variable at a time, discarding each global 3D field as soon as it has been
            regridded onto every segment. Peak memory becomes one variable instead of
            all of them, at the cost of self.ds holding lazy, un-flooded data after
            load(). Recommended for monthly data, which is 12x the annual size.
            Defaults to False, which mirrors CobaltBoundary.
        chunks (dict, optional): Passed to xr.open_dataset or xr.open_mfdataset,
            making the source dask backed. Keyed by *native* dimension names, since
            chunking happens before cobalt_rename is applied. Note that flooding and
            regridding both read across the horizontal, so splitting there costs
            more dask tasks than splitting over time. Defaults to None, which reads
            a single file without dask and a multi-file set one chunk per file.
    """

    # written alongside CobaltBoundary's 'bgc_cobalt' rather than over it
    output_name = 'bgc_cobalt_monthly'

    # time attributes required for MOM6 climatological forcing
    time_attrs = {
        'units': 'days since 0001-01-01',
        'calendar': 'noleap',
        'modulo': ' ',
        'cartesian_axis': 'T'
    }

    # Mid-month days of a 365-day year, consistent with time_attrs['units'].
    # Kept as plain floats rather than dates on purpose: to_netcdf writes
    # time_attrs['units'] as an attribute, and xarray refuses to do that on top of
    # an already-decoded time axis.
    clim_time = np.array([15.5, 45.0, 74.5, 105.0, 135.5, 166.0,
                          196.5, 227.5, 258.0, 288.5, 319.0, 349.5])

    # COBALTv2 -> v3 derived variable: (parent variable, divisor).
    # nmd/simd/femd are copies of the lg variables, so 'pmd = nmd / 20' is flattened
    # to 'nlg / 20'. That leaves every derived variable depending on a single parent,
    # which is what lets streaming mode derive them from the flooded parent in place,
    # and what lets this class build only the ones whose parent it happens to hold.
    v2_to_v3 = {
        'nmd': ('nlg', 1.0),
        'simd': ('silg', 1.0),
        'femd': ('felg', 1.0),
        'psm': ('nsm', 24.0),
        'pmd': ('nlg', 20.0),
        'plg': ('nlg', 14.0),
        'pdi': ('ndi', 40.0)
    }

    def __init__(self, fpath_cobalt_monthly, grid_file, output_dir, cache_dir,
                 segments, vars, cobalt_rename, cobalt_renamed_dims,
                 time_attrs=None, stream=False, chunks=None):
        self.chunks = chunks
        self.fpath_cobalt_monthly = fpath_cobalt_monthly
        self.grid_file = grid_file
        self.output_dir = output_dir
        self.cache_dir = cache_dir
        self.segments = segments
        self.vars = list(vars)
        self.cobalt_rename = cobalt_rename
        self.cobalt_renamed_dims = cobalt_renamed_dims
        self.time_attrs = time_attrs if time_attrs is not None \
                          else CobaltBoundaryMonthly.time_attrs
        self.stream = stream
        self.ds = None
        self.hgrid = None
        # derived variables deferred to export(); only populated in streaming mode
        self._derived = {}

    def _validate(self):
        assert 'z' in self.cobalt_rename.values(), \
            "cobalt_rename must map some key to 'z': required by regrid_tracer"
        assert 'lat' in self.cobalt_rename.values() and 'lon' in self.cobalt_rename.values(), \
            "cobalt_rename must map some keys to 'lat' and 'lon': required by assign_coords/xesmf"
        assert all(k in self.cobalt_renamed_dims for k in ['xdim', 'ydim', 'zdim']), \
            "cobalt_renamed_dims must contain keys: 'xdim', 'ydim', 'zdim'"
        if not self.vars:
            raise ValueError('vars must name at least one tracer')

    @staticmethod
    def _open(fpath, chunks=None):
        """Open a single path, a glob pattern, or a list of paths as one dataset.

        chunks is keyed by *native* dimension names, because chunking happens here,
        before cobalt_rename is applied. None reads a single file without dask, and
        leaves a multi-file read chunked one chunk per file, which is the default.
        """
        if isinstance(fpath, (list, tuple)):
            paths = list(fpath)
        elif any(c in fpath for c in '*?['):
            paths = sorted(glob(fpath))
            if not paths:
                raise FileNotFoundError(f'no files matched the pattern {fpath!r}')
        else:
            return xr.open_dataset(fpath, chunks=chunks)
        if len(paths) == 1:
            return xr.open_dataset(paths[0], chunks=chunks)
        # Handles both one-file-per-variable and one-file-per-month layouts.
        # data_vars/coords='minimal' takes the time-invariant variables from the first
        # file instead of concatenating them along time, which is what would otherwise
        # turn geolat_t/geolon_t into a 3D lat/lon that xesmf cannot build a grid from.
        return xr.open_mfdataset(paths, combine='by_coords',
                                 data_vars='minimal', coords='minimal',
                                 compat='override', chunks=chunks)

    def load(self):
        self._validate()
        ds = self._open(self.fpath_cobalt_monthly,
                        chunks=self.chunks).rename(**self.cobalt_rename)

        # lat/lon have to be coordinates to survive the ds[self.vars] subset below
        promote = [c for c in ('lat', 'lon') if c in ds.data_vars]
        if promote:
            ds = ds.set_coords(promote)
        missing = [v for v in self.vars if v not in ds.data_vars]
        if missing:
            raise KeyError(
                f'variables {missing} not found in {self.fpath_cobalt_monthly!r}')
        ds = ds[self.vars]

        # a time coordinate, not just a time dimension, is what orders the months
        if 'time' not in ds.dims or 'time' not in ds.coords:
            raise ValueError("the monthly source needs a 'time' dimension carrying a "
                             'time coordinate; map its month dimension through '
                             'cobalt_rename')
        ds = ds.sortby('time')
        if ds.sizes['time'] != 12:
            raise ValueError('expected 12 monthly climatology steps in '
                             f"{self.fpath_cobalt_monthly!r}, got {ds.sizes['time']}")

        if not self.stream:
            flooded = xr.merge((bnd.flood_missing(ds[v], **self.cobalt_renamed_dims)
                                for v in ds.data_vars))
            # load required before xesmf can recognize coordinates, and flooding
            # loses the 2D lat/lon, so put them back
            ds = flooded.load().assign_coords(lat=ds['lat'], lon=ds['lon'])

        self.ds = ds.assign_coords(time=CobaltBoundaryMonthly.clim_time)
        self.ds['time'].attrs.update(self.time_attrs)
        self.hgrid = xr.open_dataset(self.grid_file)
        return self

    def cobaltv2_to_v3(self):
        """Convert COBALTv2 variables to v3, for whichever parents are present.

        The monthly source holds only part of the tracer set, so only the derived
        variables whose parent is here can be built. The others belong to the annual
        source, and are CobaltBoundary's job.

        Eager mode adds them to self.ds. Streaming mode registers them instead, so
        that each is derived from its parent right after the parent is flooded in
        export(). That keeps the number of flood_kara calls the same either way.
        """
        available = {name: (parent, divisor)
                     for name, (parent, divisor)
                     in CobaltBoundaryMonthly.v2_to_v3.items()
                     if parent in self.ds.data_vars}

        if self.stream:
            self._derived = available
            return self

        for name, (parent, divisor) in available.items():
            self.ds[name] = self.ds[parent] if divisor == 1.0 \
                            else self.ds[parent] / divisor
        return self

    def export(self):
        segments = [
            bnd.Segment(s['id'], s['border'], self.hgrid,
                        regrid_dir=self.cache_dir,
                        output_dir=self.output_dir)
            for s in self.segments
        ]

        regridded = {seg.segstr: [] for seg in segments}

        for v in self.ds.data_vars:
            source = self.ds[v]
            if self.stream:
                source = bnd.flood_missing(
                    source, **self.cobalt_renamed_dims).load()
                source = source.assign_coords(lat=self.ds['lat'], lon=self.ds['lon'])

            # Derived in here so that the flooded parent is reused rather than
            # reflooded. regrid_tracer takes its output name from field.name, hence
            # the rename rather than carrying the name alongside.
            fields = [source]
            for name, (parent, divisor) in self._derived.items():
                if parent == v:
                    fields.append((source if divisor == 1.0
                                   else source / divisor).rename(name))

            for field in fields:
                for seg in segments:
                    # a suffix of its own, so the weight cache cannot collide with
                    # CobaltBoundary writing into the same directory
                    regridded[seg.segstr].append(
                        seg.regrid_tracer(field,
                                          regrid_suffix='cobalt_monthly',
                                          flood=False,
                                          periodic=False,
                                          write=False))

            # source and fields are rebound on the next iteration, which releases the
            # global 3D field. Only the small segment arrays are carried forward.

        for seg in segments:
            cobalt_seg = xr.merge(regridded[seg.segstr])
            for v in cobalt_seg.data_vars:
                cobalt_seg[v] = np.clip(cobalt_seg[v], 0.0, None)
            cobalt_seg = seg.add_coords(cobalt_seg)
            # 'modulo' here also stops to_netcdf from applying a gregorian calendar
            cobalt_seg['time'].attrs.update(self.time_attrs)
            seg.to_netcdf(cobalt_seg, self.output_name)
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

    # These pairs convert each dataset's native variable and dimension names into
    # the names the boundary code requires. They describe the source file, not this
    # script, so they belong in the config beside the paths. The annual and monthly
    # sources are usually different products with different conventions: a GFDL ESM
    # run uses geolat_t/geolon_t/st_ocean, a MOM6 run uses xh/yh/z_l.
    #
    # cobalt_rename maps native names to lat/lon/z/time.
    # cobalt_renamed_dims then says which dimension is x, y and z *after* that
    # rename has been applied, because flooding and regridding happen afterwards.
    # So for a GFDL file the dims stay xt_ocean/yt_ocean (they are not renamed),
    # while for a MOM6 file xh/yh have become lon/lat by that point.
    cobalt_rename = config.get(
        'cobalt_rename', {'geolat_t': 'lat', 'geolon_t': 'lon', 'st_ocean': 'z'})
    cobalt_renamed_dims = config.get(
        'cobalt_renamed_dims', dict(xdim='xt_ocean', ydim='yt_ocean', zdim='z'))

    # the monthly source falls back to the annual conventions when unspecified
    monthly_rename = config.get('cobalt_monthly_rename', cobalt_rename)
    monthly_renamed_dims = config.get('cobalt_monthly_renamed_dims',
                                      cobalt_renamed_dims)

    # Dask chunking, keyed by NATIVE dimension names because it is applied at open
    # time, before the rename. Unset means no dask for a single file, and one chunk
    # per file for a glob or list.
    cobalt_chunks = config.get('cobalt_chunks')
    monthly_chunks = config.get('cobalt_monthly_chunks')

    time0 = dtt.datetime.strptime(str(config['time0']), '%Y-%m-%d')

    # Whether to run the COBALTv2 -> v3 conversion. A source that is already v3
    # carries nmd/simd/femd and the p* variables itself, and needs no conversion.
    # This gates the parent check as well as the conversion, so the two cannot drift
    # apart: the check exists only to serve the conversion.
    apply_v2_to_v3 = config.get('cobalt_v2_to_v3', True)

    # cobalt_monthly_file is optional: without it, the annual-only path below runs
    # exactly as before. With it, cobalt_vars and cobalt_monthly_vars must both be
    # listed explicitly in the config, and the two sources are written to separate
    # files by separate runs, so neither constrains the other's grid.
    monthly_file = config.get('cobalt_monthly_file')

    if monthly_file is None:
        annual = CobaltBoundary(fpath_cobalt=config['cobalt_file'],
                                grid_file=config['grid_file'],
                                output_dir=config['output_dir'],
                                cache_dir=config['cache'],
                                segments=config['segments'],
                                cobalt_rename=cobalt_rename,
                                cobalt_renamed_dims=cobalt_renamed_dims,
                                chunks=cobalt_chunks,
                                time0=time0).load()
        if apply_v2_to_v3:
            annual = annual.cobaltv2_to_v3()
        annual.export()
    else:
        monthly_vars = list(config['cobalt_monthly_vars'])
        annual_vars = [v for v in config['cobalt_vars'] if v not in monthly_vars]

        if annual_vars:
            # checked before load() so that a bad split fails in milliseconds
            # rather than after the flooding has been paid for
            if apply_v2_to_v3:
                require_v2_to_v3_parents(annual_vars)

            annual = CobaltBoundary(fpath_cobalt=config['cobalt_file'],
                                    grid_file=config['grid_file'],
                                    output_dir=config['output_dir'],
                                    cache_dir=config['cache'],
                                    segments=config['segments'],
                                    cobalt_rename=cobalt_rename,
                                    cobalt_renamed_dims=cobalt_renamed_dims,
                                    vars=annual_vars,
                                    chunks=cobalt_chunks,
                                    time0=time0).load()
            if apply_v2_to_v3:
                annual = annual.cobaltv2_to_v3()
            annual.export()

        monthly = CobaltBoundaryMonthly(
            fpath_cobalt_monthly=monthly_file,
            grid_file=config['grid_file'],
            output_dir=config['output_dir'],
            cache_dir=config['cache'],
            segments=config['segments'],
            vars=monthly_vars,
            cobalt_rename=monthly_rename,
            cobalt_renamed_dims=monthly_renamed_dims,
            chunks=monthly_chunks,
            stream=config.get('cobalt_stream', False)).load()
        # no check needed here: this derives only the variables whose parent it holds
        if apply_v2_to_v3:
            monthly = monthly.cobaltv2_to_v3()
        monthly.export()


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
    # cobalt_renamed_dims = dict(xdim='xt_ocean', ydim='yt_ocean', zdim='z')
    # time0 = dtt.datetime.strptime(str(config['time0']), '%Y-%m-%d')

    # (CobaltBoundary(fpath_cobalt=config['cobalt_file'],
    #                 grid_file=config['grid_file'],
    #                 output_dir=config['output_dir'],
    #                 cache_dir=config['cache'],
    #                 segments=config['segments'],
    #                 cobalt_rename=cobalt_rename,
    #                 cobalt_renamed_dims=cobalt_renamed_dims,
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

    # cobalt_flooded, hgrid = load_cobalt(fpath_cobalt, grid_file, cobalt_rename, cobalt_renamed_dims)
    # cobalt_flooded = cobaltv2_to_v3(cobalt_flooded)
    # export_segments(config, cobalt_flooded)



