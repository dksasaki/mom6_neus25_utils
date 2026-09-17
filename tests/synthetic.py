"""Synthetic MOM6 inputs for testing, small enough to check by hand.

Nothing here reads the cluster, so the boundary code can be exercised on any
machine. This is test scaffolding rather than shipped code, which is why it lives
under tests/ and not in the package.

It deliberately has no pytest dependency, so it is equally usable from a plain
script or a notebook.
"""

import hashlib
import os

import numpy as np
import xarray as xr

# spherical earth, only used for the optional dx/dy/area metrics
EARTH_RADIUS = 6.371e6

# every file MOM6 reads is netCDF3, matching Segment.to_netcdf
NETCDF_FORMAT = 'NETCDF3_64BIT'

# A small but sufficient tracer set: every parent that cobaltv2_to_v3 derives from
# (nlg, silg, felg, nsm, ndi) plus two extras, so there is something left over to
# split between a monthly and an annual source.
DEFAULT_COBALT_VARS = ('nlg', 'silg', 'felg', 'nsm', 'ndi', 'nh4', 'ldon')


def _supergrid_axes(ni, nj, lon0, lat0, dlon, dlat):
    """Supergrid longitude and latitude axes for an ni x nj model grid.

    Shared by the hgrid and topography builders so the two cannot drift apart. The
    supergrid steps by half the model resolution and has an odd number of points,
    because it carries a corner at each end.

    Returns:
        tuple: (lon, lat) 1D arrays of length 2*ni+1 and 2*nj+1.
    """
    lon = lon0 + np.arange(2 * ni + 1) * (dlon / 2.0)
    lat = lat0 + np.arange(2 * nj + 1) * (dlat / 2.0)
    return lon, lat


def hgrid_dataset(ni=12, nj=10, lon0=-72.0, lat0=36.0, dlon=0.25, dlat=0.25,
                  with_metrics=True):
    """Build an ocean_hgrid.nc dataset for an ni x nj model grid.

    The result is a MOM6 supergrid: it carries every stagger of the Arakawa C grid,
    indexed by parity, so (even, even) are vorticity corners, (odd, odd) are tracer
    centres, (even, odd) are u points and (odd, even) are v points. That is why the
    supergrid holds 2*ni+1 by 2*nj+1 points and steps by half the model resolution,
    and why a boundary segment ends up 2*ni+1 points long.

    The grid is lon/lat aligned, so angle_dx is zero everywhere. Segment only uses
    angle_dx to rotate velocities, and check_angle_range merely guards against
    degrees/radians confusion, so zero is both valid and irrelevant to tracers.

    Args:
        ni (int): Model cells in the x direction. The supergrid gets 2*ni+1 points.
        nj (int): Model cells in the y direction. The supergrid gets 2*nj+1 points.
        lon0 (float): Longitude of the southwest corner, in degrees east.
        lat0 (float): Latitude of the southwest corner, in degrees north.
        dlon (float): Model cell size in longitude. Supergrid spacing is half this.
        dlat (float): Model cell size in latitude. Supergrid spacing is half this.
        with_metrics (bool): Also emit dx, dy, area and tile, so the file matches a
            real ocean_hgrid.nc. Segment needs only x, y and angle_dx. Note that these
            metrics describe *supergrid* cells, so a model cell area is the sum of its
            four sub-cells.

    Returns:
        xarray.Dataset: Dataset with the variables and dimensions of ocean_hgrid.nc.
    """
    nxp = 2 * ni + 1
    nyp = 2 * nj + 1
    lon, lat = _supergrid_axes(ni, nj, lon0, lat0, dlon, dlat)

    # x and y are 2D even though this grid is separable, as on a real curvilinear grid
    x = np.broadcast_to(lon[np.newaxis, :], (nyp, nxp)).astype('float64')
    y = np.broadcast_to(lat[:, np.newaxis], (nyp, nxp)).astype('float64')

    ds = xr.Dataset({
        'x': (('nyp', 'nxp'), x.copy(), {'units': 'degrees_east'}),
        'y': (('nyp', 'nxp'), y.copy(), {'units': 'degrees_north'}),
        'angle_dx': (('nyp', 'nxp'), np.zeros((nyp, nxp)), {'units': 'degrees'})
    })

    if with_metrics:
        dlon_rad = np.radians(dlon / 2.0)
        dlat_rad = np.radians(dlat / 2.0)
        lat_rad = np.radians(lat)

        # dx spans adjacent supergrid points in i, hence (nyp, nxp-1)
        dx = EARTH_RADIUS * np.cos(lat_rad)[:, np.newaxis] * dlon_rad
        dx = np.broadcast_to(dx, (nyp, nxp - 1)).copy()
        # dy spans adjacent supergrid points in j, hence (nyp-1, nxp)
        dy = np.full((nyp - 1, nxp), EARTH_RADIUS * dlat_rad)
        # exact spherical area of each supergrid cell, (nyp-1, nxp-1)
        band = np.sin(lat_rad[1:]) - np.sin(lat_rad[:-1])
        area = EARTH_RADIUS ** 2 * dlon_rad * band[:, np.newaxis]
        area = np.broadcast_to(area, (nyp - 1, nxp - 1)).copy()

        ds['dx'] = (('nyp', 'nx'), dx, {'units': 'm'})
        ds['dy'] = (('ny', 'nxp'), dy, {'units': 'm'})
        ds['area'] = (('ny', 'nx'), area, {'units': 'm2'})
        ds['tile'] = (('string255', ),
                      np.array(list('tile1'.ljust(255)), dtype='S1'))

    ds.attrs['history'] = 'synthetic grid from tests/synthetic.py'
    return ds


def bathymetry(lon, lat, lon_coast=-72.0, lat_ref=36.0, coast_tilt=0.2,
               shelf_width=2.0, max_depth=2000.0):
    """Analytic ocean depth in metres: zero at the surface, positive downward.

    This is the single source of truth for where the sea floor is. Sampling it onto
    both the model h-points and the coarse source grid keeps the fake topography and
    the fake forcing consistent, without either file having to know the other's grid.

    Depth is zero on the landward side of the coastline, which is also MOM6's land
    marker in ocean_topog.nc: zero water thickness.

    The coastline is deliberately tilted so the field is asymmetric in i and j. A bug
    that transposes the two indices then changes the answer instead of hiding.

    Args:
        lon: Longitudes in degrees east. Broadcast against lat.
        lat: Latitudes in degrees north. Broadcast against lon.
        lon_coast (float): Longitude of the coastline at lat_ref. Land lies west of it.
        lat_ref (float): Latitude at which the coastline sits at lon_coast.
        coast_tilt (float): Degrees of longitude the coast shifts per degree of latitude.
        shelf_width (float): Degrees of longitude over which depth ramps from 0 to max_depth.
        max_depth (float): Depth in metres reached at the seaward edge of the shelf.

    Returns:
        numpy.ndarray: Depth in metres, positive down, zero over land.
    """
    lon = np.asarray(lon, dtype='float64')
    lat = np.asarray(lat, dtype='float64')
    coast = lon_coast + coast_tilt * (lat - lat_ref)
    offshore = np.clip((lon - coast) / shelf_width, 0.0, 1.0)
    return max_depth * offshore


def topog_dataset(ni=12, nj=10, lon0=-72.0, lat0=36.0, dlon=0.25, dlat=0.25,
                  **bathy_kwargs):
    """Build an ocean_topog.nc dataset for an ni x nj model grid.

    Unlike the hgrid, this file is on the *model* grid, not the supergrid: depth is
    (nj, ni), which is why a real NEUS25 topography is 396 x 430 against a supergrid
    of 793 x 861. Depth is sampled at the model h-points, which are the (odd, odd)
    supergrid points, so the two files agree by construction given the same ni/nj.

    Matching the real file, there are no coordinate variables. All georeferencing
    lives in the hgrid, so this file alone cannot tell you where it is.

    Args:
        ni (int): Model cells in the x direction, giving the nx dimension.
        nj (int): Model cells in the y direction, giving the ny dimension.
        lon0 (float): Longitude of the southwest supergrid corner, in degrees east.
        lat0 (float): Latitude of the southwest supergrid corner, in degrees north.
        dlon (float): Model cell size in longitude.
        dlat (float): Model cell size in latitude.
        **bathy_kwargs: Passed to bathymetry(), e.g. max_depth or shelf_width.

    Returns:
        xarray.Dataset: Dataset with the variables and dimensions of ocean_topog.nc.
    """
    lon, lat = _supergrid_axes(ni, nj, lon0, lat0, dlon, dlat)
    # model h points are the (odd, odd) supergrid points, giving ni x nj of them
    lon_h = lon[1::2]
    lat_h = lat[1::2]

    depth = bathymetry(lon_h[np.newaxis, :], lat_h[:, np.newaxis], **bathy_kwargs)
    depth = np.ascontiguousarray(depth, dtype='float64')

    ds = xr.Dataset({
        'depth': (('ny', 'nx'), depth, {
            'units': 'meters',
            'standard_name': 'topographic depth at Arakawa C h-points',
            # genuinely computed, as gridtools does, rather than a placeholder string
            'sha256': hashlib.sha256(depth.tobytes()).hexdigest()
        })
    })
    # the real file carries _FillValue = NaN, which xarray writes from encoding
    ds['depth'].encoding['_FillValue'] = np.nan

    ds.attrs['grid_version'] = '0.2'
    ds.attrs['history'] = 'synthetic topography from tests/synthetic.py'
    return ds


def cobalt_dataset(vars=DEFAULT_COBALT_VARS, nmonths=1, nz=10,
                   lon_bounds=(-74.0, -67.0), lat_bounds=(34.0, 40.5), dsource=0.5,
                   source_max_depth=2400.0, depth_scale=1.25,
                   lonlat_as_coords=True, shuffle_time=False, **bathy_kwargs):
    """Build a dataset shaped like global COBALT output, on its own coarse grid.

    This stands in for the *source* of the boundary forcing, which is a different
    model from the one ocean_hgrid.nc describes. So it carries GFDL's native names,
    is coarser than the model grid, and extends well beyond it, exactly as a global
    run would. Those names are what cobalt_rename and flood_missing_rename map:
    geolat_t, geolon_t, st_ocean, xt_ocean, yt_ocean.

    The vertical axis starts at 0 at the surface and increases downward, with levels
    thickening with depth. Every level at or below the local sea floor is NaN, which
    is what real ocean output looks like and what the pipeline has to cope with:
    land columns are entirely NaN and get filled horizontally by flood_missing, while
    a water column that runs out of levels gets filled downward by regrid_tracer.

    Tracer values are deliberately distinguishable rather than realistic. Each carries
    a per-variable offset plus depth, month, latitude and longitude terms, so a test
    that picks up the wrong variable, month or transposed axis sees different numbers
    instead of a plausible-looking field.

    Args:
        vars (sequence): Tracer names to emit. The default includes every parent that
            cobaltv2_to_v3 needs, so the conversion can run against it.
        nmonths (int): Time steps. Use 1 for an annual climatology, 12 for a monthly one.
        nz (int): Number of depth levels.
        lon_bounds (tuple): (west, east) extent in degrees east. Should enclose the
            model domain, since a global source always does.
        lat_bounds (tuple): (south, north) extent in degrees north.
        dsource (float): Source grid spacing in degrees. Coarser than the model grid.
        source_max_depth (float): Deepest level, in metres. Kept just inside the
            deepest sea floor (max_depth * depth_scale) so that the bottom level is
            genuinely wet in the deep region, rather than NaN everywhere.
        depth_scale (float): The source sea floor is the analytic bathymetry times
            this. Above 1 puts the source floor below the model's, so every model
            level has data. Below 1 gives the shallow-source case that exercises the
            vertical backfill in regrid_tracer. Land stays land either way.
        lonlat_as_coords (bool): Emit geolat_t/geolon_t as coordinates. False emits
            them as data variables, which is the case that needs promoting before a
            variable subset would drop them.
        shuffle_time (bool): Reverse the time axis, so that code depending on
            chronological order has to sort rather than getting it for free.
        **bathy_kwargs: Passed to bathymetry(), which must match the topography's.

    Returns:
        xarray.Dataset: Dataset shaped like ocean_cobalt_tracers.*.nc.
    """
    vars = list(vars)
    lon = np.arange(lon_bounds[0], lon_bounds[1] + dsource / 2.0, dsource)
    lat = np.arange(lat_bounds[0], lat_bounds[1] + dsource / 2.0, dsource)
    nx, ny = lon.size, lat.size

    lon2d = np.broadcast_to(lon[np.newaxis, :], (ny, nx)).astype('float64')
    lat2d = np.broadcast_to(lat[:, np.newaxis], (ny, nx)).astype('float64')

    # levels thicken with depth, starting at 0 at the surface
    z = source_max_depth * (np.arange(nz) / max(nz - 1, 1)) ** 2

    # wet strictly above the sea floor, so a land column (floor 0) has no wet levels
    floor = bathymetry(lon2d, lat2d, **bathy_kwargs) * depth_scale
    wet = z[:, np.newaxis, np.newaxis] < floor[np.newaxis, :, :]

    # nanosecond precision up front, else xarray warns while converting
    time = np.array([np.datetime64(f'1993-{m:02d}-15', 'ns')
                     for m in range(1, nmonths + 1)])

    ds = xr.Dataset()
    for k, name in enumerate(vars):
        field = (100.0 * (k + 1)
                 + 10.0 * z[np.newaxis, :, np.newaxis, np.newaxis] / source_max_depth
                 + np.arange(nmonths)[:, np.newaxis, np.newaxis, np.newaxis]
                 + 0.5 * (lat2d - lat2d.min())[np.newaxis, np.newaxis, :, :]
                 + 0.25 * (lon2d - lon2d.min())[np.newaxis, np.newaxis, :, :])
        field = np.broadcast_to(field, (nmonths, nz, ny, nx)).astype('float64').copy()
        field[:, ~wet] = np.nan
        ds[name] = (('time', 'st_ocean', 'yt_ocean', 'xt_ocean'), field,
                    {'units': 'mol kg-1'})

    ds = ds.assign_coords(
        time=('time', time),
        st_ocean=('st_ocean', z, {'units': 'meters', 'positive': 'down',
                                  'cartesian_axis': 'Z'}),
        yt_ocean=('yt_ocean', lat, {'units': 'degrees_N'}),
        xt_ocean=('xt_ocean', lon, {'units': 'degrees_E'})
    )

    lonlat = {
        'geolon_t': (('yt_ocean', 'xt_ocean'), lon2d.copy(), {'units': 'degrees_E'}),
        'geolat_t': (('yt_ocean', 'xt_ocean'), lat2d.copy(), {'units': 'degrees_N'})
    }
    if lonlat_as_coords:
        ds = ds.assign_coords(**lonlat)
    else:
        for name, spec in lonlat.items():
            ds[name] = spec

    if shuffle_time:
        ds = ds.isel(time=slice(None, None, -1))

    ds.attrs['history'] = 'synthetic cobalt from tests/synthetic.py'
    return ds


def write_cobalt(path, split_by_month=False, **kwargs):
    """Write synthetic COBALT output and return the path or paths written.

    Args:
        path: Destination file, or destination directory when split_by_month is True.
        split_by_month (bool): Write one file per month instead of one file holding
            all of them. This is the layout that has to be opened with a glob or a
            list, and it is where open_mfdataset can wrongly concatenate the static
            geolat_t/geolon_t along time.
        **kwargs: Passed through to cobalt_dataset().

    Returns:
        The path written, or a sorted list of paths when split_by_month is True.
    """
    ds = cobalt_dataset(**kwargs)
    if not split_by_month:
        ds.to_netcdf(path, format=NETCDF_FORMAT, engine='netcdf4')
        return path

    paths = []
    for i in range(ds.sizes['time']):
        month = int(str(ds['time'].values[i])[5:7])
        fpath = os.path.join(path, f'cobalt_mon{month:02d}.nc')
        ds.isel(time=[i]).to_netcdf(fpath, format=NETCDF_FORMAT, engine='netcdf4')
        paths.append(fpath)
    return sorted(paths)


def write_hgrid(path, **kwargs):
    """Write a synthetic ocean_hgrid.nc and return its path.

    The boundary classes take file paths rather than datasets, so tests need a real
    file on disk. Pair this with pytest's tmp_path.

    Note that Segment converts angle_dx to radians in place on the dataset it is
    handed, so each Segment should be built from a freshly opened hgrid.

    Args:
        path: Destination for the netCDF file.
        **kwargs: Passed through to hgrid_dataset().

    Returns:
        The path that was written, unchanged, so it can be used inline.
    """
    hgrid_dataset(**kwargs).to_netcdf(path, format=NETCDF_FORMAT, engine='netcdf4')
    return path


def write_topog(path, **kwargs):
    """Write a synthetic ocean_topog.nc and return its path.

    Args:
        path: Destination for the netCDF file.
        **kwargs: Passed through to topog_dataset().

    Returns:
        The path that was written, unchanged, so it can be used inline.
    """
    topog_dataset(**kwargs).to_netcdf(path, format=NETCDF_FORMAT, engine='netcdf4')
    return path
