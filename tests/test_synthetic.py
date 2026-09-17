"""Tests for the synthetic inputs themselves.

Every other test leans on these builders, so it is worth pinning down their
contract: the supergrid arithmetic, the agreement between the grid files, the
positive-down depth convention, and the NaN structure below the sea floor.
"""

import numpy as np
import pytest
import xarray as xr

from synthetic import (_supergrid_axes, bathymetry, cobalt_dataset,
                       hgrid_dataset, topog_dataset, write_cobalt)


# --------------------------------------------------------------------------- grids

def test_hgrid_has_supergrid_dimensions():
    ds = hgrid_dataset(ni=12, nj=10)
    # a corner at each end, hence odd and twice the model resolution
    assert ds.sizes['nxp'] == 2 * 12 + 1
    assert ds.sizes['nyp'] == 2 * 10 + 1
    assert ds.sizes['nx'] == ds.sizes['nxp'] - 1
    assert ds.sizes['ny'] == ds.sizes['nyp'] - 1


def test_hgrid_x_and_y_are_two_dimensional():
    ds = hgrid_dataset()
    # real curvilinear grids are 2D even where a test grid would be separable
    assert ds['x'].dims == ('nyp', 'nxp')
    assert ds['y'].dims == ('nyp', 'nxp')


def test_hgrid_spacing_is_half_the_model_resolution():
    ds = hgrid_dataset(ni=12, nj=10, dlon=0.25, dlat=0.25)
    assert np.allclose(np.diff(ds['x'].values[0, :]), 0.125)
    assert np.allclose(np.diff(ds['y'].values[:, 0]), 0.125)


def test_hgrid_angle_is_zero_so_check_angle_range_passes():
    ds = hgrid_dataset()
    assert ds['angle_dx'].attrs['units'] == 'degrees'
    assert np.all(ds['angle_dx'].values == 0.0)


def test_topog_is_on_the_model_grid_not_the_supergrid():
    hg = hgrid_dataset(ni=12, nj=10)
    tp = topog_dataset(ni=12, nj=10)
    assert tp['depth'].dims == ('ny', 'nx')
    assert tp.sizes['nx'] == (hg.sizes['nxp'] - 1) // 2
    assert tp.sizes['ny'] == (hg.sizes['nyp'] - 1) // 2


def test_topog_has_no_coordinate_variables():
    # all georeferencing lives in the hgrid, matching the real ocean_topog.nc
    assert list(topog_dataset().coords) == []


def test_topog_samples_the_odd_odd_h_points():
    ni, nj = 12, 10
    lon, lat = _supergrid_axes(ni, nj, -72.0, 36.0, 0.25, 0.25)
    expected = bathymetry(lon[1::2][np.newaxis, :], lat[1::2][:, np.newaxis])
    assert np.allclose(topog_dataset(ni=ni, nj=nj)['depth'].values, expected)


def test_depth_is_positive_down_and_zero_at_the_coast():
    depth = topog_dataset()['depth'].values
    assert depth.min() == 0.0          # zero water thickness, i.e. land
    assert depth.max() > 0.0           # positive going deeper
    assert not (depth < 0).any()


def test_bathymetry_is_asymmetric_in_i_and_j():
    # a transposed-index bug has to change the answer rather than go unnoticed
    depth = topog_dataset(ni=6, nj=6)['depth'].values
    assert not np.allclose(depth, depth.T)


# ------------------------------------------------------------------------- sources

def test_cobalt_uses_native_gfdl_names_and_dim_order():
    ds = cobalt_dataset()
    assert ds['nlg'].dims == ('time', 'st_ocean', 'yt_ocean', 'xt_ocean')
    for name in ['geolat_t', 'geolon_t', 'st_ocean', 'xt_ocean', 'yt_ocean']:
        assert name in ds.coords


def test_cobalt_vertical_axis_starts_at_the_surface_and_deepens():
    z = cobalt_dataset()['st_ocean'].values
    assert z[0] == 0.0
    assert np.all(np.diff(z) > 0)
    # levels thicken with depth, as in an ocean model
    assert np.all(np.diff(np.diff(z)) > 0)


def test_cobalt_source_encloses_and_out_coarsens_the_model_domain():
    hg = hgrid_dataset()
    ds = cobalt_dataset()
    assert ds['xt_ocean'].min() < hg['x'].min()
    assert ds['xt_ocean'].max() > hg['x'].max()
    assert ds['yt_ocean'].min() < hg['y'].min()
    assert ds['yt_ocean'].max() > hg['y'].max()
    model_spacing = float(np.diff(hg['x'].values[0, :])[0]) * 2
    assert float(np.diff(ds['xt_ocean'].values)[0]) > model_spacing


def test_cobalt_is_dry_only_below_wet():
    # NaNs must form a sea floor, never a bubble in the middle of the column
    finite = np.isfinite(cobalt_dataset()['nlg'].isel(time=0).values)
    assert np.all(finite[:-1] | ~finite[1:])


def test_cobalt_land_columns_are_entirely_nan():
    ds = cobalt_dataset()
    lon2d, lat2d = np.meshgrid(ds['xt_ocean'].values, ds['yt_ocean'].values)
    land = bathymetry(lon2d, lat2d) == 0.0
    column = ds['nlg'].isel(time=0).values
    assert land.any(), 'the fake domain should contain some land'
    assert np.all(np.isnan(column[:, land]))
    assert np.isfinite(column[0][~land]).all()


def test_cobalt_source_floor_can_sit_below_the_model_floor():
    ds = cobalt_dataset(depth_scale=1.25)
    deepest_wet = ds['st_ocean'].values[
        np.isfinite(ds['nlg'].isel(time=0).values).any(axis=(1, 2))].max()
    assert deepest_wet > topog_dataset()['depth'].values.max()


def test_cobalt_values_identify_their_variable_and_month():
    ds = cobalt_dataset(nmonths=12)
    point = dict(st_ocean=0, yt_ocean=-1, xt_ocean=-1)
    assert float(ds['nlg'].isel(time=0, **point)) != float(ds['nsm'].isel(time=0, **point))
    assert float(ds['nlg'].isel(time=0, **point)) != float(ds['nlg'].isel(time=11, **point))
    assert np.nanmin(ds['nlg'].values) > 0


def test_cobalt_time_can_be_shuffled_and_sorted_back():
    ordered = cobalt_dataset(nmonths=12)['time'].values
    shuffled = cobalt_dataset(nmonths=12, shuffle_time=True)
    assert not (shuffled['time'].values == ordered).all()
    assert (shuffled.sortby('time')['time'].values == ordered).all()


# --------------------------------------------------------------------------- files

@pytest.mark.parametrize('builder,name', [(hgrid_dataset, 'hgrid'),
                                          (topog_dataset, 'topog')])
def test_written_files_are_netcdf3(tmp_path, builder, name):
    path = tmp_path / f'{name}.nc'
    builder().to_netcdf(path, format='NETCDF3_64BIT', engine='netcdf4')
    with open(path, 'rb') as f:
        assert f.read(4) == b'CDF\x02'   # 64-bit offset netCDF3


def test_split_by_month_writes_one_file_per_month(cobalt_monthly_files):
    assert len(cobalt_monthly_files) == 12
    for path in cobalt_monthly_files:
        assert xr.open_dataset(path).sizes['time'] == 1


def test_static_lonlat_gains_a_time_dim_unless_open_is_told_otherwise(tmp_path):
    """The open_mfdataset trap that CobaltBoundaryMonthly._open has to avoid.

    With geolat_t stored as a data variable, the default settings concatenate it
    along time, leaving xesmf a 3D latitude it cannot build a grid from.
    """
    d = tmp_path / 'dv'
    d.mkdir()
    paths = write_cobalt(d, split_by_month=True, nmonths=3, lonlat_as_coords=False)

    trap = xr.open_mfdataset(paths, combine='by_coords')
    assert trap['geolat_t'].dims == ('time', 'yt_ocean', 'xt_ocean')

    fixed = xr.open_mfdataset(paths, combine='by_coords', data_vars='minimal',
                              coords='minimal', compat='override')
    assert fixed['geolat_t'].dims == ('yt_ocean', 'xt_ocean')


# ------------------------------------------------------------------------ fixtures

def test_segment_accepts_the_synthetic_hgrid(hgrid_file):
    import boundary as bnd
    expected = {'south': 25, 'north': 25, 'west': 21, 'east': 21}
    for border, npts in expected.items():
        # a fresh open per segment: Segment converts angle_dx in place
        seg = bnd.Segment(1, border, xr.open_dataset(hgrid_file))
        assert len(seg.coords['lon']) == npts


def test_fake_flood_fills_land_but_keeps_empty_levels(fake_flood):
    ds = cobalt_dataset()
    flooded = fake_flood(ds['nlg'], **{'xdim': 'xt_ocean', 'ydim': 'yt_ocean',
                                       'zdim': 'st_ocean'})
    assert flooded.name == 'nlg'
    assert flooded.dims == ds['nlg'].dims
    # the surface has wet points to spread from, so it fills completely
    assert np.isfinite(flooded.isel(time=0, st_ocean=0).values).all()


def test_fake_flood_is_installed_over_the_real_one(fake_flood):
    import boundary as bnd
    assert bnd.flood_missing is fake_flood
