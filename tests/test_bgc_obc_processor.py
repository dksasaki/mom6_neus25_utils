"""Tests for CobaltBoundary and CobaltBoundaryMonthly.

The two classes are tested side by side on purpose. CobaltBoundary has been run on
real data; CobaltBoundaryMonthly has not, so several of these exist specifically to
find out whether the newer class behaves like the older one.

All of these use the fake_flood fixture, since HCtFlood is not installed here and
flooding is not what is under test.
"""

import datetime as dtt

import numpy as np
import pytest
import xarray as xr

import bgc_obc_processor as bgc
from synthetic import DEFAULT_COBALT_VARS, cobalt_dataset, write_cobalt

VARS = list(DEFAULT_COBALT_VARS)
MONTHLY_VARS = ['nlg', 'nh4']          # the rest come from the annual source
DERIVED = ['nmd', 'simd', 'femd', 'psm', 'pmd', 'plg', 'pdi']


def _monthly_kwargs(cobalt_rename, flood_missing_rename, **overrides):
    """Constructor arguments for CobaltBoundaryMonthly, with overrides applied."""
    kwargs = dict(
        fpath_cobalt='annual.nc',
        fpath_cobalt_monthly='monthly.nc',
        grid_file='hgrid.nc',
        output_dir='.',
        cache_dir='.',
        segments=[{'id': 1, 'border': 'south'}],
        vars=VARS,
        monthly_vars=MONTHLY_VARS,
        cobalt_rename=cobalt_rename,
        flood_missing_rename=flood_missing_rename,
    )
    kwargs.update(overrides)
    return kwargs


@pytest.fixture
def dirs(tmp_path):
    """Output and cache directories for a run."""
    out = tmp_path / 'out'
    cache = tmp_path / 'cache'
    out.mkdir()
    cache.mkdir()
    return str(out), str(cache)


# ------------------------------------------------------------- validation, no files

def test_monthly_vars_must_be_a_subset_of_vars(cobalt_rename, flood_missing_rename):
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, flood_missing_rename,
                          monthly_vars=['nlg', 'not_a_tracer']))
    with pytest.raises(ValueError, match='missing from vars'):
        obj.load()


def test_monthly_vars_cannot_be_empty(cobalt_rename, flood_missing_rename):
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, flood_missing_rename, monthly_vars=[]))
    with pytest.raises(ValueError, match='at least one variable'):
        obj.load()


def test_cobalt_rename_must_provide_z(cobalt_rename, flood_missing_rename):
    del cobalt_rename['st_ocean']
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, flood_missing_rename))
    with pytest.raises(AssertionError, match="'z'"):
        obj.load()


def test_flood_rename_must_provide_all_three_dims(cobalt_rename, flood_missing_rename):
    del flood_missing_rename['zdim']
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, flood_missing_rename))
    with pytest.raises(AssertionError, match='xdim'):
        obj.load()


def test_glob_matching_nothing_is_reported_with_the_pattern(
        cobalt_rename, flood_missing_rename, tmp_path):
    pattern = str(tmp_path / 'absent_*.nc')
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, flood_missing_rename,
                          fpath_cobalt_monthly=pattern))
    with pytest.raises(FileNotFoundError, match='absent_'):
        obj.load()


# --------------------------------------------------------- validation, needing files

def test_monthly_source_must_have_twelve_steps(
        cobalt_rename, flood_missing_rename, hgrid_file, tmp_path, dirs):
    out, cache = dirs
    short = str(write_cobalt(tmp_path / 'six.nc', nmonths=6))
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, flood_missing_rename,
                          fpath_cobalt_monthly=short, grid_file=hgrid_file,
                          output_dir=out, cache_dir=cache))
    with pytest.raises(ValueError, match='got 6'):
        obj.load()


def test_missing_variable_names_the_file(
        cobalt_rename, flood_missing_rename, cobalt_monthly_file, dirs):
    out, cache = dirs
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, flood_missing_rename,
                          fpath_cobalt_monthly=cobalt_monthly_file,
                          vars=VARS + ['absent_tracer'],
                          monthly_vars=['absent_tracer'],
                          output_dir=out, cache_dir=cache))
    with pytest.raises(KeyError, match='absent_tracer'):
        obj.load()


def test_monthly_source_needs_a_time_coordinate(
        cobalt_rename, flood_missing_rename, tmp_path, dirs):
    out, cache = dirs
    path = tmp_path / 'no_time_coord.nc'
    # a time dimension with no coordinate variable: the months cannot be ordered
    cobalt_dataset(nmonths=12).drop_vars('time').to_netcdf(
        path, format='NETCDF3_64BIT', engine='netcdf4')
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, flood_missing_rename,
                          fpath_cobalt_monthly=str(path),
                          output_dir=out, cache_dir=cache))
    with pytest.raises(ValueError, match='time coordinate'):
        obj.load()


# ------------------------------------------------------------------- CobaltBoundary

def _annual(fpath, hgrid_file, dirs, segments, cobalt_rename, flood_missing_rename):
    out, cache = dirs
    return bgc.CobaltBoundary(
        fpath_cobalt=fpath, grid_file=hgrid_file, output_dir=out, cache_dir=cache,
        segments=segments, cobalt_rename=cobalt_rename,
        flood_missing_rename=flood_missing_rename, vars=VARS,
        time0=dtt.datetime(1993, 1, 1))


def test_annual_load_applies_time0(fake_flood, cobalt_annual_file, hgrid_file, dirs,
                                  segments, cobalt_rename, flood_missing_rename):
    obj = _annual(cobalt_annual_file, hgrid_file, dirs, segments,
                  cobalt_rename, flood_missing_rename).load()
    assert obj.ds['time'].values[0] == np.datetime64('1993-01-01')
    assert sorted(obj.ds.data_vars) == sorted(VARS)


def test_annual_v2_to_v3_derives_every_variable(
        fake_flood, cobalt_annual_file, hgrid_file, dirs, segments,
        cobalt_rename, flood_missing_rename):
    obj = _annual(cobalt_annual_file, hgrid_file, dirs, segments,
                  cobalt_rename, flood_missing_rename).load().cobaltv2_to_v3()
    for name in DERIVED:
        assert name in obj.ds.data_vars
    # the copies and the ratios the conversion is defined by
    assert obj.ds['nmd'].equals(obj.ds['nlg'])
    assert np.allclose(obj.ds['psm'].values, obj.ds['nsm'].values / 24.0)
    assert np.allclose(obj.ds['pmd'].values, obj.ds['nlg'].values / 20.0)


def test_annual_export_writes_one_file_per_segment(
        fake_flood, cobalt_annual_file, hgrid_file, dirs, segments,
        cobalt_rename, flood_missing_rename):
    out, _ = dirs
    _annual(cobalt_annual_file, hgrid_file, dirs, segments,
            cobalt_rename, flood_missing_rename).load().export()
    for seg in segments:
        ds = xr.open_dataset(f"{out}/bgc_cobalt_{seg['id']:03d}.nc")
        assert f"nlg_segment_{seg['id']:03d}" in ds.data_vars
        assert ds.sizes['time'] == 1


# ------------------------------------------------------------ CobaltBoundaryMonthly

def _monthly(fpath, hgrid_file, dirs, segments, cobalt_rename,
             flood_missing_rename, annual=None, **overrides):
    out, cache = dirs
    return bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, flood_missing_rename,
                          fpath_cobalt=annual or fpath,
                          fpath_cobalt_monthly=fpath, grid_file=hgrid_file,
                          output_dir=out, cache_dir=cache, segments=segments,
                          **overrides))


def test_monthly_load_puts_months_on_the_modulo_axis(
        fake_flood, cobalt_monthly_file, cobalt_annual_file, hgrid_file, dirs,
        cobalt_rename, flood_missing_rename):
    obj = _monthly(cobalt_monthly_file, hgrid_file, dirs,
                   [{'id': 1, 'border': 'south'}], cobalt_rename,
                   flood_missing_rename, annual=cobalt_annual_file).load()
    assert np.allclose(obj.ds['time'].values,
                       bgc.CobaltBoundaryMonthly.clim_time)
    assert obj.ds['time'].attrs['modulo'] == ' '
    assert obj.ds['time'].attrs['calendar'] == 'noleap'


def test_monthly_load_strips_time_from_the_annual_variables(
        fake_flood, cobalt_monthly_file, cobalt_annual_file, hgrid_file, dirs,
        cobalt_rename, flood_missing_rename):
    obj = _monthly(cobalt_monthly_file, hgrid_file, dirs,
                   [{'id': 1, 'border': 'south'}], cobalt_rename,
                   flood_missing_rename, annual=cobalt_annual_file).load()
    for name in MONTHLY_VARS:
        assert 'time' in obj.ds[name].dims
    for name in set(VARS) - set(MONTHLY_VARS):
        assert 'time' not in obj.ds[name].dims
    assert obj.ds.sizes['time'] == 12


@pytest.mark.parametrize('source', ['cobalt_monthly_glob', 'cobalt_monthly_files'])
def test_monthly_accepts_a_glob_or_a_list(
        fake_flood, request, source, cobalt_annual_file, hgrid_file, dirs,
        cobalt_rename, flood_missing_rename):
    obj = _monthly(request.getfixturevalue(source), hgrid_file, dirs,
                   [{'id': 1, 'border': 'south'}], cobalt_rename,
                   flood_missing_rename, annual=cobalt_annual_file).load()
    assert obj.ds.sizes['time'] == 12


def test_monthly_export_holds_annual_fields_constant_over_the_months(
        fake_flood, cobalt_monthly_file, cobalt_annual_file, hgrid_file, dirs,
        cobalt_rename, flood_missing_rename):
    out, _ = dirs
    _monthly(cobalt_monthly_file, hgrid_file, dirs,
             [{'id': 1, 'border': 'south'}], cobalt_rename,
             flood_missing_rename, annual=cobalt_annual_file).load().export()

    ds = xr.open_dataset(f'{out}/bgc_cobalt_001.nc', decode_times=False)
    assert ds.sizes['time'] == 12

    annual_name = 'nsm_segment_001'        # from the annual source
    monthly_name = 'nlg_segment_001'       # from the monthly source
    annual = ds[annual_name].values
    assert np.allclose(annual, annual[0]), 'annual field should not vary by month'
    monthly = ds[monthly_name].values
    assert not np.allclose(monthly, monthly[0]), 'monthly field should vary by month'


def test_monthly_export_writes_a_float_modulo_time_axis(
        fake_flood, cobalt_monthly_file, cobalt_annual_file, hgrid_file, dirs,
        cobalt_rename, flood_missing_rename):
    out, _ = dirs
    _monthly(cobalt_monthly_file, hgrid_file, dirs,
             [{'id': 1, 'border': 'south'}], cobalt_rename,
             flood_missing_rename, annual=cobalt_annual_file).load().export()

    ds = xr.open_dataset(f'{out}/bgc_cobalt_001.nc', decode_times=False)
    assert ds['time'].dtype.kind == 'f'         # floats, not a decoded calendar
    assert ds['time'].attrs['modulo'] == ' '
    assert ds['time'].attrs['calendar'] == 'noleap'
    assert ds['time'].attrs['units'] == 'days since 0001-01-01'
    assert np.allclose(ds['time'].values, bgc.CobaltBoundaryMonthly.clim_time)


def test_monthly_export_output_is_finite_and_non_negative(
        fake_flood, cobalt_monthly_file, cobalt_annual_file, hgrid_file, dirs,
        cobalt_rename, flood_missing_rename):
    out, _ = dirs
    _monthly(cobalt_monthly_file, hgrid_file, dirs,
             [{'id': 1, 'border': 'south'}], cobalt_rename,
             flood_missing_rename, annual=cobalt_annual_file).load().export()

    ds = xr.open_dataset(f'{out}/bgc_cobalt_001.nc', decode_times=False)
    # tracers and their layer thicknesses only. add_coords also writes
    # lon_segment_001 and lat_segment_001, and those longitudes are negative here.
    checked = ([f'{v}_segment_001' for v in VARS]
               + [f'dz_{v}_segment_001' for v in VARS])
    for name in checked:
        values = ds[name].values
        assert np.isfinite(values).all(), f'{name} has gaps below the sea floor'
        assert (values >= 0).all(), f'{name} was not clipped at zero'


def test_monthly_v2_to_v3_matches_between_eager_and_streaming(
        fake_flood, cobalt_monthly_file, cobalt_annual_file, hgrid_file, tmp_path,
        cobalt_rename, flood_missing_rename):
    """The two code paths are claimed to be equivalent, so prove it."""
    outputs = {}
    for stream in (False, True):
        out = tmp_path / f'stream_{stream}'
        cache = tmp_path / f'cache_{stream}'
        out.mkdir()
        cache.mkdir()
        (bgc.CobaltBoundaryMonthly(
            **_monthly_kwargs(cobalt_rename, flood_missing_rename,
                              fpath_cobalt=cobalt_annual_file,
                              fpath_cobalt_monthly=cobalt_monthly_file,
                              grid_file=hgrid_file, output_dir=str(out),
                              cache_dir=str(cache), stream=stream))
            .load().cobaltv2_to_v3().export())
        outputs[stream] = xr.open_dataset(out / 'bgc_cobalt_001.nc',
                                          decode_times=False)

    eager, streamed = outputs[False], outputs[True]
    assert sorted(eager.data_vars) == sorted(streamed.data_vars)
    for name in eager.data_vars:
        np.testing.assert_allclose(
            eager[name].values, streamed[name].values,
            err_msg=f'{name} differs between eager and streaming')


def test_monthly_matches_annual_when_every_month_is_identical(
        fake_flood, cobalt_annual_file, hgrid_file, tmp_path,
        cobalt_rename, flood_missing_rename):
    """Cross-check the new class against the one that has run on real data.

    Twelve identical copies of the annual field must come back out as twelve
    identical copies of what CobaltBoundary produces from that same field.
    """
    annual_ds = xr.open_dataset(cobalt_annual_file)
    twelve = xr.concat([annual_ds.isel(time=0)] * 12, dim='time')
    twelve = twelve.assign_coords(
        time=np.array([np.datetime64(f'1993-{m:02d}-15', 'ns')
                       for m in range(1, 13)]))
    repeated = tmp_path / 'repeated.nc'
    twelve.to_netcdf(repeated, format='NETCDF3_64BIT', engine='netcdf4')

    paths = {}
    for tag in ('annual', 'monthly'):
        out = tmp_path / f'out_{tag}'
        cache = tmp_path / f'cache_{tag}'
        out.mkdir()
        cache.mkdir()
        paths[tag] = out
        if tag == 'annual':
            bgc.CobaltBoundary(
                fpath_cobalt=cobalt_annual_file, grid_file=hgrid_file,
                output_dir=str(out), cache_dir=str(cache),
                segments=[{'id': 1, 'border': 'south'}],
                cobalt_rename=cobalt_rename,
                flood_missing_rename=flood_missing_rename,
                vars=VARS, time0=dtt.datetime(1993, 1, 1)).load().export()
        else:
            bgc.CobaltBoundaryMonthly(
                **_monthly_kwargs(cobalt_rename, flood_missing_rename,
                                  fpath_cobalt=cobalt_annual_file,
                                  fpath_cobalt_monthly=str(repeated),
                                  monthly_vars=VARS,
                                  grid_file=hgrid_file, output_dir=str(out),
                                  cache_dir=str(cache))).load().export()

    old = xr.open_dataset(paths['annual'] / 'bgc_cobalt_001.nc', decode_times=False)
    new = xr.open_dataset(paths['monthly'] / 'bgc_cobalt_001.nc', decode_times=False)
    for name in VARS:
        var = f'{name}_segment_001'
        for month in range(12):
            np.testing.assert_allclose(
                new[var].isel(time=month).values, old[var].isel(time=0).values,
                err_msg=f'{var} differs from CobaltBoundary at month {month}')
