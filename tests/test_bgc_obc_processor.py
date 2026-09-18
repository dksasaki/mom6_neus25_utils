"""Tests for CobaltBoundary and CobaltBoundaryMonthly.

The two classes handle two separate sources and write two separate files, so they
are tested side by side. CobaltBoundary has been run on real data;
CobaltBoundaryMonthly has not, so several of these exist specifically to find out
whether the newer class behaves like the older one.

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
DERIVED = ['nmd', 'simd', 'femd', 'psm', 'pmd', 'plg', 'pdi']
SOUTH = [{'id': 1, 'border': 'south'}]


def _monthly_kwargs(cobalt_rename, cobalt_renamed_dims, **overrides):
    """Constructor arguments for CobaltBoundaryMonthly, with overrides applied."""
    kwargs = dict(
        fpath_cobalt_monthly='monthly.nc',
        grid_file='hgrid.nc',
        output_dir='.',
        cache_dir='.',
        segments=SOUTH,
        vars=VARS,
        cobalt_rename=cobalt_rename,
        cobalt_renamed_dims=cobalt_renamed_dims,
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

def test_vars_cannot_be_empty(cobalt_rename, cobalt_renamed_dims):
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, cobalt_renamed_dims, vars=[]))
    with pytest.raises(ValueError, match='at least one tracer'):
        obj.load()


def test_cobalt_rename_must_provide_z(cobalt_rename, cobalt_renamed_dims):
    del cobalt_rename['st_ocean']
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, cobalt_renamed_dims))
    with pytest.raises(AssertionError, match="'z'"):
        obj.load()


def test_flood_rename_must_provide_all_three_dims(cobalt_rename, cobalt_renamed_dims):
    del cobalt_renamed_dims['zdim']
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, cobalt_renamed_dims))
    with pytest.raises(AssertionError, match='xdim'):
        obj.load()


def test_glob_matching_nothing_is_reported_with_the_pattern(
        cobalt_rename, cobalt_renamed_dims, tmp_path):
    pattern = str(tmp_path / 'absent_*.nc')
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, cobalt_renamed_dims,
                          fpath_cobalt_monthly=pattern))
    with pytest.raises(FileNotFoundError, match='absent_'):
        obj.load()


# --------------------------------------------------------- validation, needing files

def test_source_must_have_twelve_steps(
        cobalt_rename, cobalt_renamed_dims, hgrid_file, tmp_path, dirs):
    out, cache = dirs
    short = str(write_cobalt(tmp_path / 'six.nc', nmonths=6))
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, cobalt_renamed_dims,
                          fpath_cobalt_monthly=short, grid_file=hgrid_file,
                          output_dir=out, cache_dir=cache))
    with pytest.raises(ValueError, match='got 6'):
        obj.load()


def test_missing_variable_names_the_file(
        cobalt_rename, cobalt_renamed_dims, cobalt_monthly_file, dirs):
    out, cache = dirs
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, cobalt_renamed_dims,
                          fpath_cobalt_monthly=cobalt_monthly_file,
                          vars=VARS + ['absent_tracer'],
                          output_dir=out, cache_dir=cache))
    with pytest.raises(KeyError, match='absent_tracer'):
        obj.load()


def test_source_needs_a_time_coordinate(
        cobalt_rename, cobalt_renamed_dims, tmp_path, dirs):
    out, cache = dirs
    path = tmp_path / 'no_time_coord.nc'
    # a time dimension with no coordinate variable: the months cannot be ordered
    cobalt_dataset(nmonths=12).drop_vars('time').to_netcdf(
        path, format='NETCDF3_64BIT', engine='netcdf4')
    obj = bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, cobalt_renamed_dims,
                          fpath_cobalt_monthly=str(path),
                          output_dir=out, cache_dir=cache))
    with pytest.raises(ValueError, match='time coordinate'):
        obj.load()


# ------------------------------------------------------------------- CobaltBoundary

def _annual(fpath, hgrid_file, dirs, segments, cobalt_rename, cobalt_renamed_dims,
            vars=None):
    out, cache = dirs
    return bgc.CobaltBoundary(
        fpath_cobalt=fpath, grid_file=hgrid_file, output_dir=out, cache_dir=cache,
        segments=segments, cobalt_rename=cobalt_rename,
        cobalt_renamed_dims=cobalt_renamed_dims, vars=vars or VARS,
        time0=dtt.datetime(1993, 1, 1))


def test_annual_load_applies_time0(fake_flood, cobalt_annual_file, hgrid_file, dirs,
                                   segments, cobalt_rename, cobalt_renamed_dims):
    obj = _annual(cobalt_annual_file, hgrid_file, dirs, segments,
                  cobalt_rename, cobalt_renamed_dims).load()
    assert obj.ds['time'].values[0] == np.datetime64('1993-01-01')
    assert sorted(obj.ds.data_vars) == sorted(VARS)


def test_annual_v2_to_v3_derives_every_variable(
        fake_flood, cobalt_annual_file, hgrid_file, dirs, segments,
        cobalt_rename, cobalt_renamed_dims):
    obj = _annual(cobalt_annual_file, hgrid_file, dirs, segments,
                  cobalt_rename, cobalt_renamed_dims).load().cobaltv2_to_v3()
    for name in DERIVED:
        assert name in obj.ds.data_vars
    # the copies and the ratios the conversion is defined by
    assert obj.ds['nmd'].equals(obj.ds['nlg'])
    assert np.allclose(obj.ds['psm'].values, obj.ds['nsm'].values / 24.0)
    assert np.allclose(obj.ds['pmd'].values, obj.ds['nlg'].values / 20.0)


def test_annual_export_writes_one_file_per_segment(
        fake_flood, cobalt_annual_file, hgrid_file, dirs, segments,
        cobalt_rename, cobalt_renamed_dims):
    out, _ = dirs
    _annual(cobalt_annual_file, hgrid_file, dirs, segments,
            cobalt_rename, cobalt_renamed_dims).load().export()
    for seg in segments:
        ds = xr.open_dataset(f"{out}/bgc_cobalt_{seg['id']:03d}.nc")
        assert f"nlg_segment_{seg['id']:03d}" in ds.data_vars
        assert ds.sizes['time'] == 1


# ------------------------------------------------------------ CobaltBoundaryMonthly

def _monthly(fpath, hgrid_file, dirs, cobalt_rename, cobalt_renamed_dims,
             segments=None, **overrides):
    out, cache = dirs
    return bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, cobalt_renamed_dims,
                          fpath_cobalt_monthly=fpath, grid_file=hgrid_file,
                          output_dir=out, cache_dir=cache,
                          segments=segments or SOUTH, **overrides))


def test_monthly_load_puts_months_on_the_modulo_axis(
        fake_flood, cobalt_monthly_file, hgrid_file, dirs,
        cobalt_rename, cobalt_renamed_dims):
    obj = _monthly(cobalt_monthly_file, hgrid_file, dirs,
                   cobalt_rename, cobalt_renamed_dims).load()
    assert np.allclose(obj.ds['time'].values, bgc.CobaltBoundaryMonthly.clim_time)
    assert obj.ds['time'].attrs['modulo'] == ' '
    assert obj.ds['time'].attrs['calendar'] == 'noleap'
    assert obj.ds.sizes['time'] == 12


def test_monthly_load_holds_only_the_monthly_tracers(
        fake_flood, cobalt_monthly_file, hgrid_file, dirs,
        cobalt_rename, cobalt_renamed_dims):
    """The annual tracers belong to CobaltBoundary and must not appear here."""
    obj = _monthly(cobalt_monthly_file, hgrid_file, dirs, cobalt_rename,
                   cobalt_renamed_dims, vars=['nlg', 'nh4']).load()
    assert sorted(obj.ds.data_vars) == ['nh4', 'nlg']
    for name in obj.ds.data_vars:
        assert 'time' in obj.ds[name].dims


def test_monthly_v2_to_v3_derives_only_what_its_parents_allow(
        fake_flood, cobalt_monthly_file, hgrid_file, dirs,
        cobalt_rename, cobalt_renamed_dims):
    """With only nlg present, the silg/felg/nsm/ndi children cannot be built."""
    obj = _monthly(cobalt_monthly_file, hgrid_file, dirs, cobalt_rename,
                   cobalt_renamed_dims, vars=['nlg', 'nh4']
                   ).load().cobaltv2_to_v3()
    built = sorted(set(obj.ds.data_vars) - {'nlg', 'nh4'})
    assert built == ['nmd', 'plg', 'pmd']
    assert np.allclose(obj.ds['plg'].values, obj.ds['nlg'].values / 14.0)


@pytest.mark.parametrize('source', ['cobalt_monthly_glob', 'cobalt_monthly_files'])
def test_monthly_accepts_a_glob_or_a_list(
        fake_flood, request, source, hgrid_file, dirs,
        cobalt_rename, cobalt_renamed_dims):
    obj = _monthly(request.getfixturevalue(source), hgrid_file, dirs,
                   cobalt_rename, cobalt_renamed_dims).load()
    assert obj.ds.sizes['time'] == 12


def test_monthly_export_writes_its_own_file_and_varies_by_month(
        fake_flood, cobalt_monthly_file, hgrid_file, dirs,
        cobalt_rename, cobalt_renamed_dims):
    out, _ = dirs
    _monthly(cobalt_monthly_file, hgrid_file, dirs,
             cobalt_rename, cobalt_renamed_dims).load().export()

    # a name of its own, so it cannot overwrite CobaltBoundary's output
    ds = xr.open_dataset(f'{out}/bgc_cobalt_monthly_001.nc', decode_times=False)
    assert ds.sizes['time'] == 12
    for name in VARS:
        values = ds[f'{name}_segment_001'].values
        assert not np.allclose(values, values[0]), f'{name} should vary by month'


def test_monthly_export_writes_a_float_modulo_time_axis(
        fake_flood, cobalt_monthly_file, hgrid_file, dirs,
        cobalt_rename, cobalt_renamed_dims):
    out, _ = dirs
    _monthly(cobalt_monthly_file, hgrid_file, dirs,
             cobalt_rename, cobalt_renamed_dims).load().export()

    ds = xr.open_dataset(f'{out}/bgc_cobalt_monthly_001.nc', decode_times=False)
    assert ds['time'].dtype.kind == 'f'         # floats, not a decoded calendar
    assert ds['time'].attrs['modulo'] == ' '
    assert ds['time'].attrs['calendar'] == 'noleap'
    assert ds['time'].attrs['units'] == 'days since 0001-01-01'
    assert np.allclose(ds['time'].values, bgc.CobaltBoundaryMonthly.clim_time)


def test_monthly_export_output_is_finite_and_non_negative(
        fake_flood, cobalt_monthly_file, hgrid_file, dirs,
        cobalt_rename, cobalt_renamed_dims):
    out, _ = dirs
    _monthly(cobalt_monthly_file, hgrid_file, dirs,
             cobalt_rename, cobalt_renamed_dims).load().export()

    ds = xr.open_dataset(f'{out}/bgc_cobalt_monthly_001.nc', decode_times=False)
    # tracers and their layer thicknesses only. add_coords also writes
    # lon_segment_001 and lat_segment_001, and those longitudes are negative here.
    checked = ([f'{v}_segment_001' for v in VARS]
               + [f'dz_{v}_segment_001' for v in VARS])
    for name in checked:
        values = ds[name].values
        assert np.isfinite(values).all(), f'{name} has gaps below the sea floor'
        assert (values >= 0).all(), f'{name} was not clipped at zero'


def test_monthly_export_works_without_the_v2_to_v3_conversion(
        fake_flood, cobalt_monthly_file, hgrid_file, dirs,
        cobalt_rename, cobalt_renamed_dims):
    """The conversion is optional, so load().export() must stand on its own.

    main() exposes that by commenting out a line, so the chain has to be valid with
    cobaltv2_to_v3 never called.
    """
    out, _ = dirs
    _monthly(cobalt_monthly_file, hgrid_file, dirs,
             cobalt_rename, cobalt_renamed_dims).load().export()

    ds = xr.open_dataset(f'{out}/bgc_cobalt_monthly_001.nc', decode_times=False)
    for name in VARS:
        assert f'{name}_segment_001' in ds.data_vars
    for name in DERIVED:
        assert f'{name}_segment_001' not in ds.data_vars


def test_require_v2_to_v3_parents_accepts_a_complete_set():
    # no exception when all five are on the annual side
    bgc.require_v2_to_v3_parents(VARS)
    with pytest.raises(ValueError, match='nlg'):
        bgc.require_v2_to_v3_parents([v for v in VARS if v != 'nlg'])


def test_monthly_chunking_does_not_change_the_result(
        fake_flood, cobalt_monthly_file, hgrid_file, tmp_path,
        cobalt_rename, cobalt_renamed_dims):
    """Chunking is a memory strategy, so it must be invisible in the output.

    The chunks are keyed by native names, since they are applied at open time and
    cobalt_rename has not run yet. This is the full-vertical, split-horizontal
    scheme: flooding and regridding both read across the horizontal, so it is the
    arrangement most likely to disagree if something is chunk-unaware.
    """
    outputs = {}
    for label, chunks in [('plain', None),
                          ('chunked', {'st_ocean': -1,
                                       'yt_ocean': 4, 'xt_ocean': 4})]:
        out = tmp_path / label
        cache = tmp_path / f'cache_{label}'
        out.mkdir()
        cache.mkdir()
        (bgc.CobaltBoundaryMonthly(
            **_monthly_kwargs(cobalt_rename, cobalt_renamed_dims,
                              fpath_cobalt_monthly=cobalt_monthly_file,
                              grid_file=hgrid_file, output_dir=str(out),
                              cache_dir=str(cache), chunks=chunks))
            .load().cobaltv2_to_v3().export())
        outputs[label] = xr.open_dataset(out / 'bgc_cobalt_monthly_001.nc',
                                         decode_times=False)

    plain, chunked = outputs['plain'], outputs['chunked']
    assert sorted(plain.data_vars) == sorted(chunked.data_vars)
    for name in plain.data_vars:
        np.testing.assert_allclose(
            plain[name].values, chunked[name].values,
            err_msg=f'{name} changed when the source was chunked')


def test_monthly_eager_and_streaming_agree(
        fake_flood, cobalt_monthly_file, hgrid_file, tmp_path,
        cobalt_rename, cobalt_renamed_dims):
    """The two code paths are claimed to be equivalent, so prove it."""
    outputs = {}
    for stream in (False, True):
        out = tmp_path / f'stream_{stream}'
        cache = tmp_path / f'cache_{stream}'
        out.mkdir()
        cache.mkdir()
        (bgc.CobaltBoundaryMonthly(
            **_monthly_kwargs(cobalt_rename, cobalt_renamed_dims,
                              fpath_cobalt_monthly=cobalt_monthly_file,
                              grid_file=hgrid_file, output_dir=str(out),
                              cache_dir=str(cache), stream=stream))
            .load().cobaltv2_to_v3().export())
        outputs[stream] = xr.open_dataset(out / 'bgc_cobalt_monthly_001.nc',
                                          decode_times=False)

    eager, streamed = outputs[False], outputs[True]
    assert sorted(eager.data_vars) == sorted(streamed.data_vars)
    for name in eager.data_vars:
        np.testing.assert_allclose(
            eager[name].values, streamed[name].values,
            err_msg=f'{name} differs between eager and streaming')


def test_main_refuses_to_split_the_v2_to_v3_parents(tmp_path, monkeypatch):
    """Moving a conversion parent to the monthly source must fail explicably.

    CobaltBoundary.cobaltv2_to_v3 reads nlg, silg, felg, nsm and ndi from its own
    dataset, so main() checks before running rather than letting a KeyError surface
    from inside the chain. The guard fires before any file is opened.
    """
    import sys
    import yaml

    config = {'boundary': {
        'output_dir': str(tmp_path), 'cache': str(tmp_path),
        'grid_file': 'unused.nc', 'time0': '1993-01-01',
        'cobalt_file': 'unused.nc', 'woa_file': 'unused.nc',
        'segments': SOUTH,
        'cobalt_monthly_file': 'unused.nc',
        'cobalt_vars': VARS,
        'cobalt_monthly_vars': ['nlg', 'nh4'],   # nlg is a conversion parent
    }}
    path = tmp_path / 'config.yaml'
    path.write_text(yaml.safe_dump(config))
    monkeypatch.setattr(sys, 'argv', ['bgc_obc_processor', '--config', str(path)])

    with pytest.raises(ValueError, match='cobaltv2_to_v3'):
        bgc.main()


def test_monthly_matches_annual_on_the_same_field(
        fake_flood, cobalt_annual_file, hgrid_file, tmp_path,
        cobalt_rename, cobalt_renamed_dims):
    """Cross-check the new class against the one that has run on real data.

    Twelve identical copies of the annual field, put through the monthly class,
    must come back out as twelve identical copies of what CobaltBoundary produces
    from that same field. Both classes see the same variables, including everything
    cobaltv2_to_v3 derives.
    """
    annual_ds = xr.open_dataset(cobalt_annual_file)
    twelve = xr.concat([annual_ds.isel(time=0)] * 12, dim='time')
    twelve = twelve.assign_coords(
        time=np.array([np.datetime64(f'1993-{m:02d}-15', 'ns')
                       for m in range(1, 13)]))
    repeated = tmp_path / 'repeated.nc'
    twelve.to_netcdf(repeated, format='NETCDF3_64BIT', engine='netcdf4')

    out_a, out_m = tmp_path / 'annual', tmp_path / 'monthly'
    for d in (out_a, out_m):
        d.mkdir()

    bgc.CobaltBoundary(
        fpath_cobalt=cobalt_annual_file, grid_file=hgrid_file,
        output_dir=str(out_a), cache_dir=str(out_a), segments=SOUTH,
        cobalt_rename=cobalt_rename, cobalt_renamed_dims=cobalt_renamed_dims,
        vars=VARS, time0=dtt.datetime(1993, 1, 1)).load().cobaltv2_to_v3().export()

    bgc.CobaltBoundaryMonthly(
        **_monthly_kwargs(cobalt_rename, cobalt_renamed_dims,
                          fpath_cobalt_monthly=str(repeated),
                          grid_file=hgrid_file, output_dir=str(out_m),
                          cache_dir=str(out_m))).load().cobaltv2_to_v3().export()

    old = xr.open_dataset(out_a / 'bgc_cobalt_001.nc', decode_times=False)
    new = xr.open_dataset(out_m / 'bgc_cobalt_monthly_001.nc', decode_times=False)
    for name in VARS + DERIVED:
        var = f'{name}_segment_001'
        for month in range(12):
            np.testing.assert_allclose(
                new[var].isel(time=month).values, old[var].isel(time=0).values,
                err_msg=f'{var} differs from CobaltBoundary at month {month}')
