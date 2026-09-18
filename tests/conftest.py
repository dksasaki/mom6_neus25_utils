"""Shared fixtures for the test suite.

pytest loads this file automatically, so nothing here is ever imported by hand.
Loading it also puts this directory at the front of sys.path, which is what lets
``from synthetic import ...`` work below and keeps path handling out of the
individual test modules.

Run the suite with ``pixi run pytest``. Calling pytest directly skips environment
activation, and esmpy then fails on a missing ESMFMKFILE.
"""

import sys
from pathlib import Path

import pytest

# The boundary modules import each other by bare name ("import boundary as bnd")
# rather than as a package, and boundary/ has no __init__.py, so that directory has
# to be importable in its own right.
BOUNDARY_DIR = (Path(__file__).resolve().parents[1]
                / 'src' / 'mom6_neus25_utils' / 'boundary')
sys.path.insert(0, str(BOUNDARY_DIR))

from synthetic import write_cobalt, write_hgrid, write_topog  # noqa: E402

# the values main() passes, kept here so the tests stay in step with production
COBALT_RENAME = {'geolat_t': 'lat', 'geolon_t': 'lon', 'st_ocean': 'z'}
COBALT_RENAMED_DIMS = dict(xdim='xt_ocean', ydim='yt_ocean', zdim='z')


@pytest.fixture
def cobalt_rename():
    """Native cobalt names mapped to the names the boundary code requires."""
    return dict(COBALT_RENAME)


@pytest.fixture
def cobalt_renamed_dims():
    """Native cobalt dimension names, as flood_missing expects them."""
    return dict(COBALT_RENAMED_DIMS)


@pytest.fixture
def segments():
    """The two open boundaries of the real NEUS25 config: south and east.

    The synthetic domain has its land wedge in the west, so these two borders are
    over open water, as they are in production.
    """
    return [{'id': 1, 'border': 'south'}, {'id': 2, 'border': 'east'}]


@pytest.fixture
def hgrid_file(tmp_path):
    """Path to a synthetic ocean_hgrid.nc for a 12 x 10 model grid."""
    return str(write_hgrid(tmp_path / 'ocean_hgrid.nc'))


@pytest.fixture
def topog_file(tmp_path):
    """Path to a synthetic ocean_topog.nc matching hgrid_file's ni and nj."""
    return str(write_topog(tmp_path / 'ocean_topog.nc'))


@pytest.fixture
def cobalt_annual_file(tmp_path):
    """One file holding one time step, as the annual climatology does."""
    return str(write_cobalt(tmp_path / 'cobalt_ann.nc', nmonths=1))


@pytest.fixture
def cobalt_monthly_file(tmp_path):
    """One file holding all twelve months."""
    return str(write_cobalt(tmp_path / 'cobalt_mon.nc', nmonths=12))


@pytest.fixture
def cobalt_monthly_files(tmp_path):
    """Twelve files, one per month: the layout that needs a glob or a list.

    This is also the layout where open_mfdataset will concatenate static variables
    along time unless told not to.
    """
    d = tmp_path / 'monthly'
    d.mkdir()
    return [str(p) for p in write_cobalt(d, split_by_month=True, nmonths=12)]


@pytest.fixture
def cobalt_monthly_glob(cobalt_monthly_files):
    """A glob pattern matching the twelve per-month files."""
    return str(Path(cobalt_monthly_files[0]).parent / 'cobalt_mon*.nc')


@pytest.fixture
def fake_flood(monkeypatch):
    """Replace flood_missing with a cheap horizontal fill.

    HCtFlood is not installed here, and flooding is not what these tests are
    checking. The stand-in keeps the properties the pipeline relies on: it returns a
    DataArray with the same name, dims and coords, and it fills land horizontally one
    level at a time. Levels that are NaN everywhere stay NaN, as with flood_kara, so
    regrid_tracer's vertical fill still has something to do.
    """
    import boundary as bnd

    def _flood(arr, xdim='lon', ydim='lat', zdim=None, **kwargs):
        filled = arr
        for dim in (xdim, ydim):
            filled = filled.ffill(dim).bfill(dim)
        return filled

    monkeypatch.setattr(bnd, 'flood_missing', _flood)
    return _flood
