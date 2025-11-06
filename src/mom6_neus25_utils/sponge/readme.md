# Nudging and Damping Setup

Scripts for generating ocean model boundary damping and interior nudging files.

## Quick Start

```bash
# Generate damping coefficients
python write_damping_tgb.py --config config.yaml

# Process nudging data
python write_nudging_data.py --config config.yaml
```

## Configuration

Create `config.yaml`:

```yaml
forecasts:
  first_year: 1993
  last_year: 1994

filesystem:
  ocean_hgrid: '/path/to/ocean_hgrid.nc'
  ocean_static: '/path/to/ocean_static.nc'
  monthly_data_nudging: '/path/to/glorys/*{year}*.nc'
  output_dir: '/path/to/output'
  tmp: '/tmp'
```

## Script Details

### Damping Generation (`write_damping_tgb.py`)

**Purpose**: Creates sponge layer damping at open boundaries to prevent wave reflection.

**Process**:
1. Extracts U, V, and T grids from hgrid
2. Calculates damping width based on grid resolution
3. Applies tanh profile: `1 - tanh((2/e)*(i-1)/(width-1))`
4. Combines south, east, and north boundaries

**Parameters** (hardcoded):
- Sponge width: 250 grid points
- Physical width: 20 km
- Damping rate: 1/(7 days)
- Variable width by boundary:
  - South: `width × 2 / dy`
  - East: `width × 4 / dx`
  - North: `width × 2 / dy / 2`

**Output**:
- `damping_tgb_uv_b.nc`: U/V velocity damping coefficients
- `damping_tgb_t_b.nc`: Tracer (T/S) damping coefficients

### Nudging Data Processing (`write_nudging_data.py`)

**Purpose**: Prepares monthly GLORYS data for interior model nudging/relaxation.

**Process**:
1. Loads GLORYS monthly averages
2. Forward-fills missing values in depth
3. Regrids to model grid (nearest neighbor)
4. Adds time bounds for monthly data
5. Backward-fills remaining gaps horizontally

**Variables Processed**:
| Variable | Grid Type | Dimensions | Description |
|----------|-----------|------------|-------------|
| thetao | Tracer | (yh, xh) | Temperature |
| so | Tracer | (yh, xh) | Salinity |
| uo | U-grid | (yh, xq) | Zonal velocity |
| vo | V-grid | (yq, xh) | Meridional velocity |

**Time Bounds**:
- Start: First day of month at 00:00
- End: Last day at 23:59:59 (except December → Jan 1 00:00)

**Output**: `nudging/nudging_monthly_{year}.nc`

## File Formats

### Damping Files
```netcdf
dimensions:
  xh, xq = longitude points
  yh, yq = latitude points
variables:
  Idamp_u(yh, xq) - U damping [s⁻¹]
  Idamp_v(yq, xh) - V damping [s⁻¹]
  Idamp(yh, xh)   - T/S damping [s⁻¹]
```

### Nudging Files
```netcdf
dimensions:
  time = unlimited
  depth = vertical levels
  xh, yh = horizontal grid
variables:
  thetao(time, depth, yh, xh) - Temperature [°C]
  so(time, depth, yh, xh)     - Salinity [PSU]
  uo(time, depth, yh, xq)     - U velocity [m/s]
  vo(time, depth, yq, xh)     - V velocity [m/s]
```

## Grid Staggering

MOM6 uses Arakawa C-grid:
- Tracers (T/S): Cell centers (xh, yh)
- U velocity: West/East faces (xq, yh)
- V velocity: South/North faces (xh, yq)

## Usage in MOM6

### Damping (MOM_input)
```
SPONGE = True
SPONGE_UV = True
SPONGE_DAMPING_FILE = "damping_tgb_t_b.nc"
SPONGE_IDAMP_VAR = "Idamp"
SPONGE_UV_DAMPING_FILE = "damping_tgb_uv_b.nc"
SPONGE_IDAMP_U_VAR = "Idamp_u"
SPONGE_IDAMP_V_VAR = "Idamp_v"
```

### Nudging (data_table)
```
"OCN", "runoff", "runoff", "./nudging/nudging_monthly_%4yr.nc", "none", 1.0
```

## Notes

- GLORYS data must be monthly averages (not daily)
- Damping strongest at boundaries, decays inward
- Nudging typically applied weakly in interior (~30-90 day timescale)
- Both scripts handle missing data via filling strategies
- Output format: NETCDF3_64BIT for compatibility