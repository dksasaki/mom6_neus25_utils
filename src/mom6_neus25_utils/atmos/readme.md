# ERA5 Meteorological Data Processor

Converts ERA5 GRIB files to NetCDF format for ocean model forcing, with proper time coordinates and derived variable calculations.

## Quick Start

```bash
# Process all years and pad (full pipeline)
python era5_processor.py --config era5_sfc.yaml

# Process single year only (no padding)
python era5_processor.py --config era5_sfc.yaml --year 2020

# Run only padding step (for SLURM workflows)
python era5_processor.py --config era5_sfc.yaml --pad-only

# Force overwrite existing files
python era5_processor.py --config era5_sfc.yaml --year 2020 --force
```

## Configuration

Create `era5_sfc.yaml`:

```yaml
first_year: 1993
last_year: 1994
era5_dir: '/path/to/era5/grib/files'
tmp_dir: '/path/to/temporary/files'      # Intermediate processing
output_dir: '/path/to/final/output'      # Final padded files
```

## Input Requirements

GRIB files must follow naming convention: `ERA5_{variable}_{year}.grib`

Example: `ERA5_t2m_2020.grib`

## Supported Variables

| Code | Description | Units |
|------|-------------|-------|
| `t2m` | 2-meter temperature | K |
| `d2m` | 2-meter dewpoint temperature | K |
| `u10` | 10-meter u-wind | m/s |
| `v10` | 10-meter v-wind | m/s |
| `sp` | Surface pressure | Pa |
| `msl` | Mean sea level pressure | Pa |
| `tp` | Total precipitation | m |
| `sf` | Snowfall | m |
| `ssrd` | Surface solar radiation down | J/m² |
| `strd` | Surface thermal radiation down | J/m² |

### Derived Variables (Calculated Automatically)
- **`sphum`**: Specific humidity (from `sp` and `d2m`)
- **`lp`**: Liquid precipitation (`tp` - `sf`)

## Processing Pipeline

### 1. Main Processing (`process_year`)
- Converts GRIB → NetCDF with proper time coordinates
- Handles stepped (forecast) and regular time data
- Reverses latitude if needed for consistency
- Saves to `tmp_dir`: `ERA5_{variable}_{year}.nc`

### 2. Derived Variables
- Calculates specific humidity using saturation vapor pressure
- Computes liquid precipitation with zero floor
- Saves to `tmp_dir`: `ERA5_sphum_{year}.nc`, `ERA5_lp_{year}.nc`

### 3. Padding Step (`pad_year`)
- Concatenates first timestep from next year
- Creates seamless time series for model forcing
- Saves to `output_dir`: `ERA5_{variable}_{year}_padded.nc`

## Output Structure

```
tmp_dir/
├── ERA5_t2m_1993.nc         # Processed base variables
├── ERA5_sphum_1993.nc       # Derived variables
└── ERA5_lp_1993.nc

output_dir/
├── ERA5_t2m_1993_padded.nc  # Final padded files
├── ERA5_sphum_1993_padded.nc
└── ERA5_lp_1993_padded.nc
```

## Key Features

### Time Handling
- Converts to numerical time: "days since 1990-01-01T00:00:00"
- Stacks time+step dimensions for forecast data
- Removes duplicates automatically

### Coordinate Standards
- CF-compliant attributes (axis labels, units)
- Automatic latitude reversal if needed
- Unlimited time dimension for model compatibility

### Error Handling
- Continues processing if individual variables fail
- Reports success/failure counts
- Validates input file availability

## Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--config` | YAML configuration file | `era5_processor.yaml` |
| `--year` | Process single year only | None (uses config range) |
| `--variables` | Specific variables to process | All variables |
| `--force` | Overwrite existing files | False |
| `--pad-only` | Skip processing, only pad | False |

## SLURM Workflow Example

```bash
# Step 1: Process each year in parallel
sbatch --array=1993-2020 process_year.sh

# Step 2: After all complete, run padding
python era5_processor.py --config era5_sfc.yaml --pad-only
```

## Physical Constants Used

- Molecular weight ratio (ε): 0.622
- Saturation pressure at 0°C: 611.2 Pa
- Saturation vapor pressure: Clausius-Clapeyron equation

## Dependencies

- `xarray`: NetCDF/GRIB handling
- `cfgrib`: GRIB engine for xarray
- `netCDF4`: Time coordinate conversion
- `pandas`: Datetime operations
- `numpy`: Numerical calculations
- `pyyaml`: Configuration parsing

## Notes

- Input GRIB files must have complete annual coverage
- Padding requires next year's data (except final year)
- NetCDF3_64BIT format used for padded files (MOM6 compatibility)
- Memory efficient: processes one variable at a time