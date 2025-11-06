# GLORYS OBC Processor

Generates ocean model boundary conditions from GLORYS reanalysis data for MOM6.

## Quick Start

```bash
# Process all years and pad
python glorys_obc_processor.py --config glorys_obc.yaml

# Process single year (development/testing)
python glorys_obc_processor.py --config glorys_obc.yaml --year 2020
```

## Configuration

Create `glorys_obc.yaml`:

```yaml
first_year: 1993
last_year: 1994
glorys_dir: '/path/to/glorys/*%d*.nc'  # %d = year placeholder
tmp_dir: '/path/to/intermediate/files'  # Intermediate processing
output_dir: '/path/to/final/output'     # Final padded files
hgrid: '/path/to/ocean_hgrid.nc'       # MOM6 grid file

segments:
  - id: 1
    border: 'south'
  - id: 2
    border: 'east'
  - id: 3
    border: 'north'
  - id: 4
    border: 'west'

variables:
  - 'thetao'  # Temperature
  - 'so'      # Salinity  
  - 'zos'     # Sea surface height
  - 'uv'      # Velocity (u,v components)
```

## Input Requirements

- **GLORYS files**: Pattern matching `glorys_dir` with year substitution
- **MOM6 grid**: `ocean_hgrid.nc` with coordinate/angle information
- **Variable names**: GLORYS standard names (thetao, so, uo, vo, zos)

## Processing Pipeline

### Step 1: Process Raw Data (`run()`)
1. Loads GLORYS NetCDF files matching pattern
2. Renames coordinates: `latitude→lat`, `longitude→lon`, `depth→z`
3. Adjusts first year's initial time to midnight
4. Regrids to boundary segments using xESMF
5. Saves to `tmp_dir`: `{variable}_{segment_id:03d}_{year}.nc`

### Step 2: Pad Files (`pad_year()`)
1. Reads files from Step 1
2. Concatenates first timestep from next year
3. Saves to `output_dir`: `{variable}_{segment_id:03d}_{year}_padded.nc`

**Note**: When using SLURM, run Step 1 in parallel, then Step 2 after all complete.

## Output Structure

```
tmp_dir/
├── thetao_001_1993.nc    # Temperature, south boundary
├── so_001_1993.nc        # Salinity, south boundary
├── uv_001_1993.nc        # Velocity, south boundary
└── zos_001_1993.nc       # SSH, south boundary

output_dir/
├── thetao_001_1993_padded.nc
├── so_001_1993_padded.nc
├── uv_001_1993_padded.nc
└── zos_001_1993_padded.nc
```

## Segment Numbering

Standard convention (configurable):
- 001: South boundary
- 002: East boundary
- 003: North boundary
- 004: West boundary

## Key Processing Details

### Velocity Handling
- Processes `uo` and `vo` together as `uv`
- Rotates from earth-relative to model grid coordinates
- Uses segment-specific coordinate rotation from hgrid

### Tracer Variables
- Direct regridding for thetao, so, zos
- No rotation needed (scalar quantities)
- Fills missing values using flood algorithm

### Time Adjustments
- First year: Floors initial time to midnight
- Padding: Adds next year's first timestep for continuity
- Format: NetCDF3_64BIT for MOM6 compatibility

## Dependencies

- `xarray`: NetCDF data handling
- `xESMF`: Regridding operations (via Segment class)
- `boundary`: Custom Segment class for boundary processing
- `pyyaml`: Configuration parsing

## Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--config` | YAML configuration file | `glorys_processor.yaml` |
| `--year` | Single year to process (testing) | None (uses config range) |

## SLURM Workflow

```bash
# Step 1: Process years in parallel
#!/bin/bash
#SBATCH --array=1993-2020
python glorys_obc_processor.py --config glorys_obc.yaml --year ${SLURM_ARRAY_TASK_ID}

# Step 2: After all complete, pad files
python glorys_obc_processor.py --config glorys_obc.yaml  # Uncomment pad_year() in main
```

## Important Notes

- The `Segment` class (from `boundary.py`) handles actual regridding
- Code currently has Step 2 (padding) commented out in `main()`
- Uncomment `glo.pad_year()` when ready to pad after processing
- Missing next year data: Padding still creates output file without concatenation
- All variables listed in config must exist in GLORYS files