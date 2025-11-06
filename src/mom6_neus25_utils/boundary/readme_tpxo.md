# Tidal OBC Processor

Generates tidal boundary conditions from TPXO9 data for MOM6 ocean models.

## Quick Start

```bash
python tidal_obc_processor.py --config tides_obc.yaml
```

## Configuration

Create `tides_obc.yaml`:

```yaml
# TPXO data directory
tpxo_dir: '/path/to/tpxo9/data'

# MOM6 grid and output
hgrid_file: '/path/to/ocean_hgrid.nc'
output_dir: '/path/to/output'

# Tidal constituents (indices 0-9 in TPXO)
constituents: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

# Start date (1 month before simulation)
start_date: '1992-12-01'

# Optional: Spatial subset
horizontal_subset:
  lon: [-180, 180]
  lat: [-90, 90]

# Boundary segments
segments:
  - id: 1
    border: 'south'
  - id: 2
    border: 'east'
  - id: 3
    border: 'north'
  - id: 4
    border: 'west'
```

## Input Requirements

TPXO9 v5a files in `tpxo_dir`:
- `h_tpxo9.v5a.nc` - Tidal elevation harmonics
- `u_tpxo9.v5a.nc` - Velocity harmonics (u and v)

## Processing Steps

1. **Load TPXO Data**
   - Elevation: Amplitude (`ha`) and phase (`hp`)
   - Velocity: u/v amplitudes (`ua`,`va`) and phases (`up`,`vp`)
   - Converts cm/s → m/s for velocities

2. **Convert to Complex Form**
   - Elevation: `h = ha * exp(-i * hp)`
   - Velocity: `u = ua * exp(-i * up)`, `v = va * exp(-i * vp)`
   - Splits into real (`hRe`,`uRe`,`vRe`) and imaginary (`hIm`,`uIm`,`vIm`) parts

3. **Regrid to Boundaries**
   - Uses `Segment` class for each boundary
   - Applies flood algorithm to fill land values
   - Creates single time point (start_date)

## Output Files

Generated in `output_dir` via Segment class:
- `tz_{segment_id:03d}.nc` - Tidal elevation harmonics
- `tu_{segment_id:03d}.nc` - Tidal u-velocity harmonics  
- `tv_{segment_id:03d}.nc` - Tidal v-velocity harmonics

Each file contains real and imaginary components for selected constituents.

## Tidal Constituents

TPXO9 constituent indices (0-9):
| Index | Name | Period (hours) |
|-------|------|----------------|
| 0 | M2 | 12.42 |
| 1 | S2 | 12.00 |
| 2 | N2 | 12.66 |
| 3 | K2 | 11.97 |
| 4 | K1 | 23.93 |
| 5 | O1 | 25.82 |
| 6 | P1 | 24.07 |
| 7 | Q1 | 26.87 |
| 8 | MM | 661.31 |
| 9 | MF | 327.86 |

## Key Processing Details

- **Complex Arithmetic**: Converts amplitude/phase to real/imaginary for MOM6
- **Unit Conversion**: TPXO velocities scaled from cm/s to m/s
- **Grid Alignment**: Different grids for h (z-points) and u/v (u/v-points)
- **Flooding**: Fills missing ocean values near boundaries

## Dependencies

- `xarray`: NetCDF handling
- `boundary`: Custom Segment class for regridding
- `numpy`: Complex number operations
- `pandas`: Date handling

## Notes

- Start date should be ~1 month before simulation start
- All processing handled by `Segment` class methods
- Output format compatible with MOM6 tidal forcing requirements