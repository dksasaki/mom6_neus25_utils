# Ocean Model Preprocessing Scripts

This repository contains scripts for generating ocean model boundary conditions and nudging data.

## Scripts

### 1. Damping File Generation (`write_damping_tgb.py`)
Generates damping files for ocean model boundaries using a tanh-based sponge layer approach.

### 2. Nudging Data Processing (`write_nudging_data.py`)
Processes monthly temperature, salinity, and velocity data for model nudging by regridding GLORYS monthly averages to the model grid.

## Requirements

- Ocean grid file (`ocean_hgrid.nc`)
- Ocean static file (`ocean_static.nc`)
- GLORYS monthly average data files
- Configuration file (YAML format)

## Configuration File Format

Create a YAML configuration file with the following structure:

```yaml
forecasts:
  first_year: 1993
  last_year: 1994 

filesystem:
  # Path to ocean_static file without mask table:
  ocean_static: "/path/to/ocean_static.nc"
  
  # Path of ocean_hgrid.nc (required for damping generation)
  ocean_hgrid: "/path/to/ocean_hgrid.nc"
  
  # Path to monthly averaged temperature and salinity data, will be used for nudging
  # (ATTENTION {year} will be replaced in the script)
  monthly_data_nudging: "/path/to/glorys/{year}*.nc"
  
  # Where to put data that will be used to force the model:
  output_dir: "/path/to/output/directory"
  
  tmp: "/path/to/tmp/directory"
```

**Note**: 
- For damping generation: only `ocean_hgrid` and `output_dir` are required
- For nudging data: `ocean_static`, `monthly_data_nudging`, `output_dir`, and `tmp` are required
- The `{year}` placeholder in `monthly_data_nudging` gets replaced with actual years from the forecast range
- GLORYS data must be monthly averages for proper nudging data generation

## Usage

### Damping Files
1. Ensure your configuration file includes the required paths (see format above)

2. Run the damping script:
   ```bash
   python write_nudging_grid_tgb.py -c write_nudging.yaml
   ```

### Nudging Data
1. Ensure GLORYS monthly average data files are available in the specified directory

2. Run the nudging script:
   ```bash
   python write_nudging_data.py -c write_nudging.yaml
   ```

## Output Files

### Damping Script
- `damping_tgb_uv_b.nc` - Damping coefficients for u and v velocity components
- `damping_tgb_t_b.nc` - Damping coefficients for tracers (temperature/salinity)

### Nudging Script
- `nudging/nudging_monthly_{year}.nc` - Monthly nudging data for each forecast year
  - Contains regridded temperature (`thetao`), salinity (`so`), u-velocity (`uo`), and v-velocity (`vo`)
  - Time bounds added for proper model integration
  - Data forward-filled in depth and backward-filled spatially

## Parameters

### Damping Script
The script uses these hardcoded parameters:
- `nsponge`: 250 grid points (sponge layer width)
- `width`: 20 km (physical width of damping zone)  
- `rate`: 1/(7 days) (damping timescale)

Damping is applied to the southern, eastern, and northern boundaries with varying intensities based on local grid spacing.

### Nudging Script
- Uses nearest neighbor regridding (`nearest_s2d`) from GLORYS to model grid
- Processes temperature, salinity, and velocity components
- Creates time bounds for monthly data extending to end of each month
- Forward fills missing values in depth dimension