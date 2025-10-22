# Ocean Model Boundary Condition Processing

This repository contains a script for processing ocean boundary conditions from GLORYS (Global Ocean Physics Reanalysis) data for MOM6 ocean models. The pipeline extracts and regrids oceanographic data along model boundaries for open boundary condition forcing.

## Script Overview

### GLORYS OBC Processor (`glorys_obc_processor.py`)
Main processing engine that orchestrates boundary condition generation:
- Loads and preprocesses GLORYS reanalysis data
- Processes multiple years of data sequentially
- Applies regridding to all specified boundary segments
- Handles velocity rotation and tracer interpolation
- Manages output file organization and naming

## Requirements

- GLORYS reanalysis data files
- MOM6 horizontal grid file (`ocean_hgrid.nc`)

## Configuration File

The pipeline uses a YAML configuration file (`glorys_obc.yaml`):

```yaml
# Time range for processing
first_year: 1993
last_year: 1993

# Data paths
glorys_dir: '/path/to/glorys/*%d*nc'  # %d gets replaced with year
output_dir: '/path/to/output/directory'
hgrid: '/path/to/ocean_hgrid.nc'

# Processing options
ncrcat_years: true  # Whether to concatenate multiple years

# Boundary segments to process
segments:
  - id: 1
    border: 'south'
  - id: 2
    border: 'east'
  - id: 3
    border: 'north'
  - id: 4
    border: 'west'

# Variables to process
variables:
  - 'thetao'  # Temperature
  - 'so'      # Salinity  
  - 'zos'     # Sea surface height
  - 'uv'      # Velocity components
```

### Configuration Parameters

- **Time Range**: `first_year` and `last_year` define the processing period
- **Data Paths**: 
  - `glorys_dir`: Pattern for GLORYS files with `%d` placeholder for year
  - `output_dir`: Directory for boundary condition output files
  - `hgrid`: MOM6 horizontal grid file path
- **Segments**: List of boundary segments with ID numbers (1-4) and border locations
- **Variables**: Oceanographic variables to extract and process

## Usage

### Basic Processing
```bash
python glorys_obc_processor.py
```

The script automatically reads configuration from `glorys_obc.yaml` and processes all specified years and segments.

### Processing Steps
1. **Load Configuration**: Reads YAML parameters and validates settings
2. **Setup Segments**: Creates boundary segment objects for each specified border
3. **Process Years**: Iterates through each year in the specified range
4. **Load GLORYS Data**: Opens and preprocesses reanalysis files for each year
5. **Regrid to Boundaries**: Interpolates data onto boundary segment coordinates
6. **Write Output**: Saves boundary condition files for each segment and variable

## Output Files

The processor generates NetCDF files for each boundary segment and variable:

### File Naming Convention
- `thetao_001_{year}.nc` - Temperature for segment 1 (south border)
- `so_002_{year}.nc` - Salinity for segment 2 (east border)
- `uv_001_{year}.nc` - Velocity components for segment 1
- `zos_003_{year}.nc` - Sea surface height for segment 3 (north border)

### Output Structure
Each file contains:
- **Coordinates**: `lon_segment_XXX`, `lat_segment_XXX`, `time`
- **Dimensions**: Segment-specific dimension names (`nx_segment_XXX`, `ny_segment_XXX`)
- **Variables**: Regridded oceanographic data with proper MOM6 naming conventions
- **Metadata**: CF-compliant attributes and coordinate information

## Key Features

### Boundary Segment Handling
- **Flexible Geometry**: Supports all four model boundaries (N/S/E/W)
- **Coordinate Rotation**: Converts earth-relative to model-relative coordinates
- **Grid Mapping**: Maps irregular GLORYS grid to boundary points using xESMF

### Data Processing
- **Conservative Regridding**: Preserves oceanographic properties during interpolation
- **Missing Data Handling**: Forward/backward filling and flooding algorithms
- **Temporal Consistency**: Proper time coordinate handling across multiple years

### Velocity Processing
- **Vector Rotation**: Rotates u/v velocities from earth-relative to model grid orientation
- **Layer Thickness**: Automatically calculates vertical layer thicknesses (`dz`)
- **Boundary-Specific**: Handles different velocity staggering for each boundary type

### Tracer Processing
- **Multi-Variable Support**: Handles temperature, salinity, sea surface height
- **Depth Interpolation**: Processes 3D fields with vertical coordinate mapping
- **Quality Control**: Validates and fills missing values

## Advanced Features

### Regridding Options
- **Method Selection**: Supports nearest neighbor, bilinear, and conservative interpolation
- **Weight Reuse**: Caches regridding weights for efficiency across multiple variables
- **Periodic Boundaries**: Handles global vs regional grid configurations

### Data Quality Control
- **Flooding Algorithms**: Uses advanced algorithms for filling data over land areas
- **Boundary Filling**: Forward/backward fill along boundaries
- **Coordinate Validation**: Checks for proper angle ranges and coordinate systems

## Data Dependencies

The pipeline requires:
- **GLORYS**: Global Ocean Physics Reanalysis from Copernicus Marine Service
- **MOM6 Grid**: Horizontal grid file (`ocean_hgrid.nc`) with coordinate information
- **Proper Directory Structure**: Organized data storage for input and output files

## Output Usage

Generated boundary condition files are ready for MOM6 model integration:
- Files follow MOM6 naming conventions and NetCDF format requirements
- Time coordinates are properly configured for model forcing
- Segment-specific dimensions match MOM6 open boundary expectations
- All variables include necessary metadata and attributes for model compatibility