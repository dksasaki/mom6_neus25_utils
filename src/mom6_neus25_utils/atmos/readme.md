# ERA5 Meteorological Data Processing

This repository contains a script for processing ERA5 reanalysis meteorological data for ocean model forcing. The pipeline converts GRIB format ERA5 data into NetCDF files with proper time coordinates and calculates derived meteorological variables.

## Script Overview

### ERA5 Processor (`era5_processor.py`)
Main processing engine for ERA5 surface meteorological data:
- Loads and processes ERA5 GRIB files for multiple variables and years
- Handles time coordinate conversion from GRIB step format to standard time
- Calculates derived variables (specific humidity, liquid precipitation)
- Manages temporal padding by concatenating adjacent years
- Outputs CF-compliant NetCDF files ready for ocean model forcing


## Configuration File

The pipeline uses a YAML configuration file (`era5_sfc.yaml`):

```yaml
# Time range for processing
first_year: 1993
last_year: 1994

# Data paths
era5_dir: '/path/to/era5/grib/files'
output_dir: '/path/to/output/directory'
```

### Configuration Parameters

- **Time Range**: `first_year` and `last_year` define the processing period
- **Data Paths**: 
  - `era5_dir`: Directory containing ERA5 GRIB files
  - `output_dir`: Directory for processed NetCDF output files

### Expected Input File Format

ERA5 GRIB files should follow the naming convention:
- `ERA5_{variable}_{year}.grib`

For example:
- `ERA5_t2m_1993.grib` (2-meter temperature for 1993)
- `ERA5_u10_1993.grib` (10-meter u-wind for 1993)

## Usage

### Basic Processing (All Years and Variables)
```bash
python era5_processor.py --config era5_sfc.yaml
```

### Process Single Year
```bash
python era5_processor.py --config era5_sfc.yaml --year 1993
```

### Process Specific Variables
```bash
python era5_processor.py --config era5_sfc.yaml --variables t2m u10 v10
```

### Force Overwrite Existing Files
```bash
python era5_processor.py --config era5_sfc.yaml --force
```

### Command Line Options
- `--config`: Path to YAML configuration file (default: `era5_sfc.yaml`)
- `--year`: Process single year instead of range from config
- `--variables`: List of specific variables to process (default: all)
- `--force`: Force overwrite of existing output files

## Supported Variables

The processor handles the following ERA5 surface variables:

| Variable | Description | ERA5 Code |
|----------|-------------|-----------|
| `t2m` | 2-meter temperature | 2m_temperature |
| `d2m` | 2-meter dewpoint temperature | 2m_dewpoint_temperature |
| `u10` | 10-meter u-component of wind | 10m_u_component_of_wind |
| `v10` | 10-meter v-component of wind | 10m_v_component_of_wind |
| `sp` | Surface pressure | surface_pressure |
| `msl` | Mean sea level pressure | mean_sea_level_pressure |
| `tp` | Total precipitation | total_precipitation |
| `sf` | Snowfall | snowfall |
| `ssrd` | Surface solar radiation downwards | surface_solar_radiation_downwards |
| `strd` | Surface thermal radiation downwards | surface_thermal_radiation_downwards |

## Processing Steps

1. **Load Configuration**: Reads YAML parameters and determines years to process
2. **Process Base Variables**: Converts each GRIB file to NetCDF with proper time coordinates
3. **Time Coordinate Handling**: 
   - Handles both regular time series and stepped forecast data
   - Concatenates with next year's data for temporal padding
   - Converts to standard time units
4. **Coordinate Adjustments**: Ensures proper latitude/longitude ordering and attributes
5. **Calculate Derived Variables**: Computes specific humidity and liquid precipitation
6. **Quality Control**: Validates data and handles missing values

## Output Files

### Base Variables
For each processed variable and year:
- `ERA5_{variable}_{year}_padded.nc` - Processed NetCDF file

Examples:
- `ERA5_t2m_1993_padded.nc` - 2-meter temperature for 1993
- `ERA5_u10_1993_padded.nc` - 10-meter u-wind for 1993

### Derived Variables
- `ERA5_sphum_{year}.nc` - Specific humidity calculated from surface pressure and dewpoint
- `ERA5_lp_{year}.nc` - Liquid precipitation (total precipitation minus snowfall)

## Key Features

### Time Coordinate Processing
- **Stepped Data Handling**: Properly processes ERA5 forecast step data by combining time and step dimensions
- **Temporal Padding**: Concatenates with adjacent year data for seamless time series
- **Standard Time Units**: Converts to "days since 1990-01-01T00:00:00" format
- **Duplicate Removal**: Automatically removes duplicate time steps

### Data Quality Control
- **Missing Value Handling**: Applies appropriate thresholds for missing data detection
- **Coordinate Validation**: Ensures proper latitude/longitude ordering
- **CF Compliance**: Adds standard coordinate attributes for ocean model compatibility

### Derived Variable Calculations

#### Specific Humidity
Calculated using the standard meteorological formula:
- Input: Surface pressure (`sp`) and 2-meter dewpoint temperature (`d2m`)
- Formula: Based on saturation vapor pressure and mixing ratio calculations
- Output: Specific humidity in kg/kg

#### Liquid Precipitation
Simple difference calculation:
- Formula: Total precipitation (`tp`) minus snowfall (`sf`)
- Quality control: Sets negative values to zero
- Output: Liquid precipitation rate

### Memory Management
- **Efficient Processing**: Processes one variable at a time to minimize memory usage
- **Resource Cleanup**: Properly closes datasets after processing
- **Error Handling**: Continues processing other variables if one fails

## Output Format

All output files are CF-compliant NetCDF files with:
- **Time Dimension**: Unlimited dimension for model compatibility
- **Coordinate Attributes**: Proper axis and units metadata
- **Standard Variables**: Following CF naming conventions where possible
- **Temporal Coverage**: Full year plus padding for smooth model integration

## Data Dependencies

The pipeline requires:
- **ERA5 Reanalysis**: Surface meteorological fields from Copernicus Climate Data Store
- **GRIB Format**: Input data must be in GRIB format as downloaded from ERA5
- **Complete Time Series**: Full annual coverage for proper temporal padding

## Integration with Ocean Models

Generated files are ready for ocean model forcing:
- Standard time coordinates compatible with MOM6 and other ocean models
- CF-compliant metadata and variable naming
- Proper spatial coordinate ordering
- Complete set of surface meteorological variables for air-sea flux calculations