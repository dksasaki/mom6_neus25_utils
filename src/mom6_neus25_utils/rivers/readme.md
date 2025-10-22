# Ocean Model River Processing Pipeline

This repository contains scripts for processing river runoff and nutrient data for ocean model forcing. The pipeline converts GloFAS (Global Flood Awareness System) discharge data and USGS river chemistry data into model-ready formats.

## Scripts Overview

### 1. GloFAS Runoff Processing (`glofas_processor.py`)
Processes daily GloFAS discharge data into model-compatible runoff fields:
- Converts m³/s discharge to kg/m²/s runoff using grid cell areas
- Maps river discharge to coastal ocean grid points
- Redistributes runoff to nearest coastal cells

### 2. Multi-year Averaging (`bgc_processor_a_runoff_ave.py`)  
Calculates multi-year monthly climatologies from processed GloFAS data:
- Reads multiple years of GloFAS runoff data
- Computes monthly averages across all years
- Saves climatological monthly means

### 3. River Biogeochemistry Processing (`bgc_processor_b.py`)
Creates nutrient concentration fields from USGS river chemistry data:
- Maps USGS station data to model grid using discharge matching
- Applies GlobalNEWS2 ratios for gap filling missing nutrients
- Uses WOA temperature for oxygen saturation calculations

### 4. Combined BGC Processing (`bgc_processor_c_combine.py`)
Integrated version combining all biogeochemistry processing steps:
- Unified workflow from USGS data to final nutrient fields
- Enhanced coordinate correction and quality control
- Monthly climatology generation with proper fractionation

## Requirements

- GloFAS discharge data files
- USGS river chemistry and discharge data
- WOA temperature climatology
- Ocean model grid files (hgrid, static, mask)

## Configuration Files

The pipeline uses two YAML configuration files:

### `runoff_glofas.yaml` - GloFAS Processing Configuration
```yaml
# Grid files for ocean model
grid_mask_file: "/path/to/ocean_mask.nc"
hgrid_file: "/path/to/ocean_hgrid.nc" 
grid_mosaic_file: "/path/to/ocean_mosaic.nc"

# GloFAS input and output
glofas_files_pattern: "/path/to/glofas_*.nc"
ldd_file: "/path/to/ldd_glofas_v4_0.nc"  # Land drainage direction
output_dir: "/path/to/output/directory"
TMPDIR: "/path/to/tmp/directory"

# Processing parameters
start_year: 1993
end_year: 1994

# Regional bounds for LDD file subsetting
latitude_range:
  start: 54.06
  end: 24.9
longitude_range:
  start: -80.06
  end: -54.9
```

### `bgc_processor.yaml` - BGC Processing Configuration
```yaml
# Multi-year averaging configuration
glofas_avg:
  path_source: "/path/to/glofas/output"    # From glofas_processor.py
  path_savemat: "/path/to/averaged/output"
  year_first: 1993
  year_final: 1994

# BGC processing configuration  
bgc_processor_combine:
  # USGS river data
  filename_chem: "/path/to/mclim_chem.nc"
  filename_discharge: "/path/to/mclim_disc.nc"
  
  # Reference datasets
  basin_file: "/path/to/GlobalNEWS2_dataset.xls"
  news_file: "/path/to/RiverNutrients_GlobalNEWS2.nc"
  grid_file: "/path/to/ocean_static.nc"
  woa_temp_pattern: "/path/to/woa18_decav_t*_04.nc"
  
  # Processed runoff (from glofas_processor.py)
  runoff_file: "/path/to/glofas_runoff_mean.nc"
  
  # Final output
  output_file: "/path/to/final_river_nutrients.nc"
```

## Usage

### Step 1: Process GloFAS Discharge Data
```bash
python glofas_processor.py
```
This reads the configuration from `runoff_glofas.yaml` and processes daily GloFAS data into annual runoff files.

### Step 2: Create Multi-year Averages
```bash
python glofas_processor.py
```
Computes monthly climatological means from the processed GloFAS data using `bgc_processor.yaml` configuration.

### Step 3: Process River Biogeochemistry
```bash
python bgc_processor_combine.py
```
Creates final nutrient concentration fields mapped to the model grid using the `bgc_processor_combine` section of `bgc_processor.yaml`.

## Output Files

### GloFAS Processing
- `glofas_runoff_{year}.nc` - Annual runoff files with monthly resolution
  - Variables: runoff (kg/m²/s), area, lat, lon coordinates
  - Runoff redistributed to coastal ocean cells

### Multi-year Averaging  
- `glofas_runoff_mean.nc` - Monthly climatological runoff
  - 12-month climatology averaged across all processed years
  - Used as input for biogeochemistry processing

### BGC Processing
- `RiverNutrients_Integrated_*.nc` - Complete nutrient climatology
  - Monthly concentrations for: DIC, ALK, NO3, NH4, dissolved/particulate organic nutrients
  - Fractionated organic matter (labile, semi-labile, refractory components)
  - Iron, silicate, and oxygen concentrations
  - Mapped to model grid with proper units (mol/m³)

## Key Features

- **Conservative regridding**: Preserves total discharge when mapping between grids
- **Coastal redistribution**: Ensures runoff reaches ocean cells rather than land
- **Gap filling**: Uses GlobalNEWS2 ratios for missing nutrient data
- **Quality control**: Coordinate corrections and data validation
- **Climatological output**: Monthly climatologies suitable for model forcing
- **Biogeochemical completeness**: All major nutrients required for ecosystem models

## Data Dependencies

The pipeline requires several external datasets:
- **GloFAS**: Global river discharge reanalysis
- **USGS**: River chemistry measurements from US monitoring stations  
- **GlobalNEWS2**: Global nutrient export model ratios
- **WOA**: World Ocean Atlas temperature climatology
- **Ocean model grids**: Horizontal grid, land mask, and static files