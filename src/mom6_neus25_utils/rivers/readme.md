# River Processing Pipeline

Processes GloFAS discharge and river chemistry data into ocean model forcing files.

## Pipeline Overview

**Input Data:**
- GloFAS Daily Discharge (m³/s)
- USGS Chemistry (concentrations)
- GlobalNEWS2 (nutrient ratios)

**Processing Chain:**
1. GloFAS Daily → Annual Runoff Files
2. Annual Runoff → Monthly Climatology
3. Monthly Climatology + USGS Chemistry + GlobalNEWS2 → Final Nutrient Fields

## Configuration Files

### `runoff_glofas.yaml`
```yaml
# Ocean model grids
grid_mask_file: '/path/to/ocean_mask.nc'
hgrid_file: '/path/to/ocean_hgrid.nc'
ldd_file: '/path/to/ldd_glofas_v4_0.nc'

# GloFAS data
glofas_files_pattern: '/path/to/glofas_{year}*.nc'
output_dir: '/path/to/output'
TMPDIR: '/tmp'

# Processing period
start_year: 1993
end_year: 1994

# Regional subset
latitude_range: {start: 54.06, end: 24.9}
longitude_range: {start: -80.06, end: -54.9}
```

### `bgc_processor.yaml`
```yaml
# Runoff averaging
glofas_avg:
  path_source: '/path/to/glofas/output'
  path_savemat: '/path/to/averaged'
  year_first: 1993
  year_final: 1994

# BGC processing
bgc_processor:
  basin_file: '/path/to/GlobalNEWS2.xls'
  news_file: '/path/to/RiverNutrients_GlobalNEWS2.nc'
  grid_file: '/path/to/ocean_static.nc'
  runoff_file: '/path/to/glofas_runoff_mean.nc'
  output_file: '/path/to/RiverNutrients_GlobalNEWS2.nc'

bgc_processor_combine:
  filename_chem: '/path/to/mclim_chem.nc'
  filename_discharge: '/path/to/mclim_disc.nc'
  woa_temp_pattern: '/path/to/woa18_t*_04.nc'
  # Plus all fields from bgc_processor
```

## Processing Steps

### Step 1: Process GloFAS Discharge
```bash
python glofas_processor.py --config runoff_glofas.yaml
```

**Process:**
1. Loads daily GloFAS discharge (m³/s)
2. Identifies coastal pour points using LDD
3. Regrids to ocean model grid (conservative)
4. Redistributes to coastal cells
5. Converts to kg/m²/s

**Output:** `glofas_runoff_{year}.nc` with monthly data

### Step 2: Create Runoff Climatology
```bash
python bgc_processor_a_runoff_ave.py --config bgc_processor.yaml
```

**Process:**
1. Reads all yearly files
2. Groups by month
3. Averages across years

**Output:** `glofas_runoff_mean.nc` (12-month climatology)

### Step 3: Process GlobalNEWS2 Nutrients (Required)
```bash
python bgc_processor_b.py --config bgc_processor.yaml
```

**Process:**
1. Filters rivers by region and discharge (Q > 100 m³/s)
2. Maps rivers to runoff grid points by matching discharge
3. Fills gaps using nearest neighbor
4. Applies nutrient fractionation

**Output:** `RiverNutrients_GlobalNEWS2.nc` (used by Step 4)

### Step 4: Process USGS Chemistry (Optional Enhancement)
```bash
python bgc_processor_c_combine.py --config bgc_processor.yaml
```

**Process:**
1. Loads USGS chemistry and discharge data
2. Maps to grid using discharge matching
3. Gap-fills with GlobalNEWS2 ratios from Step 3
4. Creates monthly climatology with enhanced chemistry

**Output:** `RiverNutrients_Integrated.nc`

**Note:** Step 3 must be completed before Step 4 as it generates required ratio files.

## Key Algorithms

### Coastal Cell Identification
- Uses Local Drainage Direction (LDD = 5) for pour points
- Expands mask iteratively to find all coastal discharge

### Discharge Mapping
- Accumulates nearest grid cells until matching river discharge
- Maximum search radius: 2° (configurable)
- Sorts rivers by size (smallest first)

### Gap Filling
- **Primary**: Monthly USGS data
- **Secondary**: Annual USGS averages
- **Tertiary**: GlobalNEWS2 ratios × DIN
- **Spatial**: Nearest neighbor interpolation

## Output Variables

### Runoff Files
| Variable | Units | Description |
|----------|-------|-------------|
| runoff | kg/m²/s | Freshwater flux |
| area | m² | Grid cell area |
| lat/lon | degrees | Coordinates |

### Nutrient Files
| Variable | Units | Fractionation |
|----------|-------|---------------|
| NO3_CONC | mol/m³ | DIN (100%) |
| LDON_CONC | mol/m³ | DON × 0.30 |
| SLDON_CONC | mol/m³ | DON × 0.35 |
| SRDON_CONC | mol/m³ | DON × 0.35 |
| PO4_CONC | mol/m³ | DIP (100%) |
| LDOP_CONC | mol/m³ | DOP × 0.30 |
| SLDOP_CONC | mol/m³ | DOP × 0.35 |
| SRDOP_CONC | mol/m³ | DOP × 0.35 |
| NDET_CONC | mol/m³ | PN (100%) |
| PDET_CONC | mol/m³ | PP × 0.30 |
| SI_CONC | mol/m³ | Silicate |
| FED_CONC | mol/m³ | 70 µmol (constant) |
| DIC_CONC | mol/m³ | DIC (USGS only) |
| ALK_CONC | eq/m³ | Alkalinity (USGS only) |

## Dependencies

```python
# Core
xarray, numpy, pandas, scipy
# Regridding
xesmf (with ESMF library)
# I/O
netCDF4, openpyxl, pyyaml
# Climatology
cftime
```

## Important Notes

1. **Coordinate Systems**: GloFAS uses 0-360° longitude, converted to -180-180°
2. **Time Handling**: Processes year ± 1 day for boundary conditions
3. **Memory**: Process years individually for large domains
4. **Phosphorus Bioavailability**: PP multiplied by 0.3 (configurable)
5. **Iron**: Set to constant 70 µM where nutrients present
6. **USGS Corrections**: Manual coordinate fixes for major rivers

## Troubleshooting

| Issue | Solution |
|-------|----------|
| No convergence in coastal mask | Increase iteration limit in `_generate_glofas_coast_mask` |
| Rivers outside domain | Check min_dist parameter (default: 2°) |
| Missing nutrients | Verify GlobalNEWS2 file has ratio calculations |
| Memory errors | Process fewer years simultaneously |
| Regridding fails | Ensure ESMF/xesmf properly installed |
