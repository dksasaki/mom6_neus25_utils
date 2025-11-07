# MOM6-NEUS25 Utils

Preprocessing toolkit for MOM6-COBALT-NEUS25 regional ocean model. Generates forcing files from ERA5 (atmosphere), GLORYS (ocean BC), TPXO (tides), GloFAS (rivers), and creates nudging/damping fields.

**Note**: This is a quick reference guide. Detailed documentation for each component is available in the README files within each subdirectory (atmos/, boundary/, rivers/, sponge/). Complete model datasets (including initial conditions) required to run the model are available at [Zenodo DOI placeholder].

## Installation

[Pixi package manager](https://pixi.sh/latest/python/tutorial/) handles all dependencies including compiled libraries.

To install:
```bash
git clone https://github.com/your-org/mom6_neus25_utils
cd mom6_neus25_utils
pixi install
pixi run install-hctdflood
pixi shell  # Activate environment
```

## Quick Start

All scripts use YAML configs from `src/mom6_neus25_utils/yaml_examples/`. Configure these YAML files first - they contain paths to input data, output directories, and processing parameters.

After `pixi shell`:
```bash
# Atmosphere
era5_processor --config era5_processor.yaml

# Ocean boundaries  
glorys_obc_processor --config glorys_obc_processor.yaml
glorys_obc_processor --config glorys_obc_processor.yaml --pad-only
tpxo_obc_processor --config tpxo_obc_processor.yaml

# Rivers
glofas_processor --config glofas_processor.yaml
bgc_processor_a_runoff_ave --config bgc_processor.yaml
bgc_processor_b --config bgc_processor.yaml
bgc_processor_c_combine --config bgc_processor.yaml

# Nudging/Damping
write_nudging_data_tgb --config write_nudging.yaml
write_nudging_grid_tgb --config write_nudging.yaml
```

## Project Structure

```
mom6_neus25_utils/
├── pyproject.toml
├── pixi.toml                # Environment management
└── src/mom6_neus25_utils/
    ├── atmos/                # Atmospheric forcing processors
    ├── boundary/             # Ocean boundary condition tools
    ├── data_download/        # Reference download scripts
    ├── rivers/               # River runoff and biogeochemistry
    ├── sponge/               # Nudging and damping generators
    └── yaml_examples/        # Configuration templates
```

## Components

| Component | Input Data | Output | Purpose | Data Source DOI |
|-----------|------------|--------|---------|-----------------|
| **Atmospheric** | ERA5 GRIB (hourly) | `ERA5_{var}_{year}_padded.nc` | Surface forcing | [10.24381/cds.adbb2d47](https://doi.org/10.24381/cds.adbb2d47) |
| **Ocean BC** | GLORYS NetCDF (daily) | `{var}_{seg:03d}_{year}_padded.nc` | Lateral boundaries | [10.48670/moi-00021](https://doi.org/10.48670/moi-00021) |
| **Tidal** | TPXO9 harmonics | `t[zuv]_{seg:03d}.nc` | Tidal forcing | [TPXO Products](https://www.tpxo.net/tpxo-products-and-registration) |
| **River Runoff** | GloFAS discharge | `glofas_runoff_{year}.nc` | Freshwater input | [10.24381/cds.a4fdd6b9](https://doi.org/10.24381/cds.a4fdd6b9) |
| **River BGC** | GlobalNEWS2/USGS | `RiverNutrients_*.nc` | Nutrient concentrations | See below |
| **Nudging** | GLORYS monthly | `nudging_monthly_{year}.nc` | Interior relaxation | [10.48670/moi-00021](https://doi.org/10.48670/moi-00021) |
| **Damping** | Ocean grid | `damping_tgb_[tuv]*.nc` | Boundary sponge | N/A (generated) |

**Note**: Scripts in `data_download/` provide reference implementations for downloading ERA5, GLORYS, and GloFAS data from their respective sources.

## Data Sources and Citations

### Primary Data Sources
- **ERA5**: ECMWF Reanalysis v5 - Hourly atmospheric data (1940-present)
  - DOI: [10.24381/cds.adbb2d47](https://doi.org/10.24381/cds.adbb2d47)
  
- **GLORYS12v1**: Global Ocean Physics Reanalysis - Daily ocean fields (1993-present)
  - DOI: [10.48670/moi-00021](https://doi.org/10.48670/moi-00021)
  
- **GloFAS**: Global Flood Awareness System - River discharge (1979-present)
  - DOI: [10.24381/cds.a4fdd6b9](https://doi.org/10.24381/cds.a4fdd6b9)

### Biogeochemical Data Sources
- **GlobalNEWS2**: Global river nutrient export model
  - DOI: [10.1016/j.envsoft.2010.08.010](https://doi.org/10.1016/j.envsoft.2010.08.010)
  
- **USGS Water Quality**: Chemistry data for US rivers
  - [USGS Water Data](https://waterdata.usgs.gov/nwis)

### Tidal Data
- **TPXO9**: Global ocean tide model
  - [TPXO Products and Registration](https://www.tpxo.net/tpxo-products-and-registration)

### Initial Conditions (separate download)
- **WOA23**: World Ocean Atlas 2023 - Temperature, Salinity, Nutrients
  - Temperature: [10.25923/54bh-1613](https://doi.org/10.25923/54bh-1613)
  - Salinity: [10.25923/70qt-9574](https://doi.org/10.25923/70qt-9574)
  - Oxygen: [10.25923/rb67-ns53](https://doi.org/10.25923/rb67-ns53)
  - Nutrients: [10.25923/39qw-7j08](https://doi.org/10.25923/39qw-7j08)
  
- **ESM4 Historical**: GFDL biogeochemistry for COBALT initialization
  - DOI: [10.22033/ESGF/CMIP6.1407](https://doi.org/10.22033/ESGF/CMIP6.1407)

## Required Files

### Grid Files (must be provided)
- `ocean_hgrid.nc` - Horizontal grid specification
- `ocean_static.nc` - Static ocean fields (bathymetry, etc.)
- `ocean_mask.nc` - Land/ocean mask

### External Data Requirements
- **ERA5**: Yearly GRIB files with hourly data
- **GLORYS**: Monthly NetCDF files with daily averages  
- **GloFAS**: Yearly files with daily discharge values
- **TPXO9**: Harmonic constituent files (elevation and velocity)
- **GlobalNEWS2**: Excel file with nutrient export data
- **USGS**: Chemistry NetCDF files (optional, enhances US rivers)
- **Nudging**: Uses monthly averages derived from GLORYS

## Processing Order

1. **Independent Processing** (can run in parallel):
   - Atmospheric forcing (ERA5)
   - Ocean boundary conditions (GLORYS)
   - Tidal harmonics (TPXO9)

2. **Sequential River Processing**:
   - GloFAS discharge → Runoff climatology → River BGC nutrients

3. **Dependent Processing** (requires GLORYS completion):
   - Nudging data generation
   - Damping coefficient calculation

## SLURM Workflow Support

Many processors support a two-step workflow for HPC environments:
1. Process multiple years in parallel (without padding)
2. Apply padding in a single job using `--pad-only` flag

Example:
```bash
# Step 1: Process years in parallel (array job)
glorys_obc_processor --config glorys_obc_processor.yaml --year ${YEAR}

# Step 2: Pad all years (single job)
glorys_obc_processor --config glorys_obc_processor.yaml --pad-only
```

## Output File Conventions

- **Temporal padding**: Most files include the first timestep of the following year
- **Segment notation**: Boundary files use `_{seg:03d}_` for segment identification
- **Variable naming**: Follows MOM6 standard names where possible

## Troubleshooting

Common issues and solutions:
- **Memory errors**: Reduce chunk sizes in YAML configurations
- **Missing dependencies**: Ensure `pixi shell` is activated
- **Grid mismatch**: Verify grid files match the NEUS25 domain specification
- **Time gaps**: Check for continuous temporal coverage in input data

## License

This preprocessing toolkit is released under LGPL-3.0 license. Individual data sources maintain their original licenses as specified above.

## Acknowledgements

This toolkit builds upon:
- [CEFI Regional MOM6](https://github.com/NOAA-GFDL/CEFI-regional-MOM6) - Regional configuration framework
- [NWA25](https://github.com/jsimkins2/nwa25) - Northwest Atlantic implementation
- MOM6 Development Team at NOAA-GFDL
- COBALT biogeochemistry model developers

