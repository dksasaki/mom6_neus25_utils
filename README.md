# MOM6-NEUS25 Utils

Preprocessing toolkit for MOM6-COBALT-NEUS25 regional ocean model. Generates forcing files from ERA5 (atmosphere), GLORYS (ocean BC), TPXO (tides), GloFAS (rivers), and creates nudging/damping fields.

**Note**: This is a quick reference guide. Detailed documentation for each component is available in the README files within each subdirectory (atmos/, boundary/, rivers/, sponge/). Other datasets (including initial conditions) required to run the model (but not these scripts) are available here.

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

All scripts use YAML configs from src/mom6_neus25_utils/yaml_examples/. Configure these YAML files first - they contain paths to input data, output directories, and processing parameters.

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

## Project structure

```
mom6_neus25_utils/
├── pyproject.toml
└── src/mom6_neus25_utils/
    ├── atmos/                # Atmospheric forcing
    ├── boundary/             # Ocean boundary conditions  
    ├── data_download/        # Reference download scripts
    ├── rivers/               # River runoff and biogeochemistry
    ├── sponge/               # Nudging and damping
    └── yaml_examples/        # Configuration templates
```

## Components

| Component | Input Data | Output | Purpose |
|-----------|------------|--------|---------|
| **Atmospheric** | ERA5 GRIB (hourly) | `ERA5_{var}_{year}_padded.nc` | Surface forcing |
| **Ocean BC** | GLORYS NetCDF (daily) | `{var}_{seg:03d}_{year}_padded.nc` | Lateral boundaries |
| **Tidal** | TPXO9 harmonics | `t[zuv]_{seg:03d}.nc` | Tidal forcing |
| **River Runoff** | GloFAS discharge | `glofas_runoff_{year}.nc` | Freshwater input |
| **River BGC** | GlobalNEWS2/USGS | `RiverNutrients_*.nc` | Nutrient concentrations |
| **Nudging** | GLORYS monthly | `nudging_monthly_{year}.nc` | Interior relaxation |
| **Damping** | Ocean grid | `damping_tgb_[tuv]*.nc` | Boundary sponge |

**Note**: Scripts in data_download/ provide reference implementations for downloading ERA5, GLORYS, and GloFAS data from their respective sources.

## Required Files

**Grid**: `ocean_hgrid.nc`, `ocean_static.nc`, `ocean_mask.nc`

**External Data**:
- ERA5: Yearly GRIB, hourly data
- GLORYS: Monthly NetCDF, daily averages  
- GloFAS: Yearly files, daily discharge
- TPXO9: Harmonic constituents
- GlobalNEWS2: Excel nutrient data
- USGS: Chemistry NetCDF (optional)
- nudged data used monthly averages of Glorys

## Processing Order

1. **Independent**: Atmospheric, Ocean BC, Tidal
2. **Sequential**: GloFAS → Runoff climatology → River BGC
3. **Dependent**: Nudging/Damping (needs GLORYS)

## Parallel Processing

```bash
#SBATCH --array=1993-2020
glorys_obc_processor --year ${SLURM_ARRAY_TASK_ID}
# Run padding after all complete
```

## Troubleshooting

- **Memory errors**: Reduce year range
- **Missing nutrients**: Check GlobalNEWS2 completed
- **Regridding fails**: Verify xesmf installation
- **Rivers misplaced**: Check coordinate corrections

## Acknowledgements

Based on [CEFI Regional MOM6](https://github.com/NOAA-GFDL/CEFI-regional-MOM6) and [NWA25](https://github.com/jsimkins2/nwa25).