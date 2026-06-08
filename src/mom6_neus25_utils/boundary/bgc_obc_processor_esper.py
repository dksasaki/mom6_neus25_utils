#!/usr/bin/env python3
"""
Ocean Nutrient Data Processing Pipeline
======================================

This script processes ocean temperature and salinity data to predict nutrient 
concentrations using neural network models from the PyESPER package.

Overview
--------
The pipeline performs the following steps:
1. Loads ocean salinity and temperature data for a specified year
2. Applies temporal averaging (monthly or seasonal)
3. Uses neural networks to predict nutrient concentrations based on T/S data
4. Saves the results to NetCDF format

Usage
-----
    python script.py [year] [nsegment] [average_type] [nutrients...]
    
Arguments:
    year         : Year of data to process (e.g., 1993)
    nsegment     : Segment identifier (e.g., '002')
    average_type : Temporal averaging method ('resample' or 'groupby')
    nutrients    : Space-separated list of nutrients to predict
                   Options: nitrate, silicate, phosphate
                   Default: nitrate (if none specified)

Examples:
    # Process nitrate only with monthly resampling
    python script.py 1993 002 resample
    
    # Process multiple nutrients with seasonal grouping
    python script.py 1993 002 groupby nitrate silicate phosphate

Required Data Structure:
    data/esper/
        ├── so_*.nc          # Salinity files
        ├── thetao_*.nc      # Temperature files
        ├── ocean_topog.nc   # Topography
        └── bgc_woa_*.nc     # World Ocean Atlas data
    
Output:
    data/output/nutrients_[segment]_[year].nc

Dependencies:
    - xarray
    - numpy
    - PyESPER (custom package for neural network predictions)

Notes:
    - The script expects specific file naming conventions for input data
    - Neural network uses Equation 8 from the ESPER methodology
    - Output uses NETCDF3_64BIT format for compatibility


"""

import sys
import xarray as xr
import glob
import numpy as np
import os.path as osp
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
import yaml

from PyESPER.nn import nn


@dataclass
class Config:
    """Configuration for nutrient processing."""
    year: int
    average_type: str
    nutrients: List[str]
    TS_bound_files_path: str
    output_path: str
    pyesper_path: str
    topog_file: str
    clim_type: str = 'season'
    nutrient_mapping: Dict[str, str] = None

    
    # Nutrient mapping
    nutrient_mapping: Dict[str, str] = None
    
    def __post_init__(self):
        """Initialize nutrient mapping after object creation."""
        if self.nutrient_mapping is None:
            self.nutrient_mapping = {
                'nitrate': 'no3',
                'silicate': 'sio4',
                'phosphate': 'po4',
                'DIC': 'dissic',
                'TA':  'talk'
            }
        
        # Validate configuration
        self._validate()

    def _validate(self):
        for nutrient in self.nutrients:
            if nutrient not in self.nutrient_mapping:
                raise ValueError(f"'{nutrient}' is not a valid nutrient. "
                               f"Valid options: {list(self.nutrient_mapping.keys())}")
    
    @classmethod
    def from_args(cls):
        """Create Config from command line arguments."""
        if len(sys.argv) < 4:
            print(__doc__)
            sys.exit(1)
        
        year = int(sys.argv[1])
        nsegment = sys.argv[2]
        average_type = sys.argv[3]
        nutrients = sys.argv[4:] if len(sys.argv) > 3 else ['nitrate']
        
        return cls(year=year, nsegment=nsegment, 
                  average_type=average_type, nutrients=nutrients)

    @classmethod
    def from_script(cls, year: int,
                         nsegment: str,
                         average_type: str,
                         nutrients: list):
        """Create Config from script arguments."""
        
        return cls(year=year, nsegment=nsegment, 
                  average_type=average_type, nutrients=nutrients)

    @classmethod
    def from_yaml(cls, config_file, config_key='esper'):
        """Create Config from a YAML configuration file."""
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        if config_key is not None:
            return cls(**config[config_key])
        else:
            return cls(**config)


class OceanDataLoader:
    """Handle loading and preprocessing of ocean data."""
    
    def __init__(self, config: Config):
        self.config = config
    
    def load_pattern(self, file_pattern: str, segment: int) -> xr.Dataset:
        """Load ocean data files matching pattern."""
        nseg = segment
        files = sorted(glob.glob(
            osp.join(self.config['TS_bound_files_path'], f"{file_pattern}*{nseg:03d}*{self.config['year']}*nc")
        ))
        
        if not files:
            raise FileNotFoundError(
                f"No files found for pattern: {file_pattern}*{nseg:03d}*{self.config['year']}*nc"
            )
        
        datasets = []
        for filepath in files:
            print(f"Loading: {filepath}")
            ds = xr.open_dataset(filepath)
            ds = ds.set_coords([
                f'lon_segment_{nseg:03d}',
                f'lat_segment_{nseg:03d}'
            ])
            datasets.append(ds)
        
        return xr.concat(datasets, dim='time')
    
    def load_all_data(self, segment) -> Tuple[xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset, xr.DataArray]:
        """Load all required datasets."""
        nseg = segment

        print(f"\nLoading data for year {self.config['year']}, segment {nseg:03d}")
        
        # Load ocean data
        ds_salt = self.load_pattern('so', nseg)
        ds_temp = self.load_pattern('thetao', nseg)
        ds = xr.merge([ds_salt, ds_temp])
        
        
        # Calculate depth coordinate
        z = ds[f'dz_thetao_segment_{nseg:03d}'].cumsum(
            dim=f'nz_segment_{nseg:03d}'
        )

        ds = ds.resample(time='1MS').mean(dim='time')
        ds.load()

        dsout = self._create_output_template(ds, nseg)

        
        return ds, dsout, z
    
    
    def _create_output_template(self, dsmonthly: xr.Dataset, nseg: int) -> xr.Dataset:
        """Create output dataset template with nutrient variables."""
        
        # Start with base template including salinity variable
        template_vars = [
            f'so_segment_{nseg:03d}',
            f'lon_segment_{nseg:03d}',
            f'lat_segment_{nseg:03d}',
            f'ny_segment_{nseg:03d}',
            f'dz_so_segment_{nseg:03d}'
        ]
        
        dsout = dsmonthly[template_vars].copy(deep=True)
        
        # Create nutrient variables from salinity template
        # First nutrient takes over the salinity variable
        if self.config['nutrients']:
            first_nutrient = self.config['nutrients'][0]
            first_var_name = self.config['nutrient_mapping'][first_nutrient]
            dsout = dsout.rename({f'so_segment_{nseg:03d}': f'{first_var_name}_segment_{nseg:03d}'})
            dsout[f'{first_var_name}_segment_{nseg:03d}'].values[:] = 0

            dsout = dsout.rename({f'dz_so_segment_{nseg:03d}': f'dz_{first_var_name}_segment_{nseg:03d}'})
            
            # Additional nutrients need to be copied
            for nutrient in self.config['nutrients'][1:]:
                var_name = self.config['nutrient_mapping'][nutrient]
                # Copy from the first nutrient variable (as template)
                dsout[f'{var_name}_segment_{nseg:03d}'] = dsout[f'{first_var_name}_segment_{nseg:03d}'].copy(deep=True)
                dsout[f'{var_name}_segment_{nseg:03d}'].values[:] = 0
                dsout[f'dz_{var_name}_segment_{nseg:03d}'] = dsout[f'dz_{first_var_name}_segment_{nseg:03d}'].copy(deep=True)
        
        return dsout


class NeuralNetworkPredictor:
    """Handle neural network predictions for nutrients."""
    
    def __init__(self, config: Config):
        self.config = config
    
    def predict(self, dsmonthly: xr.Dataset, dsout: xr.Dataset, nseg:int, nutrientlist=None) -> xr.Dataset:
        """Run neural network predictions for all time steps and nutrients."""
        dsout = dsout.copy(deep=True)

        assert (nutrientlist is None) or (type(nutrientlist) is list), \
               "nutrientlist must be None or a list"
        
        if nutrientlist is None:
            nutrientlist = self.config['nutrients']

        print(nutrientlist)
        
        # print(f"\nFitting neural network for {len(self.config['nutrients'])} nutrient(s)")
        
        for it in range(dsmonthly.time.size):
            print(f"\nProcessing time step {it+1}/{dsmonthly.time.size}")
            
            # Prepare data for this time step
            ds_time = dsmonthly.isel(time=[it])
            ds1, predictor_measurements, output_coordinates = self._prepare_inputs(ds_time, nseg)
            
            # Run predictions for each nutrient
            for nutrient in nutrientlist:
                self._predict_nutrient(nseg,
                    nutrient, ds1, predictor_measurements, 
                    output_coordinates, dsout, it, [self.config['year']]
                )
        
        return dsout
    
    def _prepare_inputs(self, ds_time: xr.Dataset, nseg:int) -> Tuple[xr.Dataset, Dict, Dict]:
        """Prepare inputs for neural network prediction."""

        
        # Rename variables for NN
        ds1 = ds_time.rename({
            f'so_segment_{nseg:03d}': 'salinity',
            f'thetao_segment_{nseg:03d}': 'temperature'
        })
        
        # Expand coordinates to 3D
        nz = ds1[f'nz_segment_{nseg:03d}'].size
        lon = ds1[f'lon_segment_{nseg:03d}'].expand_dims({f'nz_segment_{nseg:03d}': nz})
        lat = ds1[f'lat_segment_{nseg:03d}'].expand_dims({f'nz_segment_{nseg:03d}': nz})
        
        # Calculate depth
        z = ds1[f'dz_thetao_segment_{nseg:03d}'].cumsum(dim=f'nz_segment_{nseg:03d}')
        
        # Add coordinates to dataset
        ds1['depth'] = z.squeeze()
        ds1['longitude'] = lon
        ds1['latitude'] = lat
        
        # Prepare dictionaries for NN
        predictor_measurements = {
            k: list(ds1[k].values.squeeze().ravel())
            for k in ["salinity", "temperature"]
        }
        
        output_coordinates = {
            k: list(ds1[k].values.squeeze().ravel())
            for k in ["longitude", "latitude", "depth"]
        }
        
        return ds1, predictor_measurements, output_coordinates
    
    def _predict_nutrient(self, nseg: int, nutrient: str, ds1: xr.Dataset, 
                         predictor_measurements: Dict, output_coordinates: Dict,
                         dsout: xr.Dataset, time_idx: int, year:int):
        """Predict single nutrient for one time step."""
        print(f"  Predicting {nutrient}")
        
        estimates_nn, uncertainties_nn = nn(
            [nutrient],
            self.config['pyesper_path'],
            output_coordinates,
            predictor_measurements,
            EstDates=year,
            Equations=[8]
        )
        
        # Reshape and store results
        var_name = self.config['nutrient_mapping'][nutrient]
        outdata = np.reshape(estimates_nn[f'{nutrient}8'], ds1.salinity.shape)
        outdata[outdata<0] = 0
        dsout[f'{var_name}_segment_{nseg:03d}'].values[time_idx] = outdata


class NetCDFWriter:
    """Handle writing data to NetCDF files."""
    
    def __init__(self, config: Config):
        self.config = config
    
    def write(self, ds: xr.Dataset, varnames: str, nseg:int,  suffix: Optional[str] = None):
        """Write dataset to NetCDF file with proper formatting."""
        
        # Set fill values
        for var in ds:
            ds[var].encoding['_FillValue'] = 1.0e20
        
        # Construct filename
        filename = f'{varnames}_{nseg:03d}_{suffix}.nc' if suffix else f'{varnames}_{nseg:03d}.nc'
        filepath = osp.join(self.config['output_path'], filename)
        
        # Set coordinate attributes and encoding
        self._set_coordinate_encoding(ds, nseg)
        
        # Handle time encoding if needed
        self._set_time_encoding(ds)
        
        # Save to file
        print(f"\nSaving to: {filepath}")
        ds.to_netcdf(
            filepath,
            format='NETCDF3_64BIT',
            engine='netcdf4',
            unlimited_dims='time'
        )
    
    def _set_coordinate_encoding(self, ds: xr.Dataset, nseg: int):
        """Set encoding and attributes for spatial coordinates."""
        ds[f'lon_segment_{nseg:03d}'].encoding['dtype'] = 'float64'
        ds[f'lat_segment_{nseg:03d}'].encoding['dtype'] = 'float64'
        ds[f'lon_segment_{nseg:03d}'].attrs['units'] = 'degree_east'
        ds[f'lat_segment_{nseg:03d}'].attrs['units'] = 'degree_north'
    
    def _set_time_encoding(self, ds: xr.Dataset):
        """Set time encoding if needed."""
        if 'calendar' not in ds['time'].attrs and 'modulo' not in ds['time'].attrs:
            ds.time.encoding['calendar'] = 'gregorian'
            ds.time.encoding['dtype'] = 'float64'
            ds.time.encoding['_FillValue'] = 1.0e20



def read_config(config_file):
    with open(config_file, 'r') as stream:
        config = yaml.safe_load(stream)
    return config


def main():
    import argparse
    import yaml

    parser = argparse.ArgumentParser(description='ESPER estimates')
    parser.add_argument('--config', type=str, default='config.yaml')
    parser.add_argument('--year', type=int,
                        help='Single year to process')

    args = parser.parse_args()


    with open(args.config) as f:
        config = yaml.safe_load(f)
    config2 = config['boundary']
    config = config['esper']
    if args.year is not None:
        config['year'] = args.year

    data_loader = OceanDataLoader(config)
    predictor = NeuralNetworkPredictor(config)
    writer = NetCDFWriter(config)


    for segid in config2['segments']:

        ds, dsout_template1, z = data_loader.load_all_data(segment=segid['id'])

        for nutrient in config['nutrients']:

            dsout1 = predictor.predict(ds,dsout_template1,segid['id'], nutrientlist=[nutrient])
            writer.write(dsout1, nutrient, segid['id'], suffix=config.year)


if __name__ == "__main__":
    main()

    # config  = Config.from_yaml('config.yaml', config_key='esper')
    # config2 = read_config('config.yaml')['boundary']

    # # config = Config.from_script(1993,'002','resample',['nitrate'])
    # data_loader = OceanDataLoader(config)
    # predictor = NeuralNetworkPredictor(config)
    # writer = NetCDFWriter(config)


    # for segid in config2['segments']:

    #     ds, dsout_template1, z = data_loader.load_all_data(segment=segid['id'])

    #     for nutrient in config.nutrients:

    #         dsout1 = predictor.predict(ds,dsout_template1,segid['id'], nutrientlist=[nutrient])
    #         writer.write(dsout1, nutrient, segid['id'], suffix=config.year)

    #         # # Create and run processor
    #         # processor = NutrientProcessor(config)
    #         # processor.process()