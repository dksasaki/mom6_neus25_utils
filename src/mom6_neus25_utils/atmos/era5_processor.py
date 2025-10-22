import xarray as xr
import glob
import os.path as osp
import numpy as np
import netCDF4
import pandas as pd
import yaml
import argparse
import warnings
from datetime import datetime as dtt
warnings.filterwarnings('ignore')


"""
# Process all years and pad (full pipeline)
python era5_processor.py --config era5_sfc.yaml

# Process single year only (no padding)
python era5_processor.py --config era5_sfc.yaml --year 2020

# Run only padding step (for SLURM workflows)
python era5_processor.py --config era5_sfc.yaml --pad-only

# Force overwrite existing files
python era5_processor.py --config era5_sfc.yaml --year 2020 --force
"""

class ERA5Processor:
    """Process ERA5 meteorological data files for a single year."""
    
    # Constants
    variables = {
        'ssrd': 'surface_solar_radiation_downwards',
        'sp': 'surface_pressure',
        'd2m': '2m_dewpoint_temperature',
        'u10': '10m_u_component_of_wind',
        'v10': '10m_v_component_of_wind',
        'tp': 'total_precipitation',
        'msl': 'mean_sea_level_pressure',
        'sf': 'snowfall',
        'strd': 'surface_thermal_radiation_downwards',
        't2m': '2m_temperature',
    }
    
    MOLECULAR_WEIGHT_RATIO = 0.622
    SAT_PRESSURE_0C = 6.112e2  # Pa
    MISSING_VALUE_THRESHOLD = -1e9
    TIME_UNITS = 'days since 1990-01-01T00:00:00'
    
    def __init__(self, config_file, year):
        """Initialize processor with configuration and target year."""
        self.config = self._load_config(config_file)
        self.year = year
        self.output_dir = self.config.get('output_dir', '.')
        self.era5_dir = self.config.get('era5_dir', '.')
        self.tmp_dir = self.config.get('tmp_dir', '.')
        
        print(f"Initialized ERA5 processor for year {self.year}")
    

    
    def _load_config(self, config_file):
        """Load YAML configuration file."""
        with open(config_file, 'r') as file:
            return yaml.safe_load(file)

    def process_year(self, force_save=True, variables=None):
        """Process all variables for the specified year."""
        if variables is None:
            variables = list(self.variables.keys())
        
        print(f"Starting processing for year {self.year}")
        print(f"Variables to process: {variables}")
        
        success_count = 0
        
        # Process each variable
        for variable in variables:
            if self._process_variable(variable, force_save):
                success_count += 1
            else:
                print(f"Failed to process variable: {variable}")
        
        print(f"Successfully processed {success_count}/{len(variables)} variables")
        
        # Calculate derived variables if base variables were processed successfully
        derived_success = 0
        
        if self.calculate_specific_humidity(force_save):
            derived_success += 1
        
        if self.calculate_liquid_precipitation(force_save):
            derived_success += 1
        
        print(f"Successfully calculated {derived_success}/2 derived variables")
        print(f"Processing complete for year {self.year}")
        
        return success_count, derived_success

    def _get_variable_files(self, variable):
        """Get input files for a specific variable and year."""
        file_pattern = osp.join(self.era5_dir, f'ERA5_{variable}_%d.grib')
        current_year_file = glob.glob(file_pattern % self.year)
        next_year_file = glob.glob(file_pattern % (self.year + 1))
        return current_year_file, next_year_file
    
    def _file_exists(self, filepath):
        """Check if output file already exists."""
        return bool(glob.glob(filepath))
    
    def _process_variable(self, variable, force_save=True):
        """Process ERA5 data for a single variable and year."""
        current_files, next_files = self._get_variable_files(variable)
        
        if not current_files:
            raise IOError(f"No file found for variable {variable}, year {self.year}")
            return False
        
        output_file = osp.join(
            self.tmp_dir, 
            f'ERA5_{variable}_{self.year}.nc'
        )
        
        if self._file_exists(output_file) and not force_save:
            print(f'Already exists: {output_file}')
            return True
        
        print(f"Processing {variable} for year {self.year}")
        
        # Load primary dataset
        ds = xr.open_dataset(current_files[0], engine='cfgrib')
        
        # # Try to concatenate with next year's data if available
        # if next_files:
        #     with xr.open_dataset(next_files[0], engine='cfgrib') as aux:
        #         aux = aux.sel(time=slice(dtt(self.year+1,1,1),dtt(self.year+1,1,2)))
        #         ds = xr.concat([ds, aux], dim='time')
        #         print(f"  Concatenated with {self.year + 1} data")
        #         ds = ds.drop_duplicates('time')

        # Process dataset based on dimensions
        if 'step' in ds[list(ds.data_vars.keys())].dims:
            processed_ds = self._process_stepped_data(ds)
        else:
            processed_ds = self._process_regular_data(ds)
        
        # Apply final adjustments
        processed_ds = self._apply_coordinate_adjustments(processed_ds)
        
        # Save processed dataset
        self._save_dataset(processed_ds, output_file)
        print(f"  Saved: {output_file}")
        
        # Clean up
        ds.close()
        processed_ds.close()
        return True
    
    def _process_stepped_data(self, ds):
        """Process data with time-step dimensions."""
        stacked = ds.stack(time_step=('time', 'step'))
        combined_time = stacked.time + stacked.step
        datetime_index = pd.to_datetime(combined_time)
        
        stacked = stacked.assign_coords(time_step=datetime_index)
        stacked = stacked.dropna(dim='time_step')
        stacked = stacked.sel(time_step=slice(f'{self.year}-01-01', f'{self.year+1}-01-03'))
        
        # Convert to numeric time
        time_values = pd.to_datetime(stacked.time_step.values)
        numeric_time = np.array([
            netCDF4.date2num(t, self.TIME_UNITS, calendar='gregorian') 
            for t in time_values
        ])
        
        result = stacked[list(ds.data_vars.keys())].transpose('time_step', 'latitude', 'longitude')
        result = result.rename(time_step='time')
        result = result.assign_coords(time=numeric_time)
        return result
    
    def _process_regular_data(self, ds):
        """Process regular time-series data."""
        time_values = pd.to_datetime(ds.time)
        numeric_time = np.array([
            netCDF4.date2num(t, self.TIME_UNITS, calendar='gregorian') 
            for t in time_values
        ]).astype(float)
        
        result = ds[list(ds.data_vars.keys())]
        result = result.assign_coords(time=numeric_time)
        if 'step' in result.coords:
            result = result.drop('step')
        return result
    
    def _apply_coordinate_adjustments(self, ds):
        """Apply coordinate system adjustments."""

        # Flip data arrays to match reversed latitude
        if ds.latitude[0] > ds.latitude[-1]:
            print('reversing latitude')

            for var in ds.data_vars.keys():
                if 'time_step' not in var:
                    ds[var].values = ds[var].values[:, ::-1, :]

            ds = ds.assign_coords(latitude=ds.latitude[::-1])

        ds = ds.drop_duplicates(dim='time')
        
        
        # Set coordinate attributes
        ds.longitude.attrs['axis'] = "X"
        ds.latitude.attrs['axis'] = "Y"
        ds.time.attrs['axis'] = "T"
        ds.time.attrs['units'] = self.TIME_UNITS
        ds.time.attrs['calendar'] = 'gregorian'
        
        return ds
    
    def _save_dataset(self, ds, filepath):
        """Save dataset to NetCDF file."""
        ds.to_netcdf(
            filepath,
            encoding={'time': {'dtype': float}},
            unlimited_dims='time'
        )
    
    def calculate_liquid_precipitation(self, force_save=True):
        """Calculate liquid precipitation (tp - sf) for the year."""
        tmp_dir = self.tmp_dir
        output_file = osp.join(tmp_dir, f'ERA5_lp_{self.year}.nc')
        
        if self._file_exists(output_file) and not force_save:
            print(f'Already exists: {output_file}')
            return True
        
        # Check if required input files exist
        tp_file = osp.join(tmp_dir, f'ERA5_tp_{self.year}.nc')
        sf_file = osp.join(tmp_dir, f'ERA5_sf_{self.year}.nc')
        
        if not (self._file_exists(tp_file) and self._file_exists(sf_file)):
            print(f"Required input files missing for liquid precipitation calculation")
            raise IOError
        
        print(f"Calculating liquid precipitation for year {self.year}")
        
        # Load required datasets
        with xr.open_dataset(tp_file) as tp_data, xr.open_dataset(sf_file) as sf_data:
            # Calculate liquid precipitation
            lp_data = tp_data['tp'] - sf_data['sf']
            # raise IOError
            lp_data = lp_data.where(lp_data >= self.MISSING_VALUE_THRESHOLD, 0)
            
            # Convert to dataset and clean up
            lp_dataset = lp_data.to_dataset(name='lp')
            if 'valid_time' in lp_dataset.coords:
                lp_dataset = lp_dataset.drop_vars('valid_time')
            
            # Save result
            lp_dataset.to_netcdf(output_file, unlimited_dims='time')
  
        return True
    
    def calculate_specific_humidity(self, force_save=True):
        """Calculate specific humidity from pressure and dewpoint temperature for the year."""
        output_file = osp.join(self.tmp_dir, f'ERA5_sphum_{self.year}.nc')
        
        if self._file_exists(output_file) and not force_save:
            print(f'Already exists: {output_file}')
            return True

        # Check if required input files exist
        pressure_file = osp.join(self.tmp_dir, f'ERA5_sp_{self.year}.nc')
        dewpoint_file = osp.join(self.tmp_dir, f'ERA5_d2m_{self.year}.nc')
        
        if not (self._file_exists(pressure_file) and self._file_exists(dewpoint_file)):
            print(f"Required input files missing for specific humidity calculation")
            raise IOError
        
        print(f"Calculating specific humidity for year {self.year}")
        
        # Load required datasets
        with xr.open_dataset(pressure_file) as p_data, xr.open_dataset(dewpoint_file) as t_data:
            pressure = p_data['sp']  # Pa
            dewpoint_temp = t_data['d2m']  # K
            
            # Calculate specific humidity
            sat_mixing_ratio = self._saturation_mixing_ratio(pressure, dewpoint_temp)
            specific_humidity = self._specific_humidity_from_mixing_ratio(sat_mixing_ratio)
            
            # Create dataset
            sphum_dataset = specific_humidity.to_dataset(name='sphum')
            if 'valid_time' in sphum_dataset.coords:
                sphum_dataset = sphum_dataset.drop_vars('valid_time')
            
            # Prepare encoding
            all_vars = list(sphum_dataset.data_vars.keys()) + list(sphum_dataset.coords.keys())
            encodings = {v: {'_FillValue': None} for v in all_vars}
            encodings['time'].update({
                'dtype': 'float64',
                'calendar': 'gregorian',
                'units': 'hours since 1990-01-01 00:00:00'
            })
            
            # Save result
            sphum_dataset.to_netcdf(
                output_file,
                encoding=encodings,
                unlimited_dims=['time']
            )
        
        print(f"  Saved: {output_file}")
        return True
    
    def _mixing_ratio(self, partial_press, total_press):
        """Calculate mixing ratio from partial and total pressure."""
        return (self.MOLECULAR_WEIGHT_RATIO * partial_press / 
                (total_press - partial_press))
    
    def _specific_humidity_from_mixing_ratio(self, mixing_ratio):
        """Convert mixing ratio to specific humidity."""
        return mixing_ratio / (1 + mixing_ratio)
    
    def _saturation_vapor_pressure(self, temperature):
        """Calculate saturation vapor pressure from temperature."""
        return (self.SAT_PRESSURE_0C * 
                np.exp(17.67 * (temperature - 273.15) / (temperature - 29.65)))
    
    def _saturation_mixing_ratio(self, total_press, temperature):
        """Calculate saturation mixing ratio."""
        sat_vapor_press = self._saturation_vapor_pressure(temperature)
        return self._mixing_ratio(sat_vapor_press, total_press)
    
    def pad_year(self, single_year=None, variables=None):
        """Pad all variables for specified year(s) with next year overlap."""

        if variables is None:
            variables = list(self.variables.keys())
        
        if single_year:
            years_to_pad = [single_year]
        else:
            first_year = self.config.get('first_year')
            last_year = self.config.get('last_year')
            years_to_pad = range(first_year, last_year + 1)
        
        for year in years_to_pad:
            print(f"Padding year {year}")
            all_variables = list(self.variables.keys()) + ['sphum', 'lp']

            for var in all_variables:
                success = self._pad_single_variable(var, year)
                if not success:
                    print(f"  Warning: Could not pad {var} for year {year}")

    def _pad_single_variable(self, var: str, year: int):
        """Pad a single variable with overlap from next year."""
        current_file = osp.join(self.tmp_dir, f'ERA5_{var}_{year}.nc')
        next_file = osp.join(self.tmp_dir, f'ERA5_{var}_{year+1}.nc')
        output_file = osp.join(self.output_dir, f'ERA5_{var}_{year}_padded.nc')
        
        # Check if current year file exists
        if not osp.exists(current_file):
            print(f"  Current year file not found: {current_file}")
            return False
        
        print(f"  Padding {var} for year {year}")
        
        # Load current year data
        ds = xr.open_dataset(current_file)
        
        # Try to concatenate with next year if available
        if osp.exists(next_file):
            ds1 = xr.open_dataset(next_file)
            print(f"    Concatenating {ds.time[-1].values} with {ds1.time[0].values}")
            ds = xr.concat([ds, ds1.isel(time=[0])], dim='time')
            ds1.close()
        else:
            print(f"    No next year file found for {var}, saving without padding")
        
        # Save padded file
        ds.to_netcdf(
            output_file,
            format='NETCDF3_64BIT',
            engine='netcdf4',
            unlimited_dims='time'
        )
        
        print(f"    Saved: {output_file}")
        ds.close()
        return True



def main():
    """Main processing function with flexible year handling."""
    parser = argparse.ArgumentParser(description='Process ERA5 surface data')
    parser.add_argument('--config', type=str, default='era5_processor.yaml',
                        help='YAML configuration file path')
    parser.add_argument('--year', type=int,
                        help='Single year to process')
    parser.add_argument('--variables', nargs='*', 
                        help='Specific variables to process (default: all)')
    parser.add_argument('--force', action='store_true',
                        help='Force overwrite existing files')
    parser.add_argument('--pad-only', action='store_true',
                        help='Only run padding step (skip main processing)')
    

    args = parser.parse_args()
    
    # Load config to get year range if --year not specified
    with open(args.config, 'r') as file:
        config = yaml.safe_load(file)
    
    # Determine years to process
    if args.year:
        years_to_process = [args.year]
        print(f"Processing single year: {args.year}")
    else:
        first_year = config.get('first_year', 1993)
        last_year = config.get('last_year', 1994)
        years_to_process = list(range(first_year, last_year+1))
        print(f"Processing years from config: {first_year} to {last_year}")
    
    # Handle pad-only mode
    if args.pad_only:
        print("\nRunning padding step only...")
        # Create a processor with the full year range from config
        processor = ERA5Processor(args.config, config.get('first_year', 1993))
        processor.pad_year()
        return
    

    total_success = 0
    total_derived_success = 0
    failed_years = []
    
    # Process each year
    for year in years_to_process:
        print(f"\n{'='*50}")
        print(f"PROCESSING YEAR {year}")
        print(f"{'='*50}")
        
        processor = ERA5Processor(args.config, year)
        var_success, derived_success = processor.process_year(
            force_save=args.force, 
            variables=args.variables
        )
        
        total_success += var_success
        total_derived_success += derived_success
        
        if var_success == 0:
            failed_years.append(year)
            print(f"WARNING: No variables processed for year {year}")
        else:
            print(f"SUCCESS: Year {year} completed")


    # Pad all files with next year overlap (only if not processing single year)
    if not args.year:  # Only pad when processing full range, not single years
        print(f"\n{'='*50}")
        print("PADDING FILES WITH NEXT YEAR OVERLAP")
        print(f"{'='*50}")
        
        if years_to_process:
            padding_processor = ERA5Processor(args.config, years_to_process[0])
            padding_processor.pad_year()

    
    # Final summary
    print(f"\n{'='*50}")
    print(f"FINAL SUMMARY")
    print(f"{'='*50}")
    print(f"Years processed: {len(years_to_process)}")
    print(f"Total variables processed: {total_success}")
    print(f"Total derived variables calculated: {total_derived_success}")
    


if __name__ == '__main__':
    main()