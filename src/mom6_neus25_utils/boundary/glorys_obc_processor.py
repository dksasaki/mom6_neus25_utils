import os
import os.path as osp
import yaml
import xarray as xr
import warnings
from glob import glob
from typing import List, Dict, Optional
import argparse

try:
    from .boundary import Segment 
except ImportError:
    from boundary import Segment 


warnings.filterwarnings('ignore')




class GlorysOBCProcessor:

    def __init__(self, config_file):
        # expected_config_keys = {
        #     'first_year': int,
        #     'last_year': int,
        #     'output_dir': str,
        #     'glorys_dir': str,
        #     'tmp_dir': str,
        #     'hgrid': str,
        #     'ncrcat_years': bool,
        #     'segments': list,
        #     'variables': list,
        # }

        self.config = self._load_config(config_file)
    #     self.config_keys = config_keys or list(expected_config_keys.keys())
    #     self.config_keys_types = expected_config_keys
    #     self.config      = self._load_config(config_path)
        
    # def _load_config(self, config_file: str):
    #     with open(config_file, 'r') as file:
    #         config = yaml.safe_load(file)   

    #     for k in config.keys():
    #         assert k in self.config_keys, f"glorys_yaml argument" + \
    #                                  f"'{k}' must be in {self.config_keys}"
    #         assert isinstance(config[k], self.config_keys_types[k])
    #     return config   

    def _load_config(self, config_file):
        """Load YAML configuration file."""
        with open(config_file, 'r') as file:
            return yaml.safe_load(file)


    def setup(self): 

        tmp_dir = self.config['tmp_dir']

        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)

        hgrid_file = self.config['hgrid']
        self.hgrid = xr.open_dataset(hgrid_file)


        self.segments = self._create_segment()
    
    def _create_segment(self):
        segments = []
        for seg_config in self.config['segments']:
            # segment configures the path where data will be saved
            segment = Segment(seg_config['id'],
                              seg_config['border'],
                              self.hgrid,
                              output_dir=self.config['tmp_dir'])
            segments.append(segment)
        print(self.config['tmp_dir'])
        return segments

    def process_years(self):
        yi = self.config.get('first_year')
        yf = self.config.get('last_year')

        for year in range(yi, yf+1):
            print(f"Processing year: {year}")
            glorys_dir = self.config['glorys_dir'] % year          

            is_first_year = (year == yi)
            is_last_year  = (year == yf)      

            self._process_single_year(year,
                                      glorys_dir,
                                      is_first_year,
                                      is_last_year)
    def _process_single_year(self,
                             year:int,
                             glorys_dir:str,
                             is_first_year:bool = False,
                             is_last_year:bool  = False):
        # Load data
        glorys_data = self._load_glorys_data(glorys_dir,
                                                     is_first_year)
        
        # Process segments
        for segment in self.segments:
            self._process_segment(segment, glorys_data, year)

    def _load_glorys_data(self, glorys_dir: str, is_first_year: bool) -> xr.Dataset:
        """Load and preprocess GLORYS data."""
        # pattern = os.path.join(glorys_dir, '*.nc')
        files = sorted(glob(glorys_dir))

        if not files:
            raise FileNotFoundError(f"No files found matching pattern: {glorys_dir}") 

        glorys = xr.open_mfdataset(files, concat_dim='time', combine='nested')
        glorys = glorys.rename({'latitude': 'lat', 'longitude': 'lon', 'depth': 'z'})
        
        if is_first_year:
            glorys = self._adjust_first_year_time(glorys)
        
        glorys.load()
        return glorys

    def _adjust_first_year_time(self, data: xr.Dataset) -> xr.Dataset:
        """Floor first time to midnight for first year."""
        time_values = data['time'].values
        new_time = xr.concat(
            [data['time'][0].dt.floor('1d'), data['time'][1:]], 
            dim='time'
        )
        data['time'] = new_time
        return data


    def _process_segment(self, segment: Segment, data: xr.Dataset, year: int):
        """Process a single segment with the given data."""
        # Process velocity
        #saves data
        segment.regrid_velocity(data['uo'], data['vo'], suffix=year, flood=False)
        
        # Process tracers
        for tracer in ['thetao', 'so', 'zos']:
            # saves data
            segment.regrid_tracer(data[tracer], suffix=year, flood=False)

    
    def run(self):
        """Run the complete processing pipeline."""
        self.setup()
        self.process_years()


    def pad_year(self):
        """Process all variables and segments for a given year."""
        first_year = self.config['first_year']
        last_year = self.config['last_year']
        faux = self.config['output_dir']

        for year in range(first_year,last_year+1):
            for var in self.config['variables']:
                for seg in self.config['segments']:
                    fpattern = f"{var}_{seg['id']:03d}_{year}_padded.nc"
                    fpath = osp.join(faux, fpattern)
                    ds = self._pad_segment_variable(var, seg, year)


                    ds.to_netcdf(
                        fpath,
                        format='NETCDF3_64BIT',
                        engine='netcdf4',
                        # encoding=encoding,
                        unlimited_dims='time'
                    )

    def _pad_segment_variable(self, var: str, seg: dict, year: int):
        """Process a single variable-segment combination for a given year."""
        faux = self.config['tmp_dir']
        fpath = osp.join(faux, f"{var}_{seg['id']:03d}_{year}.nc")
        fpath1 = osp.join(faux, f"{var}_{seg['id']:03d}_{year+1}.nc")

        ds = xr.open_dataset(fpath)

        if osp.exists(fpath1):
            ds1 = xr.open_dataset(fpath1)
            t0 = ds.time[-1].values
            t1 = ds1.time[0].values
            print('concatenating ', var, seg, t0, t1 )
            ds = xr.concat([ds, ds1.isel(time=[0])], dim='time')
        else:
            print('nothing to concatenate')
        
        return ds


def main():

    parser = argparse.ArgumentParser(description='Process Glorys boundary data')
    parser.add_argument('--config', type=str, default='glorys_processor.yaml')
    parser.add_argument('--year', type=int,
                        help='Single year to process')
    parser.add_argument('--pad-only', action='store_true',
                        help='Only run padding step (skip main processing)')

    args = parser.parse_args()

    # with open(args.config, 'r') as file:
    #     config = yaml.safe_load(file)


    glo = GlorysOBCProcessor(args.config) 

    # Step 1: Process raw GLORYS data and create boundary condition files
    # - Loads GLORYS netCDF files from config['glorys_dir'] 
    # - Regrids velocity (uo, vo) and tracers (thetao, so, zos) to boundary segments
    # - Saves intermediate files to config['tmp_dir'] with naming: {var}_{seg_id:03d}_{year}.nc

    # Step 2: Pad boundary condition files with overlap from next year
    # - Takes each {var}_{seg_id:03d}_{year}.nc file from Step 1
    # - Adds the first timestep from {year+1} file to ensure temporal continuity
    # - Saves padded files as {var}_{seg_id:03d}_{year}_padded.nc

    # Step 1
    if not args.pad_only:
        glo.run()
    else:
    # Step 2
    # IMPORTANT: When using SLURM parallel processing, comment out this line!
    # Padding must occur AFTER all yearly files are created, since it depends
    # on having both current year and next year files available.
    # When next year is absent, padding STILL save a file as 
    # {var}_{seg_id:03d}_{year}_padded.nc
        glo.pad_year()

if __name__ == '__main__':
    main()
