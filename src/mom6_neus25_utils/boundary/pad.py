import xarray as xr
import glob
import yaml
import os.path as osp
import numpy as np

class GlorysOBCProcessor:

    def __init__(self, config_path: str, config_keys: list = None):
        expected_config_keys = {
            'first_year': int,
            'last_year': int,
            'output_dir': str,
            'glorys_dir': str,
            'hgrid': str,
            'ncrcat_years': bool,
            'segments': list,
            'variables': list,
        }

        self.config_path = config_path
        self.config_keys = config_keys or list(expected_config_keys.keys())
        self.config_keys_types = expected_config_keys
        self.config = self._load_config(config_path)
        
    def _load_config(self, config_file: str):
        with open(config_file, 'r') as file:
            config = yaml.safe_load(file)   

        for k in config.keys():
            assert k in self.config_keys, f"glorys_yaml argument" + \
                                     f"'{k}' must be in {self.config_keys}"
            assert isinstance(config[k], self.config_keys_types[k])
        return config


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
        faux = self.config['output_dir']
        fpath = osp.join(faux, f"{var}_{seg['id']:03d}_{year}.nc")
        fpath1 = osp.join(faux, f"{var}_{seg['id']:03d}_{year+1}.nc")

        ds = xr.open_dataset(fpath)

        if osp.exists(fpath1):
            ds1 = xr.open_dataset(fpath1)
            print('concatenating ', ds.time[-1].values, ds1.time[0].values)
            ds = xr.concat([ds, ds1.isel(time=[0])], dim='time')
        else:
            print('nothing to concatenate')
        
        return ds


if __name__ == '__main__':
    fpath = 'glorys_obc.yaml'
    glo = GlorysOBCProcessor(fpath)
    
    glo.pad_year()