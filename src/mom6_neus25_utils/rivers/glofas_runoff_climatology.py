from scipy.io import savemat
import xarray as xr
import os.path as osp
import yaml
import argparse

class Glofas2mat_avg:
    def __init__(self, config: str):
        config = self._load_config(config)
        self.config = config

    def _load_config(self, config_file:str):
        """Load configuration from YAML file"""
        with open(config_file, 'r') as stream:
            return yaml.safe_load(stream)
        
    def read_multiple_years(self):
        """Read glofas estimates considering a range of years"""
        ys = self.config['start_year']
        yf = self.config['end_year']

        self.averages = []
        for y in range(ys, yf+1):
            ave, ds = self._read_glofas(y)
            self.averages.append(ave)

        self.ds = ds
        # self._savemat()
        self._savenc()


    def _read_glofas(self, year: int):
        """Read individual glofas files"""
        fpath = osp.join(self.config['output_dir'], f'glofas_runoff_{year}.nc')
        print(f'reading {fpath}')
        with  xr.open_dataset(fpath) as ds:

            ave = (ds['runoff']
                    .sel(time=slice(f'{year}-01-01', f'{year+1}-01-01'))
                    .groupby('time.month')
                    .mean('time')
                    .compute()
                    )


                # Store coordinates from first file
            if not hasattr(self, 'coords'):
                self.coords = {
                    'lat': ds.lat,
                    'lon': ds.lon,
                    'area': ds.area
                }

            for c in self.coords:
                ave[c] = self.coords[c]

            return ave, ds

    def _savemat(self):
        """Save the monthly average of glofas"""


        ds = self.ds
        fpath = self.config['output_dir']
        outpath = osp.join(fpath, 'glofas_runoff_mean.mat')
        all_average = xr.concat(self.averages, dim='year').mean('year')
        savemat(outpath, self.coords)
        print(outpath + ' saved')

    def _savenc(self):
        """Save the monthly average of glofas"""


        ds = self.ds
        fpath = self.config['output_dir']
        outpath = osp.join(fpath, 'glofas_runoff_mean.nc')
        all_average = xr.concat(self.averages, dim='year').mean('year')
        all_average.to_netcdf(outpath)
        print(outpath + ' saved')


#def main():
#    glo = Glofas2mat_avg('runoff_bgc.yaml')
#    glo.read_multiple_years()

def main():
    parser = argparse.ArgumentParser(description='Average GloFAS runoff data')
    parser.add_argument('--config', type=str, required=True,
                        help='YAML configuration file path')
    args = parser.parse_args()

    glo = Glofas2mat_avg(args.config)
    glo.read_multiple_years()

if __name__ == '__main__':
    main()

