# this code has been adapted from code created by Andrew Ross - https://github.com/andrew-c-ross/nwa-shared

import datetime as dt
import numpy as np
from os import path
import pandas as pd
import xarray as xr
import yaml
import argparse
from boundary import Segment

# xarray gives a lot of unnecessary warnings
import warnings
warnings.filterwarnings('ignore')


class TidalProcessor:
    """Process TPXO tidal data for ocean model boundaries."""
    
    def __init__(self, config_file):
        """Initialize processor with configuration file."""
        self.config = self._load_config(config_file)
        self.segments = None
        
    def _load_config(self, config_file):
        """Load YAML configuration file."""
        with open(config_file, 'r') as file:
            return yaml.safe_load(file)
    
    def setup_segments(self):
        """Create boundary segments from configuration."""
        hgrid = xr.open_dataset(self.config['hgrid_file'])
        
        self.segments = []
        for seg_config in self.config['segments']:
            segment = Segment(
                seg_config['id'],
                seg_config['border'], 
                hgrid,
                output_dir=self.config['output_dir']
            )
            self.segments.append(segment)
        
        print(f"Created {len(self.segments)} boundary segments")
        return self.segments
    
    def write_tpxo(self):
        """Process TPXO tidal data and write to boundary files."""
        constituents = self.config['constituents']

        constituents = [i for i in range(0,10) if i in constituents]

        tpxo_dir = self.config['tpxo_dir']
        horizontal_subset = self.config.get('horizontal_subset', {})
        horizontal_subset = {v:slice(horizontal_subset[v][0], horizontal_subset[v][1])
                            for v in horizontal_subset}
        
        # Load TPXO elevation data
        tpxo_h = (
            xr.open_dataset(path.join(tpxo_dir, 'h_tpxo9.v5a.nc'))
            .rename({'lon_z': 'lon', 'lat_z': 'lat', 'nc': 'constituent'})
            .isel(constituent=constituents, **horizontal_subset)
        )
    


        h = tpxo_h['ha'] * np.exp(-1j * np.radians(tpxo_h['hp']))
        tpxo_h['hRe'] = np.real(h)
        tpxo_h['hIm'] = np.imag(h)
        
        # Load TPXO u-velocity data  
        tpxo_u = (
            xr.open_dataset(path.join(tpxo_dir, 'u_tpxo9.v5a.nc'))
            .rename({'lon_u': 'lon', 'lat_u': 'lat', 'nc': 'constituent'})
            .isel(constituent=constituents, **horizontal_subset)
        )
        
        tpxo_u['ua'] *= 0.01  # convert to m/s
        u = tpxo_u['ua'] * np.exp(-1j * np.radians(tpxo_u['up']))
        tpxo_u['uRe'] = np.real(u)
        tpxo_u['uIm'] = np.imag(u)
        
        # Load TPXO v-velocity data
        tpxo_v = (
            xr.open_dataset(path.join(tpxo_dir, 'u_tpxo9.v5a.nc'))
            .rename({'lon_v': 'lon', 'lat_v': 'lat', 'nc': 'constituent'})
            .isel(constituent=constituents, **horizontal_subset)
        )
        
        tpxo_v['va'] *= 0.01  # convert to m/s
        v = tpxo_v['va'] * np.exp(-1j * np.radians(tpxo_v['vp']))
        tpxo_v['vRe'] = np.real(v)
        tpxo_v['vIm'] = np.imag(v)
        
        # Create time dimension
        # The date should begin 1 month before first day of simulation
        start_date = self.config.get('start_date', '1992-12-01')
        times = xr.DataArray(
            pd.date_range(start_date, periods=1),
            dims=['time']
        )
        
        # Process each boundary segment
        for seg in self.segments:
            print(f"Processing segment {seg.num}: {seg.border}")
            
            # Process tidal elevation
            seg.regrid_tidal_elevation(
                tpxo_h[['lon', 'lat', 'hRe']],
                tpxo_h[['lon', 'lat', 'hIm']],
                times,
                flood=True
            )
            
            # Process tidal velocities
            seg.regrid_tidal_velocity(
                tpxo_u[['lon', 'lat', 'uRe']],
                tpxo_u[['lon', 'lat', 'uIm']],
                tpxo_v[['lon', 'lat', 'vRe']],
                tpxo_v[['lon', 'lat', 'vIm']],
                times,
                flood=True
            )
        
        print("Tidal processing complete!")
    
    def run(self):
        """Run the complete tidal processing pipeline."""
        self.setup_segments()
        self.write_tpxo()


def main():
    """Main function with argument parsing."""
    parser = argparse.ArgumentParser(description='Process TPXO tidal data for ocean model boundaries')
    parser.add_argument('--config', '-c', type=str, default='tides_obc.yaml',
                        help='YAML configuration file path (default: tides_obc.yaml)')
    
    args = parser.parse_args()
    
    # Initialize and run processor
    processor = TidalProcessor(args.config)
    processor.run()


if __name__ == '__main__':
    main()
