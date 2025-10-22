#!/usr/bin/env python3
"""
GloFAS runoff data processor with object-oriented design
"""
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
import os
import yaml
import argparse
import xarray as xr
import xesmf
import glob


class GloFASProcessor:
    def __init__(self, config: str):

        self.config = self._load_config(config)
        self.ocean_mask = None
        self.hgrid = None
        self.mom_coast_mask = None
        self.glofas_coast_mask = None
        self.ldd = None

    def _load_config(self, config_file:str):
        """Load configuration from YAML file"""
        with open(config_file, 'r') as stream:
            return yaml.safe_load(stream)


    def load_grid_data(self):
        """Load ocean mask and horizontal grid data"""
        self.ocean_mask = xr.open_dataset(self.config['grid_mask_file'])
        self.hgrid = xr.open_dataset(self.config['hgrid_file'])
        return self
    
    def generate_masks(self):
        """Generate coastal masks for both MOM and GloFAS grids"""
        self.mom_coast_mask = self._get_coast_mask(self.ocean_mask['mask'])
        
        glofas_subset = dict(
            lat=slice(self.config['latitude_range']['start'], 
                     self.config['latitude_range']['end']),
            lon=slice(self.config['longitude_range']['start'], 
                     self.config['longitude_range']['end'])
        )
        self.ldd = xr.open_dataset(self.config['ldd_file'])['ldd'].sel(**glofas_subset)
        self.glofas_coast_mask = self._generate_glofas_coast_mask()
        return self
    
    def _get_coast_mask(self, mask: xr.DataArray):
        """Find coastal cells using Alistair's method"""
        ocn_mask = mask.values
        cst_mask = 0 * ocn_mask
        is_ocean = ocn_mask > 0
        
        cst_mask[(is_ocean) & (np.roll(ocn_mask, 1, axis=1) == 0)] = 1
        cst_mask[(is_ocean) & (np.roll(ocn_mask, -1, axis=1) == 0)] = 1
        cst_mask[(is_ocean) & (np.roll(ocn_mask, 1, axis=0) == 0)] = 1
        cst_mask[(is_ocean) & (np.roll(ocn_mask, -1, axis=0) == 0)] = 1
        
        cst_mask[0, :] = 0
        cst_mask[:, 0] = 0
        cst_mask[-1, :] = 0
        cst_mask[:, -1] = 0
        
        return cst_mask

    def _generate_glofas_coast_mask(self):
        """Generate GloFAS coastal mask using pour point identification"""
        adjacent = np.logical_and(self.ldd == 5.0, self._expand_mask_true(np.isnan(self.ldd), 3))
        
        imax = 50
        for i in range(imax):
            npoints = int(adjacent.sum())
            adjacent = np.logical_and(self.ldd == 5.0, self._expand_mask_true(adjacent, 3))
            npoints_new = int(adjacent.sum())
            
            if npoints_new == npoints:
                print(f'Converged after {i+1} iterations')
                break
        else:
            raise Exception('Did not converge')
        
        return adjacent.values
    
    def _expand_mask_true(self, mask, window):
        """Expand true values in mask using sliding window"""
        print(mask)

        wind = sliding_window_view(mask, (window, window))
        wind_mask = wind.any(axis=(2, 3))
        final_mask = np.zeros_like(mask)
        i = int((window - 1) / 2)
        final_mask[i:-i, i:-i] = wind_mask
        return final_mask.astype('bool')
    
    def load_data(self, year:int):
        """Load GloFAS data for specified year"""
        glofas_files = self.config['glofas_files_pattern'].format(year=year)
        ds = xr.open_mfdataset(glob.glob(glofas_files), combine='by_coords')
        
        if ds.longitude[0] > 0:
            ds = ds.assign_coords(longitude=ds.longitude - 360)
        
        glofas = (
            ds
            .rename({'latitude': 'lat', 'longitude': 'lon'})
            .sel(time=slice(f'{year-1}-12-31 00:00:00', f'{year+1}-01-02 00:00:00'),
                 lon=slice(-80.1, -54.9))
            .dis24
        )
        
        return glofas
    
    def process_runoff(self, glofas:xr.Dataset, output_file:str):
        """Process runoff data and write to output file"""
        glofas_regridded = self._regrid_glofas_to_mom(glofas)
        filled_runoff = self._redistribute_to_coast(glofas_regridded)
        dataset = self._create_output_dataset(filled_runoff, glofas)
        self._write_output(dataset, output_file)
        return dataset
    
    def _regrid_glofas_to_mom(self, glofas:xr.Dataset):
        """Regrid GloFAS data to MOM grid using conservative interpolation"""
        glofas_latb = np.arange(glofas['lat'][0] + .05, glofas['lat'][-1] - .051, -.1)
        glofas_lonb = np.arange(glofas['lon'][0] - .05, glofas['lon'][-1] + .051, .1)
        
        lon = self.hgrid.x[1::2, 1::2]
        lonb = self.hgrid.x[::2, ::2]
        lat = self.hgrid.y[1::2, 1::2]
        latb = self.hgrid.y[::2, ::2]
        
        glofas_kg = self._convert_to_kg_per_m2_s(glofas)
        
        glofas_to_mom_con = self._reuse_regrid(
            {'lon': glofas.lon, 'lat': glofas.lat, 'lon_b': glofas_lonb, 'lat_b': glofas_latb},
            {'lat': lat, 'lon': lon, 'lat_b': latb, 'lon_b': lonb},
            method='conservative',
            periodic=True,
            reuse_weights=False,
            filename=os.path.join(self.config['TMPDIR'], 'glofas_to_mom.nc')
        )
        
        glofas_regridded = glofas_to_mom_con(
            glofas_kg.where(self.glofas_coast_mask[::2, ::2] > 0).fillna(0.0)
        )
        
        return glofas_regridded.rename({'nyp': 'ny', 'nxp': 'nx'}).values
    
    def _convert_to_kg_per_m2_s(self, glofas:xr.Dataset):
        """Convert m3/s to kg/m2/s"""
        distance_1deg_equator = 111000.0
        dlon = dlat = 0.1  # glofas resolution
        dx = dlon * np.cos(np.deg2rad(glofas.lat)) * distance_1deg_equator
        dy = ((glofas.lon * 0) + 1) * dlat * distance_1deg_equator
        glofas_area = dx * dy
        return glofas * 1000.0 / glofas_area
    
    def _reuse_regrid(self, *args, **kwargs):
        """Create regridder with optional weight file reuse"""
        filename = kwargs.pop('filename', None)
        reuse_weights = kwargs.pop('reuse_weights', False)
        
        if reuse_weights and filename and os.path.isfile(filename):
            return xesmf.Regridder(*args, reuse_weights=True, filename=filename, **kwargs)
        else:
            regrid = xesmf.Regridder(*args, **kwargs)
            if filename:
                regrid.to_netcdf(filename)
            return regrid
    
    def _redistribute_to_coast(self, glofas_regridded:xr.Dataset):
        """Redistribute runoff to nearest coastal cells"""
        lon = self.hgrid.x[1::2, 1::2]
        lat = self.hgrid.y[1::2, 1::2]
        
        flat_mask = self.mom_coast_mask.ravel().astype('bool')
        coast_lon = lon.values.ravel()[flat_mask]
        coast_lat = lat.values.ravel()[flat_mask]
        mom_id = np.arange(np.prod(self.mom_coast_mask.shape))
        
        coast_to_mom = self._reuse_regrid(
            {'lat': coast_lat, 'lon': coast_lon},
            {'lat': lat, 'lon': lon},
            method='nearest_s2d',
            locstream_in=True,
            reuse_weights=False,
            filename=os.path.join(self.config['TMPDIR'], 'coast_to_mom.nc')
        )
        
        coast_id = mom_id[flat_mask]
        nearest_coast = coast_to_mom(coast_id).ravel()
        
        raw = glofas_regridded.reshape([glofas_regridded.shape[0], -1])
        filled = np.zeros_like(raw)
        
        for i in coast_id:
            filled[:, i] = raw[:, nearest_coast == i].sum(axis=1)
        
        return filled.reshape(glofas_regridded.shape)
    
    def _create_output_dataset(self, filled_runoff:np.array, glofas:xr.Dataset):
        """Create xarray dataset for output"""
        lon = self.hgrid.x[1::2, 1::2]
        lat = self.hgrid.y[1::2, 1::2]

        area = ((self.hgrid.area[::2, ::2] + self.hgrid.area[1::2, 1::2]) + 
                (self.hgrid.area[1::2, ::2] + self.hgrid.area[::2, 1::2]))
        
        ds = xr.Dataset({
            'runoff': (['time', 'y', 'x'], filled_runoff),
            'area': (['y', 'x'], area.data),
            'lat': (['y', 'x'], lat.data),
            'lon': (['y', 'x'], lon.data)
        }, coords={
            'time': glofas['time'].data,
            'y': np.arange(filled_runoff.shape[1]),
            'x': np.arange(filled_runoff.shape[2])
        })
        
        self._set_attributes(ds)
        return ds

    def _set_attributes(self, ds:xr.Dataset):
        """Set variable attributes for output dataset"""
        ds['time'].attrs = {'cartesian_axis': 'T'}
        ds['x'].attrs = {'cartesian_axis': 'X'}
        ds['y'].attrs = {'cartesian_axis': 'Y'}
        ds['lat'].attrs = {'units': 'degrees_north'}
        ds['lon'].attrs = {'units': 'degrees_east'}
        ds['runoff'].attrs = {'units': 'kg m-2 s-1'}
    
    def _write_output(self, dataset:xr.Dataset, output_file:str):
        """Write dataset to NetCDF file"""
        all_vars = list(dataset.data_vars.keys()) + list(dataset.coords.keys())
        encodings = {v: {'_FillValue': None} for v in all_vars}
        
        encodings['time'].update({
            'units': 'days since 1950-01-01',
            'dtype': 'float',
            'calendar': 'gregorian'
        })
        
        dataset.to_netcdf(
            output_file,
            unlimited_dims=['time'],
            format='NETCDF3_64BIT',
            encoding=encodings,
            engine='netcdf4'
        )
        print(f'{output_file} saved successfully')

    
    def process_years(self, start_year:int, end_year:int):
        """Process multiple years of data"""
        if not os.path.exists(self.config['output_dir']):
            os.makedirs(self.config['output_dir'])
        
        results = {}
        for year in range(start_year, end_year + 1):
            print(f"Processing year {year}")
            glofas_data = self.load_data(year)
            output_file = os.path.join(self.config['output_dir'], f'glofas_runoff_{year}.nc')
            results[year] = self.process_runoff(glofas_data, output_file)
        
        return results




def main():

    
    processor = GloFASProcessor('glofas_processor.yaml')
    processor.load_grid_data().generate_masks()
    
    results = processor.process_years(
        processor.config['start_year'], 
        processor.config['end_year']
    )
    
    print(f"Processed {len(results)} years successfully")


if __name__ == '__main__':
    main()