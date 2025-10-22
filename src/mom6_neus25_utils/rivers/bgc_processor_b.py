import pandas as pd
import xarray as xr
import numpy as np
from scipy.spatial.distance import cdist
import inspect
from scipy.interpolate import griddata
import cftime
import yaml
import argparse


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

class GlobalNews2DataManager:
    def __init__(self, filename):
        self.filename = filename
        # Raw data as xarray datasets
        self.basin_ds = None
        self.hydrology_ds = None
        self.loading_ds = None

        # bio-availability of phosphorus and the fractionation of dissolved organic
        self.frac_PP=0.3
        self.load_dataset()
        
    def load_dataset(self):
        # Load Excel sheets 
        self._basin_df = pd.read_excel(self.filename, sheet_name=1)
        self._hydrology_df = pd.read_excel(self.filename, sheet_name=2)  
        self._loading_df = pd.read_excel(self.filename, sheet_name=3)
        
        # Convert to xarray and do unit conversions
        self.basin_ds = self._basin_df.to_xarray()
        self._hydrology_ds = self._hydrology_df.to_xarray()
        self._loading_ds = self._loading_df.to_xarray()

        self._get_loading_data()
        self._get_hydrology_data()

        
    def _get_hydrology_data(self):
        # Conversion factor: Mg/yr → g/s
        km3_to_m3 = 1e9
        seconds_per_year = 86400 * 365
        vol_to_time_factor = km3_to_m3 / seconds_per_year

        varb = ['Qact', 'Qnat']
        self.hydrology_ds = self._hydrology_ds[varb].copy()

        for v in varb:
            self.hydrology_ds[v+'_all'] = self.hydrology_ds[v] * vol_to_time_factor


        return self.hydrology_ds
            
    def _get_loading_data(self):
        """Convert nutrient loads from Mg/yr to mol/s using molecular weights"""

        # Conversion factor: Mg/yr → g/s
        Mg_to_g = 1e6
        seconds_per_year = 86400 * 365
        mass_to_time_factor = Mg_to_g / seconds_per_year


        # Molecular weights (g/mol)
        molecular_weights = {
            'Ld_DIN': 14,    'Ld_DON': 14,    'Ld_PN': 14,     # Nitrogen
            'Ld_DIP': 31,    'Ld_DOP': 31,    'Ld_PP': 31,     # Phosphorus
            'Ld_DSi': 28.1                                     # Silicon
        }
        
        molecular_weights['Ld_PP'] *= self.frac_PP
        
        # Apply unit conversion: Mg/yr → mol/s
        self.loading_ds = self._loading_ds[list(molecular_weights.keys())].copy()
        
        for v, mol_weight in molecular_weights.items():
            vaux = v.split('_')[-1]
            self.loading_ds[vaux+'_load'] = (
                self.loading_ds[v] * mass_to_time_factor / mol_weight
            )

    def get_basin_data(self):
        return self.basin_ds
        
    def get_hydrology_data(self):
        return self.hydrology_ds
        
    def get_loading_data(self):
        return self.loading_ds


class ModelGridDataManager:
    def __init__(self, grid_file):
        self.grid_file = grid_file
        self._load_data()
        
    def _load_data(self):
        """Load runoff data and convert from kg m-2 s-1 to m3 s-1"""
        # Water density conversion factor
        water_density = 1e3  # kg/m3
        
        with xr.open_dataset(self.grid_file) as ds:
            self.grid_ds = ds.copy(deep=True)
            
            # Convert runoff: kg/m2/s * m2 / (kg/m3) = m3/s
            self.grid_ds['Q_monthly'] = (
                self.grid_ds['runoff'] * self.grid_ds['area'] / water_density
            )
            
            # Calculate annual mean discharge
            self.grid_ds['Q_ann'] = self.grid_ds['Q_monthly'].mean(dim='month')

            
    def get_grid_data(self):
        return self.grid_ds


class RiverFilter:
    def __init__(self, grid_manager, news_manager, q_min=100):
        self.news_manager = news_manager
        self.grid_manager = grid_manager
        self.q_min = q_min
        
    def filter_rivers(self):
        grid_bounds = self._get_grid_bounds()  #(lon0,lon1,lat0,lat1)
        self._river_ds = self._apply_regional_filter(grid_bounds)
        self._calculate_diag_concentrations()

        
    def _apply_regional_filter(self, grid_bounds):
        # Get grid bounds from ModelGridDataManager

        aux1 = self.news_manager.basin_ds
        aux2 = self.news_manager.hydrology_ds
        aux3 = self.news_manager.loading_ds

        aux = xr.merge([aux1,aux2, aux3])

        aux = aux.where(
            (aux['mouth_lon'] >= grid_bounds[0]) & \
            (aux['mouth_lon'] <= grid_bounds[1]) & \
            (aux['mouth_lat'] >= grid_bounds[2]) & \
            (aux['mouth_lat'] <= grid_bounds[3]) & \
            np.isfinite(aux['Qact_all']) & \
            (aux['Qact_all']>self.q_min),
            drop=True
            )
        isort = np.argsort(aux['Qact'])        
        aux = aux.isel(index=isort.values)
        return aux

    def _get_grid_bounds(self):
        # Extract lat/lon bounds from grid_manager
        grid_ds = self.grid_manager.grid_ds
        # Use grid to filter rivers outside domain
        lon_mod_min = grid_ds['lon'].min()
        lon_mod_max = grid_ds['lon'].max()
        lat_mod_min = grid_ds['lat'].min()
        lat_mod_max = grid_ds['lat'].max()
        return (lon_mod_min, lon_mod_max, lat_mod_min, lat_mod_max)

    def _calculate_diag_concentrations(self):

        ivar = [i for i in self._river_ds.keys() if 'load' in i]

        # Loads are in moles N sec-1, Q in m3 s-1; conc in moles N m-3
        aux = self._river_ds[ivar]/self._river_ds['Qact_all']

        vdict = {}
        for v in aux.data_vars.keys():
            vdict[v] = v.replace('load', 'conc')

        aux = aux.rename(**vdict)

        aux['N_load'] = self._river_ds['DIN_load']*0
        aux['P_load'] = self._river_ds['DIP_load']*0

        aux['N_load'].values = self._river_ds['DIN_load'].values + \
                        self._river_ds['DON_load'].values + \
                        self._river_ds['PN_load'].values

        aux['P_load'].values = self._river_ds['DIP_load'].values + \
                        self._river_ds['DOP_load'].values + \
                        self._river_ds['PP_load'].values 

        self.river_ds = xr.merge([self._river_ds, aux])
            
    # Add getter method:
    def get_filtered_rivers(self):
        return self.river_ds



class RiverMapper:
    def __init__(self, grid_manager, river_filter, min_dist=2, inspect_map='n'):

        self.grid_manager = grid_manager
        self.river_filter = river_filter
        self.min_dist = min_dist
        self.inspect_map = inspect_map
        self.concentration_vectors = {}
        
        self.inner_varbs = ['DIN_conc', 'DON_conc', 'PN_conc', \
                            'DIP_conc','DOP_conc', 'PP_conc', \
                            'DSi_conc']

        # Fractionation factors from original script
        self.frac_ldon = 0.3
        self.frac_sldon = 0.35
        self.frac_srdon = 0.35
        self.frac_ldop = 0.3
        self.frac_sldop = 0.35
        self.frac_srdop = 0.35
        self.const_fed = 70.0e-6  # Iron concentration

        grid_data = self.map_rivers_to_grid()
        self.ds_final = self._create_final_variables(grid_data)

    def map_rivers_to_grid(self):
        """Main mapping algorithm"""
        self._initialize_concentration_vectors()
        glofas = self._model_ds_vec
        news2 = self.river_filter.river_ds

        
        # Main loop through each river
        num_rivers = len(news2.index)
        for k in range(num_rivers):
            print(f"Processing river {k + 1}")
            # glofas is being modified in the following function
            self._map_single_river(k, glofas, news2)

        # map data in glofas (self._model_ds_vec) onto self.model_ds_grid
        dsgrid =  self._map_vectors_to_grid()
        return dsgrid

        
    def _initialize_concentration_vectors(self):
        varbs = self.inner_varbs
                
        # Get runoff points where Q > 0
        aux = self.grid_manager.grid_ds.stack(points=['x','y'])
        aux2 =self.grid_manager.grid_ds

        #  -- vector dataset -- #
        # Initialize concentration vectors for each nutrient
        aux1= aux.where(aux['Q_ann']>0, drop=True)
        for i in varbs:
            aux1[i] = aux1['Q_ann'].copy(deep=True)*0
        self._model_ds_vec = aux1

        # -- grid dataset -- #
        # Initialize concentration vectors for each nutrient
        for i in varbs:
            aux2[i] = aux2['Q_ann'].copy(deep=True)*0
        self.model_ds_grid = aux2
        
        
    def _map_single_river(self, river_idx, glofas, news2):
        instance_name = self._get_instance_name()

        assert glofas is self._model_ds_vec, \
        f'glofas must be an alias of {instance_name}._model_ds_vec'

        """Handle mapping for one river - core of lines ~220-280"""
        i = river_idx
        dist = cdist([[news2['mouth_lon'][i], news2['mouth_lat'][i]]],
                    np.column_stack([glofas['lon'], glofas['lat']]))
        dist_sort_ind = np.argsort(dist.ravel())
        dist_sort = dist.ravel()[dist_sort_ind]

        
        if dist_sort[0] < self.min_dist:
            Qact = news2['Qact'][i].values
            Qmod = glofas['Q_ann'].values
            Q_sum1 = 0
            Q_sum2 = 0
            n = 0
            
            while Q_sum2 < Qact:
                Q_sum1 = Q_sum2
                n += 1
                Q_sum2 = Q_sum1 + Qmod[dist_sort_ind[n]]
            
            if abs(Q_sum1 - Qact) < abs(Q_sum2 - Qact):
                nrp = n - 1
                print(Q_sum1, Qact)  # a quick check for comparable flow
            else:
                nrp = n
                print(Q_sum2, Qact)  # a quick check for comparable flow
            
            # Update nutrient concentrations
            ind = dist_sort_ind[:nrp]

            varbs = ['DIN_conc', 'DON_conc', 'PN_conc',
                    'DIP_conc', 'DOP_conc', 'PP_conc',
                    'DSi_conc',]

            for v in varbs:
                glofas[v].values[ind] = news2[v].values[i]


    def _get_instance_name(self):
        """Get the variable name that references this instance"""
        frame = inspect.currentframe().f_back
        instance_name = None
        for name, obj in frame.f_locals.items():
            if obj is self:
                instance_name = name
                break
        return instance_name or 'unknown_instance'

    def _plot_river_mapping(self, river_idx, river_data, selected_indices):
        """Plotting functionality - lines ~270-290"""
        if self.inspect_map == 'y':
            # Create the inspection plot
            # Show river location, selected grid points, etc.
            pass    


    def _map_vectors_to_grid(self):
        """Map concentration vectors back to 2D grid"""
        # Find where Q_ann > 0 in the original 2D grid

        grid_2d = self.grid_manager.grid_ds

        aa = np.where(grid_2d['Q_ann'] > 0)
        
        # Map each concentration vector back to 2D
        for var in self.inner_varbs:
            # Get the vector values
            conc_vec = self._model_ds_vec[var].values

            conc_vec1 = self.fill_zeros(self._model_ds_vec['lon'],
                                        self._model_ds_vec['lat'],
                                        conc_vec)
            
            # Assign to 2D grid at runoff locations
            self.model_ds_grid[var].values[aa] = conc_vec1
        return self.model_ds_grid



    # Function to perform nearest neighbor interpolation
    def fill_zeros(self, lon, lat, conc_vec):
        # Find indices where concentration values are zero
        zero_indices = np.where(conc_vec == 0)[0]
        non_zero_indices = np.where(conc_vec > 0)[0]
        
        # Perform nearest neighbor interpolation
        filled_values = griddata((lon[non_zero_indices], lat[non_zero_indices]), 
                                conc_vec[non_zero_indices], 
                                (lon[zero_indices], lat[zero_indices]), 
                                method='nearest')
        
        # Update concentration vector with filled values
        conc_vec[zero_indices] = filled_values
        
        return conc_vec

    def _create_final_variables(self, grid_data):
        """Apply fractionation factors and create xarray Dataset - lines ~580-650"""
        
        # Extract base concentrations from grid_data
        din_conc = grid_data['DIN_conc'].values
        don_conc = grid_data['DON_conc'].values
        dip_conc = grid_data['DIP_conc'].values
        dop_conc = grid_data['DOP_conc'].values
        pn_conc = grid_data['PN_conc'].values
        pp_conc = grid_data['PP_conc'].values
        si_conc = grid_data['DSi_conc'].values
        
        # Apply fractionation to create final variables
        NO3_CONC = din_conc.copy()
        LDON_CONC = self.frac_ldon * don_conc.copy()
        SLDON_CONC = self.frac_sldon * don_conc.copy()
        SRDON_CONC = self.frac_srdon * don_conc.copy()
        PO4_CONC = dip_conc.copy()
        LDOP_CONC = self.frac_ldop * dop_conc.copy()
        SLDOP_CONC = self.frac_sldop * dop_conc.copy()
        SRDOP_CONC = self.frac_srdop * dop_conc.copy()
        NDET_CONC = pn_conc.copy()
        PDET_CONC = pp_conc.copy()
        SI_CONC = si_conc.copy()

        # Add iron concentrations
        FED_CONC = din_conc.copy()
        FEDET_CONC = din_conc.copy()
        FED_CONC[FED_CONC > 0] = self.const_fed
        FEDET_CONC[FEDET_CONC > 0] = 0.0
        
        # Get grid coordinates
        lat_mod = grid_data['lat'].values
        lon_mod = grid_data['lon'].values
        
        # Define dimensions
        nlat, nlon = lat_mod.shape
        
        # Create xarray Dataset
        ds = xr.Dataset()
        
        # Create dimensions
        ds['time'] = xr.DataArray([0], dims='time')
        ds['y'] = xr.DataArray(np.arange(nlat), dims='y')
        ds['x'] = xr.DataArray(np.arange(nlon), dims='x')
        
        # Create coordinate variables
        ds['lat'] = xr.DataArray(lat_mod, dims=('y', 'x'), attrs={'units': 'degrees north'})
        ds['lon'] = xr.DataArray(lon_mod, dims=('y', 'x'), attrs={'units': 'degrees east'})
        
        # Create concentration variables
        ds['NO3_CONC']   = xr.DataArray([NO3_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': 'DIN_CONC'})
        ds['LDON_CONC']  = xr.DataArray([LDON_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.3*DON_CONC'})
        ds['SLDON_CONC'] = xr.DataArray([SLDON_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.35*DON_CONC'})
        ds['SRDON_CONC'] = xr.DataArray([SRDON_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.35*DON_CONC'})
        ds['NDET_CONC']  = xr.DataArray([NDET_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '1.0*PN_CONC'})
        ds['PO4_CONC']   = xr.DataArray([PO4_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': 'PO4_CONC'})
        ds['LDOP_CONC']  = xr.DataArray([LDOP_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.3*DOP_CONC'})
        ds['SLDOP_CONC'] = xr.DataArray([SLDOP_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.35*DOP_CONC'})
        ds['SRDOP_CONC'] = xr.DataArray([SRDOP_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.35*DOP_CONC'})
        ds['PDET_CONC']  = xr.DataArray([PDET_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': '0.3*PP_CONC'})
        ds['FED_CONC']   = xr.DataArray([FED_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': 'FED_CONC'})
        ds['FEDET_CONC'] = xr.DataArray([FEDET_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': 'FEDET_CONC'})
        ds['SI_CONC']    = xr.DataArray([SI_CONC], dims=('time', 'y', 'x'), attrs={'units': 'mol m-3', 'long_name': 'SI_CONC'})


        return ds


    def get_final_ds(self):
        return self.ds_final

def ds_to_necdf(ds):
        t = cftime.datetime(1993, 1, 1, calendar='365_day')
        ds = ds.assign_coords(time=[t])
        ds.time.attrs['modulo'] = 'T'
        # Save to netCDF file

        for d in ds:
            print(d)
            ds[d].values#*1e-6

        ds.to_netcdf(output_file,
                format='NETCDF4',
                engine='netcdf4',
                unlimited_dims='time'
        )


class NutrientPlot:
    def __init__(self, river_mapper):
        self.river_mapper = river_mapper
        
    def plot_summary_concentrations(self):
        """Create 2x3 summary plots - lines ~420-470"""
        # Get the vector data after gap filling
        grid_data = self.river_mapper.model_ds_grid
        
        # Extract concentration vectors (where Q_ann > 0)
        aa = np.where(grid_data['Q_ann'] > 0)
        
        # Get coordinates and discharge for runoff points
        lon_mod_runoff_vec = grid_data['lon'].values[aa]
        lat_mod_runoff_vec = grid_data['lat'].values[aa]
        Q_mod_vec = grid_data['Q_ann'].values[aa]
        
        # Get concentration vectors
        din_conc_vec = grid_data['DIN_conc'].values[aa]
        don_conc_vec = grid_data['DON_conc'].values[aa] 
        pn_conc_vec = grid_data['PN_conc'].values[aa]
        dip_conc_vec = grid_data['DIP_conc'].values[aa]
        dop_conc_vec = grid_data['DOP_conc'].values[aa]
        pp_conc_vec = grid_data['PP_conc'].values[aa]
        
        # Calculate total nutrient concentrations
        totn_conc_vec = din_conc_vec + don_conc_vec + pn_conc_vec
        totp_conc_vec = dip_conc_vec + dop_conc_vec + pp_conc_vec
        
        # Scale marker size with discharge (exactly as original)
        ms_vec = np.zeros_like(Q_mod_vec)
        ms_vec[np.log10(Q_mod_vec) < 0] = 1
        ms_vec[(np.log10(Q_mod_vec) > 0) & (np.log10(Q_mod_vec) < 1)] = 2.5
        ms_vec[(np.log10(Q_mod_vec) > 1) & (np.log10(Q_mod_vec) < 2)] = 10
        ms_vec[(np.log10(Q_mod_vec) > 2) & (np.log10(Q_mod_vec) < 3)] = 25
        ms_vec[np.log10(Q_mod_vec) > 3] = 100

        # Plotting - exactly as original
        fig, axs = plt.subplots(2, 3, figsize=(18, 10))

        # Plot total nitrogen concentration
        axs[0, 0].scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=ms_vec, c=totn_conc_vec * 1e3, cmap='viridis', alpha=0.6)
        axs[0, 0].set_title('Total Nitrogen Concentration (mmoles m$^{-3}$)')
        axs[0, 0].set_xlabel('Longitude')
        axs[0, 0].set_ylabel('Latitude')
        axs[0, 0].grid(True)

        # Plot total phosphorus concentration
        axs[0, 1].scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=ms_vec, c=totp_conc_vec * 1e3, cmap='viridis', alpha=0.6)
        axs[0, 1].set_title('Total Phosphorus Concentration (mmoles m$^{-3}$)')
        axs[0, 1].set_xlabel('Longitude')
        axs[0, 1].set_ylabel('Latitude')
        axs[0, 1].grid(True)

        # Plot N:P ratio
        axs[0, 2].scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=ms_vec, c=totn_conc_vec / totp_conc_vec, cmap='viridis', alpha=0.6)
        axs[0, 2].set_title('N:P Ratio')
        axs[0, 2].set_xlabel('Longitude')
        axs[0, 2].set_ylabel('Latitude')
        axs[0, 2].grid(True)

        # Plot DIN concentration
        axs[1, 0].scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=ms_vec, c=din_conc_vec * 1e3, cmap='viridis', alpha=0.6)
        axs[1, 0].set_title('DIN Concentration (mmoles m$^{-3}$)')
        axs[1, 0].set_xlabel('Longitude')
        axs[1, 0].set_ylabel('Latitude')
        axs[1, 0].grid(True)

        # Plot DON concentration
        axs[1, 1].scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=ms_vec, c=don_conc_vec * 1e3, cmap='viridis', alpha=0.6)
        axs[1, 1].set_title('DON Concentration (mmoles m$^{-3}$)')
        axs[1, 1].set_xlabel('Longitude')
        axs[1, 1].set_ylabel('Latitude')
        axs[1, 1].grid(True)

        # Plot PN concentration
        axs[1, 2].scatter(lon_mod_runoff_vec, lat_mod_runoff_vec, s=ms_vec, c=pn_conc_vec * 1e3, cmap='viridis', alpha=0.6)
        axs[1, 2].set_title('PN Concentration (mmoles m$^{-3}$)')
        axs[1, 2].set_xlabel('Longitude')
        axs[1, 2].set_ylabel('Latitude')
        axs[1, 2].grid(True)

        plt.tight_layout()
        plt.show()

    def plot_3d_concentrations(self):
        """Create 3D scatter plots - lines ~520-620"""
        # Get final dataset with fractionated concentrations
        ds_final = self.river_mapper.get_final_ds()
        
        # Extract coordinates and concentrations
        lon_mod = ds_final['lon'].values
        lat_mod = ds_final['lat'].values
        
        # Extract all concentration arrays (remove time dimension)
        NO3_CONC = ds_final['NO3_CONC'][0].values
        LDON_CONC = ds_final['LDON_CONC'][0].values
        SLDON_CONC = ds_final['SLDON_CONC'][0].values
        SRDON_CONC = ds_final['SRDON_CONC'][0].values
        NDET_CONC = ds_final['NDET_CONC'][0].values
        PO4_CONC = ds_final['PO4_CONC'][0].values
        LDOP_CONC = ds_final['LDOP_CONC'][0].values
        SLDOP_CONC = ds_final['SLDOP_CONC'][0].values
        SRDOP_CONC = ds_final['SRDOP_CONC'][0].values
        PDET_CONC = ds_final['PDET_CONC'][0].values
        FED_CONC = ds_final['FED_CONC'][0].values
        FEDET_CONC = ds_final['FEDET_CONC'][0].values
        SI_CONC = ds_final['SI_CONC'][0].values
        
        # Define marker size
        ms = 8
        plt.close('all')
        
        # Create figure - Nitrogen compounds
        fig = plt.figure(figsize=(15, 12), constrained_layout=True)

        # Subplot 1: log10(NO3 CONC)
        ax1 = fig.add_subplot(2, 3, 1, projection='3d')
        ax1.set_title('log10(NO3 CONC)')
        sc = ax1.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(NO3_CONC.flatten()), s=ms, c=np.log10(NO3_CONC.flatten()), cmap='viridis')
        ax1.set_zlabel('log10(NO3 CONC)')
        ax1.set_xlabel('Longitude')
        ax1.set_ylabel('Latitude')
        ax1.set_zlim([-4, -1])
        plt.colorbar(sc, ax=ax1, shrink=0.8)

        # Subplot 2: log10(LDON CONC)
        ax2 = fig.add_subplot(2, 3, 2, projection='3d')
        ax2.set_title('log10(LDON CONC)')
        sc = ax2.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(LDON_CONC.flatten()), s=ms, c=np.log10(LDON_CONC.flatten()), cmap='viridis')
        ax2.set_zlabel('log10(LDON CONC)')
        ax2.set_xlabel('Longitude')
        ax2.set_ylabel('Latitude')
        ax2.set_zlim([-4, -1])
        plt.colorbar(sc, ax=ax2, shrink=0.8)

        # Subplot 3: log10(SLDON CONC)
        ax3 = fig.add_subplot(2, 3, 3, projection='3d')
        ax3.set_title('log10(SLDON CONC)')
        sc = ax3.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(SLDON_CONC.flatten()), s=ms, c=np.log10(SLDON_CONC.flatten()), cmap='viridis')
        ax3.set_zlabel('log10(SLDON CONC)')
        ax3.set_xlabel('Longitude')
        ax3.set_ylabel('Latitude')
        ax3.set_zlim([-4, -1])
        plt.colorbar(sc, ax=ax3, shrink=0.8)

        # Subplot 4: log10(SRDON CONC)
        ax4 = fig.add_subplot(2, 3, 4, projection='3d')
        ax4.set_title('log10(SRDON CONC)')
        sc = ax4.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(SRDON_CONC.flatten()), s=ms, c=np.log10(SRDON_CONC.flatten()), cmap='viridis')
        ax4.set_zlabel('log10(SRDON CONC)')
        ax4.set_xlabel('Longitude')
        ax4.set_ylabel('Latitude')
        ax4.set_zlim([-4, -1])
        plt.colorbar(sc, ax=ax4, shrink=0.8)

        # Subplot 5: log10(NDET CONC)
        ax5 = fig.add_subplot(2, 3, 5, projection='3d')
        ax5.set_title('log10(NDET CONC)')
        sc = ax5.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(NDET_CONC.flatten()), s=ms, c=np.log10(NDET_CONC.flatten()), cmap='viridis')
        ax5.set_zlabel('log10(NDET CONC)')
        ax5.set_xlabel('Longitude')
        ax5.set_ylabel('Latitude')
        ax5.set_zlim([-4, -1])
        plt.colorbar(sc, ax=ax5, shrink=0.8)

        plt.show()

    def plot_3d_phosphorus_concentrations(self):
        """Create 3D phosphorus plots - Figure 8 from original"""
        ds_final = self.river_mapper.get_final_ds()
        
        lon_mod = ds_final['lon'].values
        lat_mod = ds_final['lat'].values
        
        PO4_CONC = ds_final['PO4_CONC'][0].values
        LDOP_CONC = ds_final['LDOP_CONC'][0].values
        SLDOP_CONC = ds_final['SLDOP_CONC'][0].values
        SRDOP_CONC = ds_final['SRDOP_CONC'][0].values
        PDET_CONC = ds_final['PDET_CONC'][0].values
        
        ms = 8
        
        # Create figure 8 - Phosphorus compounds
        fig = plt.figure(figsize=(15, 12), constrained_layout=True)

        # Subplot 1: log10(PO4 CONC)
        ax6 = fig.add_subplot(2, 3, 1, projection='3d')
        sc6 = ax6.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(PO4_CONC.flatten()), s=ms, c=np.log10(PO4_CONC.flatten()), cmap='viridis')
        ax6.set_title('log10(PO4 CONC)')
        ax6.set_xlabel('Longitude')
        ax6.set_ylabel('Latitude')
        ax6.set_zlabel('log10(PO4 CONC)')
        ax6.set_zlim([-4, -2])
        fig.colorbar(sc6, ax=ax6, shrink=0.8)

        # Subplot 2: log10(LDOP CONC)
        ax7 = fig.add_subplot(2, 3, 2, projection='3d')
        sc7 = ax7.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(LDOP_CONC.flatten()), s=ms, c=np.log10(LDOP_CONC.flatten()), cmap='viridis')
        ax7.set_title('log10(LDOP CONC)')
        ax7.set_xlabel('Longitude')
        ax7.set_ylabel('Latitude')
        ax7.set_zlabel('log10(LDOP CONC)')
        ax7.set_zlim([-4, -2])
        fig.colorbar(sc7, ax=ax7, shrink=0.8)

        # Subplot 3: log10(SLDOP CONC)
        ax8 = fig.add_subplot(2, 3, 3, projection='3d')
        sc8 = ax8.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(SLDOP_CONC.flatten()), s=ms, c=np.log10(SLDOP_CONC.flatten()), cmap='viridis')
        ax8.set_title('log10(SLDOP CONC)')
        ax8.set_xlabel('Longitude')
        ax8.set_ylabel('Latitude')
        ax8.set_zlabel('log10(SLDOP CONC)')
        ax8.set_zlim([-4, -2])
        fig.colorbar(sc8, ax=ax8, shrink=0.8)

        # Subplot 4: log10(SRDOP CONC)
        ax9 = fig.add_subplot(2, 3, 4, projection='3d')
        sc9 = ax9.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(SRDOP_CONC.flatten()), s=ms, c=np.log10(SRDOP_CONC.flatten()), cmap='viridis')
        ax9.set_title('log10(SRDOP CONC)')
        ax9.set_xlabel('Longitude')
        ax9.set_ylabel('Latitude')
        ax9.set_zlabel('log10(SRDOP CONC)')
        ax9.set_zlim([-4, -2])
        fig.colorbar(sc9, ax=ax9, shrink=0.8)

        # Subplot 5: log10(PDET CONC)
        ax10 = fig.add_subplot(2, 3, 5, projection='3d')
        sc10 = ax10.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(PDET_CONC.flatten()), s=ms, c=np.log10(PDET_CONC.flatten()), cmap='viridis')
        ax10.set_title('log10(PDET CONC)')
        ax10.set_xlabel('Longitude')
        ax10.set_ylabel('Latitude')
        ax10.set_zlabel('log10(PDET CONC)')
        ax10.set_zlim([-4, -2])
        fig.colorbar(sc10, ax=ax10, shrink=0.8)

        plt.show()

    def plot_3d_other_concentrations(self):
        """Create 3D plots for iron and silicon - Figure 9 from original"""
        ds_final = self.river_mapper.get_final_ds()
        
        lon_mod = ds_final['lon'].values
        lat_mod = ds_final['lat'].values
        
        FED_CONC = ds_final['FED_CONC'][0].values
        FEDET_CONC = ds_final['FEDET_CONC'][0].values
        SI_CONC = ds_final['SI_CONC'][0].values
        
        ms = 8
        
        # Create figure 9 - Iron and Silicon
        fig = plt.figure(figsize=(15, 8), constrained_layout=True)

        # Subplot 1: log10(FED CONC)
        ax11 = fig.add_subplot(1, 3, 1, projection='3d')
        sc11 = ax11.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(FED_CONC.flatten()), s=ms, c=np.log10(FED_CONC.flatten()), cmap='viridis')
        ax11.set_title('log10(FED CONC)')
        ax11.set_xlabel('Longitude')
        ax11.set_ylabel('Latitude')
        ax11.set_zlabel('log10(FED CONC)')
        ax11.set_zlim([-5, -3])
        fig.colorbar(sc11, ax=ax11, shrink=0.8)

        # Subplot 2: log10(FEDET CONC)
        ax12 = fig.add_subplot(1, 3, 2, projection='3d')
        sc12 = ax12.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(FEDET_CONC.flatten()), s=ms, c=np.log10(FEDET_CONC.flatten()), cmap='viridis')
        ax12.set_title('log10(FEDET CONC)')
        ax12.set_xlabel('Longitude')
        ax12.set_ylabel('Latitude')
        ax12.set_zlabel('log10(FEDET CONC)')
        ax12.set_zlim([-5, -3])
        fig.colorbar(sc12, ax=ax12, shrink=0.8)

        # Subplot 3: log10(SI CONC)
        ax13 = fig.add_subplot(1, 3, 3, projection='3d')
        sc13 = ax13.scatter(lon_mod.flatten(), lat_mod.flatten(), np.log10(SI_CONC.flatten()), s=ms, c=np.log10(SI_CONC.flatten()), cmap='viridis')
        ax13.set_title('log10(SI CONC)')
        ax13.set_xlabel('Longitude')
        ax13.set_ylabel('Latitude')
        ax13.set_zlabel('log10(SI CONC)')
        ax13.set_zlim([-3, 0])
        fig.colorbar(sc13, ax=ax13, shrink=0.8)

        plt.show()

    def plot_all_3d(self):
        """Plot all three 3D figure sets"""
        self.plot_3d_concentrations()       # Figure 7 - Nitrogen
        self.plot_3d_phosphorus_concentrations()  # Figure 8 - Phosphorus  
        self.plot_3d_other_concentrations()     # Figure 9 - Iron/Silicon

def load_config(config_file:str):
    """Load configuration from YAML file"""
    with open(config_file, 'r') as stream:
        return yaml.safe_load(stream)
if __name__ =='__main__':


    parser = argparse.ArgumentParser(description='Calculate river parameters')
    parser.add_argument('--config', type=str,
                        help='YAML configuration file path')
    args = parser.parse_args()

    config = load_config(args.config)
    config = config['bgc_processor']

    basin_file          = config['basin_file']
    output_file         = config['output_file']
    news_file           = config['news_file']
    grid_file           = config['grid_file']
    woa_temp_pattern    = config['woa_temp_pattern']
    runoff_file         = config['runoff_file']
    usgs_chem_file      = config['usgs_chem_file']
    usgs_discharge_file = config['usgs_discharge_file']

    # reads and convert GlobalNEWS2 dataset into xarray dataset
    # data in GlobalNEWS2 is saved in different sheets
    g = GlobalNews2DataManager(basin_file)

    # glofas averages (calculated on mom6 grid)
    m = ModelGridDataManager(runoff_file)


    r = RiverFilter(m,g)
    r.filter_rivers()

    rivm = RiverMapper(m, r)
    ds = rivm.get_final_ds()
    ds_to_necdf(ds)



    nutplt = NutrientPlot(rivm)
    # nutplt.plot_all_3d()
