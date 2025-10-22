import xarray as xr
import numpy as np
import xesmf as xe
from bgc_processor_b import RiverMapper
from scipy.spatial.distance import cdist
from scipy.interpolate import griddata
import cftime
import yaml
import argparse

class USGSDataManager:
    def __init__(self, chem_file:str, discharge_file:str, nutrient_option=2):
        self.chem_file = chem_file
        self.discharge_file = discharge_file
        self.nutrient_option = nutrient_option
        
        # Stores monthly arrays (following script exactly)
        self.monthly_data = {}
        self.annual_data = {}
        self.station_metadata = None  # Will be xarray Dataset
        """Load both chemistry and discharge datasets"""

        self._load_chemistry_data()
        self._load_discharge_data()        

        
    def _load_chemistry_data(self):
        """Load and process chemistry NetCDF"""

        with xr.open_dataset(self.chem_file) as ds:
            # Basic chemistry (always available)
            self.monthly_data['alk'] = np.transpose(ds['alk'].values, (1, 0))
            self.monthly_data['dic'] = np.transpose(ds['dic'].values, (1, 0))
            self.monthly_data['no3'] = np.transpose(ds['no3'].values, (1, 0))
            self.monthly_data['nh4'] = np.transpose(ds['nh4'].values, (1, 0))
            self.monthly_data['din'] = self.monthly_data['no3'] + self.monthly_data['nh4']
            self.monthly_data['dip'] = np.transpose(ds['po4'].values, (1, 0))
            self.monthly_data['si'] = np.transpose(ds['sio2'].values, (1, 0))
            self.monthly_data['o2'] = np.transpose(ds['do'].values, (1, 0)) / 2.0
            self.monthly_data['don'] = np.transpose(ds['don'].values, (1, 0))
            self.monthly_data['temp'] = np.transpose(ds['temp'].values, (1, 0))
            
            # Handle nutrient options
            if self.nutrient_option == 1:
                # Set to NaN when not available
                self.monthly_data['pn'] = np.full_like(self.monthly_data['no3'], np.nan)
                self.monthly_data['dop'] = np.full_like(self.monthly_data['no3'], np.nan)
                self.monthly_data['pp'] = np.full_like(self.monthly_data['no3'], np.nan)
                
            elif self.nutrient_option == 2:
                # Calculate from total filtered/unfiltered nutrients
                tnf_monthly = np.transpose(ds['tnf'].values, (1, 0))
                tnu_monthly = np.transpose(ds['tnu'].values, (1, 0))
                tpf_monthly = np.transpose(ds['tpf'].values, (1, 0))
                tpu_monthly = np.transpose(ds['tpu'].values, (1, 0))
                
                don2_monthly = tnf_monthly - self.monthly_data['din']
                pn_monthly = tnu_monthly - tnf_monthly
                pn_monthly[pn_monthly < 0] = np.nan
                
                dop_monthly = tpf_monthly - self.monthly_data['dip']
                dop_monthly[dop_monthly < 0] = np.nan
                
                pp_monthly = (tpu_monthly - tpf_monthly) * 0.5  # frac_PP from script
                pp_monthly[pp_monthly < 0] = np.nan
                
                self.monthly_data['pn'] = pn_monthly
                self.monthly_data['dop'] = dop_monthly
                self.monthly_data['pp'] = pp_monthly
            
            # Iron (always NaN for USGS data)
            self.monthly_data['dfe'] = np.full_like(self.monthly_data['no3'], np.nan)
            self.monthly_data['pfe'] = np.full_like(self.monthly_data['no3'], np.nan)

    def _load_discharge_data(self):
        """Load discharge NetCDF"""
        with xr.open_dataset(self.discharge_file) as ds:
            self.monthly_data['Q'] = np.transpose(ds['disc'].values, (1, 0))  # m^3 sec^-1
            
            # Create station metadata as xarray Dataset
            n_stations = len(ds['river_name'].values)
            self.station_metadata = xr.Dataset({
                'river_name': (['station'], ds['river_name'].values),
                'mouth_lon': (['station'], ds['mouth_lon'].values.copy()),
                'mouth_lat': (['station'], ds['mouth_lat'].values.copy())
            }, coords={'station': range(n_stations)})
            
        # Calculate annual means for all variables
        self._calculate_annual_means()

    def _calculate_annual_means(self):
        """Calculate annual means from monthly data"""
        # Convert monthly data to xarray Dataset first
        n_stations = len(self.station_metadata.station)
        n_months = 12
        
        monthly_vars = {}
        annual_vars = {}
        
        for var, monthly_data in self.monthly_data.items():
            # Create monthly xarray DataArray
            monthly_vars[var] = (['station', 'month'], monthly_data)
            
            # Calculate annual mean
            annual_vars[var] = (['station'], np.nanmean(monthly_data, axis=1))
        
        # Store as xarray Datasets
        self.monthly_data = xr.Dataset(
            monthly_vars, 
            coords={
                'station': self.station_metadata.station,
                'month': range(n_months)
            }
        )
        
        self.annual_data = xr.Dataset(
            annual_vars,
            coords={'station': self.station_metadata.station}
        )

    def apply_coordinate_corrections(self, corrections):
        """Apply manual coordinate fixes

        corrections: dict like {'Susquehanna': {'lat': 38.5, 'lon': -77.5}}
        """

        if corrections is None:
            return
            
        station_names = self.station_metadata['river_name'].values
        
        for i, station_name in enumerate(station_names):
            if station_name in corrections:
                if 'lat' in corrections[station_name] and corrections[station_name]['lat'] is not None:
                    self.station_metadata['mouth_lat'].values[i] = corrections[station_name]['lat']
                if 'lon' in corrections[station_name] and corrections[station_name]['lon'] is not None:
                    self.station_metadata['mouth_lon'].values[i] = corrections[station_name]['lon']

    def get_monthly_data(self):
        """Return monthly datasets as xarray Dataset"""
        return self.monthly_data

    def get_annual_data(self):
        """Return annual datasets as xarray Dataset"""
        return self.annual_data
        
    def get_station_metadata(self):
        """Return station information as xarray Dataset"""
        return self.station_metadata




class GlofasRunoffManager:
    def __init__(self, runoff_file: str):
        self.runoff_file = runoff_file
        self.runoff_ds = None
        self.load_data()
        
    def load_data(self):
        """Load runoff NetCDF and convert units"""
        with xr.open_dataset(self.runoff_file) as ds:
            self.runoff_ds = ds.copy(deep=True)
            
        # Convert runoff from kg m-2 s-1 to m3 s-1
        self._convert_runoff_units()
        
    def _convert_runoff_units(self):
        """Convert runoff from kg m-2 sec-1 to m3 sec-1"""
        # Water density conversion factor
        water_density = 1000  # kg/m3
        
        # Convert runoff: kg/m2/s * m2 / (kg/m3) = m3/s
        self.runoff_ds['Q_monthly'] = (
            self.runoff_ds['runoff'] * self.runoff_ds['area'] / water_density
        )
        
        # Calculate annual mean discharge  
        self.runoff_ds['Q_ann'] = self.runoff_ds['Q_monthly'].mean(dim='month')
        
    def get_runoff_data(self):
        """Return the processed runoff dataset"""
        return self.runoff_ds
        
    def get_runoff_vectors(self):
        """Get vectors at runoff points (where Q_ann > 0) - following original script"""
        # Stack 2D data to 1D vectors first
        stacked_ds = self.runoff_ds.stack(points=['y', 'x'])
        
        # Find points where Q_ann > 0
        mask = stacked_ds['Q_ann'] > 0
        
        # Extract vector data at runoff points only
        runoff_vectors = stacked_ds.where(mask, drop=True)
        
        return runoff_vectors



class WOATemperatureManager:
    def __init__(self, woa_temp_pattern: str, grid_file: str):
        self.woa_temp_pattern = woa_temp_pattern
        self.grid_file = grid_file
        self.woa_temp_ds = None
        self.interpolated_temp_ds = None
        self.regridder = None
        
        # O2 saturation constants
        self.o2_constants = {
            'a_0': 2.00907,
            'a_1': 3.22014, 
            'a_2': 4.05010,
            'a_3': 4.94457,
            'a_4': -2.56847e-1,
            'a_5': 3.88767,
            'b_0': -6.24523e-3,
            'b_1': -7.37614e-3,
            'b_2': -1.03410e-2,
            'b_3': -8.17083e-3,
            'c_0': -4.88682e-7,
            'sal': 0  # Salinity
        }
        self.load_data()
        
    def load_data(self):
        """Load WOA multi-file dataset - lines ~48-51"""
        self.woa_temp_ds = xr.open_mfdataset(
            self.woa_temp_pattern,
            decode_times=False,
            concat_dim='time',
            combine='nested'
        )
        
    def interpolate_to_grid(self):
        """Interpolate WOA temperature onto model grid using XESMF - lines ~275-280"""
        # Load target grid
        target_grid = xr.open_dataset(self.grid_file).rename(geolon='lon', geolat='lat')
        
        # Prepare source grid (WOA)
        source_grid = self.woa_temp_ds
        
        # Create regridder
        self.regridder = xe.Regridder(
            source_grid, 
            target_grid,
            method='bilinear',  # Conservative, bilinear, or nearest_s2d
            periodic=False,
            reuse_weights=False
        )
        
        # Apply regridding to temperature data
        self.interpolated_temp_ds = self.regridder(self.woa_temp_ds[['t_an']].isel(depth=0))
        
    def calculate_o2_saturation_vectors(self, runoff_vectors):
        """Calculate O2 saturation at runoff points - lines ~340-370"""
        n_months = 12
        n_points = len(runoff_vectors.points)
        
        temp_monthly_vecs = np.zeros((n_months, n_points))
        o2sat_monthly_vecs = np.zeros((n_months, n_points))
        
        for m in range(n_months):
            # Get temperature field for this month
            temp = self.interpolated_temp_ds['t_an'][m, :, :].squeeze()
            
            # Fill missing values and apply bounds (lines ~355-358)
            temp = temp.bfill(dim='xh').ffill(dim='yh')
            temp = xr.where(temp > 40, 40, temp)
            temp = xr.where(temp < 0, 0, temp)
            
            # Extract temperature at runoff point locations
            temp_at_runoff = self._extract_temp_at_points(temp, runoff_vectors)
            temp_monthly_vecs[m, :] = temp_at_runoff
            
            # Calculate O2 saturation
            o2sat_monthly_vecs[m, :] = self._calculate_o2_saturation_formula(temp_at_runoff)
        self.o2sat_woa_monthly_vecs = o2sat_monthly_vecs

    def get_o2_saturation_vectors(self):
        return temp_monthly_vecs, o2sat_monthly_vecs
        
    def _extract_temp_at_points(self, temp_field, runoff_vectors):
        """Extract temperature values at runoff point coordinates"""
        
        # runoff_vectors is already stacked with points dimension
        # We need to get back to the 2D indices to extract from temp_field
        
        # Get the multi-index information
        points_idx = runoff_vectors.indexes['points']
        indices = np.array(points_idx.tolist())  # Convert MultiIndex to array

        y_indices = [idx[0] for idx in points_idx]  # y coordinates 
        x_indices = [idx[1] for idx in points_idx]  # x coordinates
        
        # Extract temperature at these indices
        temp_at_points = temp_field.values[y_indices,x_indices]
        
        return temp_at_points
        
    def _calculate_o2_saturation_formula(self, temp):
        """Calculate O2 saturation using Garcia & Gordon formula - lines ~360-367"""
        tt = 298.15 - temp
        tkb = 273.15 + temp
        ts = np.log(tt / tkb)
        ts2 = ts * ts
        ts3 = ts2 * ts  
        ts4 = ts3 * ts
        ts5 = ts4 * ts
        
        c = self.o2_constants
        o2sat = (1000.0 / 22391.6) * 1000 * np.exp(
            c['a_0'] + c['a_1'] * ts + c['a_2'] * ts2 + c['a_3'] * ts3 + 
            c['a_4'] * ts4 + c['a_5'] * ts5 +
            (c['b_0'] + c['b_1'] * ts + c['b_2'] * ts2 + c['b_3'] * ts3 + 
             c['c_0'] * c['sal']) * c['sal']
        )
        
        return o2sat
    
    def cleanup(self):
        """Clean up regridder to free memory"""
        if self.regridder is not None:
            self.regridder.clean_weight_file()



class GlobalNEWS2RatioManager:
    def __init__(self, news_file: str):
        self.news_file = news_file
        self.news_ds = None
        self.ratio_ds = None
        
    def load_data(self):
        """Load GlobalNEWS2 NetCDF file"""
        self.news_ds = xr.open_dataset(self.news_file)
        
    def calculate_ratio_vectors(self):
        """Calculate nutrient ratios for gap filling - lines ~400-420"""
        # Remove time dimension (select first/only time step)
        din_ann_NEWS = self.news_ds['NO3_CONC'].isel(time=0).values
        aa = np.where(din_ann_NEWS > 0)  # Get 2D indices where DIN > 0
        
        # Calculate DON (sum of LDON, SLDON, SRDON)
        temp1 = self.news_ds['LDON_CONC'].isel(time=0).values
        temp2 = self.news_ds['SLDON_CONC'].isel(time=0).values
        temp3 = self.news_ds['SRDON_CONC'].isel(time=0).values
        don_ann_NEWS = temp1 + temp2 + temp3
        
        # Calculate DOP (sum of LDOP, SLDOP, SRDOP)
        temp1 = self.news_ds['LDOP_CONC'].isel(time=0).values
        temp2 = self.news_ds['SLDOP_CONC'].isel(time=0).values
        temp3 = self.news_ds['SRDOP_CONC'].isel(time=0).values
        dop_ann_NEWS = temp1 + temp2 + temp3
        
        # Other nutrients
        pn_ann_NEWS = self.news_ds['NDET_CONC'].isel(time=0).values
        dip_ann_NEWS = self.news_ds['PO4_CONC'].isel(time=0).values
        pp_ann_NEWS = self.news_ds['PDET_CONC'].isel(time=0).values
        si_ann_NEWS = self.news_ds['SI_CONC'].isel(time=0).values
        
        # Calculate ratios AT valid indices (following original exactly)
        ratio_data = {
            'don_ratio': (['points'], don_ann_NEWS[aa] / din_ann_NEWS[aa]),
            'dip_ratio': (['points'], dip_ann_NEWS[aa] / din_ann_NEWS[aa]),
            'dop_ratio': (['points'], dop_ann_NEWS[aa] / din_ann_NEWS[aa]),
            'pn_ratio': (['points'], pn_ann_NEWS[aa] / din_ann_NEWS[aa]),
            'pp_ratio': (['points'], pp_ann_NEWS[aa] / din_ann_NEWS[aa]),
            'si_ratio': (['points'], si_ann_NEWS[aa] / din_ann_NEWS[aa])
        }
        
        # Create ratio dataset with points dimension (REMOVE DUPLICATE)
        n_points = len(din_ann_NEWS[aa])
        self.ratio_ds = xr.Dataset(
            ratio_data,
            coords={'points': range(n_points)}
        )
        
        # Store the 2D indices for potential future use
        self.valid_indices = aa
        
    def get_ratio_at_points(self, point_indices):
        """Get ratios at specific grid points for gap filling"""
        """
        point_indices: array of indices where we need ratios for gap filling
        Returns: xarray Dataset with ratios at those specific points
        """
        return self.ratio_ds.isel(points=point_indices)
    
    def get_all_ratios(self):
        """Return all calculated ratio vectors as xarray Dataset"""
        return self.ratio_ds
    
    def close(self):
        """Close the NetCDF dataset"""
        if self.news_ds is not None:
            self.news_ds.close()



class MonthlyRiverMapper:
    """Independent monthly river mapping using USGS data with GlobalNEWS2 gap filling"""
    
    def __init__(self, usgs_manager, glofas_manager, woa_manager, news_ratio_manager, 
                 q_min=0, min_dist=2.0, max_dist=2.0): #, inspect_map='n'):
        
        # Store data managers
        self.usgs_manager = usgs_manager
        self.glofas_manager = glofas_manager
        self.woa_manager = woa_manager
        self.news_ratio_manager = news_ratio_manager
        
        # Parameters
        self.q_min = q_min
        self.min_dist = min_dist
        self.max_dist = max_dist
        # self.inspect_map = inspect_map
        
        # Fractionation factors
        self.frac_ldon = 0.3
        self.frac_sldon = 0.35
        self.frac_srdon = 0.35
        self.frac_ldop = 0.3
        self.frac_sldop = 0.35
        self.frac_srdop = 0.35
        self.const_fed = 70.0e-6
        
        # Variables to process
        self.nutrient_vars = ['dic', 'alk', 'no3', 'nh4', 'din', 'don', 'pn', 
                             'dip', 'dop', 'pp', 'si', 'o2']
        
        # Results storage
        self.filtered_rivers_ds = None
        self.monthly_grid_ds = None
        self.final_ds = None
        
    def run_mapping(self):
        """Main workflow - filter, map, and process rivers"""
        print("Starting monthly river mapping...")
        
        # Step 1: Filter and sort rivers
        print("Step 1: Filtering rivers by region and discharge...")
        self._filter_and_sort_rivers()
        print(f"  → Filtered to {len(self.filtered_rivers_ds.station)} rivers")
        
        # Step 2: Map rivers to grid
        print("Step 2: Mapping rivers to runoff grid...")
        self._map_rivers_to_grid()
        print(f"  → Mapped to {len(self.monthly_grid_ds.points)} grid points")
        
        # Step 3: Create final dataset
        print("Step 3: Creating final monthly climatology...")
        self._create_final_dataset()
        print(f"  → Created dataset with {len(self.final_ds.data_vars)} variables")
        
        print("Monthly river mapping complete!")
        return self.final_ds
    
    def _filter_and_sort_rivers(self):
        """Filter rivers by region and discharge, then sort by flow"""
        # Get data from managers
        stations = self.usgs_manager.get_station_metadata()
        annual_data = self.usgs_manager.get_annual_data()
        monthly_data = self.usgs_manager.get_monthly_data()
        grid_data = self.glofas_manager.get_runoff_data()
        
        # Define regional bounds
        lon_bounds = (grid_data['lon'].min().values, grid_data['lon'].max().values)
        lat_bounds = (grid_data['lat'].min().values, grid_data['lat'].max().values)
        
        # Apply regional filter
        in_region = (
            (stations['mouth_lon'] >= lon_bounds[0]) & 
            (stations['mouth_lon'] <= lon_bounds[1]) &
            (stations['mouth_lat'] >= lat_bounds[0]) & 
            (stations['mouth_lat'] <= lat_bounds[1]) &
            np.isfinite(annual_data['Q']) & 
            (annual_data['Q'] > self.q_min)
        )
        
        # Filter datasets
        filtered_stations = stations.where(in_region, drop=True)
        filtered_annual = annual_data.where(in_region, drop=True)
        filtered_monthly = monthly_data.where(in_region, drop=True)
        
        # Sort by discharge (ascending - smallest first)
        sort_idx = np.argsort(filtered_annual['Q'].values)
        
        # Create combined filtered dataset
        self.filtered_rivers_ds = xr.Dataset({
            # Station metadata
            'river_name': filtered_stations['river_name'].isel(station=sort_idx),
            'mouth_lon': filtered_stations['mouth_lon'].isel(station=sort_idx),
            'mouth_lat': filtered_stations['mouth_lat'].isel(station=sort_idx),
            
            # Annual data
            'Q_ann': filtered_annual['Q'].isel(station=sort_idx),
            'dic_ann': filtered_annual['dic'].isel(station=sort_idx),
            'alk_ann': filtered_annual['alk'].isel(station=sort_idx),
            'din_ann': filtered_annual['din'].isel(station=sort_idx),
            'no3_ann': filtered_annual['no3'].isel(station=sort_idx),
            'nh4_ann': filtered_annual['nh4'].isel(station=sort_idx),
            'don_ann': filtered_annual['don'].isel(station=sort_idx),
            'pn_ann': filtered_annual['pn'].isel(station=sort_idx),
            'dip_ann': filtered_annual['dip'].isel(station=sort_idx),
            'dop_ann': filtered_annual['dop'].isel(station=sort_idx),
            'pp_ann': filtered_annual['pp'].isel(station=sort_idx),
            'si_ann': filtered_annual['si'].isel(station=sort_idx),
            'o2_ann': filtered_annual['o2'].isel(station=sort_idx),
            
            # Monthly data
            'Q_monthly': filtered_monthly['Q'].isel(station=sort_idx),
            'dic_monthly': filtered_monthly['dic'].isel(station=sort_idx),
            'alk_monthly': filtered_monthly['alk'].isel(station=sort_idx),
            'din_monthly': filtered_monthly['din'].isel(station=sort_idx),
            'no3_monthly': filtered_monthly['no3'].isel(station=sort_idx),
            'nh4_monthly': filtered_monthly['nh4'].isel(station=sort_idx),
            'don_monthly': filtered_monthly['don'].isel(station=sort_idx),
            'pn_monthly': filtered_monthly['pn'].isel(station=sort_idx),
            'dip_monthly': filtered_monthly['dip'].isel(station=sort_idx),
            'dop_monthly': filtered_monthly['dop'].isel(station=sort_idx),
            'pp_monthly': filtered_monthly['pp'].isel(station=sort_idx),
            'si_monthly': filtered_monthly['si'].isel(station=sort_idx),
            'o2_monthly': filtered_monthly['o2'].isel(station=sort_idx),
        })
    
    def _map_rivers_to_grid(self):
        """Map each river to runoff grid points"""
        runoff_vectors = self.glofas_manager.get_runoff_vectors()
        news_ratios = self.news_ratio_manager.get_all_ratios()
        
        # Initialize monthly concentration arrays
        n_months = 12
        n_points = len(runoff_vectors.points)
        n_rivers = len(self.filtered_rivers_ds.station)
        
        # Create monthly concentration dataset
        monthly_data = {}
        for var in self.nutrient_vars:
            monthly_data[var] = (['month', 'points'], np.zeros((n_months, n_points)))
        
        self.monthly_grid_ds = xr.Dataset(
            monthly_data,
            coords={
                'month': range(n_months),
                'points': runoff_vectors.points,
                'lon': (['points'], runoff_vectors['lon'].values),
                'lat': (['points'], runoff_vectors['lat'].values),
                'Q_ann': (['points'], runoff_vectors['Q_ann'].values)
            }
        )
        
        # Map each river
        for k in range(n_rivers):
            river_name = self.filtered_rivers_ds['river_name'][k].values
            print(f"  Processing river {k+1}/{n_rivers}: {river_name}")
            self._map_single_river(k, runoff_vectors, news_ratios)
        
        # Apply gap filling
        self._apply_gap_filling(runoff_vectors)
    
    def _map_single_river(self, river_idx, runoff_vectors, news_ratios):
        """Map single river using distance and discharge matching"""
        river = self.filtered_rivers_ds.isel(station=river_idx)
        
        # Calculate distances from river mouth to all runoff points
        distances = cdist(
            [[river['mouth_lon'].values, river['mouth_lat'].values]],
            np.column_stack([runoff_vectors['lon'].values, runoff_vectors['lat'].values])
        )[0]
        
        # Sort by distance
        dist_sort_idx = np.argsort(distances)
        
        # Check if closest point is within minimum distance
        if distances[dist_sort_idx[0]] >= self.min_dist:
            # Add plotting for outside rivers
            if self.inspect_map == 'y':
                self._plot_outside_river(river, runoff_vectors)
            return
        
        # Match discharge by accumulating closest grid points
        target_discharge = river['Q_ann'].values
        Q_sum1 = 0.0
        Q_sum2 = 0.0
        n = 0
        
        # Add max_dist check like original
        while (Q_sum2 < target_discharge and 
            n < len(dist_sort_idx) - 1 and
            distances[dist_sort_idx[n + 1]] < self.max_dist):  # ADDED MAX_DIST CHECK
            Q_sum1 = Q_sum2
            n += 1
            Q_sum2 = Q_sum1 + runoff_vectors['Q_ann'].values[dist_sort_idx[n]]
        
        # Choose optimal number of points
        if abs(Q_sum1 - target_discharge) < abs(Q_sum2 - target_discharge):
            nrp = n - 1
            matched_discharge = Q_sum1
        else:
            nrp = n
            matched_discharge = Q_sum2
        
        if nrp <= 0:
            return
            
        # Selected grid points for this river
        selected_indices = dist_sort_idx[:nrp]
        
        # # Add inspection plotting
        # if self.inspect_map == 'y':
        #     self._plot_river_mapping(river, selected_indices, runoff_vectors, matched_discharge, target_discharge)
        
        # Assign monthly concentrations
        self._assign_monthly_concentrations(river, selected_indices, news_ratios)
        
    def _assign_monthly_concentrations(self, river, selected_indices, news_ratios):
        """Assign monthly concentration values with enhanced gap filling logic"""
        for m in range(12):
            # DIC and ALK
            for var in ['dic', 'alk']:
                monthly_val = river[f'{var}_monthly'][m].values
                annual_val = river[f'{var}_ann'].values
                
                if not np.isnan(monthly_val): 
                    value = monthly_val
                elif np.isfinite(annual_val): 
                    value = annual_val
                else:
                    value = 0.0
                    
                self.monthly_grid_ds[var][m, selected_indices] = value
            
            # DIN-based processing
            din_monthly = river['din_monthly'][m].values
            din_annual = river['din_ann'].values
            
            if np.isfinite(din_monthly) or np.isfinite(din_annual):
                # Use monthly if available, otherwise annual
                din_value = din_monthly if np.isfinite(din_monthly) else din_annual
                self.monthly_grid_ds['din'][m, selected_indices] = din_value
                
                # NO3
                no3_monthly = river['no3_monthly'][m].values
                no3_annual = river['no3_ann'].values
                if np.isfinite(no3_monthly):
                    self.monthly_grid_ds['no3'][m, selected_indices] = no3_monthly
                elif np.isfinite(no3_annual):
                    self.monthly_grid_ds['no3'][m, selected_indices] = no3_annual
                
                # NH4
                nh4_monthly = river['nh4_monthly'][m].values
                nh4_annual = river['nh4_ann'].values
                if np.isfinite(nh4_monthly):
                    self.monthly_grid_ds['nh4'][m, selected_indices] = nh4_monthly
                elif np.isfinite(nh4_annual):
                    self.monthly_grid_ds['nh4'][m, selected_indices] = nh4_annual
                
                # DON with NEWS2 gap filling
                don_monthly = river['don_monthly'][m].values
                don_annual = river['don_ann'].values
                
                if np.isnan(don_monthly) and np.isnan(don_annual):
                    ratio_values = self._get_news_ratios_at_points(selected_indices, news_ratios, 'don_ratio')
                    self.monthly_grid_ds['don'][m, selected_indices] = din_value * ratio_values
                elif np.isnan(don_monthly) and not np.isnan(don_annual):
                    self.monthly_grid_ds['don'][m, selected_indices] = don_annual
                else:
                    self.monthly_grid_ds['don'][m, selected_indices] = don_monthly
                
                # PN with NEWS2 gap filling
                pn_monthly = river['pn_monthly'][m].values
                pn_annual = river['pn_ann'].values
                
                if np.isnan(pn_monthly) and np.isnan(pn_annual):
                    ratio_values = self._get_news_ratios_at_points(selected_indices, news_ratios, 'pn_ratio')
                    self.monthly_grid_ds['pn'][m, selected_indices] = din_value * ratio_values
                elif np.isnan(pn_monthly) and not np.isnan(pn_annual):
                    self.monthly_grid_ds['pn'][m, selected_indices] = pn_annual
                else:
                    self.monthly_grid_ds['pn'][m, selected_indices] = pn_monthly
                
                # DIP with NEWS2 gap filling
                dip_monthly = river['dip_monthly'][m].values
                dip_annual = river['dip_ann'].values
                
                if np.isnan(dip_monthly) and np.isnan(dip_annual):
                    ratio_values = self._get_news_ratios_at_points(selected_indices, news_ratios, 'dip_ratio')
                    self.monthly_grid_ds['dip'][m, selected_indices] = din_value * ratio_values
                elif np.isnan(dip_monthly) and not np.isnan(dip_annual):
                    self.monthly_grid_ds['dip'][m, selected_indices] = dip_annual
                else:
                    self.monthly_grid_ds['dip'][m, selected_indices] = dip_monthly
                
                # DOP with NEWS2 gap filling
                dop_monthly = river['dop_monthly'][m].values
                dop_annual = river['dop_ann'].values
                
                if np.isnan(dop_monthly) and np.isnan(dop_annual):
                    ratio_values = self._get_news_ratios_at_points(selected_indices, news_ratios, 'dop_ratio')
                    self.monthly_grid_ds['dop'][m, selected_indices] = din_value * ratio_values
                elif np.isnan(dop_monthly) and not np.isnan(dop_annual):
                    self.monthly_grid_ds['dop'][m, selected_indices] = dop_annual
                else:
                    self.monthly_grid_ds['dop'][m, selected_indices] = dop_monthly
                
                # PP with NEWS2 gap filling
                pp_monthly = river['pp_monthly'][m].values
                pp_annual = river['pp_ann'].values
                
                if np.isnan(pp_monthly) and np.isnan(pp_annual):
                    ratio_values = self._get_news_ratios_at_points(selected_indices, news_ratios, 'pp_ratio')
                    self.monthly_grid_ds['pp'][m, selected_indices] = din_value * ratio_values
                elif np.isnan(pp_monthly) and not np.isnan(pp_annual):
                    self.monthly_grid_ds['pp'][m, selected_indices] = pp_annual
                else:
                    self.monthly_grid_ds['pp'][m, selected_indices] = pp_monthly
                
                # SI with NEWS2 gap filling
                si_monthly = river['si_monthly'][m].values
                si_annual = river['si_ann'].values
                
                if np.isnan(si_monthly) and np.isnan(si_annual):
                    ratio_values = self._get_news_ratios_at_points(selected_indices, news_ratios, 'si_ratio')
                    self.monthly_grid_ds['si'][m, selected_indices] = din_value * ratio_values
                elif np.isnan(si_monthly) and not np.isnan(si_annual):
                    self.monthly_grid_ds['si'][m, selected_indices] = si_annual
                else:
                    self.monthly_grid_ds['si'][m, selected_indices] = si_monthly
            
            # O2 with WOA gap filling
            o2_monthly = river['o2_monthly'][m].values
            
            if np.isnan(o2_monthly):
                # Use WOA O2 saturation
                o2_values = self._get_woa_o2_at_points(selected_indices, m)
                self.monthly_grid_ds['o2'][m, selected_indices] = o2_values
            else:
                self.monthly_grid_ds['o2'][m, selected_indices] = o2_monthly

    def _get_woa_o2_at_points(self, selected_indices, month):
        """Get WOA O2 saturation at specific runoff points"""
        # Similar coordinate matching approach for O2 data
        if hasattr(self.woa_manager, 'o2sat_woa_monthly_vecs'):
            return self.woa_manager.o2sat_woa_monthly_vecs[month, selected_indices]
        else:
            # Calculate default O2 saturation if not available
            return np.full(len(selected_indices), 300.0)  # Default O2 saturation
            
    def _apply_gap_filling(self, runoff_vectors):
        """Fill remaining zeros with nearest neighbor interpolation"""
        print("  Applying gap filling to unmapped grid points...")
        
        for m in range(12):
            for var in self.nutrient_vars:
                conc_values = self.monthly_grid_ds[var][m, :].values
                
                # Find zero and non-zero points
                zero_idx = np.where(conc_values == 0)[0]
                nonzero_idx = np.where(conc_values > 0)[0]
                
                if len(zero_idx) > 0 and len(nonzero_idx) > 0:
                    # Nearest neighbor interpolation
                    filled_values = griddata(
                        (runoff_vectors['lon'].values[nonzero_idx], 
                         runoff_vectors['lat'].values[nonzero_idx]),
                        conc_values[nonzero_idx],
                        (runoff_vectors['lon'].values[zero_idx],
                         runoff_vectors['lat'].values[zero_idx]),
                        method='nearest'
                    )
                    
                    self.monthly_grid_ds[var][m, zero_idx] = filled_values
    
    def _create_final_dataset(self):
        """Create final monthly climatology dataset with fractionation"""
        # Get grid data
        grid_data = self.glofas_manager.get_runoff_data()
        runoff_vectors = self.glofas_manager.get_runoff_vectors()
        
        # Map vectors back to 2D grid
        nlat, nlon = grid_data['lat'].shape
        n_months = 12
        
        # Get indices where runoff > 0
        runoff_mask = grid_data['Q_ann'] > 0
        runoff_indices = np.where(runoff_mask)
        
        # Initialize final dataset
        self.final_ds = xr.Dataset(coords={
            'time': range(n_months),
            'y': range(nlat),
            'x': range(nlon),
            'lat': (['y', 'x'], grid_data['lat'].values),
            'lon': (['y', 'x'], grid_data['lon'].values)
        })
        unit_scale = 1.0 / 1000.0
        
        # Create final variables with proper units and fractionation
        for m in range(n_months):
            # Initialize monthly fields
            monthly_fields = {}
            for var in self.nutrient_vars:
                field_2d = np.zeros((nlat, nlon))
                field_2d[runoff_indices] = self.monthly_grid_ds[var][m, :].values
                monthly_fields[var] = field_2d * unit_scale
            
            # Apply fractionation and unit conversion (mmol/m³ to mol/m³)
            
            # Create variables for first month, then assign all months
            if m == 0:
                # Carbon system
                self.final_ds['DIC_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                           {'units': 'mol m-3', 'long_name': 'Dissolved Inorganic Carbon'})
                self.final_ds['ALK_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                           {'units': 'eq m-3', 'long_name': 'Total Alkalinity'})
                
                # Nitrogen species
                self.final_ds['NO3_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                           {'units': 'mol m-3', 'long_name': 'Nitrate'})
                self.final_ds['NH4_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                           {'units': 'mol m-3', 'long_name': 'Ammonium'})
                self.final_ds['LDON_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                            {'units': 'mol m-3', 'long_name': f'{self.frac_ldon}*DON_CONC'})
                self.final_ds['SLDON_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                             {'units': 'mol m-3', 'long_name': f'{self.frac_sldon}*DON_CONC'})
                self.final_ds['SRDON_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                             {'units': 'mol m-3', 'long_name': f'{self.frac_srdon}*DON_CONC'})
                self.final_ds['NDET_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                            {'units': 'mol m-3', 'long_name': 'Nitrogen Detritus'})
                
                # Phosphorus species
                self.final_ds['PO4_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                           {'units': 'mol m-3', 'long_name': 'Phosphate'})
                self.final_ds['LDOP_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                            {'units': 'mol m-3', 'long_name': f'{self.frac_ldop}*DOP_CONC'})
                self.final_ds['SLDOP_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                             {'units': 'mol m-3', 'long_name': f'{self.frac_sldop}*DOP_CONC'})
                self.final_ds['SRDOP_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                             {'units': 'mol m-3', 'long_name': f'{self.frac_srdop}*DOP_CONC'})
                self.final_ds['PDET_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                            {'units': 'mol m-3', 'long_name': 'Phosphorus Detritus'})
                
                # Other nutrients
                self.final_ds['SI_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                          {'units': 'mol m-3', 'long_name': 'Silicate'})
                self.final_ds['O2_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                          {'units': 'mol m-3', 'long_name': 'Dissolved Oxygen'})
                
                # Iron species
                self.final_ds['FED_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                           {'units': 'mol m-3', 'long_name': 'Dissolved Iron'})
                self.final_ds['FEDET_CONC'] = (['time', 'y', 'x'], np.zeros((n_months, nlat, nlon)),
                                             {'units': 'mol m-3', 'long_name': 'Iron Detritus'})
            
            # Assign monthly values
            self.final_ds['DIC_CONC'][m, :, :] = monthly_fields['dic']
            self.final_ds['ALK_CONC'][m, :, :] = monthly_fields['alk']
            self.final_ds['NO3_CONC'][m, :, :] = monthly_fields['no3']
            self.final_ds['NH4_CONC'][m, :, :] = monthly_fields['nh4']
            self.final_ds['LDON_CONC'][m, :, :] = self.frac_ldon * monthly_fields['don']
            self.final_ds['SLDON_CONC'][m, :, :] = self.frac_sldon * monthly_fields['don']
            self.final_ds['SRDON_CONC'][m, :, :] = self.frac_srdon * monthly_fields['don']
            self.final_ds['NDET_CONC'][m, :, :] = monthly_fields['pn']
            self.final_ds['PO4_CONC'][m, :, :] = monthly_fields['dip']
            self.final_ds['LDOP_CONC'][m, :, :] = self.frac_ldop * monthly_fields['dop']
            self.final_ds['SLDOP_CONC'][m, :, :] = self.frac_sldop * monthly_fields['dop']
            self.final_ds['SRDOP_CONC'][m, :, :] = self.frac_srdop * monthly_fields['dop']
            self.final_ds['PDET_CONC'][m, :, :] = monthly_fields['pp']
            self.final_ds['SI_CONC'][m, :, :] = monthly_fields['si']
            self.final_ds['O2_CONC'][m, :, :] = monthly_fields['o2']
            
            # Iron concentrations
            fed_field = monthly_fields['no3'].copy()
            fedet_field = monthly_fields['no3'].copy()
            fed_field[fed_field > 0] = self.const_fed
            fedet_field[fedet_field > 0] = 0.0
            
            self.final_ds['FED_CONC'][m, :, :] = fed_field
            self.final_ds['FEDET_CONC'][m, :, :] = fedet_field
    
    def export_to_netcdf(self, output_file):
        """Export monthly climatology to NetCDF with proper time coordinates"""
        if self.final_ds is None:
            raise ValueError("No dataset to export. Run run_mapping() first.")
        
        # Set up monthly climatology time coordinates
        dates = [
            [1993, 1, 16, 12, 0, 0], [1993, 2, 15, 0, 0, 0], [1993, 3, 16, 12, 0, 0],
            [1993, 4, 16, 0, 0, 0], [1993, 5, 16, 12, 0, 0], [1993, 6, 16, 0, 0, 0],
            [1993, 7, 16, 12, 0, 0], [1993, 8, 16, 12, 0, 0], [1993, 9, 16, 0, 0, 0],
            [1993, 10, 16, 12, 0, 0], [1993, 11, 16, 0, 0, 0], [1993, 12, 16, 12, 0, 0]
        ]
        
        time_coords = [cftime.datetime(*date, calendar='365_day') for date in dates]
        
        # Assign time coordinates
        ds_export = self.final_ds.assign_coords(time=time_coords)
        ds_export.time.attrs['cartesian_axis'] = 'T'
        ds_export.time.attrs['modulo'] = ' '  # Important for climatology
        
        # Add global attributes
        ds_export.attrs.update({
            'title': 'Monthly River Nutrient Climatology',
            'description': 'Monthly climatology of river nutrient concentrations mapped to model grid',
            'source': 'USGS river chemistry data with GlobalNEWS2 gap filling',
            'created_by': 'MonthlyRiverMapper'
        })
        
        # Export to NetCDF
        ds_export.to_netcdf(
            output_file,
            format='NETCDF4',
            engine='netcdf4',
            encoding={
                'time': {
                    'units': 'days since 0001-01-01T00:00:00',
                    'calendar': '365_day'
                }
            },
            unlimited_dims='time'
        )
        
        print(f"Monthly climatology exported to: {output_file}")
    

    def _get_news_ratios_at_points(self, selected_indices, news_ratios, ratio_name):
        """Get NEWS2 ratios at specific runoff points using coordinate matching"""
        runoff_vectors = self.glofas_manager.get_runoff_vectors()
        
        # Get coordinates of selected runoff points
        selected_lons = runoff_vectors['lon'].values[selected_indices]
        selected_lats = runoff_vectors['lat'].values[selected_indices]
        
        # Find nearest NEWS2 ratio points using nearest neighbor
        from scipy.spatial.distance import cdist
        
        # Get NEWS2 coordinates (need to reconstruct from valid_indices)
        news_manager = self.news_ratio_manager
        news_grid = news_manager.news_ds.isel(time=0)
        
        # Get coordinates where NEWS2 ratios are valid
        news_valid_coords = np.column_stack([
            news_grid['lon'].values[news_manager.valid_indices],
            news_grid['lat'].values[news_manager.valid_indices]
        ])
        
        # Find nearest NEWS2 points for each selected runoff point
        selected_coords = np.column_stack([selected_lons, selected_lats])
        distances = cdist(selected_coords, news_valid_coords)
        nearest_news_indices = np.argmin(distances, axis=1)
        
        # Get ratio values at nearest points
        ratio_values = news_ratios[ratio_name].values[nearest_news_indices]
        
        return ratio_values

    # Getter methods
    def get_filtered_rivers(self):
        return self.filtered_rivers_ds
    
    def get_monthly_grid_data(self):
        return self.monthly_grid_ds
    
    def get_final_dataset(self):
        return self.final_ds

def load_config(config_file:str):
    """Load configuration from YAML file"""
    with open(config_file, 'r') as stream:
        return yaml.safe_load(stream)


if __name__ =='__main__':



    # Define coordinate corrections
    corrections = {
        'Susquehanna': {'lat': 38.5, 'lon': -77.5},
        'Delaware': {'lat': 39.5, 'lon': -75.5},
        'Potomac': {'lat': 38.5, 'lon': -77.5},
        'Mississippi': {'lat': 29.25, 'lon': -89.25},
        'Alabama': {'lat': 30.5, 'lon': None}
    }


    parser = argparse.ArgumentParser(description='combine river bgc files')
    parser.add_argument('--config', type=str,
                        help='YAML configuration file path')
    args = parser.parse_args()


    config = load_config(args.config)
    config = config['bgc_processor_combine']


    filename_chem      = config['filename_chem']
    filename_discharge = config['filename_discharge']
    basin_file         = config['basin_file']
    news_file          = config['news_file']
    grid_file          = config['grid_file']
    woa_temp_pattern   = config['woa_temp_pattern']
    runoff_file        = config['runoff_file']
    output_file        = config['output_file']
    print(filename_chem)

    # filename_chem = 'data/input/mclim_19902022_chem.nc'
    # filename_discharge = 'data/input/mclim_19902022_disc.nc'
    # basin_file='data/input/GlobalNEWS2_RH2000Dataset-version1.0.xls'
    # news_file='data/tmp/RiverNutrients_GlobalNEWS2_plusFe_Q100_GLOFAS_NWA12.nc'
    # grid_file='data/input/19930101.ocean_static.nc'
    # woa_temp_pattern='data/input/woa18_decav_t*_04.nc'
    # runoff_file='data/output/glofas_runoff_mean.nc'
    # output_file='data/output/RiverNutrients_Integrated_NWA12_GLOFAS_RC4US1990to2022_2023_04_v2.nc'

    usgs_manager = USGSDataManager(filename_chem, filename_discharge, nutrient_option=2)
    usgs_manager.apply_coordinate_corrections(corrections)
    
    glofas_manager = GlofasRunoffManager(runoff_file)
    
    # # Get full runoff data
    runoff_data = glofas_manager.get_runoff_data()
    runoff_vectors = glofas_manager.get_runoff_vectors()

    # # interpolate woa onto mom6 grid
    woa_manager = WOATemperatureManager(woa_temp_pattern, grid_file)
    woa_manager.interpolate_to_grid()
    woa_manager.calculate_o2_saturation_vectors(runoff_vectors)  # Use with runoff vectors

    
    # # file saved created by bgc_processor.py
    news_ratio_manager = GlobalNEWS2RatioManager(news_file)
    news_ratio_manager.load_data()
    news_ratio_manager.calculate_ratio_vectors()
    
    # # Get all ratios
    ratios = news_ratio_manager.get_all_ratios()
    print(f"Available ratios: {list(ratios.data_vars.keys())}")
    # print(f"DON ratio shape: {ratios['don_ratio'].shape}")
    # print(f"Total valid points: {len(ratios.points)}")


    monthly_mapper = MonthlyRiverMapper(
         usgs_manager=usgs_manager,
         glofas_manager=glofas_manager,
         woa_manager=woa_manager,
         news_ratio_manager=news_ratio_manager,
         q_min=0,           # Minimum discharge threshold (m3/s)
         min_dist=2.0,      # Minimum distance for river assignment (degrees)
         max_dist=2.0,      # Maximum distance to search for grid points (degrees)
         # inspect_map='y'    # Set to 'y' for debugging plots
     )

     # Step 3: Run the mapping process
    print("Running mapping process...")
    final_dataset = monthly_mapper.run_mapping()
    monthly_mapper.export_to_netcdf(output_file)
    
