import cdsapi
import os
import xarray as xr
import glob 
import os.path as osp
import sys
import os
os.environ['CDSAPI_URL'] = 'https://cds.climate.copernicus.eu/api'


def download_era5(year, variables):
    c = cdsapi.Client()

    for v in variables:
        outfile = f"ERA5_{v}_{year}.grib"
        if not glob.glob(outfile):
            dataset = "reanalysis-era5-single-levels"

            request = {
                    "product_type": ["reanalysis"],
                    "variable": [variables[v]],
                    "year": [f"{year}"],
                    "month": [f"{i:02d}" for i in range(1,13)],
                    "day":   [f"{i:02d}" for i in range(1,32)],
                    "time":  [f"{i:02d}" for i in range(24)],
                    "area": [54, -80, 25,-55,],
                    "data_format": "grib",
                    #"download_format": "unarchived"
                }
            print(request)
            c.retrieve(dataset, request,outfile)           

def create_directory_if_not_exists(directory_path):
    os.makedirs(directory_path, exist_ok=True)
    if not os.path.exists(directory_path):
        print(f"Directory '{directory_path}' created successfully.")
    else:
        print(f"Directory '{directory_path}' already exists.")


def grib2nc():

    fpath = glob.glob("*grib")
    
    for f in fpath:
        outfile = osp.basename(f).replace(".grib","") + ".nc"
        if not glob.glob(outfile):
            with xr.open_dataset(f, engine="cfgrib") as ds:
                ds.to_netcdf(outfile, unlimited_dims="time")

def main(variables, start_year, end_year):
    for year in range(start_year,end_year):
        # Specify the directory path you want to check/create
        directory_path = str(year)
        # Call the function to create the directory if needed
        create_directory_if_not_exists(directory_path)
        os.chdir(directory_path)
        download_era5(year, variables)
        grib2nc()
        os.chdir("..")


if __name__ == "__main__":

    start_year = 1993
    end_year   = 1994

    # start_year = int(sys.argv[1])
    # end_year = int(sys.argv[2])
    
    variables = {
        "metss": "mean_eastward_turbulent_surface_stress",
        "sst":   "sea_surface_temperature",
        "mntss": "mean_northward_turbulent_surface_stress",
        "mslhf": "mean_surface_latent_heat_flux",
        "msshf": "mean_surface_sensible_heat_flux",
        "e"    : "evaporation",
        "tp":"total_precipitation",
        "msl": "mean_sea_level_pressure",
        "sf": "snowfall",
        "sp": "surface_pressure",
        "d2m": "2m_dewpoint_temperature",
        "ssrd": "surface_solar_radiation_downwards",
        "strd": "surface_thermal_radiation_downwards",
        "t2m": "2m_temperature",
        "u10": "10m_u_component_of_wind",
        "v10": "10m_v_component_of_wind",
    }

    main(variables, start_year, end_year)
