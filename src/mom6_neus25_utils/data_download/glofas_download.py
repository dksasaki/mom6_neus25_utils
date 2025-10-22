import cdsapi
import os
import xarray as xr
import glob 
import os.path as osp
import sys

def download_glofas(year):
    c = cdsapi.Client()


    outfile = f'glovas_{year}.grib'
    if not glob.glob(outfile):
        dataset = "cems-glofas-historical"
        request = {
                "system_version": ["version_3_1"],
                "hydrological_model": ["lisflood"],
                "product_type": ["consolidated"],
                "variable": ["river_discharge_in_the_last_24_hours"],
                "hyear": [f"{year}"],
                "hmonth": [f"{i:02d}" for i in range(1,13)],
                "hday":   [f"{i:02d}" for i in range(1,32)],
                "area": [54, -80, 25, -55,],
                "data_format": "grib",
                "download_format": "unarchived"
            }

        c.retrieve(dataset,
            request = request,
            target = outfile)
            

def create_directory_if_not_exists(directory_path):
    os.makedirs(directory_path, exist_ok=True)
    if not os.path.exists(directory_path):
        print(f"Directory '{directory_path}' created successfully.")
    else:
        print(f"Directory '{directory_path}' already exists.")


def grib2nc():

    fpath = glob.glob('*grib')
    
    for f in fpath:
        outfile = osp.basename(f).replace('.grib','') + '.nc'
        if not glob.glob(outfile):
            with xr.open_dataset(f, engine='cfgrib') as ds:
                ds.to_netcdf(outfile, unlimited_dims='time')

def main(start_year, end_year):
    for year in range(start_year,end_year):
        # Specify the directory path you want to check/create
        directory_path = str(year)
        # Call the function to create the directory if needed
        create_directory_if_not_exists(directory_path)
        os.chdir(directory_path)
        download_glofas(year)
        grib2nc()
        os.chdir('..')


if __name__ == '__main__':

    start_year = 1993
    end_year   = 1994

    # start_year = int(sys.argv[1])
    # end_year = int(sys.argv[2])


    main(start_year, end_year)

