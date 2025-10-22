import os
import pandas as pd
import copernicusmarine as cm
import datetime as dtt
import glob
import sys

def create_directory_if_not_exists(directory_path):
    os.makedirs(directory_path, exist_ok=True)
    if not os.path.exists(directory_path):
        print(f"Directory '{directory_path}' created successfully.")
    else:
        print(f"Directory '{directory_path}' already exists.")

def download_glorys(variables, start_date, end_date, test=False):
    times = pd.date_range(start=start_date, end=end_date, freq='1D')

    id="cmems_mod_glo_phy_my_0.083deg_P1D-m"
    force_dataset_version="202311"
    variables=variables
    minlon=-83
    maxlon=-54.5
    minlat=30.5
    maxlat=49.5
    mindep=0.49402499198913574
    maxdep=5727.9#5274.781

    file = lambda x: f'{id}_multi-vars_{abs(minlon):0.2f}W-{abs(maxlon):0.2f}W_' + \
                    f'{abs(minlat):0.2f}N-{abs(maxlat):0.2f}N_{mindep:0.2f}-5274.78m_{x}.nc'

    
    for t in times:
        taux = dtt.datetime.strftime(t,'%Y-%m-%d')
        if glob.glob(file(taux)):
            print(f'{taux} exists')
            continue
        else:
            pass

        print(f'download {glob.glob(file(taux))}')
        if test:
            continue
        else:
            cm.subset(
               dataset_id=id,
               dataset_version="202311",
               variables=variables,
               minimum_longitude=minlon,
               maximum_longitude=maxlon,
               minimum_latitude=minlat,
               maximum_latitude=maxlat,
               start_datetime=str(t),
               end_datetime=str(t),
               minimum_depth=mindep,
               maximum_depth=maxdep,
               username='dsasaki2',
               password='p3nU56rbpQjuCP-',
               force_download=True,
               overwrite_output_data=False
           )



if __name__ == '__main__':
    start_year = int(sys.argv[1])
    end_year   = int(sys.argv[2])
    # start_year = 1994
    # end_year   = 2020
    variables_glorys = ["bottomT", "mlotst", "siconc", "so", "thetao", "uo", "vo", "zos"]


    for year in range(start_year,end_year):
        directory_path = str(year)
        
        start_date_glorys = f"{year}-01-01T00:00:00"
        end_date_glorys = f"{year+1}-01-01T00:00:00"
        
        create_directory_if_not_exists(directory_path)
        os.chdir(directory_path)
        download_glorys(variables_glorys, start_date_glorys, end_date_glorys)
        os.chdir('..')
