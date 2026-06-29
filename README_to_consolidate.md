
# creating boundary conditions through PyESPER

``` bash
your_working_dir=/some/path    # where to clone CEFI

# if you want to start from scratch do git clone
# git clone https://github.com/dksasaki/mom6_neus25_utils.git

cd $your_working_dir

git fetch # use this if you repository already exists
git pull 
#git checkout feature_bgc_obc
git checkout feature_bgc_obc_consolidate

# getting corrections from upstream PyESPER dependency in an existing cloned mom6_neus25_utils.got
# pixi update --no-install pyesper && pixi install

pixi install



cd src/mom6_neus25_utils/boundary

# there should be a couple of bgc scripts for obc
# bgc_obc_processor_esper.py
# bgc_obc_processor.py
# both of of them use a config.yaml (see below)
```


`bgc_obc_processor.py` : allow us to get information from WOA and a climatology from a global MOM6-COBALT run that was published by Charlie Stock (2014)
bgc_obc_processor_esper.py : calculates different nutrient based on existing temperature or salinity boundary conditions (obtained through glorys_obc_processor.py)

After editing your yaml file (see below), log into a node, change the paths/dates/year according to your experiments. Than you can run the scripts with python and the should output the files where you requested. 
For the esper  script: it takes a while to run a single time, for now we are able to run a single year at a time.

Below are the instructions to use config.yaml from the command line

```
python bgc_obc_processor.py --config config.yaml
python bgc_obc_processor_esper.py --config config.yaml

#optional usage
python bgc_obc_processor_esper.py --config config.yaml --year 2003
```


You'll need a yaml file config.yaml in the boundary directory:
```yaml
boundary:
  output_dir: './outputs'
  cache: '.' 
  grid_file: '/projects/schultz/d.sasaki/experiments/v1.1_simulation/tools_and_data/data/source/ocean_hgrid.nc'  # no need to change this in explorer
  time0: 1993-01-01
  last_time: 2022-01-01
  woa_file: '/home/d.sasaki/schultz/data/woa/woa18_NEUS/woa_ann_merged.nc'   # no need to change this in explorer
  cobalt_file: '/projects/schultz/data/cobalt_global/ocean_cobalt_tracers.1988-2007.ann.nc'   # no need to change this in explorer
  TS_bound_file: '/home/d.sasaki/schultz/d.sasaki/experiments/v1.0_simulation/nwa25/INPUT/'   # no need to change this in explorer
  topog_file:    '/home/d.sasaki/schultz/d.sasaki/experiments/v1.0_simulation/nwa25/INPUT/'   # no need to change this in explorer
 
  segments:
    - id: 1
      border: 'south'
    - id: 2
      border: 'east'

esper:
  pyesper_path: '/home/d.sasaki/packages/esper_env/PyESPER'   # shouldn't need to change this in explorer
  output_path: './outputs/esper'
  TS_bound_files_path: '/projects/schultz/d.sasaki/experiments/v1.0_simulation/nwa25/INPUT/'    # no need to change this in explorer
  topog_file: '/projects/schultz/d.sasaki/experiments/v1.0_simulation/nwa25/INPUT/ocean_topog.nc'    # no need to change this in explorer
  year: 1993
  average_type: 'resample'
  clim_type: 'season'
  nutrients:  # what nutrients are outputted from T,S
    #- nitrate
    #- silicate
    #- phosphate
    - TA
    - DIC
  nutrient_mapping:
    nitrate: no3
    silicate: sio4
    phosphate: po4
    DIC: dissic
    TA: talk
```



Using SLURM job array to run ESPER for multiple years

```bash
#!/bin/bash
#SBATCH --job-name=esper_bc
#SBATCH --time=04:00:00
#SBATCH --ntasks=5
#SBATCH --partition=short 
#SBATCH --mem=16G
#SBATCH --output=esper_bc_%A_%a.out
#SBATCH --array=0-4

START_YEAR=2018 

YEAR=$((START_YEAR + SLURM_ARRAY_TASK_ID))

python bgc_obc_processor_esper.py --config config.yaml --year $YEAR
```


Using SLURM job array to merge boundary conditions together

```bash
#!/bin/bash
#SBATCH --job-name=merge_bc
#SBATCH --time=04:00:00
#SBATCH --partition=short 
#SBATCH --mem=16G
#SBATCH --ntasks=5
#SBATCH --array=0-4
#SBATCH --output=merge_bc_%A_%a.out

START_YEAR=2018 

YEAR=$((START_YEAR + SLURM_ARRAY_TASK_ID))

echo "Processing year ${YEAR}"

ncks -6 -A DIC_001_${YEAR}.nc bgc_DIC_${YEAR}.nc
ncks -6 -A DIC_002_${YEAR}.nc bgc_DIC_${YEAR}.nc

ncks -6 -A TA_001_${YEAR}.nc bgc_TA_${YEAR}.nc
ncks -6 -A TA_002_${YEAR}.nc bgc_TA_${YEAR}.nc
```

