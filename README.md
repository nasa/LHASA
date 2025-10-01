# Landslide Hazard Assessment for Situational Awareness

LHASA was created at Goddard Space Flight Center to identify the potential for rainfall-triggered landslides in real time. 

## What's new

See the [Changelog](https://github.com/nasa/LHASA/blob/master/CHANGELOG.md)

## LHASA 2.1.1

LHASA version 2 adopts machine learning to estimate the probability of landslide occurrence at a 30-arcsecond (~1 km) daily resolution. In addition, it estimates the potential exposure of human population and roads to landslide hazard and maps the basins likely to experience post-fire debris flows. 

### Real-time data availability
The latest predictions can be downloaded from https://maps.nccs.nasa.gov/download/landslides. It can also be accessed as an ArcGIS web map at https://landslides.nasa.gov/viewer. NASA provides these data on a best-effort basis, typically four times each day, but with frequent server downtime. Users requiring a fully operational system are encouraged to clone this repository and run LHASA at the desired cadence.

### Data files

LHASA requires several large data files, but not all data may be needed by all users. The contents of [static.zip](https://gpm.nasa.gov/sites/default/files/data/landslides/static.zip) are required for the global landslide forecast. The contents of [exposure.zip](https://gpm.nasa.gov/sites/default/files/data/landslides/exposure.zip) are only used for the exposure analysis. The contents of [ref_data.zip](https://gpm.nasa.gov/sites/default/files/data/landslides/ref_data.zip) are only used for the global post-fire debris flow analysis. 

### Installation

#### Linux

After cloning this repository, some setup is required prior to running LHASA. The following commands have been tested in a linux environment.

    cd LHASA
    # Set up python environment
    wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    sh Miniforge3-Linux-x86_64.sh
    conda env create -f lhasa.yml

    # Manage Earthdata connection
    # See for more info: https://disc.gsfc.nasa.gov/data-access
    touch ~/.urs_cookies
    touch ~/.netrc
    echo "machine urs.earthdata.nasa.gov login <uid> password <password>" >> ~/.netrc
    touch ~/.dodsrc
    # Manage PPS connection, which is only necessary for downloading IMERG HDF5
    # See for more info: https://registration.pps.eosdis.nasa.gov/registration/
    echo "machine jsimpsonhttps.pps.eosdis.nasa.gov login <email>  password <email>" >> ~/.netrc
    echo "HTTP.NETRC=~/.netrc" >> ~/.dodsrc
    echo "HTTP.COOKIEJAR=~/.urs_cookies" >> ~/.dodsrc
    
    # Set up directory structure
    mkdir -p nrt/hazard/tif
    mkdir -p nrt/exposure/csv
    mkdir -p fcast/hazard/tif
    mkdir -p fcast/exposure/csv
    mkdir imerg

    # Obtain required data files
    wget https://gpm.nasa.gov/sites/default/files/data/landslides/static.zip &&
    unzip static.zip &&
    rm static.zip

    wget https://gpm.nasa.gov/sites/default/files/data/landslides/exposure.zip &&
    unzip exposure.zip &&
    rm exposure.zip

    wget https://gpm.nasa.gov/sites/default/files/data/landslides/ref_data.zip &&
    unzip ref_data.zip -d pfdf/ &&
    rm ref_data.zip

    # Configure authorization for post-fire debris flow model
    python pfdf/scripts/make_netrc.py

#### Windows

On a Windows machine, use following commands within a *command prompt*:

```
@REM Set up python environment
curl -O https://github.com/conda-forge/miniforge/releases/download/25.3.0-1/Miniforge3-25.3.0-1-Windows-x86_64.exe
start /wait Miniforge3-25.3.0-1-Windows-x86_64.exe /InstallationType=JustMe /RegisterPython=0 /S /D=%UserProfile%\Miniconda3
%UserProfile%\Miniconda3\Scripts\conda.exe env create -f lhasa.yml

@REM Manage Earthdata connection
@REM See for more info: https://disc.gsfc.nasa.gov/data-access
type nul > %UserProfile%\.urs_cookies
type nul > %UserProfile%\.netrc
echo machine urs.earthdata.nasa.gov login <uid> password <password> >> %UserProfile%\.netrc
type nul > %UserProfile%\.dodsrc

@REM Manage PPS connection, which is only necessary for downloading IMERG HDF5
@REM See for more info: https://registration.pps.eosdis.nasa.gov/registration/
echo machine jsimpsonhttps.pps.eosdis.nasa.gov login <email> password <email> >> %UserProfile%\.netrc
echo HTTP.NETRC=%UserProfile%\.netrc >> %UserProfile%\.dodsrc
echo HTTP.COOKIEJAR=%UserProfile%\.urs_cookies >> %UserProfile%\.dodsrc

@REM Set up directory structure
mkdir nrt\hazard\tif
mkdir nrt\exposure\csv
mkdir fcast\hazard\tif
mkdir fcast\exposure\csv
mkdir imerg

@REM Obtain required data files
curl -O https://gpm.nasa.gov/sites/default/files/data/landslides/static.zip && ^
tar -xf static.zip && ^
del static.zip

curl -O https://gpm.nasa.gov/sites/default/files/data/landslides/exposure.zip && ^
tar -xf exposure.zip && ^
del exposure.zip

curl -O https://gpm.nasa.gov/sites/default/files/data/landslides/ref_data.zip && ^
tar -xf ref_data.zip -C pfdf/ && ^
del ref_data.zip

@REM Configure authorization for post-fire debris flow model
python pfdf\scripts\make_netrc.py
```

### Routine operation

Run [lhasa.sh](https://github.com/nasa/LHASA/blob/master/lhasa.sh) at the desired cadence, e.g. once per day. 

### Citation

Stanley T. A., J. R. Sutton, R. S. Vershel, P. M. Amatya. 2025. "Better Satellite Precipitation Algorithms Slightly Improved Landslide Hazard Assessment." Journal of Applied Meteorology and Climatology 64 (10): 1379-1394 [10.1175/jamc-d-25-0021.1](https://doi.org/10.1175/JAMC-D-25-0021.1)

Orland, E., D. Kirschbaum, and T. Stanley. 2022. "A Scalable Framework for Post Fire Debris Flow Hazard Assessment Using Satellite Precipitation Data." Geophysical Research Letters, 49 (18): [10.1029/2022gl099850](https://doi.org/10.1029/2022GL099850)

Khan, S., D. B. Kirschbaum, T. A. Stanley, P. M. Amatya, and R. A. Emberson. 2022. "Global Landslide Forecasting System for Hazard Assessment and Situational Awareness." Frontiers in Earth Science, 10: [10.3389/feart.2022.878996](https://doi.org/10.3389/feart.2022.878)

Emberson, R., D. Kirschbaum, and T. Stanley. 2020. "New global characterisation of landslide exposure." Natural Hazards and Earth System Sciences, 20 (12): 3413-3424 [10.5194/nhess-20-3413-2020](https://doi.org/10.5194/nhess-20-3413-2020)

### Model training

The software released here enables the user to run the global landslide forecast, but it does not enable the user to retrain the model on new datasets or domains. However, a demonstration workflow similar to that used in global LHASA 2.0 can be viewed [here](https://git.smce.nasa.gov/eis-freshwater/landslides/-/blob/master/brendan/Landslide-Case-Study.ipynb). This demo was created as part of the [EIS](https://freshwater.eis.smce.nasa.gov/storymap.html?story=ls) project funded by NASA. 

### Archive

A long-term archive for hazard maps from LHASA 2.0 is available at [GES-DISC](https://disc.gsfc.nasa.gov/datasets/Global_Landslide_Nowcast_2.0.0/summary). 

### Contributing

Users are encouraged to participate in this project in various ways. 

New landslide reports can be made through [Landslide Reporter](https://landslides.nasa.gov/reporter), which will enable NASA to better validate the model in the future. 

Bug reports can be made through GitHub issues, while bug fixes and feature updates are welcome through pull requests. However, it's best to contact NASA prior to embarking on a major feature, as some improvements may lie outside the scope of this project. 

Various forms of documentation are also needed. For example, a guide to installation of LHASA on Windows has already been requested. 

---

## LHASA 1.1

Although version 2 surpasses version 1 in accuracy and features, some users may prefer the simplicity of a single heuristic decision tree. Therefore, legacy code for LHASA version 1.1.1 is available [here](https://github.com/nasa/LHASA/releases/tag/v1.1.1). The R scripts are written to be easily understood, executed, and modified by potential users of this research.

Code for LHASA 1.0 is available in python at https://github.com/vightel/ojo-processing. 

### Data files

LHASA 1.1 requires the use of 2 data files, the 95th percentile rainfall and the global landslide susceptibility map. While the former is bundled with the [code release](https://github.com/nasa/LHASA/releases/tag/v1.1.1), the latter is too large and must be downloaded from https://gpm.nasa.gov/sites/default/files/downloads/global-landslide-susceptibility-map-2-27-23.tif. Users are encouraged to map susceptibility with current information on their own study areas, as well as update the rainfall threshold as needed. 

### Citation

Emberson, R., D. Kirschbaum, and T. Stanley. 2020. "New global characterisation of landslide exposure." Natural Hazards and Earth System Sciences, 20 (12): 3413-3424 [10.5194/nhess-20-3413-2020](https://doi.org/10.5194/nhess-20-3413-2020)

Kirschbaum, D., and T. Stanley. 2018. "Satellite-Based Assessment of Rainfall-Triggered Landslide Hazard for Situational Awareness." Earth's Future, 6 (3): 505-523 [10.1002/2017ef000715](https://doi.org/10.1002/2017ef000715)

Stanley, T., and D. B. Kirschbaum. 2017. "A heuristic approach to global landslide susceptibility mapping." Natural Hazards, 1-20 [10.1007/s11069-017-2757-y](https://doi.org/10.1007/s11069-017-2757-y)

### Archive

A long-term archive of hazard maps from LHASA 1.1 is available at [GES-DISC](https://disc.gsfc.nasa.gov/datasets/Global_Landslide_Nowcast_1.1/summary). An archive of exposure maps from LHASA 1.1 is available at [GES-DISC](https://disc.gsfc.nasa.gov/datasets/Global_Landslide_Exposure_Maps_1.0/summary). 