# Landslide Hazard Assessment for Situational Awareness

LHASA was created at Goddard Space Flight Center to identify the potential for rainfall-triggered landslides in real time. 

## What's new
2022-6-23 Added an option to choose between OpenDAP and file download servers for obtaining IMERG. OpenDAP enables reduced data volumes, both for transfer and storage; while file downloads from the PPS site may offer reduced latency and improved reliability. This capacity is controlled with the '-od' or '--opendap' argument in [lhasa.py](https://github.com/nasa/LHASA/blob/master/lhasa.py).

## LHASA 2.0

LHASA version 2 adopts machine learning to estimate the probability of landslide occurrence at a 30-arcsecond (~1 km) daily resolution. In addition, it estimates the potential exposure of human population and roads to landslide hazard and maps the basins likely to experience post-fire debris flows. 

### Real-time data availability
The latest predictions can be downloaded from https://maps.nccs.nasa.gov/download/landslides/nowcast. It can also be accessed as an ArcGIS web map at https://landslides.nasa.gov/viewer. NASA provides these data on a best-effort basis, typically four times each day, but with frequent server downtime. Users requiring a fully operational system are encouraged to clone this repository and run LHASA at the desired cadence.

### Data files

LHASA requires several large data files, but not all data may be needed by all users. The contents of [static.zip](https://gpm.nasa.gov/sites/default/files/data/landslides/static.zip) are required for the global landslide forecast. The contents of [exposure.zip]() are only used for the exposure analysis. The contents of [ref_data.zip](https://gpm.nasa.gov/sites/default/files/data/landslides/ref_data.zip) are only used for the global post-fire debris flow analysis. 

### Installation

After cloning this repository, some setup is required prior to running LHASA. The following commands have been tested in a linux environment. Users of Windows or other systems may be required to modify each of these steps. 

    # Set up python environment
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    sh Miniconda3-latest-Linux-x86_64.sh
    conda env create -f lhasa.yml

    # Manage Earthdata connection
    # See for more info: https://disc.gsfc.nasa.gov/data-access
    touch ~/.urs_cookies
    touch ~/.netrc
    echo "machine urs.earthdata.nasa.gov login <uid> password <password>" >> .netrc
    touch ~/.dodsrc
    # Manage PPS connection, which is only necessary for downloading IMERG HDF5
    # See for more info: https://registration.pps.eosdis.nasa.gov/registration/
    echo "machine jsimpsonhttps.pps.eosdis.nasa.gov login <email>  password <email>" >> .netrc
    echo "HTTP.NETRC=~/.netrc" >> .dodsrc
    echo "HTTP.COOKIEJAR=~/.urs_cookies" >> .dodsrc
    
    # Set up directory structure
    mkdir -p nrt/hazard/tif
    mkdir -p nrt/exposure/csv
    mkdir -p fcast/hazard/tif
    mkdir -p fcast/exposure/csv
    mkdir -p pfdf/firms

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

    # Configure post-fire debris flow model
    python pfdf/setup.py

The post-fire debris flow module uses Google Earth Engine to access Landsat imagery. Please see the [README](https://github.com/nasa/LHASA/blob/master/pfdf/README.md) for more information. 

### Routine operation

Once a month, run the following commands to build the fires database needed to run the post-fire debris flow module:

    conda activate lhasa
    python /scripts/gee_export_all.py --filepath /pfdf --gee_username username

Then run [lhasa.sh](https://github.com/nasa/LHASA/blob/master/lhasa.sh) at the desired cadence, e.g. once per day. 

### Citation

Khan, S., D. Kirschbaum, and T. Stanley. 2021. "Investigating the potential of a global precipitation forecast to inform landslide prediction." Weather and Climate Extremes, 33: 100364 [10.1016/j.wace.2021.100364](https://doi.org/10.1016/j.wace.2021.100364)

Khan, S., D. B. Kirschbaum, T. A. Stanley, P. M. Amatya and R. Emberson. 2022. "Global Landslide Forecasting System for Hazard Assessment and Situational Awareness" Frontiers in Earth Science

Emberson, R., D. Kirschbaum, and T. Stanley. 2020. "New global characterisation of landslide exposure." Natural Hazards and Earth System Sciences, 20 (12): 3413-3424 [10.5194/nhess-20-3413-2020](https://doi.org/10.5194/nhess-20-3413-2020)

Stanley, T. A., D. B. Kirschbaum, G. Benz, et al. 2021. "Data-Driven Landslide Nowcasting at the Global Scale." Frontiers in Earth Science, 9: [10.3389/feart.2021.640043](https://doi.org/10.3389/feart.2021.640043)

### Model training

The software released here enables the user to run the global landslide forecast, but it does not enable the user to retrain the model on new datasets or domains. However, a demonstration workflow similar to that used in global LHASA 2.0 can be viewed [here](https://git.mysmce.com/eis-freshwater/landslides/-/blob/master/brendan/Landslide-Case-Study.ipynb). This demo was created as part of the [EIS](https://eis.mysmce.com/) project funded by NASA. 

### Archive

No long-term archive for predictions from LHASA 2.0 has been established. Users of Google Earth Engine can find an interim archive at *users/nbiswasuw/Hazard*. 

### Contributing

Users are encouraged to participate in this project in various ways. 

New landslide reports can be made through [Landslide Reporter](https://landslides.nasa.gov/reporter), which will enable NASA to better validate the model in the future. 

Bug reports can be made through GitHub issues, while bug fixes and feature updates are welcome through pull requests. However, it's best to contact NASA prior to embarking on a major feature, as some improvements may lie outside the scope of this project. 

Various forms of documentation are also needed. For example, a guide to installation of LHASA on Windows has already been requested. 

---

## LHASA 1.1

Although version 2 surpasses version 1 in accuracy and features, some users may prefer the simplicity of a single heuristic decision tree. Therefore, LHASA 1.1 is still running and its output can be seen at https://pmm.nasa.gov/precip-apps.

Legacy code for LHASA version 1.1.1 is available [here](https://github.com/nasa/LHASA/releases/tag/v1.1.1). The R scripts are written to be easily understood, executed, and modified by potential users of this research.

Full operational code for LHASA 1.0 is available in python at https://github.com/vightel/ojo-bot. 

### Data files

LHASA 1.1 requires the use of 2 data files, the 95th percentile rainfall and the global landslide susceptibility map. While the former is bundled with the [code release](https://github.com/nasa/LHASA/releases/tag/v1.1.1), the latter is too large and must be downloaded from https://gpm.nasa.gov/sites/default/files/downloads/global-landslide-susceptibility-map-1-30-20.zip. Note that this dataset has not been updated since 2020, at which time the available data on forest loss were for 2018. Users are encouraged to map susceptibility with current information on their own study areas, as well as update the rainfall threshold as needed. 

### Citation

Emberson, R., D. Kirschbaum, and T. Stanley. 2020. "New global characterisation of landslide exposure." Natural Hazards and Earth System Sciences, 20 (12): 3413-3424 [10.5194/nhess-20-3413-2020](https://doi.org/10.5194/nhess-20-3413-2020)

Kirschbaum, D., and T. Stanley. 2018. "Satellite-Based Assessment of Rainfall-Triggered Landslide Hazard for Situational Awareness." Earth's Future, 6 (3): 505-523 [10.1002/2017ef000715](https://doi.org/10.1002/2017ef000715)

Stanley, T., and D. B. Kirschbaum. 2017. "A heuristic approach to global landslide susceptibility mapping." Natural Hazards, 1-20 [10.1007/s11069-017-2757-y](https://doi.org/10.1007/s11069-017-2757-y)

### Archive

A long-term archive for predictions from LHASA 1.1 is available at [GES-DISC](https://disc.gsfc.nasa.gov/datasets/Global_Landslide_Nowcast_1.1/summary).