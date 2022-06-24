# coding=utf-8
# Copyright © 2020, United States Government, as represented by the
# Administrator of the National Aeronautics and Space Administration.
# All rights reserved.
#
# The LHASA (Landslide Hazard Analysis for Situational Awareness) system is
# licensed under NASA OPEN SOURCE AGREEMENT VERSION 1.3;
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#        https://github.com/nasa/LHASA/blob/master/LICENSE.pdf.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Functions and script for running LHASA 
"""

import logging
import os.path
import glob
import argparse
import requests
import numpy as np
import xarray as xr
import xgboost as xgb
import pandas as pd
import affine
import rasterio
from urllib.request import urlopen
import xml.etree.ElementTree as ET
import warnings

NO_DATA = -9999.0
OPENDAP_URL = 'https://gpm1.gesdisc.eosdis.nasa.gov/opendap'
PPS_URL = 'https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/'

def build_imerg_url(start_time, run='E', version='06C', opendap=True):
    """Build URL to IMERG data"""
    if opendap:
        product = f'GPM_3IMERGHH{run}.{version[:2]}'
        url = (f'{OPENDAP_URL}/ncml/aggregation/{product}/{start_time.year}/'
            f'{product}_Aggregation_{start_time.year}{start_time.dayofyear:03}'
            '.ncml')
    else:
        end = start_time + pd.Timedelta(minutes=29, seconds=59)
        product = f'GPM_3IMERGHH{run}.{version}'
        url = (f'{PPS_URL}{"early" if run=="E" else "late"}'
            f'/{start_time.year}{start_time.month:02}/'
            f'3B-HHR-{run}.MS.MRG.3IMERG.'
            f'{start_time.year}{start_time.month:02}{start_time.day:02}-S'
            f'{start_time.hour:02}{start_time.minute:02}{start_time.second:02}'
            f'-E{end.hour:02}{end.minute:02}{end.second:02}.'
            f'{start_time.hour*60+start_time.minute:04}.V{version}.RT-H5'
            )
    return url

def download_imerg(url, path='imerg'):
    """Downloads and saves IMERG file"""
    file_name = url[url.rfind('/')+1:len(url)]
    file_path = os.path.join(path, file_name)
    if not os.path.exists(file_path):
        request = requests.get(url)
        if request.ok:
            os.makedirs(path, exist_ok=True)
            with open(file_path, 'wb') as f:
                f.write(request.content)
        else: 
            raise RuntimeError(str(request.status_code) + ': could not download ' + request.url)
    return file_path

def get_latest_imerg_year(run='E', version='06'):
    """Finds the last year IMERG data is available at GES-DISC OpenDAP"""
    url = f'{OPENDAP_URL}/ncml/aggregation/GPM_3IMERGHH{run}.{version}/catalog.xml'
    catalog = ET.fromstring(requests.get(url).content)
    name_space = {'thredds': 'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0'}
    for dataset in reversed(catalog.findall(
            'thredds:dataset/thredds:catalogRef', 
            name_space)):
        year = dataset.attrib['name']
        if year.isnumeric():
            return int(year)
    raise RuntimeError(f'cannot parse imerg catalog at {url}')

def get_latest_imerg_url(run='E', version='06'):
    """Returns a path to the latest 30-minute IMERG data at GES-DISC OpenDAP"""
    version = str(version)[:2].zfill(2)
    year = get_latest_imerg_year(run=run, version=version)
    url = (f'{OPENDAP_URL}/ncml/aggregation/GPM_3IMERGHH{run}.{version}/{year}'
        '/catalog.xml')
    catalog = ET.fromstring(requests.get(url).content)
    name_space = {
        'thredds': 
        'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0'
    }
    datasets = reversed(catalog.findall(
        'thredds:dataset/thredds:dataset', 
        name_space
    ))
    for dataset in datasets: 
        day_url = f"{OPENDAP_URL}{dataset[2].attrib['urlPath']}"
        try: # Some catalog entries are empty
            xr.open_dataset(day_url)
        except (KeyError, OSError):
            continue
        return day_url
    raise RuntimeError(f'Cannot parse imerg catalog at {url}')

def get_latest_imerg_time(run='E', version='06C', opendap=True):
    """Returns a pandas time stamp representing the latest available data"""
    if opendap:
        url = get_latest_imerg_url(run=run, version=version)
        imerg = xr.open_dataset(f'{url}')
        latest = imerg['time'].max()
        t = pd.Timestamp(
            year=int(latest.dt.year),
            month=int(latest.dt.month),
            day=int(latest.dt.day),
            hour=int(latest.dt.hour),
            minute=int(latest.dt.minute)
        )
        return t
    else:
        now = pd.Timestamp.now(tz='UTC').floor('30min')
        for i in range(6, 48):
            t = now - pd.Timedelta(hours=i/2)
            url = build_imerg_url(t, run=run, version=version, opendap=False)
            if requests.get(url).ok: 
                return t.tz_localize(None)
        raise RuntimeError(f'No IMERG data available at {url}')

def get_valid_SMAP_time(start_time):
    day = start_time.strftime('%Y-%m-%d')
    hours = pd.date_range(day +" 01:30", periods=9, freq="3h")
    return hours[hours.get_loc(start_time, method = 'nearest')]

def build_SMAP_url(start_time, version='6', minor_version='030'):
    """Build string for SMAP OpenDAP server at NSIDC
    For more information, see https://nsidc.org/data/smap/data_versions#L4"""
    return (f'https://n5eil02u.ecs.nsidc.org/opendap/SMAP/SPL4SMGP.'
        f'{int(version):03}/{start_time.year}.{start_time.month:02}.'
        f'{start_time.day:02}/SMAP_L4_SM_gph_'
        f'{start_time.year}{start_time.month:02}{start_time.day:02}'
        f'T{start_time.hour:02}{start_time.minute:02}{start_time.second:02}'
        f'_Vv{version}{minor_version}_001.h5')

def get_SMAP(start_time, variables=['Geophysical_Data_sm_profile_wetness', 'Geophysical_Data_snow_mass'], 
        load=True, version='6', minor_version='030', latitudes=slice(-60, 60), longitudes=slice(-180, 180)):
    """returns an xarray dataset representing 1 SMAP file"""
    start_time = get_valid_SMAP_time(start_time)
    url = build_SMAP_url(start_time, version, minor_version)
    try:
        smap = xr.open_dataset(url)
    except:
        logging.warning('SMAP data is not available for'
            ' the specified time and version.')
        for i in range(1, 9):
            start_time = start_time - pd.Timedelta(hours=i*3)
            url = build_SMAP_url(start_time, version, minor_version)
            #TODO check for availability without triggering warnings
            try: 
                smap = xr.open_dataset(url)
            except (KeyError, OSError):
                continue
            break
        logging.warning(f'Using latest available SMAP data: {start_time}')
    # Convert to WGS84
    smap['x'] = smap.cell_lon[1,]
    smap['y'] = smap.cell_lat[:,1]
    smap = smap.sortby('y')
    smap = smap.sel(y=latitudes, x=longitudes)
    smap = smap[variables]
    if load: 
        smap.load()
    return smap.rename(y='lat', x='lon')

def get_IMERG_precipitation(start_time: pd.Timestamp, end_time: pd.Timestamp, 
        liquid=True, load=True, run='E', version='06C', opendap=True,
        latitudes=slice(-60, 60), longitudes=slice(-180, 180)):
    """Opens IMERG data"""
    if end_time <= start_time:
        raise ValueError('End time must be later than start time')
    if opendap:
        days = pd.date_range(start_time, end_time, freq='D')
        files = [build_imerg_url(d, run=run, version=version) for d in days]
        imerg = xr.open_mfdataset(files, parallel=True)
    else:
        half_hours = pd.date_range(start_time, end_time - pd.Timedelta(minutes=30), freq='30min')
        urls = [build_imerg_url(h, run=run, version=version, opendap=False) for h in half_hours]
        files = [download_imerg(u) for u in urls]
        imerg = xr.open_mfdataset(files, group='Grid', parallel=True)
    with warnings.catch_warnings():
        warnings.filterwarnings(action='ignore', category=RuntimeWarning)
        warnings.filterwarnings(action='ignore', category=FutureWarning)
        imerg['time'] = imerg.indexes['time'].to_datetimeindex()
    imerg = imerg.sel(time=slice(start_time, end_time))
    imerg = imerg.sel(lat=latitudes, lon=longitudes)
    if liquid: 
        precipitation = imerg['precipitationCal'] * (imerg['probabilityLiquidPrecipitation'] > 0.5)
    else: 
        precipitation = imerg['precipitationCal']
    if load: 
        precipitation.load()
    return precipitation

def build_GEOS_url(run_time=None, mode='fcast'):
    """Build string for the GEOS OpenDAP server at the NCCS Data Portal"""
    if mode == 'fcast':
        if run_time == None:
            raise ValueError('GEOS FP model run time must be specified for forecast mode')
        return f'https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/fcast/'\
            f'tavg1_2d_lnd_Nx/tavg1_2d_lnd_Nx.{run_time.year}{run_time.month:02}{run_time.day:02}'\
            f'_{run_time.round("6H").hour:02}'
    if mode == 'assim':
        return 'https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_lnd_Nx'

def get_latest_GEOS_run_time(end_time=None):
    """Returns a pandas time stamp representing the latest available data file"""
    if end_time: 
        now = end_time
    else:
        now = pd.Timestamp.now(tz='UTC')
    for i in range(0, 5):
        latest = now.floor('6H') - pd.Timedelta(hours=i*6)
        url = build_GEOS_url(latest)
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings(action='ignore')
                geos = xr.open_dataset(url)
            if geos: 
                break
        except:
            latest = None
            continue
    return latest

def get_GEOS_variable(start_time: pd.Timestamp, run_time=None, 
        variable='gwetprof', mode='assim', 
        latitudes=slice(-60, 60), longitudes=slice(-180, 180)):
    """returns an xarray dataarray for a specified variable"""
    url = build_GEOS_url(run_time, mode=mode)
    with warnings.catch_warnings():
        warnings.filterwarnings(action='ignore')
        geos = xr.open_dataset(url)
    geos = geos.sel(time=start_time, method='nearest')
    geos = geos.sel(lon=longitudes, lat=latitudes)
    return geos[variable].drop('time')

def get_GEOS_run(run_time: pd.Timestamp, mode='fcast'):
    url = build_GEOS_url(run_time=run_time, mode=mode)
    return xr.open_dataset(url)

def get_GEOS_precipitation(geos: xr.Dataset, liquid=True, load=True, 
        latitudes=slice(-60, 60), longitudes=slice(-180, 180)) -> xr.DataArray:
    """returns a xarray dataarray"""
    geos = geos.sel(lon=longitudes, lat=latitudes)
    if liquid: 
        precipitation = geos['prectot'] - geos['precsno']
    else: 
        precipitation = geos['prectot']
    if load: 
        precipitation.load()
    return precipitation * 3600

def regrid(variable, template, method='nearest'):
    """Change resolution of dataarray"""
    return variable.interp_like(template, method=method)

def apply_mask(variable, mask):
    """Mask dataarray and return numpy array"""
    return variable.values[mask.values]

def fill_array(prediction, mask, start_time):
    """Convert numpy to xarray format with time and location"""
    p_values = mask.where(mask).values
    p_values[mask.values] = prediction
    filled_array = xr.DataArray(
        p_values.reshape((1,) + mask.shape),
        coords=[[start_time], mask['lat'], mask['lon']], 
        dims=['time', 'lat', 'lon'],
        name='p_landslide')
    return filled_array

def add_metadata(data_set: xr.Dataset,  run_mode='nrt'):
    """Adds metadata for compliance with CF and GES-DISC standards"""
    data_set.attrs['title'] = 'Landslide Hazard Analysis for Situational Awareness'
    data_set.attrs['institution'] = 'NASA GSFC'
    data_set.attrs['source'] = 'LHASA V2.0.0a'
    data_set.attrs['history'] = f'{pd.Timestamp.now()} File written by XArray version {xr.__version__}'
    if run_mode == 'nrt': 
        data_set.attrs['references'] = (
            'Stanley, T. A., D. B. Kirschbaum, G. Benz, et al. 2021. '
            '"Data-Driven Landslide Nowcasting at the Global Scale." '
            'Frontiers in Earth Science, 9: [10.3389/feart.2021.640043]'
        )
    else: 
        data_set.attrs['references'] = (
            'Khan, S., D. B. Kirschbaum, T. A. Stanley, P. M. Amatya and R. '
            'Emberson. 2022. "Global Landslide Forecasting System for Hazard '
            'Assessment and Situational Awareness" Frontiers in Earth Science'
        )
    data_set.attrs['comment'] = 'LHASA identifies where landslides are most probable in near real time.'
    data_set.attrs['Conventions'] = 'CF-1.8'
    data_set.attrs['ShortName'] = 'LHASA'
    data_set.attrs['LongName'] = 'Landslide Hazard Analysis for Situational Awareness (LHASA)'
    data_set.attrs['VersionID'] = '2.0.0a'
    data_set.attrs['Format'] = 'netCDF-4'
    data_set.attrs['DataSetQuality'] = 'NRT'
    data_set.attrs['IdentifierProductDOIAuthority'] = 'https://doi.org/'
    data_set.attrs['IdentifierProductDOI'] = '10.5067/nnnnn'
    data_set.attrs['ProcessingLevel'] = '4'
    data_set['lat'].attrs['long_name'] = 'Latitude (Degrees)'
    data_set['lat'].attrs['standard_name'] = 'latitude'
    data_set['lat'].attrs['units'] = 'degrees_north'
    data_set['lat'].attrs['valid_min'] = -90.0
    data_set['lat'].attrs['valid_max'] = 90.0
    data_set['lon'].attrs['long_name'] = 'Longitude (Degrees)'
    data_set['lon'].attrs['standard_name'] = 'longitude'
    data_set['lon'].attrs['units'] = 'degrees_east'
    data_set['lon'].attrs['valid_min'] = -180.0
    data_set['lon'].attrs['valid_max'] = 180.0
    data_set['time'].attrs['long_name'] = 'Time (UTC)'
    data_set['time'].attrs['standard_name'] = 'time'
    data_set['p_landslide'].attrs['long_name'] = 'Probability of Landslide Occurrence'

def save_nc(data_set: xr.Dataset, file_path: str, run_mode='nrt'):
    """Saves prediction in netcdf format"""
    add_metadata(data_set, run_mode=run_mode)
    data_set.to_netcdf(file_path, encoding = {
        'p_landslide': {'zlib': True, '_FillValue': NO_DATA},
        'lat': {'zlib': False, '_FillValue': None},
        'lon': {'zlib': False, '_FillValue': None}
    })

def save_tiff(data_array, file_path):
    """Saves prediction in geotiff format"""
    
    cell_size = 0.00833333333333333
    metadata = {
        'driver': 'GTiff', 
        'compress': 'lzw',
        'dtype': data_array.dtype, 
        'nodata': -9999.0, 
        'width': data_array['lon'].size, 
        'height': data_array['lat'].size, 
        'count': 1, 
        'crs': 'EPSG:4326', 
        'transform': affine.Affine(
            cell_size, 
            0.0, 
            data_array['lon'].min() - cell_size/2, 
            0.0, 
            -cell_size, 
            data_array['lat'].max() + cell_size/2), 
        'tiled': False, 
        'interleave': 'band'
    }
    with rasterio.open(file_path, 'w', **metadata) as dst:
        dst.write(data_array.values)

def get_model(file_path, threads=1):
    """Open trained model"""
    model = xgb.Booster()
    model.load_model(file_path)
    model.set_param('nthread', threads)
    return model

def imerg_cleanup(path: str, cache_days=0, cache_end_time=None):
    files = glob.glob(os.path.join(path, 'imerg', '*'))
    if cache_days < 1:
        for f in files:
            os.remove(f)
    else:
        cache_start_time = cache_end_time - pd.DateOffset(days=cache_days)
        df = pd.DataFrame({'file': files})
        df['time'] = df['file'].apply(lambda x: 
            pd.Timestamp(x[67:81].replace('-S', ' ')))
        excess = df['file'][~df['time'].between(cache_start_time, cache_end_time)]
        excess.map(os.remove)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='LHASA 2.0 global landslide forecast')
    parser.add_argument('-v', '--version', action='version', 
        version='LHASA version 2.0.0a')
    parser.add_argument('-p', '--path', default=os.getcwd(), 
        help='location of input files')
    parser.add_argument('-op', '--output_path', help='location of output files')
    parser.add_argument('-d', '--date', 
        help='UTC date and time, formatted as "YYYY-MM-DD HH:MM", assumes '
        'latest available if unspecified')
    parser.add_argument('-l', '--lead', type=int,  default=2, 
        help='Days of forecast to generate. 0 means only LHASA-NRT will run.')
    parser.add_argument('-rt', '--run_time', 
        help='GEOS-FP model initiation time (UTC), formatted as '
        '"YYYY-MM-DD HH", assumes latest available if unspecified')
    parser.add_argument('-o', '--overwrite', action="store_true", 
        help='assumed false if this option is not provided')
    parser.add_argument('--small', action="store_true", 
        help='shrinks output file size by masking where p_landslide < 0.01')
    parser.add_argument('-t', '--threads', type=int,  default=4, 
        help='Number of threads to use. Currently, this only affects XGBoost.')
    parser.add_argument('-f', '--format', default='nc4', 
        choices=['nc4', 'tif', 'nc4tif', 'nc', 'nctif'], 
        help='Output file formats. Both geotiff and netCDF can be chosen.'
    )
    parser.add_argument('-ex', '--exposure', action="store_true", 
        help='The memory-intensive exposure analysis is omitted unless this'
        ' argument is used')
    parser.add_argument('-N', '--north', type=float,  default=60.0, 
        help='maximum latitude (WGS84)')
    parser.add_argument('-S', '--south', type=float,  default=-60.0, 
        help='minimum latitude (WGS84)')
    parser.add_argument('-E', '--east', type=float,  default=180.0, 
        help='maximum longitude (WGS84)')
    parser.add_argument('-W', '--west', type=float,  default=-180.0, 
        help='minimum longitude (WGS84)')
    parser.add_argument('-sv', '--smap_version', default='6030',
        help='SMAP L4 major and minor version, e.g. 6030')
    parser.add_argument('-iv', '--imerg_version', default='06C',
        help='IMERG version, e.g. 06C')
    parser.add_argument('-icd', '--imerg_cache_days', type=int,  default=0, 
        help='Days of IMERG data to cache.')
    parser.add_argument('-od', '--opendap', action='store_true',
        help='Use OpenDAP to access IMERG; LHASA defaults to HDF5 download')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s', 
        datefmt='%H:%M')
    logging.info('logging started')

    if args.west > args.east:
        raise ValueError(f'West longitude ({args.west}) cannot be greater than East longitude ({args.east})')
    if args.south > args.north:
        raise ValueError(f'South latitude ({args.south}) cannot be greater than North latitude ({args.north})')

    path = args.path
    if args.output_path:
        output_path = args.output_path
    else: 
        output_path = path

    if args.date:
        forecast_start_time = pd.Timestamp(args.date)
    else:
        forecast_start_time = get_latest_imerg_time(run='E') + pd.Timedelta(minutes=30)
    if args.lead > 0:
        # GEOS is hourly data, so we may have to discard the last half hour of IMERG
        forecast_start_time = forecast_start_time.floor('H')
        if args.run_time:
            run_time = pd.Timestamp(args.run_time).floor('6H')
            if run_time > forecast_start_time:
                raise ValueError("Forecast must start later than GEOS-FP run time")
        else:
            run_time = get_latest_GEOS_run_time(forecast_start_time)
        if not run_time:
            raise RuntimeError(f'No GEOS-FP run is available for {forecast_start_time}.')
    if args.lead > 2:
        warnings.warn('Only 2 days of forecast data are considered reliable.')

    model = get_model(f'{path}/model.json', args.threads)

    static_files = [f'{path}/static/Faults.nc4',
        f'{path}/static/Lithology.nc4',
        f'{path}/static/Slope.nc4',
        f'{path}/static/mask.nc4']
    static_variables = xr.open_mfdataset(static_files, parallel=True)
    static_variables = static_variables.sel(
        lat=slice(args.south, args.north), 
        lon=slice(args.west, args.east)
    )
    static_variables.load()
    p99 = xr.open_dataarray(f'{path}/static/p99.nc4')
    p99 = p99.sel(
        lat=slice(args.south, args.north), 
        lon=slice(args.west, args.east)
    ).load()
    if args.lead > 0:
        #TODO: Replace with file at native resolution
        p99geos = xr.open_dataarray(f'{path}/static/p99GEOS_latest.nc4')
        p99geos = p99geos.sel(
            lat=slice(args.south, args.north), 
            lon=slice(args.west, args.east)
        ).load()
    logging.info('opened static variables')

    imerg_late_end_time = get_latest_imerg_time(run='L') + pd.Timedelta(minutes=30)
    precipitation_start_date = forecast_start_time - pd.DateOffset(3)
    precipitation_end_date = forecast_start_time + pd.DateOffset(args.lead)

    imerg_late = get_IMERG_precipitation(
        precipitation_start_date, 
        min(precipitation_end_date, imerg_late_end_time), 
        liquid=True, 
        load=False, 
        run='L', 
        version=args.imerg_version, 
        opendap=args.opendap,
        latitudes=slice(args.south, args.north), 
        longitudes=slice(args.west, args.east)
    )

    if imerg_late_end_time < forecast_start_time:
        imerg_early = get_IMERG_precipitation(
            max(precipitation_start_date, imerg_late_end_time), 
            forecast_start_time, 
            liquid=True, 
            load=False, 
            run='E', 
            version=args.imerg_version, 
            opendap=args.opendap,
            latitudes=slice(args.south, args.north), 
            longitudes=slice(args.west, args.east)
        )
        imerg = xr.concat([imerg_late, imerg_early], dim='time')
    else:
        imerg = imerg_late

    imerg_daily = imerg.resample({'time': '24H'}, 
        base=forecast_start_time.hour).sum().load()

    if args.lead > 0:
        geos_run_times = [run_time, run_time - pd.Timedelta(hours=6)]
        geos_runs = [get_GEOS_run(t) for t in geos_run_times]

        processed_hours = geos_runs[0]['time'].where(False) # empty array
        selected_geos_runs = []
        for geos_run in reversed(geos_runs):
            geos_run = geos_run.sel(
                time=slice(
                    forecast_start_time, 
                    precipitation_end_date
                )
            )
            new_hours = geos_run['time'].loc[~geos_run['time'].isin(processed_hours)]
            if new_hours.size > 0:
                selected_geos_runs.append(geos_run.sel(time=new_hours))
                processed_hours = processed_hours.combine_first(new_hours).dropna('time')

        selected_geos_precip = [
            get_GEOS_precipitation(
                r, 
                liquid=True,
                load=False,
                latitudes=slice(args.south, args.north), 
                longitudes=slice(args.west, args.east)
            ) for r in selected_geos_runs
        ]

        geos = xr.concat(selected_geos_precip, dim='time').sortby('time')
        geos_daily = geos.resample({'time': '24H'}, 
            base=forecast_start_time.hour).sum().load()
        daily_rain = [da for da in imerg_daily] + [da for da in geos_daily]
    else: 
        daily_rain = [da for da in imerg_daily]
        
    dates = pd.date_range(
        forecast_start_time - pd.DateOffset(1), 
        precipitation_end_date - pd.DateOffset(1)
    )
    for i, d in enumerate(daily_rain[2:]): 
        date_string = dates[i].strftime('%Y%m%dT%H%M')
        run_mode = "nrt" if i < 1 else "fcast"
        if run_mode == 'fcast':
            initialization = f"{run_time.strftime('%Y%m%dT%H%M')}+" 
        else: 
            initialization = ''
        if 'nc' in args.format:
            file_ext = args.format.strip('tif')
            nc_path = (f'{output_path}/{run_mode}/hazard/'
                f'{initialization}{date_string}.{file_ext}')
            if os.path.exists(nc_path) and not args.overwrite:
                raise FileExistsError(f'{nc_path} exists. Use the --overwrite'
                    ' argument to overwrite it.')
        else: 
            nc_path = None
        if 'tif' in args.format:
            tif_path = (f'{output_path}/{run_mode}/hazard/tif/'
                f'{initialization}{date_string}.tif')
            if os.path.exists(tif_path) and not args.overwrite:
                raise FileExistsError(f'{tif_path} exists. Use the --overwrite'
                    ' argument to overwrite it.')
        else: 
            tif_path = None

        day_before = daily_rain[i]
        yesterday = daily_rain[i + 1].interp_like(day_before, method='nearest')
        antecedent = xr.concat([yesterday, day_before], dim='time').sum('time')
        antecedent.name = 'antecedent'

        rain = d
        if run_mode == "nrt": 
            rain = rain.interp_like(p99, method='nearest') / p99
            smap = get_SMAP(
                start_time=dates[i] - pd.Timedelta(days=2, hours=1), 
                version=args.smap_version[0], 
                minor_version=args.smap_version[1:4], 
                latitudes=slice(args.south, args.north), 
                longitudes=slice(args.west, args.east)
            )
            moisture = smap['Geophysical_Data_sm_profile_wetness']
            snow = smap['Geophysical_Data_snow_mass']
            # Match variable names in GEOS
            moisture.name = 'gwetprof'
            snow.name = 'snomas'
        else: 
            moisture = get_GEOS_variable(
                start_time=dates[i] - pd.Timedelta(days=2, hours=1), 
                run_time=run_time, 
                variable='gwetprof', 
                mode='assim', 
                latitudes=slice(args.south, args.north), 
                longitudes=slice(args.west, args.east)
            )
            snow = get_GEOS_variable(
                start_time=dates[i] - pd.Timedelta(days=2, hours=1), 
                run_time=run_time, 
                variable='snomas', 
                mode='assim', 
                latitudes=slice(args.south, args.north), 
                longitudes=slice(args.west, args.east)
            )
            rain = d.interp_like(p99geos, method='nearest') / p99geos
        rain.name = 'rain'
        dynamic_variables = [rain, antecedent, moisture, snow]
        transposed = [v.transpose('lat', 'lon') for v in dynamic_variables]
        logging.info('opened ' + date_string)
        regridded = xr.merge([regrid(v, static_variables) for v in transposed])
        logging.info('interpolated')
        # It's important to retain the correct variable order
        variable_order = [
            'rain', 
            'antecedent', 
            'gwetprof', 
            'snomas', 
            'Faults', 
            'Lithology', 
            'Slope'
        ]
        variables = xr.merge([regridded, static_variables.drop('land_mask')])
        masked_values = [apply_mask(variables[v], mask=(static_variables['land_mask'] > 0)) for v in variable_order]
        logging.info('built mask')
        inputs = xgb.DMatrix(np.stack(masked_values, 1))
        prediction = model.predict(inputs)
        logging.info('made predictions')
        p_landslide = fill_array(prediction, static_variables['land_mask'] > 0, dates[i])
        # Reindex for faster loading in ArcGIS
        p_landslide = p_landslide.reindex(time=p_landslide.time, lat=p_landslide.lat[::-1])
        if args.small:
            p_landslide = p_landslide.where(p_landslide >= 0.01)
        if nc_path:
            save_nc(p_landslide.to_dataset(), nc_path)
        if tif_path:
            save_tiff(p_landslide, tif_path)
        logging.info('saved output to disk')

        if args.exposure:
            csv_path = f'{output_path}/{run_mode}/exposure/csv/{date_string}.csv'
            if os.path.exists(csv_path) and not args.overwrite:
                raise FileExistsError(f'{csv_path} exists. Use the --overwrite'
                    ' argument to overwrite it.')
            variable_files = [f'{path}/exposure/road_length.nc4', 
            f'{path}/exposure/population.nc4', f'{path}/exposure/gadm36.nc4']
            variables = xr.open_mfdataset(variable_files)
            #TODO: reduce code redundancy
            variables['l_haz'] = p_landslide > 0.1
            variables['m_haz'] = p_landslide > 0.5
            variables['h_haz'] = p_landslide > 0.9
            variables['l_haz_pp'] = variables['population'] * variables['l_haz']
            variables['m_haz_pp'] = variables['population'] * variables['m_haz']
            variables['h_haz_pp'] = variables['population'] * variables['h_haz']
            variables['l_haz_rd'] = variables['road_length'] * variables['l_haz']
            variables['m_haz_rd'] = variables['road_length'] * variables['m_haz']
            variables['h_haz_rd'] = variables['road_length'] * variables['h_haz']
            variables = variables.drop(['time', 'lat', 'lon', 'road_length', 'population']).squeeze()
            data_frame = variables.to_dataframe().dropna()
            logging.info('opened exposure variables')

            totals = data_frame.groupby('gadm_fid').sum()

            constant_totals = pd.read_csv(f'{path}/exposure/totals.csv', index_col='gadm_fid')
            # Apply a minimum to avoid divide by zero errors
            constant_totals['road_length'] = np.maximum(constant_totals['road_length'], 1e-10)
            constant_totals['population'] = np.maximum(constant_totals['population'], 1e-10)

            totals = totals.merge(constant_totals, on='gadm_fid')
            totals['l_haz_f'] = totals['l_haz'] / totals['cells']
            totals['m_haz_f'] = totals['m_haz'] / totals['cells']
            totals['h_haz_f'] = totals['h_haz'] / totals['cells']
            totals['l_haz_pp_f'] = totals['l_haz_pp'] / totals['population']
            totals['m_haz_pp_f'] = totals['m_haz_pp'] / totals['population']
            totals['h_haz_pp_f'] = totals['h_haz_pp'] / totals['population']
            totals['l_haz_rd_f'] = totals['l_haz_rd'] / totals['road_length']
            totals['m_haz_rd_f'] = totals['m_haz_rd'] / totals['road_length']
            totals['h_haz_rd_f'] = totals['h_haz_rd'] / totals['road_length']
            totals = totals.drop(columns=['population', 'road_length', 'cells'])

            admin_names = pd.read_csv(f'{path}/exposure/admin_names.csv', index_col='gadm_fid')
            admin_names = admin_names.merge(totals, on='gadm_fid', how='outer')
            admin_names.to_csv(csv_path)
            logging.info(f'saved {csv_path}')
    
    imerg_cleanup(
        path=path, 
        cache_end_time=forecast_start_time, 
        cache_days=args.imerg_cache_days
    )
