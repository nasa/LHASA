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

import argparse
import glob
import io
import logging
import os.path
import warnings
import zipfile

import numpy as np
import pandas as pd
import requests
import rioxarray
import xarray as xr
import xgboost as xgb

NO_DATA = -9999.0
PPS_URL = 'https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/'

def build_imerg_url(start_time, run="E", version="07B", liquid=True):
    """Build URL to IMERG data"""
    start_time = start_time.ceil("3h") - pd.Timedelta(minutes=30)
    end = start_time + pd.Timedelta(minutes=29, seconds=59)
    if liquid:
        extension = ".zip"
    else:
        extension = ".tif"
    url = (
        f"{PPS_URL}gis/{'early/' if run=='E' else ''}"
        f"{start_time.year}/{start_time.month:02}/"
        f"3B-HHR-{run}.MS.MRG.3IMERG."
        f"{start_time.year}{start_time.month:02}{start_time.day:02}"
        f"-S{start_time.hour:02}{start_time.minute:02}00"
        f"-E{end.hour:02}{end.minute:02}59."
        f"{start_time.hour*60+start_time.minute:04}"
        f".V{version}.1day{extension}"
    )
    return url


def download_imerg(url, path="./imerg"):
    """Downloads and saves IMERG file"""
    file_name = os.path.basename(url)
    extension = os.path.splitext(url)[1]
    if extension == ".tif":
        tif_name = file_name
    else:
        tif_name = file_name.replace(".zip", ".liquid.tif")
    file_path = os.path.join(path, tif_name)
    if not os.path.exists(file_path):
        r = requests.get(url)
        if r.ok:
            os.makedirs(path, exist_ok=True)
            if extension == ".tif":
                with open(file_path, "wb") as f:
                    f.write(r.content)
            else:
                zipped = zipfile.ZipFile(io.BytesIO(r.content))
                zipped.extract(tif_name, path)
        else:
            raise RuntimeError(f"{r.status_code}: could not download {r.url}")
    return file_path


def get_latest_imerg_time(run="E", version="07B"):
    """Returns a pandas time stamp representing the latest available data"""
    now = pd.Timestamp.now(tz="UTC").ceil("3h") + pd.Timedelta(minutes=30)
    for i in range(2, 17):
        t = now - pd.Timedelta(hours=i * 3)
        url = build_imerg_url(t, run=run, version=version)
        if requests.get(url).ok:
            url_t = pd.Timestamp(url[-43:-27].replace("-S", " "))
            return url_t.tz_localize(None)
    raise RuntimeError(f"No IMERG data available at {url}")


def get_valid_SMAP_time(start_time):
    day = start_time.strftime('%Y-%m-%d')
    hours = pd.date_range(day +" 01:30", periods=9, freq="3h")
    return hours[hours.get_indexer([start_time], method = 'nearest')[0]]

def build_smap_url(t, version='7', minor_version='030', opendap=True):
    """Build URL to SMAP data"""
    od_server = 'https://n5eil02u.ecs.nsidc.org/opendap/'
    http_server = 'https://n5eil01u.ecs.nsidc.org/'
    url = (
        f'{od_server if opendap else http_server}SMAP/'
        f'SPL4SMGP.{int(version):03}/{t.year}.{t.month:02}.{t.day:02}/'
        f'SMAP_L4_SM_gph_{t.year}{t.month:02}{t.day:02}'
        f'T{t.hour:02}{t.minute:02}00_Vv{version}{minor_version}_001.h5'
    )
    return url

def download_smap(url, path='./smap'):
    """Downloads and saves SMAP file"""
    file_name = url[url.rfind('/')+1:len(url)]
    file_path = os.path.join(path, file_name)
    if not os.path.exists(file_path):
        r = requests.get(url)
        if r.ok:
            os.makedirs(path, exist_ok=True)
            with open(file_path, 'wb') as f:
                f.write(r.content)
        else: 
            return None
    return file_path

def get_smap(
        start_time, variables=['sm_profile_wetness', 'snow_mass'], 
        load=True, version='7', minor_version='030', opendap=True, 
        latitudes=slice(-60, 60), longitudes=slice(-180, 180)
    ):
    """returns an xarray dataset representing 1 SMAP file"""
    start_time = get_valid_SMAP_time(start_time)
    url = build_smap_url(start_time, version, minor_version, opendap=opendap)
    if opendap:
        try:
            smap = xr.open_dataset(url)
        except OSError:
            logging.warning('SMAP data is not available for'
                ' the specified time and version.')
            for i in range(1, 17): # Look back 2 more days
                earlier_time = start_time - pd.Timedelta(hours=i*3)
                url = build_smap_url(earlier_time, version, minor_version, opendap=opendap)
                #TODO check for availability without triggering warnings
                try: 
                    smap = xr.open_dataset(url)
                except (KeyError, OSError):
                    smap = None
                if smap is not None:
                    break
            if smap is None:
                logging.warning(f'No  SMAP data available since: {earlier_time}')
                return None
            logging.warning(f'Using latest available SMAP data: {earlier_time}')
        # Convert to WGS84
        smap['x'] = smap.cell_lon[0,]
        smap['y'] = smap.cell_lat[:, 0]
    else:
        file_path = download_smap(url)
        if file_path is None:
            logging.warning('SMAP data is not available for'
            ' the specified time and version.')
            for i in range(1, 17):  # Look back 2 more days
                earlier_time = start_time - pd.Timedelta(hours=i*3)
                url = build_smap_url(earlier_time, version, minor_version, opendap=False)
                file_path = download_smap(url)
                if file_path is not None:
                    break
            
            if file_path is None:
                logging.warning(f'No  SMAP data available since: {earlier_time}')
                return None
            logging.warning(f'Using latest available SMAP data: {earlier_time}')
        smap_grid = xr.open_dataset(file_path)
        smap = xr.open_dataset(file_path, group='Geophysical_Data')
        # Convert to WGS84
        smap['x'] = smap_grid.cell_lon[0,].values
        smap['y'] = smap_grid.cell_lat[:, 0].values
    smap = smap[variables]
    smap = smap.sortby('y')
    smap = smap.sel(y=latitudes, x=longitudes)
    if load: 
        smap.load()
    return smap.rename(y='lat', x='lon')

def load_imerg_tiff(
    file_path: str,
    latitudes=slice(60, -60),
    longitudes=slice(-180, 180),
) -> xr.DataArray:
    """Open IMERG geotiff file"""
    imerg = xr.open_dataarray(file_path, engine="rasterio")
    imerg.name = "precipitation"
    selected = imerg.sel(x=longitudes, y=latitudes).squeeze().drop_vars("band")
    renamed = selected.rename({"x": "lon", "y": "lat"})
    return renamed


def get_IMERG_precipitation(
    start_time: pd.Timestamp,
    end_time: pd.Timestamp,
    liquid=True,
    load=True,
    run="E",
    version="07B",
    cache_dir="./imerg",
    latitudes=slice(-60, 60),
    longitudes=slice(-180, 180),
):
    """Opens IMERG data"""
    if end_time <= start_time:
        raise ValueError("End time must be later than start time")
    days = pd.date_range(start_time, end_time, freq="D", inclusive="left")
    urls = [
        build_imerg_url(
            d,
            run=run,
            version=version,
            liquid=liquid,
        )
        for d in days
    ]
    files = [download_imerg(u, path=cache_dir) for u in urls]
    arrays = [
        load_imerg_tiff(
            f,
            # geotiff y (latitude) coordinates don't match other fields
            latitudes=slice(latitudes.stop, latitudes.start),
            longitudes=longitudes,
        )
        for f in files
    ]
    imerg = xr.concat(arrays, dim="time")
    imerg["time"] = days
    precipitation = imerg / 10  # the geotiffs are stored in 0.1 mm
    if load:
        precipitation.load()
    return precipitation


def build_GEOS_url(run_time=None, mode="fcast"):
    """Build string for the GEOS OpenDAP server at the NCCS Data Portal"""
    if mode == 'fcast':
        if run_time is None:
            raise ValueError('GEOS FP model run time must be specified for forecast mode')
        return f'https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/fcast/'\
            f'tavg1_2d_lnd_Nx/tavg1_2d_lnd_Nx.{run_time.year}{run_time.month:02}{run_time.day:02}'\
            f'_{run_time.round("6h").hour:02}'
    if mode == 'assim':
        return 'https://opendap.nccs.nasa.gov/dods/GEOS-5/fp/0.25_deg/assim/tavg1_2d_lnd_Nx'

def get_latest_GEOS_run_time(end_time=None):
    """Returns a pandas time stamp representing the latest available data file"""
    if end_time: 
        now = end_time
    else:
        now = pd.Timestamp.now(tz='UTC')
    for i in range(0, 5):
        latest = now.floor('6h') - pd.Timedelta(hours=i*6)
        url = build_GEOS_url(latest)
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings(action='ignore')
                geos = xr.open_dataset(url)
            if geos: 
                break
        except OSError:
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
    return geos[variable].drop_vars('time')

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
    now = pd.Timestamp.now().strftime('%Y-%m-%dT%H:%M:%S:%fZ')
    data_set.attrs["ProductionDateTime"] = now
    data_set.attrs['title'] = 'Landslide Hazard Analysis for Situational Awareness'
    data_set.attrs['institution'] = 'NASA GSFC'
    data_set.attrs['source'] = 'LHASA V2.1'
    data_set.attrs["history"] = f"{now} File written by XArray version {xr.__version__}"
    if run_mode == 'nrt': 
        data_set.attrs['references'] = (
            'Stanley, T. A., D. B. Kirschbaum, G. Benz, et al. 2021. '
            '"Data-Driven Landslide Nowcasting at the Global Scale." '
            'Frontiers in Earth Science, 9: [10.3389/feart.2021.640043]'
        )
    else: 
        data_set.attrs['references'] = (
            'Khan, S., D. B. Kirschbaum, T. A. Stanley, P. M. Amatya, '
            'and R. A. Emberson. 2022. "Global Landslide Forecasting System '
            'for Hazard Assessment and Situational Awareness." Frontiers in '
            'Earth Science, 10: [10.3389/feart.2022.878996]'
        )
    data_set.attrs['comment'] = 'LHASA identifies where landslides are most probable in near real time.'
    data_set.attrs['Conventions'] = 'CF-1.8'
    data_set.attrs['ShortName'] = 'LHASA'
    data_set.attrs['LongName'] = 'Landslide Hazard Analysis for Situational Awareness (LHASA)'
    data_set.attrs['VersionID'] = '2.1'
    data_set.attrs['Format'] = 'netCDF-4'
    data_set.attrs['DataSetQuality'] = 'NRT'
    data_set.attrs['IdentifierProductDOIAuthority'] = 'https://doi.org/'
    data_set.attrs['IdentifierProductDOI'] = '10.5067/8VKQDQFFOTS3'
    data_set.attrs['ProcessingLevel'] = '4'
    start_date = data_set["time"].dt.strftime("%Y-%m-%d").values[0]
    data_set.attrs["RangeBeginningDate"] = start_date
    data_set.attrs["RangeBeginningTime"] = "00:00:00.000000"
    end_date = (
        data_set["time"] 
        + pd.Timedelta(days=1)
    ).dt.strftime("%Y-%m-%d").values[0]
    data_set.attrs["RangeEndingDate"] = end_date
    data_set.attrs["RangeEndingTime"] = "00:00:00.000000"
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

def save_tiff(data_array: xr.DataArray, file_path: str):
    """Saves prediction in geotiff format"""
    data_array.rio.write_nodata(NO_DATA, inplace=True)
    data_array.rio.write_crs(4326, inplace=True)
    data_array = data_array.rename(lat="latitude", lon="longitude")
    data_array.rio.to_raster(file_path, compress='zstd')


def get_model(file_path, threads=1):
    """Open trained model"""
    model = xgb.Booster()
    model.load_model(file_path)
    model.set_param('nthread', threads)
    return model

def expose(hazard: xr.DataArray, variables: xr.Dataset, 
            thresholds={'l': 0.1, 'm': 0.5, 'h': 0.9}, 
            selection = ['population', 'road_length']):
    """Totals exposed assets at flexible hazard thresholds"""
    variables = variables.copy() # To avoid mutating the assets dataset
    for k in thresholds.keys():
        variables[f'{k}_haz'] = hazard > thresholds[k]
        for v in selection:
            v_name = f"{k}_haz_{'pp' if v == 'population' else 'rd'}"
            variables[v_name] = variables[v] * variables[f'{k}_haz']
    variables = variables.drop_vars(['time', 'lat', 'lon', 'road_length', 'population']).squeeze()
    return variables

def add_ratios(totals, thresholds={'l': 0.1, 'm': 0.5, 'h': 0.9}, 
                selection = ['population', 'road_length']):
    """Calculates the fraction of each asset exposed in each county"""
    for k in thresholds.keys():
            totals[f'{k}_haz_f'] = totals[f'{k}_haz'] / totals['cells']
            for v in selection:
                v_name = f"{k}_haz_{'pp' if v == 'population' else 'rd'}_f"
                totals[v_name] = totals[v_name[:-2]] / totals[v]

def imerg_cleanup(path: str, cache_days=0, cache_end_time=None):
    files = glob.glob(os.path.join(path, 'imerg', '*'))
    if cache_days < 1:
        for f in files:
            os.remove(f)
    else:
        cache_start_time = cache_end_time - pd.DateOffset(days=cache_days)
        df = pd.DataFrame({'file': files})
        df['time'] = df['file'].apply(lambda x: 
            pd.Timestamp(x[-50:-34].replace('-S', ' ')))
        excess = df['file'][~df['time'].between(cache_start_time, cache_end_time)]
        excess.map(os.remove)

def smap_cleanup(path: str, cache_days=0, cache_end_time=None):
    files = glob.glob(os.path.join(path, 'smap', '*'))
    if cache_days < 1:
        for f in files:
            os.remove(f)
    else:
        cache_start_time = cache_end_time - pd.DateOffset(days=cache_days)
        df = pd.DataFrame({'file': files})
        df['time'] = df['file'].apply(lambda x: 
            pd.Timestamp(x[-29:-14].replace('-S', ' ')))
        excess = df['file'][~df['time'].between(cache_start_time, cache_end_time)]
        excess.map(os.remove)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='LHASA 2.0 global landslide forecast')
    parser.add_argument('-v', '--version', action='version', 
        version='LHASA version 2.1')
    parser.add_argument('-p', '--path', default=os.getcwd(), 
        help='location of input files')
    parser.add_argument('-op', '--output_path', help='location of output files')
    parser.add_argument('-d', '--date', 
        help='UTC date and time, formatted as "YYYY-MM-DD HH:MM", assumes '
        'latest available if unspecified. This represents the end time of the'
        'NRT product and the start time of the forecast product.')
    parser.add_argument('-l', '--lead', type=int,  default=2, 
        help='Days of forecast to generate. 0 means only LHASA-NRT will run.')
    parser.add_argument('-rt', '--run_time', 
        help='GEOS-FP model initiation time (UTC), formatted as '
        '"YYYY-MM-DD HH", assumes latest available if unspecified')
    parser.add_argument('-o', '--overwrite', action="store_true", 
        help='assumed false if this option is not provided')
    parser.add_argument('--small', action="store_true", 
        help='shrinks output file size by masking where p_landslide < 0.01')
    parser.add_argument(
        "-st",
        "--slope_threshold",
        type=float,
        default=10.0,
        help="Slope angle (degrees) below which predictions will not be shown",
    )
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
    parser.add_argument('-sv', '--smap_version', default='7031',
        help='SMAP L4 major and minor version, e.g. 7031')
    parser.add_argument('-iv', '--imerg_version', default='07B',
        help='IMERG version, e.g. 07B')
    parser.add_argument(
        "-icd",
        "--imerg_cache_days",
        type=int,
        default=7,
        help="Days of IMERG data to cache.",
    )
    parser.add_argument(
        "-scd",
        "--smap_cache_days",
        type=int,
        default=7,
        help="Days of SMAP data to cache.",
    )
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
        forecast_start_time = (
            get_latest_imerg_time(run='E', version=args.imerg_version)
            + pd.Timedelta(minutes=30)
        )
    if args.lead > 0:
        # GEOS is hourly data, so we may have to discard the last half hour of IMERG
        forecast_start_time = forecast_start_time.floor('h')
        if args.run_time:
            run_time = pd.Timestamp(args.run_time).floor('6h')
            if run_time > forecast_start_time:
                raise ValueError("Forecast must start later than GEOS-FP run time")
        else:
            run_time = get_latest_GEOS_run_time(forecast_start_time)
        if not run_time:
            warnings.warn(f'No GEOS-FP run is available for {forecast_start_time}.')
            args.lead = 0 # Fall back to NRT-only model run
    if args.lead > 2:
        warnings.warn('Only 2 days of forecast data are considered reliable.')

    model = get_model(f'{path}/model.json', args.threads)

    static_files = [f'{path}/static/pga.nc4',
        f'{path}/static/Lithology.nc4',
        f'{path}/static/Slope.nc4',
        f'{path}/static/mask.nc4']
    static_variables = xr.open_mfdataset(static_files, parallel=True)
    static_variables = static_variables.sel(
        lat=slice(args.south, args.north), 
        lon=slice(args.west, args.east)
    )
    p99 = xr.open_dataarray(f'{path}/static/p99.nc4')
    p99 = p99.sel(
        lat=slice(args.south, args.north), 
        lon=slice(args.west, args.east)
    ).load()
    if args.lead > 0:
        p99geos = xr.open_dataarray(f'{path}/static/p99geos.nc4')
        p99geos = p99geos.sel(
            lat=slice(args.south, args.north), 
            lon=slice(args.west, args.east)
        ).load()
    if args.exposure:
        exposure_files = [f'{path}/exposure/road_length.nc4', 
        f'{path}/exposure/population.nc4', f'{path}/exposure/gadm36.nc4']
        assets = xr.open_mfdataset(exposure_files)
        assets = assets.sel(
            lat=slice(args.north, args.south), 
            lon=slice(args.west, args.east)
        )
        constant_totals = pd.read_csv(f'{path}/exposure/totals.csv', index_col='gadm_fid')
        # Apply a minimum to avoid divide by zero errors
        constant_totals['road_length'] = np.maximum(constant_totals['road_length'], 1e-10)
        constant_totals['population'] = np.maximum(constant_totals['population'], 1e-10)
        admin_names = pd.read_csv(f'{path}/exposure/admin_names.csv', index_col='gadm_fid')
    logging.info('opened static variables')

    precipitation_start_date = forecast_start_time - pd.DateOffset(3)
    precipitation_end_date = forecast_start_time + pd.DateOffset(args.lead)

    imerg_late = get_IMERG_precipitation(
        precipitation_start_date,
        forecast_start_time - pd.DateOffset(1),
        liquid=True,
        load=False,
        run="L",
        version=args.imerg_version,
        cache_dir=os.path.join(path, "imerg"),
        latitudes=slice(args.south, args.north),
        longitudes=slice(args.west, args.east),
    )
    imerg_early = get_IMERG_precipitation(
        forecast_start_time - pd.DateOffset(1),
        forecast_start_time,
        liquid=True,
        load=False,
        run="E",
        version=args.imerg_version,
        cache_dir=os.path.join(path, "imerg"),
        latitudes=slice(args.south, args.north),
        longitudes=slice(args.west, args.east),
    )
    imerg = xr.concat([imerg_late, imerg_early], dim="time")

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
            offset=pd.Timedelta(f'{forecast_start_time.hour}h')).sum().load()
        daily_rain = [da for da in imerg] + [da for da in geos_daily]
    else: 
        daily_rain = [da for da in imerg]
        
    dates = pd.date_range(
        forecast_start_time - pd.DateOffset(1), 
        precipitation_end_date - pd.DateOffset(1)
    )
    for i, d in enumerate(daily_rain[2:(3+args.lead)]): 
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
            nc_path = ''
        if 'tif' in args.format:
            tif_path = (f'{output_path}/{run_mode}/hazard/tif/'
                f'{initialization}{date_string}.tif')
            if os.path.exists(tif_path) and not args.overwrite:
                raise FileExistsError(f'{tif_path} exists. Use the --overwrite'
                    ' argument to overwrite it.')
        else: 
            tif_path = ''

        day_before = daily_rain[i]
        yesterday = daily_rain[i + 1].interp_like(day_before, method='nearest')
        antecedent = xr.concat([yesterday, day_before], dim='time').sum('time')
        antecedent.name = 'antecedent'

        if run_mode == "nrt": 
            rain = d.sortby('lat').reindex_like(p99, method='nearest') / p99
            if args.opendap:
                smap_variables = ["Geophysical_Data_sm_profile_wetness"]
            else:
                smap_variables = ["sm_profile_wetness"]
            smap = get_smap(
                start_time=dates[i] - pd.Timedelta(days=2, hours=1), 
                variables=smap_variables, 
                version=args.smap_version[0], 
                minor_version=args.smap_version[1:4], 
                opendap=args.opendap,
                latitudes=slice(args.south, args.north), 
                longitudes=slice(args.west, args.east)
            )
            if smap is None:
                if args.lead > 0:
                    logging.warning('LHASA skipped to the forecast product.')
                    continue
                else:
                    raise SystemExit('LHASA halted because SMAP is unavailable.')
            moisture = smap[smap_variables[0]]
            # Match variable names in GEOS
            moisture.name = 'gwetprof'
        else: 
            moisture = get_GEOS_variable(
                start_time=dates[i] - pd.Timedelta(days=2, hours=1), 
                run_time=run_time, 
                variable='gwetprof', 
                mode='assim', 
                latitudes=slice(args.south, args.north), 
                longitudes=slice(args.west, args.east)
            )
            rain = d.interp_like(p99geos, method='nearest') / p99geos
        rain.name = 'rain'
        dynamic_variables = [rain, antecedent, moisture]
        transposed = [v.transpose('lat', 'lon') for v in dynamic_variables]
        logging.info('opened ' + date_string)
        regridded = xr.merge([regrid(v, static_variables) for v in transposed])
        logging.info('interpolated')
        # It's important to retain the correct variable order
        variable_order = [
            "Lithology",
            "Slope",
            "pga",
            "gwetprof",
            "antecedent",
            "rain",
        ]
        static_variables.load()
        variables = xr.merge([regridded, static_variables.drop_vars('land_mask')])
        masked_values = [apply_mask(variables[v], mask=(static_variables['land_mask'] > 0)) for v in variable_order]
        logging.info('built mask')
        inputs = xgb.DMatrix(
            np.stack(masked_values, 1),
            feature_names=variable_order,
            nthread=args.threads,
        )
        prediction = model.predict(inputs)
        logging.info('made predictions')
        p_landslide = fill_array(prediction, static_variables['land_mask'] > 0, dates[i])
        p_landslide = p_landslide.where(static_variables['Slope'] > args.slope_threshold)
        if args.small:
            warnings.warn(
                ("The --small command-line argument will be removed in a future"
                " version of LHASA. Use the slope threshold instead. E.g.:  "
                "(python lhasa.py -st 15)"),
                category=FutureWarning,
            )
            p_landslide = p_landslide.where(p_landslide >= 0.01)
        if nc_path:
            save_nc(p_landslide.to_dataset(), nc_path)
            logging.info(f'saved {nc_path}')
        if tif_path:
            save_tiff(p_landslide, tif_path)
            logging.info(f'saved {tif_path}')
        
        if args.exposure:
            csv_path = os.path.join(
                output_path, 
                run_mode, 
                'exposure', 
                'csv', 
                f'{date_string}.csv'
            )
            if os.path.exists(csv_path) and not args.overwrite:
                raise FileExistsError(f'{csv_path} exists. Use the --overwrite'
                    ' argument to overwrite it.')
            # The exposure analysis uses a lot of memory
            del dynamic_variables, transposed, regridded, variables
            del inputs
            exposed = expose(p_landslide, assets)
            data_frame = exposed.to_dataframe().dropna()
            totals = data_frame.groupby('gadm_fid').sum()
            totals = totals.merge(constant_totals, on='gadm_fid')
            add_ratios(totals)
            totals = totals.drop(columns=['population', 'road_length', 'cells'])
            logging.info('Calculated exposure levels by county')
            totals = totals[totals['l_haz'] > 0]
            named = admin_names.merge(totals, on='gadm_fid', how='inner')
            named.to_csv(csv_path)
            logging.info(f'saved {csv_path}')
    
    imerg_cleanup(
        path=path, 
        cache_end_time=forecast_start_time, 
        cache_days=args.imerg_cache_days
    )

    smap_cleanup(
        path=path, 
        cache_end_time=precipitation_start_date, 
        cache_days=args.smap_cache_days
    )
