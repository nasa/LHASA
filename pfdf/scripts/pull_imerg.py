# coding=utf-8
# Copyright © 2021, United States Government, as represented by the
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
Functions and script for global post-fire debris flow model
"""

import logging
import os
import argparse
import numpy as np
import pandas as pd
import geopandas as gpd
import requests
import xarray as xr
from scipy import stats
import xgboost as xgb
from urllib.request import urlopen
import xml.etree.ElementTree as ET
import warnings

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

def download_imerg(url):
    """Downloads and saves IMERG file"""
    file_name = url[url.rfind('/')+1:len(url)]
    file_path = 'imerg/' + file_name
    if not os.path.exists(file_path):
        request = requests.get(url)
        if request.ok:
            os.makedirs('imerg', exist_ok=True)
            with open(file_path, 'wb') as f:
                f.write(request.content)
        else: 
            raise RuntimeError(str(request.status_code) + ': could not download ' + request.url)
    return file_name

def get_latest_imerg_year(run='E', version='06'):
    """Finds the last year IMERG data is available at GES-DISC OpenDAP"""
    url = f'{OPENDAP_URL}/ncml/aggregation/GPM_3IMERGHH{run}.{version}/catalog.xml'
    with urlopen(url) as f:
        catalog = ET.parse(f).getroot()
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
    with urlopen(url) as f:
        catalog = ET.parse(f).getroot()
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
                return t
        raise RuntimeError(f'No IMERG data available at {url}')

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
        half_hours = pd.date_range(start_time, end_time, freq='30min')
        urls = [build_imerg_url(h, run=run, version=version, opendap=False) for h in half_hours]
        files = [download_imerg(u) for u in urls]
        imerg = xr.open_mfdataset(files, parallel=True)
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

def zonal_max(polygon, grid: xr.DataArray):
    """Calculate zonal maximum with xarray/numpy, not rasterstats"""
    if polygon.type == 'Polygon':
        x, y = polygon.exterior.xy
    else:
        coords = [p.exterior.xy for p in polygon.geoms]
        x = np.concatenate([a[0] for a in coords])
        y = np.concatenate([a[1] for a in coords])
    lat = xr.DataArray(y, dims='points')
    lon = xr.DataArray(x, dims='points')
    grid_values = grid.sel(lon=lon, lat=lat, method='nearest')
    return grid_values.max().values

def get_model(file_path: str, threads=1):
    """Open trained model"""
    model = xgb.Booster()
    model.load_model(file_path)
    model.set_param('nthread', threads)
    return model

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s', 
        datefmt='%H:%M')
    logging.info('logging started')
    parser = argparse.ArgumentParser(
        description='LHASA 2.0 global post-fire debris flow model')
    parser.add_argument('-v', '--version', action='version', 
        version='LHASA version 2.0.0a')
    parser.add_argument('-p', '--path', default=os.getcwd(), 
        help='location of input files')
    parser.add_argument('-op', '--output_path', help='location of output files')
    parser.add_argument('-iv', '--imerg_version', default='06C',
        help='IMERG version, e.g. 06C')
    parser.add_argument('-od', '--opendap', action='store_true',
        help='Use OpenDAP to access IMERG; LHASA defaults to HDF5 download')
    parser.add_argument('-N', '--north', type=float, default=60.0, 
        help='maximum latitude (WGS84)')
    parser.add_argument('-S', '--south', type=float, default=-60.0, 
        help='minimum latitude (WGS84)')
    parser.add_argument('-E', '--east', type=float, default=180.0, 
        help='maximum longitude (WGS84)')
    parser.add_argument('-W', '--west', type=float, default=-180.0, 
        help='minimum longitude (WGS84)')
    args = parser.parse_args()

    path = args.path
    if args.output_path:
        output_path = args.output_path
    else: 
        output_path = path

    intensity_stats = xr.open_dataset(f'{path}/ref_data/stats.nc4')
    intensity_stats = intensity_stats.sel(
        lat=slice(args.south, args.north), 
        lon=slice(args.west, args.east)
    ).load()
    accumulation_stats = xr.open_dataset(f'{path}/ref_data/daily_acc_stats.nc4')
    accumulation_stats = accumulation_stats.sel(
        lat=slice(args.south, args.north), 
        lon=slice(args.west, args.east)
    ).load()
    logging.info('opened rainfall statistics files')

    imerg_early_end_time = get_latest_imerg_time() + pd.Timedelta(minutes=30)
    imerg_late_end_time = get_latest_imerg_time(run='L') + pd.Timedelta(minutes=30)
    # These might be different from the latest IMERG time stamp if specified from the CLI
    precipitation_end_time = imerg_early_end_time
    precipitation_start_time = precipitation_end_time - pd.Timedelta(hours=24)

    imerg_late = get_IMERG_precipitation(
        precipitation_start_time, 
        min(precipitation_end_time, imerg_late_end_time), 
        liquid=False, 
        load=True, 
        run='L', 
        version=args.imerg_version, 
        latitudes=slice(args.south, args.north), 
        longitudes=slice(args.west, args.east)
    )

    if imerg_late_end_time < imerg_early_end_time:
        imerg_early = get_IMERG_precipitation(
            max(precipitation_start_time, imerg_late_end_time), 
            precipitation_end_time, 
            liquid=False, 
            load=True, 
            run='E', 
            version=args.imerg_version, 
            latitudes=slice(args.south, args.north), 
            longitudes=slice(args.west, args.east)
        )
        imerg = xr.concat([imerg_late, imerg_early], dim='time')
    else:
        imerg = imerg_late

    logging.info('opened rain file')

    rain_max = imerg.max('time')
    accumulation = imerg.mean('time') * 24
    logging.info('calculated max and sum of daily rain')

    max_perc_array = stats.norm.cdf(
        np.log(
            rain_max, 
            where=rain_max > 0, 
            out=np.zeros_like(rain_max) - 1e9
        ), # -1e9 will typically produce a quantile of zero in cases where rain==0
        loc=intensity_stats['mean'], 
        scale=intensity_stats['std']
    )
    max_perc = xr.DataArray(
        max_perc_array, 
        coords=intensity_stats['mean'].coords, 
        dims=intensity_stats['mean'].dims
    )
    acc_perc_array = stats.norm.cdf(
        #TODO: replace np.log10 with np.log when IMERG V7 is released
        np.log10(
            accumulation, 
            where=accumulation > 0, 
            out=np.zeros_like(accumulation) - 1e9
        ), # using -1e9 will typically produce a quantile of zero when rain==0
        loc=accumulation_stats['mean_acc'], 
        scale=accumulation_stats['std_acc']
    )
    acc_perc = xr.DataArray(
        acc_perc_array, 
        coords=accumulation_stats['mean_acc'].coords, 
        dims=accumulation_stats['mean_acc'].dims
    )
    logging.info('calculated rainfall percentiles')

    #TODO: Determine what the status of the basins would have been in the past
    basins = pd.read_csv(f'{path}/ref_data/basin_records.csv')
    #TODO: Implement bounding box for basins
    basins = basins[basins['Visibility']>0.95]
    logging.info(str(len(basins))+' basins sit above the visibility threshold.')
    logging.info('opened burned basin file')
    watersheds = gpd.read_file(f'{path}/ref_data/sheds.geojson')
    logging.info('opened basin geometry file')
    selected_basins = watersheds[watersheds['HYBAS_ID'].isin(basins['HYBAS_ID'])]
    # fix issue where basin_records.csv contains basins not in sheds.geojson
    basins = basins[basins['HYBAS_ID'].isin(selected_basins['HYBAS_ID'])].copy()

    basins['MaxPercentile'] = np.array([zonal_max(b, max_perc) for b in selected_basins['geometry']])
    basins['SumPercentile'] = np.array([zonal_max(b, acc_perc) for b in selected_basins['geometry']])
    basins['Slope_Burn_Count'] = basins['Slope_Burn_Abv_23']
    logging.info('calculated zonal maximum rainfall')

    model = get_model(f'{path}/ref_data/model.json')

    predictors = xgb.DMatrix(basins[[
        'MaxPercentile', 
        'SumPercentile',
        'Mean_dNBR',
        'Slope_Burn_Count',
        'SlopeBurnAreaRatio'
    ]])
    basins['p_debris_flow'] = model.predict(predictors)
    logging.info('calculated probabilities')
    basins_geom = gpd.GeoDataFrame(
        data=basins, 
        geometry=selected_basins['geometry'].values, 
        crs=selected_basins.crs
    )
    
    tstamp_f = precipitation_start_time.strftime('%Y%m%dT%H%M')
    try:
        os.mkdir(output_path)
    except FileExistsError:
        pass
    basins_geom.to_file(
        f'{output_path}/p_debris_flow{tstamp_f}.geojson', 
        driver='GeoJSON'
    )
    logging.info('saved output file')
