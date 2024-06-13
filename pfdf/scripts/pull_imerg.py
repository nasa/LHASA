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

import os
import warnings
import requests
import xml.etree.ElementTree as ET
import pandas as pd
import xarray as xr


OPENDAP_URL = "https://gpm1.gesdisc.eosdis.nasa.gov/opendap"
PPS_URL = "https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/"


def build_imerg_url(start_time, run="E", version="07B", opendap=True):
    """Build URL to IMERG data"""
    if opendap:
        product = f"GPM_3IMERGHH{run}.{version[:2]}"
        url = (
            f"{OPENDAP_URL}/ncml/aggregation/{product}/{start_time.year}/"
            f"{product}_Aggregation_{start_time.year}{start_time.dayofyear:03}"
            ".ncml"
        )
    else:
        end = start_time + pd.Timedelta(minutes=29, seconds=59)
        product = f"GPM_3IMERGHH{run}.{version}"
        url = (
            f'{PPS_URL}{"early" if run=="E" else "late"}'
            f"/{start_time.year}{start_time.month:02}/"
            f"3B-HHR-{run}.MS.MRG.3IMERG."
            f"{start_time.year}{start_time.month:02}{start_time.day:02}-S"
            f"{start_time.hour:02}{start_time.minute:02}{start_time.second:02}"
            f"-E{end.hour:02}{end.minute:02}{end.second:02}."
            f"{start_time.hour*60+start_time.minute:04}.V{version}.RT-H5"
        )
    return url


def download_imerg(url, path="imerg"):
    """Downloads and saves IMERG file"""
    file_name = url[url.rfind("/") + 1 : len(url)]
    file_path = os.path.join(path, file_name)
    if not os.path.exists(file_path):
        request = requests.get(url)
        if request.ok:
            os.makedirs(path, exist_ok=True)
            with open(file_path, "wb") as f:
                f.write(request.content)
        else:
            raise RuntimeError(
                str(request.status_code) + ": could not download " + request.url
            )
    return file_path


def get_latest_imerg_year(run="E", version="07"):
    """Finds the last year IMERG data is available at GES-DISC OpenDAP"""
    url = f"{OPENDAP_URL}/ncml/aggregation/GPM_3IMERGHH{run}.{version}/catalog.xml"
    catalog = ET.fromstring(requests.get(url).content)
    name_space = {
        "thredds": "http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0"
    }
    for dataset in reversed(
        catalog.findall("thredds:dataset/thredds:catalogRef", name_space)
    ):
        year = dataset.attrib["name"]
        if year.isnumeric():
            return int(year)
    raise RuntimeError(f"cannot parse imerg catalog at {url}")


def get_latest_imerg_url(run="E", version="07"):
    """Returns a path to the latest 30-minute IMERG data at GES-DISC OpenDAP"""
    version = str(version)[:2].zfill(2)
    year = get_latest_imerg_year(run=run, version=version)
    url = (
        f"{OPENDAP_URL}/ncml/aggregation/GPM_3IMERGHH{run}.{version}/{year}"
        "/catalog.xml"
    )
    catalog = ET.fromstring(requests.get(url).content)
    name_space = {
        "thredds": "http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0"
    }
    datasets = reversed(catalog.findall("thredds:dataset/thredds:dataset", name_space))
    for dataset in datasets:
        day_url = f"{OPENDAP_URL}{dataset[2].attrib['urlPath']}"
        try:  # Some catalog entries are empty
            xr.open_dataset(day_url)
        except (KeyError, OSError):
            continue
        return day_url
    raise RuntimeError(f"Cannot parse imerg catalog at {url}")


def get_latest_imerg_time(run="E", version="07B", opendap=True):
    """Returns a pandas time stamp representing the latest available data"""
    if opendap:
        url = get_latest_imerg_url(run=run, version=version)
        imerg = xr.open_dataset(f"{url}")
        latest = imerg["time"].max()
        t = pd.Timestamp(
            year=int(latest.dt.year),
            month=int(latest.dt.month),
            day=int(latest.dt.day),
            hour=int(latest.dt.hour),
            minute=int(latest.dt.minute),
        )
        return t
    else:
        now = pd.Timestamp.now(tz="UTC").floor("30min")
        for i in range(6, 48):
            t = now - pd.Timedelta(hours=i / 2)
            url = build_imerg_url(t, run=run, version=version, opendap=False)
            if requests.get(url).ok:
                return t.tz_localize(None)
        raise RuntimeError(f"No IMERG data available at {url}")


def get_IMERG_precipitation(
    start_time: pd.Timestamp,
    end_time: pd.Timestamp,
    liquid=True,
    load=True,
    run="E",
    version="07B",
    opendap=True,
    latitudes=slice(-60, 60),
    longitudes=slice(-180, 180),
):
    """Opens IMERG data"""
    if end_time <= start_time:
        raise ValueError("End time must be later than start time")
    if opendap:
        days = pd.date_range(start_time, end_time, freq="D")
        files = [build_imerg_url(d, run=run, version=version) for d in days]
        imerg = xr.open_mfdataset(files, parallel=True)
    else:
        half_hours = pd.date_range(
            start_time, end_time - pd.Timedelta(minutes=30), freq="30min"
        )
        urls = [
            build_imerg_url(h, run=run, version=version, opendap=False)
            for h in half_hours
        ]
        files = [download_imerg(u) for u in urls]
        imerg = xr.open_mfdataset(files, group="Grid")
    with warnings.catch_warnings():
        warnings.filterwarnings(action="ignore", category=RuntimeWarning)
        warnings.filterwarnings(action="ignore", category=FutureWarning)
        imerg["time"] = imerg.indexes["time"].to_datetimeindex()
    imerg = imerg.sel(time=slice(start_time, end_time))
    imerg = imerg.sel(lat=latitudes, lon=longitudes)
    if liquid:
        precipitation = imerg["precipitationCal"] * (
            imerg["probabilityLiquidPrecipitation"] > 0.5
        )
    else:
        precipitation = imerg["precipitationCal"]
    if load:
        precipitation.load()
    return precipitation


def pull_imerg(
    home_path,
    output_path=None,
    north=60.0,
    south=-60.0,
    east=180,
    west=-180,
    version="07B",
    opendap=False,
):
    """Loads IMERG data and runs prediction model"""
    if output_path == None:
        output_path = home_path

    imerg_early_end_time = get_latest_imerg_time(opendap=opendap) + pd.Timedelta(
        minutes=30
    )
    imerg_late_end_time = get_latest_imerg_time(
        run="L", opendap=opendap
    ) + pd.Timedelta(minutes=30)
    # These might be different from the latest IMERG time stamp if specified from the CLI
    precipitation_end_time = imerg_early_end_time
    precipitation_start_time = precipitation_end_time - pd.Timedelta(hours=24)
    imerg_late = get_IMERG_precipitation(
        precipitation_start_time,
        min(precipitation_end_time, imerg_late_end_time),
        liquid=False,
        load=True,
        run="L",
        version=version,
        opendap=opendap,
        latitudes=slice(south, north),
        longitudes=slice(west, east),
    )
    if imerg_late_end_time < imerg_early_end_time:
        imerg_early = get_IMERG_precipitation(
            max(precipitation_start_time, imerg_late_end_time),
            precipitation_end_time,
            liquid=False,
            load=True,
            run="E",
            version=version,
            opendap=opendap,
            latitudes=slice(south, north),
            longitudes=slice(west, east),
        )
        imerg = xr.concat([imerg_late, imerg_early], dim="time")
    else:
        imerg = imerg_late
    rain_max = imerg.max("time", keep_attrs=True)
    rain_max.attrs['start'] = imerg.indexes['time'][0].strftime("%Y%m%dT%H%M")
    return rain_max
