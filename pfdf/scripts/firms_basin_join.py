import glob
import geopandas as gpd
import pandas as pd
from datetime import date
import numpy as np

def read_function(file_path, bbox):
    
    if bbox is None:
        df = pd.read_csv(file_path)
    
    else:
        xmin, xmax, ymin, ymax = bbox[0], bbox[2], bbox[1], bbox[3]
        df = pd.read_csv(file_path)
        df = df[(df['longitude']>xmin) & (df['longitude']<xmax) &
                (df['latitude']>ymin) & (df['latitude']<ymax)]
        
    return df

def join(firms_path, sheds, date_range=65, count_thresh=250, pre_offset=1, bbox=None):
    all_firms = glob.glob(f'{firms_path}/*.txt')
    snpp_paths = [file for file in all_firms if "VNP" in file]
    noaa_paths = [file for file in all_firms if "VNP" not in file]
    
    snpp_paths.sort()
    noaa_paths.sort()
    
    snpp_df = pd.concat([read_function(file,bbox) for file in snpp_paths[-date_range:]],ignore_index=True)
    snpp_df['Sat'] = 'SNPP'
    
    noaa_df = pd.concat([read_function(file,bbox) for file in noaa_paths[-date_range:]],ignore_index=True)
    noaa_df['Sat'] = 'J1'
    
    all_df = pd.concat([snpp_df,noaa_df],ignore_index=True)

    gdf = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(
            all_df.longitude,
            all_df.latitude,
            crs="EPSG:4326",
        ),
        data=all_df,
    )
    
    gdf.to_crs('EPSG:6933',inplace=True)
    
    sheds = sheds.copy()
    sheds.to_crs('EPSG:6933',inplace=True)
    
    joined_gdf = gpd.sjoin(gdf,sheds,predicate='intersects',how='inner')
    
    joined_gdf['acq_date'] = pd.to_datetime(joined_gdf['acq_date'])
    
    
    # counts of daily pixels in each basin
    firms_daily_count = joined_gdf.groupby(['HYBAS_ID','acq_date']).size()

    # return all days that exceed the pixel threshold, sort by increasing basin id and increasing date
    firms_daily_count_f = firms_daily_count[firms_daily_count>count_thresh].sort_index(ascending=True).reset_index()

    # drop entries with duplicate dates for a single basin, keep the first one.
    # a.k.a the earliest data above the pixel threshold
    firms_fire_start_idx = firms_daily_count_f['HYBAS_ID'].drop_duplicates(keep='first').index
    return_df = firms_daily_count_f.loc[firms_fire_start_idx].rename({'acq_date':'FireStart',
                                                                      0:'InitialCount'},axis=1)
    return_df.set_index('HYBAS_ID',inplace=True)

    joined_gdf_copy = joined_gdf.copy().set_index('HYBAS_ID')
    joined_gdf_copy = joined_gdf_copy.loc[joined_gdf_copy.index.isin(return_df.index)] # isolate just the basins from above
    last_detection = joined_gdf_copy.groupby(joined_gdf_copy.index)['acq_date'].max().rename('FireEnd') # get last fire observation

    return_df = return_df.join(last_detection, on='HYBAS_ID') # join!
    return_df['Duration'] = (return_df['FireEnd'] - return_df['FireStart']) / np.timedelta64(1,'M')

    today = pd.to_datetime(date.today())
    return_df['DaysFromLastDetection'] = (today-return_df['FireEnd']) / np.timedelta64(1,'D')

    
    return_df['LastRunDate'] = today

    sheds.set_index('HYBAS_ID',inplace=True)
    return_df['geometry'] = sheds.loc[return_df.index]['geometry']

    run_ids = gpd.GeoDataFrame(return_df,geometry='geometry')
    run_ids = run_ids[(run_ids['DaysFromLastDetection']<200)&(run_ids['DaysFromLastDetection']>14)]

    run_ids['PreOffset'] = pre_offset
    run_ids['Notes'] = None

    return run_ids
    