# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 08:48:15 2022

@author: eorland
"""

import geopandas as gpd
from shapely.geometry import Point
import pandas as pd
import os
import helper_functions

def join(path, firms_path):
    ref_data_path = os.path.join(path, 'ref_data')
    all_files = helper_functions.get_file_list(firms_path)

    # create single df of all FIRMS data
    firms_df = pd.concat([pd.read_csv(f) for f in all_files]).reset_index() 

    # create geometry col via Shaply Point constructor
    firms_df['geometry'] = firms_df.apply(lambda x: 
                                          Point((float(x.longitude), 
                                                 float(x.latitude))), axis=1)
    # create geodataframe
    firms_gdf = gpd.GeoDataFrame(firms_df, geometry='geometry',crs="EPSG:4326")

    print('firms gdf complete')
    # set index
    firms_gdf['firms_id'] = firms_gdf.index

    # load sheds
    sheds_f = gpd.read_file(os.path.join(ref_data_path,'sheds.geojson'))
    print('load_sheds_complete')
    
    # do spatial join of FIRMS data to watersheds - takes a while
    joined_sheds = gpd.sjoin(firms_gdf,sheds_f,predicate='intersects',how='inner')
    print('join complete')
    
    f_basins, f_count, f_detect = helper_functions.get_firms_basins(joined_sheds, 100)

    helper_functions.get_all_basins(f_basins, f_count, f_detect, path)
    
    r_data = pd.read_csv(os.path.join(ref_data_path,'run_data.csv'),
                         index_col='HYBAS_ID')
    basins = r_data.index.values.tolist()
    run_dates = pd.to_datetime(r_data['FirstRunDate']).astype(str).tolist()
    detect = r_data['DaysFromFirstDetection'].values.tolist()
    assert len(basins)==len(run_dates)==len(detect)
    print('total basins for run:', len(basins))
    return basins, run_dates, detect