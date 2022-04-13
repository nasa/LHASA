# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 08:24:19 2022

@author: eorland
"""

import numpy as np
import ee
import ee_connect
import pandas as pd
import os
from datetime import datetime
import argparse

def export(basin_list,path,userid):
    '''Export workflow consolitated into a function.'''
    
    t = datetime.utcnow()
    tstamp_f = t.strftime('%Y_%m_%dT%H%M')
    
    file = 'pfdf_export_'+tstamp_f
    assetid = 'users/' + userid + '/pfdf_basin_assets/' + file
    shed_list = ee.List(basin_list)
    sheds_w_dnbr = ee.FeatureCollection(shed_list.map(get_all_data))
    
    task = ee.batch.Export.table.toAsset(collection=sheds_w_dnbr,
                                         assetId=assetid,
                                         description=file)
    
    print("Starting export for "+file+'...')
    task.start()
    
    op_str = 'projects/earthengine-legacy/operations/'
    task_id = task.status()['id']
    opid = op_str+task_id
    op = ee.data.getOperation(opid)
    state = op['metadata']['state']
    
    cols = ['File','OpID','TaskID','AssetID',
            'NumBasins', 'Status', 
            'Message','StartTime','EndTime']
    write_array = np.array([file, opid, task_id, assetid, len(basin_list),
                            state,'NA','NA','NA']).reshape(1,9)
    write_vals = pd.DataFrame(write_array,columns=cols)
    write_vals.set_index('File',inplace=True)
    
    task_records_path = os.path.join(path,'ref_data','task_records.csv')
    data = pd.read_csv(task_records_path,index_col='File')
    data_ap = data.append(write_vals,ignore_index=False)
    data_ap.to_csv(task_records_path,header=True,index=True)
    
    return task

def applyScaleFactors(image):
    opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2)
    thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0)
    
    return image.addBands(opticalBands, None, True).addBands(thermalBands, None, True);

def get_all_data(hybas_id):
    
    h_id = hybas_id
    
    h_idx = shed_list.indexOf(h_id)
    
    run_date = ee.Date(rd_list.get(h_idx)) # run date for basin
    days_from_det = ee.Number(d_list.get(h_idx)) # days before run date 
    
    feature = ee.Feature(watersheds.filter(
        ee.Filter.eq('HYBAS_ID',h_id)).first())
    geometry = feature.geometry()
    
    def lsNBR(lsImage):
    
        nbr = lsImage.normalizedDifference(['SR_B5', 'SR_B7']).toFloat()
        qa = lsImage.select(['QA_PIXEL'])
  
        return nbr.addBands([qa]).select(
                [0,1], ['nbr', 'QA_PIXEL']).copyProperties(
                                            lsImage, ['system:time_start']);
          
    def lsmask(lsImage): # updated cloud mask for LS8 Col 2
    
        quality = lsImage.select(['QA_PIXEL'])
        clear = quality.bitwiseAnd(2).eq(0).And( # dilated cloud
                quality.bitwiseAnd(4).eq(0).And( # cirrus (high confidence)
                quality.bitwiseAnd(8).eq(0).And( # cloud 
                quality.bitwiseAnd(16).eq(0).And( # cloud shadow 
                quality.bitwiseAnd(32).eq(0).And( # snow
                quality.bitwiseAnd(128).eq(0)))))) # water
                
        return lsImage.updateMask(
            clear).select([0]).copyProperties(
                lsImage, ['system:time_start'])
    
    def lscloudmask(lsImage):
    
        quality = lsImage.select(['QA_PIXEL'])
        clear = quality.bitwiseAnd(2).eq(0).And( # dilated cloud
                quality.bitwiseAnd(4).eq(0).And( # cirrus (high confidence)
                quality.bitwiseAnd(8).eq(0).And( # cloud 
                quality.bitwiseAnd(16).eq(0))))
                
        return lsImage.updateMask(
            clear).select([0]).copyProperties(
                lsImage, ['system:time_start'])
    
    ls8 = ls8SR.map(applyScaleFactors).map(lsNBR).map(lsmask)
    ls8_c = ls8SR.map(lscloudmask) # making separate col for looking at % cloudiness
    
    # start right after the initial last fire detection date
    start = run_date.advance(days_from_det.multiply(-1),'days')
    
    # end at the conclusion of available imagery
    end = ls8.limit(1, 'system:time_start', False).first().date()
        
    pre_NBR = ls8.filterDate(
                        start.advance(-12,'months'),
                        start.advance(-1,'month')).median()
    post_NBR = ls8.filterDate(start,end).median()
    
    no_band_test = ee.Algorithms.IsEqual(
                   0,post_NBR.bandNames().size()) 

    dNBR = ee.Image(ee.Algorithms.If(no_band_test,ee.Image.constant(0).rename('dnbr'),
                            ee.Image(pre_NBR.subtract(post_NBR).rename('dnbr'))))
    
    cond = ee.Image(dNBR).gt(0.1).multiply(slope.gt(23))              
    
    Mean_dNBR = ee.Image(dNBR).updateMask(cond).reduceRegion(
                       reducer=ee.Reducer.mean().combine(
                           reducer2=ee.Reducer.count(),
                           sharedInputs=True),
                       scale=30,
                       geometry=geometry,
                       maxPixels=1e13,
                       tileScale=16)
    
    slope_clipped = slope.clip(geometry)
    
    all_bands = dNBR.addBands(slope_clipped,['slope'])
    
    all_bd_reduced = ee.Image(all_bands).reduceNeighborhood(
                                              reducer=ee.Reducer.mean(),
                                              kernel=ee.Kernel.square(3)
                                              )
    
    slope_mask = ee.Image(all_bd_reduced).select('slope_mean').gt(23)
    
    slope_burn_mask = ee.Image(all_bd_reduced).select(
        'slope_mean').gt(23).multiply(
            ee.Image(all_bd_reduced).select('dnbr_mean').gt(0.4))
            
    slope_abv_23 = ee.Number(
                          ee.Image(all_bd_reduced).select(
                              'slope_mean').updateMask(
                                  slope_mask).reduceRegion(
                                      reducer=ee.Reducer.count(),
                                      scale=30,
                                      geometry=geometry,
                                      maxPixels=1e13,
                                      tileScale=16
                                      ).get('slope_mean'))
                                      
    slope_burn_abv_23 =  ee.Number(
                              ee.Image(all_bd_reduced).select(
                                  'slope_mean').updateMask(
                                  slope_burn_mask).reduceRegion(
                                      reducer=ee.Reducer.count(),
                                      scale=30,
                                      geometry=geometry,
                                      maxPixels=1e13,
                                      tileScale=16
                                      ).get('slope_mean'))
                                      
    slope_burn_area_ratio = ee.Number(
                                    slope_burn_abv_23).divide(
                                        ee.Number(slope_abv_23))
    
    pre_cloud = ls8_c.filterDate( # cloud masked imagery - pre
                        start.advance(-12,'months'),
                        start.advance(-1,'month')).median()
    
    post_cloud = ls8_c.filterDate(start,end).median() # cloud masked - post
    
    # if: no available imagery, set blank image - else: take the difference 
    cloud_diff = ee.Image(ee.Algorithms.If(no_band_test,ee.Image(),
                            ee.Image(pre_cloud.subtract(post_cloud))))
   
    # isolate imagery to high slope locations
    cloud_slope_mask = slope.gt(23)
    
    # apply mask to cloud and slope rasters
    cloud_and_slope = cloud_diff.updateMask(cloud_slope_mask).addBands(
                                slope.updateMask(cloud_slope_mask)).rename(['Clouds','Slopes'])
    
    # count number of clear AND high slope pixels and just high slope pixels
    clouds_red = ee.Image(cloud_and_slope).reduceRegion(
                       reducer=ee.Reducer.count(),
                       scale=30,
                       geometry=geometry,
                       maxPixels=1e13,
                       tileScale=16)
    
    # get ratio of pixel counts - 1 means completely clear, 0 means completely obscured
    ratio = ee.Number(clouds_red.get('Clouds')).divide(ee.Number(clouds_red.get('Slopes')))

    # last step - get LAI info
    
    pre_fire_lai = LAI.filterDate(start.advance(-12,'months'),
                                    start.advance(-1,'month')).max()\
                                    .clip(geometry)
    post_fire_lai = LAI.filterDate(start,end).max()\
                                    .clip(geometry)
    pre_lai_bool = ee.Algorithms.IsEqual(0,pre_fire_lai.bandNames().size())
    post_lai_bool = ee.Algorithms.IsEqual(0,post_fire_lai.bandNames().size())
    
    # nested if statements. if no pre-fire lai data, return all zeros
    # also if pre-fire lai but no post-fire lai, still make it all zeros
    # if both are present, the subtract post fire lai from pre fire lai 
    dLAI = ee.Image(ee.Algorithms.If(pre_lai_bool,\
                              ee.Image.constant(0),\
                              (ee.Algorithms.If(post_lai_bool,\
                                                ee.Image.constant(0),\
                                                pre_fire_lai.divide(post_fire_lai)))))
    mean_dLAI = dLAI.updateMask(dLAI.gt(1).And(slope.gt(23))).reduceRegion(
                reducer=ee.Reducer.mean(),
                crs=LAI.first().projection(),
                geometry=geometry,
                maxPixels=1e13,
                tileScale=16
                ).values().get(0)
    
    mean_dLAI_isNull = ee.Algorithms.IsEqual(mean_dLAI,None)
    mean_dLAI_final = ee.Algorithms.If(mean_dLAI_isNull,0,mean_dLAI)
    
    
    return feature.set('Mean_dNBR',Mean_dNBR.get('dnbr_mean')) \
                  .set('Mean_dNBR_cnt',Mean_dNBR.get('dnbr_count')) \
                  .set('Slope_Burn_Abv_23',slope_burn_abv_23) \
                  .set('SlopeBurnAreaRatio',slope_burn_area_ratio) \
                  .set('Visibility', ratio) \
                  .set('CountClouds', clouds_red.get('Clouds')) \
                  .set('CountSlopes', clouds_red.get('Slopes')) \
                  .set('Mean_dLAI',mean_dLAI_final) \
                  .set('StartDate',start.format())


if __name__ == '__main__':
    
    ee_connect.connect()
    parser = argparse.ArgumentParser(description='Handles directories for user.')
    parser.add_argument('--filepath', type=str, required=True,
                   help='Full file path of project folder')
    parser.add_argument('--gee_username', type=str, required=True,
                   help='Case sensitive username for GEE account.')

    args = parser.parse_args()
    main_path = args.filepath
    username = args.gee_username
    
    basin_records_path = os.path.join(main_path,'ref_data','basin_records.csv')
    # read in all basin data
    r_data = pd.read_csv(basin_records_path,index_col='HYBAS_ID')
    basins = r_data.index.values.tolist()
    run_dates = pd.to_datetime(r_data['FirstRunDate']).astype(str).tolist()
    detects = r_data['DaysFromFirstDetection'].values.tolist()
    assert len(basins)==len(run_dates)==len(detects), 'List lengths do not match.'
    
    print('num basins for export:', len(basins))
    
    # get_all_data() will look for these three lists -- 
    # these variable names are hard-coded into it.
    shed_list = ee.List(basins)
    rd_list = ee.List(run_dates)
    d_list = ee.List(detects)

    watersheds = ee.FeatureCollection(
        "WWF/HydroSHEDS/v1/Basins/hybas_12")
    elevation = ee.Image('NASA/NASADEM_HGT/001')
    ls8SR = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
    slope = ee.Terrain.slope(elevation)
    LAI = ee.ImageCollection('MODIS/006/MCD15A3H').select('Lai')
    
    task = export(basins, main_path, username)