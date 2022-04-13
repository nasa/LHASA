# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 08:28:36 2022

@author: eorland
"""
import os
import pandas as pd
import ee

def directory_setup(archive_path, model_output_path):
    """build all archive paths"""
    archive_dir_list = [os.path.join(archive_path,'firms'),
                        os.path.join(archive_path,'firms','all_firms'),
                        os.path.join(archive_path,'asset_downloads')] 

    for folder in archive_dir_list:
        os.chdir(archive_path)
        if os.path.exists(folder):
            print(os.path.basename(folder),'already exists.')
            continue
        print("Making directory",os.path.basename(folder))
        os.mkdir(folder) 
        
    model_dir_path = os.path.join(model_output_path,'model_outputs')
    os.chdir(model_output_path)
    
    if os.path.exists(model_dir_path):
            print(os.path.basename(model_dir_path),'already exists.')
            return 
    else: 
        print("Making directory",os.path.basename(model_dir_path))
        os.mkdir(model_dir_path) 
        return 

def make_records(path):
    
    basin_recs_path = os.path.join(path,'ref_data',
                                   'basin_records.csv')
    task_recs_path = os.path.join(path,'ref_data',
                                   'task_records.csv')
    run_data_path = os.path.join(path,'ref_data',
                                   'run_data.csv')
    
    if os.path.exists(basin_recs_path): # if one exists, then they all do
        print('Records files already exist.\n')
        return
    
    
    basin_records_cols = ['HYBAS_ID', 'FirstFirms', 'RecentFirms', 
                          'FirstRunDate','LastRunDate', 'DaysFromFirstDetection',
                          'DaysFromLastDetection', 'Mean_dLAI','Mean_dNBR', 'Mean_dNBR_cnt',
                          'Slope_Burn_Abv_23','SlopeBurnAreaRatio', 'Visibility', 'LastBurnUpdate']
    
    task_records_cols = ['File','OpID','TaskID','AssetID','NumBasins',
                         'Status', 'Message', 'StartTime',
                         'EndTime']
    
    basin_rec = pd.DataFrame(columns=basin_records_cols)
    task_rec = pd.DataFrame(columns=task_records_cols)
    
    print("Writing records files...\n")
    basin_rec.to_csv(basin_recs_path,index=False)
    basin_rec.to_csv(run_data_path,index=False)
    task_rec.to_csv(task_recs_path,index=False)
    
    return 
        
    
    
def gee_folder_setup(user):    
    
    ee.data.createAsset({'type': 'Folder'}, 
                            'users/'+user+'/pfdf_basin_assets')
    print('Asset folder created at users/'+user+'/pfdf_basin_assets')
    return
