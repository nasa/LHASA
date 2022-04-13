import os
import asset_functions, helper_functions
import sys


def query_and_download(path,asset_path,userid):
    task_records_path = os.path.join(path,'ref_data','task_records.csv')
    basin_records_path = os.path.join(path,'ref_data','basin_records.csv')
    asset_records_path = os.path.join(asset_path,'asset_downloads')
    dwn_list, data_update = helper_functions.update_task_list(task_records_path)

    if len(dwn_list)==0:
        print('No new completed tasks.')
    else:
        for file in dwn_list:
            print('Downloading',file+'...')
            
            try:
                asset_subpath = asset_functions.download_large_asset(file,userid,asset_records_path,
                                                                   delete_intermediate=True)
            except: 
                print('Could not download file',file)
                e = sys.exc_info()
                print('\nThe following exception occured:',e)
                continue
            try:
                helper_functions.update_dnbr_data(asset_subpath,basin_records_path)
                asset_functions.delete_asset_gee(file,asset_subpath,userid)
            except: 
                print('Downloaded but could not update asset data for file', file)
                e = sys.exc_info()
                print('The following exception occured:',e)
                continue
            
      
    data_update.to_csv(task_records_path,header=True,index=True)
    
    return