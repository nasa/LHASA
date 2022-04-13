import os
import ee
import urllib.request
import glob
import pandas as pd
   
def download_large_asset(file,userid,asset_dir,ext='csv',
                         delete_intermediate=True):
    
    file_dir_path = os.path.join(asset_dir,file)
    
    try:
        os.mkdir(file_dir_path)
        os.chdir(file_dir_path)
    
    except:
        os.chdir(file_dir_path)
        file_path = os.path.join(os.getcwd()+'/'+file+'.'+ext)
        if os.path.exists(file_path):
            print('file already exists in', file_path)
            return file_path 
    
    
    table = ee.FeatureCollection('users/'+userid+'/pfdf_basin_assets/'+file);
    table_size = table.size().getInfo()
    
    if table_size == 0:
        print('Table is empty.')
        link = table.getDownloadURL(ext)
        urllib.request.urlretrieve(link,file_path)
        return file_path
    
    interval = int(table_size / 10)
    
    if interval==0:
        interval = 1
        
    while interval > 500:
        interval = interval - 25
        
    path_list = []
    for i in range(interval,table_size+interval,interval):
        
        subset = ee.FeatureCollection(table.toList(interval,(i-interval)))
        sub_id = '_'+ str(i-interval)
        link = subset.getDownloadURL(ext)
        file_path = os.path.join(os.getcwd()+'/'+file+sub_id+'.'+ext)
        path_list.append(file_path)
        urllib.request.urlretrieve(link,file_path)

    all_filenames = [i for i in glob.glob('*.{}'.format(ext))]
    df = pd.concat([pd.read_csv(f) for f in all_filenames]).reset_index()
    df_path = os.path.join(os.getcwd()+'/'+file+'.'+ext)  
    df.to_csv(df_path,header=True,index=True)
    
    if delete_intermediate:
        for file in path_list:
            os.remove(file)
    
    return df_path


def delete_asset_gee(file,file_path,userid):
    assert os.path.exists(file_path), "File path doesn't exist."
    asset_id = 'users/'+userid+'/pfdf_basin_assets/'+file
    ee.data.deleteAsset(asset_id)
    print("Asset Deleted")
    return 

