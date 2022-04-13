import glob
from datetime import date
from time import sleep
import ee
import pandas as pd
import numpy as np
import os

def get_file_list(path):
    '''Returns a list of file names of the recent
    FIRMS data. This assumes download_firms.py has already run, 
    else the file list won't be the most recent.
    
    FIRMS files are organized by the year + day of year (doy),
    e.g., 2021046 means the 46th day of 2021. We want the last
    month of FIRMS data, so we want all files with today's doy-31
    to today's doy in them.
    '''
    
    # get todays date and year
    doy = date.today().timetuple().tm_yday
    #doy = 6
    year = date.today().year
    day_list = [] # list of str(year+doy)
    all_filenames = []
    
    # start with special case where the last month falls into
    # the previous year
    if doy < 31: 
        diff = 31 - doy # how many days remain in previous year
        for i in range(1,doy): 
            while len(str(i)) < 3: # zero pad for 3 digits total
                    i = '0'+str(i)
            string = str(year)+str(i)
            day_list.append(string)
        for i in range((366-diff),367): # 367 to include leap yr
            string = str(year-1)+str(i)
            day_list.append(string)
    else: 
        for i in range((doy-31),doy):
            if len(str(i))<3:
                while len(str(i)) < 3:
                    i = '0'+str(i)
            string = str(year)+str(i)    
            day_list.append(string)
    for day in day_list:
        extension = str(day)+'.txt'
        files = glob.glob(os.path.join(path, '*{}'.format(extension)))
        all_filenames.append(files) 
    
    # return flattened list
    return [file for file_list in all_filenames for file in file_list]

def get_firms_basins(joined_gdf,count=3,days=14):

    today = pd.to_datetime(date.today())
    
    joined_gdf['acq_date'] = pd.to_datetime(joined_gdf['acq_date'])
    
    # num days from detection
    joined_gdf['time_from_aq'] = (today-joined_gdf['acq_date']) / np.timedelta64(1,'D')
    
    # num firms pts for each basin
    firms_count = joined_gdf.groupby('HYBAS_ID').size()
    
    # most recent detection for each basin
    day_last_detected = joined_gdf.groupby('HYBAS_ID').agg({'time_from_aq': np.min})

    # basins which exceed FIRMS count threshold
    count_cond = firms_count[firms_count>count].index.tolist()
    # basins which exceed FIRMS detected threshold
    detect_cond = day_last_detected[day_last_detected>days].dropna().index.tolist()

    # unique basins for which both conditions are met
    basins_in_common = set(count_cond)&set(detect_cond)
    
    # return basins, their counts, and days from last detections
    return list(basins_in_common), firms_count.loc[basins_in_common], day_last_detected.loc[basins_in_common]

def get_all_basins(f_basins, f_count, f_detect, main_path):
    '''function which does two things: 
    
    1. Looks at the current basin records, determines which 
    had too many clouds during the last run and thus should be
    rerun. The corresponding basin records are updated.
    
    2. Looks at the recent FIRMS basin list, finds which 
    basins are already documented and marks them for rerun; adds 
    new basin entries for completely new basins.
    
    After basin records are updated or added, they are 
    are saved. Also a list of a basins to be run is compiled and 
    the record data that correspond to this list is saved in a
    a separate file. This file saves the basin info for the GEE 
    script to access, and is overwritten just to include the basin
    data for the most recent run.
    '''
    
    print('INFO FOR BASIN RECORDS: \n')
    b_data_path = os.path.join(main_path,'ref_data','basin_records.csv')
    r_data_path = os.path.join(main_path,'ref_data','run_data.csv')
    
    # read existing data
    b_data = pd.read_csv(b_data_path,index_col='HYBAS_ID')

    # first tackle existing data w/ no or low cloud cover
    nan_basins = b_data['Visibility'].index[
        b_data['Visibility'].apply(np.isnan)].values
    #print(nan_basins)
    below_thresh = b_data[b_data['Visibility']<0.95].index
    #print(below_thresh)
    
    cc_basins = list(set(np.concatenate([below_thresh,nan_basins]))) # all cc basins
    print('Basins below visibility threshold:',len(cc_basins))
    cc_run_dates = b_data.loc[cc_basins,'FirstRunDate'].values.tolist()
    b_data.loc[cc_basins,'LastRunDate'] = pd.to_datetime(date.today()) # mark rerun
    #print(b_data.loc[cc_basins])
    
    # aggregate recent firms data, convert to DF
    barray = np.array([f_basins, f_count, f_detect.values.flatten()]).T
    bdf = pd.DataFrame(barray,columns=['HYBAS_ID','FirstFirms','DaysFromFirstDetection'])
    bdf.set_index('HYBAS_ID',inplace=True)
    bdf['FirstRunDate'] = pd.to_datetime(date.today())
    
    #print(bdf.head())
    
    print('Basins for recent FIRMS run:',len(bdf))
    # update existing basin data
    intersection = list(set(bdf.index)&set(b_data.index)) # basins in both DFs
    print('Already Present in Records:',len(intersection))
    b_data_cols = ['RecentFirms','LastRunDate','DaysFromLastDetection'] # cols to write to
    recent_cols = ['FirstFirms','FirstRunDate','DaysFromFirstDetection'] # cols to use for writing

    # update data for basins already listed in record. Changes are:
    # 'FirstFirms' --> 'RecentFirms', 'FirstRunDate' --> 'LastRunDate'
    # 'DaysFromFirstDetection' --> 'DaysFromLastDetection'
    b_data.loc[intersection,b_data_cols] = bdf.loc[intersection,recent_cols].values # update!

    # add in remaining new data
    concat_data = bdf[~bdf.index.isin(intersection)] # basins not in intersection
    print('To be added:',len(concat_data))
    
    new_data = pd.concat([b_data, concat_data])
    
    print('Previous record length:',len(b_data))
    print('New record length:',len(new_data))
    
    # final step, aggregate list of all basins to be run
    all_basins = list(set(np.concatenate([np.array(f_basins),
                                np.array(cc_basins)])))
    # and save
    new_data.to_csv(b_data_path,header=True,index=True)
    
    run_data = new_data.loc[all_basins]
    run_data.to_csv(r_data_path,header=True,index=True)
    return 


def update_dnbr_data(file_path,basin_records_path):
    print('dNBR INFO UPDATE \n')
    data = pd.read_csv(file_path,index_col='HYBAS_ID') # read asset file
    b_data = pd.read_csv(basin_records_path,index_col='HYBAS_ID')
    print('# of current records:', len(b_data))

    intersection = list(set(data.index)&set(b_data.index)) 
    print('# to be updated:', len(intersection))
    
    update_cols = ['Mean_dLAI','Mean_dNBR', 'Mean_dNBR_cnt','Slope_Burn_Abv_23','SlopeBurnAreaRatio','Visibility']
    b_data.loc[intersection,update_cols] = data.loc[intersection,update_cols]
    b_data.loc[intersection,'LastBurnUpdate'] = pd.to_datetime(date.today())
    b_data.to_csv(basin_records_path,header=True,index=True)

    return

def monitor_task(task,sleep_seconds=25):
    count = 0
    state = task.status()['state']
    
    while state != 'COMPLETED':
        # handle exceptions
        if state=='FAILED':
            print(task.status()['error_message'])
            raise ValueError('Export Failed')
            
        if state=='CANCEL_REQUESTED' or state=='Canceled':
            raise ValueError('Export Canceled')
    
        # convert seconds to hrs
        runtime = ((sleep_seconds*count)/60/60)
        runtime_f = "{:.3}".format(runtime) # 3 sig figs
    
        print('Current Runtime:',runtime_f,'hrs','('+state+')' )
        state = task.status()['state']
    
        sleep(sleep_seconds)
        count += 1
    print('Export Completed (Runtime: '+runtime_f+'hrs)')
    return

def update_task_list(task_records_path,return_list=True):
    task_list = pd.read_csv(task_records_path,index_col='File')
    pend = task_list[
        (task_list['Status']=='PENDING') | (task_list['Status']=='RUNNING')]
    p_idx = pend.index
    
    download_list = []
    
    for p in p_idx:
        update = ee.data.getOperation(task_list.loc[p]['OpID'])
        meta = update['metadata']
        status = meta['state']
        start = meta['startTime']
        
        if (status == 'RUNNING') | (status == 'PENDING'):
            print(p, 'still running...\n')
            print('full status:',meta,'\n')
            task_list.loc[p,'StartTime'] = start
            continue
            
        if (status == 'CANCELED')|(status == 'CANCELLING'):
            print(p, 'canceled.')
            task_list.loc[p,'Status'] = 'CANCELED'
            task_list.loc[p,'StartTime'] = start
   
            
        if status=='FAILED':
            error = update['error']['message']
            print(p,'failed:', error)
            task_list.loc[p,'Status'] = status
            task_list.loc[p,'StartTime'] = start
            task_list.loc[p,'EndTime'] = meta['endTime']
            task_list.loc[p,'Message'] = error
            continue
            
        if status == 'SUCCEEDED':
            print(p, 'succeeded.')
            end_t = pd.to_datetime(meta['endTime'])
            start_t = pd.to_datetime(start)
            diff = (end_t - start_t) / np.timedelta64(1,'h')
            print('total runtime:', diff, 'hours')
            task_list.loc[p,'Status'] = status
            task_list.loc[p,'StartTime'] = start
            task_list.loc[p,'EndTime'] = meta['endTime']
            download_list.append(p)
            continue

    if return_list:
        return download_list, task_list.copy()
    else:
        return task_list.copy()
    
    
def get_lists(main_path):
    run_data_path = os.path.join(main_path,'ref_data','run_data.csv') 
    r_data = pd.read_csv(run_data_path,
                     index_col='HYBAS_ID')
    print('basins read in for run:', len(r_data))
    basins = r_data.index.values.tolist()
    run_dates = pd.to_datetime(r_data['FirstRunDate']).astype(str).tolist()
    detect = r_data['DaysFromFirstDetection'].values.tolist()

    assert len(basins)==len(run_dates)==len(detect)
    return basins, run_dates, detect