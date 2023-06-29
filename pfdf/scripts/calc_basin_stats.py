'''Helper Function(s)!'''
import numpy as np
from geocube.api.core import make_geocube
import pandas as pd
import stac_functions
from xrspatial import slope
import make_catchments
from datetime import date
import os

def update_basin_log(home_path,basin_info,remove=False):
    
    # open file
    df = pd.read_csv(f'{home_path}/ref_data/BasinToDoList.csv',index_col='HYBAS_ID')
    
    keys = list(basin_info.keys())
    values = list(basin_info.values())
    
    hybas_id = values[0] # always the first entry!
    replace_keys = keys[1:] # cols to update
    replace_vals = values[1:] # the values to use in each col
    
    if remove: 
        df.drop(hybas_id,axis=0,inplace=True) # drop row
    
    else: 
        df.loc[hybas_id,replace_keys] = replace_vals # update row
     
    df.to_csv(f'{home_path}/ref_data/BasinToDoList.csv',index=True,header=True)

    return
    

def basin_zonal_stats(basins,slope,dnbr):
    
    orders = basins['Order'].sort_values().unique() 
    all_merge_df = pd.DataFrame() # this will hold the data for all the basins regardless of the order
    
    for order in orders:
        
        single_order_slice = basins[basins['Order']==order].copy()
        
        # convert all catchments for that order to a raster. each catchmnt has a unique index val.
        so_cube = make_geocube(single_order_slice,measurements=['Order','Index'],like=slope)
        
        # align slope and dnbr rasters...
        slope_re = slope.rio.reproject_match(so_cube)
        slope_re = slope_re.assign_coords({"x": so_cube.x,
                                           "y": so_cube.y})

        dnbr_re = dnbr.rio.reproject_match(so_cube)
        dnbr_re = dnbr_re.assign_coords({"x": so_cube.x,
                                         "y": so_cube.y})
        
        # add these as "bands" to the raster representing the basins
        so_cube["slope"] = (slope_re.dims, slope_re.values, slope_re.attrs, slope_re.encoding)
        so_cube["dnbr"] = (dnbr_re.dims, dnbr_re.values, dnbr_re.attrs, dnbr_re.encoding)
        
        # convert to dataframe
        so_df = so_cube.to_dataframe().reset_index().dropna()
        
        # set conditions for the zonal stats 
        burn_cond = so_df['dnbr']>0.1
        hi_burn_cond = so_df['dnbr']>0.4
        slope_cond = so_df['slope']>23
        
        # do groupby based on index. 
        # because each catchment has a unique index, this summarizes all the dnbr and slope pixels that are "contained" in this index.
        mean_dnbr = so_df[burn_cond&slope_cond].groupby('Index')['dnbr'].mean().rename('Mean_dNBR')
        prophm23 = (so_df[hi_burn_cond&slope_cond].groupby('Index')['dnbr'].count() / so_df[slope_cond].groupby('Index')['dnbr'].count()).rename('PropHM23')
        slope_count = so_df[slope_cond].groupby('Index')['slope'].count().rename('SlopeCount')
        slope_median = so_df.groupby('Index')['slope'].median().rename('MedianSlope')
        
        merge = pd.concat([mean_dnbr,prophm23,slope_count,slope_median],axis=1) # index here is "Index"
        all_merge_df = pd.concat([all_merge_df,merge],axis=0)
        
    # print('Merge len:',len(all_merge_df))
    # print('Basin len:',len(basins))
    basins_j = basins.join(all_merge_df,how='inner',on='Index') # finally join back with the original data
    
    return basins_j 


def basin_dnbr_pull(basin, fire_start, fire_end, home_path, 
                    pre_offset=1, post_offset=2, 
                    cloud_cover_threshold=20,creds=None,
                    mask=True,writer=print):
        
        '''Acts as a pseudo-wrapper and runs stac_functions.calc_dnbr() 
           for a specific basin and its respective attributes. 
           Handles the various exceptions that may occur when used in 
           the context of the PFDF model. Similarly to how calc_dnbr()
           returns the creds for optional reuse, this same process occurs here.
           For any unsuccessful pull of dNBR data, this function returns None.'''
        
        hybas_id = basin['HYBAS_ID'].iloc[0]
        
        try: 
            dnbr, creds = stac_functions.calc_dnbr(basin, # this is the "generic" dnbr function
                                                   fire_start=fire_start, 
                                                   fire_end=fire_end,
                                                   pre_offset=pre_offset,
                                                   post_offset=post_offset, 
                                                   cloud_cover_threshold=cloud_cover_threshold,
                                                   creds=creds,
                                                   mask=mask)
            
            ### Condition 1 ###
            #### No imagery available for the periods specified, so calc_dnbr() returns None ###
            if dnbr is None:
                
                writer('No imagery available for listed offset periods.')
                basin_info = {'HYBAS_ID': hybas_id,'PostOffset': post_offset,'Notes':'No imagery available.'}
                update_basin_log(home_path,basin_info) # do logging
                
                return dnbr, creds # dnbr is None here
            
            # calculate visibility: count number of unmasked (i.e. not nan) pixels / all pixels 
            visibility = float(dnbr.notnull().sum() / (dnbr.shape[0]*dnbr.shape[1])) 
        
            ### Condition 2 ###
            #### Imagery IS available, but too cloudy ###
            if (visibility < 0.95) & (mask==True):
                writer('dNBR data fall below visibility threshold.')
                basin_info = {'HYBAS_ID': hybas_id,'PostOffset': post_offset,'Notes': f'Vis Level: {visibility}'}
                update_basin_log(home_path,basin_info)
                
                dnbr = None # manually set to None to indicate error
                return dnbr, creds 
            
            else: 
                return dnbr, creds # only instance in which dnbr is not None!
        
        ### Condition 3 ###
        #### Some other error is occuring, perhaps with the API response ###
        except Exception as e:
            
            writer(f'Other exception for HLS data occured: {e}')
            
            basin_info = {'HYBAS_ID': hybas_id,'PostOffset': post_offset,'Notes':f'Error: {e}'}
            update_basin_log(home_path,basin_info)
            dnbr = None
            return dnbr, creds


########################################################################
########################################################################

                        # Main Functions!
    
########################################################################
########################################################################  

def run_basins(basins, home_path, writer=print, creds=None, mask=True, 
               cloud_cover_threshold=20):
    
    basins = basins.reset_index() # remove HYBAS_IDs from index and back into DF 
    pc_stac_catalog = stac_functions.pc_catalog()
    basins_all = []
    
    for _id in basins['HYBAS_ID'].unique(): # look at each basin
                
        writer(f'Pulling STAC data for id: {_id}')

        basin = basins.loc[basins['HYBAS_ID']==_id].copy()
        
        # all original parameters
        start = str(basin['FireStart'].iloc[0])
        end = str(basin['FireEnd'].iloc[0])
        today = pd.Timestamp.now()
        
        # calc post-fire offset
        timedelta = np.ceil((today - pd.to_datetime(end)) / np.timedelta64(1,'M')) # round up to next month!

        # get original pre-fire offset
        pre_offset = float(basin['PreOffset'].iloc[0])
        
        # attempt dnbr pull
        dnbr, creds = basin_dnbr_pull(basin, 
                                      home_path=home_path, 
                                      fire_start=start, 
                                      fire_end=end,
                                      pre_offset=pre_offset,
                                      post_offset=timedelta, # months to today
                                      cloud_cover_threshold=cloud_cover_threshold,
                                      creds=creds,
                                      mask=mask) 
        
        if dnbr is not None: # if it succeeds, proceed with next steps
            
            writer(f'dNBR acquired with {timedelta} month offset.')    
            
            ### Next step: get DEM data! ###
            try: 
                dem = stac_functions.get_elevation(basin,catalog=pc_stac_catalog)
                
            except: # Rarely this will fail. Handle it below. 
                    writer('Could not retrieve DEM data. Retrying...')
                    pc_stac_catalog = stac_functions.pc_catalog() 
                    dem = stac_functions.get_elevation(basin,catalog=pc_stac_catalog)

            
            # Next step, get slope data, generate catchments!
            try: 
                slope_raster = slope(dem)
                dem.rio.to_raster("temp_elevation.tif")
                
                basins_c, branches = make_catchments.generate_catchments('temp_elevation.tif',
                                                                          basin,
                                                                          acc_thresh=100,so_filter=4)
            except Exception as e: # Again, this will rarely fail. Future work will improve on this. 
                    writer(f'Exception for slope/catchment generation occured: {e}')
                    basin_info = {'HYBAS_ID': _id,'Notes':f'Slope/Catchment Error: {e}'}
                    update_basin_log(home_path,basin_info)
                    continue
            
            # finally calc the zonal stats!
            try: basins_zonal = basin_zonal_stats(basins_c,slope_raster,dnbr)
            except Exception as e: # Again, this will rarely fail. Future work will improve on this. 
                    writer(f'Exception for zonal stats occured: {e}')
                    basin_info = {'HYBAS_ID': _id,'Notes':f'Zonal Stats Error: {e}'}
                    update_basin_log(home_path,basin_info)
                    continue
                
            # add some extra helpful info
            basins_zonal['FireStart'] = start
            basins_zonal['FireEnd'] = end
            basins_zonal['PreOffset'] = pre_offset
            basins_zonal['PostOffset'] = timedelta
            basins_zonal['dnbrCalcDate'] = pd.to_datetime(date.today()).strftime('%Y-%m-%d')
            basins_zonal['LastRunDate'] = pd.to_datetime(date.today()).strftime('%Y-%m-%d')
            # basins_zonal['dNBRCalcDiff'] = ((pd.to_datetime(date.today()) - pd.to_datetime(basins_zonal['dnbrCalcDate']))
            #                                  / np.timedelta64(1,'M'))
            basins_zonal.to_file(f'{os.path.join(home_path,"recover_files")}/{_id}_run.geojson',header=True,index=True)
            basins_all.append(basins_zonal) # add to growing list
            
            writer(f'{_id} sucessfully run and removed from list.\n')
            update_basin_log(home_path,{'HYBAS_ID':_id},remove=True)
    
    if len(basins_all):        
        basins_all = pd.concat(basins_all,ignore_index=True)
        
    return basins_all # this can be an empty list
