import os    
import shutil
import glob
import argparse
import firms_request
from firms_basin_join import join
import stac_functions
import calc_basin_stats
from shapely import wkt
from pull_imerg import pull_imerg
import geopandas as gpd
import pandas as pd
from datetime import date
import numpy as np
from fiona.errors import DriverError
import xgboost as xgb
import xarray as xr


def directory_setup(archive_path, model_output_path, home_path):
    """Create folders to store data"""
    os.makedirs(archive_path, exist_ok=True)
    model_dir = os.path.join(model_output_path, "model_outputs")
    os.makedirs(model_dir, exist_ok=True)
    os.makedirs(os.path.join(home_path,'recover_files'),exist_ok=True)
    return True


def flag_recent_failed_runs(run_data):
    '''Function that separates out data that either have never been previously run, 
       or that had prior issues which now we presume may be resolved. This includes
       data that simply had an API error, or that had visibility/imagery issues 
       where enough time has passed such that new imagery has likely been acquired.'''
    
    r_data = run_data.copy()
    # check how much time has passed since that last run
    r_data['TimeDiff'] = (pd.to_datetime(date.today()) - pd.to_datetime(r_data['LastRunDate'])) / np.timedelta64(1,'D')
    
    # isolate those with either no notes (i.e. haven't been run yet) or that haven't been run in the last 3 days
    run_idx = r_data[(r_data['Notes'].isna())|(r_data['TimeDiff']>3)].index.tolist()
    
    # Now pull data that's been recently run, and check if they had errors unrelated to avalable imagery
    recent_runs = r_data[(r_data['TimeDiff']<3)]
    recent_runs = recent_runs[~recent_runs['Notes'].isna()]
    
    # possible errors we want to skip over and run at a later date
    err1 = 'No imagery available' 
    err2 = 'NoneType' # <-- Translates to 'No imagery available'
    err3 = '0d-boolean' # <-- same as above
    err4 = 'Vis Level' # means imagery is available but it's not clear enough! 
    
    for idx, row in recent_runs.iterrows():
        # if any of the below errors are present, skip it - there hasn't been enough time to try again
        if (err1 in row['Notes']) or (err2 in row['Notes']) or (err3 in row['Notes']) or (err4 in row['Notes']):
            continue
        else: # otherwise if the note is different, we also want to rerun these!
            run_idx.append(idx)
    return run_idx # return all data ready to be rerun


def get_model(file_path: str, threads=1):
    """Open trained model"""
    model = xgb.Booster()
    model.load_model(file_path)
    model.set_param("nthread", threads)
    return model


def zonal_max(polygon, grid: xr.DataArray):
    if polygon.type == "Polygon":
        x, y = polygon.exterior.xy
    else:
        coords = [p.exterior.xy for p in polygon.geoms]
        x = np.concatenate([a[0] for a in coords])
        y = np.concatenate([a[1] for a in coords])
    lat = xr.DataArray(y, dims="points")
    lon = xr.DataArray(x, dims="points")
    grid_values = grid.sel(lon=lon, lat=lat, method="nearest")
    m = grid_values.max().values
    return m


def workflow(bbox: list, writer=print): 
    writer("Opening basin data")
    sheds = gpd.read_file(f'{home_path}/ref_data/sheds.geojson',bbox=bbox)
    writer("Opened basin data")
    writer("Downloading FIRMS")
    firms_request.download_firms(
        home_path=home_path,
        firms_path=firms_path,
        show_progress=False,
        overwrite=False,
    )
    writer('FIRMS downloaded\n')

    # Get most recent to-do list - might be empty if everything ran sucessfully
    # This list represents previous runs which failed that need a rerun now.
    try: 
        run_data = pd.read_csv(f'{home_path}/ref_data/BasinToDoList.csv',index_col='HYBAS_ID')
        run_from_scratch = False
                    
    except FileNotFoundError:
        writer('No run data file found! The data from this run will be saved to a new file.')
        run_from_scratch = True
    
    # Do basin join and get all hybas ids to run.
    new_ids = join(
        firms_path=firms_path,
        sheds=sheds,
        date_range=120,
        count_thresh=50,
        bbox=bbox,
    )
    
    if run_from_scratch: # if file not found or all run data finshed, set this most recent batch as it
        run_data = new_ids
    else: # otherwise load in the existing basins that need to run
        run_data['geometry'] = run_data.geometry.apply(wkt.loads) # instantiate geometry
        run_data = gpd.GeoDataFrame(run_data, geometry='geometry',crs='EPSG:6933')        
        # join with the new basins!
        run_data = pd.concat([run_data, new_ids[~new_ids.index.isin(run_data.index)]]) # only data that aren't already present
        
    # also get most recent set of catchments used for model runs
    try: 
        current_run_basins = gpd.read_file(f'{home_path}/ref_data/BasinsToRun.geojson')
    except (FileNotFoundError, DriverError): 
        current_run_basins = []
    
    if len(current_run_basins): # start working with the current file which houses all data for model input
        
        current_run_basins.set_index('HYBAS_ID',inplace=True)
        
        # if applicable: recover mid-run data if last run was interrupted     
        recover_files = glob.glob(os.path.join(home_path,'recover_files','*.geojson')) # list of available files
        
        if len(recover_files):
            writer('RECOVERED FILES:', recover_files)
            recovered_data = pd.concat([gpd.read_file(f) for f in recover_files]) # merge into single gdf
            recovered_data.set_index('HYBAS_ID',inplace=True)
            
            # flag the indices of the basins that got rerun
            basin_overlap = list(set(recovered_data.index)&set(current_run_basins.index)) 
            
            if basin_overlap: # update based on rerun!
                writer(f'Updating {len(basin_overlap)} IDs...')
                current_run_basins.drop(basin_overlap,inplace=True) # drop the overlapping indices from the main dataset

            # join all data, including that which were rerun
            current_run_basins = pd.concat([current_run_basins, recovered_data])
    
        # After that, look at the newly flagged data and check to see if any have already be run...
        recently_flagged_overlap = list(set(run_data.index)&set(current_run_basins.index))
    
        # If we already have data, don't run them again.     
        if recently_flagged_overlap:
            writer(f'Dropping {len(recently_flagged_overlap)} overlapping IDs...')
            run_data = run_data.drop(recently_flagged_overlap)
        
        
        # BUT: look at previous sucessful runs, and determine if it's time for dNBR reassessment.
        current_run_basins['dNBRCalcDiff'] = (pd.to_datetime(date.today()) - pd.to_datetime(current_run_basins['dnbrCalcDate'])) / np.timedelta64(1,'M')
        basins_for_rerun = current_run_basins[current_run_basins['dNBRCalcDiff']>1].copy()
        
        if len(basins_for_rerun):
            
            basins_for_rerun.index = basins_for_rerun.index.astype(int)
            sheds_c = sheds.copy().set_index('HYBAS_ID')
            rerun_sheds = sheds_c.loc[basins_for_rerun.index.unique()]
        
            all_rerun_info = gpd.GeoDataFrame(index=rerun_sheds.index.unique(),
                                              geometry=rerun_sheds.geometry,
                                              crs='EPSG:4326')
            
            cols = ['FireStart','FireEnd','PreOffset','PostOffset']
            all_rerun_info[cols] = basins_for_rerun.groupby(basins_for_rerun.index)[cols].first()
            all_rerun_info.to_crs('EPSG:6933',inplace=True)
        
            run_data = pd.concat([run_data,all_rerun_info],verify_integrity=True)
            writer(f'Adding {len(rerun_sheds)} basins for dNBR rerun...')
            
    if run_from_scratch:
        run_idx = run_data.index
    else:
        run_idx = flag_recent_failed_runs(run_data)    
        
    run_data['FireStart'] = run_data['FireStart'].apply(lambda x: pd.to_datetime(x).strftime('%Y-%m-%d'))
    run_data['FireEnd'] = run_data['FireEnd'].apply(lambda x: pd.to_datetime(x).strftime('%Y-%m-%d'))
    
    writer(f'Clearing {len(run_idx)} basins for most recent run..')
    run_data.loc[run_idx,'LastRunDate'] = pd.Timestamp.now()
    run_data.to_csv(f'{home_path}/ref_data/BasinToDoList.csv',header=True,index=True)
    
    run_data = run_data.loc[run_idx]
    
    # get API creds for HLS access
    creds = stac_functions.get_lpdaac_creds()
    
    # pull out dnbr, slope data!
    basins_all = calc_basin_stats.run_basins(basins=run_data, home_path=home_path, creds=creds)
        
    if len(basins_all): # workflow and any newly acquired data
        # set index vals for both this most recent run and all other model data
        basins_all.set_index('HYBAS_ID',inplace=True)
        
        if len(current_run_basins)==0:
            current_run_basins = basins_all
        else: 
            overlap = list(set(basins_all.index)&set(current_run_basins.index)) # flag the indices of the basins that got rerun
            if overlap: # update based on rerun!
                writer(f'Updating {len(overlap)} IDs...')                
                current_run_basins.drop(overlap,inplace=True) # drop the overlapping indices from the main dataset
            # join all data, including that which were rerun
            current_run_basins = pd.concat([current_run_basins, basins_all])
    
    if len(current_run_basins)==0:
        writer('No burned area data available. Exiting...')
        return
    
    # convert to string for easier saving
    current_run_basins['dnbrCalcDate'] = current_run_basins['dnbrCalcDate'].astype(str)
    current_run_basins['LastRunDate'] = current_run_basins['LastRunDate'].astype(str)
    current_run_basins['FireStart'] = current_run_basins['FireStart'].astype(str)
    current_run_basins['FireEnd'] = current_run_basins['FireEnd'].astype(str)
    
    dnbr_max = current_run_basins.groupby(current_run_basins.index)['Mean_dNBR'].max()  
    drop_idx = dnbr_max[dnbr_max.isna()].index
    writer(f'Dropping {len(drop_idx)} basins which no longer meet dNBR threshold...')
    current_run_basins.drop(drop_idx,axis=0,inplace=True)
    
    current_run_basins.to_file(f'{home_path}/ref_data/BasinsToRun.geojson',header=True,index=True)
    current_run_basins.to_file(f'{home_path}/ref_data/BasinsToRun_BACKUP.geojson',header=True,index=True)
    writer('all basins run\n')

    shutil.rmtree(os.path.join(home_path,'recover_files'))
    writer('Removing temp files...')
    
    
    writer('Pulling IMERG...\n')
    rain_max = pull_imerg(
        home_path=home_path,
        output_path=output_path,
        west=bbox[0],
        south=bbox[1],
        east=bbox[2],
        north=bbox[3],
    )
    
    basins = current_run_basins.reset_index()
    basins = basins.loc[basins["PropHM23"].dropna().index]
    basins["SlopeBurnAreaRatio"] = basins["PropHM23"]
    basins["MaxRain"] = np.array([zonal_max(b, rain_max) for b in basins["geometry"]])
    basins["MaxRain"] = pd.to_numeric(basins["MaxRain"])

    model = get_model(f"{home_path}/ref_data/model.json")
    cols = ["Mean_dNBR", "SlopeBurnAreaRatio", "MaxRain", "MedianSlope", "AreaSqKm"]
    predictors = xgb.DMatrix(basins[cols])

    basins["p_debris_flow"] = model.predict(predictors)
    basins_copy = basins.copy()

    basins_copy = basins_copy[
        (basins_copy["MaxRain"] >= 1)
        & (basins_copy["p_debris_flow"] >= 0.05)
        & (basins_copy["MedianSlope"] >= 10)
    ]

    out_cols = ["Mean_dNBR", "SlopeBurnAreaRatio", "MaxRain", "MedianSlope", "AreaSqKm", "p_debris_flow", "geometry"]
    print("Output shape:", basins_copy[out_cols].shape)
    json_name = os.path.join(output_path, "model_outputs", (rain_max.attrs['start'] + ".geojson"))
    basins_copy[out_cols].to_file(json_name, driver="GeoJSON")
    writer('Run finished.\n')
    return basins_copy

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Global post-fire debris flow model")
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="LHASA version 2.0.0",
    )
    parser.add_argument(
        "-p",
        "--path",
        default=os.getcwd(),
        help="file path to project folder.",
    )
    parser.add_argument(
        "-ap",
        "--archive_path",
        help=(
            "Location of output files including FIRMS data and asset downloads"
            ". Optional with default to the main directory if unspecified."
        ),
    )
    parser.add_argument(
        "-op",
        "--output_path",
        help=(
            "Location of geojson output files from model. "
            "Optional with default to the main directory if unspecified"
        ),
    )
    parser.add_argument(
        "-N",
        "--north",
        type=float,
        default=60.0,
        help="maximum latitude (WGS84)",
    )
    parser.add_argument(
        "-S",
        "--south",
        type=float,
        default=-60.0,
        help="minimum latitude (WGS84)",
    )
    parser.add_argument(
        "-E",
        "--east",
        type=float,
        default=180.0,
        help="maximum longitude (WGS84)",
    )
    parser.add_argument(
        "-W",
        "--west",
        type=float,
        default=-180.0,
        help="minimum longitude (WGS84)",
    )
    
    args = parser.parse_args()

    if args.west > args.east:
        raise ValueError(
            (
                f"West longitude ({args.west}) cannot be greater than East "
                f"longitude ({args.east})"
            )
        )
    if args.south > args.north:
        raise ValueError(
            (
                f"South latitude ({args.south}) cannot be greater than North "
                f"latitude ({args.north})"
            )
        )

    home_path = args.path
    if args.archive_path:
        firms_path = os.path.join(args.archive_path, "firms", "all_firms")
    else:
        firms_path = os.path.join(home_path, "firms", "all_firms")
    if args.output_path:
        output_path = args.output_path
    else:
        output_path = home_path
    
    directory_setup(firms_path, output_path, home_path)
    bbox = [args.west, args.south, args.east, args.north]
    workflow(bbox)
