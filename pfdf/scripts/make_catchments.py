from pysheds.grid import Grid
import geopandas as gpd
from shapely import geometry, ops
import numpy as np
from pyproj import Geod
from shapely.errors import ShapelyDeprecationWarning
import warnings

# filter warnings for now - code will need updating for Shapely 2.0+
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning) 

# global variables
geod = Geod(ellps="WGS84")

######################################################################
#                      Helper Functions
######################################################################

def calc_attributes(row,dem,profiles,so):
    '''Takes a stream segment as input and assigns various 
       elements related to its geometry and stream order. 
       Topographic metrics such as relief and slope are 
       approximations only.'''
    
    index = row.name 
    row['Index'] = row.name
    row['Length'] = geod.geometry_length(row['geometry'])
    row['Relief'] = int(np.ptp(dem.flat[profiles[index]]))
    
    row['Order'] = int(np.min(so.flat[profiles[index]]))
    # note: 'so.flat[profiles[index]]' returns a 1d array of length n
    # where all values are equal to the stream order for that profile.
    # so if it is a 1st order channel, the array has all 1's. 
    
    row['Slope'] = np.rad2deg(np.arctan2(row['Relief'],row['Length']))
    
    return row 

def make_connections(profiles, connections):
    '''creates full connections map for stream network.'''
    
    profile_list = []
    connection_list = []
    
    for i in range(0,len(profiles)): # loop through each segment profile
        
        index = i
        chain = []
        connection_map = []
        
        # create full profile and connection map
        while True:
            chain.extend(profiles[index]) # add segment
            connection_map.extend([index]) # add that segment id
            next_index = connections[index] # get idx of the connecting segment
            
            # if next idx is the same as before, then the full profile is done
            if index == next_index:  
                break
            index = next_index # otherwise go to the next connecting segment
        
        profile_list.append(chain) # add full profile
        connection_list.append(connection_map) # add full connection mapping
    
    return profile_list, connection_list

def make_profile(row,so,coords,profile_list,connection_list):
    index = row.name # stream profile index 
    row_order = row['Order']
    profile = profile_list[index] # grab full profile for segment
        
    so_list = so.flat[profile].tolist() # get stream orders along full profile
    so_list.reverse() # reversing to perform the next steps easier...
    
    final_idx = len(so_list) - so_list.index(row_order) - 1
    # ^^^ get last instance of this stream order in the full profile. 
    
    row['LocalPP'] = np.flip(coords[profile[final_idx]]) 
    # ^^^ use idx to get the appropriate coords for the pour point.
    # np.flip() puts it in the right order
    row['LocalPP_X'] = row['LocalPP'][0]
    row['LocalPP_Y'] = row['LocalPP'][1]
    
    
    row['Profile'] = profile_list[index] # add full profile
    row['Chain'] = connection_list[index] # add full chain mapping    
    
    row['Final_SO'] = so.flat[profile][-1]
    # idx of final segment along profile
    row['Final_Chain_Val'] = connection_list[index][-1] 
    return row

def make_basin(row,grid,dem,fdir,
               routing,algorithm):
    '''Calculate basin geometry via pour point 
       derived from above functions.'''
    x, y = row['LocalPP_X'], row['LocalPP_Y']
    c = grid.catchment(x=x,y=y,fdir=fdir,routing=routing,algorithm=algorithm)
    grid.clip_to(c)
    catchment_polygon = ops.unary_union([geometry.shape(shape) 
                                     for shape, value in grid.polygonize()])
    grid.viewfinder = dem.viewfinder
    
    row['BasinGeo']  = catchment_polygon
    return row

######################################################################
#                      Main Function
######################################################################

def generate_catchments(path,basin,acc_thresh=100,so_filter=3,
                        routing='d8',algorithm='iterative'):
    '''Full workflow integrating the above functions.
       Process is as follows: 
       
       1. Read and process DEM
       2. Create stream network and connection map
       3. Determine relevent pour point locations
       4. Generate all catchments.
       5. Return gdf of all catchment data and stream network.
       
       '''
    print('Reading DEM...')
    grid = Grid.from_raster(path)
    dem = grid.read_raster(path,band=1)
    pit_filled_dem = grid.fill_pits(dem)
    # Fill depressions in DEM
    flooded_dem = grid.fill_depressions(pit_filled_dem)
    
    # Resolve flats in DEM
    inflated_dem = grid.resolve_flats(flooded_dem)
    # Specify directional mapping
    dirmap = (64, 128, 1, 2, 4, 8, 16, 32)
    
    # Compute flow directions
    # -------------------------------------
    fdir = grid.flowdir(inflated_dem, dirmap=dirmap)
    
    # make flow accumulation raster
    print('Making flow accumulation map...')
    acc = grid.accumulation(fdir, dirmap=dirmap)
    # set mask
    acc_mask = (acc_thresh < acc )
    # calculate stream order
    so = grid.stream_order(fdir=fdir,mask=acc_mask)
    # update mask
    mask = acc_mask
    
    # make river network
    print('Generating stream network...')
    branches = grid.extract_river_network(fdir=fdir,mask=mask) # returns geojson
    # make main GeoDataFrame -- the indices here will match those of all profile
    # and connection lists that follow
    branch_gdf = gpd.GeoDataFrame.from_features(branches,crs='epsg:6933')
    branch_gdf['HYBAS_ID'] = int(basin['HYBAS_ID'].iloc[0])
    # generate profiles for each individual segment
    profiles, connections = grid.extract_profiles(fdir=fdir,mask=mask,include_endpoint=False)
    #print(branch_gdf)
    branch_gdf = branch_gdf.apply(lambda x: calc_attributes(x,dem,profiles,so),axis=1)
    coords = dem.coords
    
    profile_list, connection_list = make_connections(profiles, connections)
    
    # create all profile and connection lists w/ pour points
    print('Calculating pour points...')
    branch_gdf = branch_gdf.apply(lambda x: make_profile(x,so,coords,profile_list,connection_list),axis=1)
    
    branch_gdf_copy = branch_gdf.copy() # unfiltered copy to return at end
    
    # some of these sections will have duplicate orders and pour points - drop them
    unique = branch_gdf[['Order','LocalPP_X','LocalPP_Y','Final_Chain_Val']].drop_duplicates()
    
    branch_gdf = branch_gdf.loc[unique.index]
    
    # filter by stream order
    if so_filter:
        branch_gdf = branch_gdf[branch_gdf['Order']<=so_filter]
    
    print('Generating',len(branch_gdf),'catchments...')
    branch_gdf = branch_gdf.apply(lambda x: make_basin(x,grid,dem,fdir,routing,algorithm),axis=1)

    b_copy = branch_gdf.copy()
    b_copy.set_geometry('BasinGeo',inplace=True)
    b_copy.set_crs('epsg:6933',inplace=True)
    
    copy_for_area = b_copy.copy()
    copy_for_area.to_crs('epsg:6933',inplace=True)
    copy_for_area.geometry.area
    copy_for_area['AreaSqKm'] = copy_for_area.geometry.area  / 1000000

    b_copy['AreaSqKm'] = copy_for_area['AreaSqKm']
    
    basin_drop_cols = ['geometry','LocalPP'] #,'Profile','Chain']
    branch_drop_cols = ['LocalPP'] #,'Profile','Chain']
    
    b_copy.drop(basin_drop_cols,axis=1,inplace=True)
    branch_gdf_copy.drop(branch_drop_cols,axis=1,inplace=True)
    
    b_copy['Profile'] = b_copy['Profile'].astype(str)
    b_copy['Chain'] = b_copy['Chain'].astype(str)
    b_copy.rename({'BasinGeo':'geometry'},axis=1,inplace=True)
    b_copy.set_geometry('geometry',inplace=True)

    branch_gdf_copy['Profile'] = branch_gdf_copy['Profile'].astype(str)
    branch_gdf_copy['Chain'] = branch_gdf_copy['Chain'].astype(str)
    
    # final step, only filter to data within the actual basin
    basin.to_crs('EPSG:6933',inplace=True)
    
    b_copy = b_copy[b_copy.geometry.centroid.within(basin.unary_union)]
    branch_gdf_copy = branch_gdf_copy[branch_gdf_copy.geometry.centroid.within(basin.unary_union)]
    return b_copy, branch_gdf_copy
