Functions to create an updating database of burned areas following wildfires and their surface properties for use in a PFDF prediction model. Two steps are required for this process to function properly: 

1. Your environment should already have the Google Earth Engine Python API configured (installation info [here](https://developers.google.com/earth-engine/guides/python_install))

2. You should have a [NASA EarthData Account](https://urs.earthdata.nasa.gov/) and retrieve an access token for data downloads. This allows you to access FIRMS data for near real time (NRT) downloads. After making an account, go to https://nrt3.modaps.eosdis.nasa.gov/, log in by clicking on the "Profile" menu tab (upper right), select "Download Tokens" from this same menu tab, and finally "Generate a Token". Copy and paste this token into ```token.txt```, saved in ```pfdf\ref_data\```. 

After downloading, run ```setup.py``` with the following command: 

```python setup.py --filepath C:\path\to\repo\ --archive_path C:\path\to\archive\data\folder\ --model_output_path C:\path\to\model\output\folder\ --gee_username yourGEEUsername```


For example, if the repo is saved in the folder ```gitlab_repos```, the ```--filepath``` parameter should be ```C:\...\gitlab_repos\pfdf\```.

Note: the ```--archive_path``` and ```--model_output_path``` commands are OPTIONAL.  ```--archive_path``` specifies the directory for storing the FIRMS and downloaded Earth Engine files. ```--model_output_path``` is the directory for the ```.geojson``` xgboost model output files. The containing folder names and paths need to be created by the user. Otherwise, these fields can be left blank and ```setup.py``` will use the containing ```pfdf``` repo folder to house them accordingly. 

Then run:

```python request.py --filepath C:\path\to\repo\ --firms_path C:\path\to\firms\containing\folder\firms\all_firms```

```python gee_export.py --filepath C:\path\to\repo\ --gee_username yourGEEUsername --firms_path C:\path\to\firms\containing\folder\firms\all_firms```

(note: these must be run one after another and _not_ in parallel)

These two scripts download daily global FIRMS data for the past 60 days and perform a spatial join with a global set of watersheds. This is a non-trival task. After the join, burned basins are sent to GEE to calculate burn area and slope characteristics. Running these two scripts will take several minutes to complete with additional time needed for the basin calculations to be performed on GEE servers. Running ```query.py``` checks the status of the GEE operations and will download the corresponding file when it is complete. For an initial run of ~20 basins, this process can take as little as 1-2min. Run ```query.py``` with the following syntax: 

```python query.py --filepath C:\path\to\repo\ --asset_path C:\path\to\archive\data\folder\ --gee_username yourGEEUsername```

If ```--archive_path``` was specified in ```setup.py```, that path should be used for ```--asset_path```. Otherwise this field is the same as ```--filepath```. 

All summary data are saved in ```basin_records.csv``` in the ```ref_data``` folder. Full export data are saved in ```asset_downloads```


# An example workflow looks like this: 


```python C:\Users\user\Documents\pfdf\scripts\setup.py --filepath C:\Users\user\Documents\pfdf\ --gee_username myusername```

```python C:\Users\user\Documents\pfdf\scripts\request.py --filepath C:\Users\user\Documents\pfdf\ --firms_path C:\Users\user\Documents\pfdf\firms\all_firms```

```python C:\Users\user\Documents\pfdf\scripts\gee_export.py --filepath C:\Users\user\Documents\pfdf\ --gee_username myusername --firms_path C:\Users\user\Documents\pfdf\firms\all_firms```

```python C:\Users\user\Documents\pfdf\scripts\query.py --filepath C:\Users\user\Documents\pfdf\ --asset_path C:\Users\user\Documents\pfdf\ --gee_username myusername```


And for rerunning all basins at once: 
```python C:\Users\user\Documents\pfdf\scripts\gee_export_all.py --filepath C:\Users\user\Documents\pfdf\ --gee_username myusername```
