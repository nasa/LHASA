# Landslide Hazard Assessment for Situational Awareness, version 1.1
# By Thomas Stanley USRA/GESTAR 2020-4-2
# NASA Goddard Space Flight Center
# Maps the potential for landslides by identifying which locations
# exceed thresholds for a 7-day antecedent precipitation index
# and a landslide susceptibility map.

# Load packages
library(raster)
# Set working directory
setwd('C:/LHASA')

# Open antecedent precipitation threshold file
threshold <- crop(raster('ARI95.tif'), extent(-180, 180, -60, 60))
# Open the susceptibility map
susceptible <- crop(raster('global.tif'), extent(-180, 180, -60, 60))
# List antecedent precipitation files
files <- list.files(path='API', pattern='*.tif', full.names=TRUE)

# Iterate through all days in record
for(f in files){
	# Open antecedent precipitation index file for current date
	API <- raster(f)
	# Compare to API threshold
	wet <- API > threshold
	# Run decision tree model at resolution of susceptibility map
	moderate <- resample(wet, susceptible, method='ngb') & (susceptible > 2)
	high <- moderate & (susceptible > 4)
	# Save outputs
	writeRaster(moderate,filename=gsub('API/','moderate/',f), datatype='INT1U')
	writeRaster(high,filename=gsub('API/','high/',f), datatype='INT1U')
}
