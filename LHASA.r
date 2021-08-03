# Landslide Hazard Assessment for Situational Awareness, version 1.1.1
# By Thomas Stanley USRA/GESTAR 2021-8-2
# NASA Goddard Space Flight Center
# Maps the potential for landslides by identifying which locations
# exceed thresholds for a 7-day antecedent precipitation index
# and a landslide susceptibility map.

# Load packages
library(raster)
# Set working directory
setwd('C:/LHASA')

# Open antecedent rainfall threshold file
# Note that the version of this file posted at https://github.com/nasa/LHASA
# is based on the use of daily geotiff files, which are in tenths of mm.
# To use with the netcdf version of IMERG, divde by 10. 
threshold <- crop(raster('ARI95.tif'), extent(-180, 180, -60, 60))
# Open the susceptibility map
susceptible <- crop(raster('global.tif'), extent(-180, 180, -60, 60))
# List antecedent rainfall files
files <- list.files(path='ARI', pattern='*.tif', full.names=TRUE)

# Iterate through all days in record
for(f in files){
	# Open antecedent rainfall index file for current date
	ARI <- raster(f)
	# Compare to ARI threshold
	wet <- ARI > threshold
	# Run decision tree model at resolution of susceptibility map
	moderate <- resample(wet, susceptible, method='ngb') & (susceptible > 2)
	high <- moderate & (susceptible > 4)
	nowcast <- moderate + high
	# Save outputs
	writeRaster(nowcast,filename=gsub('ARI/','nowcast/',f), datatype='INT1U')
}
