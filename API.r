# 7-day ARI calculation v1.0
# 11/20/2018
# Thomas Stanley NASA GSFC/GESTAR/USRA
# Calculates a 7-day antecedent precipitation index

# Load R packages
library(raster)

# Set working directory
setwd('C:/LHASA')
# List all daily precipitation files
files <- list.files(path='IMERG', full.names=TRUE)

# Set antecedent precipitation window
API.window <- 7
# Set IDW exponent
exponent <- 2
# Calculate weights
w <- 1/seq(API.window,1)^exponent
# Iterate through every day, starting at the end of the 1st API window
for(day in API.window:length(files)){
	# Open files within window, including current day
	IMERG <- crop(stack(files[(day-API.window+1):day]), extent(-180,180,-50,50))/10
	# Calculate antecedent precipitation index
	API <- calc(w*IMERG,sum)/sum(w)
	# Save to disk
	writeRaster(API,filename=gsub('IMERG/','API/',files[day]))
}
