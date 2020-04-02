# 7-day ARI calculation v1.1
# 2020-4-2
# Thomas Stanley NASA GSFC/GESTAR/USRA
# Calculates a 7-day antecedent rainfall index

# Load R packages
library(raster)

# Set working directory
setwd('C:/LHASA')
# List all daily rainfall files
files <- list.files(path='IMERG', full.names=TRUE)

# Set antecedent rainfall window
API.window <- 7
# Set IDW exponent
exponent <- 2
# Calculate weights
w <- 1/seq(API.window, 1)^exponent
# Iterate through every day, starting at the end of the 1st API window
for(day in API.window:length(files)){
	# Open files within window, including current day
	IMERG <- crop(stack(files[(day - API.window + 1):day]), extent(-180, 180, -60, 60)) / 10
	# Calculate antecedent rainfall index
	API <- calc(w*IMERG,sum)/sum(w)
	# Save to disk
	writeRaster(API,filename=gsub('IMERG/','API/',files[day]))
}
